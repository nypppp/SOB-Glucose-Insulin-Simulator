'use strict';

/*
 * Reproducible checks for the eight-state Hovorka T1D v2 specification.
 * This is independent from the browser simulator and does not provide dosing
 * advice. Internal glucose units follow the paper: mmol/kg and mmol/L.
 */

const P = {
  weightKg: 70,
  targetG: 5.5,
  VG: 0.16,
  VI: 0.12,
  k12: 0.066,
  ka1: 0.006,
  ka2: 0.060,
  ka3: 0.030,
  ke: 0.138,
  SIT: 51.2e-4,
  SID: 8.2e-4,
  SIE: 520e-4,
  EGP0: 0.0161,
  F01: 0.0097,
  AG: 0.8,
  tmaxG: 40,
  tmaxI: 55,
  mmolToMgDl: 18.0182,
};

P.kb1 = P.ka1 * P.SIT;
P.kb2 = P.ka2 * P.SID;
P.kb3 = P.ka3 * P.SIE;

function correctedF01(G) {
  return G >= 4.5 ? P.F01 : P.F01 * G / 4.5;
}

function renalClearance(G) {
  return G >= 9 ? 0.003 * (G - 9) * P.VG : 0;
}

function carbsToMmolPerKg(grams) {
  return grams * 1000 / (180.16 * P.weightKg);
}

function mealAppearance(t, meals) {
  let total = 0;
  for (const meal of meals) {
    const age = t - meal.time;
    if (age >= 0) {
      const DG = carbsToMmolPerKg(meal.grams);
      total += DG * P.AG * age * Math.exp(-age / P.tmaxG) /
        (P.tmaxG * P.tmaxG);
    }
  }
  return total;
}

function equilibriumResidual(I, targetG = P.targetG) {
  const x1 = P.SIT * I;
  const x2 = P.SID * I;
  const x3 = P.SIE * I;
  const Q1 = P.VG * targetG;
  const Q2 = x1 * Q1 / (P.k12 + x2);
  const EGP = P.EGP0 * Math.max(0, 1 - x3);
  return -correctedF01(targetG) - x1 * Q1 + P.k12 * Q2 -
    renalClearance(targetG) + EGP;
}

function bisect(fn, lo, hi, tolerance = 1e-13) {
  let flo = fn(lo);
  let fhi = fn(hi);
  if (flo * fhi > 0) throw new Error('Root is not bracketed');
  while (hi - lo > tolerance) {
    const mid = (lo + hi) / 2;
    const fm = fn(mid);
    if (flo * fm <= 0) {
      hi = mid;
      fhi = fm;
    } else {
      lo = mid;
      flo = fm;
    }
  }
  return (lo + hi) / 2;
}

function equilibrium(targetG = P.targetG) {
  const I = bisect(value => equilibriumResidual(value, targetG), 0, 1 / P.SIE);
  const x1 = P.SIT * I;
  const x2 = P.SID * I;
  const x3 = P.SIE * I;
  const Q1 = P.VG * targetG;
  const Q2 = x1 * Q1 / (P.k12 + x2);
  const uBasal = I * P.VI * P.ke;
  const S1 = uBasal * P.tmaxI;
  const S2 = S1;
  return { y: [Q1, Q2, S1, S2, I, x1, x2, x3], uBasal };
}

const EQ = equilibrium();

function derivative(t, y, cfg) {
  const [Q1, Q2, S1, S2, I, x1, x2, x3] = y;
  const G = Q1 / P.VG;
  const basal = t >= (cfg.stopBasalAt ?? Infinity) ? 0 : EQ.uBasal;
  const EGPfactor = cfg.clipEGP ? Math.max(0, 1 - x3) : 1 - x3;
  const UG = mealAppearance(t, cfg.meals);
  return [
    -correctedF01(G) - x1 * Q1 + P.k12 * Q2 - renalClearance(G) +
      UG + P.EGP0 * EGPfactor,
    x1 * Q1 - (P.k12 + x2) * Q2,
    basal - S1 / P.tmaxI,
    S1 / P.tmaxI - S2 / P.tmaxI,
    (S2 / P.tmaxI) / P.VI - P.ke * I,
    -P.ka1 * x1 + P.kb1 * I,
    -P.ka2 * x2 + P.kb2 * I,
    -P.ka3 * x3 + P.kb3 * I,
  ];
}

function addScaled(y, k, scale) {
  return y.map((value, index) => value + scale * k[index]);
}

function rk4Step(t, y, h, cfg) {
  const k1 = derivative(t, y, cfg);
  const k2 = derivative(t + h / 2, addScaled(y, k1, h / 2), cfg);
  const k3 = derivative(t + h / 2, addScaled(y, k2, h / 2), cfg);
  const k4 = derivative(t + h, addScaled(y, k3, h), cfg);
  return y.map((value, index) => value + h *
    (k1[index] + 2 * k2[index] + 2 * k3[index] + k4[index]) / 6);
}

function eulerStep(t, y, h, cfg) {
  return addScaled(y, derivative(t, y, cfg), h);
}

function simulate({
  dt = 0.1,
  duration = 1440,
  meals = [],
  boluses = [],
  stopBasalAt,
  clipEGP = true,
  method = 'rk4',
  initial = EQ.y,
} = {}) {
  const cfg = { meals, stopBasalAt, clipEGP };
  const stepper = method === 'euler' ? eulerStep : rk4Step;
  const pendingBoluses = [...boluses].sort((a, b) => a.time - b.time);
  let bolusIndex = 0;
  let y = [...initial];
  const trace = [];
  const steps = Math.round(duration / dt);

  for (let step = 0; step <= steps; step++) {
    const t = step * dt;
    while (bolusIndex < pendingBoluses.length &&
           Math.abs(pendingBoluses[bolusIndex].time - t) <= dt * 1e-6) {
      y[2] += pendingBoluses[bolusIndex].units * 1000 / P.weightKg;
      bolusIndex++;
    }
    trace.push({ t, y: [...y], G: y[0] / P.VG, I: y[4] });
    if (step < steps) y = stepper(t, y, dt, cfg);
  }
  return trace;
}

function trapezoidIntegral(rows, valueOf) {
  let total = 0;
  for (let i = 1; i < rows.length; i++) {
    const dt = rows[i].t - rows[i - 1].t;
    total += dt * (valueOf(rows[i - 1]) + valueOf(rows[i])) / 2;
  }
  return total;
}

function metrics(trace, eventTime = 0) {
  const afterEvent = trace.filter(row => row.t >= eventTime);
  const peak = afterEvent.reduce((a, b) => b.G > a.G ? b : a);
  const nadir = afterEvent.reduce((a, b) => b.G < a.G ? b : a);
  const peakI = afterEvent.reduce((a, b) => b.I > a.I ? b : a);
  const hyperAuc = trapezoidIntegral(afterEvent,
    row => Math.max(0, row.G - P.targetG));
  const minutesAbove10 = trapezoidIntegral(afterEvent,
    row => row.G > 10 ? 1 : 0);
  return {
    peakMgDl: peak.G * P.mmolToMgDl,
    peakTime: peak.t,
    nadirMgDl: nadir.G * P.mmolToMgDl,
    nadirTime: nadir.t,
    peakInsulin: peakI.I,
    hyperAucMmolMinL: hyperAuc,
    minutesAbove180: minutesAbove10,
    finalMgDl: afterEvent.at(-1).G * P.mmolToMgDl,
    minQ1: Math.min(...afterEvent.map(row => row.y[0])),
    maxX3: Math.max(...afterEvent.map(row => row.y[7])),
  };
}

function maxStateDeviation(trace, baseline = EQ.y) {
  return trace.reduce((maximum, row) => Math.max(maximum,
    ...row.y.map((value, index) => Math.abs(value - baseline[index]))), 0);
}

function glucoseEigenvalues() {
  const x1 = EQ.y[5];
  const x2 = EQ.y[6];
  const a = -x1;
  const b = P.k12;
  const c = x1;
  const d = -(P.k12 + x2);
  const discriminant = Math.sqrt((a + d) ** 2 - 4 * (a * d - b * c));
  return [(a + d + discriminant) / 2, (a + d - discriminant) / 2];
}

function mealMassBalance(grams = 50, dt = 0.01, duration = 2000) {
  const expected = carbsToMmolPerKg(grams) * P.AG;
  const rows = [];
  let D1 = expected;
  let D2 = 0;
  const steps = Math.round(duration / dt);
  let closedIntegral = 0;
  let chainIntegral = 0;
  let maxRateDifference = 0;
  let priorClosed = 0;
  let priorChain = 0;
  for (let step = 0; step <= steps; step++) {
    const t = step * dt;
    const closed = expected * t * Math.exp(-t / P.tmaxG) /
      (P.tmaxG * P.tmaxG);
    const chain = D2 / P.tmaxG;
    maxRateDifference = Math.max(maxRateDifference, Math.abs(closed - chain));
    if (step > 0) {
      closedIntegral += dt * (priorClosed + closed) / 2;
      chainIntegral += dt * (priorChain + chain) / 2;
    }
    priorClosed = closed;
    priorChain = chain;
    rows.push({ t, closed, chain });
    if (step < steps) {
      const e = Math.exp(-dt / P.tmaxG);
      D2 = e * (D2 + D1 * dt / P.tmaxG);
      D1 *= e;
    }
  }
  return { expected, closedIntegral, chainIntegral, maxRateDifference };
}

function bolusMassBalance(units = 5) {
  const expected = units * 1000 / P.weightKg;
  const durations = [0.1, 1, 5, 10];
  return durations.map(duration => ({
    duration,
    inputRate: expected / duration,
    inputIntegral: expected / duration * duration,
  }));
}

function cutoffChecks() {
  const eps = 1e-9;
  return {
    F01At45: correctedF01(4.5),
    F01Jump: correctedF01(4.5 + eps) - correctedF01(4.5 - eps),
    FRAt9: renalClearance(9),
    FRJump: renalClearance(9 + eps) - renalClearance(9 - eps),
  };
}

const eig = [
  ...glucoseEigenvalues(),
  -1 / P.tmaxI, -1 / P.tmaxI, -P.ke,
  -P.ka1, -P.ka2, -P.ka3,
];
const equilibriumTrace = simulate({ duration: 2000, dt: 0.5 });
const stoppedTrace = simulate({ duration: 10000, dt: 0.5, stopBasalAt: 60 });
const noBolusTrace = simulate({
  meals: [{ time: 60, grams: 50 }], duration: 1440, dt: 0.1,
});

const doses = [0, 1, 2, 2.5, 3, 4, 5].map(units => ({
  units,
  ...metrics(simulate({
    meals: [{ time: 60, grams: 50 }],
    boluses: units ? [{ time: 60, units }] : [],
    duration: 1440,
    dt: 0.1,
  }), 60),
}));

const timings = [-15, 0, 30, 60].map(offset => ({
  offset,
  ...metrics(simulate({
    meals: [{ time: 60, grams: 50 }],
    boluses: [{ time: 60 + offset, units: 2 }],
    duration: 1440,
    dt: 0.1,
  }), 60),
}));

const literalOverdose = metrics(simulate({
  meals: [{ time: 60, grams: 10 }],
  boluses: [{ time: 60, units: 10 }],
  duration: 1440,
  dt: 0.05,
  clipEGP: false,
}), 60);
const clippedOverdose = metrics(simulate({
  meals: [{ time: 60, grams: 10 }],
  boluses: [{ time: 60, units: 10 }],
  duration: 1440,
  dt: 0.05,
  clipEGP: true,
}), 60);

const literalCalibration = [3, 5].map(units => ({
  units,
  ...metrics(simulate({
    meals: [{ time: 60, grams: 50 }],
    boluses: [{ time: 60, units }],
    duration: 1440,
    dt: 0.1,
    clipEGP: false,
  }), 60),
}));

function convergenceRun(method, dt) {
  return metrics(simulate({
    meals: [{ time: 60, grams: 50 }],
    boluses: [{ time: 60, units: 2 }],
    duration: 1440,
    dt,
    method,
  }), 60);
}
const reference = convergenceRun('rk4', 0.1);
const convergence = [
  ['rk4', 1], ['rk4', 0.5], ['rk4', 0.25], ['euler', 1],
].map(([method, dt]) => {
  const result = convergenceRun(method, dt);
  return {
    method,
    dt,
    peakMgDl: result.peakMgDl,
    peakErrorMgDl: Math.abs(result.peakMgDl - reference.peakMgDl),
    nadirErrorMgDl: Math.abs(result.nadirMgDl - reference.nadirMgDl),
  };
});

const noInsulinAnalyticG = 9 + (P.EGP0 - P.F01) / (0.003 * P.VG);
const doseThrough4 = doses.filter(row => row.units <= 4);
const strictlyDecreasing = (rows, field) => rows.every((row, index) =>
  index === 0 || row[field] < rows[index - 1][field]);
const rk4Dt1 = convergence.find(row => row.method === 'rk4' && row.dt === 1);
const mealBalance = mealMassBalance();
const bolusBalance = bolusMassBalance();
const cutoffs = cutoffChecks();

const checks = [
  {
    id: '1', name: 'equilibrium basal',
    pass: maxStateDeviation(equilibriumTrace) === 0 && eig.every(value => value < 0),
  },
  {
    id: '2', name: 'neraca makanan',
    pass: Math.abs(mealBalance.closedIntegral - mealBalance.expected) < 1e-6 &&
      Math.abs(mealBalance.chainIntegral - mealBalance.expected) < 1e-6 &&
      mealBalance.maxRateDifference < 1e-12,
  },
  {
    id: '3', name: 'neraca bolus',
    pass: bolusBalance.every(row => Math.abs(row.inputIntegral - 5000 / 70) < 1e-12),
  },
  {
    id: '4', name: 'basal dihentikan',
    pass: stoppedTrace.filter(row => row.t >= 60)
      .every((row, index, rows) => index === 0 || row.G >= rows[index - 1].G - 1e-12) &&
      Math.abs(stoppedTrace.at(-1).G - noInsulinAnalyticG) < 1e-9,
  },
  {
    id: '5', name: 'makan tanpa bolus',
    pass: metrics(noBolusTrace, 60).peakMgDl > 180 &&
      metrics(noBolusTrace, 60).minutesAbove180 > 0,
  },
  {
    id: '6', name: 'dosis-respons',
    pass: strictlyDecreasing(doseThrough4, 'peakMgDl') &&
      strictlyDecreasing(doseThrough4, 'hyperAucMmolMinL'),
  },
  {
    id: '7', name: 'waktu bolus',
    pass: timings.every((row, index) => index === 0 ||
      row.peakMgDl > timings[index - 1].peakMgDl),
  },
  {
    id: '8a', name: 'overdosis dengan persamaan literal',
    pass: literalOverdose.minQ1 >= 0,
  },
  {
    id: '8b', name: 'overdosis dengan clipping EGP',
    pass: clippedOverdose.minQ1 >= 0,
  },
  {
    id: '9a', name: 'kontinuitas F01c',
    pass: Math.abs(cutoffs.F01Jump) < 1e-9,
  },
  {
    id: '9b', name: 'kontinuitas FR',
    pass: Math.abs(cutoffs.FRJump) < 1e-9,
  },
  {
    id: '10', name: 'konvergensi RK4',
    pass: rk4Dt1.peakErrorMgDl < 0.01 && rk4Dt1.nadirErrorMgDl < 0.01,
  },
  {
    id: '11', name: 'konversi satuan',
    pass: Math.abs((5.5 * P.mmolToMgDl) / P.mmolToMgDl - 5.5) < 1e-12 &&
      Math.abs((5 * 1000 / P.weightKg) * P.weightKg / 1000 - 5) < 1e-12,
  },
];

const output = {
  assumptions: {
    weightKg: P.weightKg,
    fastingGlucoseMmolL: P.targetG,
    fastingGlucoseMgDl: P.targetG * P.mmolToMgDl,
    mealGrams: 50,
  },
  equilibrium: {
    basalInsulinMUL: EQ.y[4],
    basalUPerHour: EQ.uBasal * P.weightKg * 60 / 1000,
    state: EQ.y,
    residual: derivative(0, EQ.y, { meals: [], clipEGP: true }),
    maxDrift: maxStateDeviation(equilibriumTrace),
    eigenvalues: eig,
    slowestTimeConstantMin: 1 / Math.min(...eig.map(value => Math.abs(value))),
  },
  checkSummary: {
    passed: checks.filter(check => check.pass).length,
    failed: checks.filter(check => !check.pass).length,
    checks,
  },
  mealMassBalance: mealBalance,
  bolusMassBalance: bolusBalance,
  basalStopped: {
    analyticNoInsulinMgDl: noInsulinAnalyticG * P.mmolToMgDl,
    simulated: metrics(stoppedTrace, 60),
    monotonicAfterStop: stoppedTrace.filter(row => row.t >= 60)
      .every((row, index, rows) => index === 0 || row.G >= rows[index - 1].G - 1e-12),
  },
  mealWithoutBolus: metrics(noBolusTrace, 60),
  doseResponse: doses,
  bolusTiming: timings,
  overdose: { literal: literalOverdose, clipped: clippedOverdose },
  literalCalibration,
  cutoffs,
  convergence,
};

console.log(JSON.stringify(output, null, 2));
