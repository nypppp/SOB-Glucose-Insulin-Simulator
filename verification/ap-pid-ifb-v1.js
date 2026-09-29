'use strict';

/*
 * Numerical verification of the educational AP PID/PID-IFB specification.
 * Plant: eight-state Hovorka T1D v2 with non-negative EGP.
 * Controller: one-minute discrete PID and Ruiz et al. (2012) IFB estimator.
 */

const P = {
  weightKg: 70,
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

const IFB = {
  a11: 0.9802,
  a21: 0.014043,
  a31: 0.000127,
  a22: 0.98582,
  a32: 0.017889,
  a33: 0.98198,
  b1: 1.1881,
  b2: 0.0084741,
  b3: 0.00005,
  g1: 0.64935,
  g2: 0.34128,
  g3: 0.0093667,
};
IFB.gSum = IFB.g1 + IFB.g2 + IFB.g3;

const C = {
  targetMgDl: 120,
  controllerPeriod: 1,
  // Ruiz et al. starts at 150 min. The project uses 300 min after the
  // documented meal-response sweep to avoid late hypoglycemia on this plant.
  Ti: 300,
  Td: 75,
  derivativeTau: 5,
  sensorTau: 10,
  suspendMgDl: 70,
};

function f01c(G) {
  return G >= 4.5 ? P.F01 : P.F01 * G / 4.5;
}

function renal(G) {
  return G >= 9 ? 0.003 * (G - 9) * P.VG : 0;
}

function carbsToMmolKg(grams) {
  return grams * 1000 / (180.16 * P.weightKg);
}

function mealAppearance(t, meals) {
  let total = 0;
  for (const meal of meals) {
    const age = t - meal.time;
    if (age >= 0) {
      const available = carbsToMmolKg(meal.grams) * P.AG;
      total += available * age * Math.exp(-age / P.tmaxG) /
        (P.tmaxG * P.tmaxG);
    }
  }
  return total;
}

function bisect(fn, lo, hi, tolerance = 1e-13) {
  let flo = fn(lo);
  if (flo * fn(hi) > 0) throw new Error('Root not bracketed');
  while (hi - lo > tolerance) {
    const mid = (lo + hi) / 2;
    const fm = fn(mid);
    if (flo * fm <= 0) hi = mid;
    else {
      lo = mid;
      flo = fm;
    }
  }
  return (lo + hi) / 2;
}

function plantEquilibrium(targetMgDl = C.targetMgDl, sensitivityScale = 1) {
  const G = targetMgDl / P.mmolToMgDl;
  const residual = I => {
    const x1 = P.SIT * sensitivityScale * I;
    const x2 = P.SID * sensitivityScale * I;
    const x3 = P.SIE * sensitivityScale * I;
    const Q1 = P.VG * G;
    const Q2 = x1 * Q1 / (P.k12 + x2);
    return -f01c(G) - x1 * Q1 + P.k12 * Q2 - renal(G) +
      P.EGP0 * Math.max(0, 1 - x3);
  };
  const I = bisect(residual, 0, 1 / (P.SIE * sensitivityScale));
  const x1 = P.SIT * sensitivityScale * I;
  const x2 = P.SID * sensitivityScale * I;
  const x3 = P.SIE * sensitivityScale * I;
  const Q1 = P.VG * G;
  const Q2 = x1 * Q1 / (P.k12 + x2);
  const basalPlant = I * P.VI * P.ke;
  const S1 = basalPlant * P.tmaxI;
  return {
    y: [Q1, Q2, S1, S1, I, x1, x2, x3, targetMgDl],
    basalPlant,
    basalUph: basalPlant * P.weightKg * 60 / 1000,
  };
}

function ifbEquilibrium(rateUph) {
  const dosePerMinute = rateUph / 60;
  const Isc = IFB.b1 * dosePerMinute / (1 - IFB.a11);
  const Ip = (IFB.a21 * Isc + IFB.b2 * dosePerMinute) / (1 - IFB.a22);
  const Ieff = (IFB.a31 * Isc + IFB.a32 * Ip + IFB.b3 * dosePerMinute) /
    (1 - IFB.a33);
  return [Isc, Ip, Ieff];
}

function ifbStep(state, deliveredRateUph) {
  const [Isc, Ip, Ieff] = state;
  const dose = deliveredRateUph / 60;
  return [
    IFB.a11 * Isc + IFB.b1 * dose,
    IFB.a21 * Isc + IFB.a22 * Ip + IFB.b2 * dose,
    IFB.a31 * Isc + IFB.a32 * Ip + IFB.a33 * Ieff + IFB.b3 * dose,
  ];
}

function ifbValue(state) {
  return IFB.g1 * state[0] + IFB.g2 * state[1] + IFB.g3 * state[2];
}

function newController(mode, eq, overrides = {}) {
  const settings = { ...C, ...overrides };
  const dailyInsulinU = eq.basalUph * 24;
  const Kp = dailyInsulinU / 2250 * (overrides.KpScale ?? 1);
  const estimator = ifbEquilibrium(eq.basalUph);
  const feedbackAtEquilibrium = ifbValue(estimator);
  return {
    mode,
    settings,
    Kp,
    integral: mode === 'pid'
      ? eq.basalUph / (1 + IFB.gSum)
      : (eq.basalUph + feedbackAtEquilibrium) / (1 + IFB.gSum),
    derivative: 0,
    previousSG: settings.targetMgDl,
    estimator,
    deliveredRate: eq.basalUph,
    maxRate: overrides.maxRate ?? Math.max(3 * eq.basalUph, 3),
    minIntegral: 0,
    maxIntegral: 3 * eq.basalUph,
  };
}

function controllerOutput(controller, SG) {
  const c = controller;
  c.estimator = ifbStep(c.estimator, c.deliveredRate);
  const feedback = ifbValue(c.estimator);
  const error = SG - c.settings.targetMgDl;
  const rawDerivative = (SG - c.previousSG) / c.settings.controllerPeriod;
  const alpha = Math.exp(-c.settings.controllerPeriod / c.settings.derivativeTau);
  c.derivative = alpha * c.derivative + (1 - alpha) * rawDerivative;
  const candidateIntegral = Math.max(c.minIntegral, Math.min(c.maxIntegral,
    c.integral + c.Kp / c.settings.Ti * error * c.settings.controllerPeriod));

  const calculate = integral => {
    const proportional = c.Kp * error;
    const derivative = c.Kp * c.settings.Td * c.derivative;
    const pidCore = proportional + integral + derivative;
    const raw = c.mode === 'pid'
      ? (1 + IFB.gSum) * pidCore
      : (1 + IFB.gSum) * pidCore - feedback;
    return { proportional, integral, derivative, pidCore, feedback, raw };
  };

  let terms = calculate(candidateIntegral);
  let rate = Math.max(0, Math.min(c.maxRate, terms.raw));
  const saturatedHigh = terms.raw > c.maxRate && error > 0;
  const saturatedLow = terms.raw < 0 && error < 0;
  if (saturatedHigh || saturatedLow) {
    terms = calculate(c.integral);
    rate = Math.max(0, Math.min(c.maxRate, terms.raw));
  } else {
    c.integral = candidateIntegral;
  }

  let suspended = false;
  if (SG <= c.settings.suspendMgDl) {
    rate = 0;
    suspended = true;
  }
  c.previousSG = SG;
  c.deliveredRate = rate;
  return { ...terms, rate, error, suspended };
}

function derivative(t, y, cfg, rateUph) {
  const [Q1, Q2, S1, S2, I, x1, x2, x3, Gs] = y;
  const G = Q1 / P.VG;
  const GmgDl = G * P.mmolToMgDl;
  const u = rateUph * 1000 / (60 * P.weightKg);
  const sensitivity = cfg.sensitivityScale;
  return [
    -f01c(G) - x1 * Q1 + P.k12 * Q2 - renal(G) +
      mealAppearance(t, cfg.meals) + P.EGP0 * Math.max(0, 1 - x3),
    x1 * Q1 - (P.k12 + x2) * Q2,
    u - S1 / P.tmaxI,
    S1 / P.tmaxI - S2 / P.tmaxI,
    (S2 / P.tmaxI) / P.VI - P.ke * I,
    -P.ka1 * x1 + P.ka1 * P.SIT * sensitivity * I,
    -P.ka2 * x2 + P.ka2 * P.SID * sensitivity * I,
    -P.ka3 * x3 + P.ka3 * P.SIE * sensitivity * I,
    (GmgDl - Gs) / cfg.sensorTau,
  ];
}

function addScaled(y, k, scale) {
  return y.map((value, index) => value + scale * k[index]);
}

function rk4Step(t, y, h, cfg, rate) {
  const k1 = derivative(t, y, cfg, rate);
  const k2 = derivative(t + h / 2, addScaled(y, k1, h / 2), cfg, rate);
  const k3 = derivative(t + h / 2, addScaled(y, k2, h / 2), cfg, rate);
  const k4 = derivative(t + h, addScaled(y, k3, h), cfg, rate);
  return y.map((value, index) => value + h *
    (k1[index] + 2 * k2[index] + 2 * k3[index] + k4[index]) / 6);
}

function simulate({
  mode = 'pid-ifb', duration = 1440, plantDt = 0.1, meals = [],
  sensorTau = C.sensorTau, sensitivityScale = 1, controllerOverrides = {},
} = {}) {
  // Keep controller tuning and initial state nominal while perturbing plant
  // sensitivity, so robustness runs test mismatch instead of retuning the AP.
  const eq = plantEquilibrium(C.targetMgDl, 1);
  const controller = newController(mode, eq, { sensorTau, ...controllerOverrides });
  const cfg = { meals, sensorTau, sensitivityScale };
  let y = [...eq.y];
  let current = controllerOutput(controller, y[8]);
  const rows = [];
  const steps = Math.round(duration / plantDt);
  const controlEvery = Math.round(C.controllerPeriod / plantDt);

  for (let step = 0; step <= steps; step++) {
    const t = step * plantDt;
    if (step > 0 && step % controlEvery === 0) {
      current = controllerOutput(controller, y[8]);
    }
    rows.push({
      t, G: y[0] / P.VG * P.mmolToMgDl, SG: y[8], rate: current.rate,
      P: current.proportional, Ic: controller.integral, D: current.derivative,
      IFB: current.feedback, suspended: current.suspended, y: [...y],
    });
    if (step < steps) y = rk4Step(t, y, plantDt, cfg, current.rate);
  }
  return { rows, eq, controller };
}

function integrate(rows, valueOf) {
  let sum = 0;
  for (let i = 1; i < rows.length; i++) {
    const dt = rows[i].t - rows[i - 1].t;
    sum += dt * (valueOf(rows[i - 1]) + valueOf(rows[i])) / 2;
  }
  return sum;
}

function metrics(run, eventTime = 0) {
  const rows = run.rows.filter(row => row.t >= eventTime);
  const peak = rows.reduce((a, b) => b.G > a.G ? b : a);
  const nadir = rows.reduce((a, b) => b.G < a.G ? b : a);
  const totalInsulin = integrate(rows, row => row.rate / 60);
  const duration = rows.at(-1).t - rows[0].t;
  const timeIn = predicate => integrate(rows, row => predicate(row.G) ? 1 : 0);
  return {
    peakMgDl: peak.G,
    peakTime: peak.t,
    nadirMgDl: nadir.G,
    nadirTime: nadir.t,
    finalMgDl: rows.at(-1).G,
    totalInsulinU: totalInsulin,
    timeInRangePercent: 100 * timeIn(G => G >= 70 && G <= 180) / duration,
    minutesBelow70: timeIn(G => G < 70),
    minutesBelow54: timeIn(G => G < 54),
    minutesAbove180: timeIn(G => G > 180),
    maxRateUph: Math.max(...rows.map(row => row.rate)),
    maxIntegralUph: Math.max(...rows.map(row => row.Ic)),
    suspendedMinutes: integrate(rows, row => row.suspended ? 1 : 0),
  };
}

function maxEquilibriumDeviation(run) {
  return Math.max(...run.rows.flatMap(row => [
    Math.abs(row.G - C.targetMgDl),
    Math.abs(row.SG - C.targetMgDl),
    Math.abs(row.rate - run.eq.basalUph),
  ]));
}

const equilibriumPID = simulate({ mode: 'pid', duration: 720 });
const equilibriumIFB = simulate({ mode: 'pid-ifb', duration: 720 });
const mealSizes = [25, 50, 75];
const mealsPID = mealSizes.map(grams => ({
  grams,
  ...metrics(simulate({ mode: 'pid', meals: [{ time: 60, grams }] }), 60),
}));
const mealsIFB = mealSizes.map(grams => ({
  grams,
  ...metrics(simulate({ mode: 'pid-ifb', meals: [{ time: 60, grams }] }), 60),
}));
const sensorRuns = [5, 10, 15].map(sensorTau => ({
  sensorTau,
  ...metrics(simulate({
    mode: 'pid-ifb', sensorTau, meals: [{ time: 60, grams: 50 }],
  }), 60),
}));
const robustness = [0.8, 1, 1.2].map(sensitivityScale => ({
  sensitivityScale,
  ...metrics(simulate({
    mode: 'pid-ifb', sensitivityScale, meals: [{ time: 60, grams: 50 }],
  }), 60),
}));
const convergence = [0.5, 0.25, 0.1, 0.05].map(plantDt => ({
  plantDt,
  ...metrics(simulate({
    mode: 'pid-ifb', plantDt, meals: [{ time: 60, grams: 50 }],
  }), 60),
}));

const eq = plantEquilibrium();
const estimatorEq = ifbEquilibrium(eq.basalUph);
const estimatorAfter = ifbStep(estimatorEq, eq.basalUph);
const ifbEqError = Math.max(...estimatorEq.map((v, i) =>
  Math.abs(v - estimatorAfter[i])));
const ifbAtEq = ifbValue(estimatorEq);
const expectedIfbAtEq = IFB.gSum * eq.basalUph;
const conversionRoundTrip = eq.basalUph * 1000 / (60 * P.weightKg) *
  60 * P.weightKg / 1000;
const ref = convergence.at(-1);
const dt01 = convergence.find(row => row.plantDt === 0.1);
const pid50 = mealsPID.find(row => row.grams === 50);
const ifb50 = mealsIFB.find(row => row.grams === 50);

const synthetic = newController('pid-ifb', eq);
synthetic.previousSG = 120;
const syntheticTerms = controllerOutput(synthetic, 130);
const saturated = newController('pid-ifb', eq, { maxRate: eq.basalUph * 1.01 });
for (let i = 0; i < 500; i++) controllerOutput(saturated, 300);
const suspended = newController('pid-ifb', eq);
const suspendTerms = controllerOutput(suspended, 65);

const checks = [
  { id: 1, name: 'equilibrium', pass: maxEquilibriumDeviation(equilibriumIFB) < 1e-9 },
  { id: 2, name: 'bumpless start PID dan PID-IFB', pass:
    Math.abs(equilibriumPID.rows[0].rate - eq.basalUph) < 1e-9 &&
    Math.abs(equilibriumIFB.rows[0].rate - eq.basalUph) < 1e-9 },
  { id: 3, name: 'konversi satuan', pass:
    Math.abs(conversionRoundTrip - eq.basalUph) < 1e-12 },
  { id: 4, name: 'komponen PID sintetis', pass:
    syntheticTerms.proportional > 0 && syntheticTerms.derivative > 0 &&
    synthetic.integral > eq.basalUph },
  { id: 5, name: 'equilibrium estimator IFB', pass:
    ifbEqError < 1e-12 },
  { id: 6, name: 'dosis-respons makanan', pass:
    mealsIFB.every((row, i) => i === 0 || row.peakMgDl > mealsIFB[i - 1].peakMgDl) },
  { id: 7, name: 'IFB membatasi insulin aktif', pass:
    ifb50.totalInsulinU < pid50.totalInsulinU },
  { id: 8, name: 'anti-windup', pass:
    saturated.integral <= saturated.maxIntegral + 1e-12 },
  { id: 9, name: 'low-glucose suspend', pass:
    suspendTerms.rate === 0 && suspendTerms.suspended },
  { id: 10, name: 'delay sensor', pass:
    sensorRuns.every(row => Number.isFinite(row.peakMgDl) && row.nadirMgDl >= 0) },
  { id: 11, name: 'robustness sensitivitas', pass:
    robustness.every(row => row.nadirMgDl >= 0 && Number.isFinite(row.peakMgDl)) },
  { id: 12, name: 'konvergensi plant', pass:
    Math.abs(dt01.peakMgDl - ref.peakMgDl) < 0.01 &&
    Math.abs(dt01.nadirMgDl - ref.nadirMgDl) < 0.01 },
  { id: 13, name: 'metrik lengkap', pass:
    [pid50, ifb50].every(row => Object.values(row).every(Number.isFinite)) },
];

const performanceGates = [
  {
    name: 'tanpa hipoglikemia pada makan 25/50/75 g',
    pass: mealsIFB.every(row => row.minutesBelow70 === 0),
  },
  {
    name: 'tanpa hipoglikemia berat pada makan 25/50/75 g',
    pass: mealsIFB.every(row => row.minutesBelow54 === 0),
  },
  {
    name: 'time-in-range minimal 70 persen',
    pass: mealsIFB.every(row => row.timeInRangePercent >= 70),
  },
  {
    name: 'kembali dekat target pada akhir simulasi',
    pass: mealsIFB.every(row => Math.abs(row.finalMgDl - C.targetMgDl) <= 10),
  },
];

const tuningCandidates = [];
for (const KpScale of [0.25, 0.5, 0.75, 1]) {
  for (const Ti of [150, 300, 450]) {
    for (const Td of [40, 75]) {
      const scenarios = mealSizes.map(grams => metrics(simulate({
        mode: 'pid-ifb',
        meals: [{ time: 60, grams }],
        controllerOverrides: { KpScale, Ti, Td },
      }), 60));
      const score = scenarios.reduce((sum, row) => sum +
        row.minutesAbove180 + 100 * row.minutesBelow70 +
        1000 * row.minutesBelow54 + Math.abs(row.finalMgDl - C.targetMgDl), 0);
      tuningCandidates.push({ KpScale, Ti, Td, score, scenarios });
    }
  }
}
tuningCandidates.sort((a, b) => a.score - b.score);

console.log(JSON.stringify({
  assumptions: {
    targetMgDl: C.targetMgDl,
    weightKg: P.weightKg,
    controllerPeriodMin: C.controllerPeriod,
    KpMapping: 'Kp = (24 * basal U/h) / 2250',
    Kp: newController('pid-ifb', eq).Kp,
    Ti: C.Ti,
    Td: C.Td,
    sensorTau: C.sensorTau,
    suspendMgDl: C.suspendMgDl,
  },
  equilibrium: {
    basalUph: eq.basalUph,
    basalPlant_mU_kg_min: eq.basalPlant,
    maxDeviationPID: maxEquilibriumDeviation(equilibriumPID),
    maxDeviationIFB: maxEquilibriumDeviation(equilibriumIFB),
    estimatorState: estimatorEq,
    estimatorFixedPointError: ifbEqError,
    ifbAtEquilibrium: ifbAtEq,
  },
  checkSummary: {
    passed: checks.filter(check => check.pass).length,
    failed: checks.filter(check => !check.pass).length,
    checks,
  },
  performanceGates: {
    passed: performanceGates.filter(gate => gate.pass).length,
    failed: performanceGates.filter(gate => !gate.pass).length,
    gates: performanceGates,
    readyForUI: checks.every(check => check.pass) &&
      performanceGates.every(gate => gate.pass),
  },
  mealPID: mealsPID,
  mealPIDIFB: mealsIFB,
  comparison50g: { pid: pid50, pidIFB: ifb50 },
  sensorDelay: sensorRuns,
  sensitivityRobustness: robustness,
  convergence,
  safety: {
    antiWindupIntegral: saturated.integral,
    antiWindupLimit: saturated.maxIntegral,
    suspendOutput: suspendTerms.rate,
  },
  exploratoryTuningTop5: tuningCandidates.slice(0, 5),
}, null, 2));
