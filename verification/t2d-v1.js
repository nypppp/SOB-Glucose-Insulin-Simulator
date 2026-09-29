'use strict';

/*
 * Numerical verification for the educational T2D v1 plant.
 *
 * The three-state G-I-A structure follows Subramanian et al. (2024), while
 * centering each regulatory term on its fasting value is a project adaptation
 * that gives an exact reusable equilibrium for continuous simulation. It is not
 * a clinical predictor and must not be used to calculate treatment.
 */

const assert = require('node:assert/strict');

const P = {
  weightKg: 70,
  Gb: 117,                 // mg/dL; educational virtual-subject baseline
  Ib: 12,                  // µU/mL; early-T2D compensatory fasting insulin
  Ab: 20,                  // pM
  SG: 0.014,               // 1/min; Subramanian fixed population value
  a1: 6.6e-5,              // (µU/mL min)^-1; converted T2D mean
  a2: 0.16,                // effective centered coefficient; project calibration
  n1: 0.14,                // 1/min
  n2: 0.08,                // 1/min
  betaG: 3.2,              // µU/mL/min; project calibration
  betaE: 0.55,             // µU/mL/min; project incretin surrogate
  betaA: 5.0,              // pM/min; T2D mean glucagon magnitude
  K: 306,                  // mg/dL (17 mmol/L)
  hill: 1.27,
  kSuppression: 0.15,      // L/mmol; T2D mean
  AG: 0.8,
  VG: 0.16,                // L/kg
  tmaxG: 40,               // min
  incretinTau: 35,         // min; project surrogate
  mmolToMgDl: 18.0182,
};

function hill(G) {
  const ratio = Math.max(G, 0) / P.K;
  const power = ratio ** P.hill;
  return 1.5 * power / (1 + power);
}

function glucagonResponse(G) {
  return Math.exp(-P.kSuppression * Math.max(G, 0) / P.mmolToMgDl);
}

function carbsToMmolKg(grams) {
  return grams * 1000 / (180.16 * P.weightKg);
}

function mealAppearanceMgDlMin(t, meals) {
  let total = 0;
  for (const meal of meals) {
    const age = t - meal.time;
    if (age >= 0) {
      const mmolKgMin = carbsToMmolKg(meal.grams) * P.AG * age *
        Math.exp(-age / P.tmaxG) / (P.tmaxG * P.tmaxG);
      total += mmolKgMin / P.VG * P.mmolToMgDl;
    }
  }
  return total;
}

/* Dimensionless meal-driven incretin surrogate, normalized to peak at 1. */
function incretinSignal(t, meals) {
  let total = 0;
  for (const meal of meals) {
    const age = t - meal.time;
    if (age >= 0) {
      const sizeScale = meal.grams / 50;
      total += sizeScale * (age / P.incretinTau) *
        Math.exp(1 - age / P.incretinTau);
    }
  }
  return total;
}

function derivative(t, y, meals) {
  const [G, I, A] = y;
  const dG = -P.SG * (G - P.Gb) - P.a1 * (I - P.Ib) * Math.max(G, 0) +
    P.a2 * (A - P.Ab) + mealAppearanceMgDlMin(t, meals);
  const dI = -P.n1 * (I - P.Ib) + P.betaG * (hill(G) - hill(P.Gb)) +
    P.betaE * incretinSignal(t, meals) * hill(G);
  const dA = -P.n2 * (A - P.Ab) +
    P.betaA * (glucagonResponse(G) - glucagonResponse(P.Gb));
  return [dG, dI, dA];
}

function addScaled(y, k, scale) {
  return y.map((value, index) => value + scale * k[index]);
}

function rk4Step(t, y, h, meals) {
  const k1 = derivative(t, y, meals);
  const k2 = derivative(t + h / 2, addScaled(y, k1, h / 2), meals);
  const k3 = derivative(t + h / 2, addScaled(y, k2, h / 2), meals);
  const k4 = derivative(t + h, addScaled(y, k3, h), meals);
  return y.map((value, index) => value + h *
    (k1[index] + 2 * k2[index] + 2 * k3[index] + k4[index]) / 6);
}

function simulate({ duration = 1000, dt = 0.1, meals = [], assertPositive = true } = {}) {
  let y = [P.Gb, P.Ib, P.Ab];
  const rows = [{ t: 0, G: y[0], I: y[1], A: y[2] }];
  const steps = Math.round(duration / dt);
  for (let step = 0; step < steps; step++) {
    const t = step * dt;
    y = rk4Step(t, y, dt, meals);
    assert.ok(y.every(Number.isFinite), `non-finite state at ${t + dt} min`);
    if (assertPositive) {
      assert.ok(y.every(value => value >= 0), `negative state at ${t + dt} min`);
    }
    rows.push({ t: t + dt, G: y[0], I: y[1], A: y[2] });
  }
  return rows;
}

function simulateFrom(initial, { duration = 1000, dt = 0.1, meals = [] } = {}) {
  let y = [...initial];
  const rows = [{ t: 0, G: y[0], I: y[1], A: y[2] }];
  const steps = Math.round(duration / dt);
  for (let step = 0; step < steps; step++) {
    const t = step * dt;
    y = rk4Step(t, y, dt, meals);
    assert.ok(y.every(Number.isFinite));
    assert.ok(y.every(value => value >= 0));
    rows.push({ t: t + dt, G: y[0], I: y[1], A: y[2] });
  }
  return rows;
}

function metrics(rows, mealTime = 60) {
  const post = rows.filter(row => row.t >= mealTime);
  const peakG = post.reduce((best, row) => row.G > best.G ? row : best, post[0]);
  const peakI = post.reduce((best, row) => row.I > best.I ? row : best, post[0]);
  const nadirA = post.reduce((best, row) => row.A < best.A ? row : best, post[0]);
  return {
    peakG: peakG.G, peakGTime: peakG.t,
    peakI: peakI.I, peakITime: peakI.t,
    nadirA: nadirA.A, nadirATime: nadirA.t,
    finalG: rows.at(-1).G,
    finalI: rows.at(-1).I,
    finalA: rows.at(-1).A,
    minG: Math.min(...post.map(row => row.G)),
    minI: Math.min(...post.map(row => row.I)),
    minA: Math.min(...post.map(row => row.A)),
  };
}

function maxDeviation(rows) {
  return Math.max(...rows.map(row => Math.max(
    Math.abs(row.G - P.Gb), Math.abs(row.I - P.Ib), Math.abs(row.A - P.Ab))));
}

function withParameter(key, value, run) {
  const original = P[key];
  P[key] = value;
  try { return run(); } finally { P[key] = original; }
}

function determinant3(m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
    m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
    m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

function equilibriumJacobian(step = 1e-5) {
  const base = [P.Gb, P.Ib, P.Ab];
  return base.map((_, column) => {
    const plus = [...base], minus = [...base];
    plus[column] += step;
    minus[column] -= step;
    const fp = derivative(0, plus, []), fm = derivative(0, minus, []);
    return fp.map((value, row) => (value - fm[row]) / (2 * step));
  }).reduce((matrix, column, columnIndex) => {
    column.forEach((value, rowIndex) => { matrix[rowIndex][columnIndex] = value; });
    return matrix;
  }, [[], [], []]);
}

function routhHurwitz3(jacobian) {
  const a = -(jacobian[0][0] + jacobian[1][1] + jacobian[2][2]);
  const b = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0] +
    jacobian[0][0] * jacobian[2][2] - jacobian[0][2] * jacobian[2][0] +
    jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1];
  const c = -determinant3(jacobian);
  return { a, b, c, aTimesB: a * b, stable: a > 0 && b > 0 && c > 0 && a * b > c };
}

const basal = simulate();
assert.ok(maxDeviation(basal) < 1e-10, 'fasting equilibrium drifts');
const jacobian = equilibriumJacobian();
const localStability = routhHurwitz3(jacobian);
assert.ok(localStability.stable, 'fasting equilibrium is not locally stable');

for (const initial of [[150, P.Ib, P.Ab], [P.Gb, 25, P.Ab], [P.Gb, P.Ib, 35]]) {
  const final = simulateFrom(initial).at(-1);
  assert.ok(Math.abs(final.G - P.Gb) < 0.01 &&
    Math.abs(final.I - P.Ib) < 0.01 && Math.abs(final.A - P.Ab) < 0.01,
  'perturbed state does not return to fasting equilibrium');
}

const scenarios = [25, 50, 75].map(grams => ({
  grams,
  ...metrics(simulate({ meals: [{ time: 60, grams }] })),
}));
assert.ok(scenarios.every((row, index, all) => index === 0 ||
  row.peakG > all[index - 1].peakG), 'meal glucose dose-response is not monotonic');
assert.ok(scenarios.every(row => row.peakI > P.Ib), 'insulin does not respond to meals');
assert.ok(scenarios.every(row => row.nadirA < P.Ab), 'glucagon is not suppressed');
assert.ok(scenarios.every(row => Math.abs(row.finalG - P.Gb) < 0.05),
  'glucose does not return near baseline');

const stress = simulate({ meals: [{ time: 60, grams: 125 }] });
assert.ok(stress.every(row => row.G >= 0 && row.I >= 0 && row.A >= 0),
  '125 g stress scenario violates positivity');

const mealBalanceExpected = carbsToMmolKg(50) * P.AG / P.VG * P.mmolToMgDl;
let mealBalanceIntegrated = 0;
const balanceDt = 0.01;
for (let t = 0; t < 1000; t += balanceDt) {
  mealBalanceIntegrated += mealAppearanceMgDlMin(t + balanceDt / 2,
    [{ time: 0, grams: 50 }]) * balanceDt;
}
assert.ok(Math.abs(mealBalanceIntegrated - mealBalanceExpected) < 1e-6,
  'meal input mass balance fails');

const originalA1 = P.a1;
const sensitivitySweep = [0.5, 1, 2].map(scale => {
  P.a1 = originalA1 * scale;
  return { scale, ...metrics(simulate({ meals: [{ time: 60, grams: 50 }] })) };
});
P.a1 = originalA1;
assert.ok(sensitivitySweep[0].peakG > sensitivitySweep[1].peakG &&
  sensitivitySweep[1].peakG > sensitivitySweep[2].peakG,
'greater insulin sensitivity does not lower glucose peak');

const pathwayChecks = {
  noInsulinAction: withParameter('a1', 0, () =>
    metrics(simulate({ meals: [{ time: 60, grams: 50 }] }))),
  paperControlA1: withParameter('a1', 3e-4, () =>
    metrics(simulate({ meals: [{ time: 60, grams: 50 }] }))),
  noGlucagonAction: withParameter('a2', 0, () =>
    metrics(simulate({ meals: [{ time: 60, grams: 50 }] }))),
  noIncretin: withParameter('betaE', 0, () =>
    metrics(simulate({ meals: [{ time: 60, grams: 50 }] }))),
};
assert.ok(pathwayChecks.noInsulinAction.peakG > scenarios[1].peakG);
assert.ok(pathwayChecks.paperControlA1.peakG < scenarios[1].peakG);
assert.ok(pathwayChecks.noIncretin.peakI < scenarios[1].peakI);
const physiologicalAttribution = {
  insulinActionPeakEffect: pathwayChecks.noInsulinAction.peakG - scenarios[1].peakG,
  controlVsT2dA1PeakEffect: scenarios[1].peakG - pathwayChecks.paperControlA1.peakG,
  glucagonActionPeakEffect: pathwayChecks.noGlucagonAction.peakG - scenarios[1].peakG,
  incretinInsulinPeakEffect: scenarios[1].peakI - pathwayChecks.noIncretin.peakI,
};
const physiologicalGate = {
  name: 'insulin-resistance pathway materially affects meal glucose',
  thresholdMgDl: 5,
  observedMgDl: physiologicalAttribution.controlVsT2dA1PeakEffect,
  pass: physiologicalAttribution.controlVsT2dA1PeakEffect >= 5,
};

/* Exploratory calibration: the paper mean a2 is not jointly identifiable with
 * centered mean states, so test an effective project coefficient. */
const a2CalibrationSweep = [0.05, 0.08, 0.10, 0.13, 0.16, 0.20, 0.26].map(a2 =>
  withParameter('a2', a2, () => {
    const t2d = metrics(simulate({ meals: [{ time: 60, grams: 50 }] }));
    const control = withParameter('a1', 3e-4, () =>
      metrics(simulate({ meals: [{ time: 60, grams: 50 }] })));
    const stressRows = simulate({ meals: [{ time: 60, grams: 125 }], assertPositive: false });
    return {
      a2,
      t2dPeakG: t2d.peakG,
      controlPeakG: control.peakG,
      insulinSensitivitySeparation: t2d.peakG - control.peakG,
      stressMinimum: Math.min(...stressRows.flatMap(row => [row.G, row.I, row.A])),
    };
  }));

const shifted = metrics(simulate({ duration: 1300, meals: [{ time: 300, grams: 50 }] }), 300);
assert.ok(Math.abs(shifted.peakG - scenarios[1].peakG) < 1e-8,
  'meal response depends on absolute simulation time');

const convergence = [1, 0.5, 0.25, 0.1].map(dt => ({
  dt,
  ...metrics(simulate({ dt, meals: [{ time: 60, grams: 50 }] })),
}));
const reference = convergence.at(-1);
assert.ok(convergence.slice(0, -1).every(row =>
  Math.abs(row.peakG - reference.peakG) < 0.02), 'RK4 peak is not converged');

console.log(JSON.stringify({
  status: physiologicalGate.pass ? 'PASS' : 'STRUCTURAL_PASS_FUNCTIONAL_FAIL',
  model: 'T2D educational v1, equilibrium-centered G-I-A',
  basalMaxDeviation: maxDeviation(basal),
  jacobian,
  localStability,
  mealMassBalance: { expected: mealBalanceExpected, integrated: mealBalanceIntegrated },
  scenarios,
  sensitivitySweep,
  pathwayChecks,
  physiologicalAttribution,
  physiologicalGate,
  a2CalibrationSweep,
  convergence,
  readyForUI: physiologicalGate.pass,
  limitations: [
    'betaG, betaE, and incretinTau are project calibration parameters',
    'incretin is a meal-driven surrogate, not measured GLP-1/GIP',
    'clinical validation and individual prediction are out of scope',
  ],
}, null, 2));

