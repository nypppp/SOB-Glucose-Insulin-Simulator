'use strict';

/*
 * Verification of a deliberately simple educational T2D profile obtained by
 * changing parameters of the existing Normal G-X-I model. This is useful for
 * direct Normal-versus-T2D comparison; it is not a clinical disease model.
 */

const assert = require('node:assert/strict');

const NORMAL = Object.freeze({
  id: 'normal',
  p1: 0.03082,
  p2: 0.02093,
  p3: 1.062e-5,
  n: 0.3,
  p6: 0.003349,
  p5: 94,
  Gb: 92,
  Ib: 7.3,
});

const T2D_SIMPLE = Object.freeze({
  id: 't2d-simple',
  p1: 0.024,
  p2: 0.02093,
  p3: 2.655e-6,       // 25% of Normal: project severity setting
  n: 0.3,
  p6: 0.001675,       // 50% of Normal: project severity setting
  p5: 119,
  Gb: 117,
  Ib: 12,
});

const MEALS = Object.freeze({ small: 80, medium: 150, large: 220 });
const MEAL_TAU = 40;

function mealAppearance(t, meal) {
  const age = t - meal.time;
  if (age < 0) return 0;
  return meal.amplitude * age * Math.exp(-age / MEAL_TAU) / (MEAL_TAU ** 2);
}

function derivatives(y, t, tauLocal, p, meal) {
  const [G, X, I] = y;
  const D = meal ? mealAppearance(t, meal) : 0;
  return [
    -p.p1 * (G - p.Gb) - X * G + D,
    -p.p2 * X + p.p3 * (I - p.Ib),
    p.p6 * Math.max(G - p.p5, 0) * tauLocal - p.n * (I - p.Ib),
  ];
}

function addScaled(y, k, scale) {
  return y.map((value, index) => value + scale * k[index]);
}

function rk4(y, t, dt, tauLocal, p, meal) {
  const k1 = derivatives(y, t, tauLocal, p, meal);
  const k2 = derivatives(addScaled(y, k1, dt / 2), t + dt / 2, tauLocal, p, meal);
  const k3 = derivatives(addScaled(y, k2, dt / 2), t + dt / 2, tauLocal, p, meal);
  const k4 = derivatives(addScaled(y, k3, dt), t + dt, tauLocal, p, meal);
  return y.map((value, index) => value + dt *
    (k1[index] + 2 * k2[index] + 2 * k3[index] + k4[index]) / 6);
}

function simulate(p, { duration = 1000, dt = 0.1, meal = null,
  initial = [p.Gb, 0, p.Ib] } = {}) {
  let y = [...initial];
  let tauLocal = 0;
  const rows = [{ t: 0, G: y[0], X: y[1], I: y[2] }];
  for (let t = 0; t < duration - dt / 2; t += dt) {
    tauLocal = y[0] <= p.p5 ? 0 : tauLocal + dt;
    y = rk4(y, t, dt, tauLocal, p, meal);
    assert.ok(y.every(Number.isFinite), `non-finite state at ${t + dt} min`);
    assert.ok(y.every(value => value >= 0), `negative state at ${t + dt} min`);
    rows.push({ t: t + dt, G: y[0], X: y[1], I: y[2] });
  }
  return rows;
}

function summary(rows, p) {
  const peakG = rows.reduce((best, row) => row.G > best.G ? row : best);
  const peakI = rows.reduce((best, row) => row.I > best.I ? row : best);
  const final = rows.at(-1);
  return {
    peakG: peakG.G,
    peakGTime: peakG.t,
    peakI: peakI.I,
    peakITime: peakI.t,
    finalG: final.G,
    finalI: final.I,
    finalX: final.X,
    maxGExcursion: peakG.G - p.Gb,
  };
}

function maxBasalDeviation(rows, p) {
  return Math.max(...rows.map(row => Math.max(
    Math.abs(row.G - p.Gb), Math.abs(row.X), Math.abs(row.I - p.Ib))));
}

function runProfile(p, dt = 0.1) {
  return Object.fromEntries(Object.entries(MEALS).map(([size, amplitude]) => [
    size,
    summary(simulate(p, { dt, meal: { time: 60, amplitude } }), p),
  ]));
}

const basalNormal = simulate(NORMAL, { duration: 1000 });
const basalT2d = simulate(T2D_SIMPLE, { duration: 1000 });
const normal = runProfile(NORMAL);
const t2d = runProfile(T2D_SIMPLE);

const perturbations = [
  [T2D_SIMPLE.Gb + 20, 0, T2D_SIMPLE.Ib],
  [T2D_SIMPLE.Gb, 0.01, T2D_SIMPLE.Ib],
  [T2D_SIMPLE.Gb, 0, T2D_SIMPLE.Ib + 10],
].map(initial => summary(simulate(T2D_SIMPLE, { initial, duration: 1500 }), T2D_SIMPLE));

const sensitivityAudit = {
  t2dSimple: t2d.medium,
  normalInsulinSensitivity: summary(simulate(
    { ...T2D_SIMPLE, p3: NORMAL.p3 },
    { meal: { time: 60, amplitude: MEALS.medium } },
  ), { ...T2D_SIMPLE, p3: NORMAL.p3 }),
  normalBetaResponse: summary(simulate(
    { ...T2D_SIMPLE, p6: NORMAL.p6 },
    { meal: { time: 60, amplitude: MEALS.medium } },
  ), { ...T2D_SIMPLE, p6: NORMAL.p6 }),
};

assert.ok(maxBasalDeviation(basalNormal, NORMAL) < 1e-12);
assert.ok(maxBasalDeviation(basalT2d, T2D_SIMPLE) < 1e-12);
assert.ok(Object.values(t2d).every(row => row.finalG > T2D_SIMPLE.Gb - 0.01 &&
  row.finalG < T2D_SIMPLE.Gb + 0.01), 'T2D profile does not return to baseline');
assert.ok(t2d.small.peakG < t2d.medium.peakG && t2d.medium.peakG < t2d.large.peakG,
  'meal dose response is not monotonic');
assert.ok(t2d.medium.maxGExcursion > normal.medium.maxGExcursion,
  'T2D excursion is not larger than Normal under the same input');
assert.ok(perturbations.every(row => Math.abs(row.finalG - T2D_SIMPLE.Gb) < 0.01 &&
  Math.abs(row.finalI - T2D_SIMPLE.Ib) < 0.01 && Math.abs(row.finalX) < 1e-6),
  'perturbed T2D profile does not recover to equilibrium');
assert.ok(sensitivityAudit.normalInsulinSensitivity.peakG < t2d.medium.peakG,
  'restoring Normal insulin sensitivity does not reduce the glucose peak');
assert.ok(sensitivityAudit.normalBetaResponse.peakG < t2d.medium.peakG,
  'restoring Normal beta-cell response does not reduce the glucose peak');

const shifted = summary(simulate(T2D_SIMPLE, {
  duration: 1400,
  meal: { time: 400, amplitude: MEALS.medium },
}), T2D_SIMPLE);
assert.ok(Math.abs(shifted.maxGExcursion - t2d.medium.maxGExcursion) < 1e-7,
  'response depends on absolute meal time');

const convergence = [1, 0.5, 0.25, 0.1].map(dt => ({
  dt,
  ...summary(simulate(T2D_SIMPLE, {
    dt,
    meal: { time: 60, amplitude: MEALS.medium },
  }), T2D_SIMPLE),
}));
assert.ok(convergence.slice(0, -1).every(row =>
  Math.abs(row.peakG - convergence.at(-1).peakG) < 0.1));

console.log(JSON.stringify({
  status: 'PASS',
  model: 'T2D simple profile derived from Normal G-X-I equations',
  parameters: { normal: NORMAL, t2dSimple: T2D_SIMPLE },
  basalMaxDeviation: {
    normal: maxBasalDeviation(basalNormal, NORMAL),
    t2dSimple: maxBasalDeviation(basalT2d, T2D_SIMPLE),
  },
  sameInputComparison: { normal, t2dSimple: t2d },
  perturbationRecovery: perturbations,
  sensitivityAudit,
  convergence,
  limitations: [
    'p1, p3, and p6 changes are project severity settings, not a fitted T2D cohort',
    'the model has no explicit glucagon or incretin state',
    'educational comparison only; not for diagnosis, prediction, or treatment',
  ],
}, null, 2));
