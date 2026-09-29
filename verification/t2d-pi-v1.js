'use strict';

/*
 * Educational PI control of the T2D G-X-I plant.
 * Controller output is an external subcutaneous insulin infusion. This model is
 * for control-system demonstration only and must not be used for dosing.
 */

const assert = require('node:assert/strict');

const P = Object.freeze({
  weightKg: 70,
  p1: 0.024,
  p2: 0.02093,
  p3: 2.655e-6,
  n: 0.3,
  p6: 0.001675,
  p5: 119,
  Gb: 117,
  Ib: 12,
  mealTau: 40,
  VI: 0.12,       // L/kg; actuator distribution volume
  tmaxI: 55,      // min; two-compartment subcutaneous absorption
});

const C = Object.freeze({
  target: 100,
  period: 1,
  sensorTau: 10,
  kp: 0.035,      // U/h per mg/dL
  ti: 300,        // min; ki = kp/ti
  maxUph: 10,
  suspend: 75,
});

const ki = C.kp / C.ti;

function targetEquilibrium(target = C.target) {
  const G = target;
  const X = P.p1 * (P.Gb - G) / G;
  const I = P.Ib + P.p2 * X / P.p3;
  const appearance = P.n * (I - P.Ib); // µU/mL/min = mU/L/min
  const uPlant = appearance * P.VI;     // mU/kg/min
  const basalUph = uPlant * 60 * P.weightKg / 1000;
  const S = uPlant * P.tmaxI;
  return { G, X, I, S1: S, S2: S, Gs: G, uPlant, basalUph };
}

const EQ = targetEquilibrium();

function mealD(t, meals) {
  let total = 0;
  for (const meal of meals) {
    const age = t - meal.time;
    if (age >= 0) total += meal.amplitude * age * Math.exp(-age / P.mealTau) /
      (P.mealTau ** 2);
  }
  return total;
}

function rateToPlant(rateUph) {
  return rateUph * 1000 / (60 * P.weightKg);
}

function derivative(t, y, meals, rateUph) {
  const [G, X, I, S1, S2, Gs, tau] = y;
  const u = rateToPlant(rateUph);
  const tauRate = G > P.p5 ? 1 : 0;
  const endogenous = P.p6 * Math.max(G - P.p5, 0) * tau;
  return [
    -P.p1 * (G - P.Gb) - X * G + mealD(t, meals),
    -P.p2 * X + P.p3 * (I - P.Ib),
    endogenous - P.n * (I - P.Ib) + S2 / (P.tmaxI * P.VI),
    u - S1 / P.tmaxI,
    S1 / P.tmaxI - S2 / P.tmaxI,
    (G - Gs) / C.sensorTau,
    tauRate,
  ];
}

function addScaled(y, k, scale) {
  return y.map((value, index) => value + scale * k[index]);
}

function rk4(t, y, h, meals, rateUph) {
  const k1 = derivative(t, y, meals, rateUph);
  const k2 = derivative(t + h / 2, addScaled(y, k1, h / 2), meals, rateUph);
  const k3 = derivative(t + h / 2, addScaled(y, k2, h / 2), meals, rateUph);
  const k4 = derivative(t + h, addScaled(y, k3, h), meals, rateUph);
  const next = y.map((value, index) => value + h *
    (k1[index] + 2 * k2[index] + 2 * k3[index] + k4[index]) / 6);
  if (next[0] <= P.p5) next[6] = 0;
  return next;
}

function newController({ controlledEquilibrium = false } = {}) {
  return {
    integral: controlledEquilibrium ? 0 : 0,
    rate: controlledEquilibrium ? EQ.basalUph : 0,
  };
}

function controllerStep(controller, glucose) {
  const error = glucose - C.target;
  const unsaturated = EQ.basalUph + C.kp * error + ki * controller.integral;
  let rate = Math.max(0, Math.min(C.maxUph, unsaturated));
  if (glucose <= C.suspend) rate = 0;

  const saturatedHigh = unsaturated > C.maxUph && error > 0;
  const saturatedLow = unsaturated < 0 && error < 0;
  if (!saturatedHigh && !saturatedLow && glucose > C.suspend) {
    controller.integral += error * C.period;
  }
  controller.rate = rate;
  return { rate, error, unsaturated };
}

function simulate({ duration = 1440, dt = 0.1, meals = [], control = true,
  startAtTarget = false } = {}) {
  let y = startAtTarget
    ? [EQ.G, EQ.X, EQ.I, EQ.S1, EQ.S2, EQ.Gs, 0]
    : [P.Gb, 0, P.Ib, 0, 0, P.Gb, 0];
  const controller = newController({ controlledEquilibrium: startAtTarget });
  let nextControl = 0;
  let rate = startAtTarget ? EQ.basalUph : 0;
  const rows = [{ t: 0, G: y[0], X: y[1], I: y[2], Gs: y[5], rate }];

  for (let t = 0; t < duration - 1e-9;) {
    if (control && t >= nextControl - 1e-9) {
      rate = controllerStep(controller, y[5]).rate;
      nextControl += C.period;
    }
    const h = Math.min(dt, duration - t, control ? nextControl - t : dt);
    y = rk4(t, y, h, meals, control ? rate : 0);
    t += h;
    assert.ok(y.every(Number.isFinite), `non-finite state at ${t}`);
    assert.ok(y.slice(0, 6).every(value => value >= 0), `negative state at ${t}`);
    rows.push({ t, G: y[0], X: y[1], I: y[2], Gs: y[5], rate });
  }
  return rows;
}

function metrics(rows, from = 0) {
  const selected = rows.filter(row => row.t >= from);
  const peak = selected.reduce((best, row) => row.G > best.G ? row : best);
  const nadir = selected.reduce((best, row) => row.G < best.G ? row : best);
  const final = rows.at(-1);
  return {
    peakG: peak.G,
    peakTime: peak.t,
    nadirG: nadir.G,
    nadirTime: nadir.t,
    finalG: final.G,
    finalSensorG: final.Gs,
    finalRateUph: final.rate,
    maxRateUph: Math.max(...selected.map(row => row.rate)),
    minutesBelow70: selected.filter(row => row.G < 70).length *
      (selected.length > 1 ? selected[1].t - selected[0].t : 0),
  };
}

const equilibrium = simulate({ duration: 1000, startAtTarget: true });
const equilibriumDeviation = Math.max(...equilibrium.map(row => Math.max(
  Math.abs(row.G - EQ.G), Math.abs(row.Gs - EQ.G), Math.abs(row.rate - EQ.basalUph))));
assert.ok(equilibriumDeviation < 1e-10, 'controlled equilibrium drifts');

const uncontrolled = simulate({ control: false });
const controlled = simulate({ control: true });
const uncontrolledMetrics = metrics(uncontrolled);
const controlledMetrics = metrics(controlled);
assert.ok(Math.abs(uncontrolledMetrics.finalG - P.Gb) < 1e-8);
assert.ok(Math.abs(controlledMetrics.finalG - C.target) < 1,
  'PI does not bring fasting T2D glucose close to target');
assert.ok(controlledMetrics.nadirG >= 70, 'PI causes hypoglycaemia during target acquisition');

const mealScenarios = [80, 150, 220].map(amplitude => ({
  amplitude,
  ...metrics(simulate({
    control: true,
    duration: 2400,
    meals: [{ time: 360, amplitude }],
  }), 360),
}));
assert.ok(mealScenarios.every(row => row.nadirG >= 70),
  'PI meal scenario causes hypoglycaemia');
assert.ok(mealScenarios.every(row => Math.abs(row.finalG - C.target) < 2),
  'PI meal scenario does not recover');
assert.ok(mealScenarios[0].peakG < mealScenarios[1].peakG &&
  mealScenarios[1].peakG < mealScenarios[2].peakG,
  'PI meal dose response is not monotonic');

const saturationProbe = newController();
for (let i = 0; i < 1000; i++) controllerStep(saturationProbe, 400);
const integralAtSaturation = saturationProbe.integral;
for (let i = 0; i < 100; i++) controllerStep(saturationProbe, 400);
assert.equal(saturationProbe.rate, C.maxUph, 'upper saturation failed');
assert.equal(saturationProbe.integral, integralAtSaturation, 'anti-windup failed');
const suspendProbe = newController();
assert.equal(controllerStep(suspendProbe, 70).rate, 0, 'low-glucose suspend failed');

const convergence = [0.5, 0.25, 0.1].map(dt => ({
  dt,
  ...metrics(simulate({
    control: true,
    duration: 2400,
    meals: [{ time: 360, amplitude: 150 }],
    dt,
  }), 360),
}));
assert.ok(convergence.slice(0, -1).every(row =>
  Math.abs(row.peakG - convergence.at(-1).peakG) < 0.1 &&
  Math.abs(row.finalG - convergence.at(-1).finalG) < 0.1));

console.log(JSON.stringify({
  status: 'PASS',
  model: 'PI-controlled T2D G-X-I educational v1',
  targetEquilibrium: EQ,
  controller: { ...C, ki },
  equilibriumDeviation,
  uncontrolled: uncontrolledMetrics,
  controlledNoMeal: controlledMetrics,
  controlledMeals: mealScenarios,
  protectionChecks: {
    upperSaturationUph: saturationProbe.rate,
    antiWindupIntegral: saturationProbe.integral,
    suspendRateUph: suspendProbe.rate,
  },
  convergence,
  limitations: [
    'controller and actuator are project models, not a medical device algorithm',
    'gain values are numerically tuned for one virtual subject',
    'output must not be interpreted as an insulin dose recommendation',
  ],
}, null, 2));
