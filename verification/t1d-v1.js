'use strict';

/*
 * Reproducible numerical checks for the proposed five-state T1D v1 model.
 * This script is independent from the browser simulator.
 */

const P = {
  p1: 0.03082,
  p2: 0.02093,
  p3: 1.062e-5,
  Gb: 92,
  Ib: 7.3,
  mealTau: 40,
  bodyWeight: 70,
  VIperKg: 0.12,
  ke: 0.138,
  tmaxI: 55,
};
P.VI = P.VIperKg * P.bodyWeight;
P.uBasal = P.Ib * P.VI * P.ke;
P.Sbasal = P.uBasal * P.tmaxI;

function mealInput(t, meals) {
  let total = 0;
  for (const meal of meals) {
    const age = t - meal.time;
    if (age >= 0) {
      total += meal.A * age * Math.exp(-age / P.mealTau) /
        (P.mealTau * P.mealTau);
    }
  }
  return total;
}

function derivative(t, y, cfg) {
  const [G, X, I, S1, S2] = y;
  const basal = t >= (cfg.stopBasalAt ?? Infinity) ? 0 : P.uBasal;
  return [
    -P.p1 * (G - P.Gb) - X * G + mealInput(t, cfg.meals),
    -P.p2 * X + P.p3 * (I - P.Ib),
    (S2 / P.tmaxI) / P.VI - P.ke * I,
    basal - S1 / P.tmaxI,
    (S1 - S2) / P.tmaxI,
  ];
}

function rk4Step(t, y, h, cfg) {
  const k1 = derivative(t, y, cfg);
  const y2 = y.map((v, i) => v + h * k1[i] / 2);
  const k2 = derivative(t + h / 2, y2, cfg);
  const y3 = y.map((v, i) => v + h * k2[i] / 2);
  const k3 = derivative(t + h / 2, y3, cfg);
  const y4 = y.map((v, i) => v + h * k3[i]);
  const k4 = derivative(t + h, y4, cfg);
  return y.map((v, i) =>
    v + h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6);
}

function simulate({
  dt = 0.5,
  duration = 1000,
  initial = [P.Gb, 0, P.Ib, P.Sbasal, P.Sbasal],
  meals = [],
  boluses = [],
  stopBasalAt,
}) {
  const cfg = { meals, stopBasalAt };
  let y = [...initial];
  let bolusIndex = 0;
  const trace = [];
  const steps = Math.round(duration / dt);

  for (let step = 0; step <= steps; step++) {
    const t = step * dt;
    while (bolusIndex < boluses.length &&
           Math.abs(boluses[bolusIndex].time - t) < dt / 10) {
      y[3] += boluses[bolusIndex].units * 1000;
      bolusIndex++;
    }
    trace.push({ t, G: y[0], X: y[1], I: y[2], S1: y[3], S2: y[4] });
    if (step < steps) y = rk4Step(t, y, dt, cfg);
  }
  return trace;
}

function metrics(trace, eventTime = 0) {
  let maxG = trace[0], minG = trace[0], maxI = trace[0];
  for (const row of trace) {
    if (row.G > maxG.G) maxG = row;
    if (row.G < minG.G) minG = row;
    if (row.I > maxI.I) maxI = row;
  }
  let lastOutside = eventTime;
  for (const row of trace) {
    if (row.t >= eventTime &&
        (Math.abs(row.G - P.Gb) >= 1 || Math.abs(row.I - P.Ib) >= 0.1)) {
      lastOutside = row.t;
    }
  }
  const end = trace.at(-1);
  return {
    peakG: maxG.G,
    peakGTime: maxG.t,
    minG: minG.G,
    minGTime: minG.t,
    peakI: maxI.I,
    peakITime: maxI.t,
    finalG: end.G,
    finalI: end.I,
    settleTime: lastOutside < end.t ? lastOutside : null,
  };
}

function maxBaselineDeviation(trace) {
  return trace.reduce((m, r) => ({
    G: Math.max(m.G, Math.abs(r.G - P.Gb)),
    X: Math.max(m.X, Math.abs(r.X)),
    I: Math.max(m.I, Math.abs(r.I - P.Ib)),
    S1: Math.max(m.S1, Math.abs(r.S1 - P.Sbasal)),
    S2: Math.max(m.S2, Math.abs(r.S2 - P.Sbasal)),
  }), { G: 0, X: 0, I: 0, S1: 0, S2: 0 });
}

const equilibrium = simulate({});
const perturbation = simulate({
  initial: [P.Gb + 10, 0, P.Ib + 1, P.Sbasal, P.Sbasal],
});
const meal = [{ time: 60, A: 150 }];
const doseRuns = [0, 2, 4, 6, 8].map(units => ({
  units,
  metrics: metrics(simulate({
    meals: meal,
    boluses: units ? [{ time: 60, units }] : [],
  }), 60),
}));
const stopped = simulate({ duration: 2000, stopBasalAt: 60 });
const convergence = [0.5, 0.25, 0.1].map(dt => ({
  dt,
  metrics: metrics(simulate({
    dt,
    meals: meal,
    boluses: [{ time: 60, units: 4 }],
  }), 60),
}));

const xNoInsulin = -P.p3 * P.Ib / P.p2;
const analyticNoInsulinG = P.p1 * P.Gb / (P.p1 + xNoInsulin);

console.log(JSON.stringify({
  parameters: {
    VI_L: P.VI,
    tmaxI_min: P.tmaxI,
    ke_per_min: P.ke,
    basal_mU_per_min: P.uBasal,
    basal_U_per_hour: P.uBasal * 60 / 1000,
    S1b_mU: P.Sbasal,
    S2b_mU: P.Sbasal,
  },
  analytic: {
    eigenvalues: [-P.p1, -P.p2, -P.ke, -1 / P.tmaxI, -1 / P.tmaxI],
    noInsulinX: xNoInsulin,
    noInsulinG_mg_dL: analyticNoInsulinG,
  },
  equilibriumMaxDeviation: maxBaselineDeviation(equilibrium),
  perturbation: metrics(perturbation),
  mealDoseResponse: doseRuns,
  basalStopped: metrics(stopped, 60),
  timeStepConvergence: convergence,
}, null, 2));
