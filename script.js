/* =========================================================
   CONFIG — satu sumber kebenaran untuk semua parameter
   ========================================================= */
const CONFIG = {
  p1: 0.03082, p2: 0.02093, p3: 1.062e-5,
  n: 0.3, p6: 0.003349, p5: 94.0,
  Gb: 92.0, Ib: 7.3,
  dt: 0.5,                 // menit, step internal solver
  minutesPerSecond: 22,    // kecepatan simulasi (dipercepat)
  mealAmplitude: { kecil: 80, sedang: 150, besar: 220 },
  mealTau: 40,
  zones: { hipoMax: 70, normalMax: 180, chartMax: 260 },
  historyWindowMin: 480,   // jendela auto-scroll grafik (menit)
  historyStepMin: 1.5,     // jarak antar titik yang disimpan ke history
};

/* =========================================================
   MODE REGISTRY — kontrak seragam, siap upgrade
   ========================================================= */
const ModeNormal = {
  id: 'normal', label: 'Normal', sub: 'Pankreas sehat', enabled: true,
  computeIin(state, tauLocal, p) {
    const above = Math.max(state.G - p.p5, 0);
    return p.p6 * above * tauLocal;
  },
  extraControls: null,
};
const ModeDiabetes = {
  id: 'diabetes', label: 'Diabetes', sub: 'Tanpa kontrol', enabled: true,
  computeIin: null,   // T1D: engine Hovorka 8-state, bukan computeIin
  engine: 't1d',
  extraControls: null,
};
const ModeDiabetesAP = {
  id: 'diabetes_ap', label: 'Diabetes + AP', sub: 'Artificial Pancreas', enabled: true,
  computeIin: null,   // AP: plant T1D + controller PID-IFB
  engine: 'ap',
  extraControls: null,
};
const ModeDiabetesT2D = {
  id: 'diabetes_t2d', label: 'Diabetes Tipe 2', sub: 'Resistensi + sel-beta melemah', enabled: true,
  computeIin(state, tauLocal, p) {
    return p.p6 * Math.max(state.G - p.p5, 0) * tauLocal;
  },
  engine: 't2d',
  extraControls: null,
};
const ModeDiabetesT2DPI = {
  id: 'diabetes_t2d_pi', label: 'T2D + PI', sub: 'Kontrol insulin otomatis', enabled: true,
  computeIin: null,   // T2D+PI: plant G-X-I + aktuator subkutan + controller PI
  engine: 't2d_pi',
  extraControls: null,
};
const MODES = [ModeNormal, ModeDiabetes, ModeDiabetesAP, ModeDiabetesT2D, ModeDiabetesT2DPI];
let currentMode = ModeNormal;

/* =========================================================
   T1D ENGINE — Hovorka 8-state (port dari verification/t1d-v2.js,
   spec: sumber/model-normal-glucose-insulin.md Bagian 10)
   ========================================================= */
const T1D_PARAMS = {
  weightKg: 70,
  targetG: 5.5,            // mmol/L (99.1 mg/dL) — equilibrium T1D manual
  VG: 0.16, VI: 0.12,
  k12: 0.066,
  ka1: 0.006, ka2: 0.060, ka3: 0.030,
  ke: 0.138,
  SIT: 51.2e-4, SID: 8.2e-4, SIE: 520e-4,
  EGP0: 0.0161, F01: 0.0097,
  AG: 0.8, tmaxG: 40, tmaxI: 55,
  mmolToMgDl: 18.0182,
};
T1D_PARAMS.kb1 = T1D_PARAMS.ka1 * T1D_PARAMS.SIT;
T1D_PARAMS.kb2 = T1D_PARAMS.ka2 * T1D_PARAMS.SID;
T1D_PARAMS.kb3 = T1D_PARAMS.ka3 * T1D_PARAMS.SIE;

const MEAL_GRAMS = { kecil: 25, sedang: 50, besar: 75 };   // preset T1D/AP (Bagian 11.13)

function t1dF01c(G){ return G >= 4.5 ? T1D_PARAMS.F01 : T1D_PARAMS.F01 * G / 4.5; }
function t1dRenal(G){ return G >= 9 ? 0.003 * (G - 9) * T1D_PARAMS.VG : 0; }
function t1dCarbsToMmolKg(grams){ return grams * 1000 / (180.16 * T1D_PARAMS.weightKg); }

function t1dMealAppearance(t, meals){
  let total = 0;
  for(const m of meals){
    const age = t - m.time;
    if(age >= 0){
      total += t1dCarbsToMmolKg(m.grams) * T1D_PARAMS.AG * age *
        Math.exp(-age / T1D_PARAMS.tmaxG) / (T1D_PARAMS.tmaxG * T1D_PARAMS.tmaxG);
    }
  }
  return total;
}

function bisect(fn, lo, hi, tolerance = 1e-13){
  let flo = fn(lo), fhi = fn(hi);
  if(flo * fhi > 0) throw new Error('Root is not bracketed');
  while(hi - lo > tolerance){
    const mid = (lo + hi) / 2;
    const fm = fn(mid);
    if(flo * fm <= 0){ hi = mid; fhi = fm; }
    else { lo = mid; flo = fm; }
  }
  return (lo + hi) / 2;
}

/* Bagian 10.10: equilibrium 8-state via bisection — WAJIB untuk initial state */
function t1dEquilibrium(targetG = T1D_PARAMS.targetG){
  const residual = (I) => {
    const x1 = T1D_PARAMS.SIT * I, x2 = T1D_PARAMS.SID * I, x3 = T1D_PARAMS.SIE * I;
    const Q1 = T1D_PARAMS.VG * targetG;
    const Q2 = x1 * Q1 / (T1D_PARAMS.k12 + x2);
    return -t1dF01c(targetG) - x1 * Q1 + T1D_PARAMS.k12 * Q2 -
      t1dRenal(targetG) + T1D_PARAMS.EGP0 * Math.max(0, 1 - x3);
  };
  const I = bisect(residual, 0, 1 / T1D_PARAMS.SIE);
  const x1 = T1D_PARAMS.SIT * I, x2 = T1D_PARAMS.SID * I, x3 = T1D_PARAMS.SIE * I;
  const Q1 = T1D_PARAMS.VG * targetG;
  const Q2 = x1 * Q1 / (T1D_PARAMS.k12 + x2);
  const uBasal = I * T1D_PARAMS.VI * T1D_PARAMS.ke;   // mU/kg/min
  const S1 = uBasal * T1D_PARAMS.tmaxI;
  return {
    y: [Q1, Q2, S1, S1, I, x1, x2, x3],
    uBasal,
    basalUph: uBasal * T1D_PARAMS.weightKg * 60 / 1000,  // U/hour
  };
}

/* Bolus adalah loncatan massa depot S1, sekali per event, bukan laju RK4. */
function applyPendingBoluses(state){
  state.boluses = state.boluses.filter(b => {
    if(b.time > state.t + 1e-9) return true;
    state.y[2] += b.units * 1000 / T1D_PARAMS.weightKg;
    return false;
  });
}

function t1dDerivative(t, y, cfg){
  const [Q1, Q2, S1, S2, I, x1, x2, x3] = y;
  const G = Q1 / T1D_PARAMS.VG;
  const basal = cfg.basalOn ? cfg.uBasal : 0;
  const u = basal;
  const UG = t1dMealAppearance(t, cfg.meals);
  return [
    -t1dF01c(G) - x1 * Q1 + T1D_PARAMS.k12 * Q2 - t1dRenal(G) +
      UG + T1D_PARAMS.EGP0 * Math.max(0, 1 - x3),      // clipEGP (Bagian 10.5)
    x1 * Q1 - (T1D_PARAMS.k12 + x2) * Q2,
    u - S1 / T1D_PARAMS.tmaxI,
    S1 / T1D_PARAMS.tmaxI - S2 / T1D_PARAMS.tmaxI,
    (S2 / T1D_PARAMS.tmaxI) / T1D_PARAMS.VI - T1D_PARAMS.ke * I,
    -T1D_PARAMS.ka1 * x1 + T1D_PARAMS.kb1 * I,
    -T1D_PARAMS.ka2 * x2 + T1D_PARAMS.kb2 * I,
    -T1D_PARAMS.ka3 * x3 + T1D_PARAMS.kb3 * I,
  ];
}

function t1dAddScaled(y, k, scale){ return y.map((v, i) => v + scale * k[i]); }

function t1dRk4Step(t, y, h, cfg){
  const k1 = t1dDerivative(t, y, cfg);
  const k2 = t1dDerivative(t + h/2, t1dAddScaled(y, k1, h/2), cfg);
  const k3 = t1dDerivative(t + h/2, t1dAddScaled(y, k2, h/2), cfg);
  const k4 = t1dDerivative(t + h, t1dAddScaled(y, k3, h), cfg);
  return y.map((v, i) => v + h * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6);
}

/* =========================================================
   AP CONTROLLER — PID + insulin feedback Ruiz (port dari
   verification/ap-pid-ifb-v1.js, spec Bagian 11)
   ========================================================= */
const IFB_PARAMS = {
  a11: 0.9802, a21: 0.014043, a31: 0.000127,
  a22: 0.98582, a32: 0.017889, a33: 0.98198,
  b1: 1.1881, b2: 0.0084741, b3: 0.00005,
  g1: 0.64935, g2: 0.34128, g3: 0.0093667,
};
IFB_PARAMS.gSum = IFB_PARAMS.g1 + IFB_PARAMS.g2 + IFB_PARAMS.g3;

const AP_CONFIG = {
  targetMgDl: 120,
  controllerPeriod: 1,      // menit
  Ti: 300,                  // tuning proyek (dokumen 11.9), bukan 150 paper
  Td: 75,
  derivativeTau: 5,
  sensorTau: 10,
  suspendMgDl: 70,
};

/* Bagian 11.6: equilibrium estimator IFB dari basal rate */
function ifbEquilibrium(rateUph){
  const dose = rateUph / 60;   // U/jam → U/menit (dokumen 11.13 temuan #2)
  const Isc = IFB_PARAMS.b1 * dose / (1 - IFB_PARAMS.a11);
  const Ip = (IFB_PARAMS.a21 * Isc + IFB_PARAMS.b2 * dose) / (1 - IFB_PARAMS.a22);
  const Ieff = (IFB_PARAMS.a31 * Isc + IFB_PARAMS.a32 * Ip +
    IFB_PARAMS.b3 * dose) / (1 - IFB_PARAMS.a33);
  return [Isc, Ip, Ieff];
}

function ifbStep(state, deliveredRateUph){
  const [Isc, Ip, Ieff] = state;
  const dose = deliveredRateUph / 60;
  return [
    IFB_PARAMS.a11 * Isc + IFB_PARAMS.b1 * dose,
    IFB_PARAMS.a21 * Isc + IFB_PARAMS.a22 * Ip + IFB_PARAMS.b2 * dose,
    IFB_PARAMS.a31 * Isc + IFB_PARAMS.a32 * Ip + IFB_PARAMS.a33 * Ieff + IFB_PARAMS.b3 * dose,
  ];
}

function ifbValue(state){
  return IFB_PARAMS.g1 * state[0] + IFB_PARAMS.g2 * state[1] + IFB_PARAMS.g3 * state[2];
}

/* Plant equilibrium AP pada target 120 mg/dL — 9 state (y[8] = Gs sensor) */
function apPlantEquilibrium(){
  const G = AP_CONFIG.targetMgDl / T1D_PARAMS.mmolToMgDl;
  const residual = (I) => {
    const x1 = T1D_PARAMS.SIT * I, x2 = T1D_PARAMS.SID * I, x3 = T1D_PARAMS.SIE * I;
    const Q1 = T1D_PARAMS.VG * G;
    const Q2 = x1 * Q1 / (T1D_PARAMS.k12 + x2);
    return -t1dF01c(G) - x1 * Q1 + T1D_PARAMS.k12 * Q2 - t1dRenal(G) +
      T1D_PARAMS.EGP0 * Math.max(0, 1 - x3);
  };
  const I = bisect(residual, 0, 1 / T1D_PARAMS.SIE);
  const x1 = T1D_PARAMS.SIT * I, x2 = T1D_PARAMS.SID * I, x3 = T1D_PARAMS.SIE * I;
  const Q1 = T1D_PARAMS.VG * G;
  const Q2 = x1 * Q1 / (T1D_PARAMS.k12 + x2);
  const basalPlant = I * T1D_PARAMS.VI * T1D_PARAMS.ke;
  const S1 = basalPlant * T1D_PARAMS.tmaxI;
  const basalUph = basalPlant * T1D_PARAMS.weightKg * 60 / 1000;
  return {
    y: [Q1, Q2, S1, S1, I, x1, x2, x3, AP_CONFIG.targetMgDl],
    basalPlant, basalUph,
  };
}

/* Bagian 11.5-11.8: controller dengan bumpless init + anti-windup + suspend */
function newApController(eq){
  const dailyInsulinU = eq.basalUph * 24;
  const Kp = dailyInsulinU / 2250;                  // Ruiz: Kp = IDIR/2250
  const estimator = ifbEquilibrium(eq.basalUph);
  const feedbackEq = ifbValue(estimator);
  return {
    Kp,
    integral: (eq.basalUph + feedbackEq) / (1 + IFB_PARAMS.gSum),  // bumpless (11.13 #1)
    derivative: 0,
    previousSG: AP_CONFIG.targetMgDl,
    estimator,
    deliveredRate: eq.basalUph,
    maxRate: Math.max(3 * eq.basalUph, 3),
    minIntegral: 0,
    maxIntegral: 3 * eq.basalUph,
  };
}

function apControllerOutput(c, SG){
  c.estimator = ifbStep(c.estimator, c.deliveredRate);
  const feedback = ifbValue(c.estimator);
  const error = SG - AP_CONFIG.targetMgDl;
  const rawDerivative = (SG - c.previousSG) / AP_CONFIG.controllerPeriod;
  const alpha = Math.exp(-AP_CONFIG.controllerPeriod / AP_CONFIG.derivativeTau);
  c.derivative = alpha * c.derivative + (1 - alpha) * rawDerivative;

  const candidateIntegral = Math.max(c.minIntegral, Math.min(c.maxIntegral,
    c.integral + c.Kp / AP_CONFIG.Ti * error * AP_CONFIG.controllerPeriod));

  const calculate = (integral) => {
    const proportional = c.Kp * error;
    const derivative = c.Kp * AP_CONFIG.Td * c.derivative;
    const pidCore = proportional + integral + derivative;
    const raw = (1 + IFB_PARAMS.gSum) * pidCore - feedback;
    return { proportional, integral, derivative, feedback, raw };
  };

  let terms = calculate(candidateIntegral);
  let rate = Math.max(0, Math.min(c.maxRate, terms.raw));
  const saturatedHigh = terms.raw > c.maxRate && error > 0;
  const saturatedLow = terms.raw < 0 && error < 0;
  if(saturatedHigh || saturatedLow){
    terms = calculate(c.integral);      // anti-windup: integral dibekukan
    rate = Math.max(0, Math.min(c.maxRate, terms.raw));
  } else {
    c.integral = candidateIntegral;
  }

  let suspended = false;
  if(SG <= AP_CONFIG.suspendMgDl){       // low-glucose suspend (11.8)
    rate = 0;
    suspended = true;
  }
  c.previousSG = SG;
  c.deliveredRate = rate;
  return { ...terms, rate, error, suspended };
}

/* Plant AP derivative — 9 state (Gs sensor lag orde 1, Bagian 11.3) */
function apDerivative(t, y, cfg, rateUph){
  const [Q1, Q2, S1, S2, I, x1, x2, x3, Gs] = y;
  const G = Q1 / T1D_PARAMS.VG;
  const GmgDl = G * T1D_PARAMS.mmolToMgDl;
  const u = rateUph * 1000 / (60 * T1D_PARAMS.weightKg);   // U/h → mU/kg/min
  const UG = t1dMealAppearance(t, cfg.meals);
  return [
    -t1dF01c(G) - x1 * Q1 + T1D_PARAMS.k12 * Q2 - t1dRenal(G) +
      UG + T1D_PARAMS.EGP0 * Math.max(0, 1 - x3),
    x1 * Q1 - (T1D_PARAMS.k12 + x2) * Q2,
    u - S1 / T1D_PARAMS.tmaxI,
    S1 / T1D_PARAMS.tmaxI - S2 / T1D_PARAMS.tmaxI,
    (S2 / T1D_PARAMS.tmaxI) / T1D_PARAMS.VI - T1D_PARAMS.ke * I,
    -T1D_PARAMS.ka1 * x1 + T1D_PARAMS.kb1 * I,
    -T1D_PARAMS.ka2 * x2 + T1D_PARAMS.kb2 * I,
    -T1D_PARAMS.ka3 * x3 + T1D_PARAMS.kb3 * I,
    (GmgDl - Gs) / AP_CONFIG.sensorTau,
  ];
}

function apRk4Step(t, y, h, cfg, rate){
  const k1 = apDerivative(t, y, cfg, rate);
  const k2 = apDerivative(t + h/2, t1dAddScaled(y, k1, h/2), cfg, rate);
  const k3 = apDerivative(t + h/2, t1dAddScaled(y, k2, h/2), cfg, rate);
  const k4 = apDerivative(t + h, t1dAddScaled(y, k3, h), cfg, rate);
  return y.map((v, i) => v + h * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6);
}

/* =========================================================
   T2D ENGINE — persamaan G-X-I Mode Normal dengan parameter T2D
   (verification/t2d-normal-derived-v1.js, dokumen Bagian 13)
   ========================================================= */
const T2D_PARAMS = {
  p1: 0.024,
  p2: 0.02093,
  p3: 2.655e-6,             // 25% Normal; severity setting proyek
  n: 0.3,
  p6: 0.001675,             // 50% Normal; severity setting proyek
  p5: 119,
  Gb: 117,
  Ib: 12,
};
// Parameter input T2D yang tidak termasuk state diferensial G-X-I.
const T2D_UI_PARAMS = { resistanceScale: 1, mealTau: 40 };
const T2D_CONTROL_TARGET = 100; // mg/dL; setpoint kandidat controller T2D berikutnya

const T2D_RESISTANCE_MIN = 0.5, T2D_RESISTANCE_MAX = 2.0;

/* =========================================================
   T2D + PI ENGINE — G-X-I + aktuator subkutan + sensor lag
   (port dari verification/t2d-pi-v1.js, spec Bagian 14)
   ========================================================= */
const T2D_PI_ACTUATOR = {
  weightKg: 70,
  VI: 0.12,       // L/kg; volume distribusi aktuator
  tmaxI: 55,      // min; absorpsi subkutan 2-kompartemen
  mealTau: 40,
};

const T2D_PI_CONFIG = {
  target: 100,    // mg/dL (kandidat proyek, Bagian 13.5/14.3)
  period: 1,      // menit
  sensorTau: 10,  // menit
  kp: 0.035,      // U/h per mg/dL — tuning numerik proyek
  ti: 300,        // menit — ki = kp/ti
  maxUph: 10,     // U/jam — saturasi pompa
  suspend: 75,    // mg/dL — low-glucose suspend
};
function t2dPiKi(){ return T2D_PI_CONFIG.kp / T2D_PI_CONFIG.ti; }

/* Bagian 14.4: equilibrium target dihitung dari persamaan plant, bukan tebakan */
function t2dPiTargetEquilibrium(target = T2D_PI_CONFIG.target, resistanceScale = 1){
  const G = target;
  const X = T2D_PARAMS.p1 * (T2D_PARAMS.Gb - G) / G;
  const p3eff = T2D_PARAMS.p3 / resistanceScale;
  const I = T2D_PARAMS.Ib + T2D_PARAMS.p2 * X / p3eff;
  const appearance = T2D_PARAMS.n * (I - T2D_PARAMS.Ib);   // µU/mL/min
  const uPlant = appearance * T2D_PI_ACTUATOR.VI;           // mU/kg/min
  const basalUph = uPlant * 60 * T2D_PI_ACTUATOR.weightKg / 1000;
  const S = uPlant * T2D_PI_ACTUATOR.tmaxI;
  return { G, X, I, S1: S, S2: S, Gs: G, uPlant, basalUph };
}

function t2dPiRateToPlant(rateUph){
  return rateUph * 1000 / (60 * T2D_PI_ACTUATOR.weightKg);
}

/* Input makan: format amplitude (Mode Normal style) — plant G-X-I menerima D(t) */
function t2dPiMealD(t, meals){
  let total = 0;
  for(const m of meals){
    const age = t - m.startTime;
    if(age >= 0){
      total += m.A * age * Math.exp(-age / m.tauMeal) /
        (m.tauMeal ** 2);
    }
  }
  return total;
}

/* 7-state: [G, X, I, S1, S2, Gs, tau] — tau = tau lokal endogen dalam RK4 */
function t2dPiDerivative(t, y, meals, rateUph, p3eff){
  const [G, X, I, S1, S2, Gs, tau] = y;
  const u = t2dPiRateToPlant(rateUph);
  const tauRate = G > T2D_PARAMS.p5 ? 1 : 0;
  const endogenous = T2D_PARAMS.p6 * Math.max(G - T2D_PARAMS.p5, 0) * tau;
  return [
    -T2D_PARAMS.p1 * (G - T2D_PARAMS.Gb) - X * G + t2dPiMealD(t, meals),
    -T2D_PARAMS.p2 * X + p3eff * (I - T2D_PARAMS.Ib),
    endogenous - T2D_PARAMS.n * (I - T2D_PARAMS.Ib) +
      S2 / (T2D_PI_ACTUATOR.tmaxI * T2D_PI_ACTUATOR.VI),
    u - S1 / T2D_PI_ACTUATOR.tmaxI,
    S1 / T2D_PI_ACTUATOR.tmaxI - S2 / T2D_PI_ACTUATOR.tmaxI,
    (G - Gs) / T2D_PI_CONFIG.sensorTau,
    tauRate,
  ];
}

function t2dPiRk4Step(t, y, h, meals, rateUph, p3eff){
  const k1 = t2dPiDerivative(t, y, meals, rateUph, p3eff);
  const k2 = t2dPiDerivative(t + h/2, t1dAddScaled(y, k1, h/2), meals, rateUph, p3eff);
  const k3 = t2dPiDerivative(t + h/2, t1dAddScaled(y, k2, h/2), meals, rateUph, p3eff);
  const k4 = t2dPiDerivative(t + h, t1dAddScaled(y, k3, h), meals, rateUph, p3eff);
  const next = y.map((v, i) => v + h * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6);
  if(next[0] <= T2D_PARAMS.p5) next[6] = 0;   // reset tau lokal (verifier line 92)
  return next;
}

/* Bagian 14.3: PI diskrit + saturasi + anti-windup + suspend */
function t2dPiNewController(){
  return { integral: 0, rate: 0 };
}

function t2dPiControllerStep(controller, glucose, eq){
  const error = glucose - T2D_PI_CONFIG.target;
  const unsaturated = eq.basalUph + T2D_PI_CONFIG.kp * error +
    t2dPiKi() * controller.integral;
  let rate = Math.max(0, Math.min(T2D_PI_CONFIG.maxUph, unsaturated));
  if(glucose <= T2D_PI_CONFIG.suspend) rate = 0;

  const saturatedHigh = unsaturated > T2D_PI_CONFIG.maxUph && error > 0;
  const saturatedLow = unsaturated < 0 && error < 0;
  if(!saturatedHigh && !saturatedLow && glucose > T2D_PI_CONFIG.suspend){
    controller.integral += error * T2D_PI_CONFIG.period;   // conditional integration
  }
  controller.rate = rate;
  return { rate, error, unsaturated };
}

/* =========================================================
   STATE
   ========================================================= */
function initialState(){
  if(currentMode.id === 'diabetes_t2d_pi'){
    // T2D+PI: start dari baseline T2D 117 (verifier startAtTarget:false) —
    // PI menarik glukosa menuju target 100 secara dinamis (edukatif).
    const eq = t2dPiTargetEquilibrium();
    return {
      t: 0,
      engine: 't2d_pi',
      y: [T2D_PARAMS.Gb, 0, T2D_PARAMS.Ib, 0, 0, T2D_PARAMS.Gb, 0],
      piEq: eq,                            // { basalUph, ... } untuk readout bias
      piController: t2dPiNewController(),
      piRate: 0,                           // U/h — live readout
      piSuspended: false,
      piInfo: null,                        // { error, unsaturated } readout
      nextControlT: 0,
      resistanceScale: T2D_UI_PARAMS.resistanceScale,
      activeMeals: [],                     // { startTime, A, tauMeal } — Mode Normal style
      mealEvents: [],
      G: T2D_PARAMS.Gb, I: T2D_PARAMS.Ib, X: 0,
      history: [], lastHistoryT: -999,
    };
  }
  if(currentMode.id === 'diabetes_t2d'){
    // T2D: struktur G-X-I Normal, equilibrium memakai parameter T2D (Bagian 13).
    return {
      t: 0,
      engine: 't2d',
      G: T2D_PARAMS.Gb, X: 0, I: T2D_PARAMS.Ib,
      tauLocal: 0,
      resistanceScale: T2D_UI_PARAMS.resistanceScale,
      activeMeals: [],                   // format sama dengan Mode Normal
      mealEvents: [],                    // { t, size } — marker chart
      history: [], lastHistoryT: -999,
    };
  }
  if(currentMode.id === 'diabetes' || currentMode.id === 'diabetes_ap'){
    // T1D / AP: 8-state Hovorka dari equilibrium (Bagian 10.10 — tidak boleh tebakan)
    const isAp = currentMode.id === 'diabetes_ap';
    const eq = isAp ? apPlantEquilibrium() : t1dEquilibrium();
    const initialG = isAp
      ? AP_CONFIG.targetMgDl
      : T1D_PARAMS.targetG * T1D_PARAMS.mmolToMgDl;
    return {
      t: 0,
      engine: currentMode.engine,          // 't1d' | 'ap'
      y: [...eq.y],
      eq,                                   // { y, uBasal, basalUph, basalPlant? }
      basalOn: true,                        // T1D manual: toggle user
      boluses: [],                          // { time, units }
      activeMeals: [],                      // { time, grams }
      mealEvents: [],                       // { t, size } — marker chart
      // AP state (diabaikan engine t1d):
      apController: isAp ? newApController(eq) : null,
      apRate: isAp ? eq.basalUph : 0,       // U/h — live readout
      apSuspended: false,
      apCtrlInfo: null,                     // { error, P, Ic, D, IFB } readout
      nextControlT: 0,                      // tick tepat pada 0, 1, 2, ... menit
      // UI mirror — renderer membaca ini (bridge, tidak refactor):
      G: initialG, I: eq.y[4], X: 0,
      history: [], lastHistoryT: -999,
    };
  }
  return {
    t: 0, G: CONFIG.Gb, X: 0, I: CONFIG.Ib,
    tauLocal: 0,
    activeMeals: [],
    mealEvents: [],         // {t, size} — untuk tanda di grafik, tidak ikut dihapus semudah activeMeals
    history: [],           // {t, G, I, X}
    lastHistoryT: -999,
  };
}
let simState = initialState();
let lastStatusClass = null;   // zona status terakhir (live announce hanya saat berubah)
let lastSummaryT = -1;        // waktu summary chart terakhir di-update

/* =========================================================
   SIMULATION ENGINE (RK4) — mode-agnostic
   ========================================================= */
function mealD(t, activeMeals){
  let D = 0;
  for(const m of activeMeals){
    if(t >= m.startTime){
      const dtSince = t - m.startTime;
      D += m.A * dtSince * Math.exp(-dtSince / m.tauMeal) / (m.tauMeal*m.tauMeal);
    }
  }
  return D;
}

function derivatives(s, t, tauLocal, p3eff, mode, activeMeals, params = CONFIG){
  const { G, X, I } = s;
  const D = mealD(t, activeMeals);
  const dI = mode.computeIin(s, tauLocal, params) - params.n * (I - params.Ib);
  const dG = -params.p1 * (G - params.Gb) - X * G + D;
  const dX = -params.p2 * X + p3eff * (I - params.Ib);
  return { dG, dX, dI };
}

function rk4Step(state, dt, mode, params = CONFIG, p3eff = params.p3){
  const G = state.G;
  // update tau lokal (reset saat G turun <= p5)
  if(G <= params.p5){ state.tauLocal = 0; } else { state.tauLocal += dt; }

  function f(s, tt, tau){
    const d = derivatives(s, tt, tau, p3eff, mode, state.activeMeals, params);
    return { G: d.dG, X: d.dX, I: d.dI };
  }
  const s0 = { G: state.G, X: state.X, I: state.I };
  const k1 = f(s0, state.t, state.tauLocal);
  const s1 = { G: s0.G+dt/2*k1.G, X: s0.X+dt/2*k1.X, I: s0.I+dt/2*k1.I };
  const k2 = f(s1, state.t+dt/2, state.tauLocal);
  const s2 = { G: s0.G+dt/2*k2.G, X: s0.X+dt/2*k2.X, I: s0.I+dt/2*k2.I };
  const k3 = f(s2, state.t+dt/2, state.tauLocal);
  const s3 = { G: s0.G+dt*k3.G, X: s0.X+dt*k3.X, I: s0.I+dt*k3.I };
  const k4 = f(s3, state.t+dt, state.tauLocal);

  state.G = s0.G + (dt/6)*(k1.G+2*k2.G+2*k3.G+k4.G);
  state.X = s0.X + (dt/6)*(k1.X+2*k2.X+2*k3.X+k4.X);
  state.I = s0.I + (dt/6)*(k1.I+2*k2.I+2*k3.I+k4.I);
  state.t += dt;

  // buang meal yang sudah tidak relevan (>6 jam) biar array tidak numpuk
  state.activeMeals = state.activeMeals.filter(m => (state.t - m.startTime) < 360);
  // buang tanda meal yang sudah keluar dari jendela grafik
  state.mealEvents = state.mealEvents.filter(m => (state.t - m.t) < CONFIG.historyWindowMin);
}

/* =========================================================
   SCENARIOS
   ========================================================= */
function triggerMeal(size){
  if(currentMode.id !== 'normal' && currentMode.id !== 'diabetes_t2d' &&
     currentMode.id !== 'diabetes_t2d_pi'){
    // T1D/AP: gram karbohidrat (Bagian 10.6) — preset 25/50/75 g
    simState.activeMeals.push({ time: simState.t, grams: MEAL_GRAMS[size] || 50 });
  } else {
    // Normal/T2D/T2D+PI memakai event makanan G-X-I yang sama.
    simState.activeMeals.push({
      startTime: simState.t,
      A: CONFIG.mealAmplitude[size],
      tauMeal: (currentMode.id === 'diabetes_t2d' || currentMode.id === 'diabetes_t2d_pi')
        ? T2D_UI_PARAMS.mealTau : CONFIG.mealTau,
    });
  }
  simState.mealEvents.push({ t: simState.t, size });
}
function resetSimulation(){
  simState = initialState();
  pendingSimMinutes = 0;
  lastFrameTime = null;
  if(typeof setMealMenuOpen === 'function') setMealMenuOpen(false, false);
  lastStatusClass = null; // paksa render ulang status card setelah reset
  lastSummaryT = -1;
  lastRenderedKey = '';   // paksa chart render 
  lastMathRenderT = -Infinity;
  applyChartScale();      // yG max berganti per mode (260 vs 400)ulang setelah reset
  if(mainChart){
    mainChart.data.labels = [];
    mainChart.data.datasets[0].data = [];
    mainChart.data.datasets[1].data = [];
    mainChart.update('none');
  }
  if(xChart){
    xChart.data.labels = [];
    xChart.data.datasets[0].data = [];
    xChart.update('none');
  }
  if(breakdownChart){
    breakdownChart.data.datasets[0].data = [0,0,0];
    breakdownChart.update('none');
  }
  const sum = document.getElementById('main-chart-summary');
  if(sum) sum.textContent = '';
  const tbody = document.querySelector('#main-chart-data tbody');
  if(tbody) tbody.innerHTML = '';
  updateModePanels();
  updateCards();
  updateApReadout();
  updateT2dPiReadout();
  updateMathPanel();
}

/* =========================================================
   RENDERER: SVG
   ========================================================= */
let vesselPaths = [];
let pancreasPaths = [];

function initSvgRefs(){
  const vesselGroup = document.getElementById('vessel');
  const pancreasGroup = document.getElementById('pancreas');
  const gutGroup = document.getElementById('gut');
  const tissueTint = document.getElementById('tissue-tint');
  if(vesselGroup) vesselPaths = Array.from(vesselGroup.querySelectorAll('path'));
  if(pancreasGroup) pancreasPaths = Array.from(pancreasGroup.querySelectorAll('path'));
  if(gutGroup) gutPaths = Array.from(gutGroup.querySelectorAll('path'));
  tissueTintEl = tissueTint;
}
let gutPaths = [];
let tissueTintEl = null;

function glucoseToColor(G){
  const { hipoMax, normalMax, chartMax } = CONFIG.zones;
  if(G < hipoMax){
    return '#C9463D';
  } else if(G <= normalMax){
    return '#1B9C6E';
  } else {
    const t = Math.min((G - normalMax) / (chartMax - normalMax), 1);
    const r = Math.round(27 + t*(201-27));
    const g = Math.round(156 + t*(70-156));
    const b = Math.round(110 + t*(61-110));
    return `rgb(${r},${g},${b})`;
  }
}

/* Reduced motion: solver & data tetap berjalan penuh, hanya gerak visual
   (partikel & denyut jantung SMIL) yang dibatasi. Audit A5.
   Preferensi dibaca LIVE — berubah tanpa reload (motion spec). */
const reducedMotionQuery =
  window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)');
let prefersReducedMotion = reducedMotionQuery ? reducedMotionQuery.matches : false;

function applyReducedMotionState(){
  const scene = document.getElementById('scene');
  if(prefersReducedMotion){
    if(scene && typeof scene.pauseAnimations === 'function') scene.pauseAnimations();
    document.documentElement.classList.remove('motion-ready', 'motion-entered');
  } else {
    if(scene && typeof scene.unpauseAnimations === 'function' && !simPaused) scene.unpauseAnimations();
  }
}

if(reducedMotionQuery){
  const onChange = () => {
    prefersReducedMotion = reducedMotionQuery.matches;
    applyReducedMotionState();
  };
  if(reducedMotionQuery.addEventListener) reducedMotionQuery.addEventListener('change', onChange);
  else if(reducedMotionQuery.addListener) reducedMotionQuery.addListener(onChange); // browser lama
}

if(prefersReducedMotion){
  window.addEventListener('load', applyReducedMotionState);
}

/* Page entrance: sekali saat load, stagger <500ms, skip jika reduced motion */
(function initPageEntrance(){
  if(prefersReducedMotion) return;
  if(typeof document === 'undefined' || !document.documentElement) return; // non-browser guard
  document.documentElement.classList.add('motion-ready');
  requestAnimationFrame(() => {
    requestAnimationFrame(() => {
      document.documentElement.classList.add('motion-entered');
    });
  });
})();

function updateSvg(dt){
  const color = glucoseToColor(simState.G);
  for(const p of vesselPaths){ p.style.stroke = color; p.style.fill = 'none'; }

  const isT2d = currentMode.id === 'diabetes_t2d' || currentMode.id === 'diabetes_t2d_pi';
  const endocrineParams = isT2d ? T2D_PARAMS : CONFIG;
  const tauLocal = currentMode.id === 'diabetes_t2d_pi'
    ? simState.y[6] : (simState.tauLocal || 0);
  const above = Math.max(simState.G - endocrineParams.p5, 0);
  const secretionRate = endocrineParams.p6 * above * tauLocal;
  const glowStrength = Math.min(secretionRate / 3, 1);
  // Preserve the 2.5D shading with lightweight translucent activity overlays.
  // T1D/AP: pankreas tidak menyekresi — selalu idle. T2D: masih aktif (respons sel-beta).
  const pancreasIdle = (currentMode.id === 'diabetes' || currentMode.id === 'diabetes_ap');
  for(const p of pancreasPaths){
    p.style.opacity = pancreasIdle
      ? 0                               // pankreas nonaktif di T1D/AP
      : (glowStrength * 0.32).toFixed(3);
  }
  const isGxi = currentMode.id === 'normal' || isT2d;
  const D = isGxi
    ? mealD(simState.t, simState.activeMeals)
    : t1dMealAppearance(simState.t, simState.activeMeals);
  const gutGlow = Math.min(D / (isGxi ? VIZ.mealSpawnMax : 0.05), 1);
  for(const p of gutPaths) p.style.opacity = (gutGlow * 0.24).toFixed(3);

  /* Tint jaringan/otot ∝ X(t) — efek insulin aktif menyerap glukosa.
     Normal: X max ~0.009 (XFull). T1D/AP: x2 disposal ~0.005.
     T2D memakai state X yang sama dengan Normal. */
  if(tissueTintEl){
    const xFull = isGxi
      ? VIZ.XFull : 0.005;
    const Xfrac = Math.max(0, Math.min(simState.X / xFull, 1));
    tissueTintEl.style.opacity = (Xfrac * 0.22).toFixed(3);
  }

  /* Reduced motion: partikel membeku di tempat (tidak diperlambat) —
     solver, data, dan status tetap berjalan. (motion spec) */
  if(!prefersReducedMotion){
    updateFlowParticles(dt);
  }
}

/* Anatomi: koordinat visual saja; tidak digunakan oleh solver.
   Pembuluh SVG dan partikel memakai poliline yang sama. */
const TREE = { heartRoot: 19, insulinRoot: 40, nodes: {19:[212,153],40:[200,216],gut:[200,269]}, edges:[], children:{}, edgeWeight:[] };
function addAnatomyEdge(parent, child, pts, weight=1, branch=true){
  const idx = TREE.edges.length;
  TREE.nodes[parent] = pts[0];
  TREE.nodes[child] = pts[pts.length-1];
  TREE.edges.push({parent, child, pts});
  TREE.edgeWeight.push(weight);
  if(branch) (TREE.children[parent] ||= []).push(idx);
  return idx;
}
addAnatomyEdge(19, 1, [[212,153],[201,136],[200,117]], 3);
addAnatomyEdge(1, 2, [[200,117],[189,104],[187,90],[183,70]], 1);
addAnatomyEdge(1, 3, [[200,117],[213,104],[214,88],[218,70]], 1);
addAnatomyEdge(1, 4, [[200,117],[176,113],[144,119],[134,140],[126,164],[116,192]], 2);
addAnatomyEdge(4, 5, [[116,192],[110,215],[97,247],[86,280]], 1);
addAnatomyEdge(5, 6, [[86,280],[77,296],[65,311]], 1);
addAnatomyEdge(5, 7, [[86,280],[92,296],[91,310]], 1);
addAnatomyEdge(1, 8, [[200,117],[224,113],[256,119],[266,140],[274,164],[284,192]], 2);
addAnatomyEdge(8, 9, [[284,192],[290,215],[303,247],[314,280]], 1);
addAnatomyEdge(9, 10, [[314,280],[323,296],[335,311]], 1);
addAnatomyEdge(9, 11, [[314,280],[308,296],[309,310]], 1);
addAnatomyEdge(19, 12, [[212,153],[198,171],[195,197],[195,224],[200,266],[200,308]], 5);
addAnatomyEdge(12, 13, [[200,308],[174,330],[168,365],[169,399]], 2);
addAnatomyEdge(13, 14, [[169,399],[162,434],[167,478],[172,517],[172,544],[154,553]], 1);
addAnatomyEdge(13, 15, [[169,399],[178,436],[175,481],[172,517]], 1);
addAnatomyEdge(12, 16, [[200,308],[226,330],[232,365],[231,399]], 2);
addAnatomyEdge(16, 17, [[231,399],[238,434],[233,478],[228,517],[228,544],[246,553]], 1);
addAnatomyEdge(16, 18, [[231,399],[222,436],[225,481],[228,517]], 1);
const PANCREAS_VEIN = addAnatomyEdge(40, 19, [[200,216],[181,207],[170,190],[183,174],[212,153]], 1, false);
const GUT_VEIN = addAnatomyEdge('gut', 19, [[200,269],[177,251],[172,223],[170,190],[183,174],[212,153]], 1, false);
/* Depot insulin subkutan (T1D/AP): tempat injeksi/pompa — vena ke jantung.
   Posisi: sisi kanan abdomen atas, sejajar pankreas. */
TREE.nodes['subcut'] = [238, 262];
const SUBCUT_VEIN = addAnatomyEdge('subcut', 19, [[238,262],[224,243],[204,226],[194,204],[199,180],[212,153]], 1, false);
const LIVER_POS = [170,190];
const LIVER_R = 6;
const TISSUE_NODE_IDS = [4,5,8,9,13,14,15,16,17,18];
const TISSUE_NODE_POS = TISSUE_NODE_IDS.map(id => TREE.nodes[String(id)]);
// Haluskan sudut visual; pembuluh dan partikel tetap memakai titik yang sama.
function smoothAnatomyPoints(points){
  if(points.length < 3) return points;
  const result = [points[0]];
  for(let i=1;i<points.length-1;i++){
    const prev=points[i-1], p=points[i], next=points[i+1];
    const entry=[p[0]+(prev[0]-p[0])*.28,p[1]+(prev[1]-p[1])*.28];
    const exit=[p[0]+(next[0]-p[0])*.28,p[1]+(next[1]-p[1])*.28];
    result.push(entry);
    for(let j=1;j<=12;j++){
      const t=j/12, u=1-t;
      result.push([u*u*entry[0]+2*u*t*p[0]+t*t*exit[0],u*u*entry[1]+2*u*t*p[1]+t*t*exit[1]]);
    }
  }
  result.push(points[points.length-1]);
  return result;
}
TREE.edges.forEach(edge => { edge.pts = smoothAnatomyPoints(edge.pts); });
// Gambar jalur yang sama persis dengan lintasan partikel.
const anatomyVessels = document.getElementById('vessel');
TREE.edges.forEach((edge, index) => {
  const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
  path.setAttribute('d', edge.pts.map((p,i) => (i ? 'L' : 'M') + p.join(' ')).join(' '));
  path.setAttribute('fill', 'none');
  path.setAttribute('stroke-width', index >= PANCREAS_VEIN ? '1.5' : '2');
  anatomyVessels.appendChild(path);
});

const edges = TREE.edges.map(ed => {
  const pts = ed.pts;
  const cum = [0];
  for (let i=1;i<pts.length;i++){
    const dx = pts[i][0]-pts[i-1][0], dy = pts[i][1]-pts[i-1][1];
    cum.push(cum[i-1] + Math.hypot(dx,dy));
  }
  return {...ed, pts, cum, total: cum[cum.length-1] || 0.0001};
});
const children = TREE.children;
const edgeWeight = TREE.edgeWeight;

function pointOnEdge(edge, dist){
  const {pts, cum, total} = edge;
  const d = Math.max(0, Math.min(dist, total));
  let i = 1;
  while (i < cum.length && cum[i] < d) i++;
  if (i >= cum.length) i = cum.length - 1;
  const d0 = cum[i-1], d1 = cum[i];
  const t = d1 > d0 ? (d - d0)/(d1-d0) : 0;
  const p0 = pts[i-1], p1 = pts[i];
  return [ p0[0] + (p1[0]-p0[0])*t, p0[1] + (p1[1]-p0[1])*t ];
}
function pickChildEdge(nodeId){
  const opts = children[String(nodeId)];
  if (!opts || opts.length === 0) return null;
  if (opts.length === 1) return opts[0];
  let totalW = 0;
  for (const idx of opts) totalW += edgeWeight[idx];
  let r = Math.random() * totalW;
  for (const idx of opts){
    r -= edgeWeight[idx];
    if (r <= 0) return idx;
  }
  return opts[opts.length-1];
}

const svgNS = "http://www.w3.org/2000/svg";
const particleLayer = document.getElementById('particles');

/* Konstanta visual-akurasi partikel.
   Kalibrasi terukur: puncak D(t) makan besar/sedang/kecil = 2.02/1.38/0.74;
   puncak sekresi insulin ~8.7; puncak X ~0.009. */
const VIZ = {
  mealSpawnMax: 2.0,       // D(t) yang dianggap "puncak" untuk spawn gut
  basalGutFrac: 1.8,       // partikel gut tetap lahir walau basal (glukosa basal)
  insulinSecretionMax: 8.7,// sekresi puncak untuk spawn insulin
  insulinMaxCount: 28,
  glucoseMaxCount: 30,
  XFull: 0.009,            // X(t) dianggap penuh untuk tint & uptake
};

function spawnAt(p, edgeIdx){
  p.edgeIdx = edgeIdx;
  p.dist = 0;
  p.alive = true;
  p.fade = 0;               // fade-in saat lahir
}

function makeParticleSet(cls, radius, speed){
  const el = document.createElementNS(svgNS,'circle');
  el.setAttribute('r', radius);
  el.setAttribute('class', cls);
  el.style.opacity = '0';
  particleLayer.appendChild(el);
  return { el, speed, edgeIdx: null, dist: 0, alive: false, fade: 0 };
}

const bloodParticles   = Array.from({length: 34}, () => makeParticleSet('particle-blood', 2.2, 58));
const insulinParticles = Array.from({length: VIZ.insulinMaxCount}, () => makeParticleSet('particle-insulin', 2.0, 34));
const glucoseParticles = Array.from({length: VIZ.glucoseMaxCount}, () => makeParticleSet('particle-glucose', 1.6, 50));

// posisi organ, utk memarkir partikel dorman
const pancreasNode = TREE.nodes[String(TREE.insulinRoot)];
const subcutNode = TREE.nodes['subcut'];
const gutNode = TREE.nodes['gut'];
const heartNode = TREE.nodes[String(TREE.heartRoot)];
/* Asal insulin per mode: Normal = pankreas (sekresi endogen);
   T1D/AP = depot subkutan (injeksi/pompa). */
function insulinVeinForMode(){
  // T2D+PI memiliki insulin endogen dan insulin eksogen dari depot.
  if(currentMode.id === 'diabetes_t2d_pi'){
    const y = simState.y;
    const exogenous = y && y.length >= 5
      ? y[4] / (T2D_PI_ACTUATOR.tmaxI * T2D_PI_ACTUATOR.VI) : 0;
    const endogenous = Math.max(0, T2D_PARAMS.p6 *
      Math.max(simState.G - T2D_PARAMS.p5, 0) * y[6]);
    const total = exogenous + endogenous;
    return total > 0 && Math.random() < exogenous / total ? SUBCUT_VEIN : PANCREAS_VEIN;
  }
  return (currentMode.id === 'normal' || currentMode.id === 'diabetes_t2d')
    ? PANCREAS_VEIN : SUBCUT_VEIN;
}
function insulinParkNode(){
  return (currentMode.id === 'normal' || currentMode.id === 'diabetes_t2d')
    ? pancreasNode : subcutNode;
}

// scatter posisi awal darah di seluruh tree arteri (dari jantung) supaya terisi seketika
for(const p of bloodParticles){
  const e = pickChildEdge(TREE.heartRoot);
  spawnAt(p, e !== null ? e : 0);
  p.dist = Math.random() * 260;
}

// akumulator spawn-rate (partikel/detik)
const spawnAcc = { insulin: 0, gut: 0 };

function updateFlowParticles(dt){
  const dtClamped = Math.min(0.05, dt);
  const bloodColor = glucoseToColor(simState.G);

  /* --- insulin: lahir di pankreas dengan rate ∝ sekresi sesaat --- */
  const isT2d = currentMode.id === 'diabetes_t2d' || currentMode.id === 'diabetes_t2d_pi';
  const endocrineParams = isT2d ? T2D_PARAMS : CONFIG;
  const tauLocal = currentMode.id === 'diabetes_t2d_pi'
    ? simState.y[6] : (simState.tauLocal || 0);
  const above = Math.max(simState.G - endocrineParams.p5, 0);
  const secretionRate = endocrineParams.p6 * above * tauLocal;
  const secretionFrac = Math.min(secretionRate / VIZ.insulinSecretionMax, 1);
  // tingkat target jumlah partikel aktif mengikuti level plasma I absolut
  // (I/45 dari sumbu kanan grafik — insulin basal Ib juga bersirkulasi,
  // sama seperti glukosa basal), tapi dicapai lewat lahir-mati, bukan lompatan:
  const plasmaFrac = Math.max(0, Math.min(simState.I / 45, 1));
  const insulinTarget = Math.round(plasmaFrac * VIZ.insulinMaxCount);
  const insulinActive = insulinParticles.reduce((a,p)=>a + (p.alive?1:0), 0);
  let insulinRate;
  if(insulinActive < insulinTarget) insulinRate = 3 + secretionFrac * 12;   // cepat saat sekresi
  else if(insulinActive > insulinTarget) insulinRate = 0;                   // surplus: mati alami
  else insulinRate = (2 + insulinActive / 4) * Math.max(0.3, secretionFrac); // kesetimbangan
  spawnAcc.insulin += insulinRate * dtClamped;
  while(spawnAcc.insulin >= 1){
    spawnAcc.insulin -= 1;
    const p = insulinParticles.find(p=>!p.alive);
    if(p) spawnAt(p, insulinVeinForMode()); else { spawnAcc.insulin = 0; break; }
  }

  /* --- glukosa: lahir di usus dengan rate ∝ D(t) + basal --- */
  const D = currentMode.id === 'diabetes' || currentMode.id === 'diabetes_ap'
    ? t1dMealAppearance(simState.t, simState.activeMeals)
    : mealD(simState.t, simState.activeMeals);
  const mealFrac = Math.min(D / VIZ.mealSpawnMax, 1);
  const glucoseFrac = Math.max(0, Math.min(simState.G / CONFIG.zones.chartMax, 1));
  const glucoseTarget = Math.round(glucoseFrac * VIZ.glucoseMaxCount);
  const glucoseActive = glucoseParticles.reduce((a,p)=>a + (p.alive?1:0), 0);
  let gutRate;
  if(glucoseActive < glucoseTarget) gutRate = 4 + mealFrac * 16;
  else if(glucoseActive > glucoseTarget) gutRate = 0;
  else gutRate = VIZ.basalGutFrac + mealFrac * 4;
  spawnAcc.gut += gutRate * dtClamped;
  while(spawnAcc.gut >= 1){
    spawnAcc.gut -= 1;
    const p = glucoseParticles.find(p=>!p.alive);
    if(p) spawnAt(p, GUT_VEIN); else { spawnAcc.gut = 0; break; }
  }

  /* --- uptake jaringan: probabilitas terserap saat melewati node otot --- */
  const Xfrac = Math.min(simState.X / VIZ.XFull, 1);

  /* --- partikel darah: selalu bersirkulasi (sirkulasi umum) --- */
  for(const p of bloodParticles){
    advanceParticle(p, dtClamped, null);
    p.el.style.fill = bloodColor;
    if(p.fade < 1){ p.fade = 1; p.el.style.opacity = '0.9'; }
  }

  /* --- partikel insulin: clearance saat melewati area hati --- */
  for(const p of insulinParticles){
    if(!p.alive){
      park(p, insulinParkNode());
      continue;
    }
    // checkpoint ujung-edge (node) — insulin bisa dibersihkan di dekat hati:
    const cleared = advanceParticle(p, dtClamped, checkLiverNode);
    if(cleared){ killParticle(p); continue; }
    // hati juga di tengah edge vena — cek posisi kontinu (dt-independent):
    const cur = currentPos(p);
    if(cur && Math.hypot(cur[0] - LIVER_POS[0], cur[1] - LIVER_POS[1]) < LIVER_R
        && Math.random() < Math.min(0.5, 0.5 * dtClamped * 60)){
      killParticle(p); continue;
    }
    p.fade = Math.min(1, p.fade + dtClamped * 3);
    p.el.style.opacity = (0.95 * p.fade).toFixed(2);
  }

  /* --- partikel glukosa: uptake jaringan saat melewati node otot --- */
  for(const p of glucoseParticles){
    if(!p.alive){
      park(p, gutNode);
      continue;
    }
    const absorbed = advanceParticle(p, dtClamped, makeTissueCheck(Xfrac, dtClamped));
    if(absorbed){ killParticle(p); continue; }
    p.fade = Math.min(1, p.fade + dtClamped * 3);
    p.el.style.opacity = (0.85 * p.fade).toFixed(2);
  }
}

function checkLiverNode(x, y){
  // hati = area di jalur vena sekitar LIVER_POS: insulin dibersihkan (suku -n(I-Ib))
  return Math.hypot(x - LIVER_POS[0], y - LIVER_POS[1]) < LIVER_R && Math.random() < 0.35;
}

function currentPos(p){
  if(!p.alive || p.edgeIdx === null) return null;
  return pointOnEdge(edges[p.edgeIdx], p.dist);
}

function makeTissueCheck(Xfrac, dtClamped){
  const uptakeRate = 0.05 * Xfrac;   // per-detik saat berada dekat node otot
  return function(x, y){
    return uptakeRate > 0.001 && nearTissue(x, y) && Math.random() < uptakeRate * dtClamped;
  };
}

function park(p, nodePos){
  p.el.setAttribute('cx', nodePos[0].toFixed(2));
  p.el.setAttribute('cy', nodePos[1].toFixed(2));
  p.el.style.opacity = '0';
}

function killParticle(p){
  p.alive = false; p.edgeIdx = null; p.dist = 0; p.fade = 0;
  p.el.style.opacity = '0';
}

function nearTissue(x, y){
  for(const t of TISSUE_NODE_POS){
    if(Math.hypot(x - t[0], y - t[1]) < 3.2) return true;
  }
  return false;
}

/* advanceParticle: gerak sepanjang edge; callback checkpoint(x,y) boleh
   return true untuk memicu despawn (clearance/uptake). Return nilai callback
   atau false. */
function advanceParticle(p, dtClamped, checkpoint){
  if (!p.alive){
    const edgeIdx = pickChildEdge(TREE.heartRoot);
    if(edgeIdx === null) return false;
    spawnAt(p, edgeIdx);
  }
  let edge = edges[p.edgeIdx];
  p.dist += p.speed * dtClamped;

  if (p.dist > edge.total){
    // checkpoint di ujung edge (node) — insulin bisa dibersihkan, glukosa terserap
    if(checkpoint){
      const node = TREE.nodes[String(edge.child)];
      if(node && checkpoint(node[0], node[1])) return true;
    }
    const leftover = p.dist - edge.total;
    const nextEdgeIdx = pickChildEdge(edge.child);
    if (nextEdgeIdx === null){
      // leaf: darah "kembali" ke sirkulasi jantung; insulin/glukosa mati (habis terpakai/dibersihkan)
      if(p.el.classList.contains('particle-blood')){
        const ve = pickChildEdge(TREE.heartRoot);
        spawnAt(p, ve !== null ? ve : 0);
      } else {
        return true;
      }
    } else {
      p.edgeIdx = nextEdgeIdx;
      p.dist = leftover;
    }
    edge = edges[p.edgeIdx];
  }
  const [x,y] = pointOnEdge(edge, p.dist);
  p.el.setAttribute('cx', x.toFixed(2));
  p.el.setAttribute('cy', y.toFixed(2));
  return false;
}



/* =========================================================
   RENDERER: Cards
   ========================================================= */
const STATUS_ICONS = {
  'status-hipo': '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="23 18 13.5 8.5 8.5 13.5 1 6"/><polyline points="17 18 23 18 23 12"/></svg>',
  'status-normal': '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M22 11.08V12a10 10 0 1 1-5.93-9.14"/><polyline points="22 4 12 14.01 9 11.01"/></svg>',
  'status-hiper': '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="23 6 13.5 15.5 8.5 10.5 1 18"/><polyline points="17 6 23 6 23 12"/></svg>',
};
function t2dTargetStatus(G, target = T2D_CONTROL_TARGET){
  if(G < CONFIG.zones.hipoMax) return { cls: 'status-hipo', text: 'Hipoglikemia', summary: 'hipoglikemia' };
  if(Math.abs(G - target) <= 5) return { cls: 'status-normal', text: 'Target', summary: 'target' };
  if(G < target - 5) return { cls: 'status-hiper', text: 'Di bawah Target', summary: 'di bawah target' };
  if(G <= CONFIG.zones.normalMax) return { cls: 'status-hiper', text: 'Di atas Target', summary: 'di atas target' };
  return { cls: 'status-hiper', text: 'Hiperglikemia', summary: 'hiperglikemia' };
}
function updateCards(){
  document.getElementById('card-glucose').textContent = simState.G.toFixed(0);
  document.getElementById('card-insulin').textContent = simState.I.toFixed(1);
  const statusCard = document.getElementById('card-status');
  const statusText = document.getElementById('card-status-text');
  const { hipoMax, normalMax } = CONFIG.zones;
  let cls, txt;
  if(simState.G < hipoMax){ cls='status-hipo'; txt='Hipoglikemia'; }
  else if(currentMode.id === 'diabetes_t2d' || currentMode.id === 'diabetes_t2d_pi'){
    const target = currentMode.id === 'diabetes_t2d_pi' ? T2D_PI_CONFIG.target : T2D_CONTROL_TARGET;
    const status = t2dTargetStatus(simState.G, target);
    cls = status.cls;
    txt = status.text;
  }
  else if(simState.G <= normalMax){ cls='status-normal'; txt='Normal'; }
  else { cls='status-hiper'; txt='Hiperglikemia'; }

  // Update class/text/icon HANYA saat zona berubah — live region tidak spam tiap frame
  if(cls !== lastStatusClass || statusText.textContent !== txt){
    lastStatusClass = cls;
    // Pertahankan kpi-status agar kartu tetap menempati area grid yang benar.
    statusCard.className = 'kpi-card kpi-status status-card ' + cls;
    statusText.textContent = txt;
    const icon = document.getElementById('card-status-icon');
    if(icon) icon.innerHTML = STATUS_ICONS[cls] || '';
  }
}

/* =========================================================
   RENDERER: Charts (Chart.js)
   ========================================================= */
function zonesPlugin(hipoMax, normalMax, chartMax){
  return {
    id: 'zonesBg',
    beforeDatasetsDraw(chart){
      const { ctx, chartArea, scales } = chart;
      if(!chartArea) return;
      const y = scales.yG;
      const bands = [
        [0, hipoMax, 'rgba(166,43,36,0.06)'],
        [hipoMax, normalMax, 'rgba(14,122,82,0.05)'],
        [normalMax, y.max, 'rgba(143,90,16,0.07)'],
      ];
      ctx.save();
      for(const [lo,hi,color] of bands){
        const yTop = y.getPixelForValue(hi);
        const yBot = y.getPixelForValue(lo);
        ctx.fillStyle = color;
        ctx.fillRect(chartArea.left, yTop, chartArea.right-chartArea.left, yBot-yTop);
      }
      ctx.restore();
    },
    afterDatasetsDraw(chart){
      const target = currentMode.id === 'diabetes_ap' ? AP_CONFIG.targetMgDl
        : (currentMode.id === 'diabetes_t2d_pi' ? T2D_PI_CONFIG.target : null);
      if(target === null || !chart.chartArea) return;
      const {ctx, chartArea, scales} = chart;
      const y = scales.yG.getPixelForValue(target);
      ctx.save();
      ctx.strokeStyle = currentMode.id === 'diabetes_ap' ? '#0E7A52' : '#8F5A10';
      ctx.lineWidth = 1.5;
      ctx.setLineDash([7, 5]);
      ctx.beginPath(); ctx.moveTo(chartArea.left, y); ctx.lineTo(chartArea.right, y); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = ctx.strokeStyle;
      ctx.font = '600 10px sans-serif';
      ctx.textAlign = 'right';
      ctx.fillText(`Target ${target} mg/dL`, chartArea.right - 4, y - 5);
      ctx.restore();
    }
  };
}

function findNearestHistoryIndex(history, t){
  if(!history.length) return -1;
  let bestIdx = 0, bestDiff = Infinity;
  for(let i=0;i<history.length;i++){
    const diff = Math.abs(history[i].t - t);
    if(diff < bestDiff){ bestDiff = diff; bestIdx = i; }
  }
  return bestIdx;
}
function mealMarkersPlugin(){
  return {
    id: 'mealMarkers',
    afterDatasetsDraw(chart){
      const { ctx, chartArea, scales } = chart;
      if(!chartArea) return;
      const events = simState.mealEvents;
      if(!events || !events.length) return;
      const history = simState.history;
      const xScale = scales.x;
      ctx.save();
      events.forEach(ev => {
        const idx = findNearestHistoryIndex(history, ev.t);
        if(idx < 0) return;
        const x = xScale.getPixelForValue(idx);
        if(!isFinite(x) || x < chartArea.left - 1 || x > chartArea.right + 1) return;

        ctx.strokeStyle = 'rgba(143,90,16,0.85)';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([4,3]);
        ctx.beginPath();
        ctx.moveTo(x, chartArea.top);
        ctx.lineTo(x, chartArea.bottom);
        ctx.stroke();
        ctx.setLineDash([]);

        // Marker makan non-emoji: glyph garpu sederhana (path canvas)
        ctx.fillStyle = 'rgba(143,90,16,0.95)';
        ctx.strokeStyle = 'rgba(143,90,16,0.95)';
        ctx.lineWidth = 1.2;
        ctx.save();
        ctx.translate(x - 4, chartArea.top - 12);
        // tiga garis vertikal (gigi garpu) + gagang
        ctx.beginPath();
        ctx.moveTo(0, 0);  ctx.lineTo(0, 5);
        ctx.moveTo(3, 0);  ctx.lineTo(3, 5);
        ctx.moveTo(6, 0);  ctx.lineTo(6, 5);
        ctx.moveTo(3, 5);  ctx.lineTo(3, 10);
        ctx.stroke();
        ctx.restore();
      });
      ctx.restore();
    }
  };
}

let mainChart = null, xChart = null, breakdownChart = null;

function chartMaxForMode(){
  if(currentMode.id === 'normal') return CONFIG.zones.chartMax;
  if(currentMode.id === 'diabetes_t2d' || currentMode.id === 'diabetes_t2d_pi'){
    return 300;   // T2D: peak ~238 (tanpa kontrol) / ~153 (PI)
  }
  return 400;                                          // T1D/AP
}

function applyChartScale(){
  if(!mainChart) return;
  const peakG = Math.max(simState.G, ...simState.history.map(h => h.G));
  const peakI = Math.max(simState.I, ...simState.history.map(h => h.I));
  mainChart.options.scales.yG.max = Math.max(chartMaxForMode(), Math.ceil(peakG * 1.1 / 50) * 50);
  mainChart.options.scales.yI.max = Math.max(45, Math.ceil(peakI * 1.1 / 10) * 10);
}

function initCharts(){
  if(typeof Chart === 'undefined'){
    console.warn('Chart.js gagal dimuat (cek koneksi internet/CDN) — grafik dinonaktifkan, simulasi tetap berjalan.');
    return;
  }
  Chart.defaults.font.family = "'Plus Jakarta Sans', 'Inter', sans-serif";
  try{
    mainChart = new Chart(document.getElementById('mainChart'), {
      type: 'line',
      data: { labels: [], datasets: [
        { label: 'Glukosa (mg/dL) — garis penuh', data: [], borderColor:'#0284C7', backgroundColor:'transparent',
          yAxisID:'yG', tension:0.25, pointRadius:0, borderWidth:2.4 },
        { label: 'Insulin (µU/mL) — garis putus', data: [], borderColor:'#8B5CF6', backgroundColor:'transparent',
          yAxisID:'yI', tension:0.25, pointRadius:0, borderWidth:1.8, borderDash:[5,4] },
      ]},
      options: {
        animation:false, responsive:true, maintainAspectRatio:false,
        layout:{ padding:{ top:16 } },
        interaction:{ mode:'index', intersect:false },
        scales:{
          x:{ title:{display:true,text:'Waktu (menit)', font:{size:11}}, ticks:{maxTicksLimit:8, font:{size:10}} },
          yG:{ position:'left', min:0, max:CONFIG.zones.chartMax, title:{display:true,text:'Glukosa (mg/dL)', font:{size:11}}, ticks:{font:{size:10}} },
          yI:{ position:'right', min:0, max:45, grid:{drawOnChartArea:false}, title:{display:true,text:'Insulin (µU/mL)', font:{size:11}}, ticks:{font:{size:10}} },
        },
        plugins:{ legend:{ position:'bottom', labels:{ boxWidth:14, font:{size:11} } } }
      },
      plugins: [zonesPlugin(CONFIG.zones.hipoMax, CONFIG.zones.normalMax, chartMaxForMode()), mealMarkersPlugin()]
    });

    xChart = new Chart(document.getElementById('xChart'), {
      type:'line',
      data:{ labels:[], datasets:[{ label:'X(t)', data:[], borderColor:'#0284C7', pointRadius:0, borderWidth:2, tension:0.25 }]},
      options:{
        animation:false, responsive:true, maintainAspectRatio:false,
        scales:{ x:{ ticks:{maxTicksLimit:6, font:{size:10}} }, y:{ ticks:{font:{size:10}} } },
        plugins:{ legend:{ display:false } }
      }
    });

    breakdownChart = new Chart(document.getElementById('breakdownChart'), {
      type:'bar',
      data:{ labels:['-p1(G-Gb)','-X\u00b7G','D(t)'], datasets:[{
        data:[0,0,0], backgroundColor:['#475569','#A62B24','#0E7A52']
      }]},
      options:{
        animation:false, responsive:true, maintainAspectRatio:false,
        scales:{ x:{ ticks:{font:{size:10}} }, y:{ title:{display:true,text:'mg/dL per menit', font:{size:9}}, ticks:{font:{size:9}} } },
        plugins:{ legend:{ display:false } }
      }
    });
  } catch(err){
    console.error('Gagal inisialisasi Chart.js:', err);
    mainChart = null; xChart = null; breakdownChart = null;
  }
}

function pushHistory(){
  if(simState.t - simState.lastHistoryT < CONFIG.historyStepMin) return;
  simState.lastHistoryT = simState.t;
  simState.history.push({ t: simState.t, G: simState.G, I: simState.I, X: simState.X });
  const windowSteps = CONFIG.historyWindowMin / CONFIG.historyStepMin;
  if(simState.history.length > windowSteps) simState.history.shift();
}

/* Render key: chart hanya digambar ulang saat history bertambah —
   menghindari update('none') tiap frame (motion spec, performa) */
let lastRenderedKey = '';

function updateCharts(){
  if(!mainChart || !xChart) return;
  const renderKey = simState.history.length + ':' +
    (simState.history.length ? simState.history[simState.history.length-1].t.toFixed(1) : '');
  if(renderKey === lastRenderedKey) return;   // data belum berubah
  lastRenderedKey = renderKey;

  const labels = simState.history.map(h => h.t.toFixed(0));
  mainChart.data.labels = labels;
  mainChart.data.datasets[0].data = simState.history.map(h => h.G);
  mainChart.data.datasets[1].data = simState.history.map(h => h.I);
  applyChartScale();
  mainChart.update('none');

  // xChart hanya di-update saat panel math terbuka; history tetap terkumpul
  if(mathPanelOpen){
    xChart.data.labels = labels;
    xChart.data.datasets[0].data = simState.history.map(h => h.X);
    xChart.update('none');
  }

  updateChartSummary();
}

/* Summary tekstual + tabel fallback chart — aksesibilitas (audit A4) */
function updateChartSummary(){
  const hist = simState.history;
  if(!hist.length) return;
  const cur = hist[hist.length - 1];
  // Update paling banyak 1x per menit simulasi — tidak spam tiap frame
  if(Math.abs(cur.t - lastSummaryT) < 1) return;
  lastSummaryT = cur.t;

  const sum = document.getElementById('main-chart-summary');
  if(sum){
    const prev = hist.length > 1 ? hist[hist.length - 2] : cur;
    const trend = cur.G > prev.G + 0.5 ? 'naik' : (cur.G < prev.G - 0.5 ? 'turun' : 'stabil');
    const zone = (currentMode.id === 'diabetes_t2d' || currentMode.id === 'diabetes_t2d_pi')
      ? t2dTargetStatus(cur.G, currentMode.id === 'diabetes_t2d_pi' ? T2D_PI_CONFIG.target : T2D_CONTROL_TARGET).summary
      : (cur.G < CONFIG.zones.hipoMax ? 'hipoglikemia'
        : (cur.G <= CONFIG.zones.normalMax ? 'normal' : 'hiperglikemia'));
    sum.textContent = `Menit ${cur.t.toFixed(0)}: glukosa ${cur.G.toFixed(0)} mg/dL (${trend}), insulin ${cur.I.toFixed(1)} µU/mL — status ${zone}.`;
  }

  const tbody = document.querySelector('#main-chart-data tbody');
  if(tbody){
    // Tampilkan sampel data (setiap ~15 menit) agar tabel ringkas
    const rows = [];
    for(let i = 0; i < hist.length; i += 10){
      const h = hist[i];
      rows.push(`<tr><td>${h.t.toFixed(0)}</td><td>${h.G.toFixed(1)}</td><td>${h.I.toFixed(2)}</td></tr>`);
    }
    const last = hist[hist.length - 1];
    rows.push(`<tr><td>${last.t.toFixed(0)}</td><td>${last.G.toFixed(1)}</td><td>${last.I.toFixed(2)}</td></tr>`);
    tbody.innerHTML = rows.join('');
  }
}

/* =========================================================
   RENDERER: Math panel
   ========================================================= */
let mathPanelOpen = false;
let lastMathRenderT = -Infinity;
function updateMathPanel(){
  if(!mathPanelOpen) return;
  // Panel matematis tetap mengikuti simulasi tanpa menggambar ulang chart di setiap frame.
  if(Math.abs(simState.t - lastMathRenderT) < 5 && lastMathRenderT !== -Infinity) return;
  lastMathRenderT = simState.t;
  if(currentMode.id === 'diabetes_t2d_pi'){
    // T2D+PI: plant G-X-I + aktuator subkutan + controller PI (Bagian 14)
    const y = simState.y;
    if(!y || y.length < 7) return;
    const [G, X, I, S1, S2, Gs, tau] = y;
    const p3eff = T2D_PARAMS.p3 / simState.resistanceScale;
    const D = t2dPiMealD(simState.t, simState.activeMeals);
    const tP1 = -T2D_PARAMS.p1 * (G - T2D_PARAMS.Gb);
    const tX = -X * G;
    const dG = tP1 + tX + D;
    const dX = -T2D_PARAMS.p2 * X + p3eff * (I - T2D_PARAMS.Ib);
    const endogenous = T2D_PARAMS.p6 * Math.max(G - T2D_PARAMS.p5, 0) * tau;
    const extAppear = S2 / (T2D_PI_ACTUATOR.tmaxI * T2D_PI_ACTUATOR.VI);
    const dI = endogenous - T2D_PARAMS.n * (I - T2D_PARAMS.Ib) + extAppear;
    const u = t2dPiRateToPlant(simState.piRate);

    let ctrlText = '';
    if(simState.piInfo){
      const c = simState.piInfo;
      ctrlText = `
Controller PI (per 1 menit):
  e  = Gs - target = ${c.error.toFixed(2)} mg/dL
  raw = bias + Kp\u00b7e + Ki\u00b7J = ${c.unsaturated.toFixed(3)} U/h
  laju pompa = ${simState.piRate.toFixed(3)} U/h   (u = ${u.toFixed(4)} mU/kg/min)`;
    }

    document.getElementById('math-eq-text').textContent =
`T2D G-X-I + PI (Bagian 14) — t = ${simState.t.toFixed(1)} menit
p\u2083 efektif = ${p3eff.toExponential(3)} (${simState.resistanceScale.toFixed(2)}\u00d7 resistensi)

dG/dt = -p\u2081(G-Gb) - X\u00b7G + D(t)
      = ${tP1.toFixed(3)} + (${tX.toFixed(3)}) + ${D.toFixed(3)}
      = ${dG.toFixed(3)} mg/dL/min

dX/dt = -p\u2082X + p\u2083(I-Ib) = ${dX.toExponential(3)} 1/min\u00b2
dI/dt = p\u2086[G-p\u2085]\u207a\u00b7\u03c4 - n(I-Ib) + S\u2082/(tmaxI\u00b7VI)
      = ${endogenous.toFixed(3)} - ${(T2D_PARAMS.n*(I-T2D_PARAMS.Ib)).toFixed(3)} + ${extAppear.toFixed(4)}
      = ${dI.toFixed(3)} \u00b5U/mL/min
dS\u2081/dt = u - S\u2081/tmaxI,  dS\u2082/dt = S\u2081/tmaxI - S\u2082/tmaxI,  dGs/dt = (G-Gs)/\u03c4sensor${ctrlText}`;

    if(breakdownChart){
      breakdownChart.data.labels = ['-p\u2081\u00b7\u0394G','-X\u00b7G','+D(t)'];
      breakdownChart.data.datasets[0].data = [tP1, tX, D];
      breakdownChart.data.datasets[0].backgroundColor = ['#475569','#A62B24','#0E7A52'];
      breakdownChart.update('none');
    }
    return;
  }
  if(currentMode.id === 'diabetes_t2d'){
    const { G, X, I } = simState;
    const p3eff = T2D_PARAMS.p3 / simState.resistanceScale;
    const D = mealD(simState.t, simState.activeMeals);
    const tP1 = -T2D_PARAMS.p1 * (G - T2D_PARAMS.Gb);
    const tX = -X * G;
    const dG = tP1 + tX + D;
    const dX = -T2D_PARAMS.p2 * X + p3eff * (I - T2D_PARAMS.Ib);
    const dI = T2D_PARAMS.p6 * Math.max(G - T2D_PARAMS.p5, 0) * simState.tauLocal
      - T2D_PARAMS.n * (I - T2D_PARAMS.Ib);

    document.getElementById('math-eq-text').textContent =
`T2D G-X-I (struktur Mode Normal, parameter T2D) — t = ${simState.t.toFixed(1)} menit
p\u2083 efektif = ${p3eff.toExponential(3)} (${simState.resistanceScale.toFixed(2)}\u00d7 resistensi)

dG/dt = -p\u2081(G-Gb) - X\u00b7G + D(t)
      = ${tP1.toFixed(3)} + (${tX.toFixed(3)}) + ${D.toFixed(3)}
      = ${dG.toFixed(3)} mg/dL/min

dX/dt = -p\u2082X + p\u2083(I-Ib) = ${dX.toExponential(3)} 1/min\u00b2
dI/dt = p\u2086[G-p\u2085]\u207a\u00b7\u03c4 - n(I-Ib) = ${dI.toFixed(3)} \u00b5U/mL/min`;

    if(breakdownChart){
      breakdownChart.data.labels = ['-p\u2081\u00b7\u0394G','-X\u00b7G','+D(t)'];
      breakdownChart.data.datasets[0].data = [tP1, tX, D];
      breakdownChart.data.datasets[0].backgroundColor = ['#475569','#A62B24','#0E7A52'];
      breakdownChart.update('none');
    }
    return;
  }
  if(currentMode.id !== 'normal'){
    // T1D / AP: persamaan Hovorka (Bagian 10)
    const y = simState.y;
    if(!y || y.length < 8) return;
    const [Q1, Q2, S1, S2, I, x1, x2, x3] = y;
    const G = Q1 / T1D_PARAMS.VG;
    const UG = t1dMealAppearance(simState.t, simState.activeMeals);
    const u = currentMode.engine === 'ap'
      ? simState.apRate * 1000 / (60 * T1D_PARAMS.weightKg)
      : (simState.basalOn ? simState.eq.uBasal : 0);
    const termF01 = -t1dF01c(G), termX1 = -x1 * Q1, termK12 = T1D_PARAMS.k12 * Q2;
    const termFR = -t1dRenal(G), termEGP = T1D_PARAMS.EGP0 * Math.max(0, 1 - x3);
    const dQ1 = termF01 + termX1 + termK12 + termFR + UG + termEGP;

    let apText = '';
    if(currentMode.id === 'diabetes_ap' && simState.apCtrlInfo){
      const c = simState.apCtrlInfo;
      apText = `
Controller PID-IFB (per 1 menit):
  e   = SG - target = ${c.error.toFixed(2)} mg/dL
  P   = ${c.P.toFixed(4)}   Ic = ${c.Ic.toFixed(4)}   D = ${c.D.toFixed(4)} U/h
  IFB = ${c.IFB.toFixed(4)} U/h   laju pompa = ${simState.apRate.toFixed(3)} U/h
`;
    }

    document.getElementById('math-eq-text').textContent =
`Hovorka T1D (8-state) — t = ${simState.t.toFixed(1)} menit
G = Q1/VG = ${G.toFixed(3)} mmol/L (${(G * T1D_PARAMS.mmolToMgDl).toFixed(1)} mg/dL)

dQ1/dt = -F01c - x1\u00b7Q1 + k12\u00b7Q2 - FR + UG + EGP0\u00b7[1-x3]\u207a
       = ${termF01.toFixed(4)} + (${termX1.toFixed(4)}) + ${termK12.toFixed(4)} + (${termFR.toFixed(4)}) + ${UG.toFixed(4)} + ${termEGP.toFixed(4)}
       = ${dQ1.toFixed(4)} mmol/kg/min

dI/dt  = (S2/tmaxI)/VI - ke\u00b7I,  u(t) = ${u.toFixed(4)} mU/kg/min
Bolus manual: S1 ← S1 + B·1000/berat badan (mU/kg), sekali per dosis.${apText}`;

    if(breakdownChart){
      breakdownChart.data.labels = ['-F01c','-x1\u00b7Q1','+k12\u00b7Q2','-FR','+UG','+EGP'];
      const toMgDl = T1D_PARAMS.mmolToMgDl / T1D_PARAMS.VG;
      breakdownChart.data.datasets[0].data = [termF01, termX1, termK12, termFR, UG, termEGP].map(v => v * toMgDl);
      breakdownChart.data.datasets[0].backgroundColor = ['#475569','#A62B24','#0284C7','#64748B','#0E7A52','#8F5A10'];
      breakdownChart.update('none');
    }
    return;
  }
  // Mode Normal — existing
  const { G, X, I, tauLocal, t } = simState;
  const D = mealD(t, simState.activeMeals);
  const term1 = -CONFIG.p1*(G-CONFIG.Gb);
  const term2 = -X*G;
  const dG = term1+term2+D;
  const above = Math.max(G-CONFIG.p5,0);
  const secretion = CONFIG.p6*above*tauLocal;
  const clearance = -CONFIG.n*(I-CONFIG.Ib);
  const dI = secretion+clearance;

  document.getElementById('math-eq-text').textContent =
`t = ${t.toFixed(1)} menit,  \u03c4 = ${tauLocal.toFixed(1)} menit

dG/dt = -p1(G-Gb) - X\u00b7G + D(t)
      = ${term1.toFixed(3)} + ${term2.toFixed(3)} + ${D.toFixed(3)}
      = ${dG.toFixed(3)} mg/dL/menit

dI/dt = p6\u00b7[G-p5]\u207a\u00b7\u03c4 - n(I-Ib)
      = ${secretion.toFixed(3)} + (${clearance.toFixed(3)})
      = ${dI.toFixed(3)} \u00b5U/mL/menit`;

  if(breakdownChart){
    breakdownChart.data.labels = ['-p1(G-Gb)','-X\u00b7G','D(t)'];
    breakdownChart.data.datasets[0].backgroundColor = ['#475569','#A62B24','#0E7A52'];
    breakdownChart.data.datasets[0].data = [term1, term2, D];
    breakdownChart.update('none');
  }
}

/* =========================================================
   PANEL PARAMETER PER MODE
   Draft dipisahkan dari konstanta aktif. Nilai baru hanya digunakan setelah
   "Terapkan & Reset", supaya state solver lama tidak bercampur dengan model baru.
   ========================================================= */
const param = (scope, key, sym, name, min, max, step, decimals, desc, group='physiology', primary=false) =>
  ({ scope, key, sym, name, min, max, step, decimals, desc, group, primary });

const PARAMETER_SETS = {
  normal: [
    param('normal','p1','p₁','Glucose effectiveness',.01,.06,.0005,4,'Kecepatan glukosa kembali ke baseline.', 'physiology',true),
    param('normal','p3','p₃','Sensitivitas insulin',.3e-5,3e-5,.02e-5,6,'Efek insulin terhadap pengambilan glukosa.', 'physiology',true),
    param('normal','p6','p₆','Respons pankreas',.001,.008,.0001,4,'Sekresi insulin saat glukosa di atas ambang.', 'physiology',true),
    param('normal','p2','p₂','Peluruhan efek insulin',.01,.04,.0005,4,'Kecepatan efek insulin di jaringan meluruh.'),
    param('normal','n','n','Pembersihan insulin',.1,.6,.01,2,'Kecepatan insulin kembali ke baseline.'),
    param('normal','p5','p₅','Ambang sekresi pankreas',75,115,1,0,'Glukosa saat pankreas mulai meningkatkan sekresi.'),
    param('normal','Gb','Gᵇ','Glukosa basal',80,110,1,0,'Baseline glukosa puasa (mg/dL).'),
    param('normal','Ib','Iᵇ','Insulin basal',3,15,.1,1,'Baseline insulin (µU/mL).'),
    param('normal','mealTau','τ makan','Waktu penyerapan makanan',20,90,1,0,'Lebih besar berarti penyerapan makanan lebih lambat.'),
  ],
  diabetes: [
    param('t1d','targetG','Target','Target basal T1D',4.5,8,.1,1,'Target equilibrium manual (mmol/L).','physiology',true),
    param('t1d','weightKg','BB','Berat badan',40,120,1,0,'Berat badan subjek model (kg).','physiology',true),
    param('t1d','tmaxG','tmaxG','Penyerapan makanan',20,90,1,0,'Waktu penyerapan karbohidrat (menit).','physiology',true),
    param('t1d','tmaxI','tmaxI','Penyerapan insulin',30,100,1,0,'Waktu penyerapan insulin subkutan (menit).','physiology',true),
    param('t1d','SIT','SIT','Sensitivitas transport',10e-4,100e-4,1e-4,4,'Efek insulin pada transport glukosa.'),
    param('t1d','SID','SID','Sensitivitas disposal',2e-4,20e-4,.2e-4,4,'Efek insulin pada penggunaan glukosa.'),
    param('t1d','SIE','SIE','Sensitivitas EGP',100e-4,1000e-4,10e-4,4,'Efek insulin pada produksi glukosa hati.'),
    param('t1d','EGP0','EGP₀','Produksi glukosa hati',.006,.03,.0005,4,'Produksi glukosa endogen basal.'),
    param('t1d','F01','F01','Penggunaan glukosa',.004,.018,.0005,4,'Penggunaan glukosa non-insulin.'),
    param('t1d','VG','VG','Volume glukosa',.10,.25,.01,2,'Volume distribusi glukosa (L/kg).'),
    param('t1d','VI','VI','Volume insulin',.06,.20,.01,2,'Volume distribusi insulin (L/kg).'),
    param('t1d','k12','k12','Transfer antarkompartemen',.02,.12,.002,3,'Perpindahan glukosa antarkompartemen.'),
    param('t1d','ka1','ka1','Laju aksi insulin 1',.002,.02,.001,3,'Dinamika aksi insulin pertama.'),
    param('t1d','ka2','ka2','Laju aksi insulin 2',.015,.12,.005,3,'Dinamika aksi insulin kedua.'),
    param('t1d','ka3','ka3','Laju aksi insulin 3',.008,.08,.002,3,'Dinamika aksi insulin ketiga.'),
    param('t1d','ke','ke','Pembersihan insulin',.05,.25,.005,3,'Kecepatan insulin dibersihkan dari plasma.'),
    param('t1d','AG','AG','Ketersediaan makanan',.4,1,.05,2,'Fraksi karbohidrat yang masuk ke model glukosa.'),
  ],
  diabetes_ap: [],
  diabetes_t2d: [
    param('t2dUi','resistanceScale','×','Tingkat resistensi',.5,2,.05,2,'Nilai lebih tinggi menurunkan sensitivitas insulin efektif.','physiology',true),
    param('t2d','p3','p₃','Sensitivitas insulin dasar',.5e-6,8e-6,.1e-6,6,'Nilai dasar sebelum skala resistensi.','physiology',true),
    param('t2d','p6','p₆','Respons sel beta',.0005,.004,.0001,4,'Sekresi insulin endogen pada T2D.','physiology',true),
    param('t2d','Gb','Gᵇ','Glukosa basal',100,150,1,0,'Baseline glukosa puasa (mg/dL).','physiology',true),
    param('t2d','p1','p₁','Glucose effectiveness',.01,.05,.0005,4,'Penurunan glukosa tanpa insulin.'),
    param('t2d','p2','p₂','Peluruhan efek insulin',.01,.04,.0005,4,'Peluruhan efek insulin di jaringan.'),
    param('t2d','n','n','Pembersihan insulin',.1,.6,.01,2,'Kecepatan insulin kembali ke baseline.'),
    param('t2d','p5','p₅','Ambang sel beta',90,150,1,0,'Glukosa saat sekresi insulin meningkat.'),
    param('t2d','Ib','Iᵇ','Insulin basal',5,25,.5,1,'Baseline insulin (µU/mL).'),
    param('t2dUi','mealTau','τ makan','Waktu penyerapan makanan',20,90,1,0,'Lebih besar berarti makanan diserap lebih lambat.'),
  ],
  diabetes_t2d_pi: [],
};

PARAMETER_SETS.diabetes_ap = [
  ...PARAMETER_SETS.diabetes.map(d => ({...d})),
  param('ap','targetMgDl','Target','Target glukosa',90,160,1,0,'Setpoint artificial pancreas.','controller',true),
  param('ap','Ti','Ti','Waktu integral',60,600,10,0,'Respons integral controller (menit).','controller',true),
  param('ap','Td','Td','Waktu derivatif',0,150,5,0,'Respons derivatif controller (menit).','controller',true),
  param('ap','sensorTau','τ sensor','Lag sensor',1,30,1,0,'Waktu respons sensor glukosa (menit).','controller'),
  param('ap','derivativeTau','τ derivatif','Filter derivatif',1,20,1,0,'Penyaringan perubahan sensor (menit).','controller'),
  param('ap','suspendMgDl','Suspend','Batas suspend',55,90,1,0,'Pompa berhenti pada glukosa rendah (mg/dL).','controller'),
  param('ap','controllerPeriod','Δt','Periode controller',.5,5,.5,1,'Jarak antar pembaruan controller (menit).','controller'),
  param('ifb','a11','a11','Estimator IFB a11',.90,.999,.001,3,'Koefisien state estimator insulin feedback.','controller'),
  param('ifb','a21','a21','Estimator IFB a21',0,.05,.001,3,'Koefisien state estimator insulin feedback.','controller'),
  param('ifb','a31','a31','Estimator IFB a31',0,.005,.0001,4,'Koefisien state estimator insulin feedback.','controller'),
  param('ifb','a22','a22','Estimator IFB a22',.90,.999,.001,3,'Koefisien state estimator insulin feedback.','controller'),
  param('ifb','a32','a32','Estimator IFB a32',0,.05,.001,3,'Koefisien state estimator insulin feedback.','controller'),
  param('ifb','a33','a33','Estimator IFB a33',.90,.999,.001,3,'Koefisien state estimator insulin feedback.','controller'),
  param('ifb','b1','b1','Estimator IFB b1',.1,2,.01,2,'Koefisien masukan estimator insulin feedback.','controller'),
  param('ifb','b2','b2','Estimator IFB b2',0,.05,.001,3,'Koefisien masukan estimator insulin feedback.','controller'),
  param('ifb','b3','b3','Estimator IFB b3',0,.005,.0001,4,'Koefisien masukan estimator insulin feedback.','controller'),
  param('ifb','g1','g1','Bobot IFB g1',0,1,.01,2,'Bobot state estimator pertama.','controller'),
  param('ifb','g2','g2','Bobot IFB g2',0,1,.01,2,'Bobot state estimator kedua.','controller'),
  param('ifb','g3','g3','Bobot IFB g3',0,1,.01,2,'Bobot state estimator ketiga.','controller'),
];

PARAMETER_SETS.diabetes_t2d_pi = [
  ...PARAMETER_SETS.diabetes_t2d.map(d => ({...d})),
  param('pi','target','Target','Target glukosa',85,130,1,0,'Setpoint controller PI.','controller',true),
  param('pi','kp','Kp','Gain proporsional',.005,.10,.001,3,'Respons controller terhadap error glukosa.','controller',true),
  param('pi','ti','Ti','Waktu integral',60,600,10,0,'Respons integral controller (menit).','controller',true),
  param('pi','maxUph','Maks','Laju insulin maksimum',2,15,.5,1,'Batas atas sinyal kontrol (U/h).','controller'),
  param('pi','suspend','Suspend','Batas suspend',55,90,1,0,'Sinyal insulin berhenti di bawah batas ini.','controller'),
  param('pi','sensorTau','τ sensor','Lag sensor',1,30,1,0,'Waktu respons sensor glukosa (menit).','controller'),
  param('pi','period','Δt','Periode controller',.5,5,.5,1,'Jarak antar pembaruan PI (menit).','controller'),
  param('actuator','tmaxI','tmaxI','Penyerapan insulin',30,100,1,0,'Waktu penyerapan insulin subkutan (menit).','controller'),
  param('actuator','weightKg','BB','Berat badan aktuator',40,120,1,0,'Berat badan subjek model (kg).','controller'),
  param('actuator','VI','VI','Volume insulin',.06,.20,.01,2,'Volume distribusi insulin (L/kg).','controller'),
];

const parameterTarget = (scope) => ({ normal:CONFIG, t1d:T1D_PARAMS, t2d:T2D_PARAMS, t2dUi:T2D_UI_PARAMS,
  ap:AP_CONFIG, pi:T2D_PI_CONFIG, actuator:T2D_PI_ACTUATOR, ifb:IFB_PARAMS })[scope];
const parameterId = (def) => `${def.scope}.${def.key}`;
const PARAMETER_DEFAULTS = Object.fromEntries(Object.entries(PARAMETER_SETS).map(([mode, defs]) =>
  [mode, Object.fromEntries(defs.map(def => [parameterId(def), parameterTarget(def.scope)[def.key]]))]
));
const parameterDrafts = Object.fromEntries(Object.entries(PARAMETER_DEFAULTS).map(([mode, values]) => [mode, {...values}]));
let parameterDirty = false;

function recalculateDerivedParameters(){
  T1D_PARAMS.kb1 = T1D_PARAMS.ka1 * T1D_PARAMS.SIT;
  T1D_PARAMS.kb2 = T1D_PARAMS.ka2 * T1D_PARAMS.SID;
  T1D_PARAMS.kb3 = T1D_PARAMS.ka3 * T1D_PARAMS.SIE;
  IFB_PARAMS.gSum = IFB_PARAMS.g1 + IFB_PARAMS.g2 + IFB_PARAMS.g3;
}

function loadModeParameters(modeId){
  const draft = parameterDrafts[modeId];
  (PARAMETER_SETS[modeId] || []).forEach(def => { parameterTarget(def.scope)[def.key] = draft[parameterId(def)]; });
  recalculateDerivedParameters();
}

function formatParam(def, value){ return Number(value).toFixed(def.decimals); }

function buildParamSliders(){
  const container = document.getElementById('param-sliders');
  const defs = PARAMETER_SETS[currentMode.id] || [];
  const draft = parameterDrafts[currentMode.id];
  const makeRow = (def) => {
    const id = `ps-${currentMode.id}-${def.scope}-${def.key}`;
    const value = draft[parameterId(def)];
    const defaultValue = PARAMETER_DEFAULTS[currentMode.id][parameterId(def)];
    return `<div class="param-row ${value !== defaultValue ? 'is-dirty' : ''}">
      <div class="param-row-top"><label class="param-name" for="${id}"><span class="sym">${def.sym}</span> ${def.name}</label>
      <div class="param-vals"><output class="cur" id="pv-${id}" for="${id}">${formatParam(def,value)}</output> <span>(default ${formatParam(def,defaultValue)})</span></div></div>
      <input type="range" class="param-slider" id="${id}" min="${def.min}" max="${def.max}" step="${def.step}" value="${value}" aria-describedby="pd-${id}">
      <p class="param-desc" id="pd-${id}">${def.desc}</p></div>`;
  };
  const group = (name, label) => {
    const selected = defs.filter(d => d.group === name);
    if(!selected.length) return '';
    const primary = selected.filter(d => d.primary);
    const more = selected.filter(d => !d.primary);
    return `<section class="parameter-group"><h3>${label}</h3>${primary.map(makeRow).join('')}${more.length ?
      `<details class="param-more"><summary>Parameter lainnya (${more.length})</summary><div class="param-more-body">${more.map(makeRow).join('')}</div></details>` : ''}</section>`;
  };
  container.innerHTML = group('physiology','Parameter fisiologis') + group('controller','Parameter controller');
  const effective = currentMode.id === 'diabetes_t2d' || currentMode.id === 'diabetes_t2d_pi'
    ? `<p class="parameter-derived">p₃ efektif: ${(draft['t2d.p3'] / draft['t2dUi.resistanceScale']).toExponential(3)}</p>` : '';
  container.insertAdjacentHTML('beforeend', effective);
  defs.forEach(def => {
    const id = `ps-${currentMode.id}-${def.scope}-${def.key}`;
    document.getElementById(id).addEventListener('input', event => {
      draft[parameterId(def)] = Number(event.target.value);
      document.getElementById(`pv-${id}`).textContent = formatParam(def, draft[parameterId(def)]);
      parameterDirty = true;
      document.getElementById('btn-param-apply').disabled = false;
      document.getElementById('parameter-feedback').textContent = 'Perubahan belum diterapkan.';
      event.target.closest('.param-row').classList.toggle('is-dirty', draft[parameterId(def)] !== PARAMETER_DEFAULTS[currentMode.id][parameterId(def)]);
      if((def.scope === 't2d' && def.key === 'p3') || (def.scope === 't2dUi' && def.key === 'resistanceScale')) buildParamSliders();
    });
  });
  const subtitle = document.getElementById('param-subtitle');
  if(subtitle) subtitle.textContent = 'Ubah nilai, lalu terapkan untuk menghitung ulang kondisi awal.';
  document.getElementById('btn-param-apply').disabled = !parameterDirty;
}

function applyParameterDraft(){
  const targets = [CONFIG, T1D_PARAMS, T2D_PARAMS, T2D_UI_PARAMS, AP_CONFIG, T2D_PI_CONFIG, T2D_PI_ACTUATOR, IFB_PARAMS];
  const previous = targets.map(target => ({...target}));
  try {
    const draft = parameterDrafts[currentMode.id];
    if(currentMode.id === 'diabetes_t2d_pi' && draft['pi.target'] >= draft['t2d.Gb']){
      throw new Error('Target PI harus berada di bawah glukosa basal T2D agar equilibrium insulin eksternal tetap positif.');
    }
    loadModeParameters(currentMode.id);
    // Memastikan equilibrium dapat dihitung sebelum state UI diganti.
    if(currentMode.id === 'diabetes') t1dEquilibrium();
    if(currentMode.id === 'diabetes_ap') apPlantEquilibrium();
    if(currentMode.id === 'diabetes_t2d_pi') t2dPiTargetEquilibrium(T2D_PI_CONFIG.target, parameterDrafts[currentMode.id]['t2dUi.resistanceScale']);
    parameterDirty = false;
    resetSimulation();
    buildParamSliders();
    document.getElementById('parameter-feedback').textContent = 'Parameter diterapkan dan simulasi direset.';
  } catch (error) {
    targets.forEach((target, index) => Object.assign(target, previous[index]));
    recalculateDerivedParameters();
    document.getElementById('parameter-feedback').textContent = error.message || 'Parameter tidak dapat diterapkan: equilibrium tidak ditemukan.';
  }
}

function resetParamsToDefault(){
  parameterDrafts[currentMode.id] = {...PARAMETER_DEFAULTS[currentMode.id]};
  parameterDirty = true;
  buildParamSliders();
  document.getElementById('parameter-feedback').textContent = 'Nilai default siap diterapkan.';
}

/* =========================================================
   MODE TABS UI
   ========================================================= */
function buildModeTabs(){
  const container = document.getElementById('mode-tabs');
  container.innerHTML = '';
  MODES.forEach(mode => {
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'mode-tab' + (mode.id===currentMode.id?' active':'') + (!mode.enabled?' disabled':'');
    btn.innerHTML = mode.label + (mode.enabled ? '' : ' <span class="soon">Segera hadir</span>');
    btn.setAttribute('aria-pressed', mode.id===currentMode.id ? 'true' : 'false');
    if(mode.enabled){
      btn.addEventListener('click', () => {
        currentMode = mode;
        loadModeParameters(mode.id);
        parameterDirty = false;
        resetSimulation();
        buildModeTabs();
        buildParamSliders();
        updateModePanels();
      });
    } else {
      // disabled semantik: tidak fokusable via keyboard, diumumkan sebagai disabled oleh SR
      btn.disabled = true;
      btn.setAttribute('aria-disabled', 'true');
    }
    container.appendChild(btn);
  });
}

/* =========================================================
   MAIN LOOP
   ========================================================= */
let lastFrameTime = null;
let simPaused = false;
let pendingSimMinutes = 0;
/* Mode dispatch: jalankan 1 step solver sesuai mode + mirror ke simState.G/I/X */
function stepCurrentMode(dt){
  if(currentMode.id === 'normal'){
    rk4Step(simState, dt, currentMode);
    return;
  }
  if(currentMode.engine === 't2d'){
    const p3eff = T2D_PARAMS.p3 / simState.resistanceScale;
    rk4Step(simState, dt, currentMode, T2D_PARAMS, p3eff);
    return;
  }
  if(currentMode.engine === 't2d_pi'){
    // Pecah langkah pada tick controller (period 1 menit) — rate konstan sepanjang interval
    const p3eff = T2D_PARAMS.p3 / simState.resistanceScale;
    const endT = simState.t + dt;
    while(simState.t < endT - 1e-9){
      if(simState.t >= simState.nextControlT - 1e-9){
        const out = t2dPiControllerStep(simState.piController, simState.y[5], simState.piEq);
        simState.piRate = out.rate;
        simState.piSuspended = (out.rate === 0 && simState.y[5] <= T2D_PI_CONFIG.suspend);
        simState.piInfo = { error: out.error, unsaturated: out.unsaturated };
        simState.nextControlT += T2D_PI_CONFIG.period;
      }
      const h = Math.min(endT - simState.t, simState.nextControlT - simState.t);
      simState.y = t2dPiRk4Step(simState.t, simState.y, h,
        simState.activeMeals, simState.piRate, p3eff);
      simState.t += h;
    }
    // UI mirror — renderer membaca G/I/X
    simState.G = simState.y[0];
    simState.I = simState.y[2];
    simState.X = simState.y[1];
    simState.activeMeals = simState.activeMeals.filter(m =>
      (simState.t - m.startTime) < 25 * m.tauMeal);
    simState.mealEvents = simState.mealEvents.filter(m =>
      (simState.t - m.t) < CONFIG.historyWindowMin);
    return;
  }
  const cfg = {
    basalOn: simState.basalOn,
    uBasal: simState.eq.uBasal,
    boluses: simState.boluses,
    meals: simState.activeMeals,
    dt,
  };
  if(currentMode.engine === 't1d'){
    applyPendingBoluses(simState);
    simState.y = t1dRk4Step(simState.t, simState.y, dt, cfg);
    simState.t += dt;
  } else {
    // Pecah langkah pada batas tick agar laju pompa tetap sepanjang interval.
    const endT = simState.t + dt;
    while(simState.t < endT - 1e-9){
      if(simState.t >= simState.nextControlT - 1e-9){
        const out = apControllerOutput(simState.apController, simState.y[8]);
        simState.nextControlT += AP_CONFIG.controllerPeriod;
        simState.apRate = out.rate;
        simState.apSuspended = out.suspended;
        simState.apCtrlInfo = {
          error: out.error, P: out.proportional, Ic: out.integral,
          D: out.derivative, IFB: out.feedback,
        };
      }
      const h = Math.min(endT - simState.t, simState.nextControlT - simState.t);
      simState.y = apRk4Step(simState.t, simState.y, h, cfg, simState.apRate);
      simState.t += h;
    }
  }
  // Pertahankan ekor absorpsi hingga kontribusinya < 1e-9 dari dosis.
  simState.activeMeals = simState.activeMeals.filter(m => (simState.t - m.time) < 25 * T1D_PARAMS.tmaxG);
  simState.mealEvents = simState.mealEvents.filter(m => (simState.t - m.t) < CONFIG.historyWindowMin);
  // mirror ke satuan UI — renderer tidak berubah
  simState.G = simState.y[0] / T1D_PARAMS.VG * T1D_PARAMS.mmolToMgDl;
  simState.I = simState.y[4];     // mU/L ≡ µU/mL
  simState.X = simState.y[6];     // x2 — disposal effect, untuk tissue tint
}

function loop(now){
  if(lastFrameTime === null) lastFrameTime = now;
  if(simPaused){
    lastFrameTime = now;
    requestAnimationFrame(loop);
    return;
  }
  const deltaSec = Math.min((now - lastFrameTime) / 1000, 0.25);
  lastFrameTime = now;

  pendingSimMinutes += CONFIG.minutesPerSecond * deltaSec;
  const steps = Math.floor((pendingSimMinutes + 1e-9) / CONFIG.dt);
  pendingSimMinutes = Math.max(0, pendingSimMinutes - steps * CONFIG.dt);
  for(let i=0;i<steps;i++){
    stepCurrentMode(CONFIG.dt);
    pushHistory();
  }
  updateSvg(deltaSec);
  updateCards();
  updateApReadout();
  updateT2dPiReadout();
  updateCharts();
  updateMathPanel();

  requestAnimationFrame(loop);
}

/* =========================================================
   MODE PANELS (T1D manual / AP readout) + wiring
   ========================================================= */
function updateModePanels(){
  const parameters = document.querySelector('.parameter-card');
  if(parameters) parameters.hidden = false;
  document.body.classList.remove('mode-t2d');
  const pid = document.getElementById('pid-panel');
  const t1d = document.getElementById('t1d-controls');
  const ap = document.getElementById('ap-controls');
  const t2d = document.getElementById('t2d-controls');
  const pi = document.getElementById('t2d-pi-controls');
  if(!pid || !t1d || !ap) return;
  const isT1d = currentMode.id === 'diabetes';
  const isAp = currentMode.id === 'diabetes_ap';
  const isT2d = currentMode.id === 'diabetes_t2d';
  const isT2dPi = currentMode.id === 'diabetes_t2d_pi';
  pid.hidden = !(isT1d || isAp || isT2d || isT2dPi);
  t1d.hidden = !isT1d;
  ap.hidden = !isAp;
  if(t2d) t2d.hidden = !isT2d;
  if(pi) pi.hidden = !isT2dPi;
  if(isT1d && simState.eq){
    // sinkronkan toggle basal dengan state (mis. setelah reset)
    const toggle = document.getElementById('basal-toggle');
    if(toggle) toggle.checked = simState.basalOn;
  }
  if(isAp && simState.eq){
    const basalEl = document.getElementById('ap-basal');
    if(basalEl) basalEl.textContent = simState.eq.basalUph.toFixed(2) + ' U/h';
    const targetEl = document.getElementById('ap-target');
    if(targetEl) targetEl.textContent = AP_CONFIG.targetMgDl.toFixed(0) + ' mg/dL';
  }
  updateMealMenuLabels();
  const bolusFeedback = document.getElementById('bolus-feedback');
  if(bolusFeedback) bolusFeedback.textContent = '';
  const bolusInput = document.getElementById('bolus-units');
  if(bolusInput){
    bolusInput.value = '2';
    bolusInput.setCustomValidity('');
  }
  if(isT2dPi && simState.piEq){
    const biasEl = document.getElementById('pi-bias');
    if(biasEl) biasEl.textContent = simState.piEq.basalUph.toFixed(2) + ' U/h';
    const slider = document.getElementById('pi-a1-slider');
    const val = document.getElementById('pi-a1-val');
    const draftScale = parameterDrafts.diabetes_t2d_pi['t2dUi.resistanceScale'];
    if(slider) slider.value = draftScale;
    if(val) val.textContent = draftScale.toFixed(2);
    const targetEl = document.getElementById('pi-target');
    if(targetEl) targetEl.textContent = T2D_PI_CONFIG.target.toFixed(0) + ' mg/dL';
    const noteEl = document.getElementById('pi-note');
    if(noteEl) noteEl.innerHTML = `Glukosa mulai dari ${T2D_PARAMS.Gb.toFixed(0)} mg/dL tanpa kontrol; PI lalu mengarahkannya ke target ${T2D_PI_CONFIG.target.toFixed(0)} mg/dL. Makanan memakai amplitudo input model, bukan gram karbohidrat. Keluaran controller adalah <strong>sinyal kontrol model</strong>, bukan rekomendasi dosis insulin. Nilai bias tinggi merupakan konsekuensi model minimal (Bagian 14.6).`;
  }
  if(isT2d && simState.resistanceScale !== undefined){
    const slider = document.getElementById('t2d-a1-slider');
    const val = document.getElementById('t2d-a1-val');
    const draftScale = parameterDrafts.diabetes_t2d['t2dUi.resistanceScale'];
    if(slider) slider.value = draftScale;
    if(val) val.textContent = draftScale.toFixed(2);
  }
  // subtitle header mencerminkan mode
  const sub = document.querySelector('.subtitle');
  if(sub) sub.textContent = 'Model ' + (currentMode.id === 'normal'
    ? 'Bergman-Pacini \u00b7 Mode Normal (pankreas sehat)'
    : (currentMode.id === 'diabetes'
      ? 'Hovorka \u00b7 Diabetes Tipe 1 (terapi manual basal-bolus)'
      : (currentMode.id === 'diabetes_ap'
        ? 'Hovorka + PID-IFB \u00b7 Artificial Pancreas (closed-loop)'
        : (currentMode.id === 'diabetes_t2d_pi'
          ? 'G-X-I + PI \u00b7 T2D dengan kontrol insulin otomatis'
          : 'Bergman-Pacini G-X-I \u00b7 Parameter Diabetes Tipe 2'))));
}

function updateMealMenuLabels(){
  const gramsMode = currentMode.id === 'diabetes' || currentMode.id === 'diabetes_ap';
  document.querySelectorAll('#meal-options [data-size]').forEach(button => {
    const size = button.dataset.size;
    const label = size.charAt(0).toUpperCase() + size.slice(1);
    button.textContent = gramsMode
      ? `${label} · ${MEAL_GRAMS[size]} g karbohidrat`
      : `${label} · amplitudo ${CONFIG.mealAmplitude[size]}`;
  });
}

/* Bolus: tambah event { time, units } — engine memasukkan ke depot (Bagian 10.7) */
document.getElementById('btn-bolus').addEventListener('click', () => {
  if(currentMode.id !== 'diabetes') return;
  const input = document.getElementById('bolus-units');
  const feedback = document.getElementById('bolus-feedback');
  const units = Number(input.value);
  const valid = Number.isFinite(units) && units >= 0.5 && units <= 10 &&
    Math.abs(units * 2 - Math.round(units * 2)) < 1e-8;
  input.setCustomValidity(valid ? '' : 'Masukkan bolus 0,5–10 U dengan kelipatan 0,5 U.');
  if(!valid){
    if(feedback) feedback.textContent = 'Bolus belum diberikan. Masukkan nilai 0,5–10 U dengan kelipatan 0,5 U.';
    input.reportValidity();
    return;
  }
  simState.boluses.push({ time: simState.t, units });
  if(feedback) feedback.textContent = `Bolus ${units.toFixed(1)} U dimasukkan ke simulasi.`;
});

document.getElementById('bolus-units').addEventListener('input', (event) => {
  event.target.setCustomValidity('');
  const feedback = document.getElementById('bolus-feedback');
  if(feedback) feedback.textContent = '';
});

document.getElementById('basal-toggle').addEventListener('change', (e) => {
  if(currentMode.id !== 'diabetes') return;
  simState.basalOn = e.target.checked;
});

/* T2D slider resistensi — p3 efektif = p3 dasar / skala resistensi. */
document.getElementById('t2d-a1-slider').addEventListener('input', (e) => {
  if(currentMode.id !== 'diabetes_t2d') return;
  const scale = Math.max(T2D_RESISTANCE_MIN, Math.min(T2D_RESISTANCE_MAX, parseFloat(e.target.value) || 1));
  parameterDrafts.diabetes_t2d['t2dUi.resistanceScale'] = scale;
  parameterDirty = true;
  const val = document.getElementById('t2d-a1-val');
  if(val) val.textContent = scale.toFixed(2);
  buildParamSliders();
});

/* T2D+PI slider resistensi — gain PI ter-tuning utk plant nominal (disclaimer di panel) */
document.getElementById('pi-a1-slider').addEventListener('input', (e) => {
  if(currentMode.id !== 'diabetes_t2d_pi') return;
  const scale = Math.max(T2D_RESISTANCE_MIN, Math.min(T2D_RESISTANCE_MAX, parseFloat(e.target.value) || 1));
  parameterDrafts.diabetes_t2d_pi['t2dUi.resistanceScale'] = scale;
  parameterDirty = true;
  const val = document.getElementById('pi-a1-val');
  if(val) val.textContent = scale.toFixed(2);
  buildParamSliders();
});

/* T2D+PI readout live — dipanggil dari loop */
function updateT2dPiReadout(){
  if(currentMode.id !== 'diabetes_t2d_pi') return;
  const rateEl = document.getElementById('pi-rate');
  const statusEl = document.getElementById('pi-status');
  if(rateEl && Number.isFinite(simState.piRate)){
    rateEl.textContent = simState.piRate.toFixed(3) + ' U/h';
  }
  if(statusEl){
    if(simState.t === 0){
      statusEl.textContent = 'Siap — menunggu langkah PI';
      statusEl.classList.remove('suspended');
    } else if(simState.piSuspended){
      statusEl.textContent = 'Suspend (glukosa rendah)';
      statusEl.classList.add('suspended');
    } else {
      const status = t2dTargetStatus(simState.G, T2D_PI_CONFIG.target);
      statusEl.textContent = status.summary === 'target'
        ? 'Aktif — target'
        : 'Aktif — ' + status.summary;
      statusEl.classList.remove('suspended');
    }
  }
}

/* AP readout live — dipanggil dari loop via updateCards */
function updateApReadout(){
  if(currentMode.id !== 'diabetes_ap') return;
  const rateEl = document.getElementById('ap-rate');
  const statusEl = document.getElementById('ap-status');
  if(rateEl && simState.apRate !== undefined && Number.isFinite(simState.apRate)){
    rateEl.textContent = simState.apRate.toFixed(3) + ' U/h';
  }
  if(statusEl){
    if(simState.apSuspended){
      statusEl.textContent = 'Suspend (glukosa rendah)';
      statusEl.classList.add('suspended');
    } else {
      statusEl.textContent = 'Aktif';
      statusEl.classList.remove('suspended');
    }
  }
}

/* =========================================================
   INIT / EVENT WIRING
   ========================================================= */

/* ---- Meal menu: aria-expanded, Escape, click-outside, focus return ---- */
let mealMenuTimer = null;

function setMealMenuOpen(open, focusTarget){
  const menu = document.getElementById('meal-options');
  const toggle = document.getElementById('btn-meal-toggle');
  if(!menu || !toggle) return;

  if(mealMenuTimer){ clearTimeout(mealMenuTimer); mealMenuTimer = null; }

  if(open){
    menu.hidden = false;
    if(!prefersReducedMotion){
      menu.dataset.state = 'opening';
      requestAnimationFrame(() => { menu.dataset.state = 'open'; });
    }
    toggle.setAttribute('aria-expanded', 'true');
    const first = menu.querySelector('button');
    if(first && focusTarget !== false) first.focus();
  } else {
    toggle.setAttribute('aria-expanded', 'false');
    // Fokus kembali SEKARANG — penutupan visual tidak menahan keyboard
    if(focusTarget === true) toggle.focus();
    if(prefersReducedMotion){
      menu.hidden = true;
      delete menu.dataset.state;
    } else {
      // Segera keluar dari urutan fokus, lalu animasikan keluar 100ms
      menu.setAttribute('inert', '');
      menu.dataset.state = 'closing';
      mealMenuTimer = setTimeout(() => {
        menu.hidden = true;
        menu.removeAttribute('inert');
        delete menu.dataset.state;
        mealMenuTimer = null;
      }, 110);
    }
  }
}

document.getElementById('btn-meal-toggle').addEventListener('click', () => {
  const menu = document.getElementById('meal-options');
  setMealMenuOpen(menu.hidden, true);
});

document.querySelectorAll('#meal-options button').forEach(b=>{
  b.addEventListener('click', () => {
    triggerMeal(b.dataset.size);
    setMealMenuOpen(false, true); // tutup menu + fokus kembali ke toggle
  });
});

// Escape menutup menu dari dalam menu, fokus kembali ke toggle
document.getElementById('meal-options').addEventListener('keydown', (e) => {
  if(e.key === 'Escape'){
    e.stopPropagation();
    setMealMenuOpen(false, true);
  }
});

// Click di luar menutup menu
document.addEventListener('click', (e) => {
  const menu = document.getElementById('meal-options');
  const toggle = document.getElementById('btn-meal-toggle');
  if(!menu || menu.hidden) return;
  if(!menu.contains(e.target) && !toggle.contains(e.target)){
    setMealMenuOpen(false, false);
  }
});

document.getElementById('btn-reset').addEventListener('click', resetSimulation);

/* ---- Pause: aria-pressed + label dinamis + hentikan animasi SVG ---- */
const PAUSE_ICON = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><rect x="5" y="3" width="4" height="18"/><rect x="15" y="3" width="4" height="18"/></svg>';
const PLAY_ICON = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><polygon points="6 3 20 12 6 21 6 3"/></svg>';

function setSceneAnimationPaused(paused){
  const scene = document.getElementById('scene');
  if(!scene || typeof scene.pauseAnimations !== 'function') return;
  if(paused){
    scene.pauseAnimations();
  } else if(!prefersReducedMotion){
    // hanya resume jika user tidak meminta reduced-motion
    scene.unpauseAnimations();
  }
}

document.getElementById('btn-pause').addEventListener('click', () => {
  simPaused = !simPaused;
  const btn = document.getElementById('btn-pause');
  btn.innerHTML = (simPaused ? PLAY_ICON + ' Lanjutkan' : PAUSE_ICON + ' Pause');
  btn.setAttribute('aria-pressed', simPaused ? 'true' : 'false');
  btn.setAttribute('aria-label', simPaused ? 'Lanjutkan simulasi' : 'Jeda simulasi');
  btn.classList.toggle('is-paused', simPaused);
  setSceneAnimationPaused(simPaused); // jantung ikut berhenti saat simulasi pause
});
document.getElementById('btn-param-reset').addEventListener('click', resetParamsToDefault);
document.getElementById('btn-param-apply').addEventListener('click', applyParameterDraft);

/* ---- Math toggle: aria-expanded + height animation + chart sync ---- */
let mathAnimTimer = null;

document.getElementById('math-toggle').addEventListener('click', () => {
  mathPanelOpen = !mathPanelOpen;
  const body = document.getElementById('math-body');
  const toggle = document.getElementById('math-toggle');

  toggle.setAttribute('aria-expanded', mathPanelOpen ? 'true' : 'false');
  document.getElementById('math-arrow').classList.toggle('open', mathPanelOpen);

  if(mathAnimTimer){ clearTimeout(mathAnimTimer); mathAnimTimer = null; }

  if(mathPanelOpen){
    body.classList.add('open');
    body.setAttribute('aria-hidden', 'false');
    if(prefersReducedMotion){
      syncMathCharts();
    } else {
      // Buka: height 0 -> auto via animasi singkat, lalu sync chart
      body.dataset.animating = 'opening';
      body.style.height = '0px';
      body.style.opacity = '0';
      requestAnimationFrame(() => {
        body.style.height = body.scrollHeight + 'px';
        body.style.opacity = '1';
        requestAnimationFrame(() => {
          body.dataset.animating = '';
          delete body.dataset.animating;
        });
        mathAnimTimer = setTimeout(() => {
          body.style.height = 'auto';
          mathAnimTimer = null;
          syncMathCharts(); // chart diukur ulang SETELAH panel terbuka
        }, 260);
      });
    }
  } else {
    body.setAttribute('aria-hidden', 'true');
    if(prefersReducedMotion){
      body.classList.remove('open');
    } else {
      // Tutup: height current -> 0
      body.style.height = body.scrollHeight + 'px';
      body.dataset.animating = 'closing';
      requestAnimationFrame(() => {
        body.style.height = '0px';
        mathAnimTimer = setTimeout(() => {
          body.classList.remove('open');
          body.style.height = '';
          body.dataset.animating = '';
          delete body.dataset.animating;
          mathAnimTimer = null;
        }, 120);
      });
    }
  }
});

/* Sinkronkan chart di dalam math panel: data + resize */
function syncMathCharts(){
  if(!xChart || !breakdownChart) return;
  if(xChart){
    xChart.data.labels = simState.history.map(h => h.t.toFixed(0));
    xChart.data.datasets[0].data = simState.history.map(h => h.X);
    xChart.resize();
    xChart.update('none');
  }
  updateMathPanel();
}

window.addEventListener('load', () => {
  initSvgRefs();
  buildModeTabs();
  loadModeParameters(currentMode.id);
  buildParamSliders();
  initCharts();
  updateModePanels();
  requestAnimationFrame(loop);
});
