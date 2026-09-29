// Test the actual application engine, rather than a duplicate plant implementation.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const source = fs.readFileSync(path.join(__dirname, '../script.js'), 'utf8');
const context = vm.createContext({ console, assert });
vm.runInContext(source.slice(0, source.indexOf('let simState = initialState();')) + 'let simState;', context);
vm.runInContext(source.slice(source.indexOf('function mealD('), source.indexOf('/* =========================================================\n   SCENARIOS')), context);
vm.runInContext(source.slice(source.indexOf('function triggerMeal(size)'), source.indexOf('function resetSimulation()')), context);
vm.runInContext(source.slice(source.indexOf('function stepCurrentMode(dt)'), source.indexOf('function loop(now)')), context);
vm.runInContext(`
function close(actual, expected, tolerance, label){
  assert.ok(Math.abs(actual-expected) <= tolerance, label + ': ' + actual);
}
for(const dt of [.1,.5,1]){
  currentMode=ModeDiabetes; simState=initialState();
  const original=simState.y[2];
  simState.boluses.push({time:0,units:2},{time:0,units:1});
  applyPendingBoluses(simState);
  close((simState.y[2]-original)*70/1000,3,1e-12,'full dose');
  const delivered=simState.y[2]; applyPendingBoluses(simState);
  assert.equal(simState.y[2],delivered,'no duplicate bolus');
  for(let k=0;k<10;k++) stepCurrentMode(dt);
  assert.equal(simState.boluses.length,0);
}
const output=apControllerOutput, ticks=[];
apControllerOutput=(c,sg)=>{ticks.push(simState.t);return output(c,sg);};
currentMode=ModeDiabetesAP;simState=initialState();
for(let k=0;k<10;k++)stepCurrentMode(.3);
assert.deepEqual(ticks.map(t=>Math.round(t*1e6)/1e6),[0,1,2]);
apControllerOutput=output;
function scenario(mode,bolus=0,meal=true){
  currentMode=mode;simState=initialState();let min=Infinity,max=0;
  for(let k=0;k<2880;k++){
    if(k===120 && meal){simState.activeMeals.push({time:60,grams:50});
      if(bolus)simState.boluses.push({time:60,units:bolus});}
    stepCurrentMode(.5);
    assert.ok(simState.y.every(v=>Number.isFinite(v)&&v>=0));
    min=Math.min(min,simState.G);max=Math.max(max,simState.G);
  }
  return {min,max};
}
const manual=scenario(ModeDiabetes,2);
close(manual.max,182.786,.01,'T1D peak');
close(manual.min,89.558,.01,'T1D nadir');
const ap=scenario(ModeDiabetesAP);
close(ap.max,262.553,.01,'AP peak');close(ap.min,82.384,.01,'AP nadir');
for(const mode of [ModeDiabetes,ModeDiabetesAP]){
  const basal=scenario(mode,0,false);
  close(basal.max,basal.min,1e-8,'basal equilibrium');
}
function gxiMealScenario(mode,expectedPeak){
  currentMode=mode;simState=initialState();let peak=simState.G;
  for(let k=0;k<2000;k++){
    if(k===120)triggerMeal('sedang');
    stepCurrentMode(.5);
    assert.ok([simState.G,simState.X,simState.I].every(v=>Number.isFinite(v)&&v>=0));
    peak=Math.max(peak,simState.G);
  }
  close(peak,expectedPeak,.02,mode.id+' UI meal peak');
  close(simState.G,mode.id==='normal'?92:117,.01,mode.id+' return to baseline');
  return peak;
}
const normalPeak=gxiMealScenario(ModeNormal,122.54);
const t2dPeak=gxiMealScenario(ModeDiabetesT2D,157.0582);
currentMode=ModeDiabetesT2DPI;simState=initialState();
let piPeakAfterMeal=0, piNadirAfterMeal=Infinity;
const piTicks=[];
const piOutput=t2dPiControllerStep;
t2dPiControllerStep=(controller,glucose,eq)=>{
  piTicks.push(simState.t);
  return piOutput(controller,glucose,eq);
};
for(let k=0;k<4800;k++){
  if(k===720)triggerMeal('sedang');
  stepCurrentMode(.5);
  assert.ok(simState.y.slice(0,6).every(v=>Number.isFinite(v)&&v>=0));
  if(simState.t>=360){
    piPeakAfterMeal=Math.max(piPeakAfterMeal,simState.G);
    piNadirAfterMeal=Math.min(piNadirAfterMeal,simState.G);
  }
}
t2dPiControllerStep=piOutput;
assert.deepEqual(piTicks.slice(0,4),[0,1,2,3],'PI tick boundaries');
close(piPeakAfterMeal,135.9754,.03,'PI medium meal peak');
assert.ok(piNadirAfterMeal>70,'PI nadir stays above 70');
close(simState.G,98.6353,.03,'PI meal final glucose');
assert.equal(simState.activeMeals.length,0,'PI old meals are cleaned up');
assert.equal(simState.mealEvents.length,0,'PI old chart markers are cleaned up');
const piResistance=[];
for(const resistanceScale of [.5,1,2]){
  currentMode=ModeDiabetesT2DPI;simState=initialState();
  simState.resistanceScale=resistanceScale;
  simState.piEq=t2dPiTargetEquilibrium(T2D_PI_CONFIG.target,resistanceScale);
  let minG=Infinity,maxRate=0;
  for(let k=0;k<4800;k++){
    stepCurrentMode(.5);
    minG=Math.min(minG,simState.G);
    maxRate=Math.max(maxRate,simState.piRate);
  }
  piResistance.push({resistanceScale,finalG:simState.G,minG,maxRate});
  assert.ok(minG>=70,'PI resistance sweep avoids hypoglycaemia');
  close(simState.G,100,1,'PI resistance sweep reaches target');
}
const piMealResistance=[];
for(const resistanceScale of [.5,2]){
  currentMode=ModeDiabetesT2DPI;simState=initialState();
  simState.resistanceScale=resistanceScale;
  simState.piEq=t2dPiTargetEquilibrium(T2D_PI_CONFIG.target,resistanceScale);
  let minG=Infinity,maxG=0;
  for(let k=0;k<4800;k++){
    if(k===720)triggerMeal('sedang');
    stepCurrentMode(.5);
    minG=Math.min(minG,simState.G);
    if(simState.t>=360)maxG=Math.max(maxG,simState.G);
  }
  piMealResistance.push({resistanceScale,minG,maxG,finalG:simState.G});
  assert.ok(minG>=70,'PI meal resistance sweep avoids hypoglycaemia');
  close(simState.G,100,2,'PI meal resistance sweep recovers');
}
currentMode=ModeDiabetesT2DPI;simState=initialState();
let sliderChangeMinG=Infinity;
for(let k=0;k<4800;k++){
  if(k===1000){
    simState.resistanceScale=2;
    simState.piEq=t2dPiTargetEquilibrium(T2D_PI_CONFIG.target,2);
  }
  stepCurrentMode(.5);
  sliderChangeMinG=Math.min(sliderChangeMinG,simState.G);
}
assert.ok(sliderChangeMinG>=70,'live resistance change avoids hypoglycaemia');
close(simState.G,100,2,'PI recovers after live resistance change');
currentMode=ModeDiabetesT2DPI;simState=initialState();
simState.resistanceScale=2;
simState.piEq=t2dPiTargetEquilibrium(T2D_PI_CONFIG.target,2);
let reverseChangeMinG=Infinity;
for(let k=0;k<4800;k++){
  if(k===1000){
    simState.resistanceScale=.5;
    simState.piEq=t2dPiTargetEquilibrium(T2D_PI_CONFIG.target,.5);
  }
  stepCurrentMode(.5);
  reverseChangeMinG=Math.min(reverseChangeMinG,simState.G);
}
assert.ok(reverseChangeMinG>=70,'reverse live resistance change avoids hypoglycaemia');
close(simState.G,100,2,'PI recovers after reverse live resistance change');
console.log('PASS: T1D/AP regression; Normal/T2D UI meal inputs; PI tick, meal response, positivity, and event cleanup.',
  {manual,ap,normalPeak,t2dPeak,piPeakAfterMeal,piResistance,piMealResistance,
    sliderChange:{minG:sliderChangeMinG,reverseMinG:reverseChangeMinG,finalG:simState.G}});
`, context);
