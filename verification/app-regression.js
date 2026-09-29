// Test the actual application engine, rather than a duplicate plant implementation.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const source = fs.readFileSync(path.join(__dirname, '../script.js'), 'utf8');
const context = vm.createContext({ console, assert });
vm.runInContext(source.slice(0, source.indexOf('let simState = initialState();')) + 'let simState;', context);
vm.runInContext(source.slice(source.indexOf('function mealD('), source.indexOf('/* =========================================================\n   SCENARIOS')), context);
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
currentMode=ModeNormal;simState=initialState();
for(let k=0;k<2880;k++)stepCurrentMode(.5);
close(simState.G,92,1e-10,'Normal basal');
currentMode=ModeDiabetesT2D;simState=initialState();
let t2dPeak=simState.G;
for(let k=0;k<2000;k++){
  if(k===120)simState.activeMeals.push({startTime:60,A:150,tauMeal:40});
  stepCurrentMode(.5);
  assert.ok([simState.G,simState.X,simState.I].every(v=>Number.isFinite(v)&&v>=0));
  t2dPeak=Math.max(t2dPeak,simState.G);
}
close(t2dPeak,157.0582,.01,'T2D G-X-I peak');
close(simState.G,117,.01,'T2D return to baseline');
console.log('PASS: full bolus, exactly-once delivery, PID tick boundaries, T1D/AP reference responses, Normal and T2D G-X-I baselines.', {manual,ap,t2dPeak});
`, context);
