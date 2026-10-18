import {test} from 'node:test';
import assert from 'node:assert/strict';
import {defaults,single,simulate,results,derive} from '../pk-model.js';
import {parseTimes,sampleProfile,reportedSamples,analyzeNCA,trapezoid,regress,compareNCA,ncaCSV} from '../nca.js';
const close=(actual,expected,tolerance=1e-8)=>assert.ok(Math.abs(actual-expected)<=tolerance*Math.max(1,Math.abs(expected)),`${actual} ≠ ${expected}`);
const p={...defaults};
const record=(t,c)=>({t,c,status:c===null?'BLQ':'Quantifiable',event:'sample'});
const run=(config,times,lloq=1,options={},zero='post')=>{const samples=sampleProfile(config,times,lloq,zero);return {samples,nca:analyzeNCA(reportedSamples(samples),{route:config.route,...options})};};
test('ideal one-compartment IV: exact log slope, log-down AUC, extrapolation and units',()=>{
 const {nca}=run(p,[0,2,4,8,12,24,48]);
 close(nca.lambda,0.1);close(nca.halfLife,Math.log(2)/0.1);close(nca.fit.r2,1);close(nca.fit.adjustedR2,1);
 close(nca.aucLast,100000*(1-Math.exp(-4.8)));close(nca.aucInf,100000);close(nca.extrapolated,Math.exp(-4.8)*100);
 assert.deepEqual(nca.selectedTimes,[0,2,4,8,12,24,48]);close(nca.fit.span,48/nca.halfLife);
 close(nca.cmax,10000);close(nca.tmax,0);
});
test('PO missing Tmax reduces observed Cmax and shifts observed Tmax',()=>{
 const a={...p,route:'po'},theory=results(a),{nca}=run(a,[0,1,8,12,24,48]);
 assert.ok(nca.cmax<theory.global.c);assert.ok(Math.abs(nca.tmax-theory.global.t)>0.5);
 const cmp=compareNCA(nca,theory,'po');assert.ok(cmp[0].bias<0);
 assert.ok(nca.selectedTimes.every(t=>t>nca.tmax));
});
test('early two-compartment samples can have excellent R² without a terminal half-life',()=>{
 const a={...p,model:'two'},theory=results(a),{nca}=run(a,[0,0.01,0.02,0.03,0.04]);
 assert.ok(nca.fit.r2>0.99);assert.ok(Math.abs(nca.halfLife-theory.halfLife)/theory.halfLife>0.5);
 assert.ok(nca.fit.span<2);assert.ok(nca.warnings.some(w=>w.includes('terminal span')));
 assert.ok(nca.warnings.some(w=>w.includes('high R² does not')));
});
test('LLOQ removes terminal points and equality is quantifiable',()=>{
 const equal=single(p,2),samples=sampleProfile(p,[0,2,4,8],equal,'post');
 assert.equal(samples[1].status,'Quantifiable');assert.equal(samples[2].c,null);
 const nca=analyzeNCA(reportedSamples(samples));assert.equal(nca.lambda,null);assert.match(nca.reasons.lambda,/three/);
 assert.equal(nca.aucInf,null);assert.ok(nca.aucLast>0);
 assert.equal(analyzeNCA(reportedSamples(sampleProfile(p,[1,2,4],1e10))).cmax,null);
});
test('BLQ true values are inaccessible to analysis; BLQ non-null input is rejected',()=>{
 const samples=sampleProfile(p,[0,2,4,8,12,24],1000,'post'),reported=reportedSamples(samples);
 assert.ok(reported.every(s=>!('trueC' in s)));
 const before=analyzeNCA(reported);samples.forEach(s=>{s.trueC=9e99;});
 assert.deepEqual(analyzeNCA(reportedSamples(samples)),before);
 assert.throws(()=>analyzeNCA([record(0,10),{t:1,c:0.1,status:'BLQ',event:'sample'}]),/null/);
});
test('intermediate BLQ blocks AUC by default; bridging is explicit and traceable',()=>{
 const data=[record(0,100),record(1,null),record(2,50),record(3,40),record(4,30),record(8,null)];
 const strict=analyzeNCA(data);assert.equal(strict.aucLast,null);assert.equal(strict.aucInf,null);
 assert.deepEqual(strict.bridges,[{from:0,to:2,blqTimes:[1]}]);assert.match(strict.reasons.aucLast,/Intermediate BLQ/);
 const bridge=analyzeNCA(data,{blqPolicy:'bridge'});
 close(bridge.aucLast,trapezoid(data[0],data[2])+trapezoid(data[2],data[3])+trapezoid(data[3],data[4]));
 assert.ok(bridge.aucInf>bridge.aucLast);assert.ok(bridge.warnings.some(w=>w.includes('Explicit bridging')));
 assert.ok(bridge.selectedTimes.every(t=>t>1));assert.equal(bridge.lastTime,4);
});
test('manual selection uses only chosen eligible times and rejects BLQ/gap/too few/nonnegative slopes',()=>{
 const data=[record(0,100),record(1,70),record(2,50),record(4,25),record(8,10)];
 const nca=analyzeNCA(data,{method:'manual',selected:[1,2,8]});
 assert.deepEqual(nca.selectedTimes,[1,2,8]);close(nca.lambda,regress([data[1],data[2],data[4]]).lambda);
 assert.equal(analyzeNCA(data,{method:'manual',selected:[1,2]}).lambda,null);
 assert.equal(analyzeNCA(data,{method:'manual',selected:[1,2,2]}).lambda,null);
 const gap=[record(0,100),record(1,null),record(2,50),record(4,25),record(8,10)];
 assert.equal(analyzeNCA(gap,{method:'manual',selected:[0,2,4]}).lambda,null);
 assert.match(regress([record(0,10),record(1,20),record(2,30)]).reason,/negative slope/);
 assert.match(regress([record(0,10),record(1,10),record(2,10)]).reason,/negative slope/);
});
test('IV predose BLQ permits Phoenix-style back extrapolation to C0 and complete AUC',()=>{
 const times=[0,1,2,4,8,12,24];
 const samples=sampleProfile(p,times,1);
 assert.equal(samples[0].event,'pre');assert.equal(samples[0].c,null);
 const nca=analyzeNCA(reportedSamples(samples));
 close(nca.c0,10000);
 assert.equal(nca.c0Method,'log-linear-first-two');
 close(nca.aucBack,100000*(1-Math.exp(-0.1)));
 close(nca.aucLast,100000*(1-Math.exp(-2.4)));
 close(nca.aucInf,100000);
 close(nca.aucBackPct,nca.aucBack/nca.aucInf*100);
 close(nca.cmax,single(p,1));close(nca.tmax,1);
 assert.ok(!nca.selectedTimes.includes(0),'synthetic C0 must not enter lambda-z');
 const withoutTimeZero=run(p,[1,2,4,8,12,24]).nca;
 close(withoutTimeZero.aucLast,nca.aucLast);
 const post=analyzeNCA(reportedSamples(sampleProfile(p,times,1,'post')));
 assert.equal(post.c0Method,'observed-postdose');close(post.c0,10000);
 close(post.aucBack,0);close(post.aucLast,nca.aucLast);
});
test('IV C0 fallback uses first quantifiable value when early slope is nonnegative',()=>{
 const data=[{t:0,c:null,status:'BLQ',event:'pre'},record(1,10),record(2,11),record(3,8),record(4,6)];
 const nca=analyzeNCA(data);
 assert.equal(nca.c0Method,'fallback-first');close(nca.c0,10);
 close(nca.aucBack,10);assert.ok(nca.aucLast>nca.aucBack);
});
test('IV without observed 0+ and only one quantifiable postdose point cannot back-extrapolate C0',()=>{
 const data=[{t:0,c:null,status:'BLQ',event:'pre'},record(2,10)];
 const nca=analyzeNCA(data);
 assert.equal(nca.c0,null);assert.equal(nca.aucLast,null);assert.equal(nca.aucInf,null);
 assert.match(nca.reasons.c0,/two quantifiable/);
});
test('PO leading BLQ has explicit zero rule, trailing BLQ stops AUC',()=>{
 const data=[record(0,null),record(1,null),record(2,50),record(4,30),record(8,10),record(12,3),record(24,null)];
 const nca=analyzeNCA(data,{route:'po'});
 close(nca.aucLast,25+trapezoid(data[2],data[3])+trapezoid(data[3],data[4])+trapezoid(data[4],data[5]));
 assert.ok(nca.warnings.some(w=>w.includes('assigned zero')));assert.equal(nca.lastTime,12);
 assert.ok(nca.selectedTimes.every(t=>t>2));
});
test('flip-flop comparison uses apparent PO half-life; zero benchmarks have no relative bias',()=>{
 const a={...p,route:'po',ka:0.01,end:1000},theory=results(a),{nca}=run(a,[0,50,100,200,300,400,600],0.01,{method:'manual',selected:[200,300,400,600]});
 const cmp=compareNCA(nca,theory,'po');close(cmp[3].theory,Math.LN2/0.01);assert.ok(theory.apparentHalfLife>theory.halfLife);
 close(nca.halfLife,Math.LN2/0.01,1e-6);
 const iv=run(p,[0,1,2,4]).nca;assert.equal(compareNCA(iv,results(p),'iv')[1].bias,null);
 const empty=run({...p,dose:0},[0,1,2,4]).nca;assert.equal(compareNCA(empty,results({...p,dose:0}),'iv')[0].bias,null);
});
test('sampling parsing and validation',()=>{
 assert.deepEqual(parseTimes('2, 0\n1\n\n.5, 8.3e-2',24),[0,0.083,0.5,1,2]);
 for(const input of ['-1,2','1,1.0','NaN','Infinity','','1 h','1,25','0x10'])assert.throws(()=>parseTimes(input,24));
 for(const lloq of [0,-1,NaN,Infinity])assert.throws(()=>sampleProfile(p,[0,1],lloq));
 assert.throws(()=>sampleProfile({...p,dosing:'multiple'},[0,1],1),/single/);
 assert.throws(()=>sampleProfile(p,[1,1],1));assert.throws(()=>sampleProfile(p,[73],1));
 assert.throws(()=>analyzeNCA([],{thresholds:{adjustedR2:NaN,span:2,extrapolated:20}}));
});
test('linear-up/log-down trapezoids: rising, equal, zero and close concentrations',()=>{
 close(trapezoid({t:0,c:0},{t:2,c:100}),100);close(trapezoid({t:0,c:100},{t:2,c:100}),200);
 close(trapezoid({t:0,c:100},{t:2,c:50}),100/Math.log(2));close(trapezoid({t:0,c:100},{t:2,c:0}),100);
 close(trapezoid({t:0,c:1},{t:2,c:1-1e-13}),2,1e-12);
});
test('diagnostic thresholds are configurable, not a model substitute',()=>{
 const {nca}=run(p,[0,1,2,3],1,{thresholds:{adjustedR2:0.99,span:10,extrapolated:1}});
 assert.ok(nca.warnings.some(w=>w.includes('span < 10')));assert.ok(nca.warnings.some(w=>w.includes('> 1%')));
 const invalid=run(p,[0,1],1).nca;assert.equal(invalid.lambda,null);assert.equal(invalid.halfLife,null);assert.equal(invalid.aucInf,null);
});
test('CSV preserves true and reported values separately, regression selection, reasons and results',()=>{
 const {samples,nca}=run(p,[0,1,2,4,8,24],2000),comparison=compareNCA(nca,results(p),'iv');
 const text=ncaCSV(samples,nca,comparison,{lloq_ng_ml:2000,description:'A, "quoted" field'});
 assert.match(text,/Underlying theoretical concentration/);assert.match(text,/"BLQ"/);assert.match(text,/lloq_ng_ml/);assert.match(text,/lambda_z_times_h/);assert.match(text,/auc0_inf_ng_h_ml/);assert.match(text,/"A, ""quoted"" field"/);
 assert.equal(text.split('\r\n').length,samples.length+1);
});
test('NCA does not mutate PK inputs, concentration-time rows, theory or any of eight model outputs',()=>{
 for(const model of ['one','two'])for(const route of ['iv','po'])for(const dosing of ['single','multiple']){
  const a={...p,model,route,dosing},snapshot=JSON.stringify(a),before=simulate(a),theory=results(a);
  if(dosing==='single'){const {nca}=run(a,[0,1,2,4,8,12,24]);compareNCA(nca,theory,route);}
  else assert.throws(()=>sampleProfile(a,[0,1],1));
  assert.equal(JSON.stringify(a),snapshot);assert.deepEqual(simulate(a),before);assert.deepEqual(results(a),theory);
 }
});
