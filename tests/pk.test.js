import {test} from 'node:test';
import assert from 'node:assert/strict';
import {defaults,derive,single,concentration,theoreticalAUC,integrate,peak,simulate,results,validate,csv} from '../pk-model.js';
const close=(actual,expected,rel=1e-7)=>assert.ok(Math.abs(actual-expected)<=rel*Math.max(1,Math.abs(expected)),`${actual} ≠ ${expected}`);
const p={...defaults};
test('one-compartment IV: initial concentration, half-life, AUC',()=>{
 close(single(p,0),10000);close(derive(p).halfLife,6.9314718056);close(theoreticalAUC(p),100000);
 close(integrate(t=>single(p,t),0,400),100000);
});
test('one-compartment PO: origin, analytic Tmax, AUC, equal/near rates',()=>{
 const a={...p,route:'po'};close(single(a,0),0);
 close(peak(a,0,72).t,Math.log(a.ka/0.1)/(a.ka-0.1),1e-6);
 close(integrate(t=>single(a,t),0,500),theoreticalAUC(a));
 for(const ka of [0.1,0.1*(1+1e-12),0.1*(1-1e-12)])close(single({...a,ka},10),a.f*a.dose/a.v*0.1*10*Math.exp(-1)*1000);
});
test('two-compartment IV: eigenvalues, initial condition, AUC, terminal',()=>{
 const a={...p,model:'two'},d=derive(a);assert.ok(d.alpha>d.beta&&d.beta>0);
 close(single(a,0),a.dose/a.vc*1000);close(d.halfLife,Math.LN2/d.beta);
 close(integrate(t=>single(a,t),0,1000),theoreticalAUC(a));
});
// Independent RK4 integration of gut, central and peripheral amounts.
function rk4(p,t,steps=20000){
 const h=t/steps;let x=[p.f*p.dose,0,0];
 const fn=([g,c,v])=>[-p.ka*g,p.ka*g-(p.cl+p.q)/p.vc*c+p.q/p.vp*v,p.q/p.vc*c-p.q/p.vp*v];
 for(let i=0;i<steps;i++){const a=fn(x),b=fn(x.map((v,j)=>v+h*a[j]/2)),c=fn(x.map((v,j)=>v+h*b[j]/2)),d=fn(x.map((v,j)=>v+h*c[j]));x=x.map((v,j)=>v+h*(a[j]+2*b[j]+2*c[j]+d[j])/6);}
 return x[1]/p.vc*1000;
}
test('two-compartment PO: independent ODE validation, repeated roots, flip-flop, AUC',()=>{
 const a={...p,model:'two',route:'po'},d=derive(a);
 for(const ka of [1,d.alpha,d.beta,d.alpha*(1+1e-11),d.beta*(1-1e-11),0.01]){
  const b={...a,ka};close(single(b,0),0);for(const t of [0.1,2,20])close(single(b,t),rk4(b,t),2e-7);
  let area=0;for(let t=0;t<5000;t+=5)area+=integrate(x=>single(b,x),t,t+5);
  close(area,theoreticalAUC(b),1e-6);
 }
 const b={...a,ka:0.01};close(derive(b).apparentHalfLife,Math.LN2/0.01);
 close(Math.log(single(b,1000)/single(b,1010))/10,0.01,1e-6);
});
for(const model of ['one','two'])for(const route of ['iv','po'])for(const dosing of ['single','multiple'])test(`${model}/${route}/${dosing}: simulation, results and CSV`,()=>{
 const a={...p,model,route,dosing,end:72},r=results(a),rows=simulate(a);
 assert.ok(rows.length>600);assert.ok(rows.every(x=>Number.isFinite(x.c)&&x.c>=0));
 close(r.theoreticalAUC,theoreticalAUC(a));assert.ok(r.auc<=r.theoreticalAUC*(1+1e-8));
 assert.ok(csv(a,rows).includes('Time (h),Concentration (ng/mL),Event state'));
 if(dosing==='multiple'){
  for(const t of [12,24,37,48]){
   close(concentration(a,t),Array.from({length:5},(_,i)=>t>=12*i?single({...a,dosing:'single'},t-12*i):0).reduce((x,y)=>x+y));
  }
  const pre=concentration(a,24,'pre'),post=concentration(a,24,'post');
  close(post-pre,route==='iv'?a.dose/(model==='one'?a.v:a.vc)*1000:0);
  assert.ok(r.first&&r.last);assert.equal(results({...a,end:49}).last,null);
  if(route==='iv')assert.deepEqual(rows.filter(x=>x.t===24).map(x=>x.state),['pre','post']);
 }
});
test('validation rejects invalid inputs and accepts zero dose/F',()=>{
 for(const [key,value] of [['dose',-1],['cl',0],['v',-1],['end',0],['dose',NaN]])assert.ok(validate({...p,[key]:value}).length);
 for(const patch of [{model:'two',vc:0},{model:'two',vp:0},{model:'two',q:0},{route:'po',ka:0},{route:'po',f:1.1},{route:'po',f:-0.1},{dosing:'multiple',tau:0},{dosing:'multiple',n:1.5},{dosing:'multiple',n:0},{dosing:'multiple',end:48}])assert.ok(validate({...p,...patch}).length);
 assert.deepEqual(validate({...p,dose:0}),[]);assert.deepEqual(validate({...p,route:'po',f:0}),[]);
});
test('single-dose peak is independent of a truncated plot; IV interval excludes next bolus',()=>{
 const po={...p,route:'po',end:0.1};const r=results(po);
 close(r.global.t,Math.log(10)/0.9,1e-6);assert.ok(r.global.t>po.end);
 const iv={...p,dosing:'multiple',n:2,end:25};const v=results(iv);
 close(v.first.c,10000);close(v.first.t,0);close(v.last.c,concentration(iv,12));
 close(v.auc,1000*iv.dose/iv.cl*((1-Math.exp(-0.1*25))+(1-Math.exp(-0.1*13))));
});
test('zero dose and zero F produce zero exposure, including log-export source values',()=>{
 for(const a of [{...p,dose:0},{...p,model:'two',route:'po',f:0}]){
  assert.ok(simulate(a).every(r=>r.c===0));const r=results(a);close(r.global.c,0);close(r.auc,0);close(r.theoreticalAUC,0);
 }
});
