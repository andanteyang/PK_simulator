// All public concentrations are ng/mL; internal disposition coefficients are mg/L.
export const defaults = {model:'one', route:'iv', dosing:'single', dose:10, cl:0.1, v:1, vc:0.5, vp:1, q:0.3, f:0.8, ka:1, tau:12, n:5, end:72};
export function validate(p) {
  const errors=[];
  for(const [key,values] of Object.entries({model:['one','two'],route:['iv','po'],dosing:['single','multiple']})) if(!values.includes(p[key])) errors.push(`Invalid ${key}.`);
  const positive=['cl','end',...(p.model==='one'?['v']:['vc','vp','q']),...(p.route==='po'?['ka']:[]),...(p.dosing==='multiple'?['tau']:[])];
  for(const key of positive) if(!Number.isFinite(p[key])||p[key]<=0) errors.push(`${key} must be positive and finite.`);
  if(!Number.isFinite(p.dose)||p.dose<0) errors.push('Dose must be nonnegative and finite.');
  if(p.route==='po'&&(!Number.isFinite(p.f)||p.f<0||p.f>1)) errors.push('F must be between 0 and 1.');
  if(p.dosing==='multiple') {
    if(!Number.isInteger(p.n)||p.n<1||p.n>100) errors.push('Number of doses must be an integer from 1 to 100.');
    if(p.end<=(p.n-1)*p.tau) errors.push('Simulation end must be later than the last dose.');
  }
  return errors;
}
export function derive(p) {
  if(p.model==='one') { const kel=p.cl/p.v; return {kel,terminal:kel,halfLife:Math.LN2/kel,apparentHalfLife:Math.LN2/(p.route==='po'?Math.min(kel,p.ka):kel)}; }
  const k10=p.cl/p.vc,k12=p.q/p.vc,k21=p.q/p.vp;
  const sum=k10+k12+k21;
  const delta=Math.hypot(k10+k12-k21,2*Math.sqrt(k12)*Math.sqrt(k21));
  const alpha=(sum+delta)/2,beta=k10*k21/alpha;
  const a=(alpha-k21)/(alpha-beta),b=(k21-beta)/(alpha-beta);
  return {k10,k12,k21,alpha,beta,a,b,terminal:beta,halfLife:Math.LN2/beta,apparentHalfLife:Math.LN2/(p.route==='po'?Math.min(beta,p.ka):beta)};
}
// Stable convolution of exp(-r*t) with ka*exp(-ka*t), including equal rates.
export function absorption(r,ka,t) {
  const d=Math.abs(ka-r);
  return ka*Math.exp(-Math.min(r,ka)*t)*(d===0?t:-Math.expm1(-d*t)/d);
}
export function single(p,t,d=derive(p)) {
  if(t<0) return 0;
  const kernel=r=>p.route==='iv'?Math.exp(-r*t):p.f*absorption(r,p.ka,t);
  return 1000*p.dose*(p.model==='one'?kernel(d.kel)/p.v:(d.a*kernel(d.alpha)+d.b*kernel(d.beta))/p.vc);
}
export function doseTimes(p) {return Array.from({length:p.dosing==='single'?1:p.n},(_,i)=>i*p.tau);}
export function concentration(p,t,state='post',d=derive(p)) {
  return doseTimes(p).reduce((sum,td)=>sum+(td<t||(td===t&&state!=='pre')?single(p,t-td,d):0),0);
}
export function theoreticalAUC(p) {return 1000*p.dose*(p.dosing==='single'?1:p.n)*(p.route==='po'?p.f:1)/p.cl;}
export function integrate(fn,a,b,tol=1e-7,depth=22) {
  if(b<=a) return 0;
  const sim=(a,b,fa,fm,fb)=>(b-a)*(fa+4*fm+fb)/6;
  function recurse(a,b,fa,fm,fb,s,e,n) {
    const m=(a+b)/2,l=fn((a+m)/2),r=fn((m+b)/2);
    const sl=sim(a,m,fa,l,fm),sr=sim(m,b,fm,r,fb),diff=sl+sr-s;
    if(n===0||Math.abs(diff)<=15*e) return sl+sr+diff/15;
    return recurse(a,m,fa,l,fm,sl,e/2,n-1)+recurse(m,b,fm,r,fb,sr,e/2,n-1);
  }
  const fa=fn(a),fm=fn((a+b)/2),fb=fn(b);
  return recurse(a,b,fa,fm,fb,sim(a,b,fa,fm,fb),tol,depth);
}
function points(p,end=p.end) {
  const d=derive(p), times=new Set([0,end]);
  for(let i=1;i<601;i++) times.add(end*i/600);
  const fast=Math.max(d.alpha||d.kel,p.route==='po'?p.ka:0);
  for(const td of doseTimes(p)) {
    if(td>end) continue;
    times.add(td);
    for(let j=0;j<=180;j++) {const age=Math.exp(Math.log(1e-8/fast)+(Math.log(Math.max(end-td,1e-8/fast))-Math.log(1e-8/fast))*j/180);if(td+age<=end)times.add(td+age);}
  }
  return [...times].sort((a,b)=>a-b);
}
export function simulate(p) {
  const errors=validate(p);if(errors.length)throw new Error(errors.join(' '));
  const d=derive(p),events=new Set(doseTimes(p)),rows=[];
  for(const t of points(p)) {
    if(p.route==='iv'&&events.has(t)) rows.push({t,c:concentration(p,t,'pre',d),state:'pre'});
    rows.push({t,c:concentration(p,t,'post',d),state:events.has(t)?'post':'sample'});
  }
  if(rows.some(r=>!Number.isFinite(r.c)))throw new Error('Parameters exceed numerical range. Choose less extreme values.');
  return rows;
}
export function peak(p,start,end,excludeEndDose=false) {
  const d=derive(p),fn=t=>concentration(p,t,excludeEndDose&&t===end?'pre':'post',d);
  const grid=points(p,end).filter(t=>t>=start);grid.unshift(start);
  // Evaluate each smooth interval separately, retaining both sides of IV jumps.
  let best={c:fn(start),t:start};
  const take=(t,c=fn(t))=>{if(c>best.c)best={t,c};};
  for(let i=0;i<grid.length;i++) {
    take(grid[i]);
    if(p.route==='iv'&&doseTimes(p).includes(grid[i]))take(grid[i],concentration(p,grid[i],'pre',d));
    if(i>0&&i<grid.length-1&&fn(grid[i])>=fn(grid[i-1])&&fn(grid[i])>=fn(grid[i+1])&&(fn(grid[i])>fn(grid[i-1])||fn(grid[i])>fn(grid[i+1]))) {
      let a=grid[i-1],b=grid[i+1];
      // PO curves are smooth between doses; each bracket is tiny near fast phases.
      for(let j=0;j<70;j++){const x=a+(b-a)/3,y=b-(b-a)/3;if(fn(x)<fn(y))a=x;else b=y;}
      take((a+b)/2);
    }
  }
  return best;
}
export function results(p) {
  const d=derive(p),global=peak(p,0,p.dosing==='single'?Math.max(p.end,40/Math.min(d.terminal,p.route==='po'?p.ka:d.terminal)):p.end),times=doseTimes(p);
  const interval=(start)=>start+p.tau<=p.end?(()=>{const k=peak(p,start,start+p.tau,true);return {c:k.c,t:k.t-start};})():null;
  // Integrate each dose on a logarithmic grid so narrow absorption/distribution phases cannot be missed.
  let auc=0;
  for(const td of times){const end=p.end-td;const grid=points({...p,dosing:'single'},end);for(let i=1;i<grid.length;i++)auc+=integrate(t=>single(p,t,d),grid[i-1],grid[i],1e-6);}
  return {...d,global,first:p.dosing==='multiple'?interval(0):null,last:p.dosing==='multiple'?interval(times.at(-1)):null,auc,theoreticalAUC:theoreticalAUC(p)};
}
export function csv(p,rows) {
  const keys=['model','route','dosing','dose','cl',...(p.model==='one'?['v']:['vc','vp','q']),...(p.route==='po'?['f','ka']:[]),...(p.dosing==='multiple'?['tau','n']:[]),'end'];
  return ['Time (h),Concentration (ng/mL),Event state,'+keys.join(','),...rows.map(r=>[r.t,r.c,r.state,...keys.map(k=>p[k])].join(','))].join('\r\n');
}
