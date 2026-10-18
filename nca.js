// Sampling adapter reuses the unchanged PK engine. Analysis receives reported data only.
import {single,validate} from './pk-model.js';
export const defaultTimes='0, 0.083, 0.25, 0.5, 1, 2, 4, 8, 12, 24';
export const defaultThresholds={adjustedR2:0.9,span:2,extrapolated:20};
export function parseTimes(text,end){
 const parts=text.split(/[,\n]/).map(s=>s.trim()).filter(Boolean);
 if(!parts.length)throw Error('Enter at least one sampling time.');
 if(parts.length>1000)throw Error('Use at most 1,000 sampling times.');
 const times=parts.map(s=>{if(!/^(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/i.test(s))throw Error(`Invalid sampling time: ${s}. Use nonnegative numbers.`);return Number(s);});
 if(times.some(t=>!Number.isFinite(t)))throw Error('Sampling times must be finite.');
 if(new Set(times).size!==times.length)throw Error('Sampling times must not repeat (including equivalent numeric values).');
 const outside=times.filter(t=>t>end);
 if(outside.length)throw Error(`Sampling times outside simulation end (${end} h): ${outside.join(', ')} h. Extend the simulation or remove these times.`);
 return times.sort((a,b)=>a-b);
}
export function sampleProfile(p,times,lloq,zeroState='pre'){
 if(p.dosing!=='single')throw Error('NCA Sampling Mode V1 supports single dosing only.');
 const errors=validate(p);if(errors.length)throw Error(errors.join(' '));
 if(!Number.isFinite(lloq)||lloq<=0)throw Error('LLOQ must be positive and finite.');
 if(!['pre','post'].includes(zeroState))throw Error('Select a valid IV zero-time state.');
 if(!times.length||times.length>1000||new Set(times).size!==times.length||times.some(t=>!Number.isFinite(t)||t<0||t>p.end))throw Error('Sampling times must be unique, finite, nonnegative and within the simulation.');
 return [...times].sort((a,b)=>a-b).map(t=>{
  const event=p.route==='iv'&&t===0?zeroState:'sample';
  const trueC=event==='pre'?0:single(p,t);
  if(!Number.isFinite(trueC)||trueC<0)throw Error('Sample concentration exceeds numerical range.');
  return {t,trueC,c:trueC>=lloq?trueC:null,status:trueC>=lloq?'Quantifiable':'BLQ',event};
 });
}
// Explicit trust boundary: never pass trueC, LLOQ, or the model to analyzeNCA.
export function reportedSamples(samples){return samples.map(({t,c,status,event})=>({t,c,status,event}));}
export function trapezoid(a,b){
 const dt=b.t-a.t;
 if(dt<0||![a.c,b.c].every(c=>Number.isFinite(c)&&c>=0))throw Error('Invalid trapezoid endpoints.');
 if(b.c<a.c&&b.c>0){
  const relative=(a.c-b.c)/b.c;
  const logRatio=relative<0.01?Math.log1p(relative):Math.log(a.c)-Math.log(b.c);
  return dt*((a.c-b.c)/logRatio);
 }
 return dt*(a.c/2+b.c/2);
}
export function regress(points){
 if(points.length<3)return {reason:'At least three eligible quantifiable points are required.'};
 const n=points.length,x0=points[0].t;
 const xs=points.map(p=>p.t-x0),ys=points.map(p=>Math.log(p.c));
 const xm=xs.reduce((a,b)=>a+b,0)/n,ym=ys.reduce((a,b)=>a+b,0)/n;
 const sxx=xs.reduce((s,x)=>s+(x-xm)**2,0),sxy=xs.reduce((s,x,i)=>s+(x-xm)*(ys[i]-ym),0);
 const slope=sxy/sxx;
 if(!Number.isFinite(slope)||slope>=0)return {reason:'Terminal regression must have a finite negative slope.'};
 const sst=ys.reduce((s,y)=>s+(y-ym)**2,0),sse=ys.reduce((s,y,i)=>s+(y-ym-slope*(xs[i]-xm))**2,0);
 const r2=1-sse/sst,adjustedR2=1-(1-r2)*(n-1)/(n-2),lambda=-slope,halfLife=Math.LN2/lambda;
 const span=(points.at(-1).t-points[0].t)/halfLife;
 if(![r2,adjustedR2,halfLife,span].every(Number.isFinite)||halfLife<=0)return {reason:'Regression is outside numerical precision.'};
 return {lambda,halfLife,r2,adjustedR2,span,n,times:points.map(p=>p.t),anchorT:x0,logAtAnchor:ym-slope*xm,slope};
}
export function analyzeNCA(input,{route='iv',method='auto',selected=[],blqPolicy='gap',thresholds=defaultThresholds}={}){
 if(!['iv','po'].includes(route)||!['auto','manual'].includes(method)||!['gap','bridge'].includes(blqPolicy))throw Error('Invalid NCA options.');
 for(const key of Object.keys(defaultThresholds))if(!Number.isFinite(thresholds[key])||thresholds[key]<0||(key==='adjustedR2'&&thresholds[key]>1)||(key==='extrapolated'&&thresholds[key]>100))throw Error('Invalid diagnostic thresholds.');
 const data=[...input].sort((a,b)=>a.t-b.t);
 if(new Set(data.map(p=>p.t)).size!==data.length||data.some(p=>!Number.isFinite(p.t)||p.t<0||!['BLQ','Quantifiable'].includes(p.status)||!['pre','post','sample'].includes(p.event)||p.event==='pre'&&(p.t!==0||route!=='iv')||p.status==='BLQ'&&p.c!==null||p.status==='Quantifiable'&&(!Number.isFinite(p.c)||p.c<=0)))throw Error('NCA requires unique valid reported samples; BLQ concentrations must be null.');
 const exposure=data.filter(p=>p.event!=='pre'),quant=exposure.filter(p=>p.status==='Quantifiable');
 const reasons={},warnings=[];
 const out={cmax:null,tmax:null,aucLast:null,aucObserved:null,aucInf:null,extrapolated:null,lambda:null,halfLife:null,fit:null,lastTime:null,lastC:null,c0:null,c0Method:null,aucBack:null,aucBackPct:null,reasons,warnings,bridges:[],eligibleTimes:[],selectedTimes:[],method,blqPolicy,thresholds:{...thresholds}};
 if(!quant.length){for(const key of ['cmax','tmax','aucLast','lambda','halfLife','aucInf','extrapolated'])reasons[key]='No quantifiable post-dose samples.';return out;}
 const maximum=quant.reduce((a,b)=>b.c>a.c?b:a);out.cmax=maximum.c;out.tmax=maximum.t;
 const last=quant.at(-1);out.lastTime=last.t;out.lastC=last.c;
 let curve=quant.map(p=>({t:p.t,c:p.c}));
 const first=quant[0],hasOrigin=first.t===0;
 // Phoenix-style IV bolus C0: the predose 0 h BLQ is not a postdose 0+ value.
 // Estimate the dosing-time concentration from the earliest two QUANTIFIABLE
 // postdose values. Do not substitute the underlying simulator's true C(0+).
 // Estimated C0 is an integration anchor ONLY: not observed Cmax or lambda-z.
 if(route==='iv'){
  if(hasOrigin){
   out.c0=first.c;
   out.c0Method='observed-postdose';
   out.aucBack=0;
  }else if(quant.length>=2){
   const second=quant[1];
   const earlySlope=(Math.log(second.c)-Math.log(first.c))/(second.t-first.t);
   const estimated=earlySlope<0?Math.exp(Math.log(first.c)-earlySlope*first.t):first.c;
   if(Number.isFinite(estimated)&&estimated>0){
    out.c0=estimated;
    out.c0Method=earlySlope<0?'log-linear-first-two':'fallback-first';
    curve=[{t:0,c:estimated},...curve];
    out.aucBack=trapezoid(curve[0],curve[1]);
    warnings.push(earlySlope<0
     ? 'IV C0 was back-extrapolated by log-linear regression of the first two quantifiable postdose samples; 0 h predose is excluded.'
     : 'IV C0 fallback: early concentration slope was not negative, so the first quantifiable postdose concentration was used as C0.');
    if(exposure.some(p=>p.status==='BLQ'&&p.t>0&&p.t<second.t))
     warnings.push('An early postdose BLQ precedes or separates the first two quantifiable samples. The C0 estimate should be interpreted cautiously.');
   }else{
    reasons.c0='IV C0 back-extrapolation was outside numerical range.';
   }
  }else{
   reasons.c0='At least two quantifiable IV postdose points are needed to estimate C0 when 0+ is not observed.';
  }
 }
 if(route==='po'&&!hasOrigin){
  const leading=exposure.filter(p=>p.status==='BLQ'&&p.t<first.t);
  curve=[{t:0,c:0},...leading.filter(p=>p.t>0).map(p=>({t:p.t,c:0})),...curve];
  warnings.push('AUC rule: PO origin and leading BLQ samples are assigned zero for integration only; hidden concentrations are never used.');
 }
 const middle=exposure.filter(p=>p.status==='BLQ'&&p.t>first.t&&p.t<last.t);
 for(let i=1;i<quant.length;i++){
  const missing=middle.filter(p=>p.t>quant[i-1].t&&p.t<quant[i].t);
  if(missing.length)out.bridges.push({from:quant[i-1].t,to:quant[i].t,blqTimes:missing.map(p=>p.t)});
 }
 let area=0;for(let i=1;i<curve.length;i++)area+=trapezoid(curve[i-1],curve[i]);
 if(middle.length){warnings.push(`${blqPolicy==='bridge'?'Explicit bridging enabled':'Integration blocked'}: intermediate BLQ at ${middle.map(p=>p.t).join(', ')} h; quantifiable brackets ${out.bridges.map(b=>`${b.from}–${b.to} h`).join('; ')}. ${blqPolicy==='bridge'?'These intervals are integrated directly between quantifiable endpoints, ignoring BLQ values.':'No AUC is reported across these gaps.'}`);}
 if(middle.length&&blqPolicy==='gap')reasons.aucLast='Intermediate BLQ interrupts integration. Explicitly enable bridging to integrate across quantified endpoints.';
 else if(route==='iv'&&!hasOrigin&&out.c0===null){out.aucObserved=area;reasons.aucLast=reasons.c0||'Unable to determine IV C0 from reported data; AUC0-last is unavailable.';}
 else if(!Number.isFinite(area))reasons.aucLast='AUC exceeds numerical range.';
 else out.aucLast=area;
 const trailing=exposure.filter(p=>p.t>last.t&&p.status==='BLQ');
 if(trailing.length)warnings.push(`Trailing BLQ excluded: integration stops at last quantifiable sample (${last.t} h).`);
 // Regression candidates are entirely inside the final uninterrupted quantifiable run.
 const peakEnd=Math.max(...quant.filter(p=>p.c===maximum.c).map(p=>p.t));
 const priorBLQ=exposure.filter(p=>p.status==='BLQ'&&p.t<last.t).at(-1)?.t??-Infinity;
 const eligible=quant.filter(p=>(route==='po'?p.t>peakEnd:p.t>=peakEnd)&&p.t>priorBLQ);
 out.eligibleTimes=eligible.map(p=>p.t);
 let fit;
 if(method==='manual'){
  if(new Set(selected).size!==selected.length||selected.some(t=>!out.eligibleTimes.includes(t)))fit={reason:'Manual points must be unique, quantifiable, in the final uninterrupted run and after observed Tmax (IV may include Tmax).'};
  else fit=regress(eligible.filter(p=>selected.includes(p.t)));
 }else{
  // Enumerate suffixes ending at the last quantifiable point; adjusted-R² ties within 1e-4 prefer more points.
  const candidates=[];
  for(let n=3;n<=eligible.length;n++){const f=regress(eligible.slice(-n));if(!f.reason)candidates.push(f);}
  if(candidates.length){const best=Math.max(...candidates.map(f=>f.adjustedR2));fit=candidates.filter(f=>best-f.adjustedR2<=1e-4).sort((a,b)=>b.n-a.n)[0];}
  else fit={reason:'No negative-slope fit with at least three eligible terminal points. PO excludes observed Tmax; regression cannot cross BLQ gaps.'};
 }
 if(fit.reason){reasons.lambda=fit.reason;reasons.halfLife=fit.reason;}
 else{
  out.fit=fit;out.lambda=fit.lambda;out.halfLife=fit.halfLife;out.selectedTimes=fit.times;
  if(fit.adjustedR2<thresholds.adjustedR2)warnings.push(`Diagnostic: adjusted R² < ${thresholds.adjustedR2}.`);
  if(fit.span<thresholds.span)warnings.push(`Diagnostic: terminal span < ${thresholds.span} estimated half-lives.`);
 }
 if(out.aucLast===null||!out.fit){reasons.aucInf=out.aucLast===null?reasons.aucLast:reasons.lambda;reasons.extrapolated=reasons.aucInf;}
 else{
  const tail=last.c/out.lambda,total=out.aucLast+tail;
  if(!Number.isFinite(total)||total<=0){reasons.aucInf='AUC extrapolation exceeds numerical range.';reasons.extrapolated=reasons.aucInf;}
  else{out.aucInf=total;out.extrapolated=tail/total*100;if(route==='iv'&&out.aucBack!==null)out.aucBackPct=100*out.aucBack/total;if(out.extrapolated>thresholds.extrapolated)warnings.push(`Diagnostic: AUC extrapolated > ${thresholds.extrapolated}%.`);}
 }
 warnings.push('Exploratory diagnostics only; a high R² does not establish a reliable terminal phase. Thresholds are not universal regulatory acceptance criteria.');
 return out;
}
export function compareNCA(nca,theory,route){
 const bench=route==='po'?theory.apparentHalfLife:theory.halfLife;
 return [['Cmax','cmax',theory.global.c,'ng/mL'],['Tmax','tmax',theory.global.t,'h'],['AUC0-inf','aucInf',theory.theoreticalAUC,'ng·h/mL'],[route==='po'?'Apparent terminal half-life (PO benchmark)':'Disposition terminal half-life','halfLife',bench,'h']].map(([label,key,value,unit])=>({label,key,unit,theory:value,nca:nca[key],bias:Number.isFinite(value)&&value!==0&&Number.isFinite(nca[key])?100*(nca[key]-value)/value:null,reason:nca.reasons[key]||''}));
}
export function ncaCSV(samples,nca,comparison,metadata={}){
 const result={iv_c0_ng_ml:nca.c0,iv_c0_method:nca.c0Method,iv_auc_back_extrapolated_ng_h_ml:nca.aucBack,iv_auc_back_extrapolated_percent:nca.aucBackPct,observed_cmax_ng_ml:nca.cmax,observed_tmax_h:nca.tmax,auc0_last_ng_h_ml:nca.aucLast,auc_first_last_ng_h_ml:nca.aucObserved,lambda_z_per_h:nca.lambda,terminal_half_life_h:nca.halfLife,auc0_inf_ng_h_ml:nca.aucInf,auc_extrapolated_pct:nca.extrapolated,lambda_z_times_h:nca.selectedTimes.join(';'),lambda_z_n:nca.fit?.n??null,r_squared:nca.fit?.r2??null,adjusted_r_squared:nca.fit?.adjustedR2??null,terminal_span_half_lives:nca.fit?.span??null,reasons:JSON.stringify(nca.reasons),warnings:nca.warnings.join(' | '),comparison:JSON.stringify(comparison),integration_brackets:JSON.stringify(nca.bridges),...metadata};
 const header=['Sampling time (h)','Underlying theoretical concentration (ng/mL)','Reported concentration (ng/mL) or BLQ','Quantitation status','Event state','Terminal regression selection',...Object.keys(result)];
 const quote=v=>`"${String(v??'N/A').replaceAll('"','""')}"`;
 return [header,...samples.map(s=>[s.t,s.trueC,s.status==='BLQ'?'BLQ':s.c,s.status,s.event,nca.selectedTimes.includes(s.t)?'selected':'not selected',...Object.values(result)])].map(row=>row.map(quote).join(',')).join('\r\n');
}
