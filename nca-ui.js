import {parseTimes,sampleProfile,reportedSamples,analyzeNCA,compareNCA,ncaCSV} from './nca.js';
const $=id=>document.getElementById(id);
const fmt=v=>Number.isFinite(v)?new Intl.NumberFormat('en-US',{maximumSignificantDigits:5}).format(v):'N/A';
export function createSamplingMode(onChange){
 let state=null,selected=new Set();
 const controls=['nca-enabled','nca-times','nca-lloq','nca-zero','nca-method','nca-blq-policy','nca-r2','nca-span','nca-extrap'];
 for(const id of controls)$(id).addEventListener('input',onChange);
 $('nca-content').addEventListener('change',e=>{if(e.target.matches('[data-nca-time]')){const t=Number(e.target.dataset.ncaTime);if(e.target.checked)selected.add(t);else selected.delete(t);onChange();$('nca-content').querySelector(`[data-nca-time="${t}"]`)?.focus();}});
 $('nca-export').onclick=()=>{
  if(!state)return;const {samples,nca,comparison,metadata}=state;
  const url=URL.createObjectURL(new Blob([ncaCSV(samples,nca,comparison,metadata)],{type:'text/csv;charset=utf-8'}));const a=document.createElement('a');a.href=url;a.download='pk-nca-sampling.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
 };
 function refresh(p,theory){
  state=null;
  if($('tooltip').classList.contains('nca-tooltip'))$('tooltip').hidden=true;
  const enabled=$('nca-enabled').checked,single=p.dosing==='single',active=enabled&&single;
  $('nca-enabled').disabled=!single;$('nca-settings').hidden=!active;$('nca-panel').hidden=!active;$('nca-legend').hidden=true;$('nca-export').disabled=true;
  $('nca-availability').textContent=single?'Single-dose sampling and noncompartmental analysis.':'Unavailable for multiple doses. Original PK simulation remains active.';
  if(!active)return;
  $('nca-zero-setting').hidden=p.route!=='iv';
  try{
   if(!theory)throw Error('Correct the PK model configuration first.');
   const times=parseTimes($('nca-times').value,p.end),lloq=Number($('nca-lloq').value),zero=$('nca-zero').value;
   const samples=sampleProfile(p,times,lloq,zero),method=$('nca-method').value,blqPolicy=$('nca-blq-policy').value;
   const thresholds={adjustedR2:readNumber('nca-r2'),span:readNumber('nca-span'),extrapolated:readNumber('nca-extrap')};
   const nca=analyzeNCA(reportedSamples(samples),{route:p.route,method,selected:[...selected],blqPolicy,thresholds});
   const comparison=compareNCA(nca,theory,p.route);
   const metadata={lloq_ng_ml:lloq,model:p.model,route:p.route,dosing:p.dosing,dose_mg_kg:p.dose,cl_l_h_kg:p.cl,...(p.model==='one'?{v_l_kg:p.v}:{vc_l_kg:p.vc,vp_l_kg:p.vp,q_l_h_kg:p.q}),...(p.route==='po'?{f:p.f,ka_per_h:p.ka}:{}),simulation_end_h:p.end,iv_zero_state:zero,terminal_method:method,manual_requested_times_h:[...selected].sort((a,b)=>a-b).join(';'),blq_policy:blqPolicy,diagnostic_thresholds:JSON.stringify(thresholds),theoretical_disposition_half_life_h:theory.halfLife,theoretical_apparent_half_life_h:theory.apparentHalfLife};
   state={samples,nca,comparison,metadata,lloq};
   $('nca-error').hidden=true;$('nca-export').disabled=false;$('nca-legend').hidden=false;
   render(theory,p.route);
  }catch(e){$('nca-error').hidden=false;$('nca-error').textContent=e.message;$('nca-content').replaceChildren();}
 }
 function readNumber(id){return $(id).value.trim()===''?NaN:Number($(id).value);}
 function render(theory,route){
  const {samples,nca,comparison}=state,fit=nca.fit,manual=$('nca-method').value==='manual';
  const value=(key,unit)=>`<strong>${fmt(nca[key])}</strong> ${unit}${nca[key]===null?`<small>N/A: ${nca.reasons[key]}</small>`:''}`;
  const comparisonRows=comparison.map(c=>`<tr><th scope="row">${c.label}<small>${c.unit}</small></th><td>${fmt(c.nca)}${c.reason?`<small>${c.reason}</small>`:''}</td><td>${fmt(c.theory)}</td><td>${c.bias===null?'N/A':fmt(c.bias)+'%'}${c.bias===null&&c.theory===0?'<small>Theoretical value is zero.</small>':''}</td></tr>`).join('');
  const bias=comparison.find(c=>c.key==='halfLife')?.bias;
  const modelWarning=bias!==null&&bias!==undefined&&Math.abs(bias)>20?'<p class="nca-warning">Model comparison: estimated terminal half-life differs from its observable asymptotic benchmark by more than 20%. This sampling window may not cover the terminal phase, even when R² is high.</p>':'';
  $('nca-content').innerHTML=`<div class="nca-table-wrap"><table class="nca-table"><caption>NCA estimates versus full single-dose model benchmarks</caption><thead><tr><th>Parameter</th><th>NCA estimate</th><th>Theoretical model</th><th>Relative bias</th></tr></thead><tbody>${comparisonRows}</tbody></table></div>${route==='po'?`<p class="field-note">Disposition half-life: ${fmt(theory.halfLife)} h. Apparent PO terminal benchmark: ${fmt(theory.apparentHalfLife)} h${theory.apparentHalfLife>theory.halfLife?' · absorption-limited (flip-flop)':''}.</p>`:''}
  <div class="nca-metrics"><div>AUC0-last ${value('aucLast','ng·h/mL')}</div><div>λz ${value('lambda','h⁻¹')}</div><div>AUC extrapolated ${value('extrapolated','%')}</div></div>
  ${nca.aucObserved!==null?`<p class="nca-warning">Partial AUCfirst-last = ${fmt(nca.aucObserved)} ng·h/mL. This excludes the unsampled IV interval from 0 h; it is not AUC0-last.</p>`:''}
  <div class="nca-fit-stats"><span>λz points: ${fit?fit.times.map(fmt).join(', ')+' h':'N/A'}</span><span>Number of points: ${fit?.n??'N/A'}</span><span>R²: ${fmt(fit?.r2)}</span><span>Adjusted R²: ${fmt(fit?.adjustedR2)}</span><span>Terminal span: ${fmt(fit?.span)} half-lives</span></div>
  <div class="nca-warnings">${nca.warnings.map(w=>`<p>${w}</p>`).join('')}</div>${modelWarning}
  <div class="nca-table-wrap"><table class="nca-table" id="nca-samples"><caption>Reported samples · ${manual?'Select eligible rows for manual regression':'Automatic terminal selection'} · BLQ values are never assigned LLOQ</caption><thead><tr><th>λz selection</th><th>Time (h)</th><th>Reported (ng/mL)</th><th>State</th></tr></thead><tbody>${samples.map(s=>`<tr class="${nca.selectedTimes.includes(s.t)?'nca-selected':''}"><td>${manual?`<input type="checkbox" aria-label="Use ${s.t} h in terminal regression" data-nca-time="${s.t}" ${selected.has(s.t)?'checked':''} ${!nca.eligibleTimes.includes(s.t)&&!selected.has(s.t)?'disabled':''}>`:nca.selectedTimes.includes(s.t)?'Selected':'—'}</td><td>${s.t}</td><td>${s.status==='BLQ'?'<span class="nca-blq">BLQ</span>':fmt(s.c)}</td><td>${s.event==='pre'?'0 h predose · excluded':s.event==='post'?'0 h postdose (0+)':s.status}</td></tr>`).join('')}</tbody></table></div>`;
 }
 function axisValues(){return state?[state.lloq,...state.samples.filter(s=>s.c!==null).map(s=>s.c)]:[];}
 function samplePosition(s,g,log){return {x:g.x(s.t),y:s.status==='BLQ'?g.top-13:g.y(s.c)};}
 function draw(ctx,g,log){
  if(!state)return;const {samples,nca,lloq}=state;
  ctx.save();ctx.beginPath();ctx.rect(g.left-1,g.top,g.right-g.left+2,g.bottom-g.top);ctx.clip();
  ctx.strokeStyle='#bb7c30';ctx.fillStyle='#a06a26';ctx.lineWidth=1;ctx.setLineDash([5,4]);
  ctx.beginPath();ctx.moveTo(g.left,g.y(lloq));ctx.lineTo(g.right,g.y(lloq));ctx.stroke();ctx.font='10px sans-serif';ctx.textAlign='right';ctx.fillText(`LLOQ ${fmt(lloq)} ng/mL`,g.right-5,g.y(lloq)-5);
  if(nca.fit){const f=nca.fit;ctx.strokeStyle='#7c60a1';ctx.lineWidth=2;ctx.setLineDash([6,4]);ctx.beginPath();for(let i=0;i<=100;i++){const t=f.times[0]+(f.times.at(-1)-f.times[0])*i/100,c=Math.exp(f.logAtAnchor+f.slope*(t-f.anchorT));if(i===0)ctx.moveTo(g.x(t),g.y(c));else ctx.lineTo(g.x(t),g.y(c));}ctx.stroke();}
  ctx.setLineDash([]);ctx.restore();
  for(const s of samples){const pos=samplePosition(s,g,log);ctx.save();ctx.lineWidth=2;
   if(s.status==='BLQ'){ctx.strokeStyle='#b87830';ctx.beginPath();ctx.moveTo(pos.x-4,pos.y-4);ctx.lineTo(pos.x+4,pos.y+4);ctx.moveTo(pos.x-4,pos.y+4);ctx.lineTo(pos.x+4,pos.y-4);ctx.stroke();}
   else{ctx.fillStyle='white';ctx.strokeStyle=nca.selectedTimes.includes(s.t)?'#7c60a1':'#205c87';ctx.beginPath();ctx.arc(pos.x,pos.y,nca.selectedTimes.includes(s.t)?5:4,0,2*Math.PI);ctx.fill();ctx.stroke();if(nca.selectedTimes.includes(s.t)){ctx.fillStyle='#7c60a1';ctx.beginPath();ctx.arc(pos.x,pos.y,2,0,2*Math.PI);ctx.fill();}}
   ctx.restore();
  }
 }
 function tooltip(mx,my,g,log){
  if(!state)return null;let nearest=null,distance=10;
  for(const s of state.samples){const pos=samplePosition(s,g,log),d=Math.hypot(mx-pos.x,my-pos.y);if(d<distance){distance=d;nearest=s;}}
  if(!nearest)return null;
  return `${nearest.t} h · ${nearest.status==='BLQ'?'BLQ (< LLOQ; no reported concentration)':`${fmt(nearest.c)} ng/mL · Quantifiable`}${nearest.event==='pre'?' · predose (excluded)':nearest.event==='post'?' · postdose (0+)':''}${state.nca.selectedTimes.includes(nearest.t)?' · λz selected':''}`;
 }
 return {refresh,draw,tooltip,axisValues,suspend:()=>{state=null;$('nca-export').disabled=true;},caption:()=>state?'NCA: blue = samples; purple = fit; BLQ = top-rail x (no Y).':''};
}
