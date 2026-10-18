import {defaults,validate,derive,simulate,results,csv,concentration,doseTimes} from './pk-model.js';
import {createSamplingMode} from './nca-ui.js';
const $=id=>document.getElementById(id);
let p={...defaults},rows=[],stats,log=false,geometry,timer;
const sampling=createSamplingMode(()=>{clearTimeout(timer);update();});
const definitions={
 dose:['Dose','mg/kg',0,0.1,0,30],cl:['Clearance · CL','L/h/kg',0.000001,0.01,0.01,1],
 v:['Volume · V','L/kg',0.000001,0.1,0.1,5],vc:['Central volume · Vc','L/kg',0.000001,0.1],vp:['Peripheral volume · Vp','L/kg',0.000001,0.1],q:['Intercompartmental CL · Q','L/h/kg',0.000001,0.01],
 f:['Bioavailability · F','fraction',0,0.01,0,1],ka:['Absorption rate · ka','1/h',0.000001,0.05,0.01,5],
 tau:['Dosing interval · τ','h',0.000001,1],n:['Number of doses','doses',1,1],end:['Simulation end','h',0.000001,1]
};
function field(key){const [label,unit,min,step,lo,hi]=definitions[key];return `<div class="parameter"><label for="${key}">${label}<span class="unit">${unit}</span></label><div class="${lo!==undefined?'input-row':''}"><input id="${key}" name="${key}" type="number" min="${min}" ${key==='f'?'max="1"':key==='n'?'max="100"':''} step="${step}" value="${p[key]}">${lo!==undefined?`<input type="range" aria-label="${label} slider" data-for="${key}" min="${lo}" max="${hi}" step="${step}" value="${p[key]}">`:''}</div></div>`;}
function renderFields(){
 for(const key of ['model','route','dosing'])$(key).value=p[key];
 $('pk-fields').innerHTML=['dose','cl',...(p.model==='one'?['v']:['vc','vp','q']),...(p.route==='po'?['f','ka']:[])].map(field).join('');
 $('regimen-fields').innerHTML=[...(p.dosing==='multiple'?['tau','n']:[]),'end'].map(field).join('');
 $('dose-note').textContent=p.dosing==='multiple'?`Doses at 0, τ, …, (N − 1)τ. End time must exceed the last dose. Maximum 100 doses.`:'A single dose is administered at time zero.';
}
const fmt=v=>Number.isFinite(v)?new Intl.NumberFormat('en-US',{maximumSignificantDigits:5}).format(v):'N/A';
function metric(label,value,unit,extra=''){return `<div class="metric"><div class="metric-label">${label}</div><div class="metric-value">${value}</div><div class="metric-unit">${unit}</div>${extra?`<div class="metric-extra">${extra}</div>`:''}</div>`;}
function renderResults(){
 const r=stats,multiple=p.dosing==='multiple';
 let cards=[metric(multiple?'Global simulated Cmax · estimated':'Cmax · estimated',fmt(r.global.c),'ng/mL',`Tmax ${fmt(r.global.t)} h${multiple?' from first dose':''}`)];
 if(multiple){for(const [label,k]of [['First interval',r.first],['Final dose interval',r.last]])cards.push(metric(`${label} Cmax · estimated`,k?fmt(k.c):'N/A','ng/mL',k?`Tmax ${fmt(k.t)} h after dose`:'Complete interval not simulated'));
 cards.push(metric('Simulated AUC₀–end · estimated',fmt(r.auc),'ng·h/mL'));}
 cards.push(metric(multiple?'Cumulative AUC₀–∞ · theoretical':'AUC₀–∞ · theoretical',fmt(r.theoreticalAUC),'ng·h/mL'));
 cards.push(metric('Disposition terminal half-life',fmt(r.halfLife),'h'));
 if(p.route==='po')cards.push(metric('Apparent oral terminal half-life',fmt(r.apparentHalfLife),'h',p.ka<r.terminal?'Absorption-limited · flip-flop':'Disposition-limited'));
 $('results').innerHTML=cards.join('');
 $('derived').innerHTML=(p.model==='one'?[['kel',r.kel]]:[['k10',r.k10],['k12',r.k12],['k21',r.k21],['α',r.alpha],['β',r.beta]]).map(([name,v])=>`<span>${name} = ${fmt(v)} h⁻¹</span>`).join('');
}
function update(){
 const errors=validate(p);
 if(!errors.length){try{const d=derive(p);if(![d.terminal,d.halfLife,d.apparentHalfLife,...(p.model==='two'?[d.alpha,d.a,d.b]:[])].every(Number.isFinite)||d.terminal<=0)throw Error('Parameters exceed numerical range.');rows=simulate(p);stats=results(p);renderResults();}catch(e){errors.push(e.message);}}
 $('error').hidden=!errors.length;$('error').textContent=errors.join(' ');
 for(const id of ['export-csv','export-png'])$(id).disabled=!!errors.length;
 $('model-summary').textContent=`${p.model==='one'?'One':'Two'}-compartment · ${p.route==='iv'?'IV bolus':'Oral absorption'} · ${p.dosing==='single'?'Single dose':`${p.n} doses every ${p.tau} h`}`;
 if(errors.length){rows=[];stats=null;$('results').innerHTML='<p>Correct the configuration to calculate results.</p>';$('derived').innerHTML='';}
 $('plot-note').textContent=errors.length?'Invalid configuration':log?'Log scale · Zero concentrations omitted':'Analytical model · Adaptive time grid';
 sampling.refresh(p,stats);
 draw();
}
$('controls').addEventListener('input',e=>{
 const key=e.target.dataset.for||e.target.name;if(!key)return;
 if(['model','route','dosing'].includes(key)){p[key]=e.target.value;renderFields();update();return;}
 p[key]=e.target.value===''?NaN:Number(e.target.value);
 if(e.target.dataset.for)$(key).value=e.target.value;else{const slider=document.querySelector(`[data-for="${key}"]`);if(slider)slider.value=e.target.value;}
 for(const id of ['export-csv','export-png'])$(id).disabled=true;
 sampling.suspend();
 clearTimeout(timer);timer=setTimeout(update,90);
});
$('reset').onclick=()=>{clearTimeout(timer);p={...defaults};renderFields();update();};
for(const id of ['linear','log'])$(id).onclick=()=>{log=id==='log';$('linear').setAttribute('aria-pressed',String(!log));$('log').setAttribute('aria-pressed',String(log));update();};
function draw(){
 const canvas=$('chart'),rect=canvas.getBoundingClientRect(),w=rect.width,h=rect.height,dpr=devicePixelRatio||1;
 canvas.width=w*dpr;canvas.height=h*dpr;const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);ctx.fillStyle='#fff';ctx.fillRect(0,0,w,h);
 const left=69,right=w-22,top=27,bottom=h-48,pw=right-left,ph=bottom-top;
 const values=[...rows.map(r=>r.c),...sampling.axisValues()].filter(c=>c>0),max=values.reduce((a,b)=>Math.max(a,b),0)||1;
 let ymin=0,ymax=Number.isFinite(max*1.13)?max*1.13:max;
 if(log){ymax=Math.ceil(Math.log10(max))+0.15;ymin=Math.floor(Math.log10(values.reduce((a,b)=>Math.min(a,b),max)));if(ymax-ymin<1)ymin=ymax-1;}
 const x=t=>left+t/p.end*pw,y=c=>bottom-((log?Math.log10(c):c)-ymin)/(ymax-ymin)*ph;
 geometry={left,right,top,bottom,x,y,w,h};
 ctx.font='10px -apple-system, sans-serif';ctx.lineWidth=1;ctx.strokeStyle='#e9efeb';ctx.fillStyle='#90a29b';
 const ticks=log?Array.from({length:Math.min(12,Math.ceil(ymax-ymin)+1)},(_,i)=>ymin+i*Math.max(1,Math.ceil((ymax-ymin)/11))).filter(v=>v<=ymax):Array.from({length:6},(_,i)=>ymax*i/5);
 for(const v of ticks){const yy=log?y(10**v):y(v);ctx.beginPath();ctx.moveTo(left,yy);ctx.lineTo(right,yy);ctx.stroke();ctx.textAlign='right';ctx.fillText(log?`10^${v}`:fmt(v),left-12,yy+3);}
 for(let i=0;i<=6;i++){const t=p.end*i/6,xx=x(t);ctx.textAlign='center';ctx.fillText(fmt(t),xx,bottom+20);}
 ctx.save();ctx.translate(17,(top+bottom)/2);ctx.rotate(-Math.PI/2);ctx.fillStyle='#728b81';ctx.fillText('Concentration (ng/mL)',0,0);ctx.restore();ctx.textAlign='center';ctx.fillText('Time (h)',(left+right)/2,h-7);
 if(!rows.length)return;
 ctx.save();ctx.beginPath();ctx.rect(left,top,pw,ph);ctx.clip();
 ctx.setLineDash([3,5]);ctx.strokeStyle='#afc5b8';
 for(const td of doseTimes(p)){if(td>p.end)continue;ctx.beginPath();ctx.moveTo(x(td),top);ctx.lineTo(x(td),bottom);ctx.stroke();}
 ctx.setLineDash([]);
 if(!log){const gradient=ctx.createLinearGradient(0,top,0,bottom);gradient.addColorStop(0,'#1a947624');gradient.addColorStop(1,'#1a947602');ctx.beginPath();ctx.moveTo(left,bottom);for(const r of rows)ctx.lineTo(x(r.t),y(r.c));ctx.lineTo(right,bottom);ctx.closePath();ctx.fillStyle=gradient;ctx.fill();}
 ctx.strokeStyle='#16836f';ctx.lineWidth=2;ctx.lineJoin='round';ctx.beginPath();let active=false;
 for(const r of rows){if(log&&r.c<=0){active=false;continue;}if(active)ctx.lineTo(x(r.t),y(r.c));else{ctx.moveTo(x(r.t),y(r.c));active=true;}}
 ctx.stroke();ctx.restore();
 sampling.draw(ctx,geometry,log);
}
$('chart').addEventListener('pointermove',e=>{
 if(!rows.length||!geometry)return;const rect=$('chart').getBoundingClientRect(),mx=e.clientX-rect.left,my=e.clientY-rect.top,g=geometry;
 const sampleTip=sampling.tooltip(mx,my,g,log);
 if(sampleTip){const tip=$('tooltip');tip.classList.add('nca-tooltip');tip.textContent=sampleTip;tip.hidden=false;tip.style.left=`${Math.max(0,Math.min(mx+12,g.w-tip.offsetWidth-5))}px`;tip.style.top=`${Math.max(0,my-48)}px`;return;}
 if(mx<g.left||mx>g.right||my<g.top||my>g.bottom){$('tooltip').hidden=true;return;}
 let t=(mx-g.left)/(g.right-g.left)*p.end,state='post';
 const near=doseTimes(p).find(td=>Math.abs(g.x(td)-mx)<5);
 if(near!==undefined){t=near;if(p.route==='iv'){const pre=concentration(p,t,'pre'),post=concentration(p,t);if(pre>0||!log)state=Math.abs(g.y(pre)-my)<Math.abs(g.y(post)-my)?'pre':'post';}}
 const c=concentration(p,t,state),tip=$('tooltip');tip.classList.remove('nca-tooltip');tip.textContent=`${fmt(t)} h · ${fmt(c)} ng/mL${near!==undefined?` · ${state}-dose`:''}`;tip.hidden=false;tip.style.left=`${Math.max(0,Math.min(mx+12,g.w-tip.offsetWidth-5))}px`;tip.style.top=`${Math.max(0,my-48)}px`;
});
$('chart').addEventListener('pointerleave',()=>$('tooltip').hidden=true);
function download(url,name){const a=document.createElement('a');a.href=url;a.download=name;a.click();}
$('export-csv').onclick=()=>{const url=URL.createObjectURL(new Blob([csv(p,rows)],{type:'text/csv;charset=utf-8'}));download(url,'pk-simulation.csv');setTimeout(()=>URL.revokeObjectURL(url),1000);};
$('export-png').onclick=()=>{draw();const source=$('chart'),canvas=document.createElement('canvas');canvas.width=source.width;canvas.height=source.height+100*(devicePixelRatio||1);const ctx=canvas.getContext('2d'),dpr=devicePixelRatio||1;ctx.fillStyle='#fff';ctx.fillRect(0,0,canvas.width,canvas.height);ctx.drawImage(source,0,60*dpr);ctx.scale(dpr,dpr);ctx.fillStyle='#234d43';ctx.font='bold 16px sans-serif';ctx.fillText('PK Simulator · Concentration over time',22,25);ctx.font='11px sans-serif';ctx.fillText($('model-summary').textContent+` · ${log?'Log':'Linear'} scale`,22,45);ctx.fillStyle='#748d80';if(sampling.caption()){ctx.font='9px sans-serif';ctx.fillText(sampling.caption(),22,canvas.height/dpr-29);ctx.font='11px sans-serif';}ctx.fillText('Simulation and education only · Dose '+p.dose+' mg/kg · CL '+p.cl+' L/h/kg',22,canvas.height/dpr-14);download(canvas.toDataURL('image/png'),'pk-profile.png');};
new ResizeObserver(draw).observe($('chart-wrap'));
renderFields();update();
