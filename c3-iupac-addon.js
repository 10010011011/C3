/* c3-iupac-addon.js — extra IUPAC-name coverage for Chem3D (C3)
   Load AFTER the main script (see bottom of this file for the one line to add to index.html).
   It wraps Chem3DCore.fromName: names it fully understands are built here, everything else
   (and every spelling-mistake case) falls through to the original parser untouched.

   Adds (all derived from a bond graph, no per-molecule tables):
   • rings: benzene, cycloalkanes/enes, naphthalene, azulene, anthracene, phenanthrene, pyridine, furan, thiophene, pyrrole
   • any number of locants + multipliers for every suffix  (benzene-1,2,3-triol, naphthalene-1,5-diol, butane-2,3-dione …)
   • primary suffix (ane/ene/yne, with locants) + secondary suffix (ol/one/al/oic acid/amine/amide/nitrile/carboxylic acid/…)
     on main chains AND rings:  but-3-en-2-one, pent-4-en-2-one, but-2-enal, prop-2-yn-1-ol, cyclohex-2-en-1-one …
   • nested / compound substituents with the same grammar as the parent: (ethanoyloxy), (methoxycarbonyl), (2-chloroethyl),
     propan-2-yl, naphthalen-1-yl, dimethylamino, N-substituents (N-phenylethanamide, N,N-dimethylaniline) …
   • esters, anhydrides (acyclic, mixed, cyclic), acyl halides
   • retained names: catechol, resorcinol, quinol, pyrogallol, cresols, naphthols, acetanilide, cinnamic acid,
     phthalic/isophthalic/terephthalic acid, salicylic acid, succinic/maleic/… acid, acetone …  (+ small typo repair)
   Stereo descriptors such as (E)/(Z)/(R)/(S) are accepted and ignored (the 3D builder has no stereo input).            */
(function(root){
'use strict';
/* ───────────────────────── tiny molecular graph ───────────────────────── */
const MAXV={C:4,N:3,O:2,S:6,P:5,F:1,Cl:1,Br:1,I:1};
const NP=()=>{const e=new Error('np');e.np=1;return e};      // "this is not a name I can read" → fall back
const CH=m=>{const e=new Error(m);e.hard=1;return e};         // "I read it, but it is chemically impossible" (beats the old parser)
const SOFT=m=>new Error(m);                                   // "I read it, but it is odd" (old parser may still make sense of it)
const Mol=()=>({a:[],b:[]});
function smi(s){                                              // mini SMILES reader for group definitions
  const m=Mol(),st=[],rg={};let prev=-1,bo=1;
  for(let i=0;i<s.length;i++){const c=s[i];
    if(c==='(')st.push(prev);else if(c===')')prev=st.pop();else if(c==='=')bo=2;else if(c==='#')bo=3;
    else if(/\d/.test(c)){if(rg[c]!==undefined){const[a,o]=rg[c];m.b.push([a,prev,Math.max(o,bo),(m.a[a].ar&&m.a[prev].ar&&o===1&&bo===1)?1:0]);delete rg[c]}else rg[c]=[prev,bo];bo=1}
    else{let el,ar=false,q=0,h;
      if(c==='['){const e=s.indexOf(']',i),mm=s.slice(i+1,e).match(/^([A-Z][a-z]?|[cnos])(H(\d?))?([+-]+)?/);el=mm[1];ar=/[a-z]/.test(el[0]);el=el[0].toUpperCase()+el.slice(1);
        h=mm[2]?(mm[3]?+mm[3]:1):undefined;q=mm[4]?(mm[4][0]==='+'?1:-1)*mm[4].length:0;i=e}
      else if(c==='C'&&s[i+1]==='l'){el='Cl';i++}else if(c==='B'&&s[i+1]==='r'){el='Br';i++}
      else if(/[A-Z]/.test(c))el=c;else if('cnos'.includes(c)){el=c.toUpperCase();ar=true}else continue;
      const k=m.a.length;m.a.push({el,ar,q,h});if(prev>=0)m.b.push([prev,k,bo,(ar&&m.a[prev].ar&&bo===1)?1:0]);prev=k;bo=1}}
  return m}
const G=(s,ord,x)=>Object.assign({m:smi(s),att:0,ord:ord||1},x||{});   // a group: graph + attachment atom + bond order to parent
function merge(h,f){const off=h.a.length;f.a.forEach(x=>h.a.push(Object.assign({},x)));f.b.forEach(b=>b&&h.b.push([b[0]+off,b[1]+off,b[2],b[3]]));return off}
function attach(h,i,g){const off=merge(h,g.m);h.b.push([i,off+g.att,g.ord||1,0]);return off}
function used(m,i){let u=(m.a[i].h||0)+(m.a[i].res||0);for(const b of m.b)if(b&&(b[0]===i||b[1]===i))u+=b[3]?1.5:b[2];return u}
const capm=(m,i)=>(MAXV[m.a[i].el]||4)-used(m,i);
const nbrs=(m,i)=>{const r=[];m.b.forEach(b=>{if(!b)return;if(b[0]===i)r.push(b[1]);else if(b[1]===i)r.push(b[0])});return r};
function kill(m,i){m.a[i].dead=true;m.b.forEach((b,k)=>{if(b&&(b[0]===i||b[1]===i))m.b[k]=null})}

/* ───────────────────────── SMILES writer ───────────────────────── */
function toSmiles(m){
  const n=m.a.length,adj=m.a.map(()=>[]);
  m.b.forEach((b,bi)=>{if(!b)return;adj[b[0]].push([b[1],bi]);adj[b[1]].push([b[0],bi])});
  const root0=m.a.findIndex(a=>!a.dead),state=Array(n).fill(0),tree=m.a.map(()=>[]),back=[];
  (function dfs(u,pb){state[u]=1;for(const[v,bi]of adj[u]){if(bi===pb)continue;if(state[v]===0){tree[u].push([v,bi]);dfs(v,bi)}else if(state[v]===1)back.push([u,v,bi])}state[u]=2})(root0,-1);
  const open={},close={};back.forEach(([u,v,bi])=>{(open[v]=open[v]||[]).push(bi);(close[u]=close[u]||[]).push(bi)});
  const sym=bi=>{const b=m.b[bi];if(b[2]===2)return'=';if(b[2]===3)return'#';if(b[3])return'';return(m.a[b[0]].ar&&m.a[b[1]].ar)?'-':''};
  const lab=a=>{const e=a.ar?a.el.toLowerCase():a.el;if(a.q||a.h!==undefined)return'['+e+(a.h?'H'+(a.h>1?a.h:''):'')+(a.q?(a.q>0?'+':'-'):'')+']';return e};
  const dig={},free=[9,8,7,6,5,4,3,2,1];
  function emit(u){let s=lab(m.a[u]);
    (close[u]||[]).forEach(bi=>{s+=dig[bi];free.push(+dig[bi])});
    (open[u]||[]).forEach(bi=>{if(!free.length)throw CH('Too many ring closures to write');const d=free.pop();dig[bi]=String(d);s+=sym(bi)+d});
    const ch=tree[u];ch.forEach(([v,bi],k)=>{const part=sym(bi)+emit(v);s+=k<ch.length-1?'('+part+')':part});return s}
  return emit(root0)}

/* ───────────────────────── parent skeletons ───────────────────────── */
const TM={
 benzene:{L:'1 2 3 4 5 6',E:'1-2 2-3 3-4 4-5 5-6 6-1'},
 pyridine:{L:'1 2 3 4 5 6',E:'1-2 2-3 3-4 4-5 5-6 6-1',het:{1:'N'}},
 furan:{L:'1 2 3 4 5',E:'1-2 2-3 3-4 4-5 5-1',het:{1:'O'}},
 thiophene:{L:'1 2 3 4 5',E:'1-2 2-3 3-4 4-5 5-1',het:{1:'S'}},
 pyrrole:{L:'1 2 3 4 5',E:'1-2 2-3 3-4 4-5 5-1',het:{1:'N'},h:{1:1}},
 naphthalene:{L:'1 2 3 4 4a 5 6 7 8 8a',E:'1-2 2-3 3-4 4-4a 4a-5 5-6 6-7 7-8 8-8a 8a-1 4a-8a'},
 azulene:{L:'1 2 3 3a 4 5 6 7 8 8a',E:'1-2 2-3 3-3a 3a-4 4-5 5-6 6-7 7-8 8-8a 8a-1 3a-8a'},
 anthracene:{L:'1 2 3 4 4a 10 10a 5 6 7 8 8a 9 9a',E:'1-2 2-3 3-4 4-4a 4a-10 10-10a 10a-5 5-6 6-7 7-8 8-8a 8a-9 9-9a 9a-1 4a-9a 8a-10a'},
 phenanthrene:{L:'1 2 3 4 4a 4b 5 6 7 8 8a 9 10 10a',E:'1-2 2-3 3-4 4-4a 4a-4b 4b-5 5-6 6-7 7-8 8-8a 8a-9 9-10 10-10a 10a-1 4a-10a 4b-8a'}};
function tmpl(name){const t=TM[name],order=t.L.split(' '),m=Mol(),loc={};
  order.forEach(l=>{loc[l]=m.a.length;m.a.push({el:(t.het&&t.het[l])||'C',ar:true,q:0,h:t.h&&t.h[l]})});
  t.E.split(' ').forEach(e=>{const[x,y]=e.split('-');m.b.push([loc[x],loc[y],1,1])});
  return{m,loc,order,ar:true}}
function chainT(n){const m=Mol(),loc={},order=[];for(let i=1;i<=n;i++){loc[i]=m.a.length;order.push(String(i));m.a.push({el:'C',ar:false,q:0});if(i>1)m.b.push([loc[i-1],loc[i],1,0])}return{m,loc,order,n}}
function ringT(n){if(n<3)throw CH('A ring needs at least 3 carbons');const t=chainT(n);t.m.b.push([t.loc[n],t.loc[1],1,0]);t.ring=true;return t}

/* spiro[x.y] and von Baeyer bi/tri/tetra/penta/hexacyclo[a.b.c.d^p,q…] skeletons, numbered by the IUPAC rules */
const NRING={spiro:2,bicyclo:2,tricyclo:3,tetracyclo:4,pentacyclo:5,hexacyclo:6};
function polyParse(kind,desc){
  const tk=desc.split('.');
  if(kind==='spiro'){if(tk.length!==2||!tk.every(x=>/^\d+$/.test(x)))throw NP();let[x,y]=tk.map(Number);if(x>y)[x,y]=[y,x];
    if(x<2)throw CH('Each ring of a spiro compound needs at least 3 atoms — spiro[2.2]… is the smallest');return{kind,x,y,n:x+y+1}}
  if(tk.length!==NRING[kind]+1)throw CH('"'+kind+'" needs '+(NRING[kind]+1)+' numbers in the brackets');
  const main=tk.slice(0,3).map(x=>{if(!/^\d+$/.test(x))throw NP();return+x}),sec=tk.slice(3).map(x=>{const q=/^(\d+)\^(\d+),(\d+)$/.exec(x);if(!q)throw CH('Secondary bridges are written like 1^3,7 (atoms in the bridge, then the two bridgeheads)');return{d:+q[1],p:+q[2],q:+q[3]}});
  const n=main[0]+main[1]+main[2]+2+sec.reduce((s,b)=>s+b.d,0);return{kind,main,sec,n}}
function polyT(pl){
  const m=Mol(),loc={},order=[],n=pl.n;
  for(let i=1;i<=n;i++){loc[i]=i-1;order.push(String(i));m.a.push({el:'C',ar:false,q:0})}
  const bd=(i,j)=>{if(i<1||j<1||i>n||j>n)throw CH('A bridgehead number is outside 1–'+n);m.b.push([i-1,j-1,1,0])};
  if(pl.kind==='spiro'){const x=pl.x,y=pl.y;
    for(let i=1;i<x;i++)bd(i,i+1);bd(x,x+1);bd(x+1,1);bd(x+1,x+2);for(let i=x+2;i<x+y+1;i++)bd(i,i+1);bd(x+y+1,x+1)}
  else{const[a,b,c]=pl.main;
    for(let i=1;i<=a+1;i++)bd(i,i+1);                    // main ring, first bridge: C1 → second main bridgehead (a+2)
    for(let i=a+2;i<a+b+2;i++)bd(i,i+1);bd(a+b+2,1);     // second bridge back to C1
    const s=a+b+3;if(c===0)bd(1,a+2);else{bd(1,s);for(let k=0;k<c-1;k++)bd(s+k,s+k+1);bd(s+c-1,a+2)}   // main bridge, numbered from the end nearer C1
    let nx=a+b+c+3;pl.sec.forEach(({d,p,q})=>{if(d===0){bd(p,q);return}const hi=Math.max(p,q),lo=Math.min(p,q);
      bd(hi,nx);for(let k=0;k<d-1;k++)bd(nx+k,nx+k+1);bd(nx+d-1,lo);nx+=d})}   // secondary bridges: numbered from the end nearer the higher bridgehead
  return{m,loc,order,ring:true,n}}

/* ───────────────────────── vocabulary ───────────────────────── */
const MUL={di:2,tri:3,tetra:4,penta:5,hexa:6,bis:2,tris:3,tetrakis:4};
const SNUM={meth:1,eth:2,prop:3,but:4,pent:5,hex:6,hept:7,oct:8,non:9,dec:10,undec:11,dodec:12,tridec:13,tetradec:14,pentadec:15,hexadec:16,heptadec:17,octadec:18,nonadec:19,eicos:20,icos:20};
const STEMRE='undec|dodec|tridec|tetradec|pentadec|hexadec|heptadec|octadec|nonadec|eicos|icos|hept|hex|meth|eth|prop|but|pent|oct|non|dec';
const SIMPLE={fluoro:'F',chloro:'Cl',bromo:'Br',iodo:'I',hydroxy:'O',amino:'N',mercapto:'S',sulfanyl:'S',nitro:'[N+](=O)[O-]',nitroso:'N=O',cyano:'C#N',
  carboxy:'C(=O)O',formyl:'C=O',methanoyl:'C=O',carbamoyl:'C(=O)N',sulfo:'S(=O)(=O)O',vinyl:'C=C',allyl:'CC=C',phenyl:'c1ccccc1',benzyl:'Cc1ccccc1',
  benzoyl:'C(=O)c1ccccc1',anilino:'Nc1ccccc1',oxo:['O',2],keto:['O',2],imino:['N',2],isocyano:'[N+]#[C-]',carbylamino:'[N+]#[C-]',chloroformyl:'C(=O)Cl',fluoroformyl:'C(=O)F',bromoformyl:'C(=O)Br',iodoformyl:'C(=O)I',
  azido:'N=[N+]=[N-]',diazo:['[N+]=[N-]',2]};
const ALIAS={isopropyl:'propan-2-yl',isobutyl:'2-methylpropyl',secbutyl:'butan-2-yl',tertbutyl:'2-methylpropan-2-yl',neopentyl:'2,2-dimethylpropyl',
  isopentyl:'3-methylbutyl',isoamyl:'3-methylbutyl',isopropenyl:'prop-1-en-2-yl',tertpentyl:'2-methylbutan-2-yl'};
const RPAR={phenol:{s:['ol']},aniline:{s:['amine']},toluene:{p:'methyl'},anisole:{p:'methoxy'},styrene:{p:'ethenyl'},benzoicacid:{s:['carboxylicacid']},
  benzaldehyde:{s:['carbaldehyde']},benzamide:{s:['carboxamide']},benzonitrile:{s:['carbonitrile']},phenyl:{s:['yl']},benzoyl:{s:['carbonyl']}};
const FUS={naphthalen:'naphthalene',azulen:'azulene',anthracen:'anthracene',phenanthren:'phenanthrene',pyridin:'pyridine',thiophen:'thiophene',pyrrol:'pyrrole',furan:'furan'};
const EMB=new Set(['al','oicacid','amide','nitrile','oyl']);              // suffixes whose carbon is part of the chain itself
const SUFG={ol:()=>G('O'),thiol:()=>G('S'),amine:()=>G('N',1,{N:0}),one:()=>G('O',2),carboxylicacid:()=>G('C(=O)O',1,{O:2}),carbaldehyde:()=>G('C=O'),
  carboxamide:()=>G('C(=O)N',1,{N:2}),carbonitrile:()=>G('C#N'),sulfonicacid:()=>G('S(=O)(=O)O'),carbonyl:()=>G('C=O',1,{t:1})};
const HAL={fluoride:'F',chloride:'Cl',bromide:'Br',iodide:'I'};

/* ───────────────────────── parent-name reader (stem/ring + ene/yne + suffixes) ───────────────────────── */
function parseParent(t){
  const P={kind:null,unsat:[],suf:[],pre0:[]};let pos=0,m;
  if(RPAR[t]){const d=RPAR[t];P.kind='arom';P.tm='benzene';if(d.s)P.suf.push({t:d.s[0],locs:['1'],m:0,imp:1});if(d.p)P.pre0.push({loc:'1',name:d.p});return fin(P)}
  if((m=/^benzene/.exec(t))){P.kind='arom';P.tm='benzene';pos=7}
  else if((m=/^(naphthalen|azulen|anthracen|phenanthren|pyridin|thiophen|pyrrol)e?/.exec(t)||/^(furan)/.exec(t))){P.kind='arom';P.tm=FUS[m[1]];pos=m[0].length}
  else if((m=/^(spiro|bicyclo|tricyclo|tetracyclo|pentacyclo|hexacyclo)\[([0-9.^,]+)\]/.exec(t))){
    P.kind='poly';P.poly=polyParse(m[1],m[2]);pos=m[0].length;
    const sm=new RegExp('^('+STEMRE+')').exec(t.slice(pos));if(!sm)return null;pos+=sm[0].length;
    if(SNUM[sm[1]]!==P.poly.n)throw CH('['+m[2]+'] describes '+P.poly.n+' carbons, so the name needs "'+Object.keys(SNUM).find(k=>SNUM[k]===P.poly.n&&k!=='icos')+'" as its stem')}
  else if((m=new RegExp('^cyclo('+STEMRE+')').exec(t))){P.kind='ring';P.n=SNUM[m[1]];pos=m[0].length}
  else if((m=new RegExp('^('+STEMRE+')').exec(t))){P.kind='chain';P.n=SNUM[m[1]];pos=m[0].length}
  else return null;
  if(P.kind!=='arom'&&/^a(?=-?\d|di|tri|tetra)/.test(t.slice(pos)))pos++;
  if(P.kind!=='arom'){const U=/^-?(?:(\d+(?:,\d+)*)-)?(di|tri|tetra)?(an|en|yn)e?/;let u;
    while((u=U.exec(t.slice(pos)))){P.unsat.push({t:u[3],locs:u[1]?u[1].split(','):[],m:u[2]?MUL[u[2]]:0});pos+=u[0].length}}
  const S=/^-?(?:(\d+[a-z]?(?:,\d+[a-z]?)*)-)?(di|tri|tetra|penta|hexa)?(carboxylicacid|carbaldehyde|carboxamide|carbonitrile|sulfonicacid|carbonyl|oicacid|thiol|amine|amide|nitrile|oyl|one|ol|al|yl)/;let s;
  while((s=S.exec(t.slice(pos)))){P.suf.push({t:s[3],locs:s[1]?s[1].split(','):[],m:s[2]?MUL[s[2]]:0});pos+=s[0].length}
  if(pos!==t.length)return null;
  if(P.kind==='arom'&&P.suf.some(x=>EMB.has(x.t)))return null;
  return fin(P)}
function fin(P){P.hasN=P.suf.some(s=>s.t==='amine'||s.t==='amide'||s.t==='carboxamide');return P}

/* ───────────────────────── substituent / group reader ───────────────────────── */
const GC=new Map();
function grp(s){if(GC.has(s)){const v=GC.get(s);if(!v)throw NP();return v}
  let v=null;try{v=grp0(s)}catch(e){v=null}GC.set(s,v);if(!v)throw NP();return v}
const alkylize=x=>new RegExp('('+STEMRE.split('|').filter(w=>!/dec$|eicos|icos/.test(w)||w==='dec').join('|')+')$').test(x)?x+'yl':x==='phen'?'phenyl':x;
function glist(str){                                    // "dimethyl" → 2 groups, "ethylmethyl" → 2 groups
  if(!str)return[];
  try{return[grp(str)]}catch(e){}
  let m=/^(di|tri)(.+)$/.exec(str);if(m){try{const g=grp(m[2]);return Array(MUL[m[1]]).fill(g)}catch(e){}}
  for(let k=1;k<str.length;k++){let a,b;try{a=glist(str.slice(0,k));b=glist(str.slice(k))}catch(e){continue}if(a&&b&&a.length&&b.length)return a.concat(b)}
  throw NP()}
function grp0(s){
  while(s[0]==='('&&s[s.length-1]===')'&&bal(s.slice(1,-1))||s[0]==='['&&s[s.length-1]===']'&&bal(s.slice(1,-1)))s=s.slice(1,-1);
  if(!s)throw NP();
  if(s==='epoxy')return{m:smi('O'),att:0,ord:1,bridge:true};
  if(ALIAS[s])return grp(ALIAS[s]);
  let cb;if((cb=/^carbo?(meth|eth|prop|but|pent|hex|hept|oct|non|dec|phen)oxy$/.exec(s)))return grp(cb[1]+'oxycarbonyl');   // carbomethoxy, carbethoxy …
  if(SIMPLE[s]!==undefined){const d=SIMPLE[s];return Array.isArray(d)?G(d[0],d[1]):G(d)}
  let m;
  const tryc=f=>{try{return f()}catch(e){if(!e.np)throw e;return null}};
  if((m=/^(.+)oxy$/.exec(s))){const r=tryc(()=>{const g=grp(alkylize(m[1])),r=G('O');attach(r.m,0,g);return r});if(r)return r}
  if((m=/^(.+)carbonyl$/.exec(s))){const r=tryc(()=>{const g=grp(m[1]),r=G('C=O');attach(r.m,0,g);return r});if(r)return r}
  if((m=/^(.+)amino$/.exec(s))){const r=tryc(()=>{const gs=glist(m[1]);if(gs.length>2)throw NP();const r=G('N');gs.forEach(g=>attach(r.m,0,g));return r});if(r)return r}
  if((m=/^(.+)(?:sulfanyl|thio)$/.exec(s))){const r=tryc(()=>{const g=grp(m[1]),r=G('S');attach(r.m,0,g);return r});if(r)return r}
  const c=parseCompound(s);if(c.att==null)throw NP();
  const a=c.m.a[c.att];if(a.res)a.res-=1;                // the free valence becomes a real bond once the group is attached
  return{m:c.m,att:c.att,ord:1}}
function bal(s){let d=0;for(const c of s){if(c==='('||c==='[')d++;else if(c===')'||c===']'){if(--d<0)return false}}return d===0}

/* ───────────────────────── prefix list reader ───────────────────────── */
const OMP={o:2,ortho:2,m:3,meta:3,p:4,para:4};
function parsePrefixList(p,P){
  if(!p)return{items:[],pl:null};
  const memo=new Map(),LOC=/^(?:(\d+[a-z]?(?:,\d+[a-z]?)*)-|([nν](?:,[nν])*)-|(o|m|p|ortho|meta|para)-)/;
  function seg(pos){
    if(pos>=p.length)return{items:[],pl:null};
    if(memo.has(pos))return memo.get(pos);
    let res=null;
    if(p[pos]==='-'||p[pos]===','){res=seg(pos+1);memo.set(pos,res);return res}
    let locs=[],pos2=pos;const lm=LOC.exec(p.slice(pos));
    if(lm){if(lm[1])locs=lm[1].split(',');else if(lm[2]){if(!P.hasN){memo.set(pos,null);return null}locs=lm[2].split(',').map(()=>'N')}else locs=[lm[3]];pos2=pos+lm[0].length}
    if(lm&&pos2>=p.length){res={items:[],pl:locs};memo.set(pos,res);return res}
    const mm=/^(bis|tris|tetrakis|di|tri|tetra|penta|hexa)/.exec(p.slice(pos2));
    const variants=mm?[[pos2+mm[0].length,MUL[mm[0]]],[pos2,0]]:[[pos2,0]];
    if(locs.length>1&&mm)variants.reverse().reverse();
    else if(!(locs.length>1)&&false)variants.reverse();
    const vs=(locs.length>1&&mm)?variants:variants.slice().reverse().reverse();
    outer:for(const[pos3,mult]of vs){
      if(pos3>=p.length)continue;
      const cands=[];
      if(p[pos3]==='('||p[pos3]==='['){let d=0,e=-1;for(let k=pos3;k<p.length;k++){if('(['.includes(p[k]))d++;else if(')]'.includes(p[k])&&--d===0){e=k;break}}if(e>0)cands.push([p.slice(pos3+1,e),e+1])}
      for(let e=p.length;e>pos3;e--){const tx=p.slice(pos3,e);if(tx.endsWith('-')||!bal(tx))continue;cands.push([tx,e])}
      for(const[tx,e]of cands){let g;try{g=grp(tx)}catch(x){continue}
        const rest=seg(e);if(!rest)continue;
        res={items:[{locs,m:mult,g}].concat(rest.items),pl:rest.pl};break outer}}
    memo.set(pos,res);return res}
  return seg(0)}

/* ───────────────────────── builder ───────────────────────── */
function buildMol(P,items){
  const T=P.kind==='chain'?chainT(P.n):P.kind==='ring'?ringT(P.n):P.kind==='poly'?polyT(P.poly):tmpl(P.tm);
  const m=T.m,L=T.loc,N=T.order.length,isRing=!!(T.ring||T.ar),info={N:[],acidO:[],att:null};
  const pos=l=>{const k=String(l);if(!(k in L))throw CH('Position '+l+' does not exist on this '+(T.ar?'ring system':isRing?'ring':'chain'));return L[k]};
  const addG=(i,g,lab)=>{if(capm(m,i)<g.ord-1e-9)throw CH('Position '+lab+' cannot take another group (its bonds are already full)');return attach(m,i,g)};
  const rec=(g,off)=>{if(g.N!==undefined)info.N.push(off+g.N);if(g.O!==undefined)info.acidO.push(off+g.O);if(g.t===1)info.att=off+g.att};
  const bondIdx=(i,j)=>m.b.findIndex(b=>b&&((b[0]===i&&b[1]===j)||(b[0]===j&&b[1]===i)));
  const multi=(p,o)=>{if(T.ar)throw CH('ene/yne endings cannot be used on an aromatic ring system');
    const i=pos(p),np=String(+p+1);let j;
    if(np in L)j=L[np];else if(isRing&&String(p)===String(N))j=L['1'];else throw CH('A multiple bond at '+p+' needs a carbon '+np);
    const bi=bondIdx(i,j);if(bi<0)throw CH('Positions '+p+' and '+np+' are not bonded in this ring system');if(m.b[bi][2]>1)throw CH('Position '+p+' already has a multiple bond');
    if(capm(m,i)<o-1||capm(m,j)<o-1)throw CH('No room for a multiple bond at position '+p);m.b[bi][2]=o};
  const bridgeAt=(a,b)=>{const i=pos(a),j=pos(b);if(i===j)throw CH('An epoxy bridge needs two different positions');
    if(capm(m,i)<1||capm(m,j)<1)throw CH('Positions '+a+' and '+b+' cannot both carry an epoxy bridge');
    const k=m.a.length;m.a.push({el:'O',ar:false,q:0});m.b.push([i,k,1,0],[j,k,1,0])};
  const hasMulti=i=>m.b.some(b=>b&&(b[0]===i||b[1]===i)&&b[2]>1);
  const cnt=(s,what)=>{const c=s.locs.length||s.m||1;if(s.locs.length&&s.m>1&&s.locs.length!==s.m)throw CH('The number of locants does not match the multiplier in "'+what+'"');return c};
  const place=(mk,count,start,what)=>{for(let k=0;k<count;k++){let ok=false;
      for(let t=0;t<N&&!ok;t++){const i=L[T.order[(start+t)%N]],g=mk();if(capm(m,i)>=g.ord-1e-9){rec(g,attach(m,i,g));ok=true}}
      if(!ok)throw CH('There is no room for all those '+what+' groups')}};
  // 1 — double/triple bonds with explicit locants
  P.unsat.forEach(u=>{if(u.t==='an')return;cnt(u,u.t);u.locs.forEach(p=>multi(p,u.t==='en'?2:3))});
  // 2 — chain-embedded suffixes (al, oic acid, amide, nitrile, oyl), 'yl', and any suffix with locants
  const suf=P.suf.map(s=>Object.assign({},s,{locs:s.locs.slice()}));
  suf.forEach(s=>{
    if(EMB.has(s.t)){if(isRing)throw SOFT('"'+s.t+'" endings are for chains — on a ring use carbaldehyde / carboxylic acid / carboxamide / carbonitrile');
      const c=cnt(s,s.t),ls=s.locs.length?s.locs:(c===1?['1']:c===2?['1',String(N)]:null);if(!ls)throw SOFT('Give locants for "'+s.t+'"');
      ls.forEach(p=>{if(!isRing&&String(p)!=='1'&&String(p)!==String(N))throw CH('"'+s.t+'" can only sit at an end of the chain (1 or '+N+')');const i=pos(p);
        if(s.t==='al')addG(i,G('O',2),p);
        else if(s.t==='oicacid'){addG(i,G('O',2),p);const o=addG(i,G('O'),p);info.acidO.push(o)}
        else if(s.t==='amide'){addG(i,G('O',2),p);info.N.push(addG(i,G('N'),p))}
        else if(s.t==='nitrile')addG(i,G('N',3),p);
        else if(s.t==='oyl'){addG(i,G('O',2),p);m.a[i].res=(m.a[i].res||0)+1;info.att=i}});s.done=1}
    else if(s.t==='yl'){const p=s.locs[0]||'1',i=pos(p);m.a[i].res=(m.a[i].res||0)+1;if(capm(m,i)<0)throw CH('Position '+p+' has no free valence for -yl');info.att=i;s.done=1}
    else if(s.locs.length){cnt(s,s.t);s.locs.forEach(p=>{const g=SUFG[s.t](),i=pos(p);rec(g,addG(i,g,p))});s.done=1}});
  // 3 — prefixes with explicit locants (N-locants wait until the nitrogen exists)
  const deferN=[],noLoc=[];
  items.forEach(it=>{
    if(it.g.bridge){if(!it.locs.length){noLoc.push(it);return}
      if(it.locs.length%2||(it.m>1&&it.locs.length!==2*it.m))throw CH('Epoxy needs a pair of locants for each bridge, e.g. 2,3-epoxybutane');
      for(let k=0;k<it.locs.length;k+=2)bridgeAt(it.locs[k],it.locs[k+1]);return}
    if(it.locs.length){if(it.m>1&&it.locs.length!==it.m)throw CH('The number of locants does not match the multiplier');
      it.locs.forEach(l=>{if(l==='N'){deferN.push(it.g);return}addG(pos(l),it.g,l)})}
    else noLoc.push(it)});
  // 4 — remaining suffixes without locants
  suf.filter(s=>!s.done).forEach(s=>{if(!SUFG[s.t])throw NP();
    const start=(s.t==='one'&&!isRing&&N>=3)?1:0;place(SUFG[s.t],s.m||1,start,s.t)});
  // 5 — ene / yne without locants: first positions that still have valence and do not touch another multiple bond
  P.unsat.forEach(u=>{if(u.t==='an'||u.locs.length)return;const o=u.t==='en'?2:3,c=u.m||1;
    for(let k=0;k<c;k++){let ok=false;
      for(let t=1;t<=(isRing?N:N-1)&&!ok;t++){const i=pos(t),j=(String(t+1) in L)?L[t+1]:L['1'];if(hasMulti(i)||hasMulti(j))continue;try{multi(t,o);ok=true}catch(e){}}
      if(!ok)throw CH('No room for a '+(o===2?'double':'triple')+' bond')}});
  // 6 — prefixes without locants
  noLoc.forEach(it=>{
    if(it.g.bridge){for(let c=0;c<(it.m||1);c++){let ok=false;
      for(let t=0;t<N&&!ok;t++){if(!isRing&&t===N-1)break;const a=T.order[t],b=T.order[(t+1)%N];
        if(bondIdx(L[a],L[b])>=0&&capm(m,L[a])>=1&&capm(m,L[b])>=1){bridgeAt(a,b);ok=true}}
      if(!ok)throw CH('There is no room for an epoxy bridge')}return}
    place(()=>it.g,it.m||1,(!isRing&&N>=3)?1:0,'substituent')});
  // 7 — N-substituents
  deferN.forEach(g=>{if(!info.N.length)throw CH('N- locants need an amine or amide in the name');const ni=info.N[0];
    if(capm(m,ni)<1)throw CH('The nitrogen cannot take another group');attach(m,ni,g)});
  m.a.forEach((a,i)=>{if(a.dead||a.q||(a.ar&&a.el!=='C'))return;
    const rd=a.ar?m.b.filter(b=>b&&b[3]&&(b[0]===i||b[1]===i)).length:0;
    if(used(m,i)>(MAXV[a.el]||4)+(rd===3?.5:0)+1e-9)throw CH('Too many bonds on a '+a.el+' atom — check the locants')});
  return{m,att:info.att,N:info.N,acidO:info.acidO}}

function parseCompound(s){
  let first=null;
  for(let i=0;i<s.length;i++){
    const ps=s.slice(i);if(/[()]/.test(ps)||(/[\[\]]/.test(ps)&&!/^(?:spiro|(?:bi|tri|tetra|penta|hexa)cyclo)\[/.test(ps)))continue;
    const pre=s.slice(0,i);if(!bal(pre))continue;
    let P;try{P=parseParent(ps)}catch(e){if(!e.np&&!first)first=e;continue}
    if(!P)continue;
    const R=parsePrefixList(pre,P);if(!R)continue;
    try{return build2(P,R)}catch(e){if(e.np)continue;if(!first)first=e}}
  if(first)throw first;throw NP()}
function build2(P,R){
  P=Object.assign({},P,{unsat:P.unsat.map(u=>Object.assign({},u,{locs:u.locs.slice()})),suf:P.suf.map(u=>Object.assign({},u,{locs:u.locs.slice()}))});
  const items=P.pre0.map(x=>({locs:[x.loc],m:0,g:grp(x.name)}));
  const blocks=P.unsat.filter(u=>u.t!=='an').concat(P.suf);
  if(R.pl){const tgt=P.suf.find(s=>s.imp)||blocks.find(b=>!b.locs.length);
    if(!tgt)throw NP();
    const pl=R.pl.map(x=>OMP[x]?String(OMP[x]):x);if(pl.some(x=>!/^\d/.test(x)))throw NP();tgt.locs=pl}
  R.items.forEach(it=>{
    if(it.locs.length&&it.locs.every(x=>OMP[x])){
      const anchored=P.pre0.length||P.suf.some(s=>s.imp);
      if(anchored)it=Object.assign({},it,{locs:it.locs.map(x=>String(OMP[x]))});
      else if(it.m>1&&it.locs.length===1)it=Object.assign({},it,{locs:['1',String(OMP[it.locs[0]])]});
      else throw NP()}
    items.push(it)});
  if(P.kind==='arom'&&P.pre0.length===0&&!P.suf.length&&false)throw NP();
  return buildMol(P,items)}

/* ───────────────────────── whole-name layer ───────────────────────── */
const single=str=>{const r=parseCompound(str);if(r.att!=null)throw NP();return r};
function esterOf(words){
  let k=0;while(k<words.length-1&&/yl$/.test(words[k]))k++;if(!k)throw NP();
  const gs=[];words.slice(0,k).forEach(w=>glist(w).forEach(g=>gs.push(g)));
  const acid=words.slice(k).join('').replace(/oate$/,'oicacid').replace(/carboxylate$/,'carboxylicacid');
  const r=single(acid);if(!r.acidO.length||gs.length>r.acidO.length)throw NP();
  gs.forEach((g,i)=>attach(r.m,r.acidO[i],g));return r.m}
function joinAnh(r1,r2){const off=merge(r1.m,r2.m),o2=r2.acidO[0]+off,c2=nbrs(r1.m,o2)[0];kill(r1.m,o2);r1.m.b.push([r1.acidO[0],c2,1,0]);return r1.m}
function anhydrideOf(words){
  const ws=words.filter(w=>w!=='acid');if(!ws.length||ws.length>2)throw NP();
  const acid=w=>single(w.replace(/acid$/,'')+'acid');
  if(ws.length===2){const a=acid(ws[0]),b=acid(ws[1]);if(a.acidO.length!==1||b.acidO.length!==1)throw NP();return joinAnh(a,b)}
  const a=acid(ws[0]);
  if(a.acidO.length===1)return joinAnh(a,acid(ws[0]));
  if(a.acidO.length===2){const[o1,o2]=a.acidO,c1=nbrs(a.m,o1)[0],c2=nbrs(a.m,o2)[0];
    // shortest path between the two carboxyl carbons decides the ring size (+1 for the bridging O)
    const dist={[c1]:1},q=[c1];while(q.length){const u=q.shift();nbrs(a.m,u).forEach(v=>{if(v!==o1&&v!==o2&&!(v in dist)){dist[v]=dist[u]+1;q.push(v)}})}
    const ring=(dist[c2]||99)+1;if(ring<5||ring>7)throw CH('A cyclic anhydride needs the two acid groups close together (5–7 membered ring)');
    kill(a.m,o2);a.m.b.push([o1,c2,1,0]);return a.m}
  throw NP()}
function cyanideOf(body,kind){
  const g=grp(body),r=Mol(),off=merge(r,g.m),x=r.a.length;
  if(kind==='cyanide')r.a.push({el:'C',ar:false,q:0},{el:'N',ar:false,q:0});else r.a.push({el:'N',ar:false,q:1},{el:'C',ar:false,q:-1});
  r.b.push([off+g.att,x,1,0],[x,x+1,3,0]);return r}
function acylHalide(words){
  const body=words.slice(0,-1).join('');if(!/(oyl|carbonyl)$/.test(body))throw NP();
  const g=grp(body),r=Mol(),off=merge(r,g.m),x=r.a.length;r.a.push({el:HAL[words[words.length-1]],ar:false,q:0});r.b.push([off+g.att,x,1,0]);return r}

/* ── retained names & typo repair ── */
const RET=[
 [/\b(?:pyrocatechol|catechol)\b/g,'benzene-1,2-diol'],[/\bresorcinol\b/g,'benzene-1,3-diol'],[/\b(?:hydroquinone|hydroquinol|quinol)\b/g,'benzene-1,4-diol'],
 [/\bpyrogallol\b/g,'benzene-1,2,3-triol'],[/\bphloroglucinol\b/g,'benzene-1,3,5-triol'],[/\bpicric acid\b/g,'2,4,6-trinitrophenol'],
 [/\b(?:o|ortho)-?cresol\b/g,'2-methylphenol'],[/\b(?:m|meta)-?cresol\b/g,'3-methylphenol'],[/\b(?:p|para)-?cresol\b/g,'4-methylphenol'],[/\bcresol\b/g,'2-methylphenol'],
 [/\b(?:o|ortho)-?xylene\b/g,'1,2-dimethylbenzene'],[/\b(?:m|meta)-?xylene\b/g,'1,3-dimethylbenzene'],[/\b(?:p|para)-?xylene\b/g,'1,4-dimethylbenzene'],[/\bxylene\b/g,'1,2-dimethylbenzene'],
 [/\b(?:α|alpha)-?naphthol\b/g,'naphthalen-1-ol'],[/\b(?:β|beta)-?naphthol\b/g,'naphthalen-2-ol'],[/(\d[a-z]?)-?naphthol\b/g,'naphthalen-$1-ol'],[/\bnaphthol\b/g,'naphthalen-1-ol'],
 [/(\d[a-z]?)-?naphthylamine\b/g,'naphthalen-$1-amine'],[/\bnaphthylamine\b/g,'naphthalen-1-amine'],
 [/(\d[a-z]?)-?naphthoic acid\b/g,'naphthalene-$1-carboxylic acid'],[/\bnaphthoic acid\b/g,'naphthalene-1-carboxylic acid'],
 [/(\d[a-z]?)-?naphthyl\b/g,'naphthalen-$1-yl'],[/\bnaphthyl\b/g,'naphthalen-1-yl'],
 [/\bacetanilide\b/g,'ν-phenylethanamide'],[/\bcinnamaldehyde\b/g,'3-phenylprop-2-enal'],[/\bcinnamyl alcohol\b/g,'3-phenylprop-2-en-1-ol'],[/\bcinnamic acid\b/g,'3-phenylprop-2-enoic acid'],
 [/\bisophthalic\b/g,'benzene-1,3-dicarboxylic'],[/\bterephthalic\b/g,'benzene-1,4-dicarboxylic'],[/\bphthalic\b/g,'benzene-1,2-dicarboxylic'],
 [/\bacetylsalicylic acid\b/g,'2-ethanoyloxybenzoic acid'],[/\bsalicylic acid\b/g,'2-hydroxybenzoic acid'],
 [/\bacetophenone\b/g,'1-phenylethanone'],[/\bbenzophenone\b/g,'diphenylmethanone'],[/\bacetone\b/g,'propanone'],[/\bethylene oxide\b/g,'1,2-epoxyethane'],[/\bpropylene oxide\b/g,'1,2-epoxypropane'],[/oxirane/g,'-1,2-epoxyethane'],[/\baspirin\b/g,'2-ethanoyloxybenzoic acid'],[/\bbiphenyl\b/g,'phenylbenzene'],
 [/\bnorbornane\b/g,'bicyclo[2.2.1]heptane'],[/\bnorbornene\b/g,'bicyclo[2.2.1]hept-2-ene'],[/\bdecalin\b/g,'bicyclo[4.4.0]decane'],[/\bspiropentane\b/g,'spiro[2.2]pentane'],
 [/\badamantane\b/g,'tricyclo[3.3.1.1^3,7]decane'],[/\bcubane\b/g,'pentacyclo[4.2.0.0^2,5.0^3,8.0^4,7]octane'],[/\bcamphor\b/g,'1,7,7-trimethylbicyclo[2.2.1]heptan-2-one'],
 [/\bformic\b/g,'methanoic'],[/\bacetic\b/g,'ethanoic'],[/\bpropionic\b/g,'propanoic'],[/\bbutyric\b/g,'butanoic'],[/\bvaleric\b/g,'pentanoic'],[/\bcaproic\b/g,'hexanoic'],
 [/\bacrylic\b/g,'prop-2-enoic'],[/\bcrotonic\b/g,'but-2-enoic'],[/\boxalic\b/g,'ethanedioic'],[/\bmalonic\b/g,'propanedioic'],[/\bsuccinic\b/g,'butanedioic'],
 [/\bglutaric\b/g,'pentanedioic'],[/\badipic\b/g,'hexanedioic'],[/\b(?:maleic|fumaric)\b/g,'but-2-enedioic'],[/\blactic\b/g,'2-hydroxypropanoic'],[/\bpyruvic\b/g,'2-oxopropanoic'],
 [/\bformate\b/g,'methanoate'],[/\bacetate\b/g,'ethanoate'],[/\bpropionate\b/g,'propanoate'],[/\bbutyrate\b/g,'butanoate'],
 [/acetyl(?!ene)/g,'ethanoyl'],[/propionyl/g,'propanoyl'],[/butyryl/g,'butanoyl'],[/acetoxy/g,'ethanoyloxy'],[/acetamido/g,'ethanoylamino'],
 [/\bacetamide\b/g,'ethanamide'],[/\bformamide\b/g,'methanamide'],[/\bacetonitrile\b/g,'ethanenitrile'],[/\bacetaldehyde\b/g,'ethanal'],[/\bformaldehyde\b/g,'methanal'],
 [/\bpropionaldehyde\b/g,'propanal'],[/\bbutyraldehyde\b/g,'butanal']];
const KNOWN='catechol resorcinol hydroquinone quinol pyrogallol phloroglucinol cresol naphthol naphthylamine acetanilide cinnamic cinnamaldehyde isophthalic terephthalic phthalic salicylic benzoic anhydride azulene naphthalene anthracene phenanthrene benzene pyridine thiophene pyrrole furan acetophenone benzophenone succinic maleic fumaric glutaric adipic oxalic malonic acetic formic propionic butyric acetone xylene toluene aniline phenol'.split(' ');
function lev(a,b){const d=Array.from({length:a.length+1},(_,i)=>[i]);for(let j=1;j<=b.length;j++)d[0][j]=j;
  for(let i=1;i<=a.length;i++)for(let j=1;j<=b.length;j++){d[i][j]=Math.min(d[i-1][j]+1,d[i][j-1]+1,d[i-1][j-1]+(a[i-1]===b[j-1]?0:1));if(i>1&&j>1&&a[i-1]===b[j-2]&&a[i-2]===b[j-1])d[i][j]=Math.min(d[i][j],d[i-2][j-2]+1)}return d[a.length][b.length]}
const SAFE='phenyl benzyl benzoyl phenoxy anilino aniline toluene'.split(' ');   // real words that sit one letter from a retained name (phenyl ↔ phenol)
function repair(s){return s.replace(/[a-z]{5,}/g,w=>{if(KNOWN.includes(w)||SAFE.includes(w))return w;const mx=w.length>=8?2:1;let best=null,bd=9;
  KNOWN.forEach(k=>{if(k[0]!==w[0]||Math.abs(k.length-w.length)>mx)return;const d=lev(w,k);if(d<=mx&&d<bd){bd=d;best=k}});return best||w})}
function clean(x){
  x=x.replace(/[–—−]/g,'-').replace(/[’′`]/g,"'");
  const SUP={'⁰':0,'¹':1,'²':2,'³':3,'⁴':4,'⁵':5,'⁶':6,'⁷':7,'⁸':8,'⁹':9};
  x=x.replace(/[⁰¹²³⁴⁵⁶⁷⁸⁹]+(?:,[⁰¹²³⁴⁵⁶⁷⁸⁹]+)*/g,r=>'^'+r.replace(/[⁰¹²³⁴⁵⁶⁷⁸⁹]/g,c=>SUP[c]));
  x=x.replace(/(spiro|(?:bi|tri|tetra|penta|hexa)cyclo)\s*\[([^\]]*)\]/gi,(m,k,i)=>k.toLowerCase()+'['+i.replace(/\s/g,'').replace(/\^?[{(]/g,'^').replace(/[})]/g,'')+']');
  x=x.replace(/(^|[\s(\[,\-])N(?=[,\-\s'])/g,'$1ν');                 // N-locants (capital N) survive lower-casing as ν
  x=x.toLowerCase().replace(/sulph/g,'sulf').replace(/(^|[-,(\s])\d+h-/g,'$1');                // indicated hydrogen (1H-pyrrole) is implied by the ring
  x=x.replace(/\(\s*(?:\d+[a-z]?\s*,?\s*)*[ezrs](?:\s*,\s*(?:\d+[a-z]?\s*)?[ezrs])*\s*\)\s*-?/g,'').replace(/\((?:±|\+|-|dl|rac)\)-?/g,'')
   .replace(/\b(?:cis|trans)-/g,'').replace(/^(?:e|z|r|s|d|l|dl)-(?=[a-z0-9])/,'');
  x=x.replace(/\s*([-,()\[\]])\s*/g,'$1').replace(/\s+/g,' ').trim();
  x=x.replace(/\b(sec|tert|iso|neo)-(?=[a-z])/g,'$1').replace(/\bt-butyl/g,'tertbutyl');
  return x}
function expand(s){RET.forEach(([re,to])=>{s=s.replace(re,to)});return s.replace(/\s*([-,()\[\]])\s*/g,'$1').replace(/-{2,}/g,'-').replace(/^-/,'')}

function addonFromName(raw){
  let s=clean(String(raw||''));if(!/[a-z]{3}/.test(s))return null;
  s=expand(repair(s));
  const words=s.split(' ').filter(Boolean),last=words[words.length-1];let m,disp;
  try{
    if(last==='anhydride'){m=anhydrideOf(words.slice(0,-1));disp=words.join(' ')}
    else if(words.length>1&&/^(?:cyanide|isocyanide|isonitrile)$/.test(last)){m=cyanideOf(words.slice(0,-1).join(''),last==='cyanide'?'cyanide':'iso');disp=words.join(' ')}
    else if(words.length>1&&HAL[last]){m=acylHalide(words);disp=words.join(' ')}
    else if(words.length>1&&/(?:oate|carboxylate)$/.test(last)){m=esterOf(words);disp=words.join(' ')}
    else{const j=words.join('');disp=j.replace(/(ic)acid$/,'$1 acid');
      try{m=single(j).m}catch(e){if(!e.np)throw e;
        const cn=/^(.+?)(isocyanide|isonitrile|cyanide)$/.exec(j);
        if(cn)m=cyanideOf(cn[1],cn[2]==='cyanide'?'cyanide':'iso');
        else{const am=/^(.+)amine$/.exec(j);if(!am)throw e;const gs=glist(am[1]);if(!gs.length||gs.length>3)throw e;
          m=Mol();m.a.push({el:'N',ar:false,q:0});gs.forEach(g=>attach(m,0,g))}}}
    if(!m.a.some(a=>!a.dead))return null;
    const smiles=toSmiles(m).replace(/ν/g,'');
    return{src:'smiles:'+smiles,name:disp.replace(/ν/g,'N'),note:'Interpreted as '+disp.replace(/ν/g,'N')}
  }catch(e){if(e.np)return null;return{err:e.message,hard:!!e.hard}}}

/* ───────────────────────── hook into the app ───────────────────────── */
function install(Core){
  if(!Core||Core.__c3addon)return false;
  const orig=Core.fromName;
  Core.fromName=function(raw){
    const r=addonFromName(raw);
    if(r&&r.src)return r;
    if(r&&r.hard)return{err:r.err};            // name was read fine but is chemically impossible → say so, don't guess
    const o=orig.call(this,raw);
    return(o&&o.err&&r&&r.err)?{err:r.err}:o};  // if both fail, show the more specific message
  Core.fromNameOriginal=orig;Core.fromNameAddon=addonFromName;Core.__c3addon='1.0';return true}
if(typeof module!=='undefined')module.exports={install,fromName:addonFromName};
const hook=()=>root.Chem3DCore?install(root.Chem3DCore):false;
if(!hook()&&typeof document!=='undefined'){document.addEventListener('DOMContentLoaded',hook);root.addEventListener&&root.addEventListener('load',hook)}

/* ONE LINE to add to index.html, just before </body>  (same folder as index.html):
   <script src="c3-iupac-addon.js"></script>                                                                         */
})(typeof window!=='undefined'?window:globalThis);
