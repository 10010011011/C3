/* chemcore.js — everything here is DERIVED from a bond graph (no per-molecule tables).
   atoms: [{z:atomicNumber, q:formalCharge, x,y,z3}]   bonds: [[i,j,order]]            */
(function(g){
const SYM=' H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe'.split(' ');
const EN=[0,2.2,0,.98,1.57,2.04,2.55,3.04,3.44,3.98,0,.93,1.31,1.61,1.9,2.19,2.58,3.16,0,.82,1,1.36,1.54,1.63,1.66,1.55,1.83,1.88,1.91,1.9,1.65,1.81,2.01,2.18,2.55,2.96,3,.82,.95,1.22,1.33,1.6,2.16,1.9,2.2,2.28,2.2,1.93,1.69,1.78,1.96,2.05,2.1,2.66,2.6];
const GRP=[0,1,18,1,2,13,14,15,16,17,18,1,2,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
const HYB={2:'sp',3:'sp²',4:'sp³',5:'sp³d',6:'sp³d²'};
const GEO={2:['Linear'],3:['Trigonal planar','Bent'],4:['Tetrahedral','Trigonal pyramidal','Bent'],
  5:['Trigonal bipyramidal','Seesaw','T-shaped','Linear'],6:['Octahedral','Square pyramidal','Square planar','T-shaped','Linear']};
const fmtOS=v=>v>0?'+'+v:v<0?'−'+(-v):'0';

function analyze(atoms,bonds){
  const nb=atoms.map(()=>[]);
  bonds.forEach(([a,b,o])=>{nb[a].push([b,o]);nb[b].push([a,o]);});
  return atoms.map((at,i)=>{
    const Z=at.z,n=nb[i].length,sum=nb[i].reduce((s,x)=>s+x[1],0);
    // Oxidation state: formal charge + bond orders to MORE electronegative neighbours − orders to LESS electronegative ones
    let os=at.q;
    nb[i].forEach(([j,o])=>{const d=(EN[atoms[j].z]||2)-(EN[Z]||2);os+=d>.001?o:d<-.001?-o:0;});
    // VSEPR: only for s/p-block (d-block has no octet rule)
    const gp=GRP[Z],main=gp&&(gp<=2||gp>=13);
    let lp=null,sn=null,hyb='—',shape='—';
    if(main){
      const V=gp<=2?gp:gp-10;
      lp=Math.max(0,Math.floor((V-at.q-sum)/2));sn=n+lp;
      if(n<=1)shape='Terminal';
      else if(GEO[sn]&&GEO[sn][lp])shape=GEO[sn][lp];
      if(n>=1&&HYB[sn])hyb=HYB[sn];
    }
    return{i,sym:SYM[Z]||'?',z:Z,n,os,lp,sn,hyb,shape,sum};
  });
}
function formula(atoms){
  const c={};atoms.forEach(a=>{const s=SYM[a.z]||'?';c[s]=(c[s]||0)+1;});
  const hasC=!!c.C;
  const ks=Object.keys(c).sort((a,b)=>hasC&&(a==='C'||b==='C')?(a==='C'?-1:1):hasC&&(a==='H'||b==='H')?(a==='H'?-1:1):a<b?-1:1);
  const sub=n=>String(n).replace(/\d/g,d=>'₀₁₂₃₄₅₆₇₈₉'[d]);
  const q=atoms.reduce((s,a)=>s+a.q,0);
  const sup=q?(Math.abs(q)>1?Math.abs(q):'').toString().replace(/\d/g,d=>'⁰¹²³⁴⁵⁶⁷⁸⁹'[d])+(q>0?'⁺':'⁻'):'';
  return ks.map(k=>k+(c[k]>1?sub(c[k]):'')).join('')+sup;
}

/* ── SMILES parser: organic subset, aromatic rings (kekulized), brackets, hypervalent. Adds explicit H. ── */
function parseSmiles(s){
  const atoms=[],bonds=[],stack=[],ring={},V={B:[3],C:[4],N:[3,5],O:[2],F:[1],P:[3,5],S:[2,4,6],Cl:[1,3,5,7],Br:[1,3,5,7],I:[1,3,5,7]};
  let prev=-1,bo=1,bx=false;
  const link=(a,b,o,x)=>bonds.push([a,b,o,(!x&&atoms[a].ar&&atoms[b].ar)?1:0]);
  const put=(sym,q,h,org,ar)=>{const z=SYM.indexOf(sym);if(z<1)throw new Error('Unknown element '+sym);
    const k=atoms.length;atoms.push({z,q,org,sym,ar});if(prev>=0)link(prev,k,bo,bx);prev=k;bo=1;bx=false;
    for(let m=0;m<h;m++){atoms.push({z:1,q:0});bonds.push([k,atoms.length-1,1,0]);}};
  try{
    for(let i=0;i<s.length;i++){const c=s[i];
      if(c==='(')stack.push(prev);else if(c===')')prev=stack.pop();
      else if(c==='='){bo=2;bx=true}else if(c==='#'){bo=3;bx=true}else if(c==='-'){bx=true}else if('/\\:'.includes(c)){}
      else if(c==='.')prev=-1;
      else if(/\d/.test(c)){if(ring[c]){const[a,o,x]=ring[c];link(a,prev,Math.max(o,bo),x||bx);delete ring[c]}else ring[c]=[prev,bo,bx];bo=1;bx=false}
      else if(c==='['){const e=s.indexOf(']',i);if(e<0)throw new Error('Unclosed [');
        const m=s.slice(i+1,e).match(/^\d*([A-Z][a-z]?|[cnops])@*(?:H(\d*))?(?:(\++|-+)(\d*))?/);if(!m)throw new Error('Bad bracket atom');
        const sg=m[3]?(m[3][0]==='+'?1:-1):0,q=sg*(m[4]?+m[4]:(m[3]||'').length),ar=/^[a-z]$/.test(m[1]);
        put(ar?m[1].toUpperCase():m[1],q,m[2]===undefined?0:(m[2]===''?1:+m[2]),false,ar);i=e;}
      else if(/[A-Z]/.test(c)){const t=s.slice(i,i+2),sym=(t==='Cl'||t==='Br')?t:c;if(!V[sym])throw new Error('Unsupported atom '+sym);put(sym,0,0,true,false);i+=sym.length-1;}
      else if('cnosp'.includes(c))put(c.toUpperCase(),0,0,true,true);
      else if(c!==' ')throw new Error('Unexpected '+c);
    }
    // aromatic ring bonds → one valid Kekulé structure (backtracking matching)
    const inRing=b=>{const seen=new Set([b[0]]),q=[b[0]];while(q.length){const u=q.pop();for(const e of bonds){if(e===b)continue;
      const v=e[0]===u?e[1]:e[1]===u?e[0]:-1;if(v>=0&&!seen.has(v)){if(v===b[1])return true;seen.add(v);q.push(v)}}}return false};
    bonds.forEach(b=>{if(b[3]&&!inRing(b))b[3]=0});
    const ab=bonds.filter(b=>b[3]);
    if(ab.length){const deg=atoms.map(()=>0);bonds.forEach(([x,y])=>{deg[x]++;deg[y]++});
      const exo=atoms.map(()=>false);bonds.forEach(b=>{if(!b[3]&&b[2]>=2){exo[b[0]]=exo[b[1]]=true}});
      const need=atoms.map((a,k)=>!!a.ar&&!exo[k]&&(a.sym==='C'||(a.sym==='N'&&deg[k]<=2))),used=atoms.map(()=>false);
      const go=()=>{const u=need.findIndex((x,k)=>x&&!used[k]);if(u<0)return true;
        for(const b of ab){const v=b[0]===u?b[1]:b[1]===u?b[0]:-1;if(v<0||!need[v]||used[v])continue;
          used[u]=used[v]=true;b.dbl=1;if(go())return true;used[u]=used[v]=false;b.dbl=0}return false};
      if(!go())throw new Error('Cannot kekulize aromatic ring');
      ab.forEach(b=>{b[2]=b.dbl?2:1});}
    const sum=atoms.map(()=>0);bonds.forEach(([x,y,o])=>{sum[x]+=o;sum[y]+=o});const n0=atoms.length;
    for(let k=0;k<n0;k++){const a=atoms[k];if(!a.org)continue;const v=V[a.sym].find(v=>v>=sum[k]);const h=v===undefined?0:v-sum[k];
      for(let m=0;m<h;m++){atoms.push({z:1,q:0});bonds.push([k,atoms.length-1,1,0]);}}
    return{atoms:atoms.map(({z,q})=>({z,q,x:0,y:0,z3:0})),bonds:bonds.map(b=>[b[0],b[1],b[2],b[3]?1:0])};
  }catch(e){return{err:e.message}}
}

/* ── 3D: ideal VSEPR slots → tree placement → force-field relaxation (closes rings, compresses lone-pair angles) ── */
const RC={1:.31,5:.84,6:.76,7:.71,8:.66,9:.57,14:1.11,15:1.07,16:1.05,17:1.02,35:1.2,53:1.39,54:1.4};
const rc=z=>RC[z]||1.3;
const blen=(a,b,o,ar)=>(rc(a)+rc(b))*(ar?.915:([1,1,.88,.79][o]||.75));
const nz=v=>{const l=Math.hypot(v[0],v[1],v[2])||1;return[v[0]/l,v[1]/l,v[2]/l]};
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const mul=(a,k)=>[a[0]*k,a[1]*k,a[2]*k],add=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],dist=(a,b)=>Math.hypot(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
const E=a=>[Math.cos(a),Math.sin(a),0],TET=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(nz);
const SLOTS={1:[[0,0,1]],2:[[0,0,1],[0,0,-1]],3:[E(0),E(2.0944),E(4.1888)],4:TET,
  5:[[0,0,1],[0,0,-1],E(0),E(2.0944),E(4.1888)],6:[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]};
function rotAxis(k,t){const c=Math.cos(t),s=Math.sin(t);return v=>add(add(mul(v,c),mul(cross(k,v),s)),mul(k,dot(k,v)*(1-c)))}
function rotTo(u,v){const c=dot(u,v);if(c>.99999)return x=>x;
  if(c<-.99999){const p=nz(Math.abs(u[0])<.9?cross(u,[1,0,0]):cross(u,[0,1,0]));return rotAxis(p,Math.PI)}
  return rotAxis(nz(cross(u,v)),Math.acos(c))}
function relax(atoms,bonds,info,nb,asg,flatFirst){
  const n=atoms.length,P=atoms.map(a=>[a.x,a.y,a.z3]),key=(i,j)=>i<j?i*1e4+j:j*1e4+i,L=new Map(),cons=[],skip=new Set();
  bonds.forEach(([a,b,o,ar])=>{const l=blen(atoms[a].z,atoms[b].z,o,ar);L.set(key(a,b),l);cons.push([a,b,l,1]);skip.add(key(a,b))});
  nb.forEach((list,c)=>{const sn=Math.min(6,Math.max(list.length,info[c].sn||list.length)),S=SLOTS[sn],lp=info[c].lp||0;
    for(let x=0;x<list.length;x++)for(let y=x+1;y<list.length;y++){const a=list[x][0],b=list[y][0],sa=asg[c][a],sb=asg[c][b];
      if(sa===undefined||sb===undefined||!S[sa]||!S[sb])continue;
      let th=Math.acos(Math.max(-1,Math.min(1,dot(S[sa],S[sb]))));
      if(sn<=4&&lp&&th>1.8)th-=(sn===4?.0436:.035)*lp;            // −2.5° (sn4) / −2° (sn3) per lone pair
      const la=L.get(key(c,a)),lb=L.get(key(c,b));
      cons.push([a,b,Math.sqrt(Math.max(.01,la*la+lb*lb-2*la*lb*Math.cos(th))),.6]);skip.add(key(a,b));}});
  const rep=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(!skip.has(key(i,j)))rep.push([i,j,(atoms[i].z===1||atoms[j].z===1)?2.0:2.7]);
  const imp=[];nb.forEach((l,c)=>{if(info[c].sn===3&&l.length===3&&!info[c].lp)imp.push([c,l[0][0],l[1][0],l[2][0]])});
  const V=P.map(()=>[0,0,0]);
  const run=(iters,flat)=>{for(let it=0;it<iters;it++){
    const F=P.map(()=>[0,0,0]);
    const push=(i,j,d0,k,r)=>{const v=[P[i][0]-P[j][0],P[i][1]-P[j][1],P[i][2]-P[j][2]],d=Math.hypot(v[0],v[1],v[2])||1e-6;
      if(r&&d>=d0)return;const f=k*(d0-d)/d;for(let m=0;m<3;m++){F[i][m]+=f*v[m];F[j][m]-=f*v[m]}};
    cons.forEach(([i,j,d0,k])=>push(i,j,d0,k,false));rep.forEach(([i,j,d0])=>push(i,j,d0,.12,true));
    imp.forEach(([c,a,b,d])=>{const A=P[a],u=[P[b][0]-A[0],P[b][1]-A[1],P[b][2]-A[2]],v=[P[d][0]-A[0],P[d][1]-A[1],P[d][2]-A[2]],nn=nz(cross(u,v)),h=(P[c][0]-A[0])*nn[0]+(P[c][1]-A[1])*nn[1]+(P[c][2]-A[2])*nn[2];
      for(let m=0;m<3;m++){F[c][m]-=1.5*h*nn[m];[a,b,d].forEach(x=>{F[x][m]+=1.5*h*nn[m]/3})}});
    for(let i=0;i<n;i++)for(let m=0;m<3;m++){V[i][m]=.8*V[i][m]+.12*F[i][m];P[i][m]+=V[i][m]}
    if(flat)P.forEach(p=>{p[2]=0});
  }};
  if(flatFirst){P.forEach(p=>{p[2]=0});run(450,true);P.forEach((p,i)=>{p[2]=((i*7919)%13-6)/45});run(650,false)}else run(1500,false);
  atoms.forEach((a,i)=>{a.x=P[i][0];a.y=P[i][1];a.z3=P[i][2]});
}
function embed3D(atoms,bonds,info){
  const n=atoms.length,nb=atoms.map(()=>[]),asg=atoms.map(()=>({})),pos=Array(n).fill(null),dirs=Array(n),lpd=Array(n);
  bonds.forEach(([a,b,o])=>{nb[a].push([b,o]);nb[b].push([a,o]);});
  const full=i=>{const k=Math.max(1,nb[i].length),sn=Math.min(6,Math.max(k,info[i].sn||k));return SLOTS[sn].slice()};
  const frame=(i,back,t)=>{let S=full(i);if(back){S=S.map(rotTo(S[0],back));if(t)S=S.map(rotAxis(back,t));}
    const k=Math.max(1,nb[i].length);dirs[i]=S.slice(0,k);lpd[i]=S.slice(k);};
  let comp=0;
  for(;;){
    let r=-1;for(let i=0;i<n;i++)if(!pos[i]&&(r<0||nb[i].length>nb[r].length))r=i;
    if(r<0)break;
    pos[r]=[comp++*3.8,0,0];frame(r,null);
    const q=[[r,-1]];
    while(q.length){
      const[i,par]=q.shift();let slot=par<0?0:1;if(par>=0)asg[i][par]=0;
      nb[i].filter(([j])=>j!==par).forEach(([j,o])=>{
        const sl=slot++;asg[i][j]=sl;const d=dirs[i][sl];if(!d||pos[j])return;
        const ar=bonds.find(b=>(b[0]===i&&b[1]===j)||(b[0]===j&&b[1]===i))[3];
        pos[j]=add(pos[i],mul(d,blen(atoms[i].z,atoms[j].z,o,ar)));
        const back=mul(d,-1);let best=0,bs=-1;
        for(let t=0;t<12;t++){frame(j,back,t*Math.PI/6);let md=1e9;
          dirs[j].slice(1).forEach(v=>{const p=add(pos[j],mul(v,1.4));pos.forEach((a,k)=>{if(a&&k!==j)md=Math.min(md,dist(p,a))})});
          if(md>bs){bs=md;best=t}}
        frame(j,back,best*Math.PI/6);q.push([j,i]);
      });
    }
  }
  atoms.forEach((a,i)=>{const p=pos[i]||[0,0,0];a.x=p[0];a.y=p[1];a.z3=p[2];a.lps=lpd[i]||[];});
  if(n>2){const rm=new Set();let ch=true;
    while(ch){ch=false;for(let i=0;i<n;i++)if(!rm.has(i)&&nb[i].filter(([j])=>!rm.has(j)).length<=1){rm.add(i);ch=true}}
    const flat=bonds.some(b=>b[3])||atoms.some((_,i)=>!rm.has(i)&&info[i].sn===3&&nb[i].length===3);
    relax(atoms,bonds,info,nb,asg,false);}
}
function valid3D(atoms){let mn=1e9,mx=0;for(let i=0;i<atoms.length;i++)for(let j=i+1;j<atoms.length;j++){
  const d=Math.hypot(atoms[i].x-atoms[j].x,atoms[i].y-atoms[j].y,atoms[i].z3-atoms[j].z3);mn=Math.min(mn,d);mx=Math.max(mx,d);}
  return atoms.length<2?false:mn>.4&&mx>.8}
const api={analyze,formula,fmtOS,SYM,parseSmiles,embed3D,valid3D};
if(typeof module!=='undefined')module.exports=api;else g.Chem3DCore=api;
})(typeof window!=='undefined'?window:globalThis);
