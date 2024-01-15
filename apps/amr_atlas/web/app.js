// AMR Gene Atlas frontend — force-directed co-occurrence network
const CLASS_COLORS = {
  "beta-lactam": "var(--beta)", "carbapenem": "var(--carb)", "glycopeptide": "var(--glyco)",
  "tetracycline": "var(--tet)", "sulfonamide": "var(--sulf)", "aminoglycoside": "var(--amino)",
  "fluoroquinolone": "var(--fluoro)", "macrolide": "var(--macro)", "polymyxin": "var(--poly)",
};
const W = 760, H = 620;
let netData = null, nodes = [], edges = [], minLift = 1.5;

async function run() {
  const btn = document.getElementById("btnRun");
  btn.disabled = true;
  setStep("running pipeline…");
  const { job_id } = await fetch("/api/pipeline?n_genomes=140", {method:"POST"}).then(r=>r.json());
  await poll(job_id);
  setStep("");
  await load();
  btn.disabled = false;
}
function setStep(s){ document.getElementById("step").textContent = s; }
async function poll(jobId){
  for (let i=0;i<120;i++){
    const j = await fetch(`/api/jobs/${jobId}`).then(r=>r.json());
    if (j.progress && j.progress.stage) setStep(`${j.progress.stage}…`);
    if (j.status==="succeeded"||j.status==="failed") return j;
    await new Promise(r=>setTimeout(r,300));
  }
}

async function load() {
  const [net, summary] = await Promise.all([
    fetch(`/api/network?min_lift=${minLift}`).then(r=>r.json()),
    fetch(`/api/summary`).then(r=>r.json()),
  ]);
  netData = net;
  document.getElementById("rGenes").textContent = summary.genes || "—";
  document.getElementById("rGenomes").textContent = summary.genomes || "—";
  document.getElementById("rEdges").textContent = net.edges.length;
  buildGraph(net);
  drawEdgeList(net.edges);
  drawClassBars(summary.carriage_by_class);
  drawClassKey(net.nodes);
}

// ---- simple force simulation (Fruchterman-Reingold-ish) ----
function buildGraph(net) {
  const byId = {};
  nodes = net.nodes.map(n => {
    const o = {...n, x: W/2 + (Math.random()-0.5)*200, y: H/2 + (Math.random()-0.5)*200, vx:0, vy:0};
    byId[n.id] = o; return o;
  });
  edges = net.edges.map(e => ({...e, s: byId[e.source], t: byId[e.target]})).filter(e=>e.s&&e.t);

  // degree (for node sizing)
  nodes.forEach(n => n.deg = 0);
  edges.forEach(e => { e.s.deg++; e.t.deg++; });

  const iterations = 400;
  const area = W*H, kRep = Math.sqrt(area / Math.max(nodes.length,1)) * 1.25;
  for (let it=0; it<iterations; it++) {
    const t = 1 - it/iterations;
    // repulsion
    for (let i=0;i<nodes.length;i++) for (let j=i+1;j<nodes.length;j++){
      const a=nodes[i], b=nodes[j];
      let dx=a.x-b.x, dy=a.y-b.y, d=Math.hypot(dx,dy)||0.1;
      const f = (kRep*kRep)/d;
      const ux=dx/d, uy=dy/d;
      a.vx += ux*f*0.0006; a.vy += uy*f*0.0006;
      b.vx -= ux*f*0.0006; b.vy -= uy*f*0.0006;
    }
    // attraction along edges, stronger for higher lift
    edges.forEach(e => {
      let dx=e.s.x-e.t.x, dy=e.s.y-e.t.y, d=Math.hypot(dx,dy)||0.1;
      const f = (d*d)/kRep * (0.014 * Math.min(e.lift,4));
      const ux=dx/d, uy=dy/d;
      e.s.vx -= ux*f*0.001; e.s.vy -= uy*f*0.001;
      e.t.vx += ux*f*0.001; e.t.vy += uy*f*0.001;
    });
    // gravity to center + integrate
    nodes.forEach(n => {
      n.vx += (W/2 - n.x)*0.004; n.vy += (H/2 - n.y)*0.004;
      n.x += n.vx * (3*t+0.5); n.y += n.vy * (3*t+0.5);
      n.vx *= 0.82; n.vy *= 0.82;
      n.x = Math.max(45, Math.min(W-45, n.x));
      n.y = Math.max(45, Math.min(H-45, n.y));
    });
  }
  draw();
}

function draw() {
  const svg = document.getElementById("net");
  const ns = "http://www.w3.org/2000/svg";
  svg.innerHTML = "";
  const tip = document.getElementById("tip");
  const stage = document.querySelector(".net-stage");

  const maxLift = Math.max(...edges.map(e=>e.lift), 2);
  edges.forEach(e => {
    const line = document.createElementNS(ns, "line");
    line.setAttribute("x1", e.s.x); line.setAttribute("y1", e.s.y);
    line.setAttribute("x2", e.t.x); line.setAttribute("y2", e.t.y);
    line.setAttribute("class", "edge-line");
    line.setAttribute("stroke-width", 0.5 + (e.lift/maxLift)*3);
    line.setAttribute("stroke-opacity", 0.25 + (e.lift/maxLift)*0.5);
    svg.appendChild(line);
  });

  const maxDeg = Math.max(...nodes.map(n=>n.deg), 1);
  const maxPrev = Math.max(...nodes.map(n=>n.prevalence), 1);
  nodes.forEach(n => {
    const r = 5 + (n.prevalence/maxPrev)*11;
    const c = document.createElementNS(ns, "circle");
    c.setAttribute("cx", n.x); c.setAttribute("cy", n.y); c.setAttribute("r", r);
    c.setAttribute("fill", CLASS_COLORS[n.drug_class] || "var(--accent)");
    c.setAttribute("fill-opacity", n.deg>0 ? 0.92 : 0.35);
    c.setAttribute("stroke", "rgba(0,0,0,0.4)"); c.setAttribute("stroke-width", "1");
    c.setAttribute("class", "node-circle");
    c.addEventListener("mousemove", () => {
      const sr = stage.getBoundingClientRect(), svgr = svg.getBoundingClientRect();
      tip.style.left = (svgr.left - sr.left + (n.x/W)*svgr.width) + "px";
      tip.style.top  = (svgr.top - sr.top + (n.y/H)*svgr.height) + "px";
      tip.style.opacity = "1";
      tip.innerHTML = `<div class="tname">${n.id}</div><div class="tmeta">${n.drug_class} · in ${n.prevalence} genomes · ${n.deg} links</div>`;
    });
    c.addEventListener("mouseleave", () => tip.style.opacity = "0");
    svg.appendChild(c);

    const label = document.createElementNS(ns, "text");
    label.setAttribute("x", n.x); label.setAttribute("y", n.y - r - 3);
    label.setAttribute("text-anchor", "middle"); label.setAttribute("class", "node-label");
    label.textContent = n.id;
    svg.appendChild(label);
  });
}

function drawEdgeList(edgeArr) {
  const el = document.getElementById("edgeList");
  if (!edgeArr.length) { el.innerHTML = '<p class="empty">No co-occurring pairs above threshold.</p>'; return; }
  el.innerHTML = edgeArr.slice(0,8).map(e =>
    `<div class="edge"><div><div class="epair">${e.source} + ${e.target}</div>
     <div class="ecount">${e.count} shared genomes</div></div>
     <div class="elift">${e.lift}×</div></div>`
  ).join("");
}

function drawClassBars(byClass) {
  if (!byClass) return;
  const max = Math.max(...byClass.map(c=>c.carriage), 1);
  document.getElementById("classBars").innerHTML = byClass.map(c => `
    <div class="bar-row">
      <div class="bar-top"><span>${c.drug_class}</span><span>${c.carriage}</span></div>
      <div class="bar-track"><div class="bar-fill" style="width:${(c.carriage/max)*100}%;background:${CLASS_COLORS[c.drug_class]||'var(--accent)'}"></div></div>
    </div>`).join("");
}

function drawClassKey(nodeArr) {
  const classes = [...new Set(nodeArr.map(n=>n.drug_class))];
  document.getElementById("classKey").innerHTML = classes.map(c =>
    `<div class="ck"><span class="swatch" style="background:${CLASS_COLORS[c]||'var(--accent)'}"></span>${c}</div>`
  ).join("");
}

document.getElementById("liftSlider").addEventListener("input", (e) => {
  minLift = +e.target.value;
  document.getElementById("liftVal").textContent = minLift.toFixed(1);
  if (netData) load();
});
document.getElementById("btnRun").addEventListener("click", run);
document.getElementById("toPhylo").addEventListener("click", (ev)=>{ev.preventDefault(); alert("Run the Phylogeny Builder app on its own port (see README).");});

load();
