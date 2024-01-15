// Phylogeny Builder frontend
const CLADE_COLORS = {
  "Pre-1990": "var(--c0)", "Sydney-like": "var(--c1)", "Fujian-like": "var(--c2)",
  "Perth-like": "var(--c3)", "HongKong-like": "var(--c4)",
};
const CX = 380, CY = 320, R_MAX = 280;
let treeData = null, leaves = [], yearMin = 0, yearMax = 0;

// ---- radial layout: assign each leaf an angle, each node a radius by depth ----
function layout(root) {
  const allLeaves = [];
  (function collect(n){ if(!n.children) allLeaves.push(n); else n.children.forEach(collect); })(root);
  const N = allLeaves.length;
  allLeaves.forEach((lf, i) => { lf._angle = (i / N) * 2 * Math.PI; });

  let maxDepth = 0;
  (function depth(n, d){ n._depth = d; maxDepth = Math.max(maxDepth, d);
    if (n.children) n.children.forEach(c => depth(c, d + (c.branch_length || 0.01))); })(root, 0);

  (function place(n){
    if (!n.children) { n._r = R_MAX; return; }
    n.children.forEach(place);
    const angles = n.children.map(c => c._angle);
    n._angle = angles.reduce((a,b)=>a+b,0) / angles.length;
    n._r = (n._depth / maxDepth) * R_MAX;
  })(root);

  return allLeaves;
}

function polar(angle, r) { return [CX + r * Math.cos(angle - Math.PI/2), CY + r * Math.sin(angle - Math.PI/2)]; }

function render(root, cutoffYear) {
  const svg = document.getElementById("tree");
  const ns = "http://www.w3.org/2000/svg";
  svg.innerHTML = "";

  // faint concentric grid
  for (let g = 1; g <= 4; g++) {
    const c = document.createElementNS(ns, "circle");
    c.setAttribute("cx", CX); c.setAttribute("cy", CY); c.setAttribute("r", (g/4)*R_MAX);
    c.setAttribute("fill", "none"); c.setAttribute("stroke", "var(--line)"); c.setAttribute("stroke-width", "0.4");
    svg.appendChild(c);
  }

  const tip = document.getElementById("tip");
  const stage = document.querySelector(".tree-stage");

  (function draw(n){
    const [nx, ny] = polar(n._angle, n._r);
    if (n.children) {
      n.children.forEach(ch => {
        const visible = isVisible(ch, cutoffYear);
        const [cx, cy] = polar(ch._angle, ch._r);
        // radial elbow: arc along parent radius then straight out
        const [mx, my] = polar(ch._angle, n._r);
        const path = document.createElementNS(ns, "path");
        const sweep = ch._angle > n._angle ? 1 : 0;
        path.setAttribute("d", `M${nx},${ny} A${n._r},${n._r} 0 0 ${sweep} ${mx},${my} L${cx},${cy}`);
        path.setAttribute("class", "branch-line");
        path.setAttribute("opacity", visible ? "1" : "0.07");
        svg.appendChild(path);
        draw(ch);
      });
    } else {
      const visible = isVisible(n, cutoffYear);
      const color = CLADE_COLORS[n.clade] || "var(--leaf)";
      const dot = document.createElementNS(ns, "circle");
      dot.setAttribute("cx", nx); dot.setAttribute("cy", ny);
      dot.setAttribute("r", visible ? 3 : 1.2);
      dot.setAttribute("fill", visible ? color : "var(--faint)");
      dot.setAttribute("opacity", visible ? "1" : "0.12");
      dot.setAttribute("class", "leaf-dot");
      dot.addEventListener("mousemove", () => {
        const sr = stage.getBoundingClientRect();
        const svgr = svg.getBoundingClientRect();
        tip.style.left = (svgr.left - sr.left + (nx/760)*svgr.width) + "px";
        tip.style.top  = (svgr.top - sr.top + (ny/640)*svgr.height) + "px";
        tip.style.opacity = "1";
        tip.innerHTML = `<div class="tname">${n.name}</div><div class="tmeta">${n.clade||"—"} · ${n.year||"?"} · ${n.country||""}</div>`;
      });
      dot.addEventListener("mouseleave", () => tip.style.opacity = "0");
      svg.appendChild(dot);
    }
  })(root);
}

function isVisible(node, cutoffYear) {
  if (cutoffYear == null) return true;
  if (!node.children) return (node.year || 0) <= cutoffYear;
  return node.children.some(c => isVisible(c, cutoffYear));
}

// ---- data flow ----
async function poll(jobId) {
  for (let i=0;i<120;i++){
    const j = await fetch(`/api/jobs/${jobId}`).then(r=>r.json());
    if (j.progress && j.progress.stage) setStep(`${j.progress.stage}…`);
    if (j.status==="succeeded"||j.status==="failed") return j;
    await new Promise(r=>setTimeout(r,300));
  }
}
function setStep(s){ document.getElementById("step").textContent = s; }

async function run() {
  const btn = document.getElementById("btnRun");
  btn.disabled = true;
  setStep("ingesting sequences…");
  await fetch("/api/ingest?term=Influenza%20A%20virus%5BOrganism%5D&retmax=200&query_tag=default", {method:"POST"});
  setStep("building tree…");
  const { job_id } = await fetch("/api/build?query_tag=default", {method:"POST"}).then(r=>r.json());
  await poll(job_id);
  setStep("");
  await loadTree();
  btn.disabled = false;
}

async function loadTree() {
  const t = await fetch("/api/tree?query_tag=default").then(r=>r.ok?r.json():null);
  if (!t) return;
  treeData = t.tree;
  document.getElementById("rTaxa").textContent = t.n_taxa;
  document.getElementById("newick").textContent = t.newick;
  leaves = layout(treeData);
  const yrs = leaves.map(l=>l.year).filter(Boolean);
  yearMin = Math.min(...yrs); yearMax = Math.max(...yrs);
  const slider = document.getElementById("yearSlider");
  slider.min = yearMin; slider.max = yearMax; slider.value = yearMax;
  document.getElementById("yearLabel").textContent = yearMax;
  document.getElementById("timeControl").hidden = false;
  render(treeData, yearMax);
  drawCladeKey();
  await loadTimeline();
}

function drawCladeKey() {
  const present = [...new Set(leaves.map(l=>l.clade).filter(Boolean))];
  document.getElementById("cladeKey").innerHTML = present.map(c =>
    `<div class="ck"><span class="swatch" style="background:${CLADE_COLORS[c]||'var(--leaf)'}"></span>${c}</div>`
  ).join("");
}

async function loadTimeline() {
  const tl = await fetch("/api/timeline?query_tag=default").then(r=>r.json());
  const em = document.getElementById("emergence");
  const entries = Object.entries(tl.clade_emergence);
  em.innerHTML = entries.map(([clade, year]) =>
    `<div class="erow"><span class="eswatch" style="background:${CLADE_COLORS[clade]||'var(--leaf)'}"></span>
     <span class="ename">${clade}</span><span class="eyear">${year}</span></div>`
  ).join("");
}

// ---- timeline animation ----
let playing = false;
document.getElementById("yearSlider").addEventListener("input", (e) => {
  const y = +e.target.value;
  document.getElementById("yearLabel").textContent = y;
  if (treeData) render(treeData, y);
});
document.getElementById("btnPlay").addEventListener("click", () => {
  if (playing) return;
  playing = true;
  const slider = document.getElementById("yearSlider");
  let y = yearMin;
  const tick = () => {
    slider.value = y;
    document.getElementById("yearLabel").textContent = y;
    render(treeData, y);
    y++;
    if (y <= yearMax) setTimeout(tick, 220);
    else playing = false;
  };
  tick();
});
document.getElementById("btnRun").addEventListener("click", run);
document.getElementById("toAmr").addEventListener("click", (e) => {
  e.preventDefault();
  alert("Run the AMR Atlas app separately on its own port (see README).");
});

// try to load an existing tree on first paint
loadTree();
