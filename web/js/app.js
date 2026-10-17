// LocalSH Weave Lab — app state, UI bindings and the frame loop.
import { Carrier, parseSTL } from './carrier.js';
import { CpuSolver } from './solver_cpu.js';
import { GpuSolver, getDevice } from './solver_gpu.js';
import { extractLattice, chargeSummary, traceStrands, edgeStrain, orientationField } from './lattice.js';
import { View } from './view.js';
import { valenceHex, FAMILY_HEX, lutColor, misorRGB } from './colormap.js';
import { MinHeap } from './sparse.js';
import { frontState, predictGrowth, mapValue } from './frontspeed.js';
import { FRONT } from './frontspeed_table.js';
import { directionsFromNb0 } from './rosy.js';

const SQ3 = Math.sqrt(3);
const $ = (id) => document.getElementById(id);
const U_RANGE = [-0.45, 1.25];

const S = {
  manifest: [], carrier: null, solver: null, device: null,
  a: 2, valence: 6, centre: -1, centreMode: 'detected', detected: null,
  target: 3000, maxSteps: 20000, spf: 25, cgIters: 12, playing: false, busy: false,
  r: -0.08, b: 1.2, dt: 0.1, seedKind: 'wave', seedR: 2.0, seedPhi: 0, disp: true, seedMode: 'local', rosy: null, stripesRes: null, globalReq: 0,
  lat: null, summary: null, strands: null, history: [], lastExtract: -1, extractEvery: 100, shownStep: -1, scrubbed: false,
  carrierView: 'u', carrierAlpha: 1, stripeAlpha: [1, 1, 1], stripesShown: null, lathW: 0.18, lathCol: 'family', overUnder: true, msg: '',
  theme: 'auto', showValence: {}, selLath: -1, selLathKey: null, selVertex: -1, selVertexKey: -1, armPick: false,
  turntable: false, turnSpeed: 12, projection: 'persp',
  refine: 0, baseCarrier: null, side: 1, cgManual: false,
  anneal: { on: false, r2: -0.027, b2: 0.69, when: 'auto', step: 5000, active: false, startedAt: -1 },
  strainFam: -1, strainRange: 0.15, es: null, orient: null, orientMode: 'mis',
  geo: null, pred: null,                                   // geodesic distances from the centre / to the rim, and the growth prediction
};
const view = new View($('gl'));
view.attachNav($('navcube'));

// ------------------------------------------------------------------ helpers
const fmt = (x, d = 1) => (Number.isFinite(x) ? x.toFixed(d) : '—');
function setStatus(msg) { S.msg = msg; $('st-msg').textContent = msg; }
/** one-line summary shown in a section header (visible when the section is folded) */
function setSum(id, text) { const el = $('sum-' + id); if (el && el.textContent !== text) el.textContent = text; }
// the lattice edge must be resolved by the carrier (≥ 6 edges per a) and the time-step operator is conditioned by
// (a/h)² (κ ≈ 13 at 8 edges per a, ≈ 46 at 14), so the range follows the carrier resolution; refine to go denser
function aRange() { const h = S.carrier ? S.carrier.h : 0.25; return [6 * h, 14 * h]; }
function autoCgIters() { const c = S.carrier; if (!c) return 12; const epa = S.a / c.h; return Math.max(8, Math.min(40, Math.round(12 * (epa / 9) ** 2))); }
function applyCgIters() { if (S.cgManual) return; S.cgIters = autoCgIters(); $('nb-cg').value = S.cgIters; if (S.solver) S.solver.cgIters = S.cgIters; }
function aFromSlider(v) { const [lo, hi] = aRange(); return lo * Math.pow(hi / lo, v / 1000); }
function sliderFromA(a) { const [lo, hi] = aRange(); return Math.round(1000 * Math.log(a / lo) / Math.log(hi / lo)); }
function sites() { return S.carrier ? S.carrier.area / (SQ3 / 2 * S.a * S.a) : 0; }

// ------------------------------------------------------------------ carriers
async function loadManifest() {
  try { S.manifest = await (await fetch('manifest.json')).json(); } catch (e) { S.manifest = []; }
  const sel = $('sel-zoo'); sel.innerHTML = '';
  for (const m of S.manifest) { const o = document.createElement('option'); o.value = m.name; o.textContent = `${m.name}  (${m.n} v, ∫K ${m.omega_deg.toFixed(0)}°)`; sel.appendChild(o); }
  const o = document.createElement('option'); o.value = '__file__'; o.textContent = 'dropped / opened file'; o.hidden = true; sel.appendChild(o);
}

async function loadCarrierByName(name) {
  const m = S.manifest.find(x => x.name === name); if (!m) return;
  setStatus(`loading ${name}…`);
  try { const c = await Carrier.fromURL(m.file); setCarrier(c, m.a_default); $('sel-zoo').value = name; }
  catch (e) { setStatus(`load failed: ${e.message}`); console.error(e); }
}

/** opts.refine: the same surface at another resolution — keep the camera, the centre and the selection */
function setCarrier(c, aDefault, opts = {}) {
  S.carrier = c; S.baseCarrier = c.base || c; S.refine = c.level || 0; S.rosy = null; S.stripesRes = null; S.stripesShown = null;
  S.playing = false; S.lat = null; S.summary = null; S.strands = null; S.history = []; S.lastExtract = -1;
  S.extractEvery = 100 * Math.max(1, Math.round(c.n / 12000));              // the read-out is O(n log n): less often on fine carriers
  const [lo, hi] = aRange(); S.a = Math.min(hi, Math.max(lo, aDefault || c.bbox / 13));
  $('sl-a').value = sliderFromA(S.a); $('nb-a').value = S.a.toFixed(3);
  if (!opts.refine) { S.centreMode = 'detected'; S.centrePos = null; S.side = 1; view.setSide(1); renderSideChips(); }
  S.geo = null;
  view.setCarrier(c, { keepCamera: !!opts.refine });
  detect();
  applyCgIters();
  updateMeshLayer();
  newSolver(true);
  renderRefineChips();
  const q = c.meta;
  $('ro-carrier').textContent = `${c.meta.input_vertices ? `input ${c.meta.input_vertices} v → ` : ''}carrier ${c.n} v · ${c.F.length / 3} f · χ = ${c.chi} · h = ${c.h.toFixed(3)}` +
    (q.min_angle != null ? ` · min angle ${q.min_angle.toFixed(0)}° · ${q.skinny_pct.toFixed(1)}% skinny` : '') + (q.flips != null ? ` · iDT flips ${q.flips}` : '');
  const notes = (q.notes || []); $('ro-notes').hidden = !notes.length; $('ro-notes').textContent = notes.join('; ');
  $('st-mesh').textContent = `${c.name} · ${c.n} v · ${c.F.length / 3} f · χ = ${c.chi} · ∫K = ${(c.omega * 180 / Math.PI).toFixed(0)}°`;
  setSum('surface', `${c.name} · ${c.n} v`);
  updateDensityReadout(); updateCarrierColors();
  setStatus(`carrier ready`);
}

function detect() {
  const c = S.carrier; if (!c) return;
  const d = c.detectSingularity(S.a); S.detected = d;
  if (S.centreMode === 'detected' || S.centre < 0 || S.centre >= c.n) { S.centre = d.centre; S.valence = Math.max(3, Math.min(9, d.valence)); S.centreMode = 'detected'; }
  $('sl-valence').value = S.valence; $('out-valence').value = S.valence;
  const i = d.info;
  $('ro-omega').innerHTML = `${i.mode} · ${fmt(i.withinA, 0)} % of Ω within a of the centre`;
  updateCentreReadout(); drawChargeDiagram();
}
function updateCentreReadout() {
  const c = S.carrier; if (!c) return;
  $('ro-centre').textContent = `centre: vertex ${S.centre} (${S.centreMode}) · defect ${fmt(c.K[S.centre] * 180 / Math.PI)}°`;
  setSum('sing', `v${S.valence} at ${S.centre} (${S.centreMode})`);
  view.setCentre(S.centre); updateProbe();
}

// ------------------------------------------------------------------ solver
function newSolver(reseed) {
  const c = S.carrier; if (!c) return;
  if (S.solver) S.solver.dispose();
  // snapshots (float32 per vertex) for scrubbing: ~400 of them, but never more than ~300 MB on fine carriers
  const snapEvery = Math.max(10, Math.round(S.maxSteps / 400 / 10) * 10, Math.ceil(S.maxSteps * c.n * 4 / 300e6 / 10) * 10);
  const opts = { a: S.a, r: S.r, b: S.b, dt: S.dt, cgIters: S.cgIters, snapEvery, dispersion: S.disp };
  const ops = c;
  try { S.solver = S.device ? new GpuSolver(S.device, ops, opts) : new CpuSolver(ops, opts); }
  catch (e) { console.warn('GPU solver failed, using CPU', e); S.device = null; S.solver = new CpuSolver(ops, opts); }
  stopAnneal(false);
  $('st-backend').textContent = S.solver.backend === 'webgpu' ? `WebGPU${S.device.info?.architecture ? ' · ' + S.device.info.architecture : ''}` : 'CPU (JS)';
  $('st-backend').classList.toggle('gpu', S.solver.backend === 'webgpu');
  if (reseed) reseedNow();
}

function reseedNow() {
  const c = S.carrier, s = S.solver; if (!c || !s) return;
  if (S.seedMode === 'global') { globalSeed(); return; }
  const u0 = c.seedField(S.centre, S.valence, S.a, { radius: S.seedR, kind: S.seedKind, phaseDeg: S.seedPhi });
  s.setInitial(u0);
  afterReseed(`v${S.valence} at ${S.centre} (${S.centreMode})`);
}

/** no seed: smoothest 6-RoSy on the carrier → stripe phases of its three families → hexagonal u0 everywhere.
 *  Both solves run in a module Worker on the base carrier and are lifted to the refined one (global_worker.js). */
let gw = null, gwId = 0; const gwPending = new Map();
function globalWorker() {
  if (gw) return gw;
  gw = new Worker(new URL('./global_worker.js', import.meta.url), { type: 'module' });
  gw.onmessage = (e) => { const d = e.data, p = gwPending.get(d.id); if (!p) return; if (d.stage) { p.onStage(d.stage); return; } gwPending.delete(d.id); d.error ? p.reject(new Error(d.error)) : p.resolve(d); };
  gw.onerror = (e) => { for (const p of gwPending.values()) p.reject(new Error(e.message || 'worker error')); gwPending.clear(); gw = null; };
  return gw;
}
function runGlobalWorker(c, q0, onStage) {
  return new Promise((resolve, reject) => {
    const id = ++gwId; gwPending.set(id, { resolve, reject, onStage });
    const pack = (cc) => ({ n: cc.n, Fi: cc.Fi, Li: cc.Li, B: cc.boundary, V: cc.V, N: cc.N, F: cc.F });
    const base = c.base || c, key = `${c.name}|${base.n}|${c.level || 0}|${c.n}`;
    globalWorker().postMessage({ id, key, base: pack(base), fine: c.base ? pack(c) : null, q0 });
  });
}
async function globalSeed(force = false) {
  const c = S.carrier, s = S.solver; if (!c || !s) return;
  const a = S.a, req = ++S.globalReq;
  if (!force && S.stripesRes && S.stripesRes.carrier === c && S.stripesRes.a === a) { s.setInitial(S.stripesRes.u0); afterReseed(S.stripesRes.sum); return; }
  S.seedPending = true; $('btn-rosy').disabled = true;
  const base = c.base || c, progress = (stage) => { setStatus(stage); $('ro-rosy').textContent = `${stage} (base carrier ${base.n} v${c.base ? `, refined ${c.n} v` : ''})`; $('ro-stripes').textContent = 'the solver waits for the seed'; try { if (S.solver) showStep(S.solver.steps); } catch (e) { /* ignore */ } };
  progress('starting the rotation-field worker…');
  try {
    const t0 = performance.now();
    const d = await runGlobalWorker(c, 4 * Math.PI / (SQ3 * a), (stage) => { if (req === S.globalReq) progress(stage); });
    if (req !== S.globalReq || S.carrier !== c) return;                      // superseded by a newer request
    const ms = performance.now() - t0;
    const hist = {}; for (const sg of d.singular) hist[sg.valence] = (hist[sg.valence] || 0) + 1;
    const list = d.singular.map(sg => { const f = sg.face, ids = [base.Fi[3 * f], base.Fi[3 * f + 1], base.Fi[3 * f + 2]]; return { valence: sg.valence, index: sg.index, face: f, pos: cen(base.V, ids), nrm: cen(base.N, ids) }; });
    const bc = d.branchOnBase ? base : c;
    const branch = d.branch.map(b => { const f = b.face, ids = [bc.Fi[3 * f], bc.Fi[3 * f + 1], bc.Fi[3 * f + 2]]; return { ...b, pos: cen(bc.V, ids), nrm: cen(bc.N, ids) }; });
    S.rosy = { carrier: c, list, hist, theta: d.theta, nb0: d.nb0, dirs: [0, 1, 2].map(k => directionsFromNb0(c, d.theta, d.nb0, 6, k)), density: null };
    const np = branch.filter(b => b.winding > 0).length, nm = branch.length - np;
    const sum = `global field · ${d.singular.length} field singularities · ${branch.length} branch points`;
    S.stripesRes = { carrier: c, a, u0: d.u0, branch, sum, faces: d.stripeFaces || null };
    const omega6 = c.omega * 180 / Math.PI / 60, parts = Object.entries(hist).map(([v, k]) => `${k} × v${v}`).join(', ') || 'none';
    const tt = (k) => d.t[k] !== undefined ? `${(d.t[k] / 1000).toFixed(1)} s` : '—';
    $('ro-rosy').innerHTML = `rotation field (on the base carrier, ${d.nBase} v${c.base ? `, lifted to ${c.n} v` : ''}): ${d.singular.length} singular faces (${parts}) · Σ index = ${d.sumIndex}/6 vs Ω/60° = ${omega6.toFixed(1)} ${Math.abs(d.sumIndex - omega6) < 0.6 ? '✓' : '(rim turning not counted)'} · ${tt('rosy')}${c.base ? ` + lift ${tt('lift')}` : ''}`;
    $('ro-stripes').innerHTML = `stripe phases at a = ${a.toFixed(2)}: ${branch.length} branch points (+${np} / −${nm}) ≈ ${Math.round(branch.length / 2)} dislocations forced by the geometry${c.base ? ` (${d.branchBase} on the base)` : ''} · worst family match ${(d.worstMatch * 180 / Math.PI).toFixed(0)}° · ${tt('stripes')}${c.base ? ` + lift ${tt('liftStripes')}` : ''} · total ${(ms / 1000).toFixed(1)} s`;
    view.setRosy(S.rosy.dirs, 0.4 * a);
    if (S.carrierView === 'rosysing' || S.carrierView === 'stripes') updateCarrierColors();
    s.setInitial(d.u0);
    afterReseed(sum);
    setStatus(`seeded from the rotation field (${d.singular.length} singularities, ${branch.length} branch points, ${(ms / 1000).toFixed(1)} s)`);
  } catch (e) { console.error(e); setStatus('rotation field failed: ' + e.message); $('ro-rosy').textContent = 'rotation field failed: ' + e.message; }
  finally { if (req === S.globalReq) { S.seedPending = false; $('btn-rosy').disabled = false; if (S.playing) pump(); } }
}
const cen = (A, ids) => [0, 1, 2].map(k => (A[3 * ids[0] + k] + A[3 * ids[1] + k] + A[3 * ids[2] + k]) / 3);

function afterReseed(sum) {
  stopAnneal(true);
  S.history = []; S.lastExtract = -1; S.lat = null; S.summary = null; S.strands = null; S.es = null; S.orient = null; S.shownStep = -1; S.scrubbed = false;
  view.setLattice(null, null); view.setLaths(null, null, {}); view.setStrain(null); view.setOrientation(null); clearSelection(false); updateDirectionReadout();
  setSum('sing', sum);
  updateProbe(); updatePrediction();
  showStep(0); updateLatticeReadout(); drawChart();
}

const yieldTask = () => new Promise(r => setTimeout(r, 0));

async function pump() {
  // advance the solver toward the target in ~120 ms batches, yielding to the event loop (rendering) between batches
  if (S.pumping) return; S.pumping = true;
  try {
    while (S.solver && S.solver.steps < S.target) {
      if (S.seedPending) break;                                               // the global seed is still being solved: nothing to step yet
      const s = S.solver, t0 = performance.now(), n = Math.min(S.spfAuto || S.spf, S.target - s.steps);
      try { await s.step(n); } catch (e) { console.error(e); setStatus('solver error: ' + e.message); S.playing = false; S.target = s.steps; break; }
      const ms = performance.now() - t0; S.rate = n / Math.max(ms, 1e-3) * 1000;
      S.spfAuto = Math.max(4, Math.min(S.spf, Math.round(S.rate * 0.12)));
      try {
        if (s.steps - S.lastExtract >= S.extractEvery || s.steps >= S.target) { extract(s.u, true); checkAnneal(); }
        showStep(s.steps);
      } catch (e) { console.error(e); setStatus('read-out error: ' + e.message); }
      if (S.playing && s.steps >= S.target) { if (S.target >= S.maxSteps) { S.playing = false; $('btn-run').textContent = 'Run'; } else { S.target = Math.min(S.maxSteps, S.target + Math.max(S.spfAuto || S.spf, 50)); syncStepsUI(); } }
      await yieldTask();
    }
  } catch (e) { console.error(e); setStatus('solver loop error: ' + e.message); }
  finally { S.pumping = false; }
}

function tick() {
  const s = S.solver; if (!s) return;
  if (S.playing && s.steps >= S.target && S.target < S.maxSteps) { S.target = Math.min(S.maxSteps, S.target + 50); syncStepsUI(); }
  if (s.steps < S.target) pump();
  else if (S.shownStep !== S.target && !S.pumping) showStep(S.target);
  $('ro-steps').textContent = `${s.steps} computed · t = ${s.t.toFixed(0)}` + (S.rate ? ` · ${S.rate.toFixed(0)} steps/s` : '') + (Number.isFinite(s.lastRelRes) ? ` · PCG relres ${s.lastRelRes.toExponential(1)}` : '');
  $('ro-cg').textContent = s.backend === 'webgpu' ? 'fixed iterations on the GPU' : `relres ${Number.isFinite(s.lastRelRes) ? s.lastRelRes.toExponential(1) : '—'}`;
  setSum('solve', `${s.steps} / ${S.maxSteps}${S.playing || s.steps < S.target ? ' · running' : ''}`);
}

function showStep(step) {
  const s = S.solver; if (!s) return;
  const u = s.stateAt(step); if (!u) return;
  S.shownStep = step;
  if (S.carrierView === 'u') view.setFieldColors(u, 'viridis', U_RANGE[0], U_RANGE[1]);
  if (step < s.steps) { extract(u, false); S.scrubbed = true; }               // scrubbing back: lattice of that snapshot
  else if (S.scrubbed) { extract(s.u, false); S.scrubbed = false; }           // back at the live state (periodic extraction otherwise)
  $('st-run').textContent = (S.seedPending ? '· waiting for the global seed ' : '') + `· step ${step} / ${s.steps} (t = ${s.t.toFixed(0)}) · ${S.seedMode === 'global' ? 'global seed' : `v${S.valence} seed`} · a = ${S.a.toFixed(3)}` + (S.anneal.active ? ` · anneal from ${S.anneal.startedAt}` : '') + (S.playing || s.steps < S.target ? ' · running' : '');
  updateGrowthLegend(); drawChart();
}

function extract(u, record) {
  const c = S.carrier, s = S.solver; if (!c) return;
  try { S.lat = extractLattice(c, u, S.a); } catch (e) { console.error(e); S.lat = null; }
  S.summary = chargeSummary(S.lat);
  S.strands = S.lat && S.lat.tris.length ? traceStrands(S.lat.peaks.length, S.lat.tris) : null;
  S.es = S.lat && S.strands ? edgeStrain(S.lat, S.strands, S.a, null) : null;
  S.orient = S.lat ? orientationField(S.lat, c, orientStart()) : null;
  if (record) { S.lastExtract = s.steps; S.history.push([s.steps, S.lat ? S.lat.peaks.length : 0, S.summary.nInterior, S.summary.singular.length]); drawChart(); }
  resolveSelection(); view.a = S.a;
  view.setLattice(S.lat, S.summary);
  view.setSelectedVertex(S.selVertex);
  updateLaths(); updateDirectionLayers();
  updateLatticeReadout();
}

/** the regular site nearest the seed centre: the orientation arrows are made consistent outward from it */
function orientStart() {
  const lat = S.lat, g = geoFor(); if (!lat || !g) return -1;
  let best = -1, bd = Infinity;
  for (let i = 0; i < lat.peaks.length; i++) { const d = g.dist[lat.peaks[i]]; if (d < bd && lat.interior[i] && lat.valence[i] === 6) { bd = d; best = i; } }
  return best;
}
function updateDirectionLayers() {
  view.setStrain(view.visible.strain ? S.lat : null, S.strands, S.es, { family: S.strainFam, range: S.strainRange });
  view.setOrientation(view.visible.orient ? S.lat : null, S.orient, S.orientMode);
  updateDirectionReadout();
}
function updateDirectionReadout() {
  const el = $('ro-direction'), of = S.orient, es = S.es;
  const show = view.visible.orient || view.visible.strain;
  el.hidden = !show; $('row-orient').hidden = !view.visible.orient; $('row-strain').hidden = !view.visible.strain; drawOrientLegend();
  if (!show) return;
  if (!of || !es) { el.textContent = 'no lattice yet'; return; }
  const parts = [];
  if (view.visible.strain) {
    const fam = [0, 1, 2].map(f => { let s1 = 0, s2 = 0, k = 0; for (let e = 0; e < es.fam.length; e++) if (es.fam[e] === f) { s1 += es.strain[e]; s2 += es.strain[e] ** 2; k++; } return k ? `<span class="sw" style="background:${FAMILY_HEX[f]}"></span>${'ABC'[f]} ${(100 * s1 / k).toFixed(1)} ± ${(100 * Math.sqrt(Math.max(s2 / k - (s1 / k) ** 2, 0))).toFixed(1)} % (${k})` : `${'ABC'[f]} —`; });
    let unk = 0; for (const f of es.fam) if (f < 0) unk++;
    parts.push(`strain ℓ / a − 1 per lath family (mean ± sd, edges): ${fam.join(' · ')}${unk ? ` · ${unk} edges unassigned` : ''}`);
  }
  if (view.visible.orient) {
    let k = 0, s = 0, hi = 0; for (const m of of.mis) if (Number.isFinite(m)) { s += m; k++; if (m > 7.5) hi++; }
    parts.push(`orientation on ${of.count} regular sites · misorientation to the neighbours: mean ${k ? (s / k).toFixed(1) : '—'}°, ${hi} sites above 7.5° (seams) · coherence |⟨e^{6iθ}⟩| = ${of.coherence.toFixed(2)}`);
  }
  el.innerHTML = parts.join('<br>');
}

// ------------------------------------------------------------------ two-stage schedule (grow → anneal)
function checkAnneal() {
  const A = S.anneal, s = S.solver; if (!A.on || A.active || !s) return;
  let go = false;
  if (A.when === 'step') go = s.steps >= A.step;
  else { const H = S.history; if (H.length >= 4 && s.steps >= 500) { const p0 = H[H.length - 4][1], p1 = H[H.length - 1][1]; go = p1 > 10 && p1 <= p0 * 1.01; } }
  if (go) startAnneal();
}
function startAnneal() {
  const A = S.anneal, s = S.solver; if (!s || A.active) return;
  A.active = true; A.startedAt = s.steps; s.setParams({ r: A.r2, b: A.b2 });
  $('out-r').value = A.r2.toFixed(3); $('sl-r').value = A.r2; $('out-b').value = A.b2.toFixed(2); $('sl-b').value = A.b2;
  $('st-backend').classList.add('anneal');
  renderAnnealReadout(); setStatus(`annealing from step ${A.startedAt}: r = ${A.r2}, b = ${A.b2} (field kept)`); updatePrediction(); drawChart();
}
function stopAnneal(restore) {
  const A = S.anneal; if (!A.active) { renderAnnealReadout(); return; }
  A.active = false; A.startedAt = -1;
  if (restore && S.solver) S.solver.setParams({ r: S.r, b: S.b });
  $('out-r').value = S.r.toFixed(3); $('sl-r').value = S.r; $('out-b').value = S.b.toFixed(2); $('sl-b').value = S.b;
  $('st-backend').classList.remove('anneal');
  renderAnnealReadout();
}
function renderAnnealReadout() {
  const A = S.anneal, el = $('ro-anneal');
  el.innerHTML = !A.on ? 'off' : A.active ? `<span class="warn">annealing since step ${A.startedAt} at r = ${A.r2}, b = ${A.b2}</span>` : A.when === 'step' ? `armed: at step ${A.step}` : 'armed: when the peak count stops growing';
  $('nb-anneal-step').hidden = A.when !== 'step';
}

function updateLaths() {
  if (!view.visible.laths) { view.setLaths(null, null, {}); return; }
  view.setLaths(S.lat, S.strands, { width: S.lathW * S.a, overUnder: S.overUnder, colouring: S.lathCol, selected: S.selLath });
}

// ------------------------------------------------------------------ selection (laths, vertices)
// selections are remembered by carrier vertex ids so they survive re-extraction while the lattice grows
function lathKey(si) {
  const st = S.strands, lat = S.lat; if (!st || !lat || si < 0 || si >= st.strands.length) return null;
  const s = st.strands[si], e = st.edges[s.knots[Math.floor(s.knots.length / 2)]];
  return [lat.peaks[e[0]], lat.peaks[e[1]]].sort((x, y) => x - y);
}
function findLathByKey(key) {
  const st = S.strands, lat = S.lat; if (!st || !lat || !key) return -1;
  for (let si = 0; si < st.strands.length; si++) for (const e of st.strands[si].knots) {
    const ed = st.edges[e], a = lat.peaks[ed[0]], b = lat.peaks[ed[1]];
    if ((a === key[0] && b === key[1]) || (a === key[1] && b === key[0])) return si;
  }
  return -1;
}
function resolveSelection() {
  const lat = S.lat;
  S.selVertex = -1; S.selLath = -1;
  if (!lat) return;
  if (S.selVertexKey >= 0) { for (let i = 0; i < lat.peaks.length; i++) if (lat.peaks[i] === S.selVertexKey) { S.selVertex = i; break; } }
  if (S.selLathKey) S.selLath = findLathByKey(S.selLathKey);
}
function selectLath(si) {
  S.selLath = si; S.selLathKey = lathKey(si); S.selVertex = -1; S.selVertexKey = -1;
  view.setSelectedVertex(-1); updateLaths(); updateSelectionReadout();
}
function selectVertex(peak) {
  S.selVertex = peak; S.selVertexKey = S.lat ? S.lat.peaks[peak] : -1; S.selLath = -1; S.selLathKey = null;
  view.setSelectedVertex(peak); updateLaths(); updateSelectionReadout();
}
function clearSelection(redraw = true) {
  S.selLath = -1; S.selLathKey = null; S.selVertex = -1; S.selVertexKey = -1;
  if (redraw) { view.setSelectedVertex(-1); updateLaths(); }
  updateSelectionReadout();
}
function updateSelectionReadout() {
  const el = $('ro-selection'), lat = S.lat, st = S.strands;
  if (S.selLath >= 0 && st) {
    const s = st.strands[S.selLath], fam = st.family[S.selLath];
    let overs = 0; for (let i = 0; i < s.knots.length; i++) if (st.over(S.selLath, i)) overs++;
    el.innerHTML = `<span class="sw" style="background:${FAMILY_HEX[fam >= 0 ? fam % 3 : 3]}"></span>lath #${S.selLath} · family ${fam >= 0 ? 'ABC'[fam % 3] : '—'} · ${s.knots.length} knots (${overs} over / ${s.knots.length - overs} under) · ${s.closed ? 'closed' : 'open'}`;
  } else if (S.selVertex >= 0 && lat) {
    const i = S.selVertex, v = lat.valence[i];
    el.innerHTML = `<span class="sw" style="background:${valenceHex(v)}"></span>lattice vertex ${i} (OBJ v${i + 1}) · carrier vertex ${lat.peaks[i]} · valence ${v} (charge ${6 - v >= 0 ? '+' : ''}${6 - v}) · ${lat.interior[i] ? 'interior' : 'open link (front or rim)'}`;
  } else el.textContent = 'click a lath to follow its path · click a vertex to select it · Esc clears';
}

function renderValenceChips() {
  const box = $('valence-chips'), cs = S.summary, lat = S.lat;
  if (!lat || !cs.singular.length) { box.innerHTML = '<span class="hint">none yet</span>'; return; }
  const counts = {}; for (const i of cs.singular) { const v = lat.valence[i]; counts[v] = (counts[v] || 0) + 1; }
  const vals = Object.keys(counts).map(Number).sort((a, b) => a - b);
  const have = new Set(Array.from(box.querySelectorAll('.chip')).map(ch => +ch.dataset.valence));
  if (vals.length !== have.size || vals.some(v => !have.has(v))) {
    box.innerHTML = '';
    for (const v of vals) {
      const b = document.createElement('button'); b.className = 'chip' + (S.showValence[v] === false ? ' off' : ''); b.dataset.valence = v;
      b.innerHTML = `<span class="sw" style="background:${valenceHex(v)}"></span>v${v} · <span class="cnt">${counts[v]}</span>`;
      b.addEventListener('click', () => { const on = S.showValence[v] === false; S.showValence[v] = on; b.classList.toggle('off', !on); view.setValenceVisible(v, on); });
      box.appendChild(b);
    }
  } else for (const ch of box.querySelectorAll('.chip')) ch.querySelector('.cnt').textContent = counts[+ch.dataset.valence];
}

function updateLatticeReadout() {
  const lat = S.lat, cs = S.summary;
  if (!lat) { $('ro-lattice').textContent = 'no lattice yet'; $('ro-strands').textContent = ''; $('st-lat').textContent = ''; $('growth-text').textContent = 'no lattice yet — peaks · interior · singular are plotted against the step as the pattern grows'; renderValenceChips(); updateSelectionReadout(); updateGrowthLegend(); return; }
  const c = S.carrier, hist = Object.entries(cs.hist).sort((x, y) => +x[0] - +y[0]).map(([v, k]) => `v${v}: ${k}${+v !== 6 ? ` (charge ${6 - v >= 0 ? '+' : ''}${6 - v})` : ''}`).join('   ');
  $('ro-lattice').textContent = `peaks ${lat.peaks.length} · interior (closed link) ${cs.nInterior} · expected sites ≈ ${sites().toFixed(0)}\n${hist}\nΣ(6 − n) over interior peaks = ${cs.charge}    (Ω / 60° = ${(c.omega * 180 / Math.PI / 60).toFixed(2)})`;
  const st = S.strands;
  $('ro-strands').innerHTML = st ? `laths: ${st.strands.length} strands (${st.strands.filter(x => x.closed).length} closed) · families ${[0, 1, 2].map(f => `<span class="sw" style="background:${FAMILY_HEX[f]}"></span>${st.family.filter(x => x === f).length}`).join(' ')} · ${st.edges.length} knots · parity conflicts ${st.parityConflicts} · colour conflicts ${st.colourConflicts}` : '';
  $('st-lat').textContent = `· ${lat.peaks.length} peaks · ${cs.nInterior} interior · Σ(6−n) = ${cs.charge}`;
  renderValenceChips(); updateSelectionReadout(); updateGrowthLegend();
  const want = Math.round(c.omega * 180 / Math.PI / 60), ideal = cs.singular.length === 1 && cs.charge === want;
  const P = S.pred, i = P ? Math.min(P.floor.length - 1, Math.round((S.shownStep >= 0 ? S.shownStep : 0) / 50)) : -1, need = P && i >= 0 ? P.floor[i] : null;
  $('growth-text').innerHTML = `${cs.singular.length} singular · Σ(6−n) = ${cs.charge} over ${cs.nInterior} interior` +
    (need !== null ? ` · curvature inside the front needs ≥ ${need}` : '') +
    (ideal ? ` · <span class="ok">ideal: one singular peak, Σ(6−n) = round(Ω/60°) = ${want}</span>` : '') + predictionGap();
}

/** legend of the growth strip with the latest values (the chart itself has no room for it in the panel) */
function updateGrowthLegend() {
  const h = S.history.length ? S.history[S.history.length - 1] : null, step = S.shownStep >= 0 ? S.shownStep : 0;
  $('growth-legend').innerHTML = `<span class="lg c1">peaks${h ? ' ' + h[1] : ''}</span><span class="lg c2">interior${h ? ' ' + h[2] : ''}</span><span class="lg c3">singular${h ? ' ' + h[3] : ''} <small>(right axis)</small></span>` +
    `<span class="lg pred">predicted${S.pred && S.pred.coveredAt >= 0 ? ` · ideal cover ▴ ${S.pred.coveredAt.toLocaleString('en')}` : ''}</span><span class="dim">${step}${S.solver ? ' / ' + S.solver.steps : ''}</span>`;
}
/** tracked against the ideal front at the last read-out: the gap is what the surface costs the lattice */
function predictionGap() {
  const P = S.pred, H = S.history; if (!P || H.length < 2 || !(P.total > 10)) return '';
  const half = 0.5 * P.total;
  const cross = (xs, ys) => { for (let i = 1; i < xs.length; i++) if (ys[i] >= half && ys[i - 1] < half) return xs[i - 1] + (xs[i] - xs[i - 1]) * (half - ys[i - 1]) / Math.max(ys[i] - ys[i - 1], 1e-9); return -1; };
  const tp = cross(P.steps, P.peaks); if (tp <= 0) return '';
  const tt = cross(H.map(h => h[0]), H.map(h => h[1]));
  let d, txt;
  if (tt > 0) { d = tp / tt - 1; txt = `half covered at step ${Math.round(tt).toLocaleString('en')} (ideal ${Math.round(tp).toLocaleString('en')})`; }
  else {                                                                       // not there yet: compare the counts at the last read-out
    const [st, pk] = H[H.length - 1], i = Math.min(P.steps.length - 1, Math.round(st / 50)), pp = P.peaks[i];
    if (!(pp > 8)) return ''; d = (pk - pp) / pp; txt = `${pk} peaks at step ${st.toLocaleString('en')} (ideal ${pp.toFixed(0)})`;
  }
  const cls = Math.abs(d) < 0.12 ? 'ok' : 'warn';
  return ` · ${txt}: <span class="${cls}">front ${d >= 0 ? '+' : '−'}${Math.abs(100 * d).toFixed(0)} % against the ideal flat front</span>${d < -0.12 ? ' — strain or defects hold it back' : ''}`;
}

function updateCarrierColors() {
  const c = S.carrier; if (!c) return;
  const mode = S.carrierView, flat = () => view.setFlatColor(hexToRgb(getComputedStyle(document.documentElement).getPropertyValue('--mesh').trim()));
  const needGlobal = 'needs the global rotation field: 3 Singularity → global rotation field → Solve';
  if (mode !== 'stripes') view.setCarrierMode('vertex');
  if (mode === 'u') { S.shownStep = -1; updateCarrierUI(''); return; }
  if (mode === 'none') { flat(); updateCarrierUI(''); return; }
  if (mode === 'stripes') {
    const R = S.stripesRes;
    if (!R || R.carrier !== c || !R.faces) { flat(); updateCarrierUI(needGlobal); return; }
    if (S.stripesShown !== R.faces) { view.setStripes(R.faces); S.stripesShown = R.faces; }
    view.setCarrierMode('stripes');
    updateCarrierUI(`the seed's three stripe phases (before SH), drawn where each family's laths would run (band = lath width) · a band that stops is a branch point; two families stopping together make a 5–7 pair · around a field 5 or 7 the colours swap across a seam${Math.abs(R.a - S.a) > 1e-9 ? ` · <span class="warn">phases from a = ${R.a.toFixed(2)}</span>` : ''}`);
    return;
  }
  if (mode === 'rosysing') {
    const d = fieldChargeDensity();
    if (!d) { flat(); updateCarrierUI(needGlobal); return; }
    view.setFieldColors(d.vals, 'coolwarm', -d.v, d.v);
    const parts = Object.entries(S.rosy.hist).map(([v, k]) => `${k} × ${v}`).join(', ') || 'none';
    updateCarrierUI(`singularities of the rotation field (${parts}): red = 5 (+60°), blue = 7 (−60°); full colour at the singular face, fading within about two lattice edges (σ = 0.8 a along the surface)`);
    return;
  }
  const vals = new Float64Array(c.n);
  if (mode === 'K') {
    const b = c.base || c, kb = new Float64Array(b.n);                                   // defects sit on the base vertices; show their base-level density
    for (let i = 0; i < b.n; i++) kb[i] = b.boundary[i] ? 0 : b.K[i] / Math.max(b.m[i], 1e-300);
    vals.set(c.interp(kb));
  } else { const w = c.smoothedCurvature(Math.max(S.a, 2 * c.hBase)); vals.set(w); }
  const abs = Array.from(vals, Math.abs).sort((x, y) => x - y); const v = abs[Math.floor(abs.length * 0.98)] || 1;
  view.setFieldColors(vals, 'coolwarm', -v, v);
  updateCarrierUI('');
}
/** the rotation field's singular faces as charges (index: +1 for a 5, −1 for a 7, +2 for a 4 …), each spread as a
 *  Gaussian of the geodesic distance (σ = 0.8 a, bounded Dijkstra on the intrinsic edges from the face corners) and
 *  summed: one singularity is full colour at its face and fades within about two lattice edges. Cached per σ. */
function fieldChargeDensity() {
  const c = S.carrier, R = S.rosy; if (!c || !R || R.carrier !== c) return null;
  const sigma = 0.8 * S.a, key = sigma.toPrecision(3);
  if (R.density && R.density.key === key) return R.density;
  const b = c.base || c, G = c.graph, V = c.V, n = c.n, vals = new Float64Array(n), dist = new Float64Array(n).fill(Infinity), cut = 3 * sigma;
  for (const sg of R.list) {
    const ids = [b.Fi[3 * sg.face], b.Fi[3 * sg.face + 1], b.Fi[3 * sg.face + 2]], touched = [], heap = new MinHeap(256);   // base ids are kept by the refinement
    const cx = [0, 1, 2].map(k => (V[3 * ids[0] + k] + V[3 * ids[1] + k] + V[3 * ids[2] + k]) / 3);
    for (const i of ids) { const d0 = Math.hypot(V[3 * i] - cx[0], V[3 * i + 1] - cx[1], V[3 * i + 2] - cx[2]); if (d0 < dist[i]) { if (dist[i] === Infinity) touched.push(i); dist[i] = d0; heap.push(d0, i); } }
    while (heap.size) {
      const [d, i] = heap.pop(); if (d > dist[i]) continue;
      for (let k = G.indptr[i]; k < G.indptr[i + 1]; k++) {
        const j = G.indices[k], nd = d + G.data[k];
        if (nd < dist[j] && nd <= cut) { if (dist[j] === Infinity) touched.push(j); dist[j] = nd; heap.push(nd, j); }
      }
    }
    for (const i of touched) { vals[i] += sg.index * Math.exp(-0.5 * (dist[i] / sigma) ** 2); dist[i] = Infinity; }
  }
  R.density = { key, sigma, vals, v: 1 };
  return R.density;
}
/** opacity sliders: one for the vertex colourings, three (red, yellow, blue) for the stripes; plus the note line */
function updateCarrierUI(note) {
  const m = S.carrierView;
  $('row-carrier-alpha').hidden = m === 'none' || m === 'stripes';
  $('row-stripe-alpha').hidden = m !== 'stripes';
  const el = $('ro-carrier-view'); el.hidden = !note; el.innerHTML = note || '';
}
/** band half-width of the stripes view in radians of phase: the lath width w·a over the stripe wave number q0 = 4π/(√3 a) */
const stripeHalfWidth = () => Math.min(3.0, 2 * Math.PI * S.lathW / SQ3);
function hexToRgb(h) { const m = /^#?([0-9a-f]{6})$/i.exec(h); if (!m) return [0.85, 0.85, 0.85]; const n = parseInt(m[1], 16); return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255]; }

function updateDensityReadout() {
  const c = S.carrier; if (!c) return;
  const q0 = 4 * Math.PI / (SQ3 * S.a), epa = S.a / c.h, ok = epa >= 6 && epa <= 12;
  $('ro-q0').innerHTML = `q0 = 4π/(√3 a) = ${q0.toFixed(3)}`; void ok;
  const kq = 1 + (q0 * c.h) ** 2 / 48;                                           // selected k / q0 of the uncorrected operator (c = 1/24, see solver_cpu.js)
  $('ro-disp').textContent = `${epa.toFixed(1)} edges per a: uncorrected edge = ${(1 / kq).toFixed(3)} a` + (S.disp ? ` → corrected, edge = a` : ` (sites, pitch and strain are then ${(100 * (kq - 1)).toFixed(1)} % off)`);
  $('nb-sites').value = sites().toFixed(0);
  const [lo, hi] = aRange();
  $('ro-refine').textContent = `h = ${c.h.toFixed(3)} · ${c.n} v · a from ${lo.toFixed(2)} to ${hi.toFixed(2)} (${(c.area / (SQ3 / 2 * lo * lo)).toFixed(0)} to ${(c.area / (SQ3 / 2 * hi * hi)).toFixed(0)} sites)`;
  setSum('density', `a = ${S.a.toFixed(2)} · ≈ ${sites().toFixed(0)} sites${S.refine ? ` · ×${2 ** S.refine}` : ''}`);
  drawDensityDiagram();
}

function renderRefineChips() {
  const c = S.baseCarrier; if (!c) return;
  for (const ch of document.querySelectorAll('#refine-chips .chip')) {
    const L = +ch.dataset.level, n = Math.round(c.n * 4 ** L);
    ch.classList.toggle('on', L === S.refine);
    ch.title = `×${2 ** L} finer carrier: ≈ ${n >= 1e6 ? (n / 1e6).toFixed(2) + ' M' : (n / 1e3).toFixed(0) + ' k'} vertices, a down to ${(6 * c.h / 2 ** L).toFixed(3)}${L >= 3 ? ' — heavy: GPU only, read-outs take seconds' : L === 2 ? ' — needs the GPU for comfort' : ''}`;
    ch.disabled = n > 1.2e6;
  }
}
async function setRefine(level) {
  const base = S.baseCarrier; if (!base || level === S.refine) return;
  const n = Math.round(base.n * 4 ** level);
  setStatus(`refining ×${2 ** level} → ${n} vertices…`); await yieldTask(); await new Promise(r => requestAnimationFrame(r));
  const t0 = performance.now();
  try {
    const c = base.refined(level), keepValence = S.valence;                  // centre ids survive (base vertices keep their ids)
    setCarrier(c, Math.min(10 * c.h, Math.max(6 * c.h, S.a)), { refine: true });    // refining means denser: land at ≤ 10 edges per a
    if (S.valence !== keepValence) { S.valence = keepValence; $('sl-valence').value = S.valence; $('out-valence').value = S.valence; reseedNow(); }
    setStatus(`carrier ×${2 ** level}: ${c.n} v · h = ${c.h.toFixed(3)} · ${c.meta.flips} flips · ${((performance.now() - t0) / 1000).toFixed(1)} s`);
  } catch (e) { console.error(e); setStatus('refinement failed: ' + e.message); }
}
function renderSideChips() { for (const ch of document.querySelectorAll('#side-chips .chip')) ch.classList.toggle('on', +ch.dataset.side === S.side); }
function setSide(sign) {
  S.side = sign < 0 ? -1 : 1; view.setSide(S.side); renderSideChips();
  if (S.lat) { view.setLattice(S.lat, S.summary); view.setSelectedVertex(S.selVertex); updateLaths(); }
  if (S.carrier && S.centre >= 0) view.setCentre(S.centre);
}

function onDensityChange(a) {
  const [lo, hi] = aRange(); S.a = Math.min(hi, Math.max(lo, a));
  $('sl-a').value = sliderFromA(S.a); $('nb-a').value = S.a.toFixed(3);
  updateDensityReadout(); detect(); applyCgIters(); if (S.carrierView === 'Ks' || S.carrierView === 'rosysing' || S.carrierView === 'stripes') updateCarrierColors();   // field 5/7 spreads over σ = 0.8 a
  S.geo = null; S.solver.setParams({ a: S.a }); reseedNow();
  setStatus(`q0 = ${(4 * Math.PI / (SQ3 * S.a)).toFixed(3)} · ${sites().toFixed(0)} sites expected`);
}

function setCentre(idx, mode) {
  S.centre = idx; S.centreMode = mode; S.geo = null; updateCentreReadout();
  reseedNow();
  setStatus(`centre ${mode}: vertex ${idx} · defect ${(S.carrier.K[idx] * 180 / Math.PI).toFixed(1)}°`);
}

function syncStepsUI() { $('sl-steps').value = S.target; $('out-steps').value = S.target; }

function nearestCarrierVertex(p) {
  const c = S.carrier; if (!c) return -1; let best = -1, bd = Infinity;
  for (let i = 0; i < c.n; i++) { const d = (c.V[3 * i] - p.x) ** 2 + (c.V[3 * i + 1] - p.y) ** 2 + (c.V[3 * i + 2] - p.z) ** 2; if (d < bd) { bd = d; best = i; } }
  return best;
}

function setTurntable(on) {
  S.turntable = !!on; view.setTurntable(S.turntable, S.turnSpeed);
  $('nav-turn').classList.toggle('on', S.turntable); $('nav-turn').title = S.turntable ? `Turntable on · ${S.turnSpeed}°/s — click to stop (T)` : 'Turntable (T)';
  updateViewSum();
}
function setProjection(kind) {
  S.projection = kind; view.setProjection(kind);
  const ortho = kind === 'ortho';
  $('nav-proj').querySelector('.ic-persp').hidden = ortho; $('nav-proj').querySelector('.ic-ortho').hidden = !ortho;
  $('nav-proj').title = ortho ? 'Orthographic — click for perspective (P)' : 'Perspective — click for orthographic (P)';
  updateViewSum();
}
function updateViewSum() {
  const on = Array.from(document.querySelectorAll('#chips .chip.on')).map(ch => ch.dataset.layer);
  setSum('view', `${on.join(' · ')}${S.projection === 'ortho' ? ' · ortho' : ''}${S.turntable ? ' · turntable' : ''}`);
}

function setTheme(t) {
  S.theme = t;
  if (t === 'auto') document.documentElement.removeAttribute('data-theme'); else document.documentElement.setAttribute('data-theme', t);
  for (const chip of document.querySelectorAll('#theme-chips .chip')) chip.classList.toggle('on', chip.dataset.theme === t);
  try { localStorage.setItem('localsh-theme', t); } catch (e) { /* ignore */ }
}

// ------------------------------------------------------------------ layers shared by several chip rows
function setLayer(layer, on) {
  view.setVisible(layer, on);
  for (const ch of document.querySelectorAll(`.chip[data-layer="${layer}"]`)) ch.classList.toggle('on', on);
  if (layer === 'laths' && on) updateLaths();
  if (layer === 'orient' || layer === 'strain') updateDirectionLayers();
  if (layer === 'wire') updateMeshLayer();
  if (layer === 'probe') updateProbe();
  updateViewSum();
}
/** intrinsic Delaunay triangulation on the model (built when first shown: tracing the flipped edges takes ~0.1–2 s) */
function updateMeshLayer() {
  const c = S.carrier, el = $('ro-mesh');
  if (!c || !view.visible.wire) { view.setMesh(null); el.hidden = true; return; }
  const t0 = performance.now(), m = c.intrinsicEdgeLines(0);
  view.setMesh(m);
  const tot = m.nShared + m.nFlipped;
  el.hidden = false;
  el.innerHTML = `intrinsic Delaunay mesh: ${tot.toLocaleString('en')} edges, <span class="cm">${m.nFlipped} flipped (${(100 * m.nFlipped / Math.max(tot, 1)).toFixed(2)} %)</span> — each a straight line in the surface (geodesic), drawn across the faces it crosses${m.failed ? ` · ${m.failed} not traced (chords)` : ''} · the Laplacian uses these edges, so every interior cotan weight is ≥ 0 · ${(performance.now() - t0).toFixed(0)} ms`;
}
function updateProbe() {
  const c = S.carrier;
  if (!c || !view.visible.probe || S.centre < 0) { view.setProbe(null); return; }
  try { view.setProbe(c.seedProbe(S.centre, S.valence, S.a, S.seedPhi, S.seedKind)); } catch (e) { console.error(e); view.setProbe(null); }
}

// ------------------------------------------------------------------ growth prediction (ideal front, flat-lattice speed)
function geoFor() {
  const c = S.carrier; if (!c || S.centre < 0) return null;
  if (!S.geo || S.geo.centre !== S.centre || S.geo.carrier !== c) S.geo = { centre: S.centre, carrier: c, dist: c.geodesicDistance(S.centre), rim: c.rimDistance(), m: c.m, K: c.K, boundary: c.boundary };
  return S.geo;
}
function predXmax() { return Math.max(S.maxSteps, S.solver ? S.solver.steps : 0); }
function updatePrediction() {
  const geo = geoFor(); if (!geo) { S.pred = null; return; }
  try {
    S.pred = predictGrowth(geo, { a: S.a, seedR: S.seedR, seedKind: S.seedKind, valence: S.valence, r: S.r, b: S.b, dt: S.dt,
      anneal: { ...S.anneal, startedAt: S.anneal.active ? S.anneal.startedAt : -1 }, xmax: predXmax(), every: 50,
      siteScale: S.disp || !S.carrier ? 1 : 1 / (1 + (4 * Math.PI / (SQ3 * S.a) * S.carrier.h) ** 2 / 48) });
  } catch (e) { console.error(e); S.pred = null; }
  updateFrontReadout(); drawChart(); drawPhase(); drawDtChart();
}
function updateFrontReadout() {
  const f = frontState(S.r, S.b, S.dt), p = S.pred, A = S.anneal;
  const seedOk = 0.82 * S.seedR > f.RcA * 1.15;
  $('ro-seedr').innerHTML = !Number.isFinite(f.RcA) ? `<span class="warn">no growing front at (r, b) = (${S.r}, ${S.b}): ${f.label}</span>`
    : `critical radius R<sub>c</sub> = ${f.RcA.toFixed(2)} a at (r, b): a patch smaller than R<sub>c</sub> shrinks · this seed ≈ ${(0.82 * S.seedR).toFixed(2)} a <span class="${seedOk ? 'ok' : 'warn'}">${seedOk ? 'grows' : 'is marginal — enlarge it'}</span>`;
  let t = f.state === 'advance' ? `front ${f.vA100.toFixed(2)} a / 100 steps (flat, straight) · R<sub>c</sub> ${f.RcA.toFixed(2)} a` : `<span class="warn">${f.label}</span>`;
  if (f.state === 'advance' && S.r > 0) t += ` · <span class="warn">r > 0: u = 0 is unstable, patterns nucleate away from the front</span>`;
  if (p && p.coveredAt >= 0) t += ` · ideal front covers the surface by ≈ step ${p.coveredAt.toLocaleString('en')}`;
  else if (p && f.state === 'advance') t += ` · not covered within ${predXmax().toLocaleString('en')} steps`;
  if (A.on && p && p.ann) t += `\nanneal (r₂, b₂) = (${A.r2}, ${A.b2}): ${p.ann.state === 'advance' ? `front ${p.ann.vA100.toFixed(2)} a / 100 steps` : p.ann.label}${p.switchAt >= 0 && !A.active ? ` · switch predicted at ≈ step ${p.switchAt.toLocaleString('en')}` : ''}`;
  $('ro-front').innerHTML = t;
}

// ------------------------------------------------------------------ chart
const CHART_PAD = { l: 30, r: 26, t: 8, b: 16 };
const niceCeil = (x) => { const p = 10 ** Math.floor(Math.log10(Math.max(x, 1))), m = x / p; return ([1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10].find(k => m <= k + 1e-9)) * p; };
/** x range: what is computed or targeted, plus the predicted coverage (and anneal stage), not the whole 20 k budget */
function chartXmax() {
  const s = S.solver ? S.solver.steps : 0, p = S.pred;
  let want = Math.max(s, S.target, 1000);
  if (p && p.coveredAt >= 0) want = Math.max(want, p.coveredAt * 1.25, p.switchAt >= 0 ? p.switchAt + 0.5 * p.coveredAt : 0);
  return Math.min(Math.max(S.maxSteps, s), niceCeil(want * 1.02));
}
function drawChart() {
  const cv = $('chart'), ctx = cv.getContext('2d'), dpr = Math.min(window.devicePixelRatio || 1, 2);
  const W = cv.clientWidth || 360, H = cv.clientHeight || 118;
  if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(H * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); }
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
  const css = getComputedStyle(document.documentElement), v = (n) => css.getPropertyValue(n).trim();
  const dim = v('--dim'), border = v('--border'), warn = v('--warn'), cols = [v('--chart-1'), v('--chart-2'), v('--chart-3')];
  const Hh = S.history, pad = CHART_PAD, P = S.pred;
  const xmax = chartXmax();
  const inX = (st) => st <= xmax;
  const ymax = niceCeil(Math.max(10, ...Hh.map(h => h[1]), P ? P.total * 1.02 : 0));
  const ys = Math.max(4, ...Hh.map(h => h[3]), ...(P ? P.floor : [0])) + 1;            // right axis: singular count
  const X = (st) => pad.l + (W - pad.l - pad.r) * st / xmax, Y = (val) => H - pad.b - (H - pad.t - pad.b) * val / ymax, Ys = (val) => H - pad.b - (H - pad.t - pad.b) * val / ys;
  ctx.font = '10px IBM Plex Mono, monospace';
  if (S.solver && S.solver.steps > 0) { ctx.fillStyle = border; ctx.globalAlpha = 0.35; ctx.fillRect(X(0), pad.t, X(Math.min(S.solver.steps, xmax)) - X(0), H - pad.t - pad.b); ctx.globalAlpha = 1; }
  const annAt = S.anneal.active ? S.anneal.startedAt : (S.anneal.on && P && P.switchAt >= 0 ? P.switchAt : -1);
  if (annAt >= 0 && annAt <= xmax) {                                          // the anneal stage (tracked: tinted; predicted: dotted)
    const x0 = X(annAt), x1 = S.anneal.active ? X(Math.min(S.solver.steps, xmax)) : X(xmax);
    ctx.fillStyle = warn; ctx.globalAlpha = S.anneal.active ? 0.12 : 0.06; ctx.fillRect(x0, pad.t, Math.max(1, x1 - x0), H - pad.t - pad.b); ctx.globalAlpha = 1;
    ctx.strokeStyle = warn; ctx.setLineDash([2, 3]); ctx.beginPath(); ctx.moveTo(x0, pad.t); ctx.lineTo(x0, H - pad.b); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = warn; ctx.textAlign = 'left'; ctx.fillText('anneal', x0 + 3, pad.t + 9);
  }
  // axes
  ctx.strokeStyle = border; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, H - pad.b); ctx.lineTo(W - pad.r, H - pad.b); ctx.lineTo(W - pad.r, pad.t); ctx.stroke();
  ctx.fillStyle = dim; ctx.textAlign = 'right'; ctx.fillText(String(ymax), pad.l - 4, pad.t + 8); ctx.fillText('0', pad.l - 4, H - pad.b);
  ctx.fillStyle = cols[2]; ctx.textAlign = 'left'; ctx.fillText(String(ys), W - pad.r + 4, pad.t + 8); ctx.fillText('0', W - pad.r + 4, H - pad.b);
  ctx.fillStyle = dim; ctx.textAlign = 'center'; const kfmt = (st) => (st >= 1000 ? +(st / 1000).toFixed(2) + 'k' : String(Math.round(st)));
  ctx.fillText('0', X(0) + 4, H - 4); ctx.fillText(kfmt(Math.round(xmax / 2)), X(xmax / 2), H - 4); ctx.textAlign = 'right'; ctx.fillText(kfmt(xmax), X(xmax), H - 4);
  // prediction (dashed)
  if (P) {
    const line = (arr, yf, col, step) => {
      ctx.strokeStyle = col; ctx.lineWidth = 1.3; ctx.setLineDash([4, 3]); ctx.globalAlpha = 0.85; ctx.beginPath(); let first = true, prevY = 0;
      P.steps.forEach((st, i) => { if (!inX(st)) return; const x = X(st), y = yf(arr[i]); if (first) { ctx.moveTo(x, y); first = false; } else { if (step) ctx.lineTo(x, prevY); ctx.lineTo(x, y); } prevY = y; });
      ctx.stroke(); ctx.setLineDash([]); ctx.globalAlpha = 1;
    };
    line(P.peaks, Y, cols[0]); line(P.interior, Y, cols[1]); line(P.floor, Ys, cols[2], true);
    if (P.coveredAt >= 0 && P.coveredAt <= xmax) {                          // a small mark on the x axis; the number is in the legend
      const x = X(P.coveredAt); ctx.fillStyle = dim; ctx.beginPath(); ctx.moveTo(x, H - pad.b - 1); ctx.lineTo(x - 4, H - pad.b + 6); ctx.lineTo(x + 4, H - pad.b + 6); ctx.closePath(); ctx.fill();
    }
  }
  // tracked (solid)
  [0, 1, 2].forEach((k) => {
    ctx.strokeStyle = cols[k]; ctx.lineWidth = 1.6; ctx.beginPath();
    let first = true; Hh.forEach((h) => { if (!inX(h[0])) return; const x = X(h[0]), y = k === 2 ? Ys(h[3]) : Y(h[k + 1]); if (first) { ctx.moveTo(x, y); first = false; } else ctx.lineTo(x, y); }); ctx.stroke();
  });
  if (S.solver) { const x = X(Math.min(S.shownStep >= 0 ? S.shownStep : 0, xmax)); ctx.strokeStyle = dim; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, H - pad.b); ctx.stroke(); ctx.setLineDash([]); }
}

// ------------------------------------------------------------------ (r, b) map: front speed from the flat table
const PH = { r0: -0.5, r1: 0.1, b0: 0.5, b1: 2.0, pad: { l: 34, r: 8, t: 8, b: 24 } };
let phaseImg = null, phaseKey = '';
function drawPhase(live = null) {
  const cv = $('phase'); if (!cv) return;
  const ctx = cv.getContext('2d'), dpr = Math.min(window.devicePixelRatio || 1, 2);
  const W = cv.clientWidth || 372, H = Math.round(W * 0.5);
  if (cv.style.height !== H + 'px') cv.style.height = H + 'px';
  if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(H * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); phaseKey = ''; }
  const css = getComputedStyle(document.documentElement), cv_ = (n) => css.getPropertyValue(n).trim();
  const pad = PH.pad, pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
  const X = (b) => pad.l + pw * (b - PH.b0) / (PH.b1 - PH.b0), Y = (r) => pad.t + ph * (PH.r1 - r) / (PH.r1 - PH.r0);
  const toA = Math.sqrt(3) / (4 * Math.PI) * S.dt * 100;                        // q0-units → lattice edges per 100 steps
  const key = `${W}x${H}@${dpr}|${S.dt}|${cv_('--accent')}|${cv_('--panel-2')}`;
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  if (key !== phaseKey) {                                                       // the map itself, cached as an image
    const img = ctx.createImageData(cv.width, cv.height), rgb = (h) => { const m = /^#?([0-9a-f]{6})$/i.exec(h); const n = m ? parseInt(m[1], 16) : 0x888888; return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; };
    const lo = rgb(cv_('--panel-2')), hi = rgb(cv_('--accent')), pin = rgb(cv_('--warn')), bg = rgb(cv_('--panel')), ln = rgb(cv_('--border'));
    const levels = [0.25, 0.5, 0.75, 1.0];
    const spd = (px, py) => { const b = PH.b0 + (px / dpr - pad.l) / pw * (PH.b1 - PH.b0), r = PH.r1 - (py / dpr - pad.t) / ph * (PH.r1 - PH.r0); const v = mapValue(r, b); return v < -0.075 ? null : v * toA; };
    for (let py = 0; py < cv.height; py++) for (let px = 0; px < cv.width; px++) {
      const o = 4 * (py * cv.width + px), x = px / dpr, y = py / dpr;
      let c = bg;
      if (x >= pad.l && x <= pad.l + pw && y >= pad.t && y <= pad.t + ph) {
        const s = spd(px, py);
        if (s === null) { c = ((px + py) % Math.round(7 * dpr) < dpr) ? ln : bg; }                         // dies: hatched
        else if (s < 0.01 * toA) { const t = 0.32; c = [lo[0] + t * (pin[0] - lo[0]), lo[1] + t * (pin[1] - lo[1]), lo[2] + t * (pin[2] - lo[2])]; }   // pinned / retreating
        else { const t = Math.min(1, 0.12 + 0.88 * s / 1.6); c = [lo[0] + t * (hi[0] - lo[0]), lo[1] + t * (hi[1] - lo[1]), lo[2] + t * (hi[2] - lo[2])];
          const s2 = spd(px + dpr, py), s3 = spd(px, py + dpr);
          for (const L of levels) if ((s2 !== null && (s - L) * (s2 - L) < 0) || (s3 !== null && (s - L) * (s3 - L) < 0)) { c = [c[0] * 0.75, c[1] * 0.75, c[2] * 0.75]; break; } }
      }
      img.data[o] = c[0]; img.data[o + 1] = c[1]; img.data[o + 2] = c[2]; img.data[o + 3] = 255;
    }
    phaseImg = img; phaseKey = key;
  }
  ctx.putImageData(phaseImg, 0, 0);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const dim = cv_('--dim'), fg = cv_('--fg'), warn = cv_('--warn'), acc = cv_('--accent'), panel = cv_('--panel');
  ctx.font = '10px IBM Plex Mono, monospace';
  // r = 0: above it u = 0 is linearly unstable
  ctx.strokeStyle = warn; ctx.lineWidth = 1.2; ctx.setLineDash([4, 3]); ctx.beginPath(); ctx.moveTo(X(PH.b0), Y(0)); ctx.lineTo(X(PH.b1), Y(0)); ctx.stroke(); ctx.setLineDash([]);
  const tag = (t, x, y, col, align = 'left') => { ctx.textAlign = align; const w = ctx.measureText(t).width, x0 = align === 'right' ? x - w : x; ctx.fillStyle = panel; ctx.globalAlpha = 0.82; ctx.fillRect(x0 - 2, y - 9, w + 4, 12); ctx.globalAlpha = 1; ctx.fillStyle = col; ctx.fillText(t, x, y); };
  tag('r > 0: u = 0 unstable, patterns everywhere', X(PH.b0) + 4, Y(0) - 4, warn);
  // region labels
  tag('seed dies', X(0.56), Y(-0.36), dim);
  tag('pinned', X(1.98), Y(-0.44), dim, 'right');
  tag('front advances', X(1.98), Y(-0.15), fg, 'right');
  // axes
  ctx.strokeStyle = cv_('--border'); ctx.strokeRect(pad.l + 0.5, pad.t + 0.5, pw, ph);
  ctx.fillStyle = dim; ctx.textAlign = 'center';
  for (const b of [0.5, 1, 1.5, 2]) { ctx.fillText(b.toFixed(1), X(b), H - pad.b + 11); }
  ctx.fillText('b (quadratic)', pad.l + pw / 2, H - 2);
  ctx.textAlign = 'right'; for (const r of [-0.5, -0.3, -0.1, 0.1]) ctx.fillText(r.toFixed(1), pad.l - 4, Y(r) + 3);
  ctx.save(); ctx.translate(9, pad.t + ph / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; ctx.fillText('r', 0, 0); ctx.restore();
  // operating points
  const r = live && live.r !== undefined ? live.r : (S.anneal.active ? S.r : S.r), b = live && live.b !== undefined ? live.b : S.b;
  const pg = { x: X(b), y: Y(r) };
  if (S.anneal.on) {
    const pa = { x: X(S.anneal.b2), y: Y(S.anneal.r2) };
    ctx.strokeStyle = warn; ctx.lineWidth = 1.2; ctx.setLineDash([3, 2]); ctx.beginPath(); ctx.moveTo(pg.x, pg.y); ctx.lineTo(pa.x, pa.y); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = warn; ctx.strokeStyle = panel; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(pa.x, pa.y - 6); ctx.lineTo(pa.x + 6, pa.y); ctx.lineTo(pa.x, pa.y + 6); ctx.lineTo(pa.x - 6, pa.y); ctx.closePath(); ctx.stroke(); ctx.fill();
    ctx.textAlign = pa.x < pad.l + 60 ? 'left' : 'right'; ctx.fillText('anneal', pa.x + (pa.x < pad.l + 60 ? 9 : -9), pa.y + 3);
  }
  ctx.fillStyle = acc; ctx.strokeStyle = panel; ctx.lineWidth = 2.5; ctx.beginPath(); ctx.arc(pg.x, pg.y, 5.5, 0, 2 * Math.PI); ctx.stroke(); ctx.fill();
  const f = frontState(r, b, S.dt);
  const lab = f.state === 'advance' ? `${f.vA100.toFixed(2)} a/100 steps` : f.state;
  tag(lab, pg.x + (pg.x > pad.l + pw - 110 ? -9 : 9), pg.y - 8, fg, pg.x > pad.l + pw - 110 ? 'right' : 'left');
  // colour key
  tag('iso-lines 0.25 · 0.5 · 0.75 · 1 a / 100 steps', pad.l + 4, pad.t + ph - 4, dim);
}
function bindPhase() {
  const cv = $('phase'); let drag = null;
  const at = (e) => { const r = cv.getBoundingClientRect(), pad = PH.pad, pw = r.width - pad.l - pad.r, ph = r.width * 0.5 - pad.t - pad.b;
    const b = PH.b0 + (e.clientX - r.left - pad.l) / pw * (PH.b1 - PH.b0), rr = PH.r1 - (e.clientY - r.top - pad.t) / ph * (PH.r1 - PH.r0);
    return { b: Math.round(Math.max(0.5, Math.min(2, b)) * 100) / 100, r: Math.round(Math.max(-0.5, Math.min(0.1, rr)) * 200) / 200 }; };
  const near = (e, b, r) => { const rect = cv.getBoundingClientRect(), pad = PH.pad, pw = rect.width - pad.l - pad.r, ph = rect.width * 0.5 - pad.t - pad.b;
    const x = rect.left + pad.l + pw * (b - PH.b0) / (PH.b1 - PH.b0), y = rect.top + pad.t + ph * (PH.r1 - r) / (PH.r1 - PH.r0); return Math.hypot(e.clientX - x, e.clientY - y); };
  const apply = (which, p, final) => {
    if (which === 'anneal') {
      S.anneal.r2 = p.r; S.anneal.b2 = p.b; $('nb-r2').value = p.r; $('nb-b2').value = p.b;
      if (final) { if (S.anneal.active && S.solver) S.solver.setParams({ r: p.r, b: p.b }); renderAnnealReadout(); updatePrediction(); } else drawPhase();
      return;
    }
    $('sl-r').value = p.r; $('out-r').value = p.r.toFixed(3); $('sl-b').value = p.b; $('out-b').value = p.b.toFixed(2);
    if (!final) { drawPhase({ r: p.r, b: p.b }); return; }
    if (S.anneal.active) { S.anneal.r2 = p.r; S.anneal.b2 = p.b; $('nb-r2').value = p.r; $('nb-b2').value = p.b; } else { S.r = p.r; S.b = p.b; }
    if (S.solver) S.solver.setParams({ r: p.r, b: p.b }); setStatus(`r = ${p.r}, b = ${p.b} (operator rebuilt, field kept)`); renderAnnealReadout(); updatePrediction();
  };
  cv.addEventListener('pointerdown', (e) => {
    const p = at(e); drag = S.anneal.on && near(e, S.anneal.b2, S.anneal.r2) < near(e, S.b, S.r) ? 'anneal' : 'grow';
    cv.setPointerCapture(e.pointerId); apply(drag, p, false);
  });
  cv.addEventListener('pointermove', (e) => { if (drag) apply(drag, at(e), false); });
  const end = (e) => { if (!drag) return; apply(drag, at(e), true); drag = null; };
  cv.addEventListener('pointerup', end); cv.addEventListener('pointercancel', () => { drag = null; drawPhase(); });
}

// ------------------------------------------------------------------ dt: steps to cover against the front-speed error
function drawDtChart(dtLive = null) {
  const cv = $('dtchart'); if (!cv) return;
  const ctx = cv.getContext('2d'), dpr = Math.min(window.devicePixelRatio || 1, 2);
  const W = cv.clientWidth || 372, H = 78;
  if (cv.width !== Math.round(W * dpr) || cv.height !== Math.round(H * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); cv.style.height = H + 'px'; }
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
  const css = getComputedStyle(document.documentElement), cv_ = (n) => css.getPropertyValue(n).trim();
  const dim = cv_('--dim'), acc = cv_('--accent'), warn = cv_('--warn'), border = cv_('--border'), fg = cv_('--fg');
  const pad = { l: 34, r: 34, t: 8, b: 18 }, pw = W - pad.l - pad.r, ph = H - pad.t - pad.b, d0 = 0.02, d1 = 0.3;
  const X = (d) => pad.l + pw * (d - d0) / (d1 - d0);
  const rel = (d) => { const xs = FRONT.dt, ys = FRONT.dtRel[0]; for (let i = 1; i < xs.length; i++) if (d <= xs[i]) { const t = (d - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[i - 1] + t * (ys[i] - ys[i - 1]); } return ys[ys.length - 1]; };
  const dt = dtLive ?? S.dt, P = S.pred;
  // steps to cover ∝ 1 / (dt · rel(dt)) — scaled to the current prediction, or steps per 10 a of front travel
  const f = frontState(S.r, S.b, 0.1);
  let ref = null, what = '';
  if (P && P.coveredAt > 0) { ref = (d) => P.coveredAt * (S.dt * rel(S.dt)) / (d * rel(d)); what = 'steps to cover'; }
  else if (f.state === 'advance') { ref = (d) => 10 / (f.vA100 / 0.1 * d * rel(d) / rel(0.1)) * 100; what = 'steps per 10 a'; }
  ctx.font = '10px IBM Plex Mono, monospace';
  ctx.strokeStyle = border; ctx.strokeRect(pad.l + 0.5, pad.t + 0.5, pw, ph);
  if (ref) {
    const ymax = niceCeil(ref(d0) * 0.5) * 2 > ref(d0) ? niceCeil(ref(d0)) : niceCeil(ref(d0));
    const Y = (v) => pad.t + ph * (1 - Math.min(v, ymax) / ymax);
    ctx.strokeStyle = acc; ctx.lineWidth = 1.6; ctx.beginPath(); for (let i = 0; i <= 60; i++) { const d = d0 + (d1 - d0) * i / 60; const y = Y(ref(d)); if (i) ctx.lineTo(X(d), y); else ctx.moveTo(X(d), y); } ctx.stroke();
    ctx.fillStyle = acc; ctx.textAlign = 'right'; ctx.fillText(String(ymax >= 1000 ? (ymax / 1000) + 'k' : ymax), pad.l - 4, pad.t + 8); ctx.fillText('0', pad.l - 4, pad.t + ph);
    ctx.textAlign = 'left'; ctx.fillText(what, X(0.16), Math.min(pad.t + ph - 4, Y(ref(0.16)) - 6));
    var stepsNow = Math.round(ref(dt));
  }
  // error (right axis 0 … 5 %)
  const emax = 5, Ye = (e) => pad.t + ph * (1 - e / emax);
  ctx.strokeStyle = warn; ctx.lineWidth = 1.3; ctx.setLineDash([4, 3]); ctx.beginPath(); for (let i = 0; i <= 60; i++) { const d = d0 + (d1 - d0) * i / 60, e = 100 * (1 - rel(d)); if (i) ctx.lineTo(X(d), Ye(e)); else ctx.moveTo(X(d), Ye(e)); } ctx.stroke(); ctx.setLineDash([]);
  ctx.fillStyle = warn; ctx.textAlign = 'left'; ctx.fillText('5 %', W - pad.r + 4, pad.t + 8); ctx.fillText('0', W - pad.r + 4, pad.t + ph);
  ctx.textAlign = 'right'; ctx.fillText('front slower by', X(0.29), Math.min(pad.t + ph - 3, Ye(100 * (1 - rel(0.29))) + 13));
  ctx.fillStyle = dim; ctx.textAlign = 'center'; for (const d of [0.05, 0.1, 0.2, 0.3]) ctx.fillText(String(d), X(d), H - 5);
  ctx.strokeStyle = fg; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(X(dt), pad.t); ctx.lineTo(X(dt), pad.t + ph); ctx.stroke();
  const right = X(dt) < pad.l + pw * 0.6, tx = X(dt) + (right ? 5 : -5);
  ctx.textAlign = right ? 'left' : 'right'; ctx.fillStyle = fg;
  ctx.fillText(`dt ${dt.toFixed(2)}: ${typeof stepsNow === 'number' ? stepsNow.toLocaleString('en') + ' steps · ' : ''}front −${(100 * (1 - rel(dt))).toFixed(1)} %`, tx, pad.t + 10);
}

// ------------------------------------------------------------------ density: lattice edge, pitch, sites and carrier resolution at true ratio
function svgEl(tag, attrs, text) { const e = document.createElementNS('http://www.w3.org/2000/svg', tag); for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v); if (text != null) e.textContent = text; return e; }
function drawDensityDiagram() {
  const svg = $('dens-svg'), c = S.carrier; if (!svg || !c) return;
  svg.innerHTML = '';
  const css = getComputedStyle(document.documentElement), cv_ = (n) => css.getPropertyValue(n).trim();
  const acc = cv_('--accent'), knot = cv_('--knot'), edge = cv_('--mesh-edge'), mesh = cv_('--mesh'), dim = cv_('--dim'), warn = cv_('--warn'), ok = cv_('--ok'), border = cv_('--border');
  const epa = S.a / c.h, A = 58, cx = 92, cy = 78, h = A / epa;
  const defs = svgEl('defs', {}); const clip = svgEl('clipPath', { id: 'dclip' }); clip.appendChild(svgEl('circle', { cx, cy, r: 74 })); defs.appendChild(clip); svg.appendChild(defs);
  const g = svgEl('g', { 'clip-path': 'url(#dclip)' }); svg.appendChild(g);
  g.appendChild(svgEl('circle', { cx, cy, r: 74, fill: mesh, opacity: 0.55 }));
  // carrier mesh (equilateral, edge h) at the true ratio to a
  let d = ''; const R = 80, s3 = Math.sqrt(3) / 2;
  for (let k = -Math.ceil(R / (h * s3)); k <= Math.ceil(R / (h * s3)); k++) { const y = cy + k * h * s3; d += `M${cx - R} ${y.toFixed(2)}H${cx + R}`; }
  for (const ang of [60, 120]) { const t = ang * Math.PI / 180, ux = Math.cos(t), uy = Math.sin(t), nx = -uy, ny = ux;
    for (let k = -Math.ceil(R / (h * s3)); k <= Math.ceil(R / (h * s3)); k++) { const ox = cx + nx * k * h * s3, oy = cy + ny * k * h * s3; d += `M${(ox - ux * R).toFixed(2)} ${(oy - uy * R).toFixed(2)}L${(ox + ux * R).toFixed(2)} ${(oy + uy * R).toFixed(2)}`; } }
  g.appendChild(svgEl('path', { d, stroke: edge, 'stroke-width': 0.6, fill: 'none', opacity: 0.9 }));
  // lattice: centre + 6 neighbours + the next ring, edges
  const pts = []; for (let i = -3; i <= 3; i++) for (let j = -3; j <= 3; j++) { const x = cx + A * (i + j / 2), y = cy - A * j * s3; if (Math.hypot(x - cx, y - cy) < 100) pts.push([x, y]); }
  let le = ''; for (const p of pts) for (const q of pts) { const dd = Math.hypot(p[0] - q[0], p[1] - q[1]); if (dd > A * 0.99 && dd < A * 1.01 && p[0] + p[1] * 1e-3 < q[0] + q[1] * 1e-3) le += `M${p[0].toFixed(1)} ${p[1].toFixed(1)}L${q[0].toFixed(1)} ${q[1].toFixed(1)}`; }
  // two laths of one family (horizontal strands through the edge midpoints), width = lath slider × a
  const lw = Math.max(1.5, S.lathW * A);
  for (const yy of [cy - A * s3 / 2, cy + A * s3 / 2]) g.appendChild(svgEl('line', { x1: cx - 90, y1: yy, x2: cx + 90, y2: yy, stroke: FAMILY_HEX[0], 'stroke-width': lw, opacity: 0.28 }));
  g.appendChild(svgEl('path', { d: le, stroke: knot, 'stroke-width': 1.2, fill: 'none' }));
  for (const p of pts) g.appendChild(svgEl('circle', { cx: p[0], cy: p[1], r: 3, fill: knot }));
  svg.appendChild(svgEl('circle', { cx, cy, r: 74, fill: 'none', stroke: border }));
  // dimension: a along the centre → right neighbour, ticks at every carrier edge
  const ya = cy + 0; svg.appendChild(svgEl('line', { x1: cx, y1: ya, x2: cx + A, y2: ya, stroke: acc, 'stroke-width': 2.2 }));
  for (let k = 0; k <= Math.floor(epa + 1e-9); k++) svg.appendChild(svgEl('line', { x1: cx + k * h, y1: ya - 3, x2: cx + k * h, y2: ya + 3, stroke: acc, 'stroke-width': 1 }));
  svg.appendChild(svgEl('text', { x: cx + A / 2, y: ya - 6, 'text-anchor': 'middle', class: 'acc big' }, 'a'));
  // pitch between the two laths
  const px = cx - A * 0.75; svg.appendChild(svgEl('line', { x1: px, y1: cy - A * s3 / 2, x2: px, y2: cy + A * s3 / 2, stroke: FAMILY_HEX[0], 'stroke-width': 1.2, 'marker-start': '', 'marker-end': '' }));
  for (const yy of [cy - A * s3 / 2, cy + A * s3 / 2]) svg.appendChild(svgEl('line', { x1: px - 4, y1: yy, x2: px + 4, y2: yy, stroke: FAMILY_HEX[0], 'stroke-width': 1.2 }));
  svg.appendChild(svgEl('text', { x: px - 5, y: cy - 6, 'text-anchor': 'end', class: 'fg' }, 'p'));
  // numbers on the right
  const site = SQ3 / 2 * S.a * S.a, n = sites(), x0 = 186, lines = [
    ['acc big', `a = ${S.a.toFixed(3)}`],
    ['fg', `p = √3⁄2·a = ${(SQ3 / 2 * S.a).toFixed(3)}`],
    ['', `lath w = ${S.lathW.toFixed(2)}·a = ${(S.lathW * S.a).toFixed(3)}`],
    ['', `sites ≈ area / (√3⁄2·a²)`],
    ['fg', `  = ${c.area.toFixed(1)} / ${site.toFixed(2)} ≈ ${n.toFixed(0)}`],
    ['', `crossings ≈ 3 × sites = ${(3 * n).toFixed(0)}`],
  ];
  lines.forEach(([cls, t], i) => svg.appendChild(svgEl('text', { x: x0, y: 14 + i * 15, class: cls }, t)));
  // resolution gauge: carrier edges per a, 4 … 16, ok 6 … 12
  const gx = x0, gy = 118, gw = 180, G = (e) => gx + gw * (Math.min(16, Math.max(4, e)) - 4) / 12;
  svg.appendChild(svgEl('rect', { x: gx, y: gy, width: gw, height: 7, fill: warn, opacity: 0.25, rx: 2 }));
  svg.appendChild(svgEl('rect', { x: G(6), y: gy, width: G(12) - G(6), height: 7, fill: ok, opacity: 0.35 }));
  for (const t of [4, 6, 8, 10, 12, 14, 16]) svg.appendChild(svgEl('text', { x: G(t), y: gy + 18, 'text-anchor': 'middle' }, String(t)));
  const good = epa >= 6 && epa <= 12;
  svg.appendChild(svgEl('path', { d: `M${G(epa)} ${gy - 1}l-4 -6h8z`, fill: good ? ok : warn }));
  svg.appendChild(svgEl('text', { x: gx, y: gy - 9, class: good ? 'ok' : 'warn' }, `${epa.toFixed(1)} carrier edges per a${good ? '' : epa < 6 ? ' — refine' : ' — slow, coarser a ok'}`));
  svg.appendChild(svgEl('text', { x: gx + gw, y: gy + 30, 'text-anchor': 'end' }, `h = ${c.h.toFixed(3)}  (mesh drawn to scale)`));
}

// ------------------------------------------------------------------ singularity: detected curvature against the valence
function drawChargeDiagram() {
  const svg = $('charge-svg'), c = S.carrier, d = S.detected; if (!svg || !c || !d) return;
  svg.innerHTML = '';
  const css = getComputedStyle(document.documentElement), cv_ = (n) => css.getPropertyValue(n).trim();
  const acc = cv_('--accent'), fg = cv_('--fg'), border = cv_('--border');
  const x0 = 92, x1 = 360, V = (n) => x0 + (x1 - x0) * (n - 3) / 6, y = 44;
  const om = d.info.omegaDeg, nEff = 6 - om / 60;
  svg.appendChild(svgEl('text', { x: 372, y: 11, 'text-anchor': 'end', class: 'acc' }, `Ω = ∫K = ${om.toFixed(0)}° → 6 − Ω/60° = ${nEff.toFixed(2)} → v${d.valence}`));
  svg.appendChild(svgEl('line', { x1: V(3) - 12, y1: y, x2: V(9) + 12, y2: y, stroke: border, 'stroke-width': 1.5 }));
  for (const [t, yy] of [['valence n', y + 4], ['charge 6−n', y + 22], ['cone n·60°', y + 36]]) svg.appendChild(svgEl('text', { x: 0, y: yy }, t));
  for (let n = 3; n <= 9; n++) {
    const on = n === S.valence;
    svg.appendChild(svgEl('circle', { cx: V(n), cy: y, r: on ? 10 : 8, fill: valenceHex(n), stroke: on ? fg : border, 'stroke-width': on ? 2 : 1, opacity: n === 6 ? 0.6 : 1 }));
    svg.appendChild(svgEl('text', { x: V(n), y: y + 3.5, 'text-anchor': 'middle', style: 'font-weight:600; fill:#1b1d1f' }, String(n)));
    svg.appendChild(svgEl('text', { x: V(n), y: y + 22, 'text-anchor': 'middle', class: on ? 'fg' : '' }, `${6 - n > 0 ? '+' : ''}${6 - n}`));
    svg.appendChild(svgEl('text', { x: V(n), y: y + 36, 'text-anchor': 'middle', class: on ? 'fg' : '' }, `${n * 60}°`));
  }
  // the detected curvature, as a position on the same axis
  const xm = V(Math.max(2.6, Math.min(9.4, nEff)));
  svg.appendChild(svgEl('path', { d: `M${xm} ${y - 12}l-5 -9h10z`, fill: acc }));
  svg.appendChild(svgEl('line', { x1: xm, y1: y - 12, x2: xm, y2: y + 12, stroke: acc, 'stroke-width': 1, 'stroke-dasharray': '2 2' }));
}
function drawOrientLegend() {
  const svg = $('orient-legend'); if (!svg) return; svg.innerHTML = '';
  const n = 24, w = 120, base = hexToRgb(getComputedStyle(document.documentElement).getPropertyValue('--knot').trim());
  for (let i = 0; i < n; i++) { const t = i / (n - 1), rgb = S.orientMode === 'angle' ? lutColor('twilight', t) : misorRGB(Math.max(0, (15 * t - 2) / 13), base); svg.appendChild(svgEl('rect', { x: 2 + i * w / n, y: 2, width: w / n + 0.6, height: 8, fill: `rgb(${(rgb[0] * 255) | 0},${(rgb[1] * 255) | 0},${(rgb[2] * 255) | 0})` })); }
  const [l, r] = S.orientMode === 'angle' ? ['−30°', '+30°'] : ['0°', '15°+'];
  svg.appendChild(svgEl('text', { x: 2, y: 20 }, l)); svg.appendChild(svgEl('text', { x: 2 + w, y: 20, 'text-anchor': 'end' }, r));
}

// ------------------------------------------------------------------ hover cards (every ⓘ, and data-tip="<template id>")
function initTips() {
  const card = document.createElement('div'); card.id = 'tipcard'; card.setAttribute('role', 'tooltip'); card.hidden = true; document.body.appendChild(card);
  let timer = 0, owner = null;
  const place = (el) => {
    const r = el.getBoundingClientRect(), cw = card.offsetWidth, ch = card.offsetHeight;
    let x = Math.min(r.left, window.innerWidth - 8 - cw), y = r.bottom + 6;
    if (y + ch > window.innerHeight - 8) y = Math.max(8, r.top - 6 - ch);
    card.style.left = Math.max(8, x) + 'px'; card.style.top = y + 'px';
  };
  const show = (el) => {
    const tpl = el.dataset.tip ? document.getElementById(el.dataset.tip) : null;
    card.innerHTML = '';
    if (tpl) card.appendChild(tpl.content.cloneNode(true)); else card.textContent = el.dataset.tipText || '';
    for (const b of card.querySelectorAll('[data-tip-range]')) b.textContent = Math.round(100 * S.strainRange);
    card.hidden = false; owner = el; place(el);
  };
  const hide = () => { clearTimeout(timer); card.hidden = true; owner = null; };
  for (const el of document.querySelectorAll('.info, [data-tip]')) {
    if (el.title) { el.dataset.tipText = el.title; el.removeAttribute('title'); }
    if (el.classList.contains('info') && !el.getAttribute('aria-label')) el.setAttribute('aria-label', 'more information');
    const delay = el.classList.contains('info') ? 60 : 450;                  // a chip only after a pause, so clicking it does not flash the card
    el.addEventListener('mouseenter', () => { clearTimeout(timer); timer = setTimeout(() => show(el), delay); });
    el.addEventListener('mouseleave', hide);
    el.addEventListener('focus', () => { if (el.classList.contains('info')) show(el); });
    el.addEventListener('blur', hide);
    if (el.classList.contains('info')) el.addEventListener('click', (e) => { e.preventDefault(); if (owner === el) hide(); else show(el); });   // touch
  }
  window.addEventListener('keydown', (e) => { if (e.key === 'Escape') hide(); });
  $('panel').addEventListener('scroll', hide, { passive: true });
}

/** click / drag on the growth strip scrubs the shown step (same as the steps slider) */
function chartStepAt(clientX) {
  const cv = $('chart'), r = cv.getBoundingClientRect(), W = r.width, pad = CHART_PAD;
  const s = chartXmax() * (clientX - r.left - pad.l) / Math.max(1, W - pad.l - pad.r);
  return Math.max(0, Math.min(S.maxSteps, Math.round(s / 10) * 10));
}

// ------------------------------------------------------------------ export
function exportLattice() {
  const lat = S.lat, c = S.carrier; if (!lat) { setStatus('nothing to export yet'); return; }
  const lines = [`# LocalSH ${c.name}: valence ${S.valence} seed, a = ${S.a.toFixed(4)}, ${S.solver.steps} steps, r = ${S.r}, b = ${S.b}${S.anneal.active ? `, annealed from step ${S.anneal.startedAt} at r = ${S.anneal.r2}, b = ${S.anneal.b2}` : ''}, sum(6-n) interior = ${S.summary.charge}`];
  const P = lat.P; for (let i = 0; i < P.length / 3; i++) lines.push(`v ${P[3 * i].toFixed(6)} ${P[3 * i + 1].toFixed(6)} ${P[3 * i + 2].toFixed(6)}`);
  const T = lat.tris; for (let t = 0; t < T.length / 3; t++) lines.push(`f ${T[3 * t] + 1} ${T[3 * t + 1] + 1} ${T[3 * t + 2] + 1}`);
  const sing = S.summary.singular.map(i => `# singular ${i}: valence ${lat.valence[i]} charge ${6 - lat.valence[i]}`);
  const txt = lines.concat(sing).join('\n') + '\n';
  $('modal-title').textContent = `lattice OBJ — ${P.length / 3} vertices, ${T.length / 3} triangles`;
  $('modal-hint').textContent = 'Copy the text into a .obj file (the dual triangulation of the peaks; singular vertices listed as comments at the end), or use Download where the page allows it.';
  $('modal-text').value = txt;
  try { const a = $('modal-download'); if (a._url) URL.revokeObjectURL(a._url); a._url = URL.createObjectURL(new Blob([txt], { type: 'text/plain' })); a.href = a._url; a.download = `${c.name.replace(/[^A-Za-z0-9_-]+/g, '_')}_v${S.valence}_${P.length / 3}v.obj`; } catch (e) { /* ignore */ }
  $('modal').hidden = false;
}

// ------------------------------------------------------------------ bindings
function bind() {
  $('sel-zoo').addEventListener('change', (e) => { if (e.target.value !== '__file__') loadCarrierByName(e.target.value); });
  $('file-stl').addEventListener('change', (e) => { const f = e.target.files[0]; if (f) openFile(f); });
  const vp = $('viewport');
  vp.addEventListener('dragover', (e) => { e.preventDefault(); $('drop').hidden = false; });
  vp.addEventListener('dragleave', () => { $('drop').hidden = true; });
  vp.addEventListener('drop', (e) => { e.preventDefault(); $('drop').hidden = true; const f = e.dataTransfer.files[0]; if (f) openFile(f); });
  // density
  $('sl-a').addEventListener('input', (e) => { $('nb-a').value = aFromSlider(+e.target.value).toFixed(3); });
  $('sl-a').addEventListener('change', (e) => onDensityChange(aFromSlider(+e.target.value)));
  $('nb-a').addEventListener('change', (e) => onDensityChange(+e.target.value));
  $('nb-sites').addEventListener('change', (e) => { const n = +e.target.value; if (n > 3 && S.carrier) onDensityChange(Math.sqrt(S.carrier.area / (SQ3 / 2 * n))); });
  // singularity
  $('sl-valence').addEventListener('input', (e) => { $('out-valence').value = e.target.value; });
  $('sl-valence').addEventListener('change', (e) => { S.valence = +e.target.value; reseedNow(); drawChargeDiagram(); });
  $('btn-detected').addEventListener('click', () => { if (S.detected) { S.valence = Math.max(3, Math.min(9, S.detected.valence)); $('sl-valence').value = S.valence; $('out-valence').value = S.valence; setCentre(S.detected.centre, 'detected'); } });
  $('sel-seed').addEventListener('change', (e) => { S.seedKind = e.target.value; reseedNow(); });
  $('sl-seedr').addEventListener('input', (e) => { $('out-seedr').value = (+e.target.value).toFixed(1); });
  $('sl-seedr').addEventListener('change', (e) => { S.seedR = +e.target.value; reseedNow(); updateFrontReadout(); });
  $('sl-seedphi').addEventListener('input', (e) => { $('out-seedphi').value = e.target.value; });
  $('sl-seedphi').addEventListener('change', (e) => { S.seedPhi = +e.target.value; reseedNow(); });
  for (const chip of document.querySelectorAll('#seed-mode .chip')) chip.addEventListener('click', () => {
    S.seedMode = chip.dataset.mode;
    for (const ch of document.querySelectorAll('#seed-mode .chip')) ch.classList.toggle('on', ch === chip);
    $('global-block').hidden = S.seedMode !== 'global'; $('local-block').hidden = S.seedMode === 'global';
    if (S.seedMode !== 'global') view.setRosy(null);
    reseedNow();
  });
  $('btn-rosy').addEventListener('click', () => { if (S.seedMode !== 'global') { document.querySelector('#seed-mode .chip[data-mode="global"]').click(); } else globalSeed(true); });
  $('cb-disp').addEventListener('change', (e) => { S.disp = e.target.checked; if (S.solver) S.solver.setParams({ dispersion: S.disp }); updateDensityReadout(); updatePrediction(); setStatus(S.disp ? 'operator corrected for the carrier dispersion (edge = a), field kept' : 'operator uncorrected (edge ≈ a / (1 + (q0 h)²/48)), field kept'); });
  // solve
  $('sl-steps').addEventListener('input', (e) => { S.target = +e.target.value; $('out-steps').value = S.target; });
  $('btn-run').addEventListener('click', () => { S.playing = !S.playing; $('btn-run').textContent = S.playing ? 'Pause' : 'Run'; if (S.playing && S.solver && S.target <= S.solver.steps) { S.target = Math.min(S.maxSteps, S.solver.steps + S.spf); syncStepsUI(); } });
  $('btn-plus').addEventListener('click', () => { if (S.solver) { S.target = Math.min(S.maxSteps, S.solver.steps + 500); syncStepsUI(); } });
  $('btn-reseed').addEventListener('click', () => { S.playing = false; $('btn-run').textContent = 'Run'; reseedNow(); });
  $('btn-export').addEventListener('click', exportLattice);
  const paramSlider = (id, out, key, digits) => {
    $(id).addEventListener('input', (e) => { $(out).value = (+e.target.value).toFixed(digits); drawPhase(key === 'dt' ? null : { [key]: +e.target.value }); if (key === 'dt') drawDtChart(+e.target.value); });
    $(id).addEventListener('change', (e) => {
      const v = +e.target.value;
      if (S.anneal.active && (key === 'r' || key === 'b')) { S.anneal[key + '2'] = v; $('nb-' + key + '2').value = v; }   // while annealing the sliders drive the anneal stage
      else S[key] = v;
      if (S.solver) S.solver.setParams({ [key]: v }); setStatus(`${key} = ${v} (operator rebuilt, field kept)`); renderAnnealReadout(); updatePrediction();
    });
  };
  paramSlider('sl-r', 'out-r', 'r', 3); paramSlider('sl-b', 'out-b', 'b', 2); paramSlider('sl-dt', 'out-dt', 'dt', 2);
  $('nb-max').addEventListener('change', (e) => { S.maxSteps = Math.max(500, +e.target.value | 0); $('sl-steps').max = S.maxSteps; S.target = Math.min(S.target, S.maxSteps); syncStepsUI(); drawChart(); });
  $('nb-spf').addEventListener('change', (e) => { S.spf = Math.max(1, Math.min(200, +e.target.value | 0)); S.spfAuto = 0; });
  $('nb-cg').addEventListener('change', (e) => { S.cgManual = true; S.cgIters = Math.max(4, Math.min(60, +e.target.value | 0)); if (S.solver) S.solver.cgIters = S.cgIters; });
  for (const ch of document.querySelectorAll('#refine-chips .chip')) ch.addEventListener('click', () => setRefine(+ch.dataset.level));
  // anneal
  $('ck-anneal').addEventListener('change', (e) => { S.anneal.on = e.target.checked; if (!S.anneal.on) stopAnneal(true); renderAnnealReadout(); updatePrediction(); });
  $('nb-r2').addEventListener('change', (e) => { S.anneal.r2 = +e.target.value; if (S.anneal.active && S.solver) { S.solver.setParams({ r: S.anneal.r2 }); } renderAnnealReadout(); updatePrediction(); });
  $('nb-b2').addEventListener('change', (e) => { S.anneal.b2 = +e.target.value; if (S.anneal.active && S.solver) { S.solver.setParams({ b: S.anneal.b2 }); } renderAnnealReadout(); updatePrediction(); });
  $('sel-anneal-when').addEventListener('change', (e) => { S.anneal.when = e.target.value; renderAnnealReadout(); updatePrediction(); });
  $('nb-anneal-step').addEventListener('change', (e) => { S.anneal.step = Math.max(100, +e.target.value | 0); renderAnnealReadout(); updatePrediction(); });
  bindPhase();
  for (const ch of document.querySelectorAll('#strain-chips .chip')) ch.addEventListener('click', () => { S.strainFam = +ch.dataset.fam; for (const o of document.querySelectorAll('#strain-chips .chip')) o.classList.toggle('on', o === ch); updateDirectionLayers(); });
  $('nb-strain-range').addEventListener('change', (e) => { S.strainRange = Math.max(0.02, Math.min(0.6, +e.target.value / 100)); updateDirectionLayers(); });
  for (const ch of document.querySelectorAll('#orient-chips .chip')) ch.addEventListener('click', () => { S.orientMode = ch.dataset.mode; for (const o of document.querySelectorAll('#orient-chips .chip')) o.classList.toggle('on', o === ch); updateDirectionLayers(); });
  // per-direction read-outs
  for (const ch of document.querySelectorAll('#side-chips .chip')) ch.addEventListener('click', () => setSide(+ch.dataset.side));
  // view
  for (const chip of document.querySelectorAll('#chips .chip, #dens-chips .chip')) chip.addEventListener('click', () => setLayer(chip.dataset.layer, !chip.classList.contains('on')));
  $('sl-lathw').addEventListener('input', (e) => { $('out-lathw').value = (+e.target.value).toFixed(2); S.lathW = +e.target.value; view.setStripeParams({ halfWidth: stripeHalfWidth() }); updateLaths(); updateDensityReadout(); });
  $('sel-lathcol').addEventListener('change', (e) => { S.lathCol = e.target.value; updateLaths(); });
  $('ck-overunder').addEventListener('change', (e) => { S.overUnder = e.target.checked; updateLaths(); });
  $('sel-carrier').addEventListener('change', (e) => { S.carrierView = e.target.value; updateCarrierColors(); if (S.carrierView === 'u') showStep(S.shownStep >= 0 ? S.shownStep : 0); });
  const pct = (x) => `${Math.round(100 * x)} %`;
  $('sl-carrier-alpha').addEventListener('input', (e) => { S.carrierAlpha = +e.target.value; $('out-carrier-alpha').value = pct(S.carrierAlpha); view.setFieldAlpha(S.carrierAlpha); });
  for (let k = 0; k < 3; k++) $('sl-stripe-' + k).addEventListener('input', (e) => { S.stripeAlpha[k] = +e.target.value; $('out-stripe-' + k).value = pct(S.stripeAlpha[k]); view.setStripeParams({ alpha: S.stripeAlpha }); });
  // picking: click without drag
  const gl = $('gl'); let down = null;
  gl.addEventListener('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
  gl.addEventListener('pointerup', (e) => {
    if (!down) return; const moved = Math.hypot(e.clientX - down[0], e.clientY - down[1]); down = null;
    if (moved > 5 || !S.carrier) return;
    const hit = view.pick(e.clientX, e.clientY);
    if (!hit) { clearSelection(); return; }
    if (e.shiftKey || S.armPick) {                                           // centre picking: Shift-click or the armed button
      const v = hit.kind === 'carrier' ? hit.vertex : (hit.kind === 'vertex' && S.lat ? S.lat.peaks[hit.peak] : -1);
      const idx = v >= 0 ? v : nearestCarrierVertex(hit.point);
      S.armPick = false; gl.classList.remove('arm');
      if (idx >= 0) setCentre(idx, 'picked');
      return;
    }
    if (hit.kind === 'lath' && hit.strand >= 0) selectLath(hit.strand);
    else if (hit.kind === 'vertex') selectVertex(hit.peak);
    else clearSelection();
  });
  $('btn-pick').addEventListener('click', () => { S.armPick = !S.armPick; gl.classList.toggle('arm', S.armPick); setStatus(S.armPick ? 'click the surface to place the centre' : ''); });
  // navigation: view cube buttons, turntable, projection, keys
  $('nav-home').addEventListener('click', () => view.home());
  $('nav-fit').addEventListener('click', () => view.fit());
  $('nav-left').addEventListener('click', () => view.orbitBy(Math.PI / 2, 0));        // the model rolls left: its right face comes to the front
  $('nav-right').addEventListener('click', () => view.orbitBy(-Math.PI / 2, 0));
  $('nav-up').addEventListener('click', () => view.orbitBy(0, Math.PI / 2));          // the model rolls up: the camera goes under it
  $('nav-down').addEventListener('click', () => view.orbitBy(0, -Math.PI / 2));
  $('nav-turn').addEventListener('click', () => setTurntable(!S.turntable));
  $('nav-proj').addEventListener('click', () => setProjection(S.projection === 'persp' ? 'ortho' : 'persp'));
  $('sl-turn').addEventListener('input', (e) => { S.turnSpeed = +e.target.value; $('out-turn').value = S.turnSpeed; if (S.turntable) view.setTurntable(true, S.turnSpeed); });
  window.addEventListener('keydown', (e) => {
    if (e.key === 'Escape') { S.armPick = false; gl.classList.remove('arm'); clearSelection(); return; }
    if (e.metaKey || e.ctrlKey || e.altKey || /^(INPUT|TEXTAREA|SELECT)$/.test(e.target.tagName)) return;
    const k = e.key.toLowerCase();
    if (k === 'h') view.home(); else if (k === 'f') view.fit(); else if (k === 't') setTurntable(!S.turntable); else if (k === 'p') setProjection(S.projection === 'persp' ? 'ortho' : 'persp'); else return;
    e.preventDefault();
  });
  // growth strip: scrub by pointer
  const chart = $('chart'); let chartDown = false;
  const scrub = (e) => { S.target = chartStepAt(e.clientX); syncStepsUI(); };
  chart.addEventListener('pointerdown', (e) => { chartDown = true; chart.setPointerCapture(e.pointerId); scrub(e); });
  chart.addEventListener('pointermove', (e) => { if (chartDown) scrub(e); });
  chart.addEventListener('pointerup', () => { chartDown = false; }); chart.addEventListener('pointercancel', () => { chartDown = false; });
  // collapsible sections: remember what is folded
  let folds = {}; try { folds = JSON.parse(localStorage.getItem('localsh-folds') || '{}'); } catch (e) { folds = {}; }
  for (const d of document.querySelectorAll('details.sec')) {
    if (folds[d.id] === false) d.open = false;
    d.addEventListener('toggle', () => { folds[d.id] = d.open; try { localStorage.setItem('localsh-folds', JSON.stringify(folds)); } catch (e) { /* ignore */ } });
  }
  // theme
  for (const chip of document.querySelectorAll('#theme-chips .chip')) chip.addEventListener('click', () => setTheme(chip.dataset.theme));
  // modal
  $('modal-close').addEventListener('click', () => { $('modal').hidden = true; });
  $('modal').addEventListener('click', (e) => { if (e.target === $('modal')) $('modal').hidden = true; });
  $('modal-copy').addEventListener('click', async () => {
    const ta = $('modal-text');
    try { await navigator.clipboard.writeText(ta.value); $('modal-copy').textContent = 'Copied'; setTimeout(() => { $('modal-copy').textContent = 'Copy'; }, 1500); }
    catch (e) { ta.focus(); ta.select(); }
  });
  const onTheme = () => { view.applyTheme(); drawChart(); drawPhase(); drawDtChart(); drawDensityDiagram(); drawChargeDiagram(); drawOrientLegend(); updateProbe(); if (S.carrierView !== 'u') updateCarrierColors(); if (S.lat) { view.setLattice(S.lat, S.summary); view.setSelectedVertex(S.selVertex); updateLaths(); updateDirectionLayers(); } };
  try { window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', onTheme); } catch (e) { /* ignore */ }
  new MutationObserver(onTheme).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
  window.addEventListener('resize', () => { drawChart(); drawPhase(); drawDtChart(); });
}

async function openFile(file) {
  setStatus(`reading ${file.name}…`);
  try {
    const buf = await file.arrayBuffer();
    const { V, F } = parseSTL(buf);
    let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < V.length / 3; i++) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], V[3 * i + k]); hi[k] = Math.max(hi[k], V[3 * i + k]); }
    const bbox = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
    const c = Carrier.fromMesh(file.name, V, F, bbox / 104);
    setCarrier(c, bbox / 13); $('sel-zoo').value = '__file__';
    setStatus(`${file.name}: ${c.meta.notes.join('; ')} — prepare it with export_carrier.py for a remeshed carrier`);
  } catch (e) { console.error(e); setStatus(`could not read ${file.name}: ${e.message}`); }
}

// ------------------------------------------------------------------ boot
window.addEventListener('error', (e) => setStatus('error: ' + (e.message || e)));
window.addEventListener('unhandledrejection', (e) => setStatus('error: ' + (e.reason && e.reason.message || e.reason)));

async function main() {
  bind(); initTips();
  view.setFieldAlpha(S.carrierAlpha); view.setStripeParams({ alpha: S.stripeAlpha, halfWidth: stripeHalfWidth() });
  try { const t = localStorage.getItem('localsh-theme'); if (t && t !== 'auto') setTheme(t); } catch (e) { /* ignore */ }
  $('sl-steps').max = S.maxSteps; $('nb-max').value = S.maxSteps;
  if (window.top !== window.self) $('modal-download').hidden = true;      // the artifact viewer blocks page downloads
  S.device = new URLSearchParams(location.search).has('cpu') ? null : await getDevice();
  await loadManifest();
  const hash = decodeURIComponent(location.hash.slice(1));
  const wanted = hash && S.manifest.find(m => m.name === hash || m.file.toLowerCase().includes(hash.toLowerCase())) || S.manifest.find(m => m.name === 'v4Cone_sharp') || S.manifest[0];
  if (wanted) await loadCarrierByName(wanted.name); else setStatus('no carriers found — drop an .stl');
  syncStepsUI(); drawChart(); drawPhase(); drawDtChart(); updateViewSum();
  const loop = () => { tick(); view.render(); requestAnimationFrame(loop); };
  requestAnimationFrame(loop);
  setInterval(tick, 150);                                 // rAF stops in hidden tabs; keep the solver advancing anyway
}
window.LocalSH = S; S.api = { selectLath, selectVertex, clearSelection, setTheme, setCentre, updateLaths, setTurntable, setProjection, setRefine, setSide, view, startAnneal, stopAnneal, extract, reseedNow, setLayer, updatePrediction };   // for debugging from the console
main();
