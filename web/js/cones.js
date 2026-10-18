// Cone placement for the conformal size field ĝ = e^{2φ} g (carrier.sizeField): the "1 → 2 → 3" pipeline's third step.
//   1  count the windings the stripe phases are forced to have (branchPoints in stripes.js: the curl part of the three
//      wave vectors q0 n_k cannot be followed by a smooth phase; every left-over 60° of turning becomes a fork),
//   2  split that turning into what a size field can absorb (the exact / curl part) and what topology fixes (the net
//      charge Σ(6 − n) = 6χ on a closed surface, the rim's turning on an open one),
//   3  put the turning into a few cones (5s and 7s) and let ĝ be flat in between: L φ = Σ_c Ω_c δ_c − K. The lattice is
//      then regular in ĝ (constant q0 there) and has edge a·e^{−φ} in the model — no forced forks, no 5–7 pairs that
//      only compensate curvature.
// Greedy placement after Ben-Chen, Gotsman & Bunin 2008 ("Conformal flattening by curvature prescription and metric
// scaling"): a 5 where φ is lowest (the lattice would be stretched: uncompensated K > 0), a 7 where it is highest
// (squeezed), at least dmin·a from the other cones and a from the rim. Closed surfaces first get the cones their
// topology needs (#5 − #7 = 6χ, else ĝ keeps a uniform left-over curvature), then 5–7 pairs while they reduce the spread
// of φ; open surfaces add single cones while they reduce it (the rim absorbs the rest, φ = 0 there).
import { dijkstra, MinHeap } from './sparse.js';

const PI3 = Math.PI / 3;

/** mass weights of a cone's angle spread over the geodesic disc of radius R about vertex v (hat profile 1 − d/R, interior
 *  vertices only, Σ w = 1): [[i, w], …]. The lattice cannot resolve anything below a, so a cone smeared over ~a/2 is the
 *  same cone to it, and φ keeps no log spike (which would put 0.2 a edges at the cone and mislead the greedy placement). */
export function smearWeights(c, v, R, sign = 0) {
  const G = c.graph, dist = new Map([[v, 0]]), heap = new MinHeap(64); heap.push(0, v);
  while (heap.size) {
    const [d, i] = heap.pop(); if (d > dist.get(i)) continue;
    for (let k = G.indptr[i]; k < G.indptr[i + 1]; k++) { const j = G.indices[k], nd = d + G.data[k]; if (nd < R && nd < (dist.has(j) ? dist.get(j) : Infinity)) { dist.set(j, nd); heap.push(nd, j); } }
  }
  // with a sign: the cone's angle goes where the curvature of that sign sits inside the disc (an apex keeps its cone at
  // the apex, a dome spreads it like the dome) — mass weights only where there is no such curvature
  const out = []; let tot = 0;
  if (sign) { for (const [i, d] of dist) { if (c.boundary[i]) continue; const w = Math.max(sign * c.K[i], 0) * (1 - d / R); if (w > 0) { out.push([i, w]); tot += w; } } }
  if (!(tot > 0)) { out.length = 0; for (const [i, d] of dist) { if (c.boundary[i]) continue; const w = c.m[i] * (1 - d / R); out.push([i, w]); tot += w; } }
  if (!tot) return [[v, 1]];
  for (const e of out) e[1] /= tot;
  return out;
}

/** vertices that are not on the rim and at least `margin` (along the intrinsic edges) from it (all vertices if closed) */
export function interiorVertices(c, margin) {
  const src = []; for (let i = 0; i < c.n; i++) if (c.boundary[i]) src.push(i);
  const out = [];
  if (!src.length) { for (let i = 0; i < c.n; i++) out.push(i); return { cand: out, rim: null }; }
  const rim = dijkstra(c.graph, src).dist;
  for (let i = 0; i < c.n; i++) if (!c.boundary[i] && rim[i] >= margin) out.push(i);
  return { cand: out, rim };
}

/** mass-weighted mean of φ over `cand`, its extremes and the robust spread p98 − p2 (the cones' log spikes are left out) */
export function phiStats(phi, cand, m) {
  let s = 0, w = 0; for (const i of cand) { s += m[i] * phi[i]; w += m[i]; }
  const mean = w ? s / w : 0, d = cand.map(i => phi[i] - mean).sort((x, y) => x - y), q = (t) => d[Math.min(d.length - 1, Math.floor(t * d.length))];
  return { mean, lo: d[0], hi: d[d.length - 1], p2: q(0.02), p98: q(0.98), spread: q(0.98) - q(0.02) };
}

/** quantiles of the model lattice edge a·e^{−φ} (in units of a) over `cand`, φ already normalised to mean 0 */
export function edgeQuantiles(phi, cand) {
  const s = cand.map(i => Math.exp(-phi[i])).sort((x, y) => x - y), q = (t) => s[Math.min(s.length - 1, Math.floor(t * (s.length - 1)))];
  return { min: q(0), p5: q(0.05), median: q(0.5), p95: q(0.95), max: q(1) };
}

/** net charge Σ (6 − n) of a cone list, in units of 60° */
export const netCharge = (cones) => cones.reduce((s, q) => s + (6 - q.valence), 0);

/** the rotation field's own singular faces as cones: one per vertex (first corner of the face), indices summed */
export function fieldCones(c, singular) {
  const acc = new Map();
  for (const s of singular) { const v = c.Fi[3 * s.face]; acc.set(v, (acc.get(v) || 0) + s.index); }
  const out = []; for (const [v, idx] of acc) if (idx !== 0 && !c.boundary[v]) out.push({ vertex: v, valence: 6 - idx });
  return out;
}

/** geodesic discs: vertices within `radius` of v along the intrinsic edges, with their distances ([i, d]…) */
export function discAround(c, v, radius) {
  const G = c.graph, dist = new Map([[v, 0]]), heap = new MinHeap(64); heap.push(0, v);
  while (heap.size) {
    const [d, i] = heap.pop(); if (d > dist.get(i)) continue;
    for (let k = G.indptr[i]; k < G.indptr[i + 1]; k++) { const j = G.indices[k], nd = d + G.data[k]; if (nd < radius && nd < (dist.has(j) ? dist.get(j) : Infinity)) { dist.set(j, nd); heap.push(nd, j); } }
  }
  return dist;
}

/** Dijkstra from `src` that only visits vertices it improves on `far` (and writes them): the new centre's own cell */
function relaxFar(G, src, far) {
  const heap = new MinHeap(64); far[src] = 0; heap.push(0, src);
  while (heap.size) {
    const [d, i] = heap.pop(); if (d > far[i]) continue;
    for (let k = G.indptr[i]; k < G.indptr[i + 1]; k++) { const j = G.indices[k], nd = d + G.data[k]; if (nd < far[j]) { far[j] = nd; heap.push(nd, j); } }
  }
}

/** farthest-point sample of `cand` with spacing ≥ `spacing` (first point: the candidate of largest |K|); O(n log n) overall */
function sampleCentres(c, cand, spacing) {
  if (!cand.length) return [];
  let first = cand[0]; for (const i of cand) if (Math.abs(c.K[i]) > Math.abs(c.K[first])) first = i;
  const far = new Float64Array(c.n).fill(Infinity), centres = [first]; relaxFar(c.graph, first, far);
  for (let guard = 0; guard < 20000; guard++) {
    let best = -1, bd = spacing; for (const i of cand) if (far[i] > bd) { bd = far[i]; best = i; }
    if (best < 0) break;
    centres.push(best); relaxFar(c.graph, best, far);
  }
  return centres;
}

/** nearest vertex of `targets` (a Set) to `src` along the surface: Dijkstra with early exit */
function nearestOf(G, src, targets) {
  const dist = new Map([[src, 0]]), heap = new MinHeap(64); heap.push(0, src);
  while (heap.size) {
    const [d, i] = heap.pop(); if (d > dist.get(i)) continue;
    if (i !== src && targets.has(i)) return i;
    for (let k = G.indptr[i]; k < G.indptr[i + 1]; k++) { const j = G.indices[k], nd = d + G.data[k]; if (nd < (dist.has(j) ? dist.get(j) : Infinity)) { dist.set(j, nd); heap.push(nd, j); } }
  }
  return -1;
}

/**
 * placeCones(c, a, opts) on the base carrier → {cones, phi, stats, edges, history, residualDeg, net, need, closed}
 *   mode 'greedy' (quantise the curvature to the lattice scale), 'field' (the smoothest field's own singularities),
 *   'uniform' (none).
 *   greedy: the lattice sees curvature only in packets of 60° per cell, so the residual charge r = K − (cones placed so
 *   far) is swept at growing scales — discs of radius 0.75 a, 1.5 a, 3 a, … — and wherever a disc still holds ≥ 30°
 *   (plus/minus) a cone of charge round(C / 60°) (±1 … ±3, i.e. valence 5/4/3 or 7/8/9) goes to the vertex carrying most
 *   of it. Small scales catch tips and junctions, large scales what adds up over a dome or a saddle. A closed surface
 *   ends with the net charge its topology needs (#5 − #7 = 6χ). Cones are smeared over `smear`·a for φ (see smearWeights);
 *   phi is normalised so that the interior has the same area in ĝ as in the model (as many sites as a uniform lattice of edge a).
 */
export function placeCones(c, a, { mode = 'greedy', fieldSingular = null, threshold = 0.5, maxCharge = 3, margin = 1, smear = 0.5, relax = 0, log = null } = {}) {
  const { cand, rim } = interiorVertices(c, margin * a);
  const closed = !rim, need = closed ? Math.round(c.omega / PI3) : null;        // closed: #5 − #7 must be 6χ = Ω/60°
  const m = c.m, n = c.n, history = [], R = smear * a;
  let cones = [];
  const finish = () => ({ ...sizeFieldFor(c, cones, a, { smear: R, cand }), history, need, closed });
  if (mode === 'uniform') return { cones: [], phi: null, stats: null, edges: null, history, residualDeg: null, net: 0, need, closed };
  if (mode === 'field') { cones = fieldSingular ? fieldCones(c, fieldSingular) : []; return finish(); }
  // residual curvature per vertex (interior); a placed cone takes its angle out at its vertex
  const r = Float64Array.from(c.K); for (let i = 0; i < n; i++) if (c.boundary[i]) r[i] = 0;
  const inCand = new Uint8Array(n); for (const i of cand) inCand[i] = 1;
  const place = (v, q) => { cones.push({ vertex: v, valence: 6 - q }); r[v] -= q * PI3; };   // point bookkeeping: no spill into the neighbouring cells (φ alone sees the smear)
  const discCharge = (disc) => { let s = 0; for (const [i] of disc) s += r[i]; return s; };
  // the vertex of a disc that carries most of its residual (of the disc's sign), weighted towards the centre
  const carrier = (disc, sgn) => { let best = -1, bv = -Infinity; for (const [i, d] of disc) { if (!inCand[i]) continue; const v = sgn * r[i]; if (v > bv) { bv = v; best = i; } } return best; };
  // error diffusion ("dithering") of the curvature along a tour of lattice-sized cells: the charge runs up cell by
  // cell and every time it passes ±30° a cone of round(acc / 60°) is placed in the cell that contributed most since the
  // last cone; the remainder (< 30°) is carried to the next cell. So no region at any scale ≥ a keeps more than ~30° of
  // uncompensated curvature, and ĝ is flat at the scales the lattice can see.
  const centres = sampleCentres(c, cand, 0.75 * a), { source } = dijkstra(c.graph, centres), cells = new Map();
  for (const i of cand) { const k = source[i]; if (k < 0) continue; let l = cells.get(k); if (!l) cells.set(k, l = []); l.push(i); }
  const cellList = [...cells.entries()].map(([k, l]) => { let C = 0; for (const i of l) C += r[i]; return { k, l, C }; });
  // tour: always the nearest unvisited cell (along the surface), starting at the cell of largest |charge|
  const pos = new Map(cellList.map((cl, idx) => [cl.k, idx])), unvisited = new Set(cellList.map(cl => cl.k)), tour = [];
  let cur = 0; for (let j = 1; j < cellList.length; j++) if (Math.abs(cellList[j].C) > Math.abs(cellList[cur].C)) cur = j;
  for (let step = 0; step < cellList.length; step++) {
    unvisited.delete(cellList[cur].k); tour.push(cur);
    if (!unvisited.size) break;
    const nk = nearestOf(c.graph, cellList[cur].k, unvisited); if (nk < 0) break;   // disconnected rest: the tour ends there
    cur = pos.get(nk);
  }
  if (tour.length < cellList.length) for (const [k, idx] of pos) if (unvisited.has(k)) tour.push(idx);
  const carrierOf = (cl, sgn) => { let v = -1, bv = -Infinity; for (const i of cl.l) { const val = sgn * r[i]; if (val > bv) { bv = val; v = i; } } return v; };
  let acc = 0, bestCell = null, bestC = 0;
  for (const idx of tour) {
    const cl = cellList[idx]; let C = 0; for (const i of cl.l) C += r[i];
    acc += C; if (Math.abs(C) > Math.abs(bestC)) { bestC = C; bestCell = cl; }
    while (Math.abs(acc) >= threshold * PI3) {
      const sgn = Math.sign(acc), q = sgn * Math.max(1, Math.min(maxCharge, Math.round(Math.abs(acc) / PI3)));
      const target = bestCell && Math.sign(bestC) === sgn ? bestCell : cl, v = carrierOf(target, sgn);
      if (v < 0) break;
      place(v, q); acc -= q * PI3; bestCell = null; bestC = 0;
      if (log) log(`cell ${target.k}: ${6 - q} at ${v} (running charge was ${((acc + q * PI3) * 180 / Math.PI).toFixed(0)}°)`);
    }
  }
  history.push(['tour', cellList.length, cones.length, +(acc * 180 / Math.PI).toFixed(0)]);
  if (relax) {
    // Lloyd relaxation: cones of each sign share out the cells of that sign (nearest cone along the surface); each cone
    // moves to the cell of its share nearest to the |charge|-weighted centroid. Evens out the spacing the tour left.
    const cellC = cellList.map(cl => { let C = 0; for (const i of cl.l) C += c.K[i]; return C; });   // the raw charge of each cell
    const cx = cellList.map(cl => { let x = 0, y = 0, z = 0, w = 0; for (const i of cl.l) { x += m[i] * c.V[3 * i]; y += m[i] * c.V[3 * i + 1]; z += m[i] * c.V[3 * i + 2]; w += m[i]; } return [x / w, y / w, z / w]; });
    const cellOfVertex = new Int32Array(n).fill(-1); cellList.forEach((cl, idx) => { for (const i of cl.l) cellOfVertex[i] = idx; });
    for (let it = 0; it < relax; it++) {
      for (const sgn of [1, -1]) {
        const mine = cones.filter(q => Math.sign(6 - q.valence) === sgn); if (!mine.length) continue;
        const { source } = dijkstra(c.graph, mine.map(q => q.vertex)), share = new Map(mine.map(q => [q.vertex, []]));
        cellList.forEach((cl, idx) => { if (Math.sign(cellC[idx]) !== sgn) return; const owner = source[cl.k]; if (owner >= 0 && share.has(owner)) share.get(owner).push(idx); });
        for (const q of mine) {
          const idxs = share.get(q.vertex); if (!idxs || idxs.length < 2) continue;
          let x = 0, y = 0, z = 0, w = 0; for (const idx of idxs) { const wt = Math.abs(cellC[idx]); x += wt * cx[idx][0]; y += wt * cx[idx][1]; z += wt * cx[idx][2]; w += wt; }
          if (!(w > 0)) continue; x /= w; y /= w; z /= w;
          let best = -1, bd = Infinity; for (const idx of idxs) { const d = (cx[idx][0] - x) ** 2 + (cx[idx][1] - y) ** 2 + (cx[idx][2] - z) ** 2; if (d < bd) { bd = d; best = idx; } }
          if (best >= 0) { const v = carrierOf(cellList[best], sgn); if (v >= 0) q.vertex = v; }
        }
      }
    }
    // the residual bookkeeping follows the moved cones
    r.set(c.K); for (let i = 0; i < n; i++) if (c.boundary[i]) r[i] = 0; for (const q of cones) r[q.vertex] -= (6 - q.valence) * PI3;
    history.push(['relax', relax]);
    void cellOfVertex;
  }
  // a closed surface must end balanced: the needed charge goes where the residual of that sign is largest (cell scale)
  for (let guard = 0; closed && netCharge(cones) !== need && guard < 60; guard++) {
    const missing = need - netCharge(cones), sgn = Math.sign(missing), centres = sampleCentres(c, cand, 0.375 * a);
    let best = -1, bc = -Infinity; for (const v of centres) { const C = sgn * discCharge(discAround(c, v, 0.75 * a)); if (C > bc) { bc = C; best = v; } }
    if (best < 0) break;
    place(best, sgn * Math.min(Math.abs(missing), maxCharge)); history.push(['balance', netCharge(cones), need]);
    if (log) log(`balance ${netCharge(cones)}/${need}: ${6 - sgn} at ${best} (disc held ${(bc * 180 / Math.PI).toFixed(0)}°)`);
  }
  return finish();
}

/** the size field of a cone list (base carrier), normalised so that the ĝ area of the interior equals its model area
 *  (Σ m e^{2φ} = Σ m: as many sites as a uniform lattice of edge a, the edge then varies around a) → {cones, phi, stats,
 *  edges, residualDeg, net}. Cones may carry their own smear radius (`smear`, model units) and are smeared over `smear`·… */
export function sizeFieldFor(c, cones, a, { smear = 0.5 * a, cand = null } = {}) {
  const m = c.m; if (!cand) cand = interiorVertices(c, a).cand;
  const res = c.sizeField(cones, a, -1, { tol: 1e-8, smear });
  const phi = res.phi; let s2 = 0, w = 0; for (const i of cand) { s2 += m[i] * Math.exp(2 * phi[i]); w += m[i]; }
  const shift = 0.5 * Math.log(s2 / (w || 1)); for (let i = 0; i < phi.length; i++) phi[i] -= shift;
  const st = phiStats(phi, cand, m); let v = 0; for (const i of cand) v += m[i] * (phi[i] - st.mean) ** 2;
  return { cones, phi, stats: { ...st, std: Math.sqrt(v / (w || 1)) }, edges: edgeQuantiles(phi, cand), residualDeg: res.residualDeg, net: netCharge(cones) };
}

/** user edits applied to an automatic cone list: edits = [{vertex, on, valence, radius (a), angle (deg), weight}] keyed by the
 *  cone's vertex. Returns the cones for ĝ (off ones dropped, valence and smear overridden) and the annotated full list. */
export function applyConeEdits(cones, edits, a) {
  const byV = new Map((edits || []).map(e => [e.vertex, e]));
  const all = cones.map(q => { const e = byV.get(q.vertex) || {}; return { vertex: q.vertex, auto: q.valence, valence: e.valence != null ? e.valence : q.valence, on: e.on !== false, radius: e.radius != null ? e.radius : null, angle: e.angle || 0, weight: e.weight || 0 }; });
  const active = all.filter(q => q.on).map(q => ({ vertex: q.vertex, valence: q.valence, ...(q.radius != null ? { smear: q.radius * a } : {}) }));
  return { all, active };
}

/** ĝ edge lengths ℓ̂ = e^{(φi+φj)/2} ℓ per face corner, for any intrinsic triangulation (Fi, Li). A face that the
 *  per-vertex scaling would flatten (a skinny model triangle across a steep φ: the triangle inequality fails or nearly,
 *  which puts cotan weights of 1e5 into the Laplacian and blows the solver up) is scaled uniformly by e^{φ̄_f} instead,
 *  which keeps its shape. Returns the lengths; `.relaxed` counts those faces. */
export function conformalLengths(Fi, Li, phi) {
  const out = new Float64Array(Li.length), nf = Fi.length / 3; let relaxed = 0;
  for (let f = 0; f < nf; f++) {
    const i = Fi[3 * f], j = Fi[3 * f + 1], k = Fi[3 * f + 2];
    let x = Li[3 * f] * Math.exp(0.5 * (phi[i] + phi[j])), y = Li[3 * f + 1] * Math.exp(0.5 * (phi[j] + phi[k])), z = Li[3 * f + 2] * Math.exp(0.5 * (phi[k] + phi[i]));
    if (x >= 0.97 * (y + z) || y >= 0.97 * (z + x) || z >= 0.97 * (x + y)) {
      const e = Math.exp((phi[i] + phi[j] + phi[k]) / 3); x = Li[3 * f] * e; y = Li[3 * f + 1] * e; z = Li[3 * f + 2] * e; relaxed++;
    }
    out[3 * f] = x; out[3 * f + 1] = y; out[3 * f + 2] = z;
  }
  out.relaxed = relaxed;
  return out;
}

/** base-level per-vertex values → refined carrier (midpoints get the mean of their edge ends), parents from Carrier.refined */
export function interpByParents(vals, parents) {
  let cur = vals;
  for (const par of parents) {
    const n0 = cur.length, out = new Float64Array(n0 + par.length / 2); out.set(cur);
    for (let i = 0; i < par.length / 2; i++) out[n0 + i] = 0.5 * (cur[par[2 * i]] + cur[par[2 * i + 1]]);
    cur = out;
  }
  return cur;
}
