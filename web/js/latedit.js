// Connectivity edits of a frozen lattice (smooth.js · freeze). The field does not know them: they live on the lattice
// until Unfreeze.
//   flipEdge · collapseEdge · splitEdge   one edge, named by its two sites (collapse and split also on a naked edge of the rim)
//   healCluster · healAll                 the 5s and 7s round a site taken down to their net charge by a short sequence
//                                         of those three, found by search (5-7-5 → 5, 7-5-7 → 7, 5-7 + 7-5 → nothing);
//                                         at the rim they are taken out altogether: rim sites have no valence to keep
//   toRim                                 a 5–7 pair walked to the nearest rim (flips glide it, collapses and split + flip
//                                         change its row) and taken out there
//   patchSite · patchRim                  the fan of a rim site completed to six triangles; the rim grown until the
//                                         lattice covers the carrier, its new sites regular, on the carrier or — past the
//                                         boundary — on its continuation (F.ext: the quadric each such site lies on)
//   snapshot · restore                    undo
// Valences: a flip takes one from each end of the edge and gives one to each opposite site; a collapse leaves
// val(u) + val(v) − 4 on the merged site and takes one from each opposite site; a split makes a 4 and gives one to each
// opposite site. On a naked edge there is one opposite site: collapsing the edge takes one from it (a 7 behind the rim
// becomes a 6), splitting it gives it one (a 5 becomes a 6); the rim gains or loses a site. A new or merged site starts at the middle of its edge (a patched one where the regular lattice would
// continue), on the carrier, and its edges ask for the mean edge length around (F.restFix): relax (smooth.js) spreads
// the change. After an edit the lattice is its positions as frozen with the new connectivity, to be relaxed again.
import { project, rhombi, rhombusClass, fitExt, liftExt, relax } from './smooth.js';

const KEY = 2097152;                                   // an edge i < j is the number i · 2^21 + j (smooth.js · tables reads F.restFix by it)
export const edgeKey = (i, j) => (i < j ? i * KEY + j : j * KEY + i);
const sq = (x) => x * x;
const dist = (X, i, j) => Math.hypot(X[3 * i] - X[3 * j], X[3 * i + 1] - X[3 * j + 1], X[3 * i + 2] - X[3 * j + 2]);

/** valence (triangles at the site), interior (closed link), rhombi and their classes from the triangles; drops what was cached */
export function derive(F) {
  const T = F.tris, n = F.peaks.length, valence = new Int32Array(n), interior = new Uint8Array(n), cnt = new Map(), bad = new Uint8Array(n);
  for (let k = 0; k < T.length; k++) valence[T[k]]++;
  for (let t = 0; t < T.length; t += 3) for (let k = 0; k < 3; k++) { const key = edgeKey(T[t + k], T[t + (k + 1) % 3]); cnt.set(key, (cnt.get(key) || 0) + 1); }
  for (const [key, m] of cnt) if (m !== 2) { bad[Math.floor(key / KEY)] = 1; bad[key % KEY] = 1; }
  for (let i = 0; i < n; i++) interior[i] = !bad[i] && valence[i] >= 3 ? 1 : 0;
  F.valence = valence; F.interior = interior; F.Q = rhombi(F); F.cls = rhombusClass(F, F.Q); F._tb = null; F._nbr = null;
  return F;
}
/** site → its neighbours; kept on the lattice until the next edit */
export function neighbours(F) {
  if (F._nbr) return F._nbr;
  const T = F.tris, nb = Array.from({ length: F.peaks.length }, () => []);
  for (let t = 0; t < T.length; t += 3) for (let k = 0; k < 3; k++) { const i = T[t + k], j = T[t + (k + 1) % 3]; if (!nb[i].includes(j)) nb[i].push(j); if (!nb[j].includes(i)) nb[j].push(i); }
  return (F._nbr = nb);
}
/** every edge of a flat triangle list, key → {u, v, p, tp, q, tq, n}: (u, v, p) is triangle tp and (v, u, q) triangle tq
 *  (q = −1 on a rim edge), n the triangles on it; bad: more than two, or two that run the same way */
function edges(T) {
  const m = new Map();
  for (let t = 0; t < T.length; t += 3) {
    if (T[t] < 0) continue;
    for (let k = 0; k < 3; k++) {
      const a = T[t + k], b = T[t + (k + 1) % 3], o = T[t + (k + 2) % 3], key = edgeKey(a, b), e = m.get(key);
      if (!e) m.set(key, { u: a, v: b, p: o, tp: t, q: -1, tq: -1, n: 1, bad: false });
      else { e.n++; if (e.n === 2 && e.u === b) { e.q = o; e.tq = t; } else e.bad = true; }
    }
  }
  return m;
}

// ------------------------------------------------------------------ the lattice open for edits
function open(F) {
  const n = F.peaks.length, val = new Array(n).fill(0); for (let k = 0; k < F.tris.length; k++) val[F.tris[k]]++;
  return { F, n0: n, n, T: Array.from(F.tris), P: Array.from(F.P), N: Array.from(F.nrm), A: Array.from(F.anchor), val, touched: new Set(), dead: new Set(), mark: new Set(), fix: new Map(F.restFix || []), ext: new Map(F.ext || []), gen: new Map(), _em: null, _adj: null };
}
const em = (M) => M._em || (M._adj = null, M._em = edges(M.T));
function nbrs(M, i) {
  if (!M._adj || !M._em) { const A = new Map(), add = (x, y) => { let s = A.get(x); if (!s) A.set(x, s = new Set()); s.add(y); }; for (const e of em(M).values()) { add(e.u, e.v); add(e.v, e.u); } M._adj = A; }
  return M._adj.get(i) || new Set();
}
const onRim = (M, i) => { for (const x of nbrs(M, i)) if (em(M).get(edgeKey(i, x)).n === 1) return true; return false; };
/** the length an edge asks for: the mean set by an edit, else its length as frozen */
function restOf(M, i, j) {
  const f = M.fix.get(edgeKey(i, j)); if (f !== undefined) return f;
  return i < M.n0 && j < M.n0 && !M.touched.has(i) && !M.touched.has(j) ? dist(M.F.P0, i, j) : dist(M.P, i, j);
}
/** shape of the triangle i j k, 1 when equilateral, negative when it faces away from the surface normal at its sites */
function triQ(M, i, j, k) {
  const P = M.P, N = M.N, ux = P[3 * j] - P[3 * i], uy = P[3 * j + 1] - P[3 * i + 1], uz = P[3 * j + 2] - P[3 * i + 2], vx = P[3 * k] - P[3 * i], vy = P[3 * k + 1] - P[3 * i + 1], vz = P[3 * k + 2] - P[3 * i + 2];
  const cx = uy * vz - uz * vy, cy = uz * vx - ux * vz, cz = ux * vy - uy * vx, wx = vx - ux, wy = vy - uy, wz = vz - uz;
  const l2 = ux * ux + uy * uy + uz * uz + vx * vx + vy * vy + vz * vz + wx * wx + wy * wy + wz * wz, side = cx * (N[3 * i] + N[3 * j] + N[3 * k]) + cy * (N[3 * i + 1] + N[3 * j + 1] + N[3 * k + 1]) + cz * (N[3 * i + 2] + N[3 * j + 2] + N[3 * k + 2]);
  return (side < 0 ? -1 : 1) * 2 * Math.sqrt(3) * Math.hypot(cx, cy, cz) / (l2 || 1);
}

function mFlip(M, a, b, geo) {
  const E = em(M), e = E.get(edgeKey(a, b)); if (!e) return 'no lattice edge there';
  if (e.n !== 2 || e.bad) return 'an edge of the rim cannot be flipped';
  const { u, v, p, q, tp, tq } = e, T = M.T;
  if (p === q || E.has(edgeKey(p, q))) return 'the other diagonal is an edge already';
  if (M.val[u] - 1 < (onRim(M, u) ? 2 : 3) || M.val[v] - 1 < (onRim(M, v) ? 2 : 3)) return 'a site would be left with too few triangles';   // inside three, on the rim two
  if (geo && !(triQ(M, u, q, p) > 0.05 && triQ(M, q, v, p) > 0.05)) return 'the two triangles do not form a convex quad: the flipped edge would fold';
  const L = (restOf(M, u, p) + restOf(M, p, v) + restOf(M, v, q) + restOf(M, q, u)) / 4;
  T[tp] = u; T[tp + 1] = q; T[tp + 2] = p; T[tq] = q; T[tq + 1] = v; T[tq + 2] = p;
  M.val[u]--; M.val[v]--; M.val[p]++; M.val[q]++;
  M.fix.delete(edgeKey(u, v)); M.fix.set(edgeKey(p, q), L); M._em = null; for (const x of [u, v, p, q]) M.mark.add(x);
  return { sites: [u, v, p, q] };
}
function mCollapse(M, a, b) {
  const E = em(M), e = E.get(edgeKey(a, b)); if (!e) return 'no lattice edge there';
  if (e.bad || e.n > 2) return 'this edge has more than two triangles';
  const { u, v, p, q, tp, tq } = e, T = M.T, naked = e.n === 1, Nu = new Set(nbrs(M, u)), Nv = new Set(nbrs(M, v));
  let common = 0; for (const x of Nu) if (Nv.has(x)) common++;
  if (common !== (naked ? 1 : 2)) return 'collapsing this edge would pinch the lattice';
  const ru = onRim(M, u), rv = onRim(M, v); if (!naked && ru && rv) return 'both ends are on the rim and the edge is not: collapsing it would pinch the lattice';
  const merged = M.val[u] + M.val[v] - (naked ? 2 : 4), few = (x) => M.val[x] - 1 < (onRim(M, x) ? 1 : 3);
  if (few(p) || (!naked && few(q)) || merged < (ru || rv ? 1 : 3)) return 'a site would be left with two triangles';
  // the site that stays: the one on the rim when only one is, else the lower number (search() keeps the same)
  const keep = ru !== rv ? (rv ? v : u) : Math.min(u, v), gone = keep === u ? v : u;
  let L = 0, m = 0; for (const x of Nu) if (x !== v) { L += restOf(M, u, x); m++; } for (const x of Nv) if (x !== u) { L += restOf(M, v, x); m++; } L /= m;
  if (ru === rv) for (let k = 0; k < 3; k++) M.P[3 * keep + k] = 0.5 * (M.P[3 * u + k] + M.P[3 * v + k]);   // both inside, or both on the rim: the middle
  T[tp] = T[tp + 1] = T[tp + 2] = -1; if (!naked) T[tq] = T[tq + 1] = T[tq + 2] = -1;
  for (let k = 0; k < T.length; k++) if (T[k] === gone) T[k] = keep;
  M.val[keep] = merged; M.val[gone] = 0; M.val[p]--; if (!naked) M.val[q]--;
  for (const x of Nu) M.fix.delete(edgeKey(u, x)); for (const x of Nv) M.fix.delete(edgeKey(v, x));
  for (const x of new Set([...Nu, ...Nv])) if (x !== u && x !== v) M.fix.set(edgeKey(keep, x), L);
  if (M.ext.has(gone) && !M.ext.has(keep)) M.ext.set(keep, M.ext.get(gone));
  M.ext.delete(gone); M.dead.add(gone); M.touched.add(keep); M._em = null; for (const x of naked ? [keep, p] : [keep, p, q]) M.mark.add(x);
  return { sites: naked ? [keep, p] : [keep, p, q], gone };
}
function mSplit(M, a, b) {
  const E = em(M), e = E.get(edgeKey(a, b)); if (!e) return 'no lattice edge there';
  if (e.bad || e.n > 2) return 'this edge has more than two triangles';
  const { u, v, p, q, tp, tq } = e, T = M.T, P = M.P, N = M.N, m = M.n++;
  let L = restOf(M, u, v) + restOf(M, u, p) + restOf(M, p, v), cnt = 3; if (q >= 0) { L += restOf(M, v, q) + restOf(M, q, u); cnt = 5; } L /= cnt;
  const nl = Math.hypot(N[3 * u] + N[3 * v], N[3 * u + 1] + N[3 * v + 1], N[3 * u + 2] + N[3 * v + 2]) || 1;
  for (let k = 0; k < 3; k++) { P.push(0.5 * (P[3 * u + k] + P[3 * v + k])); N.push((N[3 * u + k] + N[3 * v + k]) / nl); }
  M.A.push(M.A[u]); if (M.ext.has(u) && M.ext.has(v)) M.ext.set(m, M.ext.get(u));
  T[tp] = u; T[tp + 1] = m; T[tp + 2] = p; T.push(m, v, p); M.val[p]++; M.val[m] = 2;
  if (q >= 0) { T[tq] = v; T[tq + 1] = m; T[tq + 2] = q; T.push(m, u, q); M.val[q]++; M.val[m] = 4; }
  M.fix.delete(edgeKey(u, v)); for (const x of [u, v, p, q]) if (x >= 0) { M.fix.set(edgeKey(m, x), L); M.mark.add(x); }
  M.mark.add(m); M._em = null;
  return { sites: [u, v, p, q].filter(x => x >= 0), born: m };
}
/** The edits written back: dead sites leave (the others keep their order, new ones follow), new and merged sites go
 *  onto the carrier, and every site returns to its position as frozen (sweeps 0), to be relaxed with the new
 *  connectivity. Returns old index → new index (−1: gone). */
function close(c, M) {
  const F = M.F, n = M.n, remap = new Int32Array(n).fill(-1); let k = 0;
  for (let i = 0; i < n; i++) if (!M.dead.has(i)) remap[i] = k++;
  const peaks = new Int32Array(k), P0 = new Float64Array(3 * k), nrm0 = new Float64Array(3 * k), anchor0 = new Int32Array(k), shift = new Float64Array(k), out = new Float64Array(7);
  for (let i = 0; i < n; i++) {
    const j = remap[i]; if (j < 0) continue;
    if (i < M.n0 && !M.touched.has(i)) { for (let d = 0; d < 3; d++) { P0[3 * j + d] = F.P0[3 * i + d]; nrm0[3 * j + d] = F.nrm0[3 * i + d]; } anchor0[j] = F.anchor0[i]; peaks[j] = F.peaks[i]; shift[j] = F.shift ? F.shift[i] : 0; continue; }
    const E = M.ext.get(i);
    if (E) { liftExt(E, M.P[3 * i], M.P[3 * i + 1], M.P[3 * i + 2], out); out[6] = M.A[i]; } else project(c, M.A[i], M.P[3 * i], M.P[3 * i + 1], M.P[3 * i + 2], out);
    for (let d = 0; d < 3; d++) { P0[3 * j + d] = out[d]; nrm0[3 * j + d] = out[3 + d]; }
    anchor0[j] = out[6]; peaks[j] = i < M.n0 ? F.peaks[i] : out[6]; shift[j] = i < M.n0 && F.shift ? F.shift[i] : 0;
  }
  const tr = []; for (let t = 0; t < M.T.length; t += 3) if (M.T[t] >= 0) tr.push(remap[M.T[t]], remap[M.T[t + 1]], remap[M.T[t + 2]]);
  const fix = new Map(); for (const [key, L] of M.fix) { const i = remap[Math.floor(key / KEY)], j = remap[key % KEY]; if (i >= 0 && j >= 0) fix.set(edgeKey(i, j), L); }
  const ext = new Map(); for (const [i, E] of M.ext) if (remap[i] >= 0) ext.set(remap[i], E);
  Object.assign(F, { peaks, P0, nrm0, anchor0, shift, tris: Int32Array.from(tr), restFix: fix, ext, P: Float64Array.from(P0), nrm: Float64Array.from(nrm0), anchor: Int32Array.from(anchor0), sweeps: 0 });
  derive(F); return remap;
}
function single(c, F, fn, a, b) {
  const M = open(F), r = fn(M, a, b, true); if (typeof r === 'string') return { ok: false, why: r };
  const before = r.sites.map(i => F.valence[i]), after = r.sites.map(i => M.val[i]);
  if (r.gone !== undefined) before[0] = `${F.valence[a]}+${F.valence[b]}`;
  if (r.born !== undefined) { before.push('·'); after.push(M.val[r.born]); }
  close(c, M); F.edits = (F.edits || 0) + 1;
  return { ok: true, before, after };
}
/** one edge, a–b. Each returns {ok, why} or {ok, before, after}: the valences of the sites it changed
 *  (flip, split: the two ends, then the two opposite sites, then the new site; collapse: the merged site, then the two opposite) */
export const flipEdge = (c, F, a, b) => single(c, F, mFlip, a, b);
export const collapseEdge = (c, F, a, b) => single(c, F, mCollapse, a, b);
export const splitEdge = (c, F, a, b) => single(c, F, mSplit, a, b);

// ------------------------------------------------------------------ a cluster of 5s and 7s down to its net charge
/** the singular sites that hang together with `seed`: each within `reach` edges of another (2: neighbours, or one site between) */
export function cluster(F, seed, reach = 2) {
  const nb = neighbours(F), sing = (i) => F.interior[i] && F.valence[i] !== 6;
  if (seed < 0 || !sing(seed)) return [];
  const seen = new Set([seed]), stack = [seed];
  while (stack.length) {
    const i = stack.pop(), near = new Set([i]); let front = [i];
    for (let d = 0; d < reach; d++) { const nf = []; for (const x of front) for (const y of nb[x]) if (!near.has(y)) { near.add(y); nf.push(y); } front = nf; }
    for (const j of near) if (sing(j) && !seen.has(j)) { seen.add(j); stack.push(j); }
  }
  return [...seen].sort((x, y) => x - y);
}
const energy = (F, S) => { let E = 0, q = 0; for (const i of S) { E += sq(F.valence[i] - 6); q += 6 - F.valence[i]; } return { E, q }; };
/** The shortest sequences of flips, collapses and splits (at most maxDepth, each one bettering a 5 or a 7) that leave
 *  the least Σ (valence − 6)² on the cluster and the `rings` rings round it. Only valences are followed here, on the
 *  triangles at those sites; new sites take the numbers n, n + 1, …; a collapse keeps the site on the rim when only one
 *  is, else the lower number.
 *  rim (true): rim sites count for nothing in the sum (they have no valence to keep), and a naked edge can be split
 *  (the site behind it +1) or collapsed (−1): a 5 or 7 by the rim is taken out there. Without it they are left alone.
 *  neutral: no step may raise the sum, which leaves few steps to try at each turn — a 5–7 pair is walked (flip by
 *  flip, the sum unchanged) to the 5 or 7 that takes it in, in sequences too long for the plain search.
 *  walk (the distance of every site to the rim, rimDistance): also returned, as `walk`, the sequences that leave the
 *  sum as it is and lower Σ distance over the 5s and 7s the most — one step of a pair towards the rim — and as
 *  `glide` those among them made of flips alone (no site added or removed).
 *  hold: sites no step may touch; forbid: an edge the first step may not flip. */
function search(F, sites, { maxDepth = 3, budget = 60000, rings = 1, neutral = false, rim = true, walk = null, hold = null, forbid = -1 } = {}) {
  const n = F.peaks.length, T = F.tris, nb = neighbours(F), core = new Set(sites);
  { let front = [...sites]; for (let r = 0; r < rings; r++) { const nf = []; for (const s of front) for (const x of nb[s]) if (!core.has(x)) { core.add(x); nf.push(x); } front = nf; } }
  const loc = []; for (let t = 0; t < T.length; t += 3) if (core.has(T[t]) || core.has(T[t + 1]) || core.has(T[t + 2])) loc.push(T[t], T[t + 1], T[t + 2]);
  // in: the sites whose triangles are all here (the core); wt0: 1 for an interior site, 0 for one of the rim
  const inC = new Uint8Array(n + maxDepth), wt0 = new Uint8Array(n), val0 = new Int8Array(n + maxDepth); let E0 = 0, q0 = 0, P0 = 0, atRim = false;
  for (const i of core) {
    val0[i] = F.valence[i];
    if (F.interior[i]) { inC[i] = 1; wt0[i] = 1; E0 += sq(F.valence[i] - 6); q0 += 6 - F.valence[i]; if (walk && F.valence[i] !== 6) P0 += walk[i]; }
    else if (rim && F.valence[i] > 0) { inC[i] = 1; atRim = true; }
  }
  for (let k = 0; k < maxDepth; k++) inC[n + k] = 1;
  if (hold) for (const i of hold) if (inC[i]) { inC[i] = 0; if (wt0[i]) { E0 -= sq(F.valence[i] - 6); q0 -= 6 - F.valence[i]; if (walk && F.valence[i] !== 6) P0 -= walk[i]; } }
  const target = atRim ? 0 : Math.abs(q0), best = { E: E0, depth: 0, seqs: [] }, bestW = { phi: P0, depth: 0, seqs: [] }, bestG = { phi: P0, depth: 0, seqs: [] }, seq = [];
  const note = (B, P2, D) => { if (P2 < B.phi - 1e-9) { B.phi = P2; B.depth = D; B.seqs = [seq.slice()]; } else if (B.depth === D && Math.abs(P2 - B.phi) < 1e-9 && B.seqs.length < 20) B.seqs.push(seq.slice()); };
  const d = (v, s) => sq(v + s - 6) - sq(v - 6), sg = (v) => (v !== 6 ? 1 : 0); let nodes = 0;
  // a node: its triangles and valences, the new sites so far (born; mask: which of them are on the rim; dn: their distance to the rim), the sum E and the distance sum Phi
  const dfs = (Tl, val, born, mask, dn, E, Phi, depth, D, last) => {
    const E_ = edges(Tl); let adj = null;
    const W = (i) => (i < n ? wt0[i] : (mask >> (i - n)) & 1 ? 0 : 1), DR = (i) => (walk ? (i < n ? walk[i] : dn[i - n]) : 0), dP = (i, v, s) => (walk && W(i) ? DR(i) * (sg(v + s) - sg(v)) : 0);
    const common = (u, v) => {
      if (!adj) { adj = new Map(); const add = (x, y) => { let s = adj.get(x); if (!s) adj.set(x, s = new Set()); s.add(y); }; for (const e of E_.values()) { add(e.u, e.v); add(e.v, e.u); } }
      let cm = 0; const B = adj.get(v); for (const x of adj.get(u)) if (B.has(x)) cm++; return cm;
    };
    const without = (gone, keep, t1, t2) => { const T2 = []; for (let t = 0; t < Tl.length; t += 3) { if (t === t1 || t === t2) continue; for (let k = 0; k < 3; k++) T2.push(Tl[t + k] === gone ? keep : Tl[t + k]); } return T2; };
    const visit = (op, a, b, E2, P2, make) => {
      if ((neutral && E2 > E) || ++nodes > budget) return;
      seq.push({ op, a, b });
      if (depth + 1 === D) {
        if (E2 < best.E) { best.E = E2; best.depth = D; best.seqs = [seq.slice()]; } else if (E2 === best.E && E2 < E0 && best.depth === D && best.seqs.length < 40) best.seqs.push(seq.slice());
        if (walk && E2 === E0) { note(bestW, P2, D); if (seq.every(o => o.op === 'flip')) note(bestG, P2, D); }
      } else if (E2 <= E0 + 5) { const ch = make(); dfs(ch.T, ch.val, ch.born, ch.mask, ch.dn, E2, P2, depth + 1, D, ch.made === undefined ? -1 : ch.made); }
      seq.pop();
    };
    for (const e of E_.values()) {
      if (nodes > budget) return;
      if (e.bad) continue;
      const { u, v, p, tp } = e; if (!inC[u] || !inC[v] || !inC[p]) continue;
      const vu = val[u], vv = val[v], vp = val[p], wu = W(u), wv = W(v), wp = W(p);
      if (e.n === 1) {                                     // a naked edge: both ends in the core, so it is an edge of the rim itself
        if (wu || wv || !wp) continue;
        if (born < maxDepth && vp < 6) {
          const id = n + born;
          visit('split', u, v, E + d(vp, 1), Phi + dP(p, vp, 1), () => { const T2 = Tl.slice(), v2 = val.slice(); T2[tp] = u; T2[tp + 1] = id; T2[tp + 2] = p; T2.push(id, v, p); v2[id] = 2; v2[p]++; return { T: T2, val: v2, born: born + 1, mask: mask | (1 << born), dn }; });
        }
        if (vp > 6 && common(u, v) === 1) {
          const keep = Math.min(u, v), gone = Math.max(u, v);
          visit('collapse', u, v, E + d(vp, -1), Phi + dP(p, vp, -1), () => { const v2 = val.slice(); v2[keep] = vu + vv - 2; v2[gone] = 0; v2[p]--; return { T: without(gone, keep, tp, -1), val: v2, born, mask, dn }; });
        }
        continue;
      }
      if (e.n !== 2) continue;
      const w = e.q, tq = e.tq; if (!inC[w]) continue;
      const vw = val[w], ww = W(w);
      // flip: an interior site keeps four triangles, a rim site two; not the edge the step before made (that would undo it)
      if (((wu && vu > 6) || (wv && vv > 6) || (wp && vp < 6) || (ww && vw < 6)) && vu > (wu ? 4 : 2) && vv > (wv ? 4 : 2) && edgeKey(u, v) !== last && !E_.has(edgeKey(p, w)))
        visit('flip', u, v, E + wu * d(vu, -1) + wv * d(vv, -1) + wp * d(vp, 1) + ww * d(vw, 1), Phi + dP(u, vu, -1) + dP(v, vv, -1) + dP(p, vp, 1) + dP(w, vw, 1),
          () => { const T2 = Tl.slice(), v2 = val.slice(); T2[tp] = u; T2[tp + 1] = w; T2[tp + 2] = p; T2[tq] = w; T2[tq + 1] = v; T2[tq + 2] = p; v2[u]--; v2[v]--; v2[p]++; v2[w]++; return { T: T2, val: v2, born, mask, dn, made: edgeKey(p, w) }; });
      // collapse (not a chord between two rim sites)
      if ((wu || wv) && vp > (wp ? 4 : 2) && vw > (ww ? 4 : 2)) {
        const m = vu + vv - 4, keep = wu !== wv ? (wu ? v : u) : Math.min(u, v), gone = keep === u ? v : u, vg = val[gone], side = wp * d(vp, -1) + ww * d(vw, -1), sideP = dP(p, vp, -1) + dP(w, vw, -1); let E2 = NaN, P2 = Phi;
        if (wu && wv) { if (m >= 4 && m <= 8 && ((wp && vp > 6) || (ww && vw > 6) || sq(m - 6) < sq(vu - 6) + sq(vv - 6))) { E2 = E - sq(vu - 6) - sq(vv - 6) + sq(m - 6) + side; if (walk) P2 = Phi - DR(gone) * sg(vg) + DR(keep) * (sg(m) - sg(val[keep])) + sideP; } }
        else if (m >= 2 && (vg !== 6 || (wp && vp > 6) || (ww && vw > 6))) { E2 = E - sq(vg - 6) + side; if (walk) P2 = Phi - DR(gone) * sg(vg) + sideP; }   // an interior site merged into the rim
        if (E2 === E2 && common(u, v) === 2) visit('collapse', u, v, E2, P2, () => { const v2 = val.slice(); v2[keep] = m; v2[gone] = 0; v2[p]--; v2[w]--; return { T: without(gone, keep, tp, tq), val: v2, born, mask, dn }; });
      }
      // split: a new interior site of valence 4
      if (born < maxDepth && ((wp && vp < 6) || (ww && vw < 6))) {
        const id = n + born, dnew = Math.min(DR(u), DR(v));
        visit('split', u, v, E + 4 + wp * d(vp, 1) + ww * d(vw, 1), Phi + dnew + dP(p, vp, 1) + dP(w, vw, 1),
          () => { const T2 = Tl.slice(), v2 = val.slice(); T2[tp] = u; T2[tp + 1] = id; T2[tp + 2] = p; T2[tq] = v; T2[tq + 1] = id; T2[tq + 2] = w; T2.push(id, v, p, id, u, w); v2[id] = 4; v2[p]++; v2[w]++; let d2 = dn; if (walk) { d2 = dn.slice(); d2[born] = dnew; } return { T: T2, val: v2, born: born + 1, mask, dn: d2 }; });
      }
    }
  };
  for (let D = 1; D <= maxDepth && best.E > target && nodes <= budget; D++) { dfs(loc, val0, 0, 0, new Float64Array(maxDepth), E0, P0, 0, D, forbid); if (walk && best.E === E0 && bestW.seqs.length) break; }
  return { E0, target, E: best.E, seqs: best.seqs, walk: bestW.seqs, glide: bestG.seqs, nodes, atRim };
}
const OPS = { flip: (M, a, b) => mFlip(M, a, b, false), collapse: mCollapse, split: mSplit };
const lonePair = (F, S) => S.length === 2 && F.valence[S[0]] + F.valence[S[1]] === 12 && Math.abs(F.valence[S[0]] - 6) === 1;
/** among sequences that reach the same valences, the one that leaves the best-shaped triangles: applied and written back.
 *  Returns {M, seq, remap} or null when none of them can be carried out. */
function applyBest(c, F, seqs) {
  let bestM = null, bw = -Infinity, bs = null;
  for (const sq_ of seqs) {
    const M = open(F); let bad = false; for (const o of sq_) if (typeof OPS[o.op](M, o.a, o.b) === 'string') { bad = true; break; }
    if (bad) continue;
    let w = Infinity; for (let t = 0; t < M.T.length; t += 3) if (M.T[t] >= 0 && (M.mark.has(M.T[t]) || M.mark.has(M.T[t + 1]) || M.mark.has(M.T[t + 2]))) w = Math.min(w, triQ(M, M.T[t], M.T[t + 1], M.T[t + 2]));
    if (w > bw) { bw = w; bestM = M; bs = sq_; }
  }
  return bestM ? { M: bestM, seq: bs, remap: close(c, bestM) } : null;
}
/** the 5s and 7s that are left where an edit was (its marked sites and the sites it started from), in the new numbering */
function leftOf(F, M, remap, sites, was = null) {
  const left = new Set(), own = new Set(sites);
  for (const i of [...M.mark, ...sites]) {
    const j = remap[i]; if (j < 0 || !F.interior[j] || F.valence[j] === 6) continue;
    if (was && !own.has(i) && i < was.length && was[i] === F.valence[j]) continue;      // was: the valences before — a 5 or 7 the edit touched and left as it was stood by
    left.add(j);
  }
  return [...left].sort((x, y) => x - y);
}
const byRim = (F, S) => { const nb = neighbours(F); return S.some(i => nb[i].some(j => !F.interior[j] && F.valence[j] > 0)); };
/**
 * healCluster(carrier, F, site): the cluster of 5s and 7s at `site` (cluster()) is taken as far down to its net charge
 * as short sequences go — a 5-7-5 to one 5, a 7-5-7 to one 7, two opposite 5-7 pairs side by side to nothing. Among the
 * sequences that reach the same valences, the one that leaves the best-shaped triangles is applied.
 * opts.rim (true): next to the rim the cluster is taken out altogether where a short sequence does it — the naked edge
 * behind a 5 split, the one behind a 7 collapsed. false leaves the rim alone.
 * opts.reach (4): a lone 5–7 pair, or a cluster at its net charge already, is taken together with the 5s and 7s within
 * that many edges, and the pair walked to them (search, neutral) — 2 leaves it alone.
 * Returns {ok, why} or {ok, before: [valences], after: [valences left], ops: ['flip', …]}.
 */
export function healCluster(c, F, site, opts = {}) {
  let sites = cluster(F, site); if (!sites.length) return { ok: false, why: 'click a 5 or a 7' };
  const reach = opts.reach ?? 4, rim = opts.rim ?? true, atRim = (S) => rim && byRim(F, S);
  const settled = (S) => { const en = energy(F, S); return !atRim(S) && (en.E === Math.abs(en.q) || lonePair(F, S)); };
  if (reach > 2 && settled(sites)) { const wide = cluster(F, site, reach); if (wide.length > sites.length) sites = wide; }
  const before = sites.map(i => F.valence[i]), ops = []; let why = '';
  for (let round = 0; round < 12 && sites.length; round++) {
    const en = energy(F, sites), near = atRim(sites);
    if (!near && en.E === Math.abs(en.q)) { why = 'this is the net charge of the cluster already: nothing to take out'; break; }
    if (!near && lonePair(F, sites)) { why = 'a lone 5–7 pair is a dislocation: no edit round it takes it out (to rim walks it out; flip moves it)'; break; }
    let s = search(F, sites, { ...opts, rim, walk: null, rings: near ? 2 : 1 });      // by the rim two rings: the naked edges that take a 5 or 7 in lie a site further out
    if (!s.seqs.length) s = search(F, sites, { maxDepth: 6, budget: opts.wideBudget ?? 30000, rings: 2, neutral: true, rim });
    const r = s.seqs.length ? applyBest(c, F, s.seqs) : null;
    if (!r) { why = 'no short sequence of flips, collapses and splits lowers this cluster'; break; }
    for (const o of r.seq) ops.push(o.op);
    const left = leftOf(F, r.M, r.remap, sites), next = new Set(); for (const j of left) for (const x of cluster(F, j)) next.add(x);
    sites = [...next].sort((x, y) => x - y);
  }
  if (!ops.length) return { ok: false, why };
  F.edits = (F.edits || 0) + 1;
  return { ok: true, before, after: sites.map(i => F.valence[i]), ops };
}
/** every site's distance to the rim, in edges; null when the lattice has no rim (a closed surface) */
export function rimDistance(F) {
  const n = F.peaks.length, nb = neighbours(F), D = new Float64Array(n).fill(-1); let front = [];
  for (let i = 0; i < n; i++) if (!F.interior[i] && F.valence[i] > 0) { D[i] = 0; front.push(i); }
  if (!front.length) return null;
  for (let d = 1; front.length; d++) { const nf = []; for (const i of front) for (const j of nb[i]) if (D[j] < 0) { D[j] = d; nf.push(j); } front = nf; }
  for (let i = 0; i < n; i++) if (D[i] < 0) D[i] = 1e6;
  return D;
}
/** In a walk of many edits the lattice is not left to the end: after each step the sites within `rings` of it are
 *  relaxed a few sweeps and stay there (their positions as frozen follow; their edges go on asking for the lengths they
 *  asked for before), so the next edit meets triangles of a sound shape. */
function settleNear(c, F, at, { rings = 3, sweeps = 8 } = {}) {
  const nb = neighbours(F), n = F.peaks.length, move = new Uint8Array(n); let front = at.filter(i => i >= 0 && i < n); for (const i of front) move[i] = 1;
  for (let r = 0; r < rings; r++) { const nf = []; for (const i of front) for (const j of nb[i]) if (!move[j]) { move[j] = 1; nf.push(j); } front = nf; }
  if (!F.restFix) F.restFix = new Map();
  for (let i = 0; i < n; i++) if (move[i]) for (const j of nb[i]) { const key = edgeKey(i, j); if (!F.restFix.has(key)) F.restFix.set(key, dist(F.P0, i, j)); }
  F._tb = null; relax(c, F, { sweeps, from: 0, mu: 0.3, move });
  for (let i = 0; i < n; i++) if (move[i]) { for (let d = 0; d < 3; d++) { F.P0[3 * i + d] = F.P[3 * i + d]; F.nrm0[3 * i + d] = F.nrm[3 * i + d]; } F.anchor0[i] = F.anchor[i]; }
  F.sweeps = 0;
}
/** `sites` walked to the rim and taken out there, as far as it goes: {ops, left: [sites left], dist, why}. Each step, in
 *  this order: a sequence that takes all of them out (within three rows of the rim, up to three edits); else the one or
 *  two edits that leave Σ (valence − 6)² as it is and bring them closer to the rim — flips alone if they do (a glide
 *  keeps every site), else with a collapse or a split; else a sequence that at least lowers the sum. */
function walkOut(c, F, start, { maxSteps = 200, hold = null, forbid = -1, until = Infinity } = {}) {
  let sites = start.slice(), why = '', dist = 0; const ops = [];
  for (let step = 0; step < maxSteps && sites.length; step++) {
    if (Date.now() > until) { why = 'stopped here: the walk was taking too long'; break; }
    const D = rimDistance(F), Ew = energy(F, sites).E, o = { rim: true, hold, forbid: step ? -1 : forbid }; dist = Math.min(...sites.map(i => D[i])); let seqs = null, part = null;
    if (dist <= 3) { const s = search(F, sites, { ...o, maxDepth: 3, budget: 15000, rings: 2 }); if (s.seqs.length) { if (s.E <= s.E0 - Ew) seqs = s.seqs; else part = s.seqs; } }
    if (!seqs) { const s = search(F, sites, { ...o, maxDepth: 2, budget: 6000, walk: D }); if (s.seqs.length && s.E <= s.E0 - Ew) seqs = s.seqs; else if (s.glide.length) seqs = s.glide; else if (s.walk.length) seqs = s.walk; else if (!part && s.seqs.length) part = s.seqs; }   // a glide before a change of row
    const was = F.valence, r = seqs || part ? applyBest(c, F, seqs || part) : null;
    if (!r) { why = sites.length === 1 ? `a lone ${F.valence[sites[0]]}, ${dist} row${dist === 1 ? '' : 's'} from the rim` : `no step takes it closer: ${dist} row${dist === 1 ? '' : 's'} from the rim`; break; }
    for (const q of r.seq) ops.push(q.op);
    if (hold) hold = new Set([...hold].map(i => r.remap[i]).filter(i => i >= 0));
    sites = leftOf(F, r.M, r.remap, sites, was);
    settleNear(c, F, [...r.M.mark].map(i => r.remap[i]));
  }
  if (sites.length) { const D = rimDistance(F); dist = Math.min(...sites.map(i => D[i])); }
  return { ops, left: sites, dist: sites.length ? dist : 0, why, hold };
}
/** A lone 5 or 7 moved one site towards the rim by one flip, which leaves a site of the same valence nearer the rim and
 *  a 5–7 pair beside it; the pair is then walked out. (A 7: an edge at it is flipped, its other end becomes a 5, the two
 *  opposite sites 7s. A 5: an edge opposite it is flipped, the ends become 5s, the far site a 7 — or nothing, when the
 *  far site is on the rim: then two 5s are left, both a row nearer.) Returns {ops, left: the 5s / 7s it leaves} or null. */
function hopOut(c, F, x, until = Infinity) {
  const D = rimDistance(F), sign = F.valence[x] - 6, E = em(open(F)), cand = [], plain = (i) => F.interior[i] && F.valence[i] === 6;
  for (const e of E.values()) {
    if (e.n !== 2 || e.bad || E.has(edgeKey(e.p, e.q))) continue; const { u, v, p, q } = e;
    if (sign > 0) {                                        // a 7 at one end of the edge
      if ((u !== x && v !== x) || !plain(u === x ? v : u) || !plain(p) || !plain(q)) continue;
      for (const k of [p, q]) if (D[k] < D[x]) cand.push({ a: u, b: v, stay: k, pair: [k === p ? q : p, u === x ? v : u], d: D[k] });
    } else {                                               // a 5 opposite the edge
      if ((p !== x && q !== x) || !plain(u) || !plain(v)) continue; const far = p === x ? q : p;
      if (plain(far)) { for (const k of [u, v]) if (D[k] < D[x]) cand.push({ a: u, b: v, stay: k, pair: [k === u ? v : u, far], d: D[k] }); }
      else if (!F.interior[far] && D[u] < D[x] && D[v] < D[x]) cand.push({ a: u, b: v, both: [u, v], d: D[u] + D[v] - 0.5 });   // the far site on the rim: no 7, two 5s
    }
  }
  cand.sort((s, t) => s.d - t.d);
  for (const k of cand.slice(0, 5)) {
    if (Date.now() > until) break;
    const snap = snapshot(F), M = open(F); if (typeof mFlip(M, k.a, k.b, false) === 'string') continue;
    const remap = close(c, M); settleNear(c, F, [...M.mark].map(i => remap[i]));
    if (k.both) { const w = walkOut(c, F, k.both.map(i => remap[i]), { until }); if (w.left.length <= 1) return { ops: ['flip', ...w.ops], left: w.left }; }
    else {
      const stay = remap[k.stay], w = walkOut(c, F, k.pair.map(i => remap[i]), { hold: new Set([stay]), forbid: edgeKey(remap[k.pair[0]], stay), until });
      if (!w.left.length) return { ops: ['flip', ...w.ops], left: [...w.hold].filter(i => F.interior[i] && F.valence[i] !== 6) };
    }
    restore(F, snap);
  }
  return null;
}
/**
 * toRim(carrier, F, site): the 5 or 7 at `site` is walked to the nearest rim and taken out there.
 * With a neighbour of the other kind it goes as a 5–7 pair: a flip glides the pair one site along its row, a collapse
 * or a split + flip takes it to the next row, each step chosen to leave Σ (valence − 6)² as it is and to bring the
 * pair closest to the rim; at the rim the naked edge behind the 5 is split and the one behind the 7 collapsed.
 * Alone it does not walk (moving it one site leaves a pair behind): within two rows of the rim it is taken out
 * directly, further in it hops — one flip moves it a site outwards and leaves a pair, the pair is walked out, and so
 * on. Other 5s and 7s on the way are not touched unless the pair can join with them.
 * It gives up after `ms` milliseconds and keeps the furthest it got.
 * Returns {ok, why} or {ok, before, after: [valences left], ops, out (nothing left), dist (rows left to the rim), why}.
 */
export function toRim(c, F, site, { maxSteps = 200, hops = 16, ms = 1500 } = {}) {
  const sing = (i) => F.interior[i] && F.valence[i] !== 6; if (site < 0 || !sing(site)) return { ok: false, why: 'click a 5 or a 7' };
  const D0 = rimDistance(F); if (!D0) return { ok: false, why: 'this lattice has no rim: on a closed surface a 5–7 pair can only be walked to another' };
  const sign = F.valence[site] - 6, mates = neighbours(F)[site].filter(j => sing(j) && (F.valence[j] - 6) * sign < 0).sort((x, y) => D0[x] - D0[y]);
  const edits = F.edits || 0, start = snapshot(F), until = Date.now() + ms; let best = null, bestSnap = null;
  const better = (r) => !best || (r.out && !best.out) || (r.out === best.out && (r.dist < best.dist || (r.dist === best.dist && r.ops.length < best.ops.length)));
  for (const walker of [...mates.map(j => [site, j]), [site]]) {
    if (best && Date.now() > until) break;
    restore(F, start); const before = walker.map(i => F.valence[i]); let w = walkOut(c, F, walker, { maxSteps, until }), ops = w.ops.slice();
    for (let h = 0; h < hops && w.left.length === 1 && Date.now() <= until; h++) {  // what is left is one 5 or 7: hop by hop
      const hop = hopOut(c, F, w.left[0], until); if (!hop) break;
      ops.push(...hop.ops); w = walkOut(c, F, hop.left, { maxSteps, until }); ops.push(...w.ops);
    }
    const r = { ok: ops.length > 0, before, after: w.left.map(i => F.valence[i]), ops, out: !w.left.length, dist: w.dist, why: w.why };
    if (r.ok && better(r)) { best = r; bestSnap = snapshot(F); if (r.out) break; }
    if (!r.ok && !best) best = r;
  }
  if (!best || !best.ok) { restore(F, start); return { ok: false, why: best && best.why ? `it does not move: ${best.why}` : 'no edit moves it towards the rim' }; }
  restore(F, bestSnap); F.edits = edits + 1;
  return best;
}
/** every cluster that holds more than its net charge, one after the other (the rim is left alone: opts.rim false). Returns {ok, clusters, ops, singular: [before, after]} */
export function healAll(c, F, opts = {}) {
  const count = () => { let s = 0; for (let i = 0; i < F.valence.length; i++) if (F.interior[i] && F.valence[i] !== 6) s++; return s; };
  const failed = new Set(), s0 = count(), edits = F.edits || 0; let clusters = 0, ops = 0;
  for (let it = 0; it < 2000; it++) {
    const n = F.peaks.length, seen = new Uint8Array(n); let pick = -1;
    for (let i = 0; i < n && pick < 0; i++) {
      if (seen[i] || !F.interior[i] || F.valence[i] === 6) continue;
      const cl = cluster(F, i), en = energy(F, cl); for (const x of cl) seen[x] = 1;
      if (cl.length < 2 || en.E === Math.abs(en.q) || lonePair(F, cl) || failed.has(cl.map(x => F.peaks[x]).join(','))) continue;
      pick = i;
    }
    if (pick < 0) break;
    const key = cluster(F, pick).map(x => F.peaks[x]).join(','), r = healCluster(c, F, pick, { budget: 20000, wideBudget: 10000, reach: 2, rim: false, ...opts });
    if (r.ok) { clusters++; ops += r.ops.length; } else failed.add(key);
  }
  F.edits = edits + (clusters ? 1 : 0);
  return clusters ? { ok: true, clusters, ops, singular: [s0, count()] } : { ok: false, why: 'no cluster holds more than its net charge, or none gives way to a short edit' };
}

// ------------------------------------------------------------------ the rim: fans completed to six triangles
function rimOf(M) {
  const next = new Map(), prev = new Map(), multi = new Set();
  for (const e of em(M).values()) {
    if (e.n !== 1) { if (e.bad) { multi.add(e.u); multi.add(e.v); } continue; }
    if (next.has(e.u)) multi.add(e.u); else next.set(e.u, e.v);
    if (prev.has(e.v)) multi.add(e.v); else prev.set(e.v, e.u);
  }
  return { next, prev, multi };
}
/** the sites round v, counter-clockwise from the rim neighbour its rim edge runs to (w) to the one it comes from */
function fan(M, v, w) {
  const E = em(M), ring = [w]; let x = w;
  for (let k = 0; k < 12; k++) { const e = E.get(edgeKey(v, x)); if (!e) return null; const y = e.u === v ? e.p : e.q; if (y < 0) break; ring.push(y); x = y; }
  return ring;
}
/** The fan of the rim site v completed to six triangles. With t triangles there, 6 − t are missing: one (the notch is
 *  closed by the triangle of v and its two rim neighbours) or more (5 − t new sites). The new sites continue the
 *  lattice: round a regular site every neighbour offset is the one before minus the one before that, taken in the
 *  tangent plane at v from both ends of the gap and blended. Each goes onto the carrier; one that falls past the
 *  boundary (up to o.beyond edge lengths) goes onto the carrier continued there (fitExt). The fan is refused when a
 *  site lands on another, when a triangle is of the wrong size or sense, and — unless o.force — when the carrier ends
 *  less than o.gap edge lengths outside both rim edges at v (there is nothing left to cover). */
function mFill(c, M, v, R, out, o) {
  if (!R.next.has(v) || !R.prev.has(v) || R.multi.has(v)) return 'not a site of the rim';
  const w = R.next.get(v), u = R.prev.get(v), ring = fan(M, v, w);
  if (!ring || ring[ring.length - 1] !== u) return 'the triangles at this site are not one fan';
  const t = ring.length - 1, g = 6 - t; if (g < 1) return 'this site has six triangles already';
  const P = M.P, px = P[3 * v], py = P[3 * v + 1], pz = P[3 * v + 2], nl = Math.hypot(M.N[3 * v], M.N[3 * v + 1], M.N[3 * v + 2]) || 1, nx = M.N[3 * v] / nl, ny = M.N[3 * v + 1] / nl, nz = M.N[3 * v + 2] / nl;
  let L = 0; for (const x of ring) L += restOf(M, v, x); L /= ring.length;
  const okLen = (l) => l > 0.45 * L && l < 1.7 * L;
  const off = (x, y, z) => { project(c, M.A[v], x, y, z, out); const ex = x - out[0], ey = y - out[1], ez = z - out[2], en = ex * out[3] + ey * out[4] + ez * out[5]; return Math.hypot(ex - en * out[3], ey - en * out[4], ez - en * out[5]); };   // how far past the boundary a point lies (0 over the carrier)
  // is there carrier left uncovered beside v? Points o.gap edge lengths outside the two rim edges at v are tried, near v and at the middle of each
  const outside = (b, opp, tau) => {
    const ex = P[3 * b] - px, ey = P[3 * b + 1] - py, ez = P[3 * b + 2] - pz, e2 = ex * ex + ey * ey + ez * ez || 1, mx = px + tau * ex, my = py + tau * ey, mz = pz + tau * ez;
    let dx = mx - P[3 * opp], dy = my - P[3 * opp + 1], dz = mz - P[3 * opp + 2]; const de = (dx * ex + dy * ey + dz * ez) / e2; dx -= de * ex; dy -= de * ey; dz -= de * ez; const s = o.gap * L / (Math.hypot(dx, dy, dz) || 1);
    return off(mx + s * dx, my + s * dy, mz + s * dz) <= 0.03 * L;
  };
  const open_ = () => outside(u, ring[t - 1], 0.5) || outside(w, ring[1], 0.5) || outside(u, ring[t - 1], 0.12) || outside(w, ring[1], 0.12);
  if (g === 1) {
    if (u === w || em(M).has(edgeKey(u, w))) return 'no room for a triangle here';
    if (!okLen(dist(P, u, w)) || !(triQ(M, v, u, w) > 0.15)) return 'the gap at this site is not one triangle wide';
    if (!o.force && !open_()) return 'nothing of the surface is left to cover here';
    M.T.push(v, u, w); M.val[v]++; M.val[u]++; M.val[w]++; M.fix.set(edgeKey(u, w), L); for (const x of [v, u, w]) M.mark.add(x);
    return { sites: 0, tris: 1, beyond: 0, touched: [v, u, w] };
  }
  if (!o.force && !open_()) return 'nothing of the surface is left to cover here';
  const D = ring.map(x => {                               // offsets in the tangent plane at v, each as long as its edge
    let dx = P[3 * x] - px, dy = P[3 * x + 1] - py, dz = P[3 * x + 2] - pz; const l = Math.hypot(dx, dy, dz), dn = dx * nx + dy * ny + dz * nz;
    dx -= dn * nx; dy -= dn * ny; dz -= dn * nz; const s = l / (Math.hypot(dx, dy, dz) || 1); return [dx * s, dy * s, dz * s];
  });
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], f = [], b = [];
  let a1 = D[t], a0 = D[t - 1]; for (let j = 1; j < g; j++) { const d = sub(a1, a0); f.push(d); a0 = a1; a1 = d; }
  a1 = D[0]; a0 = D[1]; for (let j = g - 1; j >= 1; j--) { const d = sub(a1, a0); b[j - 1] = d; a0 = a1; a1 = d; }
  const X = [], NN = [], AA = [], EE = []; let beyond = 0;
  for (let j = 1; j < g; j++) {
    const al = g > 2 ? (j - 1) / (g - 2) : 0.5, x = px + (1 - al) * f[j - 1][0] + al * b[j - 1][0], y = py + (1 - al) * f[j - 1][1] + al * b[j - 1][1], z = pz + (1 - al) * f[j - 1][2] + al * b[j - 1][2];
    const d = off(x, y, z); let E = null;
    if (d > 0.12 * L) {                                    // past the boundary: on the carrier continued from there
      if (d > o.beyond * L) return 'the next site would lie too far past the boundary of the surface';
      const an = out[6]; E = fitExt(c, an, out[0], out[1], out[2], out[3], out[4], out[5], 2 * L); liftExt(E, x, y, z, out); out[6] = an; beyond++;
    }
    for (let i = 0; i < M.n; i++) if (!M.dead.has(i) && Math.hypot(P[3 * i] - out[0], P[3 * i + 1] - out[1], P[3 * i + 2] - out[2]) < 0.55 * L) return 'another site is in the way';
    for (let i = 0; i < X.length; i += 3) if (Math.hypot(X[i] - out[0], X[i + 1] - out[1], X[i + 2] - out[2]) < 0.55 * L) return 'the gap at this site is too narrow for six triangles';
    X.push(out[0], out[1], out[2]); NN.push(out[3], out[4], out[5]); AA.push(out[6]); EE.push(E);
  }
  const first = M.n, seq = [u]; for (let j = 0; j < g - 1; j++) seq.push(first + j); seq.push(w);
  P.push(...X); M.N.push(...NN); M.A.push(...AA);                           // tried in place, taken back when a triangle is wrong
  let bad = false;
  for (let k = 0; k < g && !bad; k++) bad = !okLen(dist(P, seq[k], seq[k + 1])) || !okLen(dist(P, v, seq[k])) || !okLen(dist(P, v, seq[k + 1])) || !(triQ(M, v, seq[k], seq[k + 1]) > 0.15);
  if (bad) { P.length = 3 * first; M.N.length = 3 * first; M.A.length = first; return 'the gap at this site does not take six regular triangles'; }
  M.n += g - 1;
  for (let k = 0; k < g; k++) { M.T.push(v, seq[k], seq[k + 1]); M.fix.set(edgeKey(seq[k], seq[k + 1]), L); M.fix.set(edgeKey(v, seq[k]), L); M.fix.set(edgeKey(v, seq[k + 1]), L); }
  M.val[v] += g; M.val[u]++; M.val[w]++;
  for (let j = 0; j < g - 1; j++) { M.val[first + j] = 2; M.gen.set(first + j, (M.gen.get(v) || 0) + 1); if (EE[j]) M.ext.set(first + j, EE[j]); }
  for (const x of seq) M.mark.add(x); M.mark.add(v);
  return { sites: g - 1, tris: g, beyond, touched: [v, ...seq] };
}
function fillPasses(c, M, want, maxPasses, o) {
  let sites = 0, tris = 0, beyond = 0, why = ''; const blocked = new Set(), out = new Float64Array(7);
  for (let pass = 0; pass < maxPasses; pass++) {
    M._em = null; const R = rimOf(M), cand = [];
    for (const v of R.next.keys()) { const g = 6 - M.val[v]; if (want(v, g) && !blocked.has(8 * v + g)) cand.push([g, v]); }
    cand.sort((x, y) => x[0] - y[0] || x[1] - y[1]);                          // the notches first
    const dirty = new Set(); let done = 0;
    for (const [g, v] of cand) {
      if (dirty.has(v)) continue;
      const r = mFill(c, M, v, R, out, o); if (typeof r === 'string') { blocked.add(8 * v + g); why = why || r; continue; }
      for (const x of r.touched) dirty.add(x); sites += r.sites; tris += r.tris; beyond += r.beyond; done++;
    }
    if (!done) break;
  }
  M._em = null; return { sites, tris, beyond, why };
}
/** the fan of one rim site completed to six triangles, whether or not there is carrier left to cover (new sites up to
 *  `beyond` edge lengths past the boundary). Returns {ok, why} or {ok, sites, tris, beyond} (what was added; beyond: the sites past the boundary) */
export function patchSite(c, F, site, { beyond = 1.5 } = {}) {
  const M = open(F), r = fillPasses(c, M, (v) => v === site, 1, { beyond, gap: 0, force: true });
  if (!r.tris) return { ok: false, why: r.why || 'not a site of the rim' };
  close(c, M); F.edits = (F.edits || 0) + 1; return { ok: true, sites: r.sites, tris: r.tris, beyond: r.beyond };
}
/** The rim grown until the lattice covers the carrier: the fan of every rim site is completed to six triangles as long
 *  as the carrier goes on more than `gap` edge lengths outside a rim edge there — the notches of a staircase rim first,
 *  then whole rows — at most `gens` rows deep in one call (again for more). maxGap 2 fills notches only (rim sites with
 *  four or five triangles). The new sites are regular; those past the boundary lie on the carrier continued there (at
 *  most `beyond` edge lengths out). */
export function patchRim(c, F, { maxGap = 5, gens = 2, beyond = 1.2, gap = 0.15 } = {}) {
  const M = open(F), r = fillPasses(c, M, (v, g) => g >= 1 && g <= maxGap && (M.gen.get(v) || 0) < gens, 4000, { beyond, gap, force: false });
  if (!r.tris) return { ok: false, why: r.why ? `nothing to patch at the rim (${r.why})` : 'nothing to patch at the rim' };
  close(c, M); F.edits = (F.edits || 0) + 1; return { ok: true, sites: r.sites, tris: r.tris, beyond: r.beyond };
}

// ------------------------------------------------------------------ picking, undo
/** the lattice edge nearest to a point (within maxDist), as {a, b, d} */
export function nearestEdge(F, x, y, z, maxDist = Infinity) {
  const T = F.tris, P = F.P; let best = null, bd = maxDist;
  for (let t = 0; t < T.length; t += 3) for (let k = 0; k < 3; k++) {
    const a = T[t + k], b = T[t + (k + 1) % 3]; if (a > b && F.interior[a] && F.interior[b]) continue;      // an inner edge is seen from both sides: once is enough
    const ax = P[3 * a], ay = P[3 * a + 1], az = P[3 * a + 2], ex = P[3 * b] - ax, ey = P[3 * b + 1] - ay, ez = P[3 * b + 2] - az;
    const s = Math.max(0, Math.min(1, ((x - ax) * ex + (y - ay) * ey + (z - az) * ez) / (ex * ex + ey * ey + ez * ez || 1))), d = Math.hypot(ax + s * ex - x, ay + s * ey - y, az + s * ez - z);
    if (d < bd) { bd = d; best = { a, b, d }; }
  }
  return best;
}
/** the site nearest to a point among those `keep` lets through (within maxDist), or −1 */
export function nearestSite(F, x, y, z, keep = null, maxDist = Infinity) {
  const P = F.P; let best = -1, bd = maxDist;
  for (let i = 0; i < F.peaks.length; i++) { if (keep && !keep(i)) continue; const d = Math.hypot(P[3 * i] - x, P[3 * i + 1] - y, P[3 * i + 2] - z); if (d < bd) { bd = d; best = i; } }
  return best;
}
export const isRimSite = (F, i) => !F.interior[i] && F.valence[i] > 0;
/** what an edit changes, to put back with restore */
export function snapshot(F) {
  return { peaks: F.peaks.slice(), P0: F.P0.slice(), nrm0: F.nrm0.slice(), anchor0: F.anchor0.slice(), shift: F.shift ? F.shift.slice() : null, tris: F.tris.slice(), restFix: F.restFix ? new Map(F.restFix) : null, ext: F.ext ? new Map(F.ext) : null, edits: F.edits || 0 };
}
export function restore(F, s) {
  Object.assign(F, { peaks: s.peaks.slice(), P0: s.P0.slice(), nrm0: s.nrm0.slice(), anchor0: s.anchor0.slice(), shift: s.shift ? s.shift.slice() : null, tris: s.tris.slice(), restFix: s.restFix ? new Map(s.restFix) : null, ext: s.ext ? new Map(s.ext) : null, edits: s.edits });
  F.P = Float64Array.from(F.P0); F.nrm = Float64Array.from(F.nrm0); F.anchor = Int32Array.from(F.anchor0); F.sweeps = 0;
  return derive(F);
}
