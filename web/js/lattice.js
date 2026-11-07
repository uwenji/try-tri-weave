// Lattice read-out: peaks of u → geodesic Voronoi dual on the carrier → triangulation among peaks → valences
// (closed-link test), and the Kagome strands (laths) of that triangulation with crossing parity and 3-colouring.
import { dijkstra } from './sparse.js';

export function localMaxima(u, G) {
  const n = G.n, out = [];
  for (let i = 0; i < n; i++) {
    let ok = true;
    for (let k = G.indptr[i]; k < G.indptr[i + 1]; k++) if (u[G.indices[k]] >= u[i]) { ok = false; break; }
    if (ok) out.push(i);
  }
  return out;
}

/** returns null when there is no pattern yet; else {peaks, P, tris, valence, interior, label, threshold} */
export function extractLattice(carrier, u, a, { relThreshold = 0.4, voidRel = 0.15 } = {}) {
  const c = carrier, n = c.n;
  let umax = -Infinity; for (let i = 0; i < n; i++) if (u[i] > umax) umax = u[i];
  if (!(umax > 1e-6)) return null;
  const thr = relThreshold * umax;
  const pk = localMaxima(u, c.graph).filter(i => u[i] > thr);
  if (pk.length < 3) return null;
  // patterned region: rms of u smoothed over ~a third of the lattice edge
  const rad = 0.35 * a, rhs = new Float64Array(n); for (let i = 0; i < n; i++) rhs[i] = c.m[i] * u[i] * u[i];
  const amp = c.heatSolve(rhs, rad * rad);
  const voidMask = new Uint8Array(n); for (let i = 0; i < n; i++) voidMask[i] = (Math.sqrt(Math.max(amp[i], 0)) < voidRel * umax || c.boundary[i]) ? 1 : 0;
  // geodesic Voronoi labels
  const { source } = dijkstra(c.graph, pk);
  const pos = new Int32Array(n).fill(-1); pk.forEach((p, i) => { pos[p] = i; });
  const label = new Int32Array(n); for (let i = 0; i < n; i++) label[i] = source[i] >= 0 ? pos[source[i]] : -1;
  // dual triangles: intrinsic faces with three different labels, none void
  const Fi = c.Fi, nf = Fi.length / 3, triSet = new Map();
  for (let f = 0; f < nf; f++) {
    const i0 = Fi[3 * f], i1 = Fi[3 * f + 1], i2 = Fi[3 * f + 2];
    const l0 = label[i0], l1 = label[i1], l2 = label[i2];
    if (l0 < 0 || l1 < 0 || l2 < 0 || l0 === l1 || l1 === l2 || l0 === l2) continue;
    if (voidMask[i0] || voidMask[i1] || voidMask[i2]) continue;
    const s = [l0, l1, l2].sort((x, y) => x - y);
    triSet.set(s[0] * 1e10 + s[1] * 1e5 + s[2], s);
  }
  const tris = Array.from(triSet.values());
  const npk = pk.length, valence = new Int32Array(npk), interior = new Uint8Array(npk);
  for (const t of tris) { valence[t[0]]++; valence[t[1]]++; valence[t[2]]++; }
  // closed link: every dual edge at a peak is shared by exactly two dual triangles
  const edgeCount = new Map();
  for (const [p, q, r] of tris) for (const [x, y] of [[p, q], [q, r], [p, r]]) { const k = x * 1e5 + y; edgeCount.set(k, (edgeCount.get(k) || 0) + 1); }
  const bad = new Uint8Array(npk);
  for (const [k, cnt] of edgeCount) if (cnt !== 2) { bad[Math.floor(k / 1e5)] = 1; bad[k % 1e5] = 1; }
  for (let i = 0; i < npk; i++) interior[i] = (!bad[i] && valence[i] >= 3) ? 1 : 0;
  // positions and orientation (like the carrier normals)
  const P = new Float64Array(npk * 3); pk.forEach((p, i) => { P[3 * i] = c.V[3 * p]; P[3 * i + 1] = c.V[3 * p + 1]; P[3 * i + 2] = c.V[3 * p + 2]; });
  const T = new Int32Array(tris.length * 3);
  tris.forEach((t, k) => {
    const [p, q, r] = t;
    const ax = P[3 * q] - P[3 * p], ay = P[3 * q + 1] - P[3 * p + 1], az = P[3 * q + 2] - P[3 * p + 2];
    const bx = P[3 * r] - P[3 * p], by = P[3 * r + 1] - P[3 * p + 1], bz = P[3 * r + 2] - P[3 * p + 2];
    const nx = ay * bz - az * by, ny = az * bx - ax * bz, nz = ax * by - ay * bx;
    const N = c.N, v = pk[p], w = pk[q], x = pk[r];
    const mnx = N[3 * v] + N[3 * w] + N[3 * x], mny = N[3 * v + 1] + N[3 * w + 1] + N[3 * x + 1], mnz = N[3 * v + 2] + N[3 * w + 2] + N[3 * x + 2];
    const flip = nx * mnx + ny * mny + nz * mnz < 0;
    T[3 * k] = p; T[3 * k + 1] = flip ? r : q; T[3 * k + 2] = flip ? q : r;
  });
  return { peaks: Int32Array.from(pk), P, tris: T, valence, interior, label, threshold: thr, voidMask };
}

export function chargeSummary(lat) {
  if (!lat) return { hist: {}, charge: 0, singular: [], nInterior: 0 };
  const hist = {}, singular = []; let charge = 0, nInterior = 0;
  for (let i = 0; i < lat.valence.length; i++) {
    if (!lat.interior[i]) continue;
    const v = lat.valence[i]; hist[v] = (hist[v] || 0) + 1; charge += 6 - v; nInterior++;
    if (v !== 6) singular.push(i);
  }
  return { hist, charge, singular, nInterior };
}

/** Kagome strands of a triangulation: each lath passes through the midpoints of lattice edges, taking in every
 *  triangle the mid-segment that wraps the vertex it did *not* wrap in the previous triangle (a straight line).
 *  Returns {edges: [[a,b],...], strands: [{knots: [edgeIds], closed}], family: Int32Array, over: (strand, i) → bool,
 *           parityConflicts, colourConflicts}. */
export function traceStrands(npk, tris) {
  const nt = tris.length / 3, edgeId = new Map(), edges = [], triEdges = new Int32Array(nt * 3), edgeTris = [];
  const key = (a, b) => (a < b ? a * 1e5 + b : b * 1e5 + a);
  for (let t = 0; t < nt; t++) for (let k = 0; k < 3; k++) {
    const a = tris[3 * t + k], b = tris[3 * t + (k + 1) % 3], kk = key(a, b);
    let e = edgeId.get(kk);
    if (e === undefined) { e = edges.length; edgeId.set(kk, e); edges.push([Math.min(a, b), Math.max(a, b)]); edgeTris.push([]); }
    triEdges[3 * t + k] = e;                     // edge k of triangle t joins corners k and k+1
    edgeTris[e].push(t);
  }
  const cornerOf = (t, v) => tris[3 * t] === v ? 0 : tris[3 * t + 1] === v ? 1 : tris[3 * t + 2] === v ? 2 : -1;
  // segment (t, k) wraps corner k and joins the midpoints of edges (k-1,k) and (k,k+1): edge ids triEdges[t,(k+2)%3] and triEdges[t,k]
  const segEdges = (t, k) => [triEdges[3 * t + (k + 2) % 3], triEdges[3 * t + k]];
  const used = new Uint8Array(nt * 3), strands = [];
  const otherTri = (e, t) => { const l = edgeTris[e]; return l.length === 2 ? (l[0] === t ? l[1] : l[0]) : -1; };
  const advance = (t, k, exitEdge) => {
    // leave triangle t through exitEdge (one of the two edges at corner k); in the next triangle wrap the other vertex of that edge
    const t2 = otherTri(exitEdge, t); if (t2 < 0) return null;
    const v = tris[3 * t + k], [ea, eb] = edges[exitEdge], w = ea === v ? eb : ea;
    const k2 = cornerOf(t2, w); if (k2 < 0) return null;
    const [e1, e2] = segEdges(t2, k2); const nextExit = e1 === exitEdge ? e2 : e1;
    return { t: t2, k: k2, entry: exitEdge, exit: nextExit };
  };
  for (let t0 = 0; t0 < nt; t0++) for (let k0 = 0; k0 < 3; k0++) {
    if (used[3 * t0 + k0]) continue;
    const [eA, eB] = segEdges(t0, k0);
    // forward: exit through eB; backward: exit through eA
    const fwd = [], bwd = []; let closed = false;
    used[3 * t0 + k0] = 1;
    let cur = { t: t0, k: k0, exit: eB };
    for (;;) {
      const nx = advance(cur.t, cur.k, cur.exit); if (!nx) break;
      if (nx.t === t0 && nx.k === k0) { closed = true; break; }
      if (used[3 * nx.t + nx.k]) break;
      used[3 * nx.t + nx.k] = 1; fwd.push(nx); cur = nx;
    }
    if (!closed) {
      cur = { t: t0, k: k0, exit: eA };
      for (;;) {
        const nx = advance(cur.t, cur.k, cur.exit); if (!nx) break;
        if (used[3 * nx.t + nx.k]) break;
        used[3 * nx.t + nx.k] = 1; bwd.push(nx); cur = nx;
      }
    }
    // knots in order: backward chain reversed … start … forward chain
    const knots = [];
    for (let i = bwd.length - 1; i >= 0; i--) knots.push(bwd[i].exit);       // the far edge of each backward segment
    knots.push(eA, eB);
    for (const s of fwd) knots.push(s.exit);
    if (closed) knots.pop();                                                 // last exit equals the first knot
    strands.push({ knots, closed, segs: [...bwd.slice().reverse().map(s => [s.t, s.k]), [t0, k0], ...fwd.map(s => [s.t, s.k])] });
  }
  // crossings: each interior lattice edge is a knot shared by two strands
  const atEdge = edges.map(() => []);
  strands.forEach((s, si) => s.knots.forEach((e, i) => atEdge[e].push([si, i])));
  // 3-colouring (families): strands crossing at a knot get different colours (greedy, BFS)
  const ns = strands.length, family = new Int32Array(ns).fill(-1), adj = Array.from({ length: ns }, () => new Set());
  for (const l of atEdge) if (l.length === 2 && l[0][0] !== l[1][0]) { adj[l[0][0]].add(l[1][0]); adj[l[1][0]].add(l[0][0]); }
  let colourConflicts = 0;
  const order = [...Array(ns).keys()].sort((a, b) => strands[b].knots.length - strands[a].knots.length);
  for (const root of order) {
    if (family[root] >= 0) continue;
    const queue = [root]; family[root] = 0;
    while (queue.length) {
      const s = queue.shift();
      for (const q of adj[s]) {
        if (family[q] >= 0) continue;
        const usedC = new Set(); for (const z of adj[q]) if (family[z] >= 0) usedC.add(family[z]);
        let col = [0, 1, 2].find(c => !usedC.has(c));
        if (col === undefined) { col = (family[s] + 1) % 3; colourConflicts++; }
        family[q] = col; queue.push(q);
      }
    }
  }
  // min-conflict sweeps: recolour a strand when another family has fewer conflicting crossings
  for (let sweep = 0; sweep < 6; sweep++) {
    let changed = 0;
    for (let s = 0; s < ns; s++) {
      const cnt = [0, 0, 0]; for (const q of adj[s]) if (family[q] >= 0) cnt[family[q]]++;
      let best = family[s]; for (let c = 0; c < 3; c++) if (cnt[c] < cnt[best]) best = c;
      if (best !== family[s]) { family[s] = best; changed++; }
    }
    if (!changed) break;
  }
  colourConflicts = 0;
  for (let s = 0; s < ns; s++) for (const q of adj[s]) if (q > s && family[q] === family[s]) colourConflicts++;
  // over/under parity: over(s, i) = base[s] XOR (i & 1); at a knot the two strands must differ
  const base = new Int8Array(ns).fill(-1); let parityConflicts = 0;
  const padj = Array.from({ length: ns }, () => []);
  for (const l of atEdge) if (l.length === 2 && l[0][0] !== l[1][0]) {
    const [[sa, ia], [sb, ib]] = l; const rel = 1 ^ ((ia + ib) & 1);   // base[sa] XOR base[sb] must equal rel
    padj[sa].push([sb, rel]); padj[sb].push([sa, rel]);
  }
  for (let root = 0; root < ns; root++) {
    if (base[root] >= 0) continue; base[root] = 0; const queue = [root];
    while (queue.length) {
      const s = queue.shift();
      for (const [q, rel] of padj[s]) {
        const want = base[s] ^ rel;
        if (base[q] < 0) { base[q] = want; queue.push(q); } else if (base[q] !== want) parityConflicts++;
      }
    }
  }
  for (const s of strands) if (s.closed && (s.knots.length & 1)) parityConflicts++;
  const over = (si, i) => ((base[si] ^ (i & 1)) === 1);
  return { edges, strands, family, over, parityConflicts: parityConflicts / 2 | 0, colourConflicts };
}

// ------------------------------------------------------------------ per-direction read-outs
/** Direction family of every lattice edge: the two laths crossing at its midpoint knot belong to two families, the
 *  edge runs along the third. Returns Int8Array over strandsData.edges (−1 where the colouring is inconsistent). */
export function edgeFamilies(strandsData) {
  const { edges, strands, family } = strandsData, nE = edges.length, fam = new Int8Array(nE).fill(-1), seen = new Int8Array(nE).fill(-1);
  for (let si = 0; si < strands.length; si++) {
    const f = family[si]; if (f < 0) continue;
    for (const e of strands[si].knots) {
      if (seen[e] < 0) seen[e] = f;
      else if (seen[e] !== f && fam[e] < 0) fam[e] = 3 - seen[e] - f;
      else if (seen[e] === f) fam[e] = -2;                                 // same family twice: conflict
    }
  }
  for (let e = 0; e < nE; e++) if (fam[e] === -2) fam[e] = -1;
  return fam;
}

/** Lattice edge lengths (chord between the peaks) and strain relative to the target edge: a, or a·e^{−φ} at the
 *  edge midpoint when a size field is active. Returns {len, strain, fam, target} over strandsData.edges. */
export function edgeStrain(lat, strandsData, a, phi = null) {
  const { edges } = strandsData, P = lat.P, pk = lat.peaks, nE = edges.length;
  const len = new Float64Array(nE), strain = new Float64Array(nE), target = new Float64Array(nE);
  for (let e = 0; e < nE; e++) {
    const [i, j] = edges[e];
    len[e] = Math.hypot(P[3 * i] - P[3 * j], P[3 * i + 1] - P[3 * j + 1], P[3 * i + 2] - P[3 * j + 2]);
    target[e] = phi ? a * Math.exp(-0.5 * (phi[pk[i]] + phi[pk[j]])) : a;
    strain[e] = len[e] / target[e] - 1;
  }
  return { len, strain, target, fam: edgeFamilies(strandsData) };
}

/** Lattice orientation at every interior 6-valent peak: θ = arg(Σ_k e^{6iα_k}) / 6 ∈ (−30°, 30°], α_k = angles of
 *  the six neighbour directions in the tangent plane, measured from the global x axis projected onto the surface
 *  (y where the normal is nearly ±x). Exact on flat regions; on curved ones the reference turns smoothly, grain
 *  boundaries still show as jumps. Returns {theta (deg, NaN elsewhere), coherence = |⟨e^{6iθ}⟩|, count}. */
export function orientationField(lat, carrier, startPeak = -1) {
  const { peaks, P, tris, valence, interior } = lat, npk = peaks.length, N = carrier.N;
  const nbr = Array.from({ length: npk }, () => []);
  for (let t = 0; t < tris.length / 3; t++) { const p = tris[3 * t], q = tris[3 * t + 1], r = tris[3 * t + 2]; nbr[p].push(q, r); nbr[q].push(p, r); nbr[r].push(p, q); }
  const theta = new Float64Array(npk).fill(NaN), dir = new Float64Array(3 * npk); let cx = 0, cy = 0, count = 0;
  for (let i = 0; i < npk; i++) {
    if (!interior[i] || valence[i] !== 6) continue;
    const v = peaks[i], nx = N[3 * v], ny = N[3 * v + 1], nz = N[3 * v + 2];
    let ex, ey, ez;                                                           // reference direction in the tangent plane
    if (Math.abs(nx) < 0.9) { ex = 1 - nx * nx; ey = -nx * ny; ez = -nx * nz; } else { ex = -ny * nx; ey = 1 - ny * ny; ez = -ny * nz; }
    const el = Math.hypot(ex, ey, ez) || 1; ex /= el; ey /= el; ez /= el;
    const fx = ny * ez - nz * ey, fy = nz * ex - nx * ez, fz = nx * ey - ny * ex;   // e2 = N × e1
    let sx = 0, sy = 0;
    for (const j of new Set(nbr[i])) {
      const dx = P[3 * j] - P[3 * i], dy = P[3 * j + 1] - P[3 * i + 1], dz = P[3 * j + 2] - P[3 * i + 2];
      const al = Math.atan2(dx * fx + dy * fy + dz * fz, dx * ex + dy * ey + dz * ez);
      sx += Math.cos(6 * al); sy += Math.sin(6 * al);
    }
    const th = Math.atan2(sy, sx) / 6; theta[i] = th * 180 / Math.PI;
    dir[3 * i] = Math.cos(th) * ex + Math.sin(th) * fx; dir[3 * i + 1] = Math.cos(th) * ey + Math.sin(th) * fy; dir[3 * i + 2] = Math.cos(th) * ez + Math.sin(th) * fz;
    cx += Math.cos(6 * th); cy += Math.sin(6 * th); count++;
  }
  // misorientation: mean |Δθ| (mod 60°) to the regular neighbours, comparing directions after projecting the
  // neighbour's onto this site's tangent plane (a discrete parallel transport over one lattice edge) — 0 inside a
  // grain on any curved surface, large along seams, unlike θ itself, whose reference (global x) drifts with the surface
  const mis = new Float64Array(npk).fill(NaN);
  for (let i = 0; i < npk; i++) {
    if (!Number.isFinite(theta[i])) continue;
    const v = peaks[i], nx = N[3 * v], ny = N[3 * v + 1], nz = N[3 * v + 2];
    const dx = dir[3 * i], dy = dir[3 * i + 1], dz = dir[3 * i + 2], gx = ny * dz - nz * dy, gy = nz * dx - nx * dz, gz = nx * dy - ny * dx;
    let s = 0, k = 0;
    for (const j of new Set(nbr[i])) {
      if (!Number.isFinite(theta[j])) continue;
      const ux = dir[3 * j], uy = dir[3 * j + 1], uz = dir[3 * j + 2];
      let d = Math.atan2(ux * gx + uy * gy + uz * gz, ux * dx + uy * dy + uz * dz) * 180 / Math.PI;
      d = ((d % 60) + 90) % 60 - 30; s += Math.abs(d); k++;
    }
    if (k) mis[i] = s / k;
  }
  // arrows: a lattice direction is defined only mod 60°, so one arrow per site is a choice. Make it consistent: start at
  // `startPeak` (the regular site nearest the seed centre) with the direction closest to the reference, then walk the
  // regular sites breadth-first and give each the one of its six directions closest to its parent's arrow transported
  // into its tangent plane. Inside a grain the arrows are then parallel; around a disclination (charge not a multiple
  // of 6) a cut where they jump by 60° is unavoidable — the misorientation colour stays quiet there.
  const arrow = new Float64Array(3 * npk), seenA = new Uint8Array(npk);
  const choose = (j, tx, ty, tz) => {
    const v = peaks[j], nx = N[3 * v], ny = N[3 * v + 1], nz = N[3 * v + 2], dn = tx * nx + ty * ny + tz * nz;
    tx -= dn * nx; ty -= dn * ny; tz -= dn * nz;
    const dx = dir[3 * j], dy = dir[3 * j + 1], dz = dir[3 * j + 2], gx = ny * dz - nz * dy, gy = nz * dx - nx * dz, gz = nx * dy - ny * dx;
    const k = Math.round(Math.atan2(tx * gx + ty * gy + tz * gz, tx * dx + ty * dy + tz * dz) / (Math.PI / 3)) * Math.PI / 3;
    arrow[3 * j] = Math.cos(k) * dx + Math.sin(k) * gx; arrow[3 * j + 1] = Math.cos(k) * dy + Math.sin(k) * gy; arrow[3 * j + 2] = Math.cos(k) * dz + Math.sin(k) * gz;
  };
  const starts = Number.isFinite(theta[startPeak]) ? [startPeak] : [];
  for (let i = 0; i < npk; i++) if (Number.isFinite(theta[i])) starts.push(i);
  for (const s0 of starts) {
    if (seenA[s0]) continue;
    seenA[s0] = 1; arrow[3 * s0] = dir[3 * s0]; arrow[3 * s0 + 1] = dir[3 * s0 + 1]; arrow[3 * s0 + 2] = dir[3 * s0 + 2];
    const queue = [s0];
    for (let q = 0; q < queue.length; q++) {
      const i = queue[q];
      for (const j of nbr[i]) { if (seenA[j] || !Number.isFinite(theta[j])) continue; seenA[j] = 1; choose(j, arrow[3 * i], arrow[3 * i + 1], arrow[3 * i + 2]); queue.push(j); }
    }
  }
  return { theta, dir, arrow, mis, coherence: count ? Math.hypot(cx, cy) / count : 0, count };
}
