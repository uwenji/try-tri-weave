// Mesh repair for the browser import — the steps of carrier.py Carrier.__init__ before the remesh, plus what a
// user's mesh needs before a corner table can be built on it:
//   weld (rel 1e-6 of the bbox diagonal) · degenerate and duplicate faces · T-junction cracks (repair_cracks) ·
//   consistent orientation · non-manifold edges cut, non-manifold vertices split · tiny loose parts dropped ·
//   planar caps (strip_planar_caps) · small holes (fill_small_holes, perimeter < 3 h).
// Corner convention (same as remesh.js): corner c of face f = ⌊c/3⌋ sits at vertex F[c]; the edge opposite c runs
// F[next(c)] → F[prev(c)]; O[c] is the corner across that edge (−1 on the boundary).

export const nxt = (c) => (c % 3 === 2 ? c - 2 : c + 1);
export const prv = (c) => (c % 3 === 0 ? c + 2 : c - 1);

export function bounds(P) {
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < P.length; i += 3) for (let k = 0; k < 3; k++) { const v = P[i + k]; if (v < lo[k]) lo[k] = v; if (v > hi[k]) hi[k] = v; }
  return { lo, hi, diag: Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) };
}

// ------------------------------------------------------------------ weld (hash on quantised coordinates)
/** P: positions, F: triangle indices or null (soup). Points closer than relTol·diag (same rounding cell, as in
 *  carrier.py weld) become one vertex; faces with a repeated vertex are dropped. */
export function weld(P, F, relTol = 1e-6) {
  const np = P.length / 3, { diag } = bounds(P), q = relTol * (diag || 1);
  let size = 1; while (size < 2 * np) size <<= 1; const mask = size - 1;
  const tab = new Int32Array(size).fill(-1), kx = new Int32Array(np), ky = new Int32Array(np), kz = new Int32Array(np);
  const remap = new Int32Array(np), out = new Float64Array(3 * np); let nv = 0;
  for (let i = 0; i < np; i++) {
    const x = Math.round(P[3 * i] / q), y = Math.round(P[3 * i + 1] / q), z = Math.round(P[3 * i + 2] / q);
    let hsh = (Math.imul(x, 73856093) ^ Math.imul(y, 19349663) ^ Math.imul(z, 83492791)) & mask;
    for (;;) {
      const id = tab[hsh];
      if (id < 0) { tab[hsh] = nv; kx[nv] = x; ky[nv] = y; kz[nv] = z; out[3 * nv] = P[3 * i]; out[3 * nv + 1] = P[3 * i + 1]; out[3 * nv + 2] = P[3 * i + 2]; remap[i] = nv++; break; }
      if (kx[id] === x && ky[id] === y && kz[id] === z) { remap[i] = id; break; }
      hsh = (hsh + 1) & mask;
    }
  }
  const m = F ? F.length / 3 : np / 3, G = new Int32Array(3 * m); let nf = 0;
  for (let f = 0; f < m; f++) {
    const a = remap[F ? F[3 * f] : 3 * f], b = remap[F ? F[3 * f + 1] : 3 * f + 1], c = remap[F ? F[3 * f + 2] : 3 * f + 2];
    if (a === b || b === c || a === c) continue;
    G[3 * nf] = a; G[3 * nf + 1] = b; G[3 * nf + 2] = c; nf++;
  }
  return { V: out.slice(0, 3 * nv), F: G.slice(0, 3 * nf) };
}

// ------------------------------------------------------------------ face / edge bookkeeping
/** drop faces that repeat another face's vertex set (either orientation) */
export function dedupeFaces(F, n) {
  const m = F.length / 3, cnt = new Int32Array(n + 1);
  const lo = new Int32Array(m), mid = new Int32Array(m), hi = new Int32Array(m);
  for (let f = 0; f < m; f++) {
    let a = F[3 * f], b = F[3 * f + 1], c = F[3 * f + 2], t;
    if (a > b) { t = a; a = b; b = t; } if (b > c) { t = b; b = c; c = t; } if (a > b) { t = a; a = b; b = t; }
    lo[f] = a; mid[f] = b; hi[f] = c; cnt[a + 1]++;
  }
  for (let i = 0; i < n; i++) cnt[i + 1] += cnt[i];
  const pos = cnt.slice(0, n), ord = new Int32Array(m); for (let f = 0; f < m; f++) ord[pos[lo[f]]++] = f;
  const keep = new Uint8Array(m).fill(1); let removed = 0;
  for (let v = 0; v < n; v++) {
    const s = cnt[v], e = cnt[v + 1];
    for (let i = s; i < e; i++) { const f = ord[i]; if (!keep[f]) continue;
      for (let j = i + 1; j < e; j++) { const g = ord[j]; if (keep[g] && mid[g] === mid[f] && hi[g] === hi[f]) { keep[g] = 0; removed++; } } }
  }
  if (!removed) return { F, removed };
  const G = new Int32Array(3 * (m - removed)); let k = 0;
  for (let f = 0; f < m; f++) if (keep[f]) { G[k++] = F[3 * f]; G[k++] = F[3 * f + 1]; G[k++] = F[3 * f + 2]; }
  return { F: G, removed };
}

/** Group the corners by the undirected edge opposite them. Returns {O, same, nonManifold, boundary}:
 *  strict = true links only pairs of opposite direction (a consistently oriented manifold); strict = false links any
 *  pair and flags pairs of the same direction in `same` (for the orientation pass). Edges with more than two faces
 *  are paired greedily (opposite directions first) and the rest is cut. */
export function edgeMates(F, n, strict = true) {
  const nc = F.length, O = new Int32Array(nc).fill(-1), same = strict ? null : new Uint8Array(nc);
  const cnt = new Int32Array(n + 1), A = new Int32Array(nc), Bv = new Int32Array(nc);
  for (let c = 0; c < nc; c++) { const a = F[nxt(c)], b = F[prv(c)]; A[c] = a; Bv[c] = b; cnt[Math.min(a, b) + 1]++; }
  for (let i = 0; i < n; i++) cnt[i + 1] += cnt[i];
  const pos = cnt.slice(0, n), ord = new Int32Array(nc);
  for (let c = 0; c < nc; c++) ord[pos[Math.min(A[c], Bv[c])]++] = c;
  let nonManifold = 0, boundary = 0; const grp = [];
  for (let v = 0; v < n; v++) {
    const s = cnt[v], e = cnt[v + 1]; if (e - s === 0) continue;
    // sort this bucket's corners by the other end (buckets are small: insertion sort)
    const key = (c) => Math.max(A[c], Bv[c]);
    if (e - s > 48) { const sub = Array.from(ord.subarray(s, e)).sort((x, y) => key(x) - key(y)); ord.set(sub, s); }
    else for (let i = s + 1; i < e; i++) { const c = ord[i], kc = key(c); let j = i - 1; while (j >= s && key(ord[j]) > kc) { ord[j + 1] = ord[j]; j--; } ord[j + 1] = c; }
    let i = s;
    while (i < e) {
      let j = i + 1; const k0 = key(ord[i]); while (j < e && key(ord[j]) === k0) j++;
      const g = j - i;
      if (g === 1) boundary++;
      else if (g === 2) {
        const c1 = ord[i], c2 = ord[i + 1], opp = A[c1] === Bv[c2];
        if (opp || !strict) { O[c1] = c2; O[c2] = c1; if (!opp && same) { same[c1] = 1; same[c2] = 1; } }
        else boundary += 2;
      } else {
        nonManifold++; grp.length = 0; for (let t = i; t < j; t++) grp.push(ord[t]);
        const used = new Uint8Array(g);
        for (let x = 0; x < g; x++) { if (used[x]) continue;
          for (let y = x + 1; y < g; y++) if (!used[y] && A[grp[x]] === Bv[grp[y]]) { used[x] = used[y] = 1; O[grp[x]] = grp[y]; O[grp[y]] = grp[x]; break; } }
      }
      i = j;
    }
  }
  return { O, same, nonManifold, boundary };
}

/** consistent orientation over manifold edges (BFS per connected part); conflicts (Möbius-like) are left as cuts.
 *  With V, a non-orientable mesh gets a short cut: the BFS leaves the conflicts where its fronts met (a zigzag through
 *  a row of faces, a staircase), and single faces are then turned over while that shortens the cut (each turn moves the
 *  cut across the face; it settles as a locally shortest edge path). Returns {F, flipped, conflicts, cut: [a, b, …]
 *  the cut edges, cutLength0 / cutLength before / after the shortening}. */
export function orientFaces(F, n, V = null) {
  const m = F.length / 3, { O, same } = edgeMates(F, n, false);
  const flip = new Uint8Array(m), seen = new Uint8Array(m), stack = new Int32Array(m); let conflicts = 0, flipped = 0;
  for (let r = 0; r < m; r++) {
    if (seen[r]) continue; seen[r] = 1; let sp = 0; stack[sp++] = r;
    while (sp) {
      const f = stack[--sp];
      for (let s = 0; s < 3; s++) {
        const c = 3 * f + s, o = O[c]; if (o < 0) continue;
        const g = (o / 3) | 0, want = flip[f] ^ same[c];
        if (!seen[g]) { seen[g] = 1; flip[g] = want; stack[sp++] = g; }
        else if (flip[g] !== want) conflicts++;
      }
    }
  }
  const bad = (c) => { const o = O[c]; return o >= 0 && (flip[(c / 3) | 0] ^ same[c]) !== flip[(o / 3) | 0]; };
  const len = (c) => { const a = 3 * F[nxt(c)], b = 3 * F[prv(c)]; return Math.hypot(V[a] - V[b], V[a + 1] - V[b + 1], V[a + 2] - V[b + 2]); };
  let cutLength0 = 0, cutLength = 0;
  if (conflicts && V) {
    for (let c = 0; c < 3 * m; c++) if (O[c] > c && bad(c)) cutLength0 += len(c);
    const inQ = new Uint8Array(m), queue = [];
    for (let c = 0; c < 3 * m; c++) if (O[c] >= 0 && bad(c)) { const f = (c / 3) | 0; if (!inQ[f]) { inQ[f] = 1; queue.push(f); } }
    const tol = 1e-9 * (cutLength0 || 1); let guard = 0;
    while (queue.length && guard++ < 200 * m) {
      const f = queue.pop(); inQ[f] = 0; let gain = 0;
      for (let s = 0; s < 3; s++) { const c = 3 * f + s; if (O[c] >= 0) gain += bad(c) ? -len(c) : len(c); }
      if (gain >= -tol) continue;
      flip[f] ^= 1;
      for (let s = 0; s < 3; s++) { const o = O[3 * f + s]; if (o < 0) continue; const g = (o / 3) | 0; if (!inQ[g]) { inQ[g] = 1; queue.push(g); } }
    }
  }
  const cut = [];
  for (let c = 0; c < 3 * m; c++) if (O[c] > c && bad(c)) { cut.push(F[nxt(c)], F[prv(c)]); if (V) cutLength += len(c); }
  // the global sign (outward / convex side) is applied to the final carrier by carrier.js orientSide
  const G = Int32Array.from(F);
  for (let f = 0; f < m; f++) if (flip[f]) { const t = G[3 * f + 1]; G[3 * f + 1] = G[3 * f + 2]; G[3 * f + 2] = t; flipped++; }
  return { F: G, flipped, conflicts: V ? cut.length / 2 : conflicts >> 1, cut, cutLength0, cutLength };
}

/** The cut edges of orientFaces as polylines: [{P: Float64Array (positions along the path), closed}], or null when the
 *  cut branches (a vertex on three or more cut edges) — then it is not re-glued. */
export function cutPaths(V, cut) {
  const adj = new Map();
  for (let k = 0; k < cut.length; k += 2) for (const [x, y] of [[cut[k], cut[k + 1]], [cut[k + 1], cut[k]]]) { let l = adj.get(x); if (!l) adj.set(x, l = []); l.push(y); }
  for (const l of adj.values()) if (l.length > 2) return null;
  const used = new Set(), out = [], key = (a, b) => (a < b ? a + ',' + b : b + ',' + a);
  const walk = (v0) => {
    const pts = [v0]; let prev = -1, cur = v0;
    for (;;) { const nx = adj.get(cur).find(w => w !== prev && !used.has(key(cur, w))); if (nx === undefined) break; used.add(key(cur, nx)); pts.push(nx); prev = cur; cur = nx; if (nx === v0) break; }
    return pts;
  };
  for (const [v, l] of adj) if (l.length === 1 && !used.has(key(v, l[0]))) out.push({ ids: walk(v), closed: false });
  for (const [v, l] of adj) if (l.some(w => !used.has(key(v, w)))) out.push({ ids: walk(v), closed: true });
  return out.map(p => ({ closed: p.closed, P: Float64Array.from(p.ids.flatMap(i => [V[3 * i], V[3 * i + 1], V[3 * i + 2]])) }));
}

/**
 * Re-glue the orientation cut of a non-orientable surface after the remesh (orientFaces cut it so that the corner
 * table and the remesher see an orientable disc; the remesher samples the two sides of the cut independently).
 * paths: cutPaths (open paths only); V, F: the remeshed mesh (consistently oriented, the cut two boundary runs per path);
 * h: the target edge. The boundary vertices on a path (they were projected onto it) get its arc length s; the two runs
 * are merged into one sampling — samples of the two sides closer than h/2 (mutual nearest) become one vertex at their
 * mean s, every other sample is inserted into the other side's boundary edge (a fan in that edge's face) — and the
 * two sides are welded. Returns {V, F, seam: vertices on the glued seam, moved: per vertex 1 if on the seam} or
 * {error} (then the caller keeps the cut).
 */
export function zipCut(V, F, paths, h) {
  const n = V.length / 3, m = F.length / 3, { O } = edgeMates(F, n, true);
  const nextB = new Int32Array(n).fill(-1), isB = new Uint8Array(n), faceOf = new Map();
  for (let c = 0; c < 3 * m; c++) if (O[c] < 0) { const a = F[nxt(c)], b = F[prv(c)]; nextB[a] = b; isB[a] = isB[b] = 1; faceOf.set(a * 4294967296 + b, (c / 3) | 0); }
  // closest point of every boundary vertex on every path: (path, s, d)
  const segs = paths.map(p => { const k = p.P.length / 3 - 1, S = new Float64Array(k + 1); for (let i = 0; i < k; i++) S[i + 1] = S[i] + Math.hypot(p.P[3 * i + 3] - p.P[3 * i], p.P[3 * i + 4] - p.P[3 * i + 1], p.P[3 * i + 5] - p.P[3 * i + 2]); return { k, S, P: p.P }; });
  const onPath = new Int32Array(n).fill(-1), sOf = new Float64Array(n), tol = 1e-6 * h;
  for (let v = 0; v < n; v++) {
    if (!isB[v]) continue;
    const x = V[3 * v], y = V[3 * v + 1], z = V[3 * v + 2]; let bd = Infinity, bp = -1, bs = 0;
    segs.forEach((sg, pi) => {
      for (let i = 0; i < sg.k; i++) {
        const P = sg.P, ax = P[3 * i], ay = P[3 * i + 1], az = P[3 * i + 2], dx = P[3 * i + 3] - ax, dy = P[3 * i + 4] - ay, dz = P[3 * i + 5] - az, L2 = dx * dx + dy * dy + dz * dz;
        let t = L2 > 0 ? ((x - ax) * dx + (y - ay) * dy + (z - az) * dz) / L2 : 0; t = t < 0 ? 0 : t > 1 ? 1 : t;
        const d = Math.hypot(ax + t * dx - x, ay + t * dy - y, az + t * dz - z);
        if (d < bd) { bd = d; bp = pi; bs = sg.S[i] + t * (sg.S[i + 1] - sg.S[i]); }
      }
    });
    if (bd < tol) { onPath[v] = bp; sOf[v] = bs; }
  }
  // runs of consecutive on-path boundary vertices, per boundary loop
  const runs = [], seenL = new Uint8Array(n);
  for (let v0 = 0; v0 < n; v0++) {
    if (!isB[v0] || seenL[v0]) continue;
    const loop = []; let v = v0, guard = 0; while (!seenL[v] && guard++ < n) { seenL[v] = 1; loop.push(v); v = nextB[v]; if (v < 0) break; }
    const L = loop.length, st = loop.findIndex(u => onPath[u] < 0);
    if (st < 0) return { error: 'a boundary loop lies entirely on the cut (a closed cut is not re-glued)' };
    let cur = null;
    for (let i = 1; i <= L; i++) {
      const u = loop[(st + i) % L];
      if (onPath[u] < 0) { cur = null; continue; }
      if (!cur || onPath[u] !== cur.path) { cur = { path: onPath[u], v: [] }; runs.push(cur); }
      cur.v.push(u);
    }
  }
  const remap = new Int32Array(n); for (let i = 0; i < n; i++) remap[i] = i;
  const Vout = Float64Array.from(V), addF = [], dropF = new Uint8Array(m), seamV = new Uint8Array(n);
  for (let pi = 0; pi < paths.length; pi++) {
    const sg = segs[pi], Stot = sg.S[sg.k], mine = runs.filter(r => r.path === pi && r.v.length >= 2);
    if (mine.length !== 2) return { error: `cut path ${pi}: ${mine.length} boundary runs (2 expected)` };
    const [RA, RB] = mine.map(r => { const s = r.v.map(u => sOf[u]), inc = s[s.length - 1] >= s[0]; return { v: r.v, s, inc }; });
    for (const R of [RA, RB]) {
      const lo = Math.min(R.s[0], R.s[R.s.length - 1]), hi = Math.max(R.s[0], R.s[R.s.length - 1]);
      if (lo > 1e-6 * Stot + tol || hi < Stot - 1e-6 * Stot - tol) return { error: `cut path ${pi}: a side does not span the cut` };
      for (let i = 1; i < R.s.length; i++) if ((R.s[i] - R.s[i - 1]) * (R.inc ? 1 : -1) < -tol) return { error: `cut path ${pi}: a side is not monotone along the cut` };
    }
    const sortedOf = (R) => R.v.map((u, i) => ({ u, s: R.s[i] })).sort((p, q) => p.s - q.s);
    const A = sortedOf(RA), B = sortedOf(RB);
    // mutual nearest pairs closer than h/2 → one vertex
    const nearest = (X, Y) => { const out = new Int32Array(X.length); let j = 0; for (let i = 0; i < X.length; i++) { while (j + 1 < Y.length && Math.abs(Y[j + 1].s - X[i].s) <= Math.abs(Y[j].s - X[i].s)) j++; out[i] = j; } return out; };
    const nAB = nearest(A, B), nBA = nearest(B, A), pairA = new Int32Array(A.length).fill(-1), pairB = new Int32Array(B.length).fill(-1);
    for (let i = 0; i < A.length; i++) { const j = nAB[i]; if (nBA[j] === i && Math.abs(A[i].s - B[j].s) < 0.5 * h) { pairA[i] = j; pairB[j] = i; } }
    if (pairA[0] < 0 || pairA[A.length - 1] < 0) return { error: `cut path ${pi}: the ends of the two sides do not meet` };
    const at = (s) => { let i = 0; while (i < sg.k - 1 && sg.S[i + 1] < s) i++; const t = (s - sg.S[i]) / Math.max(sg.S[i + 1] - sg.S[i], 1e-300), P = sg.P; return [0, 1, 2].map(d => P[3 * i + d] + t * (P[3 * i + 3 + d] - P[3 * i + d])); };
    for (let i = 0; i < A.length; i++) if (pairA[i] >= 0) {
      const a = A[i].u, b = B[pairA[i]].u, s = 0.5 * (A[i].s + B[pairA[i]].s), p = at(s);
      remap[b] = a; Vout.set(p, 3 * a); A[i].s = B[pairA[i]].s = s; seamV[a] = 1;
    } else seamV[A[i].u] = 1;
    for (let j = 0; j < B.length; j++) if (pairB[j] < 0) seamV[B[j].u] = 1;
    // every unpaired sample into the other side's boundary edge: a fan in that edge's face
    const insert = (X, Y) => {                                              // samples of X without a partner → edges of Y
      for (let j = 0; j + 1 < Y.length; j++) {
        const lo = Y[j].s, hi = Y[j + 1].s, ins = X.filter((x, i) => (X === A ? pairA[i] : pairB[i]) < 0 && x.s > lo && x.s < hi);
        if (!ins.length) continue;
        const p = Y[j].u, q = Y[j + 1].u, k1 = p * 4294967296 + q, k2 = q * 4294967296 + p;
        const f = faceOf.has(k1) ? faceOf.get(k1) : faceOf.has(k2) ? faceOf.get(k2) : -1;
        if (f < 0) return `no boundary face between two consecutive samples`;
        const a = faceOf.has(k1) ? p : q, b = a === p ? q : p;             // the face runs a → b along the edge
        const r = [F[3 * f], F[3 * f + 1], F[3 * f + 2]].find(w => w !== a && w !== b);
        const seq = ins.map(x => x.u); if (a !== p) seq.reverse();         // from a to b
        const chain = [a, ...seq, b];
        for (let t = 0; t + 1 < chain.length; t++) addF.push(chain[t], chain[t + 1], r);
        dropF[f] = 1;
      }
      return null;
    };
    const e1 = insert(A, B) || insert(B, A); if (e1) return { error: e1 };
  }
  // assemble: the kept faces and the fans, welded; compact the vertices
  const G = []; for (let f = 0; f < m; f++) if (!dropF[f]) G.push(F[3 * f], F[3 * f + 1], F[3 * f + 2]);
  G.push(...addF);
  const used = new Uint8Array(n); for (let k = 0; k < G.length; k++) { G[k] = remap[G[k]]; used[G[k]] = 1; }
  const idx = new Int32Array(n).fill(-1); let nv = 0; for (let i = 0; i < n; i++) if (used[i]) idx[i] = nv++;
  const V2 = new Float64Array(3 * nv), seam = new Uint8Array(nv); for (let i = 0; i < n; i++) if (idx[i] >= 0) { V2.set(Vout.subarray(3 * i, 3 * i + 3), 3 * idx[i]); seam[idx[i]] = seamV[i]; }
  const F2 = new Int32Array(G.length); let bad = 0;
  for (let f = 0; f < G.length / 3; f++) { for (let k = 0; k < 3; k++) F2[3 * f + k] = idx[G[3 * f + k]]; if (F2[3 * f] === F2[3 * f + 1] || F2[3 * f + 1] === F2[3 * f + 2] || F2[3 * f] === F2[3 * f + 2]) bad++; }
  if (bad) return { error: `${bad} degenerate faces after the weld` };
  let nSeam = 0; for (const x of seam) nSeam += x;
  return { V: V2, F: F2, seam, nSeam, idx };
}

/** strict corner table; vertices whose corners form more than one fan (non-manifold vertices, also those left by
 *  cut edges) are split, one copy per fan. Returns {V, F, O, splitVerts, nonManifoldEdges}. */
export function manifoldTable(V, F) {
  let n = V.length / 3; const { O, nonManifold } = edgeMates(F, n, true);
  const nc = F.length, seen = new Uint8Array(nc), owner = new Uint8Array(n); const extra = []; let split = 0;
  const G = Int32Array.from(F); const fan = [];
  for (let c0 = 0; c0 < nc; c0++) {
    if (seen[c0]) continue; const v = G[c0];
    // collect the fan of c0: swing both ways
    fan.length = 0; let c = c0, guard = 0;
    do { fan.push(c); seen[c] = 1; const o = O[prv(c)]; if (o < 0) { c = -1; break; } c = prv(o); if (++guard > nc) break; } while (c !== c0);
    if (c === -1) { let o = O[nxt(c0)]; while (o >= 0) { const d = nxt(o); if (seen[d]) break; fan.push(d); seen[d] = 1; o = O[nxt(d)]; } }
    if (!owner[v]) { owner[v] = 1; continue; }
    const nv = n + extra.length / 3; extra.push(V[3 * v], V[3 * v + 1], V[3 * v + 2]); split++;
    for (const d of fan) G[d] = nv;
  }
  if (split) { const W = new Float64Array(V.length + extra.length); W.set(V); W.set(extra, V.length); V = W; }
  return { V, F: G, O, splitVerts: split, nonManifoldEdges: nonManifold };
}

export function faceAreas(V, F) {
  const m = F.length / 3, A = new Float64Array(m);
  for (let f = 0; f < m; f++) {
    const a = 3 * F[3 * f], b = 3 * F[3 * f + 1], c = 3 * F[3 * f + 2];
    const ux = V[b] - V[a], uy = V[b + 1] - V[a + 1], uz = V[b + 2] - V[a + 2], wx = V[c] - V[a], wy = V[c + 1] - V[a + 1], wz = V[c + 2] - V[a + 2];
    A[f] = 0.5 * Math.hypot(uy * wz - uz * wy, uz * wx - ux * wz, ux * wy - uy * wx);
  }
  return A;
}

export function faceNormals(V, F) {
  const m = F.length / 3, N = new Float64Array(3 * m);
  for (let f = 0; f < m; f++) {
    const a = 3 * F[3 * f], b = 3 * F[3 * f + 1], c = 3 * F[3 * f + 2];
    const ux = V[b] - V[a], uy = V[b + 1] - V[a + 1], uz = V[b + 2] - V[a + 2], wx = V[c] - V[a], wy = V[c + 1] - V[a + 1], wz = V[c + 2] - V[a + 2];
    const nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx, l = Math.hypot(nx, ny, nz) || 1;
    N[3 * f] = nx / l; N[3 * f + 1] = ny / l; N[3 * f + 2] = nz / l;
  }
  return N;
}

function uf(m) {
  const p = new Int32Array(m); for (let i = 0; i < m; i++) p[i] = i;
  const find = (x) => { while (p[x] !== x) { p[x] = p[p[x]]; x = p[x]; } return x; };
  return { p, find, union: (a, b) => { a = find(a); b = find(b); if (a !== b) p[a < b ? b : a] = a < b ? a : b; } };
}

/** keep the faces with keep[f]; unused vertices are dropped */
export function subsetFaces(V, F, keep) {
  const n = V.length / 3, m = F.length / 3, map = new Int32Array(n).fill(-1); let nv = 0, nf = 0;
  for (let f = 0; f < m; f++) if (keep[f]) { nf++; for (let s = 0; s < 3; s++) { const v = F[3 * f + s]; if (map[v] < 0) map[v] = nv++; } }
  const W = new Float64Array(3 * nv), G = new Int32Array(3 * nf);
  for (let v = 0; v < n; v++) if (map[v] >= 0) { W[3 * map[v]] = V[3 * v]; W[3 * map[v] + 1] = V[3 * v + 1]; W[3 * map[v] + 2] = V[3 * v + 2]; }
  let k = 0; for (let f = 0; f < m; f++) if (keep[f]) for (let s = 0; s < 3; s++) G[k++] = map[F[3 * f + s]];
  return { V: W, F: G };
}

/** connected parts by the corner table; drops parts with area < minArea. Returns {V, F, dropped, parts} */
export function dropSmallParts(V, F, O, minArea) {
  const m = F.length / 3, U = uf(m);
  for (let c = 0; c < 3 * m; c++) if (O[c] >= 0) U.union((c / 3) | 0, (O[c] / 3) | 0);
  const A = faceAreas(V, F), area = new Float64Array(m);
  for (let f = 0; f < m; f++) area[U.find(f)] += A[f];
  let parts = 0, dropped = 0; const keep = new Uint8Array(m);
  for (let f = 0; f < m; f++) if (U.find(f) === f) parts++;
  const big = new Set(); let best = -1, bestA = -1;
  for (let f = 0; f < m; f++) if (U.find(f) === f) { if (area[f] >= minArea) big.add(f); if (area[f] > bestA) { bestA = area[f]; best = f; } }
  if (!big.size) big.add(best);
  for (let f = 0; f < m; f++) { const r = U.find(f); if (big.has(r)) keep[f] = 1; }
  dropped = parts - big.size;
  if (!dropped) return { V, F, dropped, parts };
  return { ...subsetFaces(V, F, keep), dropped, parts };
}

// ------------------------------------------------------------------ repair_cracks (T-junctions)
/** A boundary vertex lying on a boundary edge it does not belong to (closer than relTol × the edge length,
 *  strictly between its ends) splits the face owning that edge, so both sides of the seam share edges again.
 *  One split per face per round, repeated until nothing changes (carrier.py repair_cracks). */
export function repairCracks(V, F, relTol = 0.02, maxRounds = 64) {
  const n = V.length / 3; let splits = 0;
  for (let round = 0; round < maxRounds; round++) {
    const m = F.length / 3, { O } = edgeMates(F, n, false);
    const bc = []; for (let c = 0; c < 3 * m; c++) if (O[c] < 0) bc.push(c);
    if (!bc.length) break;
    const isB = new Uint8Array(n); let lsum = 0;
    for (const c of bc) { const a = F[nxt(c)], b = F[prv(c)]; isB[a] = isB[b] = 1; lsum += Math.hypot(V[3 * a] - V[3 * b], V[3 * a + 1] - V[3 * b + 1], V[3 * a + 2] - V[3 * b + 2]); }
    const s = Math.max(lsum / bc.length, 1e-12), grid = new Map();
    const key = (x, y, z) => x * 73856093 + y * 19349663 + z * 83492791;     // hashed (collisions only cost time)
    const cell = (v) => [Math.floor(V[3 * v] / s), Math.floor(V[3 * v + 1] / s), Math.floor(V[3 * v + 2] / s)];
    for (let v = 0; v < n; v++) if (isB[v]) { const [x, y, z] = cell(v), k = key(x, y, z); let l = grid.get(k); if (!l) grid.set(k, l = []); l.push(v); }
    const stamp = new Int32Array(n).fill(-1), faceSplit = new Map();
    for (let e = 0; e < bc.length; e++) {
      const c = bc[e], f = (c / 3) | 0; if (faceSplit.has(f)) continue;
      const a = F[nxt(c)], b = F[prv(c)];
      const dx = V[3 * b] - V[3 * a], dy = V[3 * b + 1] - V[3 * a + 1], dz = V[3 * b + 2] - V[3 * a + 2], L2 = dx * dx + dy * dy + dz * dz, L = Math.sqrt(L2);
      if (L <= 0) continue;
      const tol = relTol * L, R = Math.ceil(tol / s) + 1, steps = Math.ceil(L / s) + 1;
      let best = null;
      for (let t = 0; t <= steps; t++) {
        const px = V[3 * a] + dx * t / steps, py = V[3 * a + 1] + dy * t / steps, pz = V[3 * a + 2] + dz * t / steps;
        const cx = Math.floor(px / s), cy = Math.floor(py / s), cz = Math.floor(pz / s);
        for (let i = -R; i <= R; i++) for (let j = -R; j <= R; j++) for (let k = -R; k <= R; k++) {
          const l = grid.get(key(cx + i, cy + j, cz + k)); if (!l) continue;
          for (const v of l) {
            if (stamp[v] === e || v === a || v === b) continue; stamp[v] = e;
            const qx = V[3 * v] - V[3 * a], qy = V[3 * v + 1] - V[3 * a + 1], qz = V[3 * v + 2] - V[3 * a + 2];
            const tt = (qx * dx + qy * dy + qz * dz) / L2; if (!(tt > 1e-3 && tt < 1 - 1e-3)) continue;
            const d = Math.hypot(qx - tt * dx, qy - tt * dy, qz - tt * dz);
            if (d < tol && (!best || tt < best[0])) best = [tt, v];
          }
        }
      }
      if (best) faceSplit.set(f, [c, best[1]]);
    }
    if (!faceSplit.size) break;
    const G = new Int32Array(F.length + 3 * faceSplit.size); G.set(F); let k = F.length;
    for (const [f, [c, v]] of faceSplit) {
      const r = F[c], p = F[nxt(c)], q = F[prv(c)];                           // face (r, p, q), v on edge p→q
      G[3 * f] = r; G[3 * f + 1] = p; G[3 * f + 2] = v; G[k++] = r; G[k++] = v; G[k++] = q;
    }
    F = G; splits += faceSplit.size;
  }
  return { F, splits };
}

// ------------------------------------------------------------------ strip_planar_caps
/** Remove planar face groups whose whole rim is a crease (the flat base of a closed Rhino cone / tube); groups
 *  smaller than minFrac of the area, or touching an open edge, are surface, not caps. */
export function stripPlanarCaps(V, F, O, creaseDeg = 30, coplanarDeg = 1, minFrac = 0.01) {
  const m = F.length / 3, N = faceNormals(V, F), A = faceAreas(V, F), U = uf(m);
  const cCo = Math.cos(coplanarDeg * Math.PI / 180), cCr = Math.cos(creaseDeg * Math.PI / 180);
  const dotc = (c) => { const f = (c / 3) | 0, g = (O[c] / 3) | 0; return N[3 * f] * N[3 * g] + N[3 * f + 1] * N[3 * g + 1] + N[3 * f + 2] * N[3 * g + 2]; };
  for (let c = 0; c < 3 * m; c++) if (O[c] > c && dotc(c) > cCo) U.union((c / 3) | 0, (O[c] / 3) | 0);
  const n = V.length / 3, bv = new Uint8Array(n); for (let c = 0; c < 3 * m; c++) if (O[c] < 0) { bv[F[nxt(c)]] = 1; bv[F[prv(c)]] = 1; }
  let total = 0; for (let f = 0; f < m; f++) total += A[f];
  const size = new Int32Array(m), area = new Float64Array(m), touches = new Uint8Array(m), rimOk = new Uint8Array(m).fill(1), rimAny = new Uint8Array(m);
  for (let f = 0; f < m; f++) { const r = U.find(f); size[r]++; area[r] += A[f]; if (bv[F[3 * f]] || bv[F[3 * f + 1]] || bv[F[3 * f + 2]]) touches[r] = 1; }
  for (let c = 0; c < 3 * m; c++) {
    if (O[c] < 0) continue; const rf = U.find((c / 3) | 0), rg = U.find((O[c] / 3) | 0); if (rf === rg) continue;
    rimAny[rf] = 1; if (!(dotc(c) < cCr)) rimOk[rf] = 0;
  }
  const keep = new Uint8Array(m).fill(1); let removed = 0, removedArea = 0;
  for (let f = 0; f < m; f++) {
    const r = U.find(f);
    if (size[r] < 2 || size[r] >= m - 2 || area[r] < minFrac * total || touches[r] || !rimAny[r] || !rimOk[r]) continue;
    keep[f] = 0; removed++; removedArea += A[f];
  }
  // caps are the small flat ends of a modelled cone or tube: a box (all of it flat) is not stripped
  if (!removed || removed === m || removedArea > 0.6 * total) return { V, F, removed: 0 };
  return { ...subsetFaces(V, F, keep), removed };
}

// ------------------------------------------------------------------ fill_small_holes
/** boundary loops (corner table walk) with their perimeter */
export function boundaryLoops(V, F, O) {
  const nc = F.length, n = V.length / 3, out = new Int32Array(n).fill(-1);
  for (let c = 0; c < nc; c++) if (O[c] < 0 && F[c] >= 0) out[F[nxt(c)]] = c;      // boundary half-edge leaving vertex F[next c]
  const used = new Uint8Array(nc), loops = [];
  for (let c0 = 0; c0 < nc; c0++) {
    if (O[c0] >= 0 || used[c0] || F[c0] < 0) continue;
    const edges = []; let c = c0, per = 0, guard = 0;
    while (c >= 0 && !used[c] && guard++ < nc) {
      used[c] = 1; edges.push(c); const a = F[nxt(c)], b = F[prv(c)];
      per += Math.hypot(V[3 * a] - V[3 * b], V[3 * a + 1] - V[3 * b + 1], V[3 * a + 2] - V[3 * b + 2]);
      c = out[b];
    }
    loops.push({ edges, perimeter: per });
  }
  return loops.sort((x, y) => y.perimeter - x.perimeter);
}

/** close loops shorter than maxPerimeter with a fan from the loop centroid (carrier.py fill_small_holes) */
export function fillSmallHoles(V, F, O, maxPerimeter) {
  const loops = boundaryLoops(V, F, O).filter(l => l.perimeter < maxPerimeter && l.edges.length >= 3);
  if (!loops.length) return { V, F, filled: 0 };
  const n0 = V.length / 3, W = new Float64Array(V.length + 3 * loops.length); W.set(V);
  const add = [];
  loops.forEach((l, i) => {
    const cv = n0 + i; let x = 0, y = 0, z = 0;
    for (const c of l.edges) { const a = F[nxt(c)]; x += V[3 * a]; y += V[3 * a + 1]; z += V[3 * a + 2]; }
    W[3 * cv] = x / l.edges.length; W[3 * cv + 1] = y / l.edges.length; W[3 * cv + 2] = z / l.edges.length;
    for (const c of l.edges) add.push(F[prv(c)], F[nxt(c)], cv);                      // opposite to the face's own edge
  });
  const G = new Int32Array(F.length + add.length); G.set(F); G.set(add, F.length);
  return { V: W, F: G, filled: loops.length };
}
