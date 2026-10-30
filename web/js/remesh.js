// Isotropic remeshing in the browser (Botsch & Kobbelt 2004, the scheme gpytoolbox remesh_botsch runs in
// carrier.py): per iteration split edges > 4/3 h, collapse edges < 4/5 h (if no new edge exceeds 4/3 h), flip
// towards valence 6 (4 on the boundary), tangential relaxation, projection back onto the input surface.
//
// Features are kept as in carrier.py: the boundary and creases (dihedral > creaseDeg) are feature chains; their
// vertices only slide along the chain and are projected onto the input chain; chain ends, junctions, chain corners
// (turning > cornerDeg) and sharp apices (|K| > sharpDeg, the cone tips) are fixed. A fixed apex keeps the whole
// angle defect at one regular-valence vertex (carrier.py snaps the nearest vertex back onto it; same outcome).
//
// Data structure: corner table on typed arrays (corner c of face ⌊c/3⌋ at vertex C[c]; the edge opposite c runs
// C[next c] → C[prev c]; O[c] the corner across it, −1 on the boundary; FL[c] the feature chain of that edge or −1).
// Closest points: AABB tree over the input triangles (with a normal-compatibility filter so a vertex of a thin part
// does not jump to the other side) and over the input feature segments (filtered by chain).
import { nxt, prv, edgeMates, faceNormals } from './meshrepair.js';

// ------------------------------------------------------------------ AABB tree
function buildBVH(blo, bhi, n, leafSize = 6) {
  const idx = new Int32Array(n), cen = new Float64Array(3 * n);
  for (let i = 0; i < n; i++) { idx[i] = i; for (let k = 0; k < 3; k++) cen[3 * i + k] = 0.5 * (blo[3 * i + k] + bhi[3 * i + k]); }
  let cap = Math.max(16, 4 * Math.ceil(n / leafSize));
  let nlo = new Float64Array(3 * cap), nhi = new Float64Array(3 * cap), ns = new Int32Array(cap), nc = new Int32Array(cap), nn = 1;
  const grow = () => {
    const c2 = cap * 2, a = new Float64Array(3 * c2), b = new Float64Array(3 * c2), s = new Int32Array(c2), c = new Int32Array(c2);
    a.set(nlo); b.set(nhi); s.set(ns); c.set(nc); nlo = a; nhi = b; ns = s; nc = c; cap = c2;
  };
  const st = [0, 0, n];
  while (st.length) {
    const e = st.pop(), s = st.pop(), node = st.pop();
    let x0 = Infinity, y0 = Infinity, z0 = Infinity, x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
    let cx0 = Infinity, cy0 = Infinity, cz0 = Infinity, cx1 = -Infinity, cy1 = -Infinity, cz1 = -Infinity;
    for (let i = s; i < e; i++) {
      const p = idx[i], q = 3 * p;
      if (blo[q] < x0) x0 = blo[q]; if (blo[q + 1] < y0) y0 = blo[q + 1]; if (blo[q + 2] < z0) z0 = blo[q + 2];
      if (bhi[q] > x1) x1 = bhi[q]; if (bhi[q + 1] > y1) y1 = bhi[q + 1]; if (bhi[q + 2] > z1) z1 = bhi[q + 2];
      if (cen[q] < cx0) cx0 = cen[q]; if (cen[q + 1] < cy0) cy0 = cen[q + 1]; if (cen[q + 2] < cz0) cz0 = cen[q + 2];
      if (cen[q] > cx1) cx1 = cen[q]; if (cen[q + 1] > cy1) cy1 = cen[q + 1]; if (cen[q + 2] > cz1) cz1 = cen[q + 2];
    }
    nlo[3 * node] = x0; nlo[3 * node + 1] = y0; nlo[3 * node + 2] = z0; nhi[3 * node] = x1; nhi[3 * node + 1] = y1; nhi[3 * node + 2] = z1;
    const ex = cx1 - cx0, ey = cy1 - cy0, ez = cz1 - cz0;
    if (e - s <= leafSize || Math.max(ex, ey, ez) <= 0) { ns[node] = s; nc[node] = e - s; continue; }
    const ax = ex >= ey && ex >= ez ? 0 : ey >= ez ? 1 : 2, mid = ax === 0 ? cx0 + ex / 2 : ax === 1 ? cy0 + ey / 2 : cz0 + ez / 2;
    let i = s, j = e - 1;
    while (i <= j) { if (cen[3 * idx[i] + ax] < mid) i++; else { const t = idx[i]; idx[i] = idx[j]; idx[j] = t; j--; } }
    let m = i; if (m === s || m === e) m = (s + e) >> 1;
    while (nn + 2 > cap) grow();
    const L = nn; nn += 2; ns[node] = L; nc[node] = 0;
    st.push(L + 1, m, e, L, s, m);
  }
  return { idx, nlo, nhi, ns, nc, nn };
}

const STACK = new Int32Array(1024);
/** best primitive by `test(prim) → squared distance or Infinity` among boxes closer than best (init maxD2) */
function bvhQuery(T, px, py, pz, maxD2, test) {
  const { idx, nlo, nhi, ns, nc } = T;
  const bd = (k) => {
    let d = 0, t;
    t = nlo[3 * k] - px; if (t > 0) d += t * t; else { t = px - nhi[3 * k]; if (t > 0) d += t * t; }
    t = nlo[3 * k + 1] - py; if (t > 0) d += t * t; else { t = py - nhi[3 * k + 1]; if (t > 0) d += t * t; }
    t = nlo[3 * k + 2] - pz; if (t > 0) d += t * t; else { t = pz - nhi[3 * k + 2]; if (t > 0) d += t * t; }
    return d;
  };
  let best = maxD2, bestP = -1, sp = 0; STACK[sp++] = 0;
  while (sp) {
    const k = STACK[--sp]; if (bd(k) >= best) continue;
    if (nc[k] > 0) {
      for (let i = ns[k], e = ns[k] + nc[k]; i < e; i++) { const d = test(idx[i], best); if (d < best) { best = d; bestP = idx[i]; } }
    } else {
      const L = ns[k], R = L + 1, dl = bd(L), dr = bd(R);
      if (sp + 2 >= STACK.length) throw new Error('bvh stack');
      if (dl < dr) { if (dr < best) STACK[sp++] = R; if (dl < best) STACK[sp++] = L; }
      else { if (dl < best) STACK[sp++] = L; if (dr < best) STACK[sp++] = R; }
    }
  }
  return { prim: bestP, d2: best };
}

// closest point on triangle abc to p (Ericson, Real-Time Collision Detection 5.1.5); writes Q, returns squared distance
const Q = new Float64Array(3);
function closestOnTri(px, py, pz, V, ia, ib, ic) {
  const ax = V[ia], ay = V[ia + 1], az = V[ia + 2], abx = V[ib] - ax, aby = V[ib + 1] - ay, abz = V[ib + 2] - az;
  const acx = V[ic] - ax, acy = V[ic + 1] - ay, acz = V[ic + 2] - az, apx = px - ax, apy = py - ay, apz = pz - az;
  const d1 = abx * apx + aby * apy + abz * apz, d2 = acx * apx + acy * apy + acz * apz;
  let qx, qy, qz;
  if (d1 <= 0 && d2 <= 0) { qx = ax; qy = ay; qz = az; }
  else {
    const bpx = px - V[ib], bpy = py - V[ib + 1], bpz = pz - V[ib + 2];
    const d3 = abx * bpx + aby * bpy + abz * bpz, d4 = acx * bpx + acy * bpy + acz * bpz;
    if (d3 >= 0 && d4 <= d3) { qx = V[ib]; qy = V[ib + 1]; qz = V[ib + 2]; }
    else {
      const vc = d1 * d4 - d3 * d2;
      if (vc <= 0 && d1 >= 0 && d3 <= 0) { const v = d1 / Math.max(d1 - d3, 1e-300); qx = ax + v * abx; qy = ay + v * aby; qz = az + v * abz; }
      else {
        const cpx = px - V[ic], cpy = py - V[ic + 1], cpz = pz - V[ic + 2];
        const d5 = abx * cpx + aby * cpy + abz * cpz, d6 = acx * cpx + acy * cpy + acz * cpz;
        if (d6 >= 0 && d5 <= d6) { qx = V[ic]; qy = V[ic + 1]; qz = V[ic + 2]; }
        else {
          const vb = d5 * d2 - d1 * d6;
          if (vb <= 0 && d2 >= 0 && d6 <= 0) { const w = d2 / Math.max(d2 - d6, 1e-300); qx = ax + w * acx; qy = ay + w * acy; qz = az + w * acz; }
          else {
            const va = d3 * d6 - d5 * d4;
            if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) {
              const w = (d4 - d3) / Math.max((d4 - d3) + (d5 - d6), 1e-300);
              qx = V[ib] + w * (V[ic] - V[ib]); qy = V[ib + 1] + w * (V[ic + 1] - V[ib + 1]); qz = V[ib + 2] + w * (V[ic + 2] - V[ib + 2]);
            } else {
              const den = va + vb + vc;
              if (!(den > 1e-300)) { qx = ax; qy = ay; qz = az; }
              else { const v = vb / den, w = vc / den; qx = ax + abx * v + acx * w; qy = ay + aby * v + acy * w; qz = az + abz * v + acz * w; }
            }
          }
        }
      }
    }
  }
  Q[0] = qx; Q[1] = qy; Q[2] = qz;
  return (px - qx) ** 2 + (py - qy) ** 2 + (pz - qz) ** 2;
}

function closestOnSeg(px, py, pz, V, ia, ib) {
  const dx = V[ib] - V[ia], dy = V[ib + 1] - V[ia + 1], dz = V[ib + 2] - V[ia + 2], L2 = dx * dx + dy * dy + dz * dz;
  let t = L2 > 0 ? ((px - V[ia]) * dx + (py - V[ia + 1]) * dy + (pz - V[ia + 2]) * dz) / L2 : 0; t = t < 0 ? 0 : t > 1 ? 1 : t;
  Q[0] = V[ia] + t * dx; Q[1] = V[ia + 1] + t * dy; Q[2] = V[ia + 2] + t * dz;
  return (px - Q[0]) ** 2 + (py - Q[1]) ** 2 + (pz - Q[2]) ** 2;
}

// ------------------------------------------------------------------ input features
/** feature chains of the input: boundary + creases, split at ends / junctions / corners; short crease chains
 *  (< minCreaseLen) dropped; apices |K| > sharpDeg fixed. Returns per-corner chain id, per-vertex rank / chain. */
/*  fixed / soft (per vertex, optional): vertices forced to be fixed (the ends of a non-orientable mesh's cut on the rim) /
 *  never made corners by the turning test (the cut itself: its zigzag along the input edges is not a feature). */
export function detectFeatures(V, F, O, { creaseDeg = 40, sharpDeg = 20, cornerDeg = 45, minCreaseLen = 0, window = 0, fixed = null, soft = null } = {}) {
  const n = V.length / 3, nc = F.length, N = faceNormals(V, F), cosCr = Math.cos(creaseDeg * Math.PI / 180);
  const flag = new Uint8Array(nc), VB = new Uint8Array(n), degen = new Uint8Array(nc / 3);
  for (let f = 0; f < nc / 3; f++) {                                        // needles / zero-area faces have no normal:
    const a = 3 * F[3 * f], b = 3 * F[3 * f + 1], d = 3 * F[3 * f + 2];       // no crease next to them
    const ux = V[b] - V[a], uy = V[b + 1] - V[a + 1], uz = V[b + 2] - V[a + 2], wx = V[d] - V[a], wy = V[d + 1] - V[a + 1], wz = V[d + 2] - V[a + 2];
    const A2 = (uy * wz - uz * wy) ** 2 + (uz * wx - ux * wz) ** 2 + (ux * wy - uy * wx) ** 2;
    const L2 = Math.max(ux * ux + uy * uy + uz * uz, wx * wx + wy * wy + wz * wz, (wx - ux) ** 2 + (wy - uy) ** 2 + (wz - uz) ** 2);
    if (A2 < 1e-12 * L2 * L2) degen[f] = 1;
  }
  for (let c = 0; c < nc; c++) {
    const o = O[c];
    if (o < 0) { flag[c] = 1; VB[F[nxt(c)]] = 1; VB[F[prv(c)]] = 1; continue; }
    if (o < c) continue;
    const f = (c / 3) | 0, g = (o / 3) | 0;
    if (degen[f] || degen[g]) continue;
    if (N[3 * f] * N[3 * g] + N[3 * f + 1] * N[3 * g + 1] + N[3 * f + 2] * N[3 * g + 2] < cosCr) { flag[c] = 1; flag[o] = 1; }
  }
  // angle sums (apex detection, corner valence targets)
  const ang = new Float64Array(n);
  for (let c = 0; c < nc; c++) {
    const a = 3 * F[c], b = 3 * F[nxt(c)], d = 3 * F[prv(c)];
    const ux = V[b] - V[a], uy = V[b + 1] - V[a + 1], uz = V[b + 2] - V[a + 2], wx = V[d] - V[a], wy = V[d + 1] - V[a + 1], wz = V[d + 2] - V[a + 2];
    const cs = (ux * wx + uy * wy + uz * wz) / Math.max(Math.hypot(ux, uy, uz) * Math.hypot(wx, wy, wz), 1e-300);
    ang[F[c]] += Math.acos(Math.max(-1, Math.min(1, cs)));
  }
  const apex = new Uint8Array(n); let nApex = 0;
  for (let v = 0; v < n; v++) if (!VB[v] && Math.abs(2 * Math.PI - ang[v]) > sharpDeg * Math.PI / 180) { apex[v] = 1; nApex++; }
  let chains, chainOf, rank, VL, dropped = 0, final = !(minCreaseLen > 0), nTurn = 0;
  const cosCorner = Math.cos(cornerDeg * Math.PI / 180);
  for (;;) {
    // feature edges (one representative corner each) and per-vertex incidence
    const edges = []; for (let c = 0; c < nc; c++) if (flag[c] && (O[c] < 0 || c < O[c])) edges.push(c);
    const deg = new Int32Array(n + 1); for (const c of edges) { deg[F[nxt(c)] + 1]++; deg[F[prv(c)] + 1]++; }
    for (let i = 0; i < n; i++) deg[i + 1] += deg[i];
    const pos = deg.slice(0, n), inc = new Int32Array(2 * edges.length);
    edges.forEach((c, e) => { inc[pos[F[nxt(c)]]++] = e; inc[pos[F[prv(c)]]++] = e; });
    const other = (e, v) => { const c = edges[e], a = F[nxt(c)]; return a === v ? F[prv(c)] : a; };
    rank = new Uint8Array(n);
    for (let v = 0; v < n; v++) { const k = deg[v + 1] - deg[v]; rank[v] = apex[v] || (k > 0 && k !== 2) || (fixed && fixed[v] && k > 0) ? 3 : k === 2 ? 2 : 0; }
    // walk chains between fixed vertices (and closed loops without any); a chain keeps its vertex sequence
    const walkAll = () => {
      chainOf = new Int32Array(edges.length).fill(-1); chains = [];
      const walk = (e0, v0) => {
        const id = chains.length, ch = { edges: [], verts: [v0], length: 0, boundary: false, closed: false }; chains.push(ch);
        let e = e0, v = v0;
        for (;;) {
          chainOf[e] = id; ch.edges.push(e); if (O[edges[e]] < 0) ch.boundary = true;
          const w = other(e, v); ch.length += Math.hypot(V[3 * w] - V[3 * v], V[3 * w + 1] - V[3 * v + 1], V[3 * w + 2] - V[3 * v + 2]);
          if (w === v0) { ch.closed = true; break; }
          ch.verts.push(w);
          if (rank[w] !== 2) break;
          const e1 = inc[deg[w]] === e ? inc[deg[w] + 1] : inc[deg[w]];
          if (chainOf[e1] >= 0) break;
          e = e1; v = w;
        }
      };
      for (let v = 0; v < n; v++) if (rank[v] === 3) for (let i = deg[v]; i < deg[v + 1]; i++) if (chainOf[inc[i]] < 0) walk(inc[i], v);
      for (let e = 0; e < edges.length; e++) if (chainOf[e] < 0) walk(e, F[nxt(edges[e])]);
    };
    walkAll();
    // corners on a chain: turning between the chord a window behind and a window ahead (so the zigzag of a dense
    // crease does not read as corners), one per window (the sharpest)
    let added = 0;
    for (const ch of chains) {
      const vs = ch.verts, m = vs.length; if (m < 3 && !ch.closed) continue;
      const S = new Float64Array(m + 1);
      for (let i = 1; i <= m; i++) { const a = vs[i - 1], b = vs[i % m]; if (i === m && !ch.closed) break; S[i] = S[i - 1] + Math.hypot(V[3 * b] - V[3 * a], V[3 * b + 1] - V[3 * a + 1], V[3 * b + 2] - V[3 * a + 2]); }
      const Ltot = ch.closed ? S[m] : S[m - 1];
      const at = (i) => vs[((i % m) + m) % m];
      const arc = (i, j) => { const d = Math.abs(S[((i % m) + m) % m] - S[((j % m) + m) % m]); return ch.closed ? Math.min(d, Ltot - d) : d; };
      const cand = [];
      for (let i = ch.closed ? 0 : 1; i < (ch.closed ? m : m - 1); i++) {
        let j = i - 1, k = i + 1;
        while ((ch.closed || j > 0) && arc(i, j) < window && i - j < m - 1) j--;
        while ((ch.closed || k < m - 1) && arc(i, k) < window && k - i < m - 1) k++;
        const v = at(i), x = at(j), y = at(k);
        if (soft && soft[v]) continue;
        const ax = V[3 * v] - V[3 * x], ay = V[3 * v + 1] - V[3 * x + 1], az = V[3 * v + 2] - V[3 * x + 2];
        const bx = V[3 * y] - V[3 * v], by = V[3 * y + 1] - V[3 * v + 1], bz = V[3 * y + 2] - V[3 * v + 2];
        const cs = (ax * bx + ay * by + az * bz) / Math.max(Math.hypot(ax, ay, az) * Math.hypot(bx, by, bz), 1e-300);
        if (cs < cosCorner) cand.push([cs, i]);
      }
      cand.sort((p, q) => p[0] - q[0]); const acc = [];
      for (const [, i] of cand) if (acc.every(j => arc(i, j) >= window)) { acc.push(i); if (rank[at(i)] === 2) { rank[at(i)] = 3; added++; } }
    }
    nTurn += added;
    if (added) walkAll();
    if (final) {
      VL = new Int32Array(n).fill(-1);
      const FL = new Int32Array(nc).fill(-1);
      edges.forEach((c, e) => { FL[c] = chainOf[e]; if (O[c] >= 0) FL[O[c]] = chainOf[e]; const a = F[nxt(c)], b = F[prv(c)]; if (rank[a] === 2) VL[a] = chainOf[e]; if (rank[b] === 2) VL[b] = chainOf[e]; });
      const tgt = new Uint8Array(n), imp = new Float64Array(n);
      for (let v = 0; v < n; v++) if (rank[v] === 3) {
        tgt[v] = Math.max(VB[v] ? 2 : 3, Math.min(12, Math.round(ang[v] / (Math.PI / 3)) + (VB[v] ? 1 : 0)));
        imp[v] = Math.abs((VB[v] ? Math.PI : 2 * Math.PI) - ang[v]) + (deg[v + 1] - deg[v]);
      }
      return { FL, rank, VL, VB, tgt, imp, chains, nApex, dropped, nFixed: rank.reduce((s, r) => s + (r === 3), 0) };
    }
    // drop short crease chains and redo once
    for (const ch of chains) if (!ch.boundary && ch.length < minCreaseLen) { dropped++; for (const e of ch.edges) { const c = edges[e]; flag[c] = 0; if (O[c] >= 0) flag[O[c]] = 0; } }
    final = true;                                                           // rebuild once (or just emit when nothing was dropped)
  }
}

// ------------------------------------------------------------------ the remesher
/** V0, F0: manifold, consistently oriented input (O0 its corner table, optional). Returns {V, F, stats}. */
export function remeshIsotropic(V0, F0, h, opts = {}) {
  const iters = opts.iters ?? 10, progress = opts.progress || (() => {});
  const lo2 = (0.8 * h) ** 2, hi2 = (4 / 3 * h) ** 2, hard2 = (0.5 * h) ** 2;
  const n0 = V0.length / 3, m0 = F0.length / 3;
  const O0 = opts.O || edgeMates(F0, n0, true).O;
  const feat = detectFeatures(V0, F0, O0, { creaseDeg: opts.creaseDeg ?? 40, sharpDeg: opts.sharpDeg ?? 20, cornerDeg: opts.cornerDeg ?? 45, minCreaseLen: opts.minCreaseLen ?? 1.5 * h, window: 0.5 * h, fixed: opts.fixed || null, soft: opts.soft || null });

  // --- projection targets: input triangles (with unit normals) and input feature segments (with chain id)
  const TN = faceNormals(V0, F0);
  const tlo = new Float64Array(3 * m0), thi = new Float64Array(3 * m0);
  for (let f = 0; f < m0; f++) for (let k = 0; k < 3; k++) {
    const a = V0[3 * F0[3 * f] + k], b = V0[3 * F0[3 * f + 1] + k], c = V0[3 * F0[3 * f + 2] + k];
    tlo[3 * f + k] = Math.min(a, b, c); thi[3 * f + k] = Math.max(a, b, c);
  }
  const triT = buildBVH(tlo, thi, m0);
  const segA = [], segB = [], segC = [];
  for (let c = 0; c < 3 * m0; c++) if (feat.FL[c] >= 0 && (O0[c] < 0 || c < O0[c])) { segA.push(F0[nxt(c)]); segB.push(F0[prv(c)]); segC.push(feat.FL[c]); }
  const ns0 = segA.length, slo = new Float64Array(3 * ns0), shi = new Float64Array(3 * ns0);
  for (let s = 0; s < ns0; s++) for (let k = 0; k < 3; k++) { const a = V0[3 * segA[s] + k], b = V0[3 * segB[s] + k]; slo[3 * s + k] = Math.min(a, b); shi[3 * s + k] = Math.max(a, b); }
  const segT = ns0 ? buildBVH(slo, shi, ns0) : null;

  // --- working mesh
  let vcap = Math.max(64, Math.ceil(n0 * 1.2)), fcap = Math.max(64, Math.ceil(m0 * 1.2));
  let P = new Float64Array(3 * vcap), VC = new Int32Array(vcap).fill(-1), RK = new Uint8Array(vcap), VL = new Int32Array(vcap).fill(-1);
  let VB = new Uint8Array(vcap), VAL = new Int32Array(vcap), TG = new Uint8Array(vcap), MK = new Int32Array(vcap), MK2 = new Int32Array(vcap), IM = new Float64Array(vcap);
  let C = new Int32Array(3 * fcap).fill(-1), O = new Int32Array(3 * fcap).fill(-1), FL = new Int32Array(3 * fcap).fill(-1);
  P.set(V0); RK.set(feat.rank); VL.set(feat.VL); VB.set(feat.VB); TG.set(feat.tgt); IM.set(feat.imp); C.set(F0); O.set(O0); FL.set(feat.FL);
  let nv = n0, nf = m0, stamp = 0, stamp2 = 0;
  for (let c = 0; c < 3 * nf; c++) { VC[C[c]] = c; VAL[C[c]]++; }
  for (let v = 0; v < nv; v++) if (VB[v]) VAL[v]++;                          // boundary: faces + 1 neighbours

  const growV = (need) => {
    if (need <= vcap) return; const c2 = Math.max(need, vcap * 2);
    const g = (A, k, fill) => { const B = new A.constructor(k * c2); if (fill !== undefined) B.fill(fill); B.set(A); return B; };
    P = g(P, 3); VC = g(VC, 1, -1); RK = g(RK, 1); VL = g(VL, 1, -1); VB = g(VB, 1); VAL = g(VAL, 1); TG = g(TG, 1); MK = g(MK, 1); MK2 = g(MK2, 1); IM = g(IM, 1); vcap = c2;
  };
  const growF = (need) => {
    if (need <= fcap) return; const c2 = Math.max(need, fcap * 2);
    const g = (A) => { const B = new Int32Array(3 * c2).fill(-1); B.set(A); return B; };
    C = g(C); O = g(O); FL = g(FL); fcap = c2;
  };
  const d2v = (a, b) => (P[3 * a] - P[3 * b]) ** 2 + (P[3 * a + 1] - P[3 * b + 1]) ** 2 + (P[3 * a + 2] - P[3 * b + 2]) ** 2;
  const cornerOf = (f, v) => (C[3 * f] === v ? 3 * f : C[3 * f + 1] === v ? 3 * f + 1 : C[3 * f + 2] === v ? 3 * f + 2 : -1);
  const ringInto = (v, out) => {
    out.length = 0; const s = VC[v]; let c = s, guard = 0;
    for (;;) { out.push(c); const o = O[prv(c)]; if (o < 0) { c = -1; break; } c = prv(o); if (c === s) break; if (++guard > 1e6) throw new Error('ring loop'); }
    if (c === -1) { let o = O[nxt(s)]; while (o >= 0) { const d = nxt(o); out.push(d); o = O[nxt(d)]; if (++guard > 1e6) throw new Error('ring loop'); } }
    return out;
  };
  const target = (v) => (RK[v] === 3 ? TG[v] : VB[v] ? 4 : 6);
  const R1 = [], R2 = [];

  // --- split
  function splitEdge(c) {
    growV(nv + 1); growF(nf + 2);
    const n = nxt(c), p = prv(c), a = C[c], b = C[n], d = C[p], o = O[c], line = FL[c];
    const m = nv++;
    P[3 * m] = 0.5 * (P[3 * b] + P[3 * d]); P[3 * m + 1] = 0.5 * (P[3 * b + 1] + P[3 * d + 1]); P[3 * m + 2] = 0.5 * (P[3 * b + 2] + P[3 * d + 2]);
    RK[m] = line >= 0 ? 2 : 0; VL[m] = line; VB[m] = o < 0 ? 1 : 0; TG[m] = 0; IM[m] = 0;
    const f2 = nf++, r0 = 3 * f2, r1 = r0 + 1, r2 = r0 + 2, xn = O[n], fln = FL[n];
    C[r0] = a; C[r1] = m; C[r2] = d; C[p] = m;
    O[n] = r2; O[r2] = n; FL[n] = -1; FL[r2] = -1;
    O[r1] = xn; if (xn >= 0) O[xn] = r1; FL[r1] = fln;
    FL[c] = line; FL[r0] = line;
    if (o >= 0) {
      const no = nxt(o), po = prv(o), e = C[o], yno = O[no], flno = FL[no];
      const g2 = nf++, q0 = 3 * g2, q1 = q0 + 1, q2 = q0 + 2;
      C[q0] = e; C[q1] = m; C[q2] = b; C[po] = m;
      O[c] = q0; O[q0] = c; FL[q0] = line;
      O[r0] = o; O[o] = r0; FL[o] = line;
      O[no] = q2; O[q2] = no; FL[no] = -1; FL[q2] = -1;
      O[q1] = yno; if (yno >= 0) O[yno] = q1; FL[q1] = flno;
      VAL[e]++; VAL[m] = 4;
    } else { O[c] = -1; O[r0] = -1; VAL[m] = 3; }
    VAL[a]++; VC[m] = r1; VC[d] = r2; VC[b] = n;
  }
  function splitPass() {
    let total = 0;
    for (let round = 0; round < 30; round++) {
      let cnt = 0;
      for (let f = 0; f < nf; f++) {
        if (C[3 * f] < 0) continue;
        for (let guard = 0; guard < 64; guard++) {
          let best = -1, bl = hi2;
          for (let s = 0; s < 3; s++) { const c = 3 * f + s, l = d2v(C[nxt(c)], C[prv(c)]); if (l > bl) { bl = l; best = c; } }
          if (best < 0) break; splitEdge(best); cnt++;
        }
      }
      total += cnt; if (!cnt) break;
    }
    return total;
  }

  // --- collapse
  const faceOK = (q, px, py, pz) => {                                       // corner q moves to p: no flip, not degenerate
    const v = 3 * C[q], x = 3 * C[nxt(q)], y = 3 * C[prv(q)];
    const ux = P[x] - P[v], uy = P[x + 1] - P[v + 1], uz = P[x + 2] - P[v + 2], wx = P[y] - P[v], wy = P[y + 1] - P[v + 1], wz = P[y + 2] - P[v + 2];
    const n0x = uy * wz - uz * wy, n0y = uz * wx - ux * wz, n0z = ux * wy - uy * wx;
    const sx = P[x] - px, sy = P[x + 1] - py, sz = P[x + 2] - pz, tx = P[y] - px, ty = P[y + 1] - py, tz = P[y + 2] - pz;
    const n1x = sy * tz - sz * ty, n1y = sz * tx - sx * tz, n1z = sx * ty - sy * tx;
    const l1 = Math.hypot(n1x, n1y, n1z); if (l1 < 1e-10 * hi2) return false;
    return n0x * n1x + n0y * n1y + n0z * n1z > 0.3 * Math.hypot(n0x, n0y, n0z) * l1;
  };
  function tryCollapse(c) {
    const cu = nxt(c), cw = prv(c), u = C[cu], w = C[cw], fe = FL[c] >= 0, ru = RK[u], rw = RK[w];
    let r, k, mid = false;
    if (ru === 0 && rw === 0) { r = u; k = w; mid = true; }
    else if (ru === 0) { r = u; k = w; }
    else if (rw === 0) { r = w; k = u; }
    else if (ru === 3 && rw === 3) {                                        // two fixed vertices closer than h/2: merge
      if (d2v(u, w) > hard2) return false;                                  // into the more important one
      if (IM[u] <= IM[w]) { r = u; k = w; } else { r = w; k = u; }
    } else if (fe) {
      if (ru === 2 && rw === 2) { if (VL[u] !== VL[w]) return false; r = u; k = w; mid = true; }
      else if (ru === 2) { r = u; k = w; } else { r = w; k = u; }
    } else {                                                                // a feature vertex next to another chain's fixed
      if (d2v(u, w) > hard2 || (ru === 2 && rw === 2)) return false;        // vertex (closer than h/2) joins it
      if (ru === 2) { r = u; k = w; } else { r = w; k = u; }
    }
    const o = O[c], a = C[c], e = o >= 0 ? C[o] : -1;
    if (o >= 0 && VB[u] && VB[w]) return false;
    if (VAL[a] - 1 < (VB[a] ? 2 : 3)) return false;
    if (e >= 0 && VAL[e] - 1 < (VB[e] ? 2 : 3)) return false;
    ringInto(k, R1); stamp++;
    for (const q of R1) { MK[C[nxt(q)]] = stamp; MK[C[prv(q)]] = stamp; }
    ringInto(r, R2); stamp2++; let common = 0;
    for (const q of R2) for (let t = 0; t < 2; t++) {
      const x = t ? C[prv(q)] : C[nxt(q)]; if (MK2[x] === stamp2) continue; MK2[x] = stamp2;
      if (x !== k && MK[x] === stamp) { if (x !== a && x !== e) return false; common++; }
    }
    if (common !== (o >= 0 ? 2 : 1)) return false;
    if (VAL[k] + VAL[r] - common - 2 < (VB[k] || VB[r] ? 2 : 3)) return false;
    let px, py, pz;
    if (mid) { px = 0.5 * (P[3 * u] + P[3 * w]); py = 0.5 * (P[3 * u + 1] + P[3 * w + 1]); pz = 0.5 * (P[3 * u + 2] + P[3 * w + 2]); }
    else { px = P[3 * k]; py = P[3 * k + 1]; pz = P[3 * k + 2]; }
    const f = (c / 3) | 0, g = o >= 0 ? (o / 3) | 0 : -1;
    const far = (x) => (P[3 * x] - px) ** 2 + (P[3 * x + 1] - py) ** 2 + (P[3 * x + 2] - pz) ** 2 > hi2;
    for (const q of R2) { const ff = (q / 3) | 0; if (ff === f || ff === g) continue; if (far(C[nxt(q)]) || far(C[prv(q)]) || !faceOK(q, px, py, pz)) return false; }
    if (mid) for (const q of R1) { const ff = (q / 3) | 0; if (ff === f || ff === g) continue; if (far(C[nxt(q)]) || far(C[prv(q)]) || !faceOK(q, px, py, pz)) return false; }
    // --- do it
    const cr = u === r ? cu : cw, ck = cr === cu ? cw : cu, xr = O[ck], xk = O[cr], fl = FL[ck] >= 0 ? FL[ck] : FL[cr];
    if (xr >= 0) { O[xr] = xk; FL[xr] = fl; } if (xk >= 0) { O[xk] = xr; FL[xk] = fl; }
    let yr = -1, yk = -1;
    if (o >= 0) {
      const o1 = nxt(o), o2 = prv(o), gr = C[o1] === r ? o1 : o2, gk = gr === o1 ? o2 : o1, fl2 = FL[gk] >= 0 ? FL[gk] : FL[gr];
      yr = O[gk]; yk = O[gr];
      if (yr >= 0) { O[yr] = yk; FL[yr] = fl2; } if (yk >= 0) { O[yk] = yr; FL[yk] = fl2; }
    }
    for (const q of R2) C[q] = k;
    for (const ff of [f, g]) if (ff >= 0) for (let s = 0; s < 3; s++) { C[3 * ff + s] = -1; O[3 * ff + s] = -1; FL[3 * ff + s] = -1; }
    VC[r] = -1; P[3 * k] = px; P[3 * k + 1] = py; P[3 * k + 2] = pz; VB[k] |= VB[r];
    VAL[k] = VAL[k] + VAL[r] - common - 2; VAL[a]--; if (e >= 0) VAL[e]--;
    const fa = ((xk >= 0 ? xk : xr) / 3) | 0; VC[k] = cornerOf(fa, k); VC[a] = cornerOf(fa, a);
    if (o >= 0) { const fb = ((yk >= 0 ? yk : yr) / 3) | 0; VC[e] = cornerOf(fb, e); }
    return true;
  }
  function collapsePass() {
    let total = 0;
    for (let round = 0; round < 8; round++) {
      let cnt = 0;
      for (let f = 0; f < nf; f++) {
        for (let s = 0; s < 3; s++) {
          if (C[3 * f] < 0) break;
          const c = 3 * f + s; if (d2v(C[nxt(c)], C[prv(c)]) < lo2 && tryCollapse(c)) { cnt++; s = -1; }
        }
      }
      total += cnt; if (cnt <= 0.002 * nf) break;
    }
    return total;
  }

  // --- flips
  const nrm = (a, b, d, out) => {
    const ux = P[3 * b] - P[3 * a], uy = P[3 * b + 1] - P[3 * a + 1], uz = P[3 * b + 2] - P[3 * a + 2];
    const wx = P[3 * d] - P[3 * a], wy = P[3 * d + 1] - P[3 * a + 1], wz = P[3 * d + 2] - P[3 * a + 2];
    out[0] = uy * wz - uz * wy; out[1] = uz * wx - ux * wz; out[2] = ux * wy - uy * wx; return Math.hypot(out[0], out[1], out[2]);
  };
  const N1 = [0, 0, 0], N2 = [0, 0, 0], N3 = [0, 0, 0], N4 = [0, 0, 0];
  function tryFlip(c) {
    const o = O[c]; if (o < 0 || FL[c] >= 0) return false;
    const n = nxt(c), p = prv(c), no = nxt(o), po = prv(o), a = C[c], b = C[n], d = C[p], e = C[o];
    if (a === e) return false;
    const ta = target(a), tb = target(b), td = target(d), te = target(e);
    const before = (VAL[a] - ta) ** 2 + (VAL[b] - tb) ** 2 + (VAL[d] - td) ** 2 + (VAL[e] - te) ** 2;          // squared deviation
    const after = (VAL[a] + 1 - ta) ** 2 + (VAL[b] - 1 - tb) ** 2 + (VAL[d] - 1 - td) ** 2 + (VAL[e] + 1 - te) ** 2; // (as pmp-library)
    if (after >= before) return false;
    if (VAL[b] - 1 < (VB[b] ? 2 : 3) || VAL[d] - 1 < (VB[d] ? 2 : 3)) return false;
    ringInto(a, R1); for (const q of R1) if (C[nxt(q)] === e || C[prv(q)] === e) return false;
    const lf = nrm(a, b, d, N1), lg = nrm(e, d, b, N2), l1 = nrm(a, b, e, N3), l2 = nrm(e, d, a, N4);
    if (l1 < 1e-10 * hi2 || l2 < 1e-10 * hi2 || lf <= 0 || lg <= 0) return false;
    const sx = N1[0] / lf + N2[0] / lg, sy = N1[1] / lf + N2[1] / lg, sz = N1[2] / lf + N2[2] / lg;
    if (N3[0] * sx + N3[1] * sy + N3[2] * sz <= 0 || N4[0] * sx + N4[1] * sy + N4[2] * sz <= 0) return false;
    const cb = (N1[0] * N2[0] + N1[1] * N2[1] + N1[2] * N2[2]) / (lf * lg), ca = (N3[0] * N4[0] + N3[1] * N4[1] + N3[2] * N4[2]) / (l1 * l2);
    if (ca < Math.min(cb, 0.5)) return false;
    const x = O[n], y = O[no], fx = FL[n], fy = FL[no];
    C[p] = e; C[po] = a;
    O[c] = y; if (y >= 0) O[y] = c; FL[c] = fy;
    O[o] = x; if (x >= 0) O[x] = o; FL[o] = fx;
    O[n] = no; O[no] = n; FL[n] = -1; FL[no] = -1;
    VC[b] = n; VC[d] = no; VC[a] = c; VC[e] = o;
    VAL[a]++; VAL[e]++; VAL[b]--; VAL[d]--;
    return true;
  }
  function flipPass() {
    let cnt = 0;
    for (let c = 0; c < 3 * nf; c++) { if (C[c] < 0) continue; const o = O[c]; if (o > c && tryFlip(c)) cnt++; }
    return cnt;
  }

  // --- normals, relaxation, projection
  let NV = new Float64Array(0), AR = new Float64Array(0);
  function vertexNormals() {
    if (NV.length < 3 * nv) { NV = new Float64Array(3 * vcap); AR = new Float64Array(vcap); }
    NV.fill(0, 0, 3 * nv); AR.fill(0, 0, nv);
    for (let f = 0; f < nf; f++) {
      if (C[3 * f] < 0) continue; const l = nrm(C[3 * f], C[3 * f + 1], C[3 * f + 2], N1);
      for (let s = 0; s < 3; s++) { const v = C[3 * f + s]; NV[3 * v] += N1[0]; NV[3 * v + 1] += N1[1]; NV[3 * v + 2] += N1[2]; AR[v] += l / 6; }
    }
    for (let v = 0; v < nv; v++) { const l = Math.hypot(NV[3 * v], NV[3 * v + 1], NV[3 * v + 2]) || 1; NV[3 * v] /= l; NV[3 * v + 1] /= l; NV[3 * v + 2] /= l; }
  }
  function relax() {
    vertexNormals();
    const acc = new Float64Array(3 * nv), ws = new Float64Array(nv), fn1 = new Int32Array(nv).fill(-1), fn2 = new Int32Array(nv).fill(-1), bad = new Uint8Array(nv);
    const addF = (v, x) => { if (fn1[v] < 0) fn1[v] = x; else if (fn1[v] !== x) { if (fn2[v] < 0) fn2[v] = x; else if (fn2[v] !== x) bad[v] = 1; } };
    for (let q = 0; q < 3 * nf; q++) {
      const v = C[q]; if (v < 0) continue; const x = C[nxt(q)], y = C[prv(q)];
      if (RK[v] === 0) {
        acc[3 * v] += AR[x] * P[3 * x] + AR[y] * P[3 * y]; acc[3 * v + 1] += AR[x] * P[3 * x + 1] + AR[y] * P[3 * y + 1]; acc[3 * v + 2] += AR[x] * P[3 * x + 2] + AR[y] * P[3 * y + 2];
        ws[v] += AR[x] + AR[y];
      }
      if (FL[q] >= 0) { addF(x, y); addF(y, x); }
    }
    const NP = P.slice(0, 3 * nv);
    for (let v = 0; v < nv; v++) {
      if (VC[v] < 0) continue;
      if (RK[v] === 0 && ws[v] > 0) {
        const qx = acc[3 * v] / ws[v], qy = acc[3 * v + 1] / ws[v], qz = acc[3 * v + 2] / ws[v];
        const nx = NV[3 * v], ny = NV[3 * v + 1], nz = NV[3 * v + 2];
        const dn = nx * (P[3 * v] - qx) + ny * (P[3 * v + 1] - qy) + nz * (P[3 * v + 2] - qz);
        NP[3 * v] = qx + dn * nx; NP[3 * v + 1] = qy + dn * ny; NP[3 * v + 2] = qz + dn * nz;
      } else if (RK[v] === 2 && !bad[v] && fn1[v] >= 0 && fn2[v] >= 0) {
        const a = fn1[v], b = fn2[v];
        let tx = P[3 * b] - P[3 * a], ty = P[3 * b + 1] - P[3 * a + 1], tz = P[3 * b + 2] - P[3 * a + 2]; const tl = Math.hypot(tx, ty, tz); if (!(tl > 0)) continue;
        tx /= tl; ty /= tl; tz /= tl;
        const mx = 0.5 * (P[3 * a] + P[3 * b]) - P[3 * v], my = 0.5 * (P[3 * a + 1] + P[3 * b + 1]) - P[3 * v + 1], mz = 0.5 * (P[3 * a + 2] + P[3 * b + 2]) - P[3 * v + 2];
        const s = mx * tx + my * ty + mz * tz;
        NP[3 * v] = P[3 * v] + s * tx; NP[3 * v + 1] = P[3 * v + 1] + s * ty; NP[3 * v + 2] = P[3 * v + 2] + s * tz;
      }
    }
    P.set(NP);
  }
  let segChain = -1;
  const triTestFiltered = (f, best) => {
    const d = closestOnTri(PX, PY, PZ, V0, 3 * F0[3 * f], 3 * F0[3 * f + 1], 3 * F0[3 * f + 2]);
    if (d >= best) return Infinity;
    return TN[3 * f] * NX + TN[3 * f + 1] * NY + TN[3 * f + 2] * NZ > 0 ? d : Infinity;
  };
  const triTest = (f) => closestOnTri(PX, PY, PZ, V0, 3 * F0[3 * f], 3 * F0[3 * f + 1], 3 * F0[3 * f + 2]);
  const segTest = (s) => (segC[s] === segChain ? closestOnSeg(PX, PY, PZ, V0, 3 * segA[s], 3 * segB[s]) : Infinity);
  let PX = 0, PY = 0, PZ = 0, NX = 0, NY = 0, NZ = 0, maxDev = 0;
  function project() {
    vertexNormals(); maxDev = 0;
    for (let v = 0; v < nv; v++) {
      if (VC[v] < 0 || RK[v] === 3) continue;
      PX = P[3 * v]; PY = P[3 * v + 1]; PZ = P[3 * v + 2];
      if (RK[v] === 2 && segT) {
        segChain = VL[v]; const r = bvhQuery(segT, PX, PY, PZ, Infinity, segTest);
        if (r.prim >= 0) { closestOnSeg(PX, PY, PZ, V0, 3 * segA[r.prim], 3 * segB[r.prim]); P[3 * v] = Q[0]; P[3 * v + 1] = Q[1]; P[3 * v + 2] = Q[2]; maxDev = Math.max(maxDev, Math.sqrt(r.d2)); continue; }
      }
      NX = NV[3 * v]; NY = NV[3 * v + 1]; NZ = NV[3 * v + 2];
      const r0 = bvhQuery(triT, PX, PY, PZ, Infinity, triTest);
      const lim = (Math.sqrt(r0.d2) + 0.5 * h) ** 2;
      const r1 = bvhQuery(triT, PX, PY, PZ, lim, triTestFiltered);
      const pick = r1.prim >= 0 ? r1.prim : r0.prim; if (pick < 0) continue;
      const d = closestOnTri(PX, PY, PZ, V0, 3 * F0[3 * pick], 3 * F0[3 * pick + 1], 3 * F0[3 * pick + 2]);
      P[3 * v] = Q[0]; P[3 * v + 1] = Q[1]; P[3 * v + 2] = Q[2]; maxDev = Math.max(maxDev, Math.sqrt(d));
    }
  }

  // --- compaction (drop dead faces / vertices)
  function compact() {
    const fmap = new Int32Array(nf).fill(-1), vmap = new Int32Array(nv).fill(-1); let mf = 0, mv = 0;
    for (let f = 0; f < nf; f++) if (C[3 * f] >= 0) fmap[f] = mf++;
    for (let v = 0; v < nv; v++) if (VC[v] >= 0) vmap[v] = mv++;
    const C2 = new Int32Array(3 * fcap).fill(-1), O2 = new Int32Array(3 * fcap).fill(-1), FL2 = new Int32Array(3 * fcap).fill(-1);
    for (let f = 0; f < nf; f++) {
      const g = fmap[f]; if (g < 0) continue;
      for (let s = 0; s < 3; s++) {
        const c = 3 * f + s, o = O[c]; C2[3 * g + s] = vmap[C[c]]; FL2[3 * g + s] = FL[c];
        O2[3 * g + s] = o >= 0 ? 3 * fmap[(o / 3) | 0] + (o % 3) : -1;
      }
    }
    for (let v = 0; v < nv; v++) {
      const w = vmap[v]; if (w < 0) continue;
      P[3 * w] = P[3 * v]; P[3 * w + 1] = P[3 * v + 1]; P[3 * w + 2] = P[3 * v + 2];
      VC[w] = 3 * fmap[(VC[v] / 3) | 0] + (VC[v] % 3); RK[w] = RK[v]; VL[w] = VL[v]; VB[w] = VB[v]; VAL[w] = VAL[v]; TG[w] = TG[v]; IM[w] = IM[v];
    }
    VC.fill(-1, mv);
    C = C2; O = O2; FL = FL2; nf = mf; nv = mv;
  }

  // --- main loop
  const stats = { input: { n: n0, nf: m0 }, chains: feat.chains.length, fixed: feat.nFixed, apices: feat.nApex, droppedCreases: feat.dropped, iters: [] };
  for (let it = 0; it < iters; it++) {
    const t0 = Date.now();
    const ns = splitPass(); compact();
    const ncl = collapsePass(); compact();
    const nfl = flipPass();
    relax(); project();
    stats.iters.push({ split: ns, collapse: ncl, flip: nfl, n: nv, nf, maxDev, ms: Date.now() - t0 });
    progress((it + 1) / iters, `remesh ${it + 1}/${iters}: ${nv} vertices`);
  }
  compact();
  return { V: P.slice(0, 3 * nv), F: C.slice(0, 3 * nf), rank: RK.slice(0, nv), stats };
}

// ------------------------------------------------------------------ smoothing a patch (the re-glued seam)
/** Tangential smoothing of the vertices flagged in `verts` (the re-glued seam of a non-orientable carrier and its
 *  neighbourhood), each sweep followed by the projection back onto the input surface (V0, F0). The faces around the
 *  seam come with both windings, so a vertex's normal is its face normals turned to agree with its largest face.
 *  V is changed in place. */
export function relaxRegion(V, F, verts, V0, F0, { iters = 8, step = 0.5 } = {}) {
  const m = F.length / 3, m0 = F0.length / 3;
  const tlo = new Float64Array(3 * m0), thi = new Float64Array(3 * m0);
  for (let f = 0; f < m0; f++) for (let k = 0; k < 3; k++) {
    const a = V0[3 * F0[3 * f] + k], b = V0[3 * F0[3 * f + 1] + k], c = V0[3 * F0[3 * f + 2] + k];
    tlo[3 * f + k] = Math.min(a, b, c); thi[3 * f + k] = Math.max(a, b, c);
  }
  const T = buildBVH(tlo, thi, m0), nb = new Map(), fc = new Map();
  for (let f = 0; f < m; f++) for (let k = 0; k < 3; k++) {
    const v = F[3 * f + k]; if (!verts[v]) continue;
    let s = nb.get(v); if (!s) { nb.set(v, s = new Set()); fc.set(v, []); }
    s.add(F[3 * f + (k + 1) % 3]); s.add(F[3 * f + (k + 2) % 3]); fc.get(v).push(f);
  }
  const fn = (f) => {
    const a = 3 * F[3 * f], b = 3 * F[3 * f + 1], c = 3 * F[3 * f + 2];
    const ux = V[b] - V[a], uy = V[b + 1] - V[a + 1], uz = V[b + 2] - V[a + 2], wx = V[c] - V[a], wy = V[c + 1] - V[a + 1], wz = V[c + 2] - V[a + 2];
    return [uy * wz - uz * wy, uz * wx - ux * wz, ux * wy - uy * wx];
  };
  let moved = 0;
  for (let it = 0; it < iters; it++) {
    const next = [];
    for (const [v, s] of nb) {
      let cx = 0, cy = 0, cz = 0; for (const w of s) { cx += V[3 * w]; cy += V[3 * w + 1]; cz += V[3 * w + 2]; }
      cx /= s.size; cy /= s.size; cz /= s.size;
      const ns = fc.get(v).map(fn); let ref = ns[0], rl = 0; for (const q of ns) { const l = Math.hypot(...q); if (l > rl) { rl = l; ref = q; } }
      let nx = 0, ny = 0, nz = 0; for (const q of ns) { const sg = q[0] * ref[0] + q[1] * ref[1] + q[2] * ref[2] < 0 ? -1 : 1; nx += sg * q[0]; ny += sg * q[1]; nz += sg * q[2]; }
      const nl = Math.hypot(nx, ny, nz) || 1; nx /= nl; ny /= nl; nz /= nl;
      let dx = cx - V[3 * v], dy = cy - V[3 * v + 1], dz = cz - V[3 * v + 2]; const dn = dx * nx + dy * ny + dz * nz;
      dx -= dn * nx; dy -= dn * ny; dz -= dn * nz;
      next.push(v, V[3 * v] + step * dx, V[3 * v + 1] + step * dy, V[3 * v + 2] + step * dz);
    }
    for (let k = 0; k < next.length; k += 4) {
      const v = next[k], px = next[k + 1], py = next[k + 2], pz = next[k + 3];
      const r = bvhQuery(T, px, py, pz, Infinity, (f) => closestOnTri(px, py, pz, V0, 3 * F0[3 * f], 3 * F0[3 * f + 1], 3 * F0[3 * f + 2]));
      if (r.prim < 0) continue;
      closestOnTri(px, py, pz, V0, 3 * F0[3 * r.prim], 3 * F0[3 * r.prim + 1], 3 * F0[3 * r.prim + 2]);
      V[3 * v] = Q[0]; V[3 * v + 1] = Q[1]; V[3 * v + 2] = Q[2]; moved++;
    }
  }
  return { moved: nb.size };
}
