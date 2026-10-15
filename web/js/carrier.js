// Carrier mesh for the browser: loads prepared .lshc carriers (remeshed + intrinsic Delaunay in Python), or builds
// a carrier from a dropped STL (weld, orient, Loop subdivision, intrinsic Delaunay flips — no remesher here), then
// cotan Laplacian / mass / defects, singularity detection, heat-method distance and the n-fold lattice seed.
import { csrFromCOO, spmv, csrDiagonal, pcg, dijkstra } from './sparse.js';

export const SQ3 = Math.sqrt(3);

// ------------------------------------------------------------------ .lshc loader
export function parseLSHC(buf) {
  const dv = new DataView(buf);
  const magic = String.fromCharCode(dv.getUint8(0), dv.getUint8(1), dv.getUint8(2), dv.getUint8(3));
  if (magic !== 'LSHC') throw new Error('not an .lshc carrier');
  const hlen = dv.getUint32(4, true);
  const header = JSON.parse(new TextDecoder().decode(new Uint8Array(buf, 8, hlen)));
  const base = 8 + hlen;
  const out = { name: header.name, meta: header.meta };
  const ctor = { float32: Float32Array, uint32: Uint32Array, uint8: Uint8Array, int32: Int32Array };
  for (const [k, a] of Object.entries(header.arrays)) {
    const C = ctor[a.dtype]; out[k] = new C(buf.slice(base + a.offset, base + a.offset + a.bytes));
  }
  return out;
}

// ------------------------------------------------------------------ STL / mesh utilities
export function parseSTL(buf) {
  const u8 = new Uint8Array(buf);
  const head = new TextDecoder().decode(u8.subarray(0, Math.min(80, u8.length)));
  let tri;
  const n = buf.byteLength >= 84 ? new DataView(buf).getUint32(80, true) : 0;
  if (buf.byteLength === 84 + n * 50 && n > 0) {
    tri = new Float64Array(n * 9); const dv = new DataView(buf);
    for (let i = 0; i < n; i++) for (let k = 0; k < 9; k++) tri[i * 9 + k] = dv.getFloat32(84 + i * 50 + 12 + k * 4, true);
  } else {
    const txt = new TextDecoder().decode(u8); const pts = [];
    for (const line of txt.split('\n')) { const s = line.trim().split(/\s+/); if (s.length === 4 && s[0] === 'vertex') pts.push(+s[1], +s[2], +s[3]); }
    tri = Float64Array.from(pts);
  }
  const nv = tri.length / 3;
  const F = new Uint32Array(nv); for (let i = 0; i < nv; i++) F[i] = i;
  return weld(tri, F);
}

export function weld(V, F, relTol = 1e-6) {
  const nv = V.length / 3; let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < nv; i++) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], V[3 * i + k]); hi[k] = Math.max(hi[k], V[3 * i + k]); }
  const diag = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;
  const q = relTol * diag, map = new Map(), remap = new Int32Array(nv), out = [];
  for (let i = 0; i < nv; i++) {
    const key = Math.round(V[3 * i] / q) + ',' + Math.round(V[3 * i + 1] / q) + ',' + Math.round(V[3 * i + 2] / q);
    let id = map.get(key);
    if (id === undefined) { id = out.length / 3; map.set(key, id); out.push(V[3 * i], V[3 * i + 1], V[3 * i + 2]); }
    remap[i] = id;
  }
  const F2 = [];
  for (let f = 0; f < F.length / 3; f++) {
    const a = remap[F[3 * f]], b = remap[F[3 * f + 1]], c = remap[F[3 * f + 2]];
    if (a !== b && b !== c && a !== c) F2.push(a, b, c);
  }
  return { V: Float64Array.from(out), F: Uint32Array.from(F2) };
}

function edgeKey(a, b) { return a < b ? a * 4294967296 + b : b * 4294967296 + a; }

export function edgeTable(F, n) {
  // unique undirected edges with face counts; boundary vertex flags
  const nf = F.length / 3, map = new Map();
  for (let f = 0; f < nf; f++) for (let s = 0; s < 3; s++) {
    const a = F[3 * f + s], b = F[3 * f + (s + 1) % 3], k = edgeKey(a, b);
    const e = map.get(k); if (e) e.count++; else map.set(k, { a: Math.min(a, b), b: Math.max(a, b), count: 1 });
  }
  const boundary = new Uint8Array(n); let ne = 0;
  for (const e of map.values()) { ne++; if (e.count === 1) { boundary[e.a] = 1; boundary[e.b] = 1; } }
  return { map, ne, boundary };
}

export function orientConsistently(V, F) {
  const nf = F.length / 3, n = V.length / 3;
  const sideMap = new Map();   // undirected edge → faces
  for (let f = 0; f < nf; f++) for (let s = 0; s < 3; s++) {
    const k = edgeKey(F[3 * f + s], F[3 * f + (s + 1) % 3]); let l = sideMap.get(k); if (!l) sideMap.set(k, l = []); l.push(f);
  }
  const visited = new Uint8Array(nf); const Fo = Uint32Array.from(F);
  const hasDirected = (f, a, b) => (Fo[3 * f] === a && Fo[3 * f + 1] === b) || (Fo[3 * f + 1] === a && Fo[3 * f + 2] === b) || (Fo[3 * f + 2] === a && Fo[3 * f] === b);
  for (let root = 0; root < nf; root++) {
    if (visited[root]) continue; visited[root] = 1; const stack = [root];
    while (stack.length) {
      const f = stack.pop();
      for (let s = 0; s < 3; s++) {
        const a = Fo[3 * f + s], b = Fo[3 * f + (s + 1) % 3];
        for (const g of sideMap.get(edgeKey(a, b))) {
          if (g === f || visited[g]) continue;
          if (hasDirected(g, a, b)) { const t = Fo[3 * g + 1]; Fo[3 * g + 1] = Fo[3 * g + 2]; Fo[3 * g + 2] = t; }
          visited[g] = 1; stack.push(g);
        }
      }
    }
  }
  // global sign: closed → outward (signed volume); open → +z, else away from the centroid
  const { boundary } = edgeTable(Fo, n); let closed = true; for (let i = 0; i < n; i++) if (boundary[i]) { closed = false; break; }
  let vol = 0, sz = 0, tot = 0, cdot = 0; const cx = [0, 0, 0]; for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) cx[k] += V[3 * i + k] / n;
  for (let f = 0; f < nf; f++) {
    const a = Fo[3 * f], b = Fo[3 * f + 1], c = Fo[3 * f + 2];
    const ux = V[3 * b] - V[3 * a], uy = V[3 * b + 1] - V[3 * a + 1], uz = V[3 * b + 2] - V[3 * a + 2];
    const wx = V[3 * c] - V[3 * a], wy = V[3 * c + 1] - V[3 * a + 1], wz = V[3 * c + 2] - V[3 * a + 2];
    const nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx;
    vol += (V[3 * a] * nx + V[3 * a + 1] * ny + V[3 * a + 2] * nz) / 6; sz += nz / 2; tot += Math.hypot(nx, ny, nz) / 2;
    const mx = (V[3 * a] + V[3 * b] + V[3 * c]) / 3 - cx[0], my = (V[3 * a + 1] + V[3 * b + 1] + V[3 * c + 1]) / 3 - cx[1], mz = (V[3 * a + 2] + V[3 * b + 2] + V[3 * c + 2]) / 3 - cx[2];
    cdot += nx * mx + ny * my + nz * mz;
  }
  return orientSide(V, Fo);
}

/** Global sign of a consistently oriented mesh (same rule as carrier.py): closed → outward (signed volume); open →
 *  the convex side (normals against the mean-curvature vector: tubes, branching tubes, domes face outward), else +z
 *  (saddles, sheets, where the total mean curvature is ~0). */
export function orientSide(V, F) {
  const n = V.length / 3, nf = F.length / 3, Fo = Uint32Array.from(F);
  const { boundary } = edgeTable(Fo, n); let closed = true; for (let i = 0; i < n; i++) if (boundary[i]) { closed = false; break; }
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], V[3 * i + k]); hi[k] = Math.max(hi[k], V[3 * i + k]); }
  const bbox = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;
  const N = new Float64Array(3 * n), LP = new Float64Array(3 * n);          // area-weighted normals; cotan Laplacian of the positions
  let vol = 0, sz = 0, tot = 0;
  for (let f = 0; f < nf; f++) {
    const I = [Fo[3 * f], Fo[3 * f + 1], Fo[3 * f + 2]], a = I[0], b = I[1], c = I[2];
    const ux = V[3 * b] - V[3 * a], uy = V[3 * b + 1] - V[3 * a + 1], uz = V[3 * b + 2] - V[3 * a + 2];
    const wx = V[3 * c] - V[3 * a], wy = V[3 * c + 1] - V[3 * a + 1], wz = V[3 * c + 2] - V[3 * a + 2];
    const nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx;
    vol += (V[3 * a] * nx + V[3 * a + 1] * ny + V[3 * a + 2] * nz) / 6; sz += nz / 2; tot += Math.hypot(nx, ny, nz) / 2;
    for (const i of I) { N[3 * i] += nx; N[3 * i + 1] += ny; N[3 * i + 2] += nz; }
    for (let k = 0; k < 3; k++) {                                                // corner k: cot of its angle weights the opposite edge
      const i = I[k], j = I[(k + 1) % 3], l = I[(k + 2) % 3];
      const e1 = [V[3 * j] - V[3 * i], V[3 * j + 1] - V[3 * i + 1], V[3 * j + 2] - V[3 * i + 2]], e2 = [V[3 * l] - V[3 * i], V[3 * l + 1] - V[3 * i + 1], V[3 * l + 2] - V[3 * i + 2]];
      const cr = Math.hypot(e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]);
      const cot = cr > 1e-300 ? (e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2]) / cr : 0, w = 0.5 * Math.max(-10, Math.min(10, cot));
      for (let d = 0; d < 3; d++) { const diff = V[3 * j + d] - V[3 * l + d]; LP[3 * j + d] += w * diff; LP[3 * l + d] -= w * diff; }
    }
  }
  let conv = 0;
  for (let i = 0; i < n; i++) { const ln = Math.hypot(N[3 * i], N[3 * i + 1], N[3 * i + 2]) || 1; conv += (N[3 * i] * LP[3 * i] + N[3 * i + 1] * LP[3 * i + 1] + N[3 * i + 2] * LP[3 * i + 2]) / ln; }
  const s = conv / bbox;                                                         // ≈ 4π for a sphere, 2π for a tube, ~0 for a saddle
  const flip = closed ? vol < 0 : (Math.abs(s) > 0.3 ? s < 0 : (Math.abs(sz) > 0.05 * tot ? sz < 0 : false));
  if (flip) for (let f = 0; f < nf; f++) { const t = Fo[3 * f + 1]; Fo[3 * f + 1] = Fo[3 * f + 2]; Fo[3 * f + 2] = t; }
  return Fo;
}

export function loopSubdivide(V, F) {
  const n = V.length / 3, nf = F.length / 3;
  const nbrs = Array.from({ length: n }, () => new Set()), edgeFaces = new Map();
  for (let f = 0; f < nf; f++) {
    const a = F[3 * f], b = F[3 * f + 1], c = F[3 * f + 2];
    nbrs[a].add(b).add(c); nbrs[b].add(a).add(c); nbrs[c].add(a).add(b);
    for (const [u, v] of [[a, b], [b, c], [c, a]]) { const k = edgeKey(u, v); let l = edgeFaces.get(k); if (!l) edgeFaces.set(k, l = []); l.push(f); }
  }
  const bnd = new Uint8Array(n);
  for (const [k, l] of edgeFaces) if (l.length === 1) { const a = Math.floor(k / 4294967296), b = k % 4294967296; bnd[a] = 1; bnd[b] = 1; }
  const NV = new Float64Array(V.length);
  for (let i = 0; i < n; i++) {
    const k = nbrs[i].size;
    if (bnd[i] || k === 0) { NV[3 * i] = V[3 * i]; NV[3 * i + 1] = V[3 * i + 1]; NV[3 * i + 2] = V[3 * i + 2]; continue; }
    const t = 3 / 8 + 0.25 * Math.cos(2 * Math.PI / k), beta = (5 / 8 - t * t) / k;
    let sx = 0, sy = 0, sz = 0; for (const j of nbrs[i]) { sx += V[3 * j]; sy += V[3 * j + 1]; sz += V[3 * j + 2]; }
    NV[3 * i] = (1 - k * beta) * V[3 * i] + beta * sx; NV[3 * i + 1] = (1 - k * beta) * V[3 * i + 1] + beta * sy; NV[3 * i + 2] = (1 - k * beta) * V[3 * i + 2] + beta * sz;
  }
  const verts = Array.from(NV), eid = new Map();
  for (const [k, l] of edgeFaces) {
    const u = Math.floor(k / 4294967296), v = k % 4294967296;
    const opp = []; for (const f of l) for (let s = 0; s < 3; s++) { const w = F[3 * f + s]; if (w !== u && w !== v) opp.push(w); }
    let p;
    if (opp.length === 2) p = [0, 1, 2].map(d => 3 / 8 * (V[3 * u + d] + V[3 * v + d]) + 1 / 8 * (V[3 * opp[0] + d] + V[3 * opp[1] + d]));
    else p = [0, 1, 2].map(d => 0.5 * (V[3 * u + d] + V[3 * v + d]));
    eid.set(k, verts.length / 3); verts.push(...p);
  }
  const NF = [];
  for (let f = 0; f < nf; f++) {
    const a = F[3 * f], b = F[3 * f + 1], c = F[3 * f + 2];
    const ab = eid.get(edgeKey(a, b)), bc = eid.get(edgeKey(b, c)), ca = eid.get(edgeKey(c, a));
    NF.push(a, ab, ca, b, bc, ab, c, ca, bc, ab, bc, ca);
  }
  return { V: Float64Array.from(verts), F: Uint32Array.from(NF) };
}

/** Midpoint (linear) subdivision: every triangle → 4 similar ones. The surface is unchanged, old vertex ids are
 *  kept, new vertices (edge midpoints, appended) are flat (defect 0) and lie on the boundary iff their edge does.
 *  Returns the parent edge of every new vertex so base-level fields can be interpolated. */
export function midpointSubdivide(V, F) {
  const n = V.length / 3, nf = F.length / 3, eid = new Map(), count = [], verts = Array.from(V), parents = [];
  const mid = (u, v) => {
    const k = edgeKey(u, v); let id = eid.get(k);
    if (id === undefined) {
      id = verts.length / 3; eid.set(k, id); count.push(0); parents.push(u, v);
      verts.push(0.5 * (V[3 * u] + V[3 * v]), 0.5 * (V[3 * u + 1] + V[3 * v + 1]), 0.5 * (V[3 * u + 2] + V[3 * v + 2]));
    }
    count[id - n]++; return id;
  };
  const NF = new Uint32Array(nf * 12);
  for (let f = 0; f < nf; f++) {
    const a = F[3 * f], b = F[3 * f + 1], c = F[3 * f + 2], ab = mid(a, b), bc = mid(b, c), ca = mid(c, a);
    NF.set([a, ab, ca, b, bc, ab, c, ca, bc, ab, bc, ca], 12 * f);
  }
  const newBoundary = Uint8Array.from(count, (k) => (k === 1 ? 1 : 0));
  return { V: Float64Array.from(verts), F: NF, parents: Uint32Array.from(parents), newBoundary };
}

// ------------------------------------------------------------------ intrinsic Delaunay (edge flips on lengths)
function oppAngle(L, f, s) {
  const a = L[3 * f + s], b = L[3 * f + (s + 1) % 3], c = L[3 * f + (s + 2) % 3];
  return Math.acos(Math.max(-1, Math.min(1, (b * b + c * c - a * a) / Math.max(2 * b * c, 1e-300))));
}

export function flipToDelaunay(V, F) {
  const nf = F.length / 3, n = V.length / 3;
  const Fi = Int32Array.from(F), L = new Float64Array(nf * 3), G = new Int32Array(nf * 6).fill(-1);
  for (let f = 0; f < nf; f++) for (let s = 0; s < 3; s++) {
    const a = Fi[3 * f + s], b = Fi[3 * f + (s + 1) % 3];
    L[3 * f + s] = Math.hypot(V[3 * a] - V[3 * b], V[3 * a + 1] - V[3 * b + 1], V[3 * a + 2] - V[3 * b + 2]);
  }
  const dir = new Map();
  for (let f = 0; f < nf; f++) for (let s = 0; s < 3; s++) dir.set(Fi[3 * f + s] * 4294967296 + Fi[3 * f + (s + 1) % 3], f * 3 + s);
  for (let f = 0; f < nf; f++) for (let s = 0; s < 3; s++) {
    const m = dir.get(Fi[3 * f + (s + 1) % 3] * 4294967296 + Fi[3 * f + s]);
    if (m !== undefined) { G[2 * (3 * f + s)] = Math.floor(m / 3); G[2 * (3 * f + s) + 1] = m % 3; }
  }
  const deg = new Int32Array(n); for (let k = 0; k < Fi.length; k++) deg[Fi[k]]++;
  const isDelaunay = (f, s) => { const g = G[2 * (3 * f + s)], t = G[2 * (3 * f + s) + 1]; return g < 0 || oppAngle(L, f, s) + oppAngle(L, g, t) <= Math.PI + 1e-6; };
  const flip = (f, s) => {
    const g = G[2 * (3 * f + s)], t = G[2 * (3 * f + s) + 1];
    if (g < 0 || g === f) return false;
    const s1 = (s + 1) % 3, s2 = (s + 2) % 3, t1 = (t + 1) % 3, t2 = (t + 2) % 3;
    const v0 = Fi[3 * f + s], v1 = Fi[3 * f + s1], v2 = Fi[3 * f + s2], v3 = Fi[3 * g + t2];
    if (v2 === v3 || v0 === v1 || deg[v0] <= 3 || deg[v1] <= 3) return false;
    const d = L[3 * f + s], l12 = L[3 * f + s1], l20 = L[3 * f + s2], l03 = L[3 * g + t1], l31 = L[3 * g + t2];
    const x2 = (l20 * l20 - l12 * l12 + d * d) / (2 * d), y2 = Math.sqrt(Math.max(l20 * l20 - x2 * x2, 0));
    const x3 = (l03 * l03 - l31 * l31 + d * d) / (2 * d), y3 = -Math.sqrt(Math.max(l03 * l03 - x3 * x3, 0));
    const lnew = Math.hypot(x2 - x3, y2 - y3);
    if (lnew <= 1e-12) return false;
    const nf1 = [G[2 * (3 * f + s1)], G[2 * (3 * f + s1) + 1]], nf2 = [G[2 * (3 * f + s2)], G[2 * (3 * f + s2) + 1]];
    const ng1 = [G[2 * (3 * g + t1)], G[2 * (3 * g + t1) + 1]], ng2 = [G[2 * (3 * g + t2)], G[2 * (3 * g + t2) + 1]];
    Fi[3 * f] = v2; Fi[3 * f + 1] = v0; Fi[3 * f + 2] = v3; L[3 * f] = l20; L[3 * f + 1] = l03; L[3 * f + 2] = lnew;
    Fi[3 * g] = v3; Fi[3 * g + 1] = v1; Fi[3 * g + 2] = v2; L[3 * g] = l31; L[3 * g + 1] = l12; L[3 * g + 2] = lnew;
    const setG = (f_, s_, p) => { G[2 * (3 * f_ + s_)] = p[0]; G[2 * (3 * f_ + s_) + 1] = p[1]; };
    setG(f, 0, nf2); setG(f, 1, ng1); setG(f, 2, [g, 2]); setG(g, 0, ng2); setG(g, 1, nf1); setG(g, 2, [f, 2]);
    for (const [nb, nw] of [[nf2, [f, 0]], [ng1, [f, 1]], [ng2, [g, 0]], [nf1, [g, 1]]]) if (nb[0] >= 0) setG(nb[0], nb[1], nw);
    deg[v0]--; deg[v1]--; deg[v2]++; deg[v3]++;
    return true;
  };
  const queue = [];
  for (let f = 0; f < nf; f++) for (let s = 0; s < 3; s++) if (G[2 * (3 * f + s)] > f) queue.push(f * 3 + s);
  let flips = 0; const maxFlips = 50 * nf + 1000;
  while (queue.length && flips < maxFlips) {
    const q = queue.pop(), f = Math.floor(q / 3), s = q % 3;
    if (isDelaunay(f, s) || !flip(f, s)) continue;
    flips++; const g = G[2 * (3 * f + 2)];
    queue.push(f * 3, f * 3 + 1, g * 3, g * 3 + 1);
  }
  return { Fi: Uint32Array.from(Fi), Li: Float32Array.from(L), flips };
}

// ------------------------------------------------------------------ operators from intrinsic faces
export function cotanLaplacian(Fi, Li, n) {
  const nf = Fi.length / 3, rows = [], cols = [], vals = [], m = new Float64Array(n);
  for (let f = 0; f < nf; f++) {
    const a = Li[3 * f], b = Li[3 * f + 1], c = Li[3 * f + 2], s = 0.5 * (a + b + c);
    const area = Math.sqrt(Math.max(s * (s - a) * (s - b) * (s - c), 0));
    for (let k = 0; k < 3; k++) {
      const x = Li[3 * f + k], y = Li[3 * f + (k + 1) % 3], z = Li[3 * f + (k + 2) % 3];
      let cot = (y * y + z * z - x * x) / Math.max(4 * area, 1e-300); cot = Math.max(-1e6, Math.min(1e6, cot));
      const w = 0.5 * cot, i = Fi[3 * f + k], j = Fi[3 * f + (k + 1) % 3];
      rows.push(i, j, i, j); cols.push(j, i, i, j); vals.push(-w, -w, w, w);
      m[Fi[3 * f + k]] += area / 3;
    }
  }
  return { L: csrFromCOO(n, rows, cols, vals), m };
}

export function angleDefects(V, F, n) {
  const angsum = new Float64Array(n), nf = F.length / 3;
  for (let f = 0; f < nf; f++) for (let k = 0; k < 3; k++) {
    const a = F[3 * f + k], b = F[3 * f + (k + 1) % 3], c = F[3 * f + (k + 2) % 3];
    const ux = V[3 * b] - V[3 * a], uy = V[3 * b + 1] - V[3 * a + 1], uz = V[3 * b + 2] - V[3 * a + 2];
    const wx = V[3 * c] - V[3 * a], wy = V[3 * c + 1] - V[3 * a + 1], wz = V[3 * c + 2] - V[3 * a + 2];
    const d = (ux * wx + uy * wy + uz * wz) / Math.max(Math.hypot(ux, uy, uz) * Math.hypot(wx, wy, wz), 1e-300);
    angsum[a] += Math.acos(Math.max(-1, Math.min(1, d)));
  }
  const { boundary } = edgeTable(F, n);
  const K = new Float64Array(n); for (let i = 0; i < n; i++) K[i] = (boundary[i] ? Math.PI : 2 * Math.PI) - angsum[i];
  return { K, boundary };
}

export function vertexNormals(V, F, n) {
  const N = new Float64Array(3 * n), nf = F.length / 3;
  for (let f = 0; f < nf; f++) {
    const a = F[3 * f], b = F[3 * f + 1], c = F[3 * f + 2];
    const ux = V[3 * b] - V[3 * a], uy = V[3 * b + 1] - V[3 * a + 1], uz = V[3 * b + 2] - V[3 * a + 2];
    const wx = V[3 * c] - V[3 * a], wy = V[3 * c + 1] - V[3 * a + 1], wz = V[3 * c + 2] - V[3 * a + 2];
    const nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx;
    for (const i of [a, b, c]) { N[3 * i] += nx; N[3 * i + 1] += ny; N[3 * i + 2] += nz; }
  }
  for (let i = 0; i < n; i++) { const l = Math.hypot(N[3 * i], N[3 * i + 1], N[3 * i + 2]) || 1; N[3 * i] /= l; N[3 * i + 1] /= l; N[3 * i + 2] /= l; }
  return N;
}

// ------------------------------------------------------------------ the carrier
export class Carrier {
  /** data: {name, meta, V, F, Fi, Li, K, B} (from .lshc or buildFromMesh) */
  constructor(data) {
    this.name = data.name; this.meta = data.meta || {};
    const n = this.n = data.V.length / 3;
    this.V = Float64Array.from(data.V); this.F = Uint32Array.from(data.F);
    this.Fi = Uint32Array.from(data.Fi); this.Li = Float64Array.from(data.Li);
    // side rule applied on load too, so carriers exported with the older +z rule (branching tubes) face outward
    if (!data.sideChecked) this.F = orientSide(this.V, this.F);
    // centre the mesh
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], this.V[3 * i + k]); hi[k] = Math.max(hi[k], this.V[3 * i + k]); }
    const cen = [0, 1, 2].map(k => 0.5 * (lo[k] + hi[k]));
    for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) this.V[3 * i + k] -= cen[k];
    this.lo = lo.map((v, k) => v - cen[k]); this.hi = hi.map((v, k) => v - cen[k]);
    this.bbox = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
    this.K = Float64Array.from(data.K); this.boundary = Uint8Array.from(data.B);
    const { L, m } = cotanLaplacian(this.Fi, this.Li, n);
    this.L = L; this.m = m; this.Ldiag = csrDiagonal(L);
    this.area = m.reduce((s, v) => s + v, 0);
    this.N = vertexNormals(this.V, this.F, n);
    let omega = 0; for (let i = 0; i < n; i++) if (!this.boundary[i]) omega += this.K[i];
    this.omega = omega;
    // Euler characteristic from the extrinsic faces
    this.chi = n - edgeTable(this.F, n).ne + this.F.length / 3;
    // intrinsic edge graph (lengths) for Dijkstra / peaks / Voronoi
    const rows = [], cols = [], vals = [], seen = new Set();
    for (let f = 0; f < this.Fi.length / 3; f++) for (let s = 0; s < 3; s++) {
      const a = this.Fi[3 * f + s], b = this.Fi[3 * f + (s + 1) % 3], k = edgeKey(a, b);
      if (seen.has(k)) continue; seen.add(k);
      rows.push(a, b); cols.push(b, a); vals.push(this.Li[3 * f + s], this.Li[3 * f + s]);
    }
    this.graph = csrFromCOO(n, rows, cols, vals);
    this.h = this.meta.h || medianEdge(this.V, this.F);
    this.hBase = this.meta.hBase || this.h; this.level = this.meta.level || 0;
    this.base = null; this.parents = [];                 // set by refined(): the level-0 carrier and the midpoint parents per level
    this._heatCache = new Map();
  }

  /** A finer copy of this carrier by `levels` midpoint subdivisions: same surface, same vertex ids for the base
   *  vertices, same defects, intrinsic Delaunay rebuilt. Detection, smoothing and the seed chart stay on the base
   *  carrier (fast) and are interpolated, because the surface is identical. Level 0 returns the base itself. */
  refined(levels) {
    const base = this.base || this;
    if (levels <= 0) return base;
    let V = base.V, F = base.F, B = base.boundary; const parents = [];
    for (let l = 0; l < levels; l++) {
      const r = midpointSubdivide(V, F);
      const B1 = new Uint8Array(r.V.length / 3); B1.set(B); B1.set(r.newBoundary, V.length / 3);
      V = r.V; F = r.F; B = B1; parents.push(r.parents);
    }
    const n = V.length / 3, K = new Float64Array(n); K.set(base.K);
    const { Fi, Li, flips } = flipToDelaunay(V, F);
    const c = new Carrier({ name: base.name, meta: { ...base.meta, h: base.h / 2 ** levels, hBase: base.h, level: levels, flips, n, nf: F.length / 3 }, V, F, Fi, Li, K, B, sideChecked: true });
    c.base = base; c.parents = parents;
    return c;
  }

  /** interpolate a base-level per-vertex field to this (refined) carrier: midpoints get the mean of their edge ends */
  interp(vals) {
    if (!this.base) return vals;
    let cur = vals;
    for (const par of this.parents) {
      const n0 = cur.length, out = new Float64Array(n0 + par.length / 2); out.set(cur);
      for (let i = 0; i < par.length / 2; i++) out[n0 + i] = 0.5 * (cur[par[2 * i]] + cur[par[2 * i + 1]]);
      cur = out;
    }
    return cur;
  }

  static async fromURL(url) {
    const res = await fetch(url); if (!res.ok) throw new Error(`fetch ${url}: ${res.status}`);
    return new Carrier(parseLSHC(await res.arrayBuffer()));
  }

  /** Build from a raw triangle soup: weld → orient → Loop-subdivide until the median edge ≤ h → intrinsic Delaunay. */
  static fromMesh(name, V, F, targetH, maxVerts = 40000) {
    let mesh = weld(V, F); mesh.F = orientConsistently(mesh.V, mesh.F);
    let levels = 0;
    while (medianEdge(mesh.V, mesh.F) > targetH && mesh.V.length / 3 * 4 <= maxVerts && levels < 5) { mesh = loopSubdivide(mesh.V, mesh.F); levels++; }
    const n = mesh.V.length / 3;
    const { Fi, Li, flips } = flipToDelaunay(mesh.V, mesh.F);
    const { K, boundary } = angleDefects(mesh.V, mesh.F, n);
    const notes = [`no remesher in the browser: Loop ×${levels} + ${flips} Delaunay flips`];
    return new Carrier({ name, meta: { h: medianEdge(mesh.V, mesh.F), flips, loop: levels, notes, input_vertices: V.length / 3 }, V: mesh.V, F: mesh.F, Fi, Li, K, B: boundary });
  }

  // --- solves on the carrier ------------------------------------------------------------------
  /** (M + t L) x = rhs, Jacobi-PCG */
  heatSolve(rhs, t, x0 = null) {
    const { L, m, n } = this, tmp = new Float64Array(n);
    const apply = (x, out) => { spmv(L, x, tmp); for (let i = 0; i < n; i++) out[i] = m[i] * x[i] + t * tmp[i]; };
    const diag = new Float64Array(n); for (let i = 0; i < n; i++) diag[i] = m[i] + t * this.Ldiag[i];
    return pcg(apply, diag, rhs, x0, { maxIter: 800, tol: 1e-7 }).x;
  }

  smoothedCurvature(radius) {
    if (this.base) return this.interp(this.base.smoothedCurvature(radius));
    const key = radius.toPrecision(3);
    if (!this._heatCache.has(key)) {
      const src = new Float64Array(this.n); for (let i = 0; i < this.n; i++) src[i] = this.boundary[i] ? 0 : this.K[i];
      this._heatCache.set(key, this.heatSolve(src, 0.5 * radius * radius));
    }
    return this._heatCache.get(key);
  }

  /** heat-method geodesic distance from vertex c (Crane–Weischedel–Wardetzky), intrinsic L, extrinsic gradients */
  heatDistance(c) {
    const { n, V, F, L, m } = this;
    const t = 4 * this.h * this.h;
    const delta = new Float64Array(n); delta[c] = 1;
    const u = this.heatSolve(delta, t);
    const nf = F.length / 3, div = new Float64Array(n);
    for (let f = 0; f < nf; f++) {
      const a = F[3 * f], b = F[3 * f + 1], cc = F[3 * f + 2];
      const pa = [V[3 * a], V[3 * a + 1], V[3 * a + 2]], pb = [V[3 * b], V[3 * b + 1], V[3 * b + 2]], pc = [V[3 * cc], V[3 * cc + 1], V[3 * cc + 2]];
      const e0 = sub(pc, pb), e1 = sub(pa, pc), e2 = sub(pb, pa);           // edge opposite a, b, c
      const nrm = cross(e2, sub(pc, pa)); const A2 = norm(nrm); if (A2 < 1e-300) continue;
      const nn = nrm.map(v => v / A2);
      // ∇u = (1/2A) Σ u_i (N × e_i)
      let g = [0, 0, 0];
      for (const [ui, e] of [[u[a], e0], [u[b], e1], [u[cc], e2]]) { const ce = cross(nn, e); g = [g[0] + ui * ce[0], g[1] + ui * ce[1], g[2] + ui * ce[2]]; }
      const gl = norm(g); if (gl < 1e-300) continue;
      const X = g.map(v => -v / gl);
      // divergence: Σ over corners ½ (cot θ1 (e1·X) + cot θ2 (e2·X))
      const corners = [[a, pa, pb, pc], [b, pb, pc, pa], [cc, pc, pa, pb]];
      for (const [i, pi, pj, pk] of corners) {
        const eij = sub(pj, pi), eik = sub(pk, pi);
        const cotk = cotAt(pk, pi, pj), cotj = cotAt(pj, pi, pk);          // angle at k opposite edge ij, angle at j opposite edge ik
        div[i] += 0.5 * (cotk * dotv(eij, X) + cotj * dotv(eik, X));
      }
    }
    // Poisson L φ = div (L psd, singular): zero-mean right-hand side, zero-mean iterates
    let mean = 0; for (let i = 0; i < n; i++) mean += div[i]; mean /= n; for (let i = 0; i < n; i++) div[i] -= mean;
    const diag = new Float64Array(n); for (let i = 0; i < n; i++) diag[i] = Math.max(this.Ldiag[i], 1e-12);
    const proj = (v) => { let s = 0; for (let i = 0; i < n; i++) s += v[i]; s /= n; for (let i = 0; i < n; i++) v[i] -= s; };
    const phi = pcg((x, out) => spmv(L, x, out), diag, div, null, { maxIter: 1500, tol: 1e-7, projectMean: proj }).x;
    const p0 = phi[c]; const d = new Float64Array(n); let msum = 0;
    for (let i = 0; i < n; i++) { d[i] = phi[i] - p0; msum += d[i]; }
    if (msum < 0) for (let i = 0; i < n; i++) d[i] = -d[i];                // sign convention of div vs. positive L
    for (let i = 0; i < n; i++) d[i] = Math.max(d[i], 0);
    return d;
  }

  /** Geodesic distance from vertex c by fast marching on the intrinsic triangulation (Dijkstra order, triangle
   *  updates with the unfolded virtual source; edge update where the source is not visible through the opposite edge).
   *  Exact on a developable patch (cone) seen from its apex; unlike the heat method it keeps the distance right up
   *  to the rim (the heat method's Neumann condition compressed the v4 cone's rim from 6.1 a to 5.0 a). Cached. */
  geodesicDistance(c) {
    if (this.base) return this.interp(this.base.geodesicDistance(c));
    this._geoCache = this._geoCache || new Map();
    const key = Array.isArray(c) ? 'set:' + c.length + ':' + c[0] : c;      // c: one vertex, or an array of sources (the rim)
    if (this._geoCache.has(key)) return this._geoCache.get(key);
    const { n, Fi, Li } = this, nf = Fi.length / 3;
    const vf = Array.from({ length: n }, () => []);
    for (let f = 0; f < nf; f++) for (let s = 0; s < 3; s++) vf[Fi[3 * f + s]].push(3 * f + s);
    const d = new Float64Array(n).fill(Infinity), done = new Uint8Array(n);
    const heap = [], push = (v, key) => { heap.push([key, v]); let i = heap.length - 1; while (i > 0) { const p = (i - 1) >> 1; if (heap[p][0] <= heap[i][0]) break; [heap[p], heap[i]] = [heap[i], heap[p]]; i = p; } };
    const pop = () => { const top = heap[0], last = heap.pop(); if (heap.length) { heap[0] = last; let i = 0; for (;;) { const l = 2 * i + 1, r = l + 1; let m = i; if (l < heap.length && heap[l][0] < heap[m][0]) m = l; if (r < heap.length && heap[r][0] < heap[m][0]) m = r; if (m === i) break; [heap[m], heap[i]] = [heap[i], heap[m]]; i = m; } } return top; };
    // candidate distance at k from known i, j in triangle with edge lengths lij, lik, ljk
    const update = (di, dj, lij, lik, ljk) => {
      let best = Math.min(di + lik, dj + ljk);
      if (!Number.isFinite(di) || !Number.isFinite(dj)) return best;
      const xk = (lik * lik - ljk * ljk + lij * lij) / (2 * lij), yk = Math.sqrt(Math.max(lik * lik - xk * xk, 0));
      const xs = (di * di - dj * dj + lij * lij) / (2 * lij), ys2 = di * di - xs * xs;
      if (ys2 < 0) return best;
      const ys = -Math.sqrt(ys2);                                             // virtual source on the far side of edge ij
      const t = -ys / (yk - ys), xc = xs + t * (xk - xs);                     // where the segment source → k crosses y = 0
      if (xc >= 0 && xc <= lij) best = Math.min(best, Math.hypot(xk - xs, yk - ys));
      return best;
    };
    for (const src of (Array.isArray(c) ? c : [c])) { d[src] = 0; push(src, 0); }
    while (heap.length) {
      const [key, v] = pop(); if (done[v] || key > d[v]) continue; done[v] = 1;
      for (const corner of vf[v]) {
        const f = (corner / 3) | 0, s = corner % 3;
        const v1 = Fi[3 * f + (s + 1) % 3], v2 = Fi[3 * f + (s + 2) % 3];
        const l01 = Li[3 * f + s], l12 = Li[3 * f + (s + 1) % 3], l20 = Li[3 * f + (s + 2) % 3];   // v→v1, v1→v2, v2→v
        for (const [k, j, lvk, lvj, ljk] of [[v1, v2, l01, l20, l12], [v2, v1, l20, l01, l12]]) {
          if (done[k]) continue;
          const cand = done[j] ? update(d[v], d[j], lvj, lvk, ljk) : d[v] + lvk;
          if (cand < d[k] - 1e-12) { d[k] = cand; push(k, cand); }
        }
      }
    }
    this._geoCache.set(key, d);
    return d;
  }
  /** The intrinsic Delaunay triangulation drawn on the surface. Edges shared with the extrinsic mesh are straight;
   *  every flipped edge is a geodesic across several extrinsic faces, traced from one end: its direction at that
   *  vertex is read off the intrinsic corner angles (measured from an edge both triangulations share), then the
   *  straight line is unfolded face by face for the edge's intrinsic length.
   *  Returns { shared: Float32Array (segment pairs), flipped: Float32Array (segment pairs), nShared, nFlipped, failed }. */
  intrinsicEdgeLines(lift = 0) {
    if (this._iel && this._iel.lift === lift) return this._iel.out;
    const { n, V, F, Fi, Li, N } = this, nf = F.length / 3, nfi = Fi.length / 3;
    const ext = new Set(), adj = new Map();
    for (let f = 0; f < nf; f++) for (let s = 0; s < 3; s++) {
      const a = F[3 * f + s], b = F[3 * f + (s + 1) % 3], k = edgeKey(a, b); ext.add(k);
      let l = adj.get(k); if (!l) adj.set(k, l = []); l.push(f);
    }
    const P = (i) => [V[3 * i], V[3 * i + 1], V[3 * i + 2]];
    const dist3 = (i, j) => Math.hypot(V[3 * i] - V[3 * j], V[3 * i + 1] - V[3 * j + 1], V[3 * i + 2] - V[3 * j + 2]);
    const lifted = (p, nrm) => { const l = Math.hypot(nrm[0], nrm[1], nrm[2]) || 1; return [p[0] + lift * nrm[0] / l, p[1] + lift * nrm[1] / l, p[2] + lift * nrm[2] / l]; };
    // flipped edges first, so the vertex fans below are only built where a trace starts
    const need = new Set();
    for (let f = 0; f < nfi; f++) for (let s = 0; s < 3; s++) { const i = Fi[3 * f + s], j = Fi[3 * f + (s + 1) % 3]; if (!ext.has(edgeKey(i, j))) { need.add(i); need.add(j); } }
    // fan around a vertex: outgoing neighbours in ccw order with cumulative angles, for a triangulation
    const fans = (T, angleAt) => {
      const nT = T.length / 3, next = new Map(), isB = new Map();
      for (let f = 0; f < nT; f++) for (let s = 0; s < 3; s++) {
        const i = T[3 * f + s]; if (!need.has(i)) continue;
        const a = T[3 * f + (s + 1) % 3], b = T[3 * f + (s + 2) % 3];
        if (!next.has(i)) { next.set(i, new Map()); isB.set(i, new Set()); }
        next.get(i).set(a, { b, f, s, ang: angleAt(f, s) }); isB.get(i).add(b);
      }
      return (i) => {
        const nx = next.get(i); if (!nx || !nx.size) return null;
        const ib = isB.get(i);
        let start = nx.keys().next().value;
        for (const a of nx.keys()) if (!ib.has(a)) { start = a; break; }       // boundary vertex: begin at the rim edge
        const nb = [], cum = [], face = [], ang = []; let a = start, total = 0, guard = 0;
        while (guard++ <= nx.size) {
          const e = nx.get(a); nb.push(a); cum.push(total); if (!e) { face.push(-1); ang.push(0); break; }
          face.push(e.f); ang.push(e.ang); total += e.ang; a = e.b; if (a === start) break;
        }
        return { nb, cum, face, ang, total, closed: a === start };
      };
    };
    const extAngle = (f, s) => { const i = F[3 * f + s], a = F[3 * f + (s + 1) % 3], b = F[3 * f + (s + 2) % 3];
      const u = [V[3 * a] - V[3 * i], V[3 * a + 1] - V[3 * i + 1], V[3 * a + 2] - V[3 * i + 2]], w = [V[3 * b] - V[3 * i], V[3 * b + 1] - V[3 * i + 1], V[3 * b + 2] - V[3 * i + 2]];
      return Math.atan2(Math.hypot(u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]), u[0] * w[0] + u[1] * w[1] + u[2] * w[2]); };
    const intAngle = (f, s) => { const x = Li[3 * f + s], y = Li[3 * f + (s + 2) % 3], z = Li[3 * f + (s + 1) % 3];
      return Math.acos(Math.max(-1, Math.min(1, (x * x + y * y - z * z) / (2 * x * y)))); };
    const extFan = fans(F, extAngle), intFan = fans(Fi, intAngle);
    const shared = new Float32Array(nfi * 3 * 6), flipped = []; let nShared = 0, nFlipped = 0, failed = 0;
    const putShared = (i, j) => { const o = 6 * nShared++, sc = lift; for (let d = 0; d < 3; d++) { shared[o + d] = V[3 * i + d] + sc * N[3 * i + d]; shared[o + 3 + d] = V[3 * j + d] + sc * N[3 * j + d]; } };
    const seen = new Set();
    const cross2 = (a, b) => a[0] * b[1] - a[1] * b[0];
    for (let f = 0; f < nfi; f++) for (let s = 0; s < 3; s++) {
      const i0 = Fi[3 * f + s], j0 = Fi[3 * f + (s + 1) % 3], k = edgeKey(i0, j0);
      if (seen.has(k)) continue; seen.add(k);
      if (ext.has(k)) { putShared(i0, j0); continue; }
      const len = Li[3 * f + s];
      let pts = null;
      for (const [i, j] of [[i0, j0], [j0, i0]]) {                            // trace from one end, else from the other
        pts = traceFrom(i, j, len); if (pts) break;
      }
      if (!pts) { failed++; pts = [lifted(P(i0), [N[3 * i0], N[3 * i0 + 1], N[3 * i0 + 2]]), lifted(P(j0), [N[3 * j0], N[3 * j0 + 1], N[3 * j0 + 2]])]; }
      for (let q = 0; q + 1 < pts.length; q++) flipped.push(...pts[q], ...pts[q + 1]);
      nFlipped++;
    }
    function traceFrom(i, j, len) {
      const fi = intFan(i), fe = extFan(i); if (!fi || !fe) return null;
      const ti = fi.nb.indexOf(j); if (ti < 0) return null;
      let ref = -1; for (const w of fi.nb) if (fe.nb.includes(w)) { ref = w; break; }
      if (ref < 0) return null;
      let th = fe.cum[fe.nb.indexOf(ref)] + (fi.cum[ti] - fi.cum[fi.nb.indexOf(ref)]);
      if (fe.closed) { th %= fe.total; if (th < 0) th += fe.total; }
      // extrinsic sector containing th
      let k = -1; for (let q = 0; q < fe.face.length; q++) if (fe.face[q] >= 0 && th >= fe.cum[q] - 1e-12 && th <= fe.cum[q] + fe.ang[q] + 1e-12) { k = q; break; }
      if (k < 0) return null;
      // 2-D frame of face (i, a, b): i at 0, a on +x
      const f0 = fe.face[k]; let sI = 0; while (F[3 * f0 + sI] !== i) sI++;
      const va = F[3 * f0 + (sI + 1) % 3], vb = F[3 * f0 + (sI + 2) % 3];
      const lia = dist3(i, va), lib = dist3(i, vb), al = fe.ang[k];
      const pos = new Map([[i, [0, 0]], [va, [lia, 0]], [vb, [lib * Math.cos(al), lib * Math.sin(al)]]]);
      const psi = th - fe.cum[k], d = [Math.cos(psi), Math.sin(psi)];
      const out = [lifted(P(i), [N[3 * i], N[3 * i + 1], N[3 * i + 2]])];
      let face = f0, enter = [i, i], tPrev = 0;
      for (let guard = 0; guard < 400; guard++) {
        const tri = [F[3 * face], F[3 * face + 1], F[3 * face + 2]];
        if (tri.includes(j) && j !== i) {
          const pj = pos.get(j); const tj = Math.hypot(pj[0], pj[1]);
          if (Math.abs(tj - len) < 0.05 * len + 1e-9 && Math.abs(cross2(d, pj)) < 0.05 * len + 1e-9) { out.push(lifted(P(j), [N[3 * j], N[3 * j + 1], N[3 * j + 2]])); return out; }
        }
        let best = null;
        for (let e = 0; e < 3; e++) {
          const p = tri[e], q = tri[(e + 1) % 3];
          if ((p === enter[0] && q === enter[1]) || (p === enter[1] && q === enter[0])) continue;
          if (p === i || q === i) { if (guard === 0) continue; }
          const A = pos.get(p), B = pos.get(q), E = [B[0] - A[0], B[1] - A[1]];
          const den = cross2(d, E); if (Math.abs(den) < 1e-15) continue;
          const t = cross2(A, E) / den, u = cross2(A, d) / den;
          if (u < -1e-9 || u > 1 + 1e-9 || t <= tPrev + 1e-12) continue;
          if (!best || t < best.t) best = { t, u: Math.max(0, Math.min(1, u)), p, q };
        }
        if (!best || best.t > len * 1.05) return null;
        const pp = P(best.p), qq = P(best.q), u = best.u;
        const nrm = [(1 - u) * N[3 * best.p] + u * N[3 * best.q], (1 - u) * N[3 * best.p + 1] + u * N[3 * best.q + 1], (1 - u) * N[3 * best.p + 2] + u * N[3 * best.q + 2]];
        out.push(lifted([(1 - u) * pp[0] + u * qq[0], (1 - u) * pp[1] + u * qq[1], (1 - u) * pp[2] + u * qq[2]], nrm));
        // unfold the face across (p, q)
        const fl = adj.get(edgeKey(best.p, best.q)); const nfce = fl && fl.length === 2 ? (fl[0] === face ? fl[1] : fl[0]) : -1;
        if (nfce < 0) return null;
        const third = [F[3 * nfce], F[3 * nfce + 1], F[3 * nfce + 2]].find(v => v !== best.p && v !== best.q);
        {
          const A = pos.get(best.p), B = pos.get(best.q), lAB = Math.hypot(B[0] - A[0], B[1] - A[1]);
          const lAr = dist3(best.p, third), lBr = dist3(best.q, third);
          const x = (lAr * lAr - lBr * lBr + lAB * lAB) / (2 * lAB), y = Math.sqrt(Math.max(lAr * lAr - x * x, 0));
          const ex = [(B[0] - A[0]) / lAB, (B[1] - A[1]) / lAB], ey = [-ex[1], ex[0]];
          const old = tri.find(v => v !== best.p && v !== best.q), O = pos.get(old);
          const side = Math.sign((O[0] - A[0]) * ey[0] + (O[1] - A[1]) * ey[1]) || 1;     // the new vertex goes on the other side
          pos.set(third, [A[0] + x * ex[0] - side * y * ey[0], A[1] + x * ex[1] - side * y * ey[1]]);
        }
        enter = [best.p, best.q]; tPrev = best.t; face = nfce;
      }
      return null;
    }
    const out = { shared: shared.slice(0, 6 * nShared), flipped: Float32Array.from(flipped), nShared, nFlipped, failed };
    this._iel = { lift, out };
    return out;
  }

  /** The seed's first ring on the surface, as a scale bar: the centre, its n neighbours at lattice distance a (in the
   *  folded chart of the seed) and the geodesic spokes to them, plus the circle ρ = a. For the "lattice edge" probe. */
  seedProbe(c, valence, a, phaseDeg = 0, kind = 'wave') {
    const { px, py, rho } = this.logMap(c), { n, V, F } = this, nf = F.length / 3;
    const cand = [];
    for (let f = 0; f < nf; f++) {
      const i = F[3 * f], j = F[3 * f + 1], k = F[3 * f + 2];
      if (rho[i] > 1.4 * a && rho[j] > 1.4 * a && rho[k] > 1.4 * a) continue;
      const ax = px[j] - px[i], ay = py[j] - py[i], bx = px[k] - px[i], by = py[k] - py[i];
      const ca = ax * by - ay * bx;                                          // chart area (signed, doubled)
      const ux = V[3 * j] - V[3 * i], uy = V[3 * j + 1] - V[3 * i + 1], uz = V[3 * j + 2] - V[3 * i + 2];
      const wx = V[3 * k] - V[3 * i], wy = V[3 * k + 1] - V[3 * i + 1], wz = V[3 * k + 2] - V[3 * i + 2];
      const sa = Math.hypot(uy * wz - uz * wy, uz * wx - ux * wz, ux * wy - uy * wx);
      if (!(Math.abs(ca) > 0.2 * sa && Math.abs(ca) < 5 * sa)) continue;    // the branch cut of the chart: skip
      cand.push(f);
    }
    const locate = (X, Y) => {
      for (const f of cand) {
        const i = F[3 * f], j = F[3 * f + 1], k = F[3 * f + 2];
        const ax = px[j] - px[i], ay = py[j] - py[i], bx = px[k] - px[i], by = py[k] - py[i], qx = X - px[i], qy = Y - py[i];
        const den = ax * by - ay * bx; if (Math.abs(den) < 1e-300) continue;
        const u = (qx * by - qy * bx) / den, v = (ax * qy - ay * qx) / den;
        if (u >= -1e-6 && v >= -1e-6 && u + v <= 1 + 1e-6) {
          const w = 1 - u - v;
          return [0, 1, 2].map(d => w * V[3 * i + d] + u * V[3 * j + d] + v * V[3 * k + d]).concat([f]);
        }
      }
      return null;
    };
    const ph0 = phaseDeg * Math.PI / 180, sites = [], spokes = [];
    for (let k = 0; k < valence; k++) {
      // wave seed: lattice points of the folded chart ζ = ρ e^{iθ n/6} at |ζ| = a, arg ζ = k·60° → θ = (k·60° − φ0)·6/n
      const th = kind === 'spots' ? ph0 + 2 * Math.PI * k / valence : (k * Math.PI / 3 - ph0) * 6 / valence;
      const poly = [];
      for (let t = 0; t <= 16; t++) { const r = a * t / 16, p = locate(r * Math.cos(th), r * Math.sin(th)); if (p) poly.push(p.slice(0, 3)); }
      if (poly.length) { spokes.push(poly); sites.push(poly[poly.length - 1]); }
    }
    // ρ = a isoline (marching triangles)
    const ring = [];
    for (let f = 0; f < nf; f++) {
      const ids = [F[3 * f], F[3 * f + 1], F[3 * f + 2]], r = ids.map(i => rho[i] - a), pts = [];
      for (let e = 0; e < 3; e++) { const p = ids[e], q = ids[(e + 1) % 3], rp = r[e], rq = r[(e + 1) % 3];
        if ((rp < 0) !== (rq < 0)) { const t = rp / (rp - rq); pts.push([0, 1, 2].map(d => V[3 * p + d] + t * (V[3 * q + d] - V[3 * p + d]))); } }
      if (pts.length === 2) ring.push(pts[0], pts[1]);
    }
    return { sites, spokes, ring, centre: [V[3 * c], V[3 * c + 1], V[3 * c + 2]] };
  }

  /** geodesic distance to the rim (Infinity everywhere on a closed surface) */
  rimDistance() {
    if (this.base) return this.interp(this.base.rimDistance());
    const src = []; for (let i = 0; i < this.n; i++) if (this.boundary[i]) src.push(i);
    return src.length ? this.geodesicDistance(src) : new Float64Array(this.n).fill(Infinity);
  }

  interiorCentre() {
    if (this.base) return this.base.interiorCentre();                         // base vertex ids are valid here
    const src = []; for (let i = 0; i < this.n; i++) if (this.boundary[i]) src.push(i);
    if (src.length) { const { dist } = dijkstra(this.graph, src); let best = 0; for (let i = 1; i < this.n; i++) if (dist[i] > dist[best]) best = i; return best; }
    let best = 0, bd = Infinity; for (let i = 0; i < this.n; i++) { const d = Math.hypot(this.V[3 * i], this.V[3 * i + 1], this.V[3 * i + 2]); if (d < bd) { bd = d; best = i; } }
    return best;
  }

  /** one ideal singularity: {valence, centre, info} */
  detectSingularity(a, flatDeg = 15) {
    if (this.base) return this.base.detectSingularity(a, flatDeg);           // same surface, same ids, 4^level cheaper
    const omegaDeg = this.omega * 180 / Math.PI, charge = Math.round(this.omega / (Math.PI / 3));
    const info = { omegaDeg, charge };
    let n, c;
    if (Math.abs(omegaDeg) < flatDeg) { c = this.interiorCentre(); n = 6; info.mode = 'flat: seed at the point farthest from the boundary'; }
    else {
      n = Math.max(2, Math.min(12, 6 - charge));
      const w = this.smoothedCurvature(Math.max(a, 2 * this.h)), sg = Math.sign(this.omega);
      c = -1; let best = -Infinity;
      for (let i = 0; i < this.n; i++) if (!this.boundary[i] && sg * w[i] > best) { best = sg * w[i]; c = i; }
      info.mode = 'curvature maximum';
    }
    const { dist } = dijkstra(this.graph, [c]);
    let inside = 0; for (let i = 0; i < this.n; i++) if (dist[i] <= a && !this.boundary[i]) inside += this.K[i];
    info.withinA = Math.abs(this.omega) > 1e-9 ? 100 * inside / this.omega : 0;
    info.defectAtCentre = this.K[c] * 180 / Math.PI;
    return { valence: n, centre: c, info };
  }

  /** per-vertex 1-ring in face order with cumulative corner angles, rescaled to 2π at interior vertices
   *  (the intrinsic tangent space of geometry-central). Built once. */
  _rings() {
    if (this._ring) return this._ring;
    const { n, V, F } = this, nf = F.length / 3;
    const outA = Array.from({ length: n }, () => []), outB = Array.from({ length: n }, () => []), outAng = Array.from({ length: n }, () => []);
    const ang = (v, a, b) => {
      const ux = V[3 * a] - V[3 * v], uy = V[3 * a + 1] - V[3 * v + 1], uz = V[3 * a + 2] - V[3 * v + 2];
      const wx = V[3 * b] - V[3 * v], wy = V[3 * b + 1] - V[3 * v + 1], wz = V[3 * b + 2] - V[3 * v + 2];
      return Math.acos(Math.max(-1, Math.min(1, (ux * wx + uy * wy + uz * wz) / Math.max(Math.hypot(ux, uy, uz) * Math.hypot(wx, wy, wz), 1e-300))));
    };
    for (let f = 0; f < nf; f++) for (let k = 0; k < 3; k++) {
      const v = F[3 * f + k], a = F[3 * f + (k + 1) % 3], b = F[3 * f + (k + 2) % 3];
      outA[v].push(a); outB[v].push(b); outAng[v].push(ang(v, a, b));
    }
    const ring = new Array(n);
    for (let v = 0; v < n; v++) {
      const A = outA[v], B = outB[v], G = outAng[v], m = A.length;
      if (!m) { ring[v] = { nb: [], cum: [], total: 0 }; continue; }
      let start = 0; const bset = new Set(B);
      for (let i = 0; i < m; i++) if (!bset.has(A[i])) { start = i; break; }     // boundary: begin at the open side
      const nb = [A[start]], cum = [0]; let cur = start, total = 0, guard = 0;
      while (guard++ <= m) {
        total += G[cur]; const nxt = B[cur]; const i2 = A.indexOf(nxt);
        if (i2 < 0 || nxt === A[start]) { if (i2 >= 0 && nxt === A[start]) break; nb.push(nxt); cum.push(total); break; }
        nb.push(nxt); cum.push(total); cur = i2;
      }
      const closed = !this.boundary[v] && nb.length === m;
      const scale = closed && total > 1e-9 ? 2 * Math.PI / total : 1;
      ring[v] = { nb, cum: cum.map(x => x * scale), total: total * scale };
    }
    this._ring = ring; return ring;
  }

  /** geodesic polar chart around c: ρ from the heat method; θ by the discrete exponential map propagated upwind
   *  on ρ with intrinsic corner angles (exact on a cone from its apex, where the angles are rescaled to 2π). */
  logMap(c) {
    if (this.base) {                                                          // chart on the base carrier, interpolated (it is smooth at the base scale)
      const b = this.base.logMap(c), rho = this.interp(b.rho), px = this.interp(b.px), py = this.interp(b.py);
      const theta = new Float64Array(this.n); for (let i = 0; i < this.n; i++) theta[i] = Math.atan2(py[i], px[i]);
      return { rho, theta, px, py };
    }
    if (this._logCache && this._logCache.c === c) return this._logCache.map;
    const map = this._logMapCompute(c); this._logCache = { c, map }; return map;
  }
  _logMapCompute(c) {
    const rho = this.heatDistance(c), n = this.n, ring = this._rings(), V = this.V;
    const px = new Float64Array(n), py = new Float64Array(n), rot = new Float64Array(n), done = new Uint8Array(n);
    const angIn = (v, w) => { const r = ring[v], i = r.nb.indexOf(w); return i < 0 ? NaN : r.cum[i]; };
    const order = Array.from({ length: n }, (_, i) => i).sort((i, j) => rho[i] - rho[j]);
    done[c] = 1; rot[c] = 0;
    for (const j of order) {
      if (j === c) continue;
      // parent: processed neighbour with the smallest ρ
      let p = -1, best = Infinity; const G = this.graph;
      for (let k = G.indptr[j]; k < G.indptr[j + 1]; k++) { const q = G.indices[k]; if (done[q] && rho[q] < best && Number.isFinite(angIn(q, j))) { best = rho[q]; p = q; } }
      if (p < 0) { done[j] = 1; continue; }
      const len = Math.hypot(V[3 * j] - V[3 * p], V[3 * j + 1] - V[3 * p + 1], V[3 * j + 2] - V[3 * p + 2]);
      const dir = rot[p] + angIn(p, j);
      px[j] = px[p] + len * Math.cos(dir); py[j] = py[p] + len * Math.sin(dir);
      const back = angIn(j, p);
      rot[j] = Number.isFinite(back) ? Math.atan2(py[p] - py[j], px[p] - px[j]) - back : 0;
      done[j] = 1;
    }
    const theta = new Float64Array(n);
    for (let i = 0; i < n; i++) theta[i] = Math.atan2(py[i], px[i]);
    return { rho, theta, px, py };
  }

  // --- conformal size field --------------------------------------------------------------------------------
  /** Size field φ with  L φ = Σ_c Ω_c δ_c − K  on the interior (K = angle defects, Ω_c = (6 − n_c)·60° at the cone
   *  vertices): ĝ = e^{2φ} g is flat except at the cones, so a lattice that is regular in ĝ has edge a·e^{−φ} in the
   *  model and needs no curvature-induced dislocations. Open surface: φ = 0 on the rim (the rim absorbs the
   *  residual Σ Ω_c − ∫K). Closed surface: the residual is spread as uniform curvature (Gauss–Bonnet needs it to be
   *  0) and φ has zero mean. φ is then shifted so that its mean on the ring ρ ≈ a around `centre` is 0 (edge a there).
   *  cones: [{vertex, valence}]. Returns {phi, residualDeg, min, max, iters}. Computed on the base carrier. */
  sizeField(cones, a, centre) {
    if (this.base) { const r = this.base.sizeField(cones, a, centre); return { ...r, phi: this.interp(r.phi) }; }
    const { n, L, K, boundary, m } = this;
    const rhs = new Float64Array(n);
    for (let i = 0; i < n; i++) rhs[i] = boundary[i] ? 0 : -K[i];
    const open = boundary.some(b => b);
    let sumOmega = 0, skipped = 0;
    for (const c of cones) {
      if (open && boundary[c.vertex]) { skipped++; continue; }              // a cone on the rim is the rim's business
      const om = (6 - c.valence) * Math.PI / 3; rhs[c.vertex] += om; sumOmega += om;
    }
    let residual = sumOmega - this.omega;                                 // what the cones do not account for
    let phi, iters;
    if (open) {
      // Dirichlet φ = 0 on the rim: boundary rows are the identity, their unknowns stay 0 (rhs 0, zero start)
      const tmp = new Float64Array(n);
      const apply = (x, out) => { spmv(L, x, tmp); for (let i = 0; i < n; i++) out[i] = boundary[i] ? x[i] : tmp[i]; };
      const diag = new Float64Array(n); for (let i = 0; i < n; i++) diag[i] = boundary[i] ? 1 : Math.max(this.Ldiag[i], 1e-12);
      const res = pcg(apply, diag, rhs, null, { maxIter: 3000, tol: 1e-8 }); phi = res.x; iters = res.iters;
    } else {
      for (let i = 0; i < n; i++) rhs[i] -= residual * m[i] / this.area;    // make Σ rhs = 0 (L is singular on a closed surface)
      const diag = new Float64Array(n); for (let i = 0; i < n; i++) diag[i] = Math.max(this.Ldiag[i], 1e-12);
      const proj = (v) => { let s = 0; for (let i = 0; i < n; i++) s += v[i]; s /= n; for (let i = 0; i < n; i++) v[i] -= s; };
      const res = pcg((x, out) => spmv(L, x, out), diag, rhs, null, { maxIter: 3000, tol: 1e-8, projectMean: proj }); phi = res.x; iters = res.iters;
    }
    // normalise: mean 0 on the ring ρ ∈ [0.8 a, 1.2 a] around the centre (chart ρ = heat distance)
    if (centre >= 0 && a > 0) {
      const { rho } = this.logMap(centre); let s = 0, k = 0;
      for (let i = 0; i < n; i++) if (rho[i] > 0.8 * a && rho[i] < 1.2 * a) { s += phi[i]; k++; }
      if (k) { s /= k; for (let i = 0; i < n; i++) phi[i] -= s; }
    }
    let lo = Infinity, hi = -Infinity; for (let i = 0; i < n; i++) { if (phi[i] < lo) lo = phi[i]; if (phi[i] > hi) hi = phi[i]; }
    return { phi, residualDeg: residual * 180 / Math.PI, open, min: lo, max: hi, iters, skipped };
  }

  /** Operators of the conformally rescaled metric ĝ = e^{2φ} g: intrinsic lengths ℓ̂ = e^{(φi+φj)/2} ℓ (discrete
   *  conformal equivalence), cotan Laplacian and lumped mass rebuilt on the same intrinsic triangulation. What the
   *  solvers read: {n, L, m, Ldiag, area}; plus the resolution check in ĝ (edges per lattice edge a). */
  conformalOps(phi, a) {
    const { Fi, Li, n } = this, nf = Fi.length / 3, Lh = new Float64Array(Li.length);
    let under = 0, bad = 0, maxE = 0;
    for (let f = 0; f < nf; f++) {
      for (let s = 0; s < 3; s++) {
        const i = Fi[3 * f + s], j = Fi[3 * f + (s + 1) % 3];
        const l = Li[3 * f + s] * Math.exp(0.5 * (phi[i] + phi[j])); Lh[3 * f + s] = l;
        if (l > a / 6) under++; if (l > maxE) maxE = l;
      }
      const x = Lh[3 * f], y = Lh[3 * f + 1], z = Lh[3 * f + 2]; if (x + y <= z || y + z <= x || z + x <= y) bad++;
    }
    const { L, m } = cotanLaplacian(Fi, Lh, n);
    let area = 0; for (let i = 0; i < n; i++) area += m[i];
    return { n, L, m, Ldiag: csrDiagonal(L), area, underResolvedPct: 100 * under / (3 * nf), badTriangles: bad, maxEdgeInGhat: maxE, phi };
  }

  /** ρ̂ = ∫₀^ρ e^{φ} dρ along the chart ray to each vertex within `limit` of c (the geodesic radius in ĝ, for a
   *  seed that is regular in ĝ); φ is sampled by nearest neighbour in chart coordinates. Outside `limit` ρ̂ = ρ e^{φ}. */
  _rhoHat(c, phi, limit, samples = 24) {
    const { rho, px, py } = this.logMap(c), n = this.n, out = new Float64Array(n);
    const cell = Math.max(this.h, 1e-9), grid = new Map(), key = (x, y) => Math.floor(x / cell) * 1048576 + Math.floor(y / cell);
    for (let i = 0; i < n; i++) if (rho[i] < 1.3 * limit) { const k = key(px[i], py[i]); let l = grid.get(k); if (!l) { l = []; grid.set(k, l); } l.push(i); }
    const nearest = (x, y) => {
      let best = c, bd = Infinity; const gx = Math.floor(x / cell), gy = Math.floor(y / cell);
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) { const l = grid.get((gx + dx) * 1048576 + gy + dy); if (!l) continue; for (const i of l) { const d = (px[i] - x) ** 2 + (py[i] - y) ** 2; if (d < bd) { bd = d; best = i; } } }
      return best;
    };
    for (let i = 0; i < n; i++) {
      if (!(rho[i] < limit)) { out[i] = rho[i] * Math.exp(phi[i]); continue; }
      let s = 0; for (let k = 0; k < samples; k++) { const t = (k + 0.5) / samples; s += Math.exp(phi[nearest(t * px[i], t * py[i])]); }
      out[i] = rho[i] * s / samples;
    }
    return out;
  }

  seedField(c, valence, a, { radius = 2.0, amplitude = 1.0, kind = 'wave', phaseDeg = 0, phi = null } = {}) {
    const { theta } = this.logMap(c), u = new Float64Array(this.n), q0 = 4 * Math.PI / (SQ3 * a);
    const rho = phi ? this._rhoHat(c, phi, 1.6 * radius * a) : this.logMap(c).rho;   // with a size field the seed is regular in ĝ
    const ph0 = phaseDeg * Math.PI / 180;
    if (kind === 'spots') {
      const sig = 0.28 * a;
      for (let i = 0; i < this.n; i++) {
        const x = rho[i] * Math.cos(theta[i]), y = rho[i] * Math.sin(theta[i]);
        let v = Math.exp(-rho[i] * rho[i] / (2 * sig * sig));
        for (let k = 0; k < valence; k++) { const tk = ph0 + 2 * Math.PI * k / valence, dx = x - a * Math.cos(tk), dy = y - a * Math.sin(tk); v += Math.exp(-(dx * dx + dy * dy) / (2 * sig * sig)); }
        u[i] = amplitude * v;
      }
      return u;
    }
    for (let i = 0; i < this.n; i++) {
      const env = Math.exp(-Math.pow(rho[i] / Math.max(radius * a, 1e-9), 4));
      if (env < 1e-6) continue;
      const ph = theta[i] * valence / 6 + ph0, zx = rho[i] * Math.cos(ph), zy = rho[i] * Math.sin(ph);
      let s = 0;
      for (let j = 0; j < 3; j++) { const ang = 2 * Math.PI * j / 3 + Math.PI / 2; s += Math.cos(q0 * (Math.cos(ang) * zx + Math.sin(ang) * zy)); }
      u[i] = amplitude * (s / 3) * env;
    }
    return u;
  }
}

// ------------------------------------------------------------------ small vector helpers
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dotv = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const norm = (a) => Math.hypot(a[0], a[1], a[2]);
function cotAt(p, q, r) { // cot of the angle at p between (q−p) and (r−p)
  const u = sub(q, p), v = sub(r, p); const c = cross(u, v); const s = norm(c); return s < 1e-300 ? 0 : dotv(u, v) / s;
}
export function medianEdge(V, F) {
  const nf = F.length / 3, l = new Float64Array(nf * 3);
  for (let f = 0; f < nf; f++) for (let s = 0; s < 3; s++) { const a = F[3 * f + s], b = F[3 * f + (s + 1) % 3]; l[3 * f + s] = Math.hypot(V[3 * a] - V[3 * b], V[3 * a + 1] - V[3 * b + 1], V[3 * a + 2] - V[3 * b + 2]); }
  return l.sort()[Math.floor(l.length / 2)];
}
