// Globally smoothest N-RoSy field on the carrier (Knöppel, Crane, Pinkall, Schröder 2013), intrinsic version:
// vertex frames from the intrinsic triangulation (ring of neighbours, corner angles from the edge lengths, rescaled to
// 2π at interior vertices), parallel transport across an edge by matching the edge direction in both frames, cotan
// weights, and the complex connection Laplacian L_ij = −w_ij e^{−iN r_ij}. The field is the eigenvector of the smallest
// eigenvalue of L z = λ M z (inverse iteration with Jacobi-PCG on the real 2n embedding). Singularities are read per
// face from the winding of the phase relative to the connection; index m/N, valence 6 − m for N = 6.
// Works on any metric given as (Fi, lengths): the model metric (Li) or the size-field metric ĝ (rescaled lengths).
import { csrFromCOO, pcg } from './sparse.js';

const TWO_PI = 2 * Math.PI;
export const wrapPi = (x) => x - TWO_PI * Math.floor((x + Math.PI) / TWO_PI);

/** Rings of a triangulation (the intrinsic Fi or the extrinsic F) that need not be consistently oriented: along the
 *  seam of a non-orientable carrier (a Möbius band) the faces around a vertex come with both windings. Per vertex the
 *  neighbours in the order of the faces around it, walked through the shared edges whatever the faces' winding, with
 *  cumulative corner angles (rescaled to 2π at interior vertices). Where the faces agree, the walk follows their
 *  winding (against it with `reverse`) and starts where the oriented rings always started — the open side of a
 *  boundary vertex, else the first corner — so an orientable carrier gets exactly the rings it had. Where they disagree,
 *  the walk takes the direction in which most faces run counter-clockwise about orient.N (orient = {V, N}; without it,
 *  the direction of most windings). CA: the corner angles (3 per face); CL (optional): the length of each corner's side
 *  to the next ring vertex, for h. Returns rings [{nb, ang, total, h, sg?}] with rings.eps: Int8Array(3 nf), +1 where
 *  the ring of F[c] passes face ⌊c/3⌋ along its (reversed) winding, −1 against; and per ring sg[i] = σ(v, nb[i]) = ±1,
 *  whether the two rings turn the same way across that edge. eps and sg only where some value is −1 (else null /
 *  absent): a parallel transport across a σ = −1 edge is a reflection. */
export function buildRings(F, n, boundary, CA, CL = null, reverse = false, orient = null) {
  const nf = F.length / 3, nc = 3 * nf, nx3 = (c) => (c % 3 === 2 ? c - 2 : c + 1), pv3 = (c) => (c % 3 === 0 ? c + 2 : c - 1);
  const off = new Int32Array(n + 1); for (let c = 0; c < nc; c++) off[F[c] + 1]++;
  for (let i = 0; i < n; i++) off[i + 1] += off[i];
  const fill = off.slice(0, n), cor = new Int32Array(nc); for (let c = 0; c < nc; c++) cor[fill[F[c]]++] = c;   // corners per vertex, face order
  const eps = new Int8Array(nc).fill(1), rings = new Array(n), A = [], B = [], used = [];
  let anyMixed = false;
  for (let v = 0; v < n; v++) {
    const s0 = off[v], m = off[v + 1] - s0;
    if (!m) { rings[v] = { nb: [], ang: [], total: 0, h: 0 }; continue; }
    A.length = B.length = used.length = 0;
    for (let i = 0; i < m; i++) { const c = cor[s0 + i]; A.push(reverse ? F[pv3(c)] : F[nx3(c)]); B.push(reverse ? F[nx3(c)] : F[pv3(c)]); used.push(0); }
    const count = (w) => { let k = 0; for (let i = 0; i < m; i++) { if (A[i] === w) k++; if (B[i] === w) k++; } return k; };
    let st = -1, d = 1;                                                       // boundary: begin at the open side
    for (let i = 0; i < m; i++) if (count(A[i]) === 1) { st = i; break; }
    if (st < 0) for (let i = 0; i < m; i++) if (count(B[i]) === 1) { st = i; d = -1; break; }
    if (st < 0) st = 0;
    let nb = [d > 0 ? A[st] : B[st]], ang = [0], cur = st, total = 0, guard = 0, mixed = false; used[st] = 1;
    while (guard++ <= m) {
      const c = cor[s0 + cur]; total += CA[c]; eps[c] = d; if (d < 0) mixed = true;
      const ex = d > 0 ? B[cur] : A[cur];
      if (ex === nb[0]) break;                                                // back at the start: closed
      let k2 = -1, d2 = 1;
      for (let i = 0; i < m; i++) { if (used[i]) continue; if (A[i] === ex) { k2 = i; d2 = 1; break; } if (B[i] === ex) { k2 = i; d2 = -1; break; } }
      nb.push(ex); ang.push(total);
      if (k2 < 0) break;
      used[k2] = 1; cur = k2; d = d2;
    }
    if (mixed) {                                                             // faces of both windings: pick the direction
      anyMixed = true; let score = 0;
      for (let i = 0; i < m; i++) {
        const c = cor[s0 + i]; let w = 1;
        if (orient) {
          const { V, N } = orient, f = c - c % 3, a = 3 * F[f], b = 3 * F[f + 1], e = 3 * F[f + 2];
          const ux = V[b] - V[a], uy = V[b + 1] - V[a + 1], uz = V[b + 2] - V[a + 2], wx = V[e] - V[a], wy = V[e + 1] - V[a + 1], wz = V[e + 2] - V[a + 2];
          const dn = (uy * wz - uz * wy) * N[3 * v] + (uz * wx - ux * wz) * N[3 * v + 1] + (ux * wy - uy * wx) * N[3 * v + 2];
          w = (reverse ? -dn : dn) < 0 ? -1 : 1;
        }
        score += eps[c] * w;
      }
      if (score < 0) {
        const closedWalk = nb.length === m && !boundary[v];
        if (closedWalk) { nb = [nb[0], ...nb.slice(1).reverse()]; ang = [0, ...ang.slice(1).reverse().map(x => total - x)]; }
        else { nb = nb.reverse(); ang = ang.reverse().map(x => total - x); }
        for (let i = 0; i < m; i++) eps[cor[s0 + i]] *= -1;
      }
    }
    const closed = !boundary[v] && nb.length === m, scale = closed && total > 1e-9 ? TWO_PI / total : 1;
    let hs = 0; if (CL) for (let i = 0; i < m; i++) hs += CL[cor[s0 + i]];
    rings[v] = { nb, ang: ang.map(x => x * scale), total: total * scale, h: hs / m };
  }
  rings.eps = null;
  if (!anyMixed) { let neg = false; for (let c = 0; c < nc && !neg; c++) if (eps[c] < 0) neg = true; if (!neg) return rings; }
  // σ across every ring edge, from a face holding both ends: σ = eps(v's corner) · eps(w's corner)
  rings.eps = eps;
  for (let v = 0; v < n; v++) {
    const r = rings[v], s0 = off[v], m = off[v + 1] - s0; let sg = null;
    for (let i = 0; i < r.nb.length; i++) {
      const w = r.nb[i]; let s = 1;
      for (let j = 0; j < m; j++) {
        const c = cor[s0 + j], f = c - c % 3;
        const k = F[f] === w ? f : F[f + 1] === w ? f + 1 : F[f + 2] === w ? f + 2 : -1;
        if (k >= 0) { s = eps[c] * eps[k]; break; }
      }
      if (s < 0) { if (!sg) sg = new Int8Array(r.nb.length).fill(1); sg[i] = -1; }
    }
    if (sg) r.sg = sg;
  }
  return rings;
}

/** Rings and frames on an intrinsic triangulation. Returns {nb, ang, total, h} per vertex: ordered neighbours,
 *  cumulative angle of each (rescaled so the total is 2π at interior vertices), and the mean incident edge length;
 *  on a non-orientable carrier also sg (per neighbour) and frames.eps (buildRings). orient = {V, N} (optional). */
export function intrinsicFrames(Fi, Len, n, boundary, reverse = false, orient = null) {
  const nf = Fi.length / 3, CA = new Float64Array(3 * nf), CL = new Float64Array(3 * nf);
  for (let f = 0; f < nf; f++) {
    const l = [Len[3 * f], Len[3 * f + 1], Len[3 * f + 2]];
    for (let k = 0; k < 3; k++) {
      // ring order follows the face orientation; `reverse` when Fi is wound opposite to F (Python-exported carriers)
      const la = reverse ? l[(k + 2) % 3] : l[k], lb = reverse ? l[k] : l[(k + 2) % 3], lo = l[(k + 1) % 3];   // sides v→a, b→v, opposite
      const cs = Math.max(-1, Math.min(1, (la * la + lb * lb - lo * lo) / (2 * la * lb)));
      CA[3 * f + k] = Math.acos(cs); CL[3 * f + k] = la;
    }
  }
  return buildRings(Fi, n, boundary, CA, CL, reverse, orient);
}

/** parallel transport of a tangent angle along the edge p → q of the frames: T(α) = s·α + r — s = −1 (a reflection)
 *  across the seam of a non-orientable carrier, where the two rings turn opposite ways. null if q is not in p's ring. */
export function transport(frames, p, q) {
  const fp = frames[p], fq = frames[q], ip = fp.nb.indexOf(q), iq = fq.nb.indexOf(p);
  if (ip < 0 || iq < 0) return null;
  const s = fp.sg ? fp.sg[ip] : 1;
  return { s, r: s > 0 ? fq.ang[iq] + Math.PI - fp.ang[ip] : fq.ang[iq] + Math.PI + fp.ang[ip] };
}

/** true when the intrinsic faces are wound clockwise about the carrier normals (sampled over faces), so that the rings
 *  built from them must be reversed to run counter-clockwise about N — the convention the index signs and the 3-D
 *  directions assume. (Python-exported carriers have Fi wound opposite to F.) */
export function fiReversed(c) {
  const { Fi, V, N } = c, nf = Fi.length / 3, step = Math.max(1, Math.floor(nf / 500)); let s = 0;
  for (let f = 0; f < nf; f += step) {
    const a = Fi[3 * f], b = Fi[3 * f + 1], cc = Fi[3 * f + 2];
    const ux = V[3 * b] - V[3 * a], uy = V[3 * b + 1] - V[3 * a + 1], uz = V[3 * b + 2] - V[3 * a + 2];
    const wx = V[3 * cc] - V[3 * a], wy = V[3 * cc + 1] - V[3 * a + 1], wz = V[3 * cc + 2] - V[3 * a + 2];
    const nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx;
    s += Math.sign(nx * (N[3 * a] + N[3 * b] + N[3 * cc]) + ny * (N[3 * a + 1] + N[3 * b + 1] + N[3 * cc + 1]) + nz * (N[3 * a + 2] + N[3 * b + 2] + N[3 * cc + 2]));
  }
  return s < 0;
}

/** cotan weights (clamped ≥ wmin · median) and lumped mass from lengths; returns {W: Map('i,j'→w) as CSR-ready arrays, m} */
function cotanFromLengths(Fi, Len, n) {
  const nf = Fi.length / 3, m = new Float64Array(n), wmap = new Map();
  const key = (i, j) => i < j ? i * n + j : j * n + i;
  for (let f = 0; f < nf; f++) {
    const l = [Len[3 * f], Len[3 * f + 1], Len[3 * f + 2]];
    const s = 0.5 * (l[0] + l[1] + l[2]), area = Math.sqrt(Math.max(s * (s - l[0]) * (s - l[1]) * (s - l[2]), 0));
    for (let k = 0; k < 3; k++) {
      const v = Fi[3 * f + k], a = Fi[3 * f + (k + 1) % 3], b = Fi[3 * f + (k + 2) % 3];
      m[v] += area / 3;
      // cot of the angle at v, opposite edge (a, b)
      const la = l[k], lb = l[(k + 2) % 3], lo = l[(k + 1) % 3];
      const cs = (la * la + lb * lb - lo * lo) / (2 * la * lb), sn = Math.sqrt(Math.max(1 - cs * cs, 1e-12));
      const kk = key(a, b); wmap.set(kk, (wmap.get(kk) || 0) + 0.5 * cs / sn);
    }
  }
  return { wmap, m };
}

/**
 * Smoothest N-RoSy. opts: {N = 6, Len = carrier.Li (or ĝ lengths), outer = 25, cgTol = 1e-6}.
 * Returns {N, frames, re, im, theta (angle of one branch in the vertex frame), mag (|z|/max), idx (Int8 per face),
 *          singular: [{face, index, valence}], sumIndex, energy, iters}
 */
export function smoothestRoSy(c, { N = 6, Len = null, outer = 20, cgTol = 1e-5, cgMax = 300, log = null, constraints = null, x0 = null, mbo = null } = {}) {
  const { Fi, n, boundary } = c, L = Len || c.Li, reversed = fiReversed(c);
  const frames = intrinsicFrames(Fi, L, n, boundary, reversed, c.V && c.N ? { V: c.V, N: c.N } : null);
  const { wmap, m } = cotanFromLengths(Fi, L, n);
  // directed edge data from the rings: transport r_ij for every (i → j) in i's ring
  const rows = [], cols = [], vre = [], vim = [], diag = new Float64Array(n), refl = [];
  let wmed = 0; { const ws = Array.from(wmap.values()).sort((a, b) => a - b); wmed = ws[ws.length >> 1] || 1; }
  const wmin = 1e-3 * wmed;
  const angIn = (v, w) => { const f = frames[v], i = f.nb.indexOf(w); return i < 0 ? NaN : f.ang[i]; };
  const rOf = (i, j) => angIn(j, i) + Math.PI - angIn(i, j);                        // vector angle t_i → t_i + r_ij
  for (let i = 0; i < n; i++) {
    const f = frames[i];
    for (let q = 0; q < f.nb.length; q++) {
      const j = f.nb[q]; const r = rOf(i, j); if (!Number.isFinite(r)) continue;
      let w = wmap.get(i < j ? i * n + j : j * n + i) || 0; if (!(w > wmin)) w = wmin;
      if (f.sg && f.sg[q] < 0) {
        // the seam of a non-orientable carrier: the transport is a reflection, t_j = ρ − t_i, so the field arrives
        // conjugated — the energy w|z_j − e^{iNρ} z̄_i|² couples (re, im) by the symmetric block −w [[cos Nρ, sin Nρ], [sin Nρ, −cos Nρ]]
        const rho = angIn(j, i) + Math.PI + f.ang[q];
        refl.push(i, j, w * Math.cos(N * rho), w * Math.sin(N * rho)); diag[i] += w; continue;
      }
      // L_ij = −w e^{−iN r_ij}
      rows.push(i); cols.push(j); vre.push(-w * Math.cos(N * r)); vim.push(w * Math.sin(N * r)); diag[i] += w;
    }
  }
  for (let i = 0; i < n; i++) { rows.push(i); cols.push(i); vre.push(diag[i]); vim.push(0); }
  const Ar = csrFromCOO(n, rows, cols, vre), Ai = csrFromCOO(n, rows, cols, vim), nR = refl.length / 4;
  let mTot = 0, dTot = 0; for (let i = 0; i < n; i++) { mTot += m[i]; dTot += diag[i]; }
  const eps = 1e-6 * dTot / mTot;                                                  // tiny shift keeps A positive definite
  // real embedding: [re; im] → [Ar re − Ai im; Ai re + Ar im] (+ the reflection blocks) + eps M. Ar and Ai come from the
  // same COO lists, so they share one pattern: one pass reads it once (the four row sums kept apart, as four passes would)
  const ip = Ar.indptr, ix = Ar.indices, dr = Ar.data, di = Ai.data, em = new Float64Array(n);
  for (let i = 0; i < n; i++) em[i] = eps * m[i];
  const applyA = (x, out) => {
    for (let i = 0; i < n; i++) {
      let s1 = 0, s2 = 0, s3 = 0, s4 = 0;
      for (let k = ip[i], e = ip[i + 1]; k < e; k++) { const j = ix[k], xr = x[j], xi = x[n + j], ar = dr[k], ai = di[k]; s1 += ar * xr; s2 += ai * xi; s3 += ai * xr; s4 += ar * xi; }
      if (nR) { out[i] = s1 - s2; out[n + i] = s3 + s4; }
      else { out[i] = s1 - s2 + em[i] * x[i]; out[n + i] = s3 + s4 + em[i] * x[n + i]; }
    }
    if (!nR) return;
    for (let k = 0; k < nR; k++) {
      const i = refl[4 * k], j = refl[4 * k + 1], wc = refl[4 * k + 2], ws = refl[4 * k + 3];
      out[i] -= wc * x[j] + ws * x[n + j]; out[n + i] -= ws * x[j] - wc * x[n + j];
    }
    for (let i = 0; i < n; i++) { out[i] += em[i] * x[i]; out[n + i] += em[i] * x[n + i]; }
  };
  applyA.par = { kind: 'rosy', n, ip, ix, dr, di, em, refl: nR ? refl : null };   // the same operator for the solver threads (parallel.js)
  const jac = new Float64Array(2 * n); for (let i = 0; i < n; i++) jac[i] = jac[n + i] = diag[i] + eps * m[i];
  const normalise = (v) => { let s = 0; for (let i = 0; i < n; i++) s += m[i] * (v[i] * v[i] + v[n + i] * v[n + i]); s = Math.sqrt(s); for (let i = 0; i < 2 * n; i++) v[i] /= s; };
  let z = new Float64Array(2 * n), iters = 0, lam = 0;
  const cons = (constraints || []).filter(q => q.weight > 0 && q.vertex >= 0 && q.vertex < n && Number.isFinite(q.theta));
  if (cons.length) {
    // alignment constraints (Knöppel et al. 2013): minimise z*Az + Σ w_v |z_v − z̄_v|², z̄_v = e^{iNθ_v} → (A + W) z = W z̄, one
    // linear solve instead of the eigenproblem; the weight grows from a nudge (0.3: 4× the vertex's own coupling) to a
    // firm pin (1: 600×); a small weight mostly picks the global rotation of the free field
    const W = new Float64Array(n), rhs = new Float64Array(2 * n);
    for (const q of cons) { const wq = Math.min(1, q.weight), w = 30 * wq * wq / (1.05 - wq) * diag[q.vertex]; W[q.vertex] += w; rhs[q.vertex] += w * Math.cos(N * q.theta); rhs[n + q.vertex] += w * Math.sin(N * q.theta); }
    const applyC = (x, out) => { applyA(x, out); for (let i = 0; i < n; i++) { out[i] += W[i] * x[i]; out[n + i] += W[i] * x[n + i]; } };
    const jacC = Float64Array.from(jac); for (let i = 0; i < n; i++) { jacC[i] += W[i]; jacC[n + i] += W[i]; }
    const res = pcg(applyC, jacC, rhs, x0, { maxIter: 6 * cgMax, tol: cgTol }); z = res.x; iters = res.iters;
    if (log) log(`rosy with ${cons.length} alignment constraints: ${iters} PCG iterations, relres ${res.relres.toExponential(2)}`);
  } else {
  // inverse iteration
  { let s = 12345; const rnd = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff - 0.5; }; for (let i = 0; i < 2 * n; i++) z[i] = rnd(); }
  normalise(z);
  const rhs = new Float64Array(2 * n), Az = new Float64Array(2 * n);
  for (let it = 0; it < outer; it++) {
    for (let i = 0; i < n; i++) { rhs[i] = m[i] * z[i]; rhs[n + i] = m[i] * z[n + i]; }
    const res = pcg(applyA, jac, rhs, z, { maxIter: cgMax, tol: cgTol }); iters += res.iters;
    const znew = res.x; normalise(znew);
    // Rayleigh quotient and change
    applyA(znew, Az); let num = 0, ch = 0;
    for (let i = 0; i < 2 * n; i++) { num += znew[i] * Az[i]; ch += (znew[i] - z[i]) ** 2; }
    lam = num; z = znew; if (log) log(`rosy it ${it}: λ = ${lam.toExponential(3)}, Δ = ${Math.sqrt(ch).toExponential(2)}`);
    if (it > 2 && Math.sqrt(ch) < 1e-5) break;
  }
  }
  // Unit-length relaxation (diffusion-generated motion, after Merriman–Bence–Osher). The eigenvector is the relaxed
  // problem: |z| may shrink where the field turns, and that parks the singular cores away from the curvature that should
  // hold them (cylinder: 1.3–1.8 a inside the flat caps, 52 forced branch points at the creases). Alternating "normalise to
  // |z| = 1" with a short implicit diffusion (M + τ A) z' = M z moves the cores to the unit-length optimum, where the
  // curvature screens them (cylinder: all 12 on the creases, 0–4 branch points). mbo = {tau: core radius² (model units), iters}
  let mboInfo = null;
  if (mbo && mbo.iters > 0 && !cons.length) {
    const tau = mbo.tau, rhsB = new Float64Array(2 * n), jacB = new Float64Array(2 * n), tmpB = new Float64Array(2 * n);
    for (let i = 0; i < n; i++) jacB[i] = jacB[n + i] = m[i] + tau * (diag[i] + eps * m[i]);
    const applyB = nR ? (x, out) => { applyA(x, tmpB); for (let i = 0; i < n; i++) { out[i] = m[i] * x[i] + tau * tmpB[i]; out[n + i] = m[i] * x[n + i] + tau * tmpB[n + i]; } }
      : (x, out) => {                                                         // the same sums as applyA, then M x + τ (A x)
        for (let i = 0; i < n; i++) {
          let s1 = 0, s2 = 0, s3 = 0, s4 = 0;
          for (let k = ip[i], e = ip[i + 1]; k < e; k++) { const j = ix[k], xr = x[j], xi = x[n + j], ar = dr[k], ai = di[k]; s1 += ar * xr; s2 += ai * xi; s3 += ai * xr; s4 += ar * xi; }
          out[i] = m[i] * x[i] + tau * (s1 - s2 + em[i] * x[i]); out[n + i] = m[i] * x[n + i] + tau * (s3 + s4 + em[i] * x[n + i]);
        }
      };
    applyB.par = { kind: 'mbo', n, ip, ix, dr, di, em, refl: nR ? refl : null, m, tau };
    for (let it = 0; it < mbo.iters; it++) {
      for (let i = 0; i < n; i++) { const r = Math.hypot(z[i], z[n + i]) || 1; rhsB[i] = m[i] * z[i] / r; rhsB[n + i] = m[i] * z[n + i] / r; }
      const res = pcg(applyB, jacB, rhsB, z, { maxIter: cgMax, tol: cgTol }); iters += res.iters; z = res.x;
    }
    normalise(z); mboInfo = { tau, iters: mbo.iters };
    if (log) log(`rosy: ${mbo.iters} unit-length (MBO) steps, core ${Math.sqrt(tau).toPrecision(3)}`);
  }
  const re = z.subarray(0, n), im = z.subarray(n, 2 * n);
  const theta = new Float64Array(n), mag = new Float64Array(n); let mx = 0;
  for (let i = 0; i < n; i++) { theta[i] = Math.atan2(im[i], re[i]) / N; mag[i] = Math.hypot(re[i], im[i]); if (mag[i] > mx) mx = mag[i]; }
  for (let i = 0; i < n; i++) mag[i] /= mx || 1;
  // per-face index
  const nf = Fi.length / 3, idx = new Int8Array(nf), singular = [], fe = frames.eps;
  for (let f = 0; f < nf; f++) {
    const a = Fi[3 * f], b = Fi[3 * f + 1], cc = Fi[3 * f + 2]; let s = 0, h = 0;
    if (fe && (fe[3 * f] < 0 || fe[3 * f + 1] < 0 || fe[3 * f + 2] < 0)) {
      // a ring passes this face against its winding (non-orientable seam): measure every angle of the three frames
      // in the face's own turning sense (θ → e θ, ring angles → e · angle, e = ±1 per corner), then as below
      const e = (v) => fe[Fi[3 * f] === v ? 3 * f : Fi[3 * f + 1] === v ? 3 * f + 1 : 3 * f + 2];
      for (const [p, q] of (reversed ? [[a, cc], [cc, b], [b, a]] : [[a, b], [b, cc], [cc, a]])) {
        const ep = e(p), eq = e(q), r = eq * angIn(q, p) + Math.PI - ep * angIn(p, q); if (!Number.isFinite(r)) { s = NaN; break; }
        h += r; s += wrapPi(N * eq * theta[q] - N * ep * theta[p] - N * r);
      }
    } else
    for (const [p, q] of (reversed ? [[a, cc], [cc, b], [b, a]] : [[a, b], [b, cc], [cc, a]])) {   // loop in the ring orientation
      const r = rOf(p, q); if (!Number.isFinite(r)) { s = NaN; break; }
      h += r; s += wrapPi(N * theta[q] - N * theta[p] - N * r);
    }
    if (!Number.isFinite(s)) continue;
    const mI = Math.round((s + N * wrapPi(h)) / TWO_PI); idx[f] = mI;
    if (mI !== 0) singular.push({ face: f, index: mI, valence: 6 - mI });
  }
  let sumIndex = 0; for (const sg of singular) sumIndex += sg.index;
  return { N, frames, re, im, theta, mag, idx, singular, sumIndex, energy: lam, iters, m, reversed, z, constrained: cons.length, mbo: mboInfo };
}

/** 3-D unit directions of branch k at every vertex for display, from the angle field and each vertex's first ring
 *  neighbour (the frame's zero direction, projected to the tangent plane of the vertex normal; rings run CCW about N) */
export function directionsFromNb0(c, theta, nb0, N, k = 0) {
  const { V, N: Nrm, n } = c, out = new Float64Array(3 * n), step = TWO_PI / N;
  for (let i = 0; i < n; i++) {
    const j = nb0[i]; if (j < 0) continue;
    let ex = V[3 * j] - V[3 * i], ey = V[3 * j + 1] - V[3 * i + 1], ez = V[3 * j + 2] - V[3 * i + 2];
    const nx = Nrm[3 * i], ny = Nrm[3 * i + 1], nz = Nrm[3 * i + 2]; const d = ex * nx + ey * ny + ez * nz; ex -= d * nx; ey -= d * ny; ez -= d * nz;
    const l = Math.hypot(ex, ey, ez) || 1; ex /= l; ey /= l; ez /= l;
    const px = ny * ez - nz * ey, py = nz * ex - nx * ez, pz = nx * ey - ny * ex;            // n × e (ring direction)
    const t = theta[i] + k * step, cs = Math.cos(t), sn = Math.sin(t);
    out[3 * i] = cs * ex + sn * px; out[3 * i + 1] = cs * ey + sn * py; out[3 * i + 2] = cs * ez + sn * pz;
  }
  return out;
}
export function nb0Of(frames) { const out = new Int32Array(frames.length); for (let i = 0; i < frames.length; i++) out[i] = frames[i].nb.length ? frames[i].nb[0] : -1; return out; }
export function rosyDirections3D(c, rosy, k = 0) { return directionsFromNb0(c, rosy.theta, nb0Of(rosy.frames), rosy.N, k); }

/** nearest carrier vertex of a singular face (for cone lists) */
export function singularVertices(c, rosy) {
  const { Fi, V } = c;
  return rosy.singular.map(s => {
    const f = s.face, ids = [Fi[3 * f], Fi[3 * f + 1], Fi[3 * f + 2]];
    const cx = (V[3 * ids[0]] + V[3 * ids[1]] + V[3 * ids[2]]) / 3, cy = (V[3 * ids[0] + 1] + V[3 * ids[1] + 1] + V[3 * ids[2] + 1]) / 3, cz = (V[3 * ids[0] + 2] + V[3 * ids[1] + 2] + V[3 * ids[2] + 2]) / 3;
    return { vertex: ids[0], valence: s.valence, index: s.index, face: f, pos: [cx, cy, cz] };
  });
}

// ------------------------------------------------------------------ lifting a field from the base carrier to a refined one
/** transport angle for the directed edge p → q from frames (NaN if q is not in p's ring) */
export function transportAngle(frames, p, q) {
  const fp = frames[p], fq = frames[q], ip = fp.nb.indexOf(q), iq = fq.nb.indexOf(p);
  return ip < 0 || iq < 0 ? NaN : fq.ang[iq] + Math.PI - fp.ang[ip];
}

/**
 * Lift an N-RoSy solved on the base carrier (vertices 0 … nBase−1 keep their ids after midpoint subdivision) to the
 * refined carrier: base vertices get their angle re-expressed in the refined frame (the refined ring's zero direction
 * differs from the base ring's), new vertices are filled layer by layer by transporting from already-assigned
 * neighbours, then a few Jacobi sweeps smooth the new vertices. Returns {theta, frames}.
 */
export function liftField(thetaBase, framesBase, nBase, cFine, N = 6, { sweeps = 3, Len = null } = {}) {
  const { Fi, Li, n, boundary, V } = cFine;                                   // Len: the fine carrier's lengths in ĝ, when the field lives there
  const reversed = fiReversed(cFine), frames = intrinsicFrames(Fi, Len || Li, n, boundary, reversed, cFine.V && cFine.N ? { V, N: cFine.N } : null);
  const theta = new Float64Array(n).fill(NaN), done = new Uint8Array(n);
  for (let i = 0; i < nBase; i++) {
    const j0 = framesBase[i].nb[0]; if (j0 === undefined || !frames[i].nb.length) continue;
    const dx = V[3 * j0] - V[3 * i], dy = V[3 * j0 + 1] - V[3 * i + 1], dz = V[3 * j0 + 2] - V[3 * i + 2];
    let best = -2, alpha = 0;
    for (let q = 0; q < frames[i].nb.length; q++) {
      const m = frames[i].nb[q], ex = V[3 * m] - V[3 * i], ey = V[3 * m + 1] - V[3 * i + 1], ez = V[3 * m + 2] - V[3 * i + 2];
      const cs = (dx * ex + dy * ey + dz * ez) / (Math.hypot(dx, dy, dz) * Math.hypot(ex, ey, ez) || 1);
      if (cs > best) { best = cs; alpha = frames[i].ang[q]; }
    }
    theta[i] = thetaBase[i] + alpha; done[i] = 1;
  }
  // layers outward from the assigned set
  let frontier = []; for (let i = 0; i < n; i++) if (done[i]) frontier.push(i);
  const assign = (m, onlyDone) => {
    let re = 0, im = 0; const f = frames[m];
    for (let q = 0; q < f.nb.length; q++) { const k = f.nb[q]; if (onlyDone && !done[k]) continue; if (Number.isNaN(theta[k])) continue;
      const t = transport(frames, k, m); if (!t || !Number.isFinite(t.r)) continue; const a = N * (t.s * theta[k] + t.r); re += Math.cos(a); im += Math.sin(a); }
    return re * re + im * im > 1e-30 ? Math.atan2(im, re) / N : NaN;
  };
  while (frontier.length) {
    const next = new Set();
    for (const i of frontier) for (const m of frames[i].nb) if (!done[m]) next.add(m);
    const layer = Array.from(next), vals = layer.map(m => assign(m, true));
    layer.forEach((m, k) => { theta[m] = vals[k]; done[m] = 1; });
    frontier = layer.filter(m => Number.isFinite(theta[m]));
    if (!frontier.length) break;
  }
  for (let s = 0; s < sweeps; s++) {
    const nt = Float64Array.from(theta);
    for (let m = nBase; m < n; m++) { const v = assign(m, false); if (Number.isFinite(v)) nt[m] = v; }
    theta.set(nt);
  }
  for (let i = 0; i < n; i++) if (Number.isNaN(theta[i])) theta[i] = 0;
  return { theta, frames, reversed };
}
