// Stripe phases for the three families of a 6-RoSy field (after Knöppel, Crane, Pinkall, Schröder 2015, "Stripe
// patterns on surfaces"), as one rank-3 bundle: at each vertex three complex unknowns ψ^(k) (one per stripe family,
// wave vector ⊥ lattice branch k), coupled along every edge to the family the branch is carried to by the connection
// (the three families permute around RoSy singularities, and a family may arrive reversed — then ψ is conjugated).
// Energy Σ_edges w |ψ_j^(k') − R_ij ψ_i^(k)|² with R = rotation by the phase increment ω_ij = ∫ q0 n·dl (times a
// conjugation when reversed), minimised under Σ m|ψ|² = 1 → smallest eigenvector, inverse iteration with Jacobi-PCG.
// Branch points (faces around which a family's phase winds) are where the stripes must end or start: the dislocations
// the geometry forces. u0 = ⅓ Σ_k Re ψ̂^(k) is a hexagonal initial condition for SH over the whole surface — no seed.
// The solve can run on the base carrier and be lifted to a refined one (liftStripes): ψ is a slow modulation of the
// carrier wave, so transporting it across the fine edges reproduces it.
import { csrFromCOO, pcg } from './sparse.js';
import { wrapPi, intrinsicFrames } from './rosy.js';

const TWO_PI = 2 * Math.PI;

/** per-face lengths → edge lengths, cotan weights and lumped mass */
function meshData(Fi, Len, n) {
  const nf = Fi.length / 3, key = (i, j) => i < j ? i * n + j : j * n + i;
  const lenOf = new Map(), wmap = new Map(), m = new Float64Array(n);
  for (let f = 0; f < nf; f++) {
    const l = [Len[3 * f], Len[3 * f + 1], Len[3 * f + 2]];
    const s = 0.5 * (l[0] + l[1] + l[2]), area = Math.sqrt(Math.max(s * (s - l[0]) * (s - l[1]) * (s - l[2]), 0));
    for (let k = 0; k < 3; k++) {
      const v = Fi[3 * f + k], a = Fi[3 * f + (k + 1) % 3], b = Fi[3 * f + (k + 2) % 3];
      m[v] += area / 3; lenOf.set(key(v, a), l[k]);
      const la = l[k], lb = l[(k + 2) % 3], lo = l[(k + 1) % 3];
      const cs = (la * la + lb * lb - lo * lo) / (2 * la * lb), sn = Math.sqrt(Math.max(1 - cs * cs, 1e-12));
      const kk = key(a, b); wmap.set(kk, (wmap.get(kk) || 0) + 0.5 * cs / sn);
    }
  }
  return { lenOf, wmap, m, key };
}

/** edge coupling of the rank-3 bundle: edge(p, q, k) → {kp, rev, om, mis} for family k at p carried to q */
export function stripeCoupling(frames, theta, N, q0, lenOf, key) {
  const step = TWO_PI / N;
  const angIn = (v, w) => { const f = frames[v], i = f.nb.indexOf(w); return i < 0 ? NaN : f.ang[i]; };
  const nu = (v, k) => theta[v] + k * step + Math.PI / 2;
  return (p, q, k) => {
    const ap = angIn(p, q), aq = angIn(q, p); if (!Number.isFinite(ap) || !Number.isFinite(aq)) return null;
    const r = aq + Math.PI - ap, target = nu(p, k) + r;
    let best = null;
    for (let kp = 0; kp < 3; kp++) for (const s of [1, -1]) {
      const d = Math.abs(wrapPi(nu(q, kp) + (s < 0 ? Math.PI : 0) - target));
      if (!best || d < best.d) best = { kp, s, d };
    }
    const l = lenOf.get(key(p, q)) || 0;
    const om = q0 * l * (Math.cos(nu(p, k) - ap) + best.s * Math.cos(nu(q, best.kp) - aq - Math.PI)) / 2;
    return { kp: best.kp, rev: best.s < 0, om, mis: best.d };
  };
}

export function u0FromPsi(psi, n) {
  const u0 = new Float64Array(n);
  for (let i = 0; i < n; i++) { let s = 0; for (let k = 0; k < 3; k++) { const re = psi[6 * i + 2 * k], im = psi[6 * i + 2 * k + 1], a = Math.hypot(re, im); s += a > 1e-12 ? re / a : 0; } u0[i] = s / 3; }
  return u0;
}

/** faces around which a family's phase winds (families that permute around a RoSy singularity are skipped) */
export function branchPoints(Fi, edge, psi, reversed = false) {
  const nf = Fi.length / 3, arg = (v, k) => Math.atan2(psi[6 * v + 2 * k + 1], psi[6 * v + 2 * k]), out = [];
  for (let f = 0; f < nf; f++) {
    const vs = reversed ? [Fi[3 * f], Fi[3 * f + 2], Fi[3 * f + 1]] : [Fi[3 * f], Fi[3 * f + 1], Fi[3 * f + 2]];   // loop in the ring orientation
    for (let k0 = 0; k0 < 3; k0++) {
      let k = k0, s = 1, Phi = arg(vs[0], k0), sum = 0, ok = true;
      for (let e = 0; e < 3; e++) {
        const p = vs[e], q = vs[(e + 1) % 3]; const ed = edge(p, q, k); if (!ed) { ok = false; break; }
        const s2 = ed.rev ? -s : s, Phi2 = s2 * arg(q, ed.kp);
        sum += wrapPi(Phi2 - Phi - s * ed.om); Phi = Phi2; k = ed.kp; s = s2;
      }
      if (!ok || k !== k0 || s !== 1) continue;
      const wnd = Math.round(sum / TWO_PI); if (wnd !== 0) out.push({ face: f, family: k0, winding: wnd });
    }
  }
  return out;
}

/** stripes(c, rosy, q0, {Len}) → {psi: Float64Array(6n) [re,im per family], u0, branch, iters, energy, worstMatch} */
export function stripes(c, rosy, q0, { Len = null, outer = 12, cgTol = 1e-5, cgMax = 200, log = null } = {}) {   // moderately loose on purpose: the eigenvalue is near-degenerate on flat regions and SH relaxes the result anyway
  const { Fi, n } = c, L = Len || c.Li, fr = rosy.frames, th = rosy.theta, N = rosy.N;
  const { lenOf, wmap, m, key } = meshData(Fi, L, n);
  let wmed = 0; { const ws = Array.from(wmap.values()).sort((a, b) => a - b); wmed = ws[ws.length >> 1] || 1; }
  const wmin = 1e-3 * wmed;
  const edge = stripeCoupling(fr, th, N, q0, lenOf, key);
  const rows = [], cols = [], vals = [];
  const add = (r, c2, v) => { rows.push(r); cols.push(c2); vals.push(v); };
  let worst = 0;
  for (let i = 0; i < n; i++) for (const j of fr[i].nb) {
    if (j < i) continue;
    let w = wmap.get(key(i, j)) || 0; if (!(w > wmin)) w = wmin;
    for (let k = 0; k < 3; k++) {
      const e = edge(i, j, k); if (!e) continue; if (e.mis > worst) worst = e.mis;
      const c1 = Math.cos(e.om), s1 = Math.sin(e.om);
      const R = e.rev ? [c1, -s1, -s1, -c1] : [c1, -s1, s1, c1];
      const I = 6 * i + 2 * k, J = 6 * j + 2 * e.kp;
      add(I, I, w); add(I + 1, I + 1, w); add(J, J, w); add(J + 1, J + 1, w);
      add(J, I, -w * R[0]); add(J, I + 1, -w * R[1]); add(J + 1, I, -w * R[2]); add(J + 1, I + 1, -w * R[3]);
      add(I, J, -w * R[0]); add(I + 1, J, -w * R[1]); add(I, J + 1, -w * R[2]); add(I + 1, J + 1, -w * R[3]);
    }
  }
  const A = csrFromCOO(6 * n, rows, cols, vals);
  const diag = new Float64Array(6 * n); for (let i = 0; i < 6 * n; i++) for (let k = A.indptr[i]; k < A.indptr[i + 1]; k++) if (A.indices[k] === i) diag[i] += A.data[k];
  let mTot = 0, dTot = 0; for (let i = 0; i < n; i++) { mTot += m[i]; dTot += diag[6 * i]; }
  const eps = 1e-6 * dTot / mTot;
  const applyA = (x, out) => { for (let i = 0; i < 6 * n; i++) { let s = 0; for (let k = A.indptr[i]; k < A.indptr[i + 1]; k++) s += A.data[k] * x[A.indices[k]]; out[i] = s + eps * m[i / 6 | 0] * x[i]; } };
  const jac = new Float64Array(6 * n); for (let i = 0; i < 6 * n; i++) jac[i] = diag[i] + eps * m[i / 6 | 0];
  const normalise = (v) => { let s = 0; for (let i = 0; i < 6 * n; i++) s += m[i / 6 | 0] * v[i] * v[i]; s = Math.sqrt(s) || 1; for (let i = 0; i < 6 * n; i++) v[i] /= s; };
  let z = new Float64Array(6 * n); { let s = 777; const rnd = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff - 0.5; }; for (let i = 0; i < 6 * n; i++) z[i] = rnd(); } normalise(z);
  const rhs = new Float64Array(6 * n); let iters = 0, lam = 0;
  for (let it = 0; it < outer; it++) {
    for (let i = 0; i < 6 * n; i++) rhs[i] = m[i / 6 | 0] * z[i];
    const res = pcg(applyA, jac, rhs, z, { maxIter: cgMax, tol: cgTol }); iters += res.iters;
    const zn = res.x; normalise(zn); const Az = new Float64Array(6 * n); applyA(zn, Az); let num = 0, ch = 0;
    for (let i = 0; i < 6 * n; i++) { num += zn[i] * Az[i]; ch += (zn[i] - z[i]) ** 2; }
    lam = num; z = zn; if (log) log(`stripes it ${it}: λ = ${lam.toExponential(3)}, Δ = ${Math.sqrt(ch).toExponential(2)}`);
    if (it > 1 && Math.sqrt(ch) < 3e-4) break;
  }
  return { psi: z, u0: u0FromPsi(z, n), branch: branchPoints(Fi, edge, z, !!rosy.reversed), iters, energy: lam, worstMatch: worst, edge };
}

/**
 * Lift stripe phases solved on the base carrier to a refined carrier whose rotation field (thetaFine, framesFine) has
 * already been lifted: base vertices keep ψ (same family labelling), new vertices are filled layer by layer by
 * transporting ψ across the fine edges (with the fine coupling), then a few sweeps smooth the new vertices.
 * Returns {psi, u0, branch, edge}.
 */
export function liftStripes(psiBase, nBase, cFine, thetaFine, framesFine, q0, N = 6, { sweeps = 2, reversed = false } = {}) {
  const { Fi, Li, n } = cFine;
  const { lenOf, key } = meshData(Fi, Li, n);
  const edge = stripeCoupling(framesFine, thetaFine, N, q0, lenOf, key);
  const psi = new Float64Array(6 * n), done = new Uint8Array(n);
  psi.set(psiBase.subarray(0, 6 * nBase)); for (let i = 0; i < nBase; i++) done[i] = 1;
  const assign = (m, onlyDone) => {
    const acc = new Float64Array(6); let any = false;
    for (const k of framesFine[m].nb) {
      for (let kk = 0; kk < 3; kk++) {
        const e = edge(k, m, kk); if (!e) continue;
        const re = psi[6 * k + 2 * kk], im = psi[6 * k + 2 * kk + 1]; if (re === 0 && im === 0) continue;
        const c1 = Math.cos(e.om), s1 = Math.sin(e.om); let tr = c1 * re - s1 * im, ti = s1 * re + c1 * im; if (e.rev) ti = -ti;
        acc[2 * e.kp] += tr; acc[2 * e.kp + 1] += ti; any = true;
      }
    }
    return any ? acc : null;
  };
  let frontier = []; for (let i = 0; i < nBase; i++) frontier.push(i);
  while (frontier.length) {
    const next = new Set();
    for (const i of frontier) for (const m of framesFine[i].nb) if (!done[m]) next.add(m);
    const layer = Array.from(next), vals = layer.map(m => assign(m, true));
    layer.forEach((m, k) => { if (vals[k]) psi.set(vals[k], 6 * m); done[m] = 1; });
    frontier = layer;
  }
  for (let s = 0; s < sweeps; s++) {
    const np = Float64Array.from(psi);
    for (let m = nBase; m < n; m++) { const v = assign(m, false); if (v) np.set(v, 6 * m); }
    psi.set(np);
  }
  return { psi, u0: u0FromPsi(psi, n), branch: branchPoints(Fi, edge, psi, reversed), edge };
}
