// Stripe phases for the three families of a 6-RoSy field (after Knöppel, Crane, Pinkall, Schröder 2015, "Stripe
// patterns on surfaces"), as one rank-3 bundle: at each vertex three complex unknowns ψ^(k) (one per stripe family,
// wave vector ⊥ lattice branch k), coupled along every edge to the family the branch is carried to by the connection
// (the three families permute around RoSy singularities, and a family may arrive reversed — then ψ is conjugated).
// Energy Σ_edges w |ψ_j^(k') − R_ij ψ_i^(k)|² with R = rotation by the phase increment ω_ij = ∫ q0 n·dl (times a
// conjugation when reversed), minimised under Σ m|ψ|² = 1 → smallest eigenvector, inverse iteration with Jacobi-PCG.
// Branch points (faces around which a family's phase winds) are where the stripes must end or start: the dislocations
// the geometry forces. u0 = ⅓ Σ_k Re ψ̂^(k) is a hexagonal initial condition for SH over the whole surface — no seed —
// with the triad sign fixed in u0FromPsi (third moment Σ m u0³ < 0 → negate u0 and ψ).
// The solve can run on the base carrier and be lifted to a refined one (liftStripes): ψ is a slow modulation of the
// carrier wave, so transporting it across the fine edges reproduces it. stripeFaces prepares the phases for display.
import { csrFromCOO, pcg, MinHeap } from './sparse.js';
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
  const sgIn = (v, w) => { const f = frames[v]; if (!f.sg) return 1; const i = f.nb.indexOf(w); return i < 0 ? 1 : f.sg[i]; };
  const nu = (v, k) => theta[v] + k * step + Math.PI / 2;
  return (p, q, k) => {
    const ap = angIn(p, q), aq = angIn(q, p); if (!Number.isFinite(ap) || !Number.isFinite(aq)) return null;
    // across the seam of a non-orientable carrier the transport is a reflection: ν ↦ aq + π − (ν − ap)
    const r = aq + Math.PI - ap, target = sgIn(p, q) > 0 ? nu(p, k) + r : aq + Math.PI - (nu(p, k) - ap);
    let bk = -1, bs = 1, bd = 0;
    for (let kp = 0; kp < 3; kp++) for (let s = 1; s >= -1; s -= 2) {
      const d = Math.abs(wrapPi(nu(q, kp) + (s < 0 ? Math.PI : 0) - target));
      if (bk < 0 || d < bd) { bk = kp; bs = s; bd = d; }
    }
    const l = lenOf.get(key(p, q)) || 0;
    const om = q0 * l * (Math.cos(nu(p, k) - ap) + bs * Math.cos(nu(q, bk) - aq - Math.PI)) / 2;
    return { kp: bk, rev: bs < 0, om, mis: bd };
  };
}

/** Seed for SH: u0 = ⅓ Σ_k cos φ^(k) (the three stripe waves added). The stripe solve does not impose the triad
 *  Φ = φ⁽⁰⁾ − φ⁽¹⁾ + φ⁽²⁾ ≡ 0; where it lands on Φ = π the sum is a honeycomb that SH23 (b > 0) erases. The third moment
 *  Σ m u0³ is > 0 for the lattice and < 0 for the honeycomb; negating u0 and ψ adds π to every phase (Φ → Φ + π).
 *  ψ is negated in place so that ψ and u0 stay consistent. m: lumped vertex areas (uniform weights if omitted). */
export function u0FromPsi(psi, n, m = null) {                     // psi: 6 reals per vertex = (Re, Im) of ψ⁽⁰⁾, ψ⁽¹⁾, ψ⁽²⁾; m: vertex areas
  const u0 = new Float64Array(n); let m3 = 0;
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let k = 0; k < 3; k++) {
      const re = psi[6 * i + 2 * k], im = psi[6 * i + 2 * k + 1]; // ψᵢ⁽ᵏ⁾
      const a = Math.hypot(re, im);                               // |ψᵢ⁽ᵏ⁾|
      s += a > 1e-12 ? re / a : 0;                                // Re ψ̂ = cos φᵢ⁽ᵏ⁾
    }
    u0[i] = s / 3;                                                // average of the three
    m3 += (m ? m[i] : 1) * u0[i] ** 3;                            // third moment: > 0 lattice, < 0 honeycomb
  }
  if (m3 < 0) {                                                   // Φ = π: negate, which adds π to every phase
    for (let i = 0; i < n; i++) u0[i] = -u0[i];
    for (let i = 0; i < 6 * n; i++) psi[i] = -psi[i];             // keep ψ consistent with u0
  }
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

/** vertex farthest (in edges) from all `sources` over the rings nb: the root of the colour walk in stripeFaces, so that
 *  the seams start at the field singularities and run away from it */
export function farthestVertex(nb, n, sources) {
  if (!sources.length) return 0;
  const dist = new Int32Array(n).fill(-1), queue = new Int32Array(n); let head = 0, tail = 0, last = sources[0];
  for (const v of sources) if (v >= 0 && v < n && dist[v] < 0) { dist[v] = 0; queue[tail++] = v; }
  while (head < tail) { const p = queue[head++]; last = p; for (const q of nb[p]) if (dist[q] < 0) { dist[q] = dist[p] + 1; queue[tail++] = q; } }
  return last;
}

/** Stripe phases per render face for the view (view.js · setStripes). The three families get colour slots by a walk
 *  outward from `root` in order of distance (Dijkstra over the intrinsic rings, 3-D edge lengths from V), carrying each
 *  slot along the edge coupling (family k at p arrives as family k' at q, possibly reversed = conjugated). Around a field
 *  5 or 7 no consistent labelling exists: the slots swap across a seam where the walk closes round the singularity
 *  (the side away from the root). Each render face takes its slots from a corner that agrees with at least one other
 *  corner (so a face on the seam goes with the side holding two of its corners), transports them to the other corners,
 *  and unwraps the phases inside the face (expected increment ω plus the wrapped remainder) for linear interpolation.
 *  F: render faces (carrier.F); nb: intrinsic rings; edge: stripeCoupling; psi: 6 reals per vertex; V: positions.
 *  Returns Float32Array(9 nf): face f, corner j, slot s → 9f + 3j + s. */
export function stripeFaces(F, n, nb, edge, psi, root = 0, V = null) {
  const lab = new Int8Array(3 * n).fill(-1), sgn = new Int8Array(3 * n), dist = new Float64Array(n).fill(Infinity), heap = new MinHeap(1024);
  const arg = (v, k) => Math.atan2(psi[6 * v + 2 * k + 1], psi[6 * v + 2 * k]);
  const len = (p, q) => V ? Math.hypot(V[3 * q] - V[3 * p], V[3 * q + 1] - V[3 * p + 1], V[3 * q + 2] - V[3 * p + 2]) : 1;
  const parent = new Int32Array(n).fill(-1), done = new Uint8Array(n);
  const label = (p, q) => {                                            // slots of q carried from p; null if not a permutation
    const kk = [0, 0, 0], ss = [1, 1, 1];
    for (let s = 0; s < 3; s++) { const e = edge(p, q, lab[3 * p + s]); if (!e) return null; kk[s] = e.kp; ss[s] = e.rev ? -sgn[3 * p + s] : sgn[3 * p + s]; }
    return kk[0] === kk[1] || kk[1] === kk[2] || kk[0] === kk[2] ? null : [kk, ss];
  };
  let next = 0, r0 = root >= 0 && root < n ? root : 0;
  for (;;) {
    lab[3 * r0] = 0; lab[3 * r0 + 1] = 1; lab[3 * r0 + 2] = 2; sgn[3 * r0] = sgn[3 * r0 + 1] = sgn[3 * r0 + 2] = 1; dist[r0] = 0; heap.push(0, r0);
    while (heap.size) {
      const [d, p] = heap.pop(); if (done[p] || d > dist[p]) continue;
      if (parent[p] >= 0) { const L = label(parent[p], p); if (!L) { done[p] = 0; dist[p] = Infinity; continue; } for (let s = 0; s < 3; s++) { lab[3 * p + s] = L[0][s]; sgn[3 * p + s] = L[1][s]; } }
      done[p] = 1;
      for (const q of nb[p]) { if (done[q]) continue; const nd = d + len(p, q); if (nd < dist[q]) { dist[q] = nd; parent[q] = p; heap.push(nd, q); } }
    }
    while (next < n && done[next]) next++;
    if (next >= n) break;
    r0 = next; parent[r0] = -1;                                         // another component, or a vertex no neighbour could label
  }
  // does the walk's labelling at v agree with slot transport from u (all three slots, signs included)?
  const agrees = (u, v) => { for (let s = 0; s < 3; s++) { const e = edge(u, v, lab[3 * u + s]); if (!e || e.kp !== lab[3 * v + s] || (e.rev ? -sgn[3 * u + s] : sgn[3 * u + s]) !== sgn[3 * v + s]) return false; } return true; };
  const nf = F.length / 3, out = new Float32Array(9 * nf), cs = [0, 0, 0];
  for (let f = 0; f < nf; f++) {
    cs[0] = F[3 * f]; cs[1] = F[3 * f + 1]; cs[2] = F[3 * f + 2];
    let j0 = 0; if (!agrees(cs[0], cs[1]) && !agrees(cs[0], cs[2]) && agrees(cs[1], cs[2])) j0 = 1;   // corner 0 alone on the other side of the seam
    const a = cs[j0];
    for (let s = 0; s < 3; s++) {
      const k = lab[3 * a + s], sg = sgn[3 * a + s], pa = sg * arg(a, k);
      for (let j = 0; j < 3; j++) {
        if (j === j0) { out[9 * f + 3 * j + s] = pa; continue; }
        const v = cs[j], e = edge(a, v, k);
        let kv, sv, om = 0;
        if (e) { kv = e.kp; sv = e.rev ? -sg : sg; om = sg * e.om; }
        else { kv = lab[3 * v + s]; sv = sgn[3 * v + s]; }            // a render edge that is not intrinsic (flipped): the walk's slots at v
        out[9 * f + 3 * j + s] = pa + om + wrapPi(sv * arg(v, kv) - pa - om);
      }
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
  const ap = A.indptr, ax = A.indices, ad = A.data, m6 = new Float64Array(6 * n), em6 = new Float64Array(6 * n);
  for (let i = 0; i < 6 * n; i++) { m6[i] = m[i / 6 | 0]; em6[i] = eps * m6[i]; }   // the vertex mass per unknown, and eps times it
  // the operator in 2×2 blocks: the rows (Re, Im) of a family at a vertex, 6i + 2k and 6i + 2k + 1, couple to a
  // neighbour's family as one 2×2 block (the diagonal is w·I), so one column index per block instead of four. Within a
  // row the products are added in the same column order as the CSR row (an absent entry adds 0), so the same numbers;
  // applyA also returns dot(x, out) for the CG
  const nq = 3 * n, bp = new Int32Array(nq + 1); let nbk = 0;
  for (let q = 0; q < nq; q++) {                                             // count the column pairs of each row pair
    let k0 = ap[2 * q], k1 = ap[2 * q + 1]; const e0 = ap[2 * q + 1], e1 = ap[2 * q + 2];
    while (k0 < e0 || k1 < e1) {
      const J = Math.min(k0 < e0 ? ax[k0] : Infinity, k1 < e1 ? ax[k1] : Infinity) & ~1;
      while (k0 < e0 && ax[k0] <= J + 1) k0++; while (k1 < e1 && ax[k1] <= J + 1) k1++; nbk++;
    }
    bp[q + 1] = nbk;
  }
  const bc = new Int32Array(nbk), bv = new Float64Array(4 * nbk);
  for (let q = 0, b = 0; q < nq; q++) {
    let k0 = ap[2 * q], k1 = ap[2 * q + 1]; const e0 = ap[2 * q + 1], e1 = ap[2 * q + 2];
    while (k0 < e0 || k1 < e1) {
      const J = Math.min(k0 < e0 ? ax[k0] : Infinity, k1 < e1 ? ax[k1] : Infinity) & ~1;
      bc[b] = J;
      while (k0 < e0 && ax[k0] <= J + 1) { bv[4 * b + (ax[k0] - J)] = ad[k0]; k0++; }
      while (k1 < e1 && ax[k1] <= J + 1) { bv[4 * b + 2 + (ax[k1] - J)] = ad[k1]; k1++; }
      b++;
    }
  }
  const applyA = (x, out) => {
    let d = 0;
    for (let q = 0; q < nq; q++) {
      let s0 = 0, s1 = 0;
      for (let b = bp[q], e = bp[q + 1]; b < e; b++) { const J = bc[b], x0 = x[J], x1 = x[J + 1], v = 4 * b; s0 += bv[v] * x0; s0 += bv[v + 1] * x1; s1 += bv[v + 2] * x0; s1 += bv[v + 3] * x1; }
      const r = 2 * q, o0 = s0 + em6[r] * x[r], o1 = s1 + em6[r + 1] * x[r + 1];
      out[r] = o0; out[r + 1] = o1; d += x[r] * o0; d += x[r + 1] * o1;
    }
    return d;
  };
  applyA.par = { kind: 'bsr6', nq, bp, bc, bv, em6 };                            // the same operator for the solver threads (parallel.js)
  const jac = new Float64Array(6 * n); for (let i = 0; i < 6 * n; i++) jac[i] = diag[i] + em6[i];
  const normalise = (v) => { let s = 0; for (let i = 0; i < 6 * n; i++) s += m6[i] * v[i] * v[i]; s = Math.sqrt(s) || 1; for (let i = 0; i < 6 * n; i++) v[i] /= s; };
  let z = new Float64Array(6 * n); { let s = 777; const rnd = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff - 0.5; }; for (let i = 0; i < 6 * n; i++) z[i] = rnd(); } normalise(z);
  const rhs = new Float64Array(6 * n), Az = new Float64Array(6 * n); let iters = 0, lam = 0;
  for (let it = 0; it < outer; it++) {
    for (let i = 0; i < 6 * n; i++) rhs[i] = m6[i] * z[i];
    const res = pcg(applyA, jac, rhs, z, { maxIter: cgMax, tol: cgTol }); iters += res.iters;
    const zn = res.x; normalise(zn); applyA(zn, Az); let num = 0, ch = 0;
    for (let i = 0; i < 6 * n; i++) { num += zn[i] * Az[i]; ch += (zn[i] - z[i]) ** 2; }
    lam = num; z = zn; if (log) log(`stripes it ${it}: λ = ${lam.toExponential(3)}, Δ = ${Math.sqrt(ch).toExponential(2)}`);
    if (it > 1 && Math.sqrt(ch) < 3e-4) break;
  }
  return { psi: z, u0: u0FromPsi(z, n, m), branch: branchPoints(Fi, edge, z, !!rosy.reversed), iters, energy: lam, worstMatch: worst, edge };
}

/**
 * Lift stripe phases solved on the base carrier to a refined carrier whose rotation field (thetaFine, framesFine) has
 * already been lifted: base vertices keep ψ (same family labelling), new vertices are filled layer by layer by
 * transporting ψ across the fine edges (with the fine coupling), then a few sweeps smooth the new vertices.
 * Returns {psi, u0, branch, edge}.
 */
export function liftStripes(psiBase, nBase, cFine, thetaFine, framesFine, q0, N = 6, { sweeps = 2, reversed = false, Len = null } = {}) {
  const { Fi, Li, n } = cFine;                                                // Len: the fine carrier's lengths in ĝ (q0 is constant in ĝ)
  const { lenOf, key, m: area } = meshData(Fi, Len || Li, n);
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
  return { psi, u0: u0FromPsi(psi, n, area), branch: branchPoints(Fi, edge, psi, reversed), edge };
}
