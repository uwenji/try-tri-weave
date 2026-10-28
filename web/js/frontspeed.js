// Front-speed model of the subcritical SH23 pattern and the growth prediction drawn over the tracked read-outs.
//   flat table v(r, b)  →  circular front dR/dt = (v − D/(q0 R)) / q0  →  covered area A(R) on the carrier (geodesic
//   distance from the centre)  →  predicted peaks = A(R) / (√3⁄2 · (s a)²), interior = A(R − 0.85 a) / site, with s the
//   realised lattice edge / a: 1 with the carrier dispersion correction (solver_cpu.js setParams), else ≈ 1/(1 + (q0 h)²/48).
import { FRONT, BOUND } from './frontspeed_table.js';

const SQ3 = Math.sqrt(3);
export const SITE_SCALE = 1;                       // realised edge / a with the dispersion correction on (zoo, closed shapes: peaks / sites = 0.99); p.siteScale overrides
const INTERIOR_LAG = 0.85;                         // interior peaks (closed link) trail the outermost peaks by ≈ 0.85 a

export function dtFactor(dt) {
  const xs = FRONT.dt, ys = FRONT.dtRel[0];
  const at = (x) => { if (x <= xs[0]) return ys[0]; for (let i = 1; i < xs.length; i++) if (x <= xs[i]) { const t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[i - 1] + t * (ys[i] - ys[i - 1]); } const n = xs.length; return ys[n - 1] + (x - xs[n - 1]) * (ys[n - 1] - ys[n - 2]) / (xs[n - 1] - xs[n - 2]); };
  return at(dt) / at(0.1);                         // the table was measured at dt = 0.1
}

// monotone cubic (Fritsch–Carlson, as Matlab's pchip) through (xs, ys), constant outside [xs[0], xs[n − 1]]
function pchip(xs, ys) {
  const n = xs.length, h = [], d = [], m = new Array(n);
  for (let i = 0; i < n - 1; i++) { h.push(xs[i + 1] - xs[i]); d.push((ys[i + 1] - ys[i]) / h[i]); }
  m[0] = d[0]; m[n - 1] = d[n - 2];
  for (let i = 1; i < n - 1; i++) {
    if (d[i - 1] * d[i] <= 0) { m[i] = 0; continue; }
    const w1 = 2 * h[i] + h[i - 1], w2 = h[i] + 2 * h[i - 1]; m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i]);
  }
  return (x) => {
    if (x <= xs[0]) return ys[0]; if (x >= xs[n - 1]) return ys[n - 1];
    let i = 0; while (i < n - 2 && xs[i + 1] < x) i++;
    const t = (x - xs[i]) / h[i], t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * ys[i] + (t3 - 2 * t2 + t) * h[i] * m[i] + (3 * t2 - 2 * t3) * ys[i + 1] + (t3 - t2) * h[i] * m[i + 1];
  };
}

/** the two measured edges (frontspeed_table.js BOUND): below rDie(b) the seed dies, above rGrow(b) the front advances */
export const rDie = pchip(BOUND.b, BOUND.die), rGrow = pchip(BOUND.b, BOUND.grow);

// above the growth edge: each column of the table as v against x = r − r_grow(b) (0.01 at the edge), and between two
// columns the speed at equal distance from the edge — the speed follows the edge instead of being smeared across it
const PROF = FRONT.b.map((bj, j) => {
  const g = rGrow(bj), xs = [0], vs = [0.01];
  FRONT.r.forEach((ri, i) => { const v = FRONT.v[i][j]; if (v !== null && v > 0.01 && ri - g > 1e-3) { xs.push(ri - g); vs.push(v); } });
  return { xs, vs };
});
function profAt({ xs, vs }, x) {
  const n = xs.length; if (n < 2) return vs[0];
  let i = 0; while (i < n - 2 && xs[i + 1] < x) i++;
  return vs[i] + (x - xs[i]) * (vs[i + 1] - vs[i]) / (xs[i + 1] - xs[i]);                // linear; past the last row, extrapolated
}
const clampB = (b) => Math.max(BOUND.b[0], Math.min(BOUND.b[BOUND.b.length - 1], b));

/** front speed (q0 units, dt = 0.1) against r at a fixed b, for r ≥ rGrow(b) — one column of the map */
export function speedColumn(b) {
  const B = FRONT.b, bc = clampB(b);
  let j = 0; while (j < B.length - 2 && B[j + 1] < bc) j++;
  const t = (bc - B[j]) / (B[j + 1] - B[j]), g = rGrow(bc), P = PROF[j], Q = PROF[j + 1];
  return (r) => Math.max(0.01, (1 - t) * profAt(P, r - g) + t * profAt(Q, r - g));
}

/** flat front speed (q0 units, dt = 0.1) at (r, b): null where the seeded patch dies, 0 where it is pinned, else ≥ 0.01 */
export function tableSpeed(r, b) {
  const bc = clampB(b);
  if (r < rDie(bc)) return null;
  if (r < rGrow(bc)) return 0;
  return speedColumn(bc)(r);
}

/** front state at (r, b, dt): v in q0 units, in lattice edges per 100 steps, D, critical radius in a, and a verdict */
export function frontState(r, b, dt = 0.1) {
  const v0 = tableSpeed(r, b);
  if (v0 === null) return { state: 'dies', v: null, vA100: 0, D: 0, RcA: Infinity, label: 'the seed dies out (u = 0 wins)' };
  const v = v0 * dtFactor(dt), D = 0.8 + 3.4 * Math.max(v, 0);
  const toA = SQ3 / (4 * Math.PI);                 // q0-units of length → lattice edges a  (q0 a = 4π/√3)
  const vA100 = v * toA * dt * 100;
  const RcA = v > 0.01 ? D / v * toA : Infinity;
  let state = 'advance', label = 'the front advances';
  if (v0 < 0.01) { state = 'pinned'; label = 'pinned: the patch survives but does not grow'; }
  if (r > 0) label += ' · r > 0: u = 0 is unstable, patterns also appear away from the front';
  return { state, v, vA100, D, RcA, label };
}

/** Prediction for the growth strip.
 *  geo = { dist (Float64Array, geodesic distance from the centre), m (vertex areas), K (angle defects), boundary }
 *  p   = { a, seedR, seedKind, valence, r, b, dt, xmax, every, siteScale } */
export function predictGrowth(geo, p) {
  const { dist, m, K, boundary, rim } = geo, n = dist.length, a = p.a;
  const order = Array.from({ length: n }, (_, i) => i).sort((i, j) => dist[i] - dist[j]);
  const cumA = new Float64Array(n), cumI = new Float64Array(n), cumK = new Float64Array(n); let sa = 0, si = 0, sk = 0;
  for (let k = 0; k < n; k++) {
    const i = order[k]; sa += m[i]; if (!boundary[i]) sk += K[i];
    if (!rim || rim[i] >= INTERIOR_LAG * a) si += m[i];                     // peaks next to the rim have an open link: not interior
    cumA[k] = sa; cumI[k] = si; cumK[k] = sk;
  }
  const total = sa, Rmax = dist[order[n - 1]];
  const at = (arr, R) => { if (R <= 0) return 0; let lo = 0, hi = n; while (lo < hi) { const mid = (lo + hi) >> 1; if (dist[order[mid]] < R) lo = mid + 1; else hi = mid; } return lo ? arr[lo - 1] : 0; };
  const q0 = 4 * Math.PI / (SQ3 * a), site = SQ3 / 2 * ((p.siteScale ?? SITE_SCALE) * a) ** 2;
  const seedCharge = 6 - p.valence;
  const grow = frontState(p.r, p.b, p.dt);
  let R = (p.seedKind === 'spots' ? 1.3 : 0.82 * p.seedR) * a;
  const every = p.every || 50, xmax = p.xmax;
  const out = { steps: [], peaks: [], interior: [], floor: [], R: [], coveredAt: -1, total: total / site, totalInterior: si / site, grow, Rmax };
  let fullAt = -1; const st = grow;
  const sub = 10, dstep = every / sub;
  for (let s = 0; s <= xmax; s += every) {
    if (s > 0) {
      for (let k = 0; k < sub; k++) {
        if (st.v === null) { R = Math.max(0, R - 0.002 * a * dstep); continue; }       // dying: shrink (only qualitative)
        if (st.state === 'pinned') continue;
        const vq = st.v > 0 ? st.v - st.D / Math.max(R * q0, 1e-6) : st.v;   // q0 units; below R_c a patch shrinks
        R = Math.max(0, R + vq / q0 * p.dt * dstep);
      }
    }
    const A = Math.min(total, at(cumA, R)), Ai = at(cumI, R - INTERIOR_LAG * a);
    const omega = at(cumK, Math.max(R - INTERIOR_LAG * a, 0)) * 180 / Math.PI;
    const need = Math.round(Math.abs(omega / 60 - seedCharge));
    out.steps.push(s); out.peaks.push(A / site); out.interior.push(Ai / site); out.R.push(R / a);
    out.floor.push((seedCharge !== 0 ? 1 : 0) + need);
    if (fullAt < 0 && A >= 0.995 * total) { fullAt = s; out.coveredAt = s; }
  }
  return out;
}
