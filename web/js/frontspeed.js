// Front-speed model of the subcritical SH23 pattern and the growth prediction drawn over the tracked read-outs.
//   flat table v(r, b)  →  circular front dR/dt = (v − D/(q0 R)) / q0  →  covered area A(R) on the carrier (geodesic
//   distance from the centre)  →  predicted peaks = A(R) / (√3⁄2 · (s a)²), interior = A(R − 0.85 a) / site, with s the
//   realised lattice edge / a: 1 with the carrier dispersion correction (solver_cpu.js setParams), else ≈ 1/(1 + (q0 h)²/48).
import { FRONT } from './frontspeed_table.js';

const SQ3 = Math.sqrt(3);
export const SITE_SCALE = 1;                       // realised edge / a with the dispersion correction on (zoo, closed shapes: peaks / sites = 0.99); p.siteScale overrides
const INTERIOR_LAG = 0.85;                         // interior peaks (closed link) trail the outermost peaks by ≈ 0.85 a

function dtFactor(dt) {
  const xs = FRONT.dt, ys = FRONT.dtRel[0];
  const at = (x) => { if (x <= xs[0]) return ys[0]; for (let i = 1; i < xs.length; i++) if (x <= xs[i]) { const t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[i - 1] + t * (ys[i] - ys[i - 1]); } const n = xs.length; return ys[n - 1] + (x - xs[n - 1]) * (ys[n - 1] - ys[n - 2]) / (xs[n - 1] - xs[n - 2]); };
  return at(dt) / at(0.1);                         // the table was measured at dt = 0.1
}

/** bilinear look-up of the flat front speed (q0 units) at (r, b); null when the seeded patch dies there */
export function tableSpeed(r, b) {
  const R = FRONT.r, B = FRONT.b, V = FRONT.v;
  const rc = Math.max(R[0], Math.min(R[R.length - 1], r)), bc = Math.max(B[0], Math.min(B[B.length - 1], b));
  let i = 0; while (i < R.length - 2 && R[i + 1] < rc) i++;
  let j = 0; while (j < B.length - 2 && B[j + 1] < bc) j++;
  const tr = (rc - R[i]) / (R[i + 1] - R[i]), tb = (bc - B[j]) / (B[j + 1] - B[j]);
  const c = [[V[i][j], V[i][j + 1]], [V[i + 1][j], V[i + 1][j + 1]]], w = [[(1 - tr) * (1 - tb), (1 - tr) * tb], [tr * (1 - tb), tr * tb]];
  let s = 0, ws = 0, wn = 0;
  for (let p = 0; p < 2; p++) for (let q = 0; q < 2; q++) { if (c[p][q] === null) { wn += w[p][q]; continue; } s += w[p][q] * c[p][q]; ws += w[p][q]; }
  if (wn > 0.5 || ws === 0) return null;
  return s / ws;
}

/** smooth field for drawing the map: v where the patch lives, −0.15 where it dies (bilinear, no holes) */
export function mapValue(r, b) {
  const R = FRONT.r, B = FRONT.b, V = FRONT.v;
  const rc = Math.max(R[0], Math.min(R[R.length - 1], r)), bc = Math.max(B[0], Math.min(B[B.length - 1], b));
  let i = 0; while (i < R.length - 2 && R[i + 1] < rc) i++;
  let j = 0; while (j < B.length - 2 && B[j + 1] < bc) j++;
  const tr = (rc - R[i]) / (R[i + 1] - R[i]), tb = (bc - B[j]) / (B[j + 1] - B[j]), g = (p, q) => (V[p][q] === null ? -0.15 : V[p][q]);
  return (1 - tr) * ((1 - tb) * g(i, j) + tb * g(i, j + 1)) + tr * ((1 - tb) * g(i + 1, j) + tb * g(i + 1, j + 1));
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
  if (Math.abs(v0) < 0.01) { state = 'pinned'; label = 'pinned: the patch survives but does not grow'; }
  else if (v0 < 0) { state = 'retreat'; label = 'the front retreats'; }
  if (r > 0) label += ' · r > 0: u = 0 is unstable, patterns also appear away from the front';
  return { state, v, vA100, D, RcA, label };
}

/** Prediction for the growth strip.
 *  geo = { dist (Float64Array, geodesic distance from the centre), m (vertex areas), K (angle defects), boundary }
 *  p   = { a, seedR, seedKind, valence, r, b, dt, anneal: {on, r2, b2, when: 'auto'|'step', step, startedAt}, xmax, every, siteScale } */
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
  const grow = frontState(p.r, p.b, p.dt), ann = p.anneal && p.anneal.on ? frontState(p.anneal.r2, p.anneal.b2, p.dt) : null;
  let R = (p.seedKind === 'spots' ? 1.3 : 0.82 * p.seedR) * a;
  const every = p.every || 50, xmax = p.xmax;
  const out = { steps: [], peaks: [], interior: [], floor: [], R: [], coveredAt: -1, switchAt: -1, total: total / site, totalInterior: si / site, grow, ann, Rmax };
  let switchAt = ann && p.anneal.startedAt >= 0 ? p.anneal.startedAt : (ann && p.anneal.when === 'step' ? p.anneal.step : -1);
  let fullAt = -1, st = grow;
  const sub = 10, dstep = every / sub;
  for (let s = 0; s <= xmax; s += every) {
    if (s > 0) {
      for (let k = 0; k < sub; k++) {
        const stepNow = s - every + k * dstep;
        st = ann && switchAt >= 0 && stepNow >= switchAt ? ann : grow;
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
    if (fullAt < 0 && A >= 0.995 * total) { fullAt = s; out.coveredAt = s; if (ann && switchAt < 0 && p.anneal.when === 'auto') switchAt = s + 300; }
  }
  out.switchAt = switchAt;
  return out;
}
