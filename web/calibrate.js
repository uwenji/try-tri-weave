// Equilibrium strain of the SH23 hexagonal ground state (after ComFiT's conf_strain_to_equilibrium):
// relax the pattern on a periodic rectangular cell of the triangular lattice, evaluate the free-energy density
// F = ⟨ ½ ((1 + Δ/q0²) u)² − ½ r u² − (b/3) u³ + ¼ u⁴ ⟩ as a function of a uniform strain ε of the cell (period
// a (1+ε) × √3 a (1+ε)), and minimise over ε. The realised lattice edge of SH23 with operator q0 = 4π/(√3 a) is
// a (1+ε*); to get an edge of exactly a, run the operator at a_op = a / (1+ε*).
//
//   import { equilibriumStrain } from './calibrate.js';
//   const { eps, F, samples } = equilibriumStrain(-0.08, 1.2);      // eps ≈ equilibrium strain, F its energy density
//
// Units: q0 = 1 (a = 4π/√3). Spectral semi-implicit stepping on an Nx × Ny grid (power-of-two FFT).
const SQ3 = Math.sqrt(3);

function fft1(re, im, inv) {                                  // in-place iterative radix-2 FFT
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit;
    if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = 2 * Math.PI / len * (inv ? 1 : -1), wr = Math.cos(ang), wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let cr = 1, ci = 0;
      for (let k = 0; k < len / 2; k++) {
        const a = i + k, b = a + len / 2;
        const tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
        re[b] = re[a] - tr; im[b] = im[a] - ti; re[a] += tr; im[a] += ti;
        const nr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = nr;
      }
    }
  }
  if (inv) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
}
function fft2(re, im, nx, ny, inv) {                           // re, im: nx*ny row-major (x slow, y fast)
  const r = new Float64Array(Math.max(nx, ny)), m = new Float64Array(Math.max(nx, ny));
  for (let x = 0; x < nx; x++) {                              // rows (along y)
    const rr = re.subarray(x * ny, (x + 1) * ny), mm = im.subarray(x * ny, (x + 1) * ny); fft1(rr, mm, inv);
  }
  const rc = r.subarray(0, nx), mc = m.subarray(0, nx);
  for (let y = 0; y < ny; y++) {                              // columns (along x)
    for (let x = 0; x < nx; x++) { rc[x] = re[x * ny + y]; mc[x] = im[x * ny + y]; }
    fft1(rc, mc, inv);
    for (let x = 0; x < nx; x++) { re[x * ny + y] = rc[x]; im[x * ny + y] = mc[x]; }
  }
}

/** Free-energy density of the relaxed SH23 pattern on the cell a(1+eps) × √3 a(1+eps). */
export function relaxedEnergy(r, b, eps, { nx = 16, ny = 32, steps = 1500, dt = 0.5, q0 = 1 } = {}) {
  const a = 4 * Math.PI / (SQ3 * q0), Lx = a * (1 + eps), Ly = SQ3 * a * (1 + eps), n = nx * ny;
  // k² on the grid and the linear symbol
  const k2 = new Float64Array(n), lin = new Float64Array(n);
  for (let x = 0; x < nx; x++) {
    const kx = 2 * Math.PI * (x < nx / 2 ? x : x - nx) / Lx;
    for (let y = 0; y < ny; y++) { const ky = 2 * Math.PI * (y < ny / 2 ? y : y - ny) / Ly; const kk = kx * kx + ky * ky; k2[x * ny + y] = kk; const s = 1 - kk / (q0 * q0); lin[x * ny + y] = s * s; }
  }
  // initial condition: ideal hexagons with the cell's own wavevectors (so the pattern is periodic on the cell for any eps)
  const u = new Float64Array(n), G = [[0, 1], [-SQ3 / 2, -0.5], [SQ3 / 2, -0.5]];
  const qx = 2 * Math.PI * 2 / Lx, qy = 2 * Math.PI * 2 / Ly;   // hex wavevectors on this cell: q1 = (0, 2·2π/Ly) etc.
  const q = [[0, 2 * 2 * Math.PI / Ly], [-2 * Math.PI / Lx, -2 * Math.PI / Ly], [2 * Math.PI / Lx, -2 * Math.PI / Ly]]; void qx; void qy; void G;
  for (let x = 0; x < nx; x++) for (let y = 0; y < ny; y++) {
    const px = Lx * x / nx, py = Ly * y / ny; let s = 0;
    for (const [gx, gy] of q) s += Math.cos(gx * px + gy * py);
    u[x * ny + y] = 0.5 * s;
  }
  const re = new Float64Array(n), im = new Float64Array(n), nr = new Float64Array(n), ni = new Float64Array(n);
  re.set(u); fft2(re, im, nx, ny, false);
  for (let it = 0; it < steps; it++) {
    for (let i = 0; i < n; i++) { const v = u[i]; nr[i] = b * v * v - v * v * v; ni[i] = 0; }
    fft2(nr, ni, nx, ny, false);
    for (let i = 0; i < n; i++) { const d = 1 + dt * (lin[i] - r); re[i] = (re[i] + dt * nr[i]) / d; im[i] = (im[i] + dt * ni[i]) / d; }
    nr.set(re); ni.set(im); fft2(nr, ni, nx, ny, true); u.set(nr);
  }
  // energy density: ½ (L u)² − ½ r u² − b/3 u³ + ¼ u⁴,  L u = ifft((1 − k²/q0²) û)
  nr.set(re); ni.set(im); for (let i = 0; i < n; i++) { const s = 1 - k2[i] / (q0 * q0); nr[i] *= s; ni[i] *= s; }
  fft2(nr, ni, nx, ny, true);
  let F = 0, amp = 0;
  for (let i = 0; i < n; i++) { const v = u[i], l = nr[i]; F += 0.5 * l * l - 0.5 * r * v * v - b / 3 * v * v * v + 0.25 * v * v * v * v; amp = Math.max(amp, Math.abs(v)); }
  return { F: F / n, amp };
}

/** Golden-section minimisation of the relaxed energy density over the strain; returns {eps, F, samples}. */
export function equilibriumStrain(r, b, { lo = -0.08, hi = 0.08, tol = 2e-4, opts = {} } = {}) {
  const g = (Math.sqrt(5) - 1) / 2, samples = [];
  const f = (e) => { const { F, amp } = relaxedEnergy(r, b, e, opts); samples.push([e, F, amp]); return amp < 1e-6 ? 1e9 : F; };
  let a0 = lo, b0 = hi, c = b0 - g * (b0 - a0), d = a0 + g * (b0 - a0), fc = f(c), fd = f(d);
  while (b0 - a0 > tol) {
    if (fc < fd) { b0 = d; d = c; fd = fc; c = b0 - g * (b0 - a0); fc = f(c); }
    else { a0 = c; c = d; fc = fd; d = a0 + g * (b0 - a0); fd = f(d); }
  }
  const eps = 0.5 * (a0 + b0);
  return { eps, F: f(eps), samples: samples.sort((p, q) => p[0] - q[0]) };
}
