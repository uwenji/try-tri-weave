// Semi-implicit SH23 stepping on the CPU (typed arrays), matrix-free Jacobi-PCG.
//   [(1 − dt r) M + dt S M⁻¹ S] u⁺ = M (u + dt (b u² − u³)),   S = M − L / q0²
// The GPU solver (solver_gpu.js) implements the same interface with WGSL kernels.
import { spmv, dot } from './sparse.js';

export const SQ3 = Math.sqrt(3);

export class CpuSolver {
  constructor(carrier, { a, r = -0.08, b = 1.2, dt = 0.1, cgIters = 12, snapEvery = 20, dispersion = true }) {
    this.c = carrier; this.n = carrier.n; this.backend = 'cpu';
    this.cgIters = cgIters; this.snapEvery = snapEvery;
    this.u = new Float64Array(this.n); this.steps = 0; this.snapshots = [new Float32Array(this.n)];
    this.tmp = new Float64Array(this.n); this.t1 = new Float64Array(this.n); this.t2 = new Float64Array(this.n);
    this.rhs = new Float64Array(this.n); this.res = new Float64Array(this.n); this.z = new Float64Array(this.n);
    this.p = new Float64Array(this.n); this.Ap = new Float64Array(this.n);
    this.setParams({ a, r, b, dt, dispersion });
    this.lastRelRes = 0;
  }

  setParams({ a = this.a, r = this.r, b = this.b, dt = this.dt, dispersion = this.dispersion ?? true } = {}) {
    this.a = a; this.r = r; this.b = b; this.dt = dt; this.dispersion = dispersion; this.q0 = 4 * Math.PI / (SQ3 * a);
    // Dispersion correction: the cotan Laplacian's symbol on a mesh of edge h is λ(k) ≈ k² (1 − c k²h²) — c = 1/16 on an
    // equilateral grid, measured c ≈ 1/24 (h = median edge) on the isotropic remeshed carriers — so the operator (1 + L/q0²)²
    // selects k = q0 (1 + c (q0 h)²/2): an edge of 0.981 a at 8 edges per a, 0.966 a at 6 (flat disk, test_dispersion.mjs).
    // Replacing q0² by λ(q0) makes the selected wavenumber q0 itself (edge = a ± 0.5 %). The continuum SH23 ground state
    // itself sits at q0 to +0.3 % (experiments/calibrate.js, ComFiT-style strain minimisation).
    this.disp = dispersion && this.c.h ? 1 - (this.q0 * this.c.h) ** 2 / 24 : 1;
    this.iq2 = 1 / (this.q0 * this.q0 * this.disp);
    // Jacobi diagonal of A = (1 − dt r) M + dt S M⁻¹ S
    const { L, m, Ldiag } = this.c, n = this.n, diag = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let s = 0;
      for (let k = L.indptr[i]; k < L.indptr[i + 1]; k++) {
        const j = L.indices[k]; const Sij = (j === i ? m[i] : 0) - L.data[k] * this.iq2; s += Sij * Sij / m[j];
      }
      diag[i] = (1 - dt * r) * m[i] + dt * s;
    }
    this.diag = diag;
  }

  /** y = S x = m∘x − L x / q0² */
  applyS(x, y) { spmv(this.c.L, x, this.tmp); const m = this.c.m; for (let i = 0; i < this.n; i++) y[i] = m[i] * x[i] - this.iq2 * this.tmp[i]; }

  applyA(x, out) {
    const m = this.c.m, n = this.n, { dt, r } = this;
    this.applyS(x, this.t1); for (let i = 0; i < n; i++) this.t1[i] /= m[i];
    this.applyS(this.t1, this.t2);
    for (let i = 0; i < n; i++) out[i] = (1 - dt * r) * m[i] * x[i] + dt * this.t2[i];
  }

  setInitial(u0) { this.u = Float64Array.from(u0); this.steps = 0; this.snapshots = [Float32Array.from(u0)]; }

  step(nsteps) {
    const { n, u, rhs, res, z, p, Ap, diag } = this; const m = this.c.m, { dt, b } = this;
    for (let s = 0; s < nsteps; s++) {
      for (let i = 0; i < n; i++) { const v = u[i]; rhs[i] = m[i] * (v + dt * (b * v * v - v * v * v)); }
      // warm-started PCG from the current u
      this.applyA(u, Ap);
      for (let i = 0; i < n; i++) { res[i] = rhs[i] - Ap[i]; z[i] = res[i] / diag[i]; p[i] = z[i]; }
      let rz = dot(res, z);
      for (let it = 0; it < this.cgIters; it++) {
        this.applyA(p, Ap);
        const pAp = dot(p, Ap); if (!(pAp > 0)) break;
        const alpha = rz / pAp;
        for (let i = 0; i < n; i++) { u[i] += alpha * p[i]; res[i] -= alpha * Ap[i]; z[i] = res[i] / diag[i]; }
        const rzNew = dot(res, z), beta = rzNew / rz; rz = rzNew;
        if (rz < 1e-30) break;
        for (let i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
      }
      this.lastRelRes = Math.sqrt(Math.max(rz, 0)) / (Math.sqrt(dot(rhs, rhs)) + 1e-30);
      this.steps++;
      if (this.steps % this.snapEvery === 0) this.snapshots.push(Float32Array.from(u));
    }
    return Promise.resolve(this.u);
  }

  /** field at a recorded step (nearest snapshot at or below), or the live state */
  stateAt(step) {
    if (step >= this.steps) return this.u;
    const k = Math.min(Math.floor(step / this.snapEvery), this.snapshots.length - 1);
    return this.snapshots[k];
  }

  get t() { return this.steps * this.dt; }
  async sync() { return this.u; }
  dispose() {}
}
