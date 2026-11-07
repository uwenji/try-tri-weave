// Semi-implicit SH23 stepping on WebGPU: the same matrix-free Jacobi-PCG as solver_cpu.js, in WGSL compute
// kernels. Per time step: K_S, K_A, K_init, then cgIters × (K_S, K_A+dot, K_update+dot, K_p) — all scalars
// (α, β) stay on the GPU: every workgroup re-reduces the small array of per-workgroup partial sums; on large
// carriers (> 256 workgroups) a one-workgroup K_fold kernel sums the partials first, so the re-reduction does not
// cost nwg² reads. Readback (the field for display, snapshots for scrubbing) happens once per batch of steps.
export const SQ3 = Math.sqrt(3);
const WG = 256;

const WGSL = /* wgsl */`
struct Params { n: u32, nwg: u32, rzOld: u32, rzNew: u32, dt: f32, r: f32, b: f32, iq2: f32, omdr: f32, fold: u32, pad1: f32, pad2: f32 };
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> m: array<f32>;
@group(0) @binding(2) var<storage, read> indptr: array<u32>;
@group(0) @binding(3) var<storage, read> indices: array<u32>;
@group(0) @binding(4) var<storage, read> data: array<f32>;
@group(0) @binding(5) var<storage, read> x: array<f32>;
@group(0) @binding(6) var<storage, read_write> t1: array<f32>;
@group(0) @binding(7) var<storage, read_write> y: array<f32>;
@group(0) @binding(8) var<storage, read_write> part: array<f32>;

var<workgroup> sh: array<f32, ${WG}>;

fn lap_x(i: u32) -> f32 {         // L applied to x
  var s = 0.0;
  for (var k = indptr[i]; k < indptr[i + 1u]; k++) { s += data[k] * x[indices[k]]; }
  return s;
}
fn lap_rw(i: u32) -> f32 {        // L applied to t1 (read_write buffer)
  var s = 0.0;
  for (var k = indptr[i]; k < indptr[i + 1u]; k++) { s += data[k] * t1[indices[k]]; }
  return s;
}
fn reduce_part(lid: u32, v: f32, wg: u32) {
  sh[lid] = v; workgroupBarrier();
  for (var s = ${WG / 2}u; s > 0u; s = s >> 1u) { if (lid < s) { sh[lid] += sh[lid + s]; } workgroupBarrier(); }
  if (lid == 0u) { part[P.nwg + wg] = sh[0]; }   // slot 1: pAp
}

// t1 = S x / m = x − iq2 · (L x) / m
@compute @workgroup_size(${WG}) fn k_S(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x; if (i >= P.n) { return; }
  t1[i] = x[i] - P.iq2 * lap_x(i) / m[i];
}
// y = A x = (1 − dt r) m x + dt (m t1 − iq2 L t1);  part[wg] = Σ x·y
@compute @workgroup_size(${WG}) fn k_A(@builtin(global_invocation_id) g: vec3<u32>, @builtin(local_invocation_id) l: vec3<u32>, @builtin(workgroup_id) w: vec3<u32>) {
  let i = g.x; var prod = 0.0;
  if (i < P.n) {
    let yi = P.omdr * m[i] * x[i] + P.dt * (m[i] * t1[i] - P.iq2 * lap_rw(i));
    y[i] = yi; prod = x[i] * yi;
  }
  reduce_part(l.x, prod, w.x);
}
`;

const WGSL2 = /* wgsl */`
struct Params { n: u32, nwg: u32, rzOld: u32, rzNew: u32, dt: f32, r: f32, b: f32, iq2: f32, omdr: f32, fold: u32, pad1: f32, pad2: f32 };
@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read_write> u: array<f32>;
@group(0) @binding(2) var<storage, read> Ap: array<f32>;
@group(0) @binding(3) var<storage, read_write> res: array<f32>;
@group(0) @binding(4) var<storage, read_write> z: array<f32>;
@group(0) @binding(5) var<storage, read_write> p: array<f32>;
@group(0) @binding(6) var<storage, read> diag: array<f32>;
@group(0) @binding(7) var<storage, read> m: array<f32>;
@group(0) @binding(8) var<storage, read_write> part: array<f32>;

var<workgroup> sh: array<f32, ${WG}>;
var<workgroup> tot: array<f32, 2>;

fn reduce_part(lid: u32, v: f32, slot: u32, wg: u32) {
  sh[lid] = v; workgroupBarrier();
  for (var s = ${WG / 2}u; s > 0u; s = s >> 1u) { if (lid < s) { sh[lid] += sh[lid + s]; } workgroupBarrier(); }
  if (lid == 0u) { part[slot * P.nwg + wg] = sh[0]; }
}
// every workgroup sums the nwg partials of a slot → tot[k]  (or reads the folded total when k_fold ran)
fn sum_slot(lid: u32, slot: u32, k: u32) {
  if (P.fold == 1u) { if (lid == 0u) { tot[k] = part[3u * P.nwg + slot]; } workgroupBarrier(); return; }
  var s = 0.0;
  for (var j = lid; j < P.nwg; j += ${WG}u) { s += part[slot * P.nwg + j]; }
  sh[lid] = s; workgroupBarrier();
  for (var st = ${WG / 2}u; st > 0u; st = st >> 1u) { if (lid < st) { sh[lid] += sh[lid + st]; } workgroupBarrier(); }
  if (lid == 0u) { tot[k] = sh[0]; }
  workgroupBarrier();
}
fn fold_slot(lid: u32, slot: u32) {
  var s = 0.0;
  for (var j = lid; j < P.nwg; j += ${WG}u) { s += part[slot * P.nwg + j]; }
  sh[lid] = s; workgroupBarrier();
  for (var st = ${WG / 2}u; st > 0u; st = st >> 1u) { if (lid < st) { sh[lid] += sh[lid + st]; } workgroupBarrier(); }
  if (lid == 0u) { part[3u * P.nwg + slot] = sh[0]; }
  workgroupBarrier();
}
// one workgroup: totals of the three slots → part[3 nwg + slot]
@compute @workgroup_size(${WG}) fn k_fold(@builtin(local_invocation_id) l: vec3<u32>) { fold_slot(l.x, 0u); fold_slot(l.x, 1u); fold_slot(l.x, 2u); }

// res = rhs(u) − Ap; z = res/diag; p = z; part[rzOld] = Σ res·z      (Ap = A u on entry)
@compute @workgroup_size(${WG}) fn k_init(@builtin(global_invocation_id) g: vec3<u32>, @builtin(local_invocation_id) l: vec3<u32>, @builtin(workgroup_id) w: vec3<u32>) {
  let i = g.x; var prod = 0.0;
  if (i < P.n) {
    let v = u[i];
    let rhs = m[i] * (v + P.dt * (P.b * v * v - v * v * v));
    let ri = rhs - Ap[i]; let zi = ri / diag[i];
    res[i] = ri; z[i] = zi; p[i] = zi; prod = ri * zi;
  }
  reduce_part(l.x, prod, P.rzOld, w.x);
}
// α = rz/pAp; u += α p; res −= α Ap; z = res/diag; part[rzNew] = Σ res·z     (slot 1 holds pAp)
@compute @workgroup_size(${WG}) fn k_update(@builtin(global_invocation_id) g: vec3<u32>, @builtin(local_invocation_id) l: vec3<u32>, @builtin(workgroup_id) w: vec3<u32>) {
  sum_slot(l.x, P.rzOld, 0u); sum_slot(l.x, 1u, 1u);
  let pAp = tot[1]; var alpha = 0.0; if (pAp > 0.0) { alpha = tot[0] / pAp; }
  let i = g.x; var prod = 0.0;
  if (i < P.n) {
    let ui = u[i] + alpha * p[i]; let ri = res[i] - alpha * Ap[i]; let zi = ri / diag[i];
    u[i] = ui; res[i] = ri; z[i] = zi; prod = ri * zi;
  }
  reduce_part(l.x, prod, P.rzNew, w.x);
}
// β = rzNew/rzOld; p = z + β p
@compute @workgroup_size(${WG}) fn k_p(@builtin(global_invocation_id) g: vec3<u32>, @builtin(local_invocation_id) l: vec3<u32>) {
  sum_slot(l.x, P.rzNew, 0u); sum_slot(l.x, P.rzOld, 1u);
  var beta = 0.0; if (tot[1] > 1e-30) { beta = tot[0] / tot[1]; }
  let i = g.x; if (i >= P.n) { return; }
  p[i] = z[i] + beta * p[i];
}
`;

export async function getDevice() {
  if (!navigator.gpu) return null;
  try {
    const adapter = await navigator.gpu.requestAdapter({ powerPreference: 'high-performance' });
    if (!adapter) return null;
    const device = await adapter.requestDevice();
    device.info = adapter.info || {};
    return device;
  } catch (e) { console.warn('WebGPU unavailable:', e); return null; }
}

export class GpuSolver {
  constructor(device, carrier, { a, r = -0.08, b = 1.2, dt = 0.1, cgIters = 12, snapEvery = 20, dispersion = true }) {
    this.dev = device; this.c = carrier; this.n = carrier.n; this.backend = 'webgpu'; this.cgIters = cgIters; this.snapEvery = snapEvery;
    this.nwg = Math.ceil(this.n / WG); this.steps = 0; this.lastRelRes = NaN;
    this.fold = this.nwg > 256;                                                 // large carriers: fold the partials once per reduction
    const dev = device, n = this.n;
    const mk = (size, usage, label) => dev.createBuffer({ size: Math.max(16, Math.ceil(size / 4) * 4), usage, label });
    const ST = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC;
    const L = carrier.L;
    this.bufM = mk(n * 4, ST, 'm'); dev.queue.writeBuffer(this.bufM, 0, Float32Array.from(carrier.m));
    this.bufIndptr = mk((n + 1) * 4, ST, 'indptr'); dev.queue.writeBuffer(this.bufIndptr, 0, Uint32Array.from(L.indptr));
    this.bufIndices = mk(L.indices.length * 4, ST, 'indices'); dev.queue.writeBuffer(this.bufIndices, 0, Uint32Array.from(L.indices));
    this.bufData = mk(L.data.length * 4, ST, 'data'); dev.queue.writeBuffer(this.bufData, 0, Float32Array.from(L.data));
    for (const nm of ['u', 't1', 'Ap', 'res', 'z', 'p', 'diag']) this['buf' + nm] = mk(n * 4, ST, nm);
    this.bufPart = mk((3 * this.nwg + 4) * 4, ST, 'part');
    this.paramA = dev.createBuffer({ size: 48, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    this.paramB = dev.createBuffer({ size: 48, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
    this.staging = []; this.uHost = new Float32Array(n); this.snapshots = [new Float32Array(n)];
    const mod1 = dev.createShaderModule({ code: WGSL }), mod2 = dev.createShaderModule({ code: WGSL2 });
    // explicit layouts: every binding present, so bind groups can be shared between kernels that use subsets
    const lay = (types) => dev.createBindGroupLayout({ entries: types.map((t, i) => ({ binding: i, visibility: GPUShaderStage.COMPUTE, buffer: { type: t } })) });
    const RO = 'read-only-storage', RW = 'storage';
    this.bgl1 = lay(['uniform', RO, RO, RO, RO, RO, RW, RW, RW]);
    this.bgl2 = lay(['uniform', RW, RO, RW, RW, RW, RO, RO, RW]);
    const pipe = (mod, entry, bgl) => dev.createComputePipeline({ layout: dev.createPipelineLayout({ bindGroupLayouts: [bgl] }), compute: { module: mod, entryPoint: entry } });
    this.pS = pipe(mod1, 'k_S', this.bgl1); this.pA = pipe(mod1, 'k_A', this.bgl1);
    this.pInit = pipe(mod2, 'k_init', this.bgl2); this.pUpd = pipe(mod2, 'k_update', this.bgl2); this.pP = pipe(mod2, 'k_p', this.bgl2); this.pFold = pipe(mod2, 'k_fold', this.bgl2);
    const bg = (pl, entries) => dev.createBindGroup({ layout: pl === this.pS || pl === this.pA ? this.bgl1 : this.bgl2, entries: entries.map((buf, i) => ({ binding: i, resource: { buffer: buf } })) });
    const common1 = (x, param) => [param, this.bufM, this.bufIndptr, this.bufIndices, this.bufData, x, this.buft1, this.bufAp, this.bufPart];
    this.bgS_u = bg(this.pS, common1(this.bufu, this.paramA)); this.bgS_p = bg(this.pS, common1(this.bufp, this.paramA));
    this.bgA_u = bg(this.pA, common1(this.bufu, this.paramA)); this.bgA_p = bg(this.pA, common1(this.bufp, this.paramA));
    const common2 = (param) => [param, this.bufu, this.bufAp, this.bufres, this.bufz, this.bufp, this.bufdiag, this.bufM, this.bufPart];
    this.bgInit = bg(this.pInit, common2(this.paramA));
    this.bgUpdA = bg(this.pUpd, common2(this.paramA)); this.bgUpdB = bg(this.pUpd, common2(this.paramB));
    this.bgPA = bg(this.pP, common2(this.paramA)); this.bgPB = bg(this.pP, common2(this.paramB));
    this.setParams({ a, r, b, dt, dispersion });
  }

  setParams({ a = this.a, r = this.r, b = this.b, dt = this.dt, dispersion = this.dispersion ?? true } = {}) {
    this.a = a; this.r = r; this.b = b; this.dt = dt; this.dispersion = dispersion; this.q0 = 4 * Math.PI / (SQ3 * a);
    this.disp = dispersion && this.c.h ? 1 - (this.q0 * this.c.h) ** 2 / 24 : 1;  // mesh dispersion correction, see solver_cpu.js
    this.iq2 = 1 / (this.q0 * this.q0 * this.disp);
    const { L, m } = this.c, n = this.n, diag = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      let s = 0;
      for (let k = L.indptr[i]; k < L.indptr[i + 1]; k++) { const j = L.indices[k]; const Sij = (j === i ? m[i] : 0) - L.data[k] * this.iq2; s += Sij * Sij / m[j]; }
      diag[i] = (1 - dt * r) * m[i] + dt * s;
    }
    this.dev.queue.writeBuffer(this.bufdiag, 0, diag);
    const pk = (rzOld, rzNew) => { const b32 = new ArrayBuffer(48), u32 = new Uint32Array(b32), f32 = new Float32Array(b32); u32[0] = n; u32[1] = this.nwg; u32[2] = rzOld; u32[3] = rzNew; f32[4] = dt; f32[5] = r; f32[6] = b; f32[7] = this.iq2; f32[8] = 1 - dt * r; u32[9] = this.fold ? 1 : 0; return b32; };
    // slot 0 and 2 alternate as rzOld / rzNew; slot 1 holds pAp
    this.dev.queue.writeBuffer(this.paramA, 0, pk(0, 2)); this.dev.queue.writeBuffer(this.paramB, 0, pk(2, 0));
  }

  setInitial(u0) {
    this.uHost = Float32Array.from(u0); this.dev.queue.writeBuffer(this.bufu, 0, this.uHost);
    this.steps = 0; this.snapshots = [Float32Array.from(this.uHost)];
  }
  /** overwrite the field without touching the step count (heal.js: nucleated sites, pins) */
  setU(u) { this.uHost = Float32Array.from(u); this.dev.queue.writeBuffer(this.bufu, 0, this.uHost); }

  _encodeStep(pass) {
    const d = (pipe, bgp) => { pass.setPipeline(pipe); pass.setBindGroup(0, bgp); pass.dispatchWorkgroups(this.nwg); };
    d(this.pS, this.bgS_u); d(this.pA, this.bgA_u); d(this.pInit, this.bgInit);
    const fold = () => { if (this.fold) { pass.setPipeline(this.pFold); pass.setBindGroup(0, this.bgInit); pass.dispatchWorkgroups(1); } };
    for (let it = 0; it < this.cgIters; it++) {
      const even = (it & 1) === 0;                                   // rzOld in slot 0 on even iterations
      d(this.pS, this.bgS_p); d(this.pA, this.bgA_p); fold();
      d(this.pUpd, even ? this.bgUpdA : this.bgUpdB); fold(); d(this.pP, even ? this.bgPA : this.bgPB);
    }
  }

  _stagingBuffer(k) {
    while (this.staging.length <= k) this.staging.push(this.dev.createBuffer({ size: this.n * 4, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST }));
    return this.staging[k];
  }

  /** advance nsteps; resolves when uHost holds the new field (and snapshots are recorded) */
  async step(nsteps) {
    if (this.disposed) throw new Error('solver disposed');
    this.inflight = (this.inflight || 0) + 1;
    try { return await this._step(nsteps); } finally { this.inflight--; if (this.disposed && !this.inflight) this._destroy(); }
  }
  async _step(nsteps) {
    const enc = this.dev.createCommandEncoder();
    const snaps = []; let done = 0, k = 0;
    while (done < nsteps) {
      const toSnap = this.snapEvery - (this.steps % this.snapEvery);
      const chunk = Math.min(toSnap, nsteps - done);
      const pass = enc.beginComputePass();
      for (let s = 0; s < chunk; s++) this._encodeStep(pass);
      pass.end();
      this.steps += chunk; done += chunk;
      if (this.steps % this.snapEvery === 0 || done === nsteps) {
        const st = this._stagingBuffer(k++); enc.copyBufferToBuffer(this.bufu, 0, st, 0, this.n * 4);
        snaps.push({ st, isSnap: this.steps % this.snapEvery === 0 });
      }
    }
    this.dev.queue.submit([enc.finish()]);
    for (const { st, isSnap } of snaps) {
      await st.mapAsync(GPUMapMode.READ);
      const arr = new Float32Array(st.getMappedRange().slice(0)); st.unmap();
      if (isSnap) this.snapshots.push(arr);
      this.uHost = arr;
    }
    return this.uHost;
  }

  get u() { return this.uHost; }
  stateAt(step) {
    if (step >= this.steps) return this.uHost;
    const k = Math.min(Math.floor(step / this.snapEvery), this.snapshots.length - 1);
    return this.snapshots[k];
  }
  get t() { return this.steps * this.dt; }
  async sync() { return this.uHost; }
  /** release the GPU buffers — after a step still in flight has read its result (destroying a mapped staging buffer
   *  rejects that step, which used to clamp the app's target to the old solver's count when a carrier or seed changed) */
  dispose() { this.disposed = true; if (!this.inflight) this._destroy(); }
  _destroy() { for (const b of [this.bufM, this.bufIndptr, this.bufIndices, this.bufData, this.bufu, this.buft1, this.bufAp, this.bufres, this.bufz, this.bufp, this.bufdiag, this.bufPart, this.paramA, this.paramB, ...this.staging]) try { b.destroy(); } catch (e) { /* ignore */ } }
}
