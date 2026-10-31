// Multi-threaded CG for the solves of the global pipeline: the rotation field and its unit-length (MBO) steps
// (rosy.js) and the stripe phases (stripes.js). The worker that runs the pipeline (global_worker.js) coordinates; a
// pool of helper workers (pcg_worker.js) computes the rows and vector entries of every CG iteration. Operators and
// vectors live in one shared WebAssembly memory and every thread runs the WebAssembly kernels of pcg_kernels.c on it.
//
// Same numbers as the single-thread pcg (sparse.js, pcg_wasm.js), bit for bit: every row and every vector entry is
// computed by the same operations in the same order, and the three dot products of an iteration are summed by the
// coordinator in index order — the helpers write the products, the coordinator adds them up chunk by chunk as the
// chunks complete, so the (inherently sequential) sums run alongside the helpers' work instead of after it.
//
// Needs a cross-origin isolated page: COOP + COEP headers (multi_threading/serve.py, or Cloudflare Pages with _headers)
// or the service worker coi-sw.js (GitHub Pages). Without it startPool() returns null and every solve runs in one
// thread: the same results, slower. A solve smaller than MIN_N unknowns or an operator without a kernel description
// (`applyA.par` unset) stays single-threaded too, and so does every solve of a size where the pool measured itself
// slower than one thread in the current run (a browser at low priority, see MARGIN).
import { setParallelPCG, dot } from './sparse.js';
import { MT } from './pcg_wasm_bin.js';
import { wasmBytes, groupRefl, wasmOn } from './pcg_wasm.js';

// control block: Int32 slots shared by all threads
const GEN = 0, KERNEL = 1, NEXT = 2, NCHUNK = 3, DONE = 4, READY = 5, ERR = 6, CTL = 8;
const K_SPMV = 1, K_UPDATE = 2, K_PUPDATE = 3, K_END = 9;
const ELEMS = 2048, ROWS = 256;                 // vector entries / operator rows (row pairs for the stripes) per chunk
const MAX_CHUNKS = 1 << 16, SPIN = 1 << 15;     // helpers spin this many checks before they sleep on Atomics.wait
const PAGES0 = 1024, PAGES_MAX = 16384;          // the shared memory: 64 MB to start, at most 1 GB (pcg_mt.wasm's limit)
export const MIN_N = 8192;                       // smaller systems: the single thread is as fast
// helper threads by default. Every CG phase waits for its slowest thread, so more threads than performance cores
// lose: Bunny ×2 on an 11-core M-series Mac (Chrome 152) took 1.57 s in one thread (WebAssembly), 1.22 / 1.09 / 1.06 /
// 1.11 / 1.30 s with 2 / 3 / 4 / 5 / 8 helpers
export const HELPERS = 4;
// The pool is not always faster: while the browser runs at low priority (another app in front, a battery saver) its
// threads land on slow cores and every phase waits for the slowest (Bunny ×2 in a Chrome behind the Claude app: one
// thread 6.2 s, five 11.7 s). So the first solve of each size in a run measures: after WARM iterations on the pool it
// times PROBE more on the pool, then PROBE in the coordinator alone (the fastest of them counts). Both run the same
// operations and the same sums in the same order (the chunks of a phase only split the index range), so the numbers
// do not change. The pool stays for that size when it is faster than the coordinator alone by MARGIN; otherwise the
// rest of the solve runs alone and the later solves of that size in the one-thread WebAssembly (pcg_wasm.js). Alone,
// the coordinator runs the pool's unfused passes, a little slower than those fused kernels (Bunny ×2, whole pipeline:
// 4–13 %), so a pool that beats it by MARGIN is about as fast as one thread or faster. At low priority the speed of the
// pool comes and goes, so a pool that stays keeps proving itself: after every WINDOW iterations on it (the start of a
// solve included) their time over WINDOW iterations alone, minus 1, is added to a debt that never goes below 0. A debt
// over DEBT — the pool lost two whole windows of iterations alone, net of its gains since it last broke even — steps
// the rest aside as above; one slow window (a helper preempted for a moment) does not. The 4 test models ×2 with 2
// helpers in a 4-core VM: 0.5 to 0.66, and windows of 1 to 2.3 now and then; Bunny ×2 in a Chrome behind the Claude
// app: 1.2 to 15, and a pool measured at 0.4 that turned 2.4× slower than one thread within the same solve.
// pool.force (tests): 'pool' never measures, 'one' steps aside after the measurement, 'late' after the first window,
// 'alone' runs every solve alone.
const WARM = 1, PROBE = 4, WINDOW = 16, MARGIN = 0.85, DEBT = 2;

/** an event on a browser Worker (addEventListener, data in e.data) or a Node worker_threads Worker (on, data itself) */
const on = (w, type, f) => (w.addEventListener ? w.addEventListener(type, (e) => f(type === 'message' ? e.data : e)) : w.on(type, f));

/** one chunk of a phase, on the WebAssembly kernels (helpers and coordinator) */
function runChunk(ex, kernel, op, V, scal, c) {
  if (kernel === K_SPMV) {
    const lo = ROWS * c, hi = Math.min(op.units, lo + ROWS);
    if (op.kind === 'bsr6') ex.rows_bsr6(op.bp, op.bc, op.bv, op.em6, V.p, V.Ap, V.pr1, lo, hi);
    else ex.rows_rosy(op.n, op.ip, op.ix, op.dr, op.di, op.em, op.rp, op.rj, op.rc, op.rs, op.m, op.tau, op.mbo, V.p, V.Ap, V.pr1, lo, hi);
  } else {
    const lo = ELEMS * c, hi = Math.min(V.N, lo + ELEMS);
    if (kernel === K_UPDATE) ex.cg_update(V.x, V.p, V.r, V.Ap, V.dg, V.z, V.pr2, V.pr3, scal[0], lo, hi);
    else ex.cg_pupdate(V.p, V.z, scal[1], lo, hi);
  }
}

// ------------------------------------------------------------------ the pool (coordinator side)
export class Pool {
  /** threads: helper workers (the coordinator is one more thread); spawn(): a Worker running pcg_worker.js.
   *  Await pool.ready before the first solve: a helper started from a worker loads through that worker's event loop
   *  (Chrome), and a solve keeps the coordinator busy without yielding. */
  constructor(threads, spawn) {
    this.K = threads; this.broken = false; this.op = null; this.V = null; this.rz = 0; this.rr = 0;
    this.ctl = new Int32Array(new SharedArrayBuffer(4 * CTL));
    this.cdone = new Int32Array(new SharedArrayBuffer(4 * MAX_CHUNKS));
    this.scal = new Float64Array(new SharedArrayBuffer(16));
    this.memory = new WebAssembly.Memory({ initial: PAGES0, maximum: PAGES_MAX, shared: true });
    this.gen = 1; this.cache = new WeakMap(); this.refl = new WeakMap(); this.vecs = new Map(); this.ops = new WeakMap(); this.pick = new Map(); this.force = null;
    this.workers = [];
    for (let k = 0; k < threads; k++) this.workers.push(spawn());
    this.ready = (async () => {
      const module = await WebAssembly.compile(wasmBytes(MT));
      const inst = await WebAssembly.instantiate(module, { env: { memory: this.memory } });
      this.ex = inst.exports; this.base = this.top = (this.ex.__heap_base.value + 15) & ~15;
      await Promise.all(this.workers.map((w) => new Promise((res, rej) => {
        on(w, 'message', (d) => { if (d && d.type === 'ready') res(); });
        on(w, 'error', (e) => rej(new Error((e && e.message) || 'a solver thread did not start')));
        w.postMessage({ type: 'init', ctl: this.ctl, cdone: this.cdone, scal: this.scal, module, memory: this.memory });
      })));
    })();
  }
  terminate() { for (const w of this.workers) w.terminate(); this.workers = []; this.broken = true; }
  /** a new pipeline run: the operators of the last one are dropped (the memory is reused) */
  reset() { this.top = this.base; this.gen++; this.vecs.clear(); this.pick.clear(); }
  /** what the solves of this run use so far: 'threads', 'one' (one thread measured faster) or 'off' */
  mode() {
    if (this.broken) return 'off';
    let pool = 0, one = 0; for (const p of this.pick.values()) { if (p.use === 'pool') pool++; else if (p.use) one++; }
    return one && !pool ? 'one' : 'threads';
  }

  alloc(nbytes) {                                    // between solves only: the helpers are idle while the memory grows
    const p = (this.top + 15) & ~15, end = p + Math.max(8, nbytes), have = this.memory.buffer.byteLength;
    if (end > have) this.memory.grow(Math.ceil((end - have) / 65536));
    this.top = end; return p;
  }
  /** a typed array in the shared memory, copied once per run */
  put(a) {
    if (!a) return 0;
    const c = this.cache.get(a); if (c && c.gen === this.gen) return c.p;
    const p = this.alloc(a.byteLength);
    new Uint8Array(this.memory.buffer, p, a.byteLength).set(new Uint8Array(a.buffer, a.byteOffset, a.byteLength));
    this.cache.set(a, { gen: this.gen, p }); return p;
  }
  /** the operator (applyA.par) as pointers into the shared memory */
  bind(par) {
    const c = this.ops.get(par); if (c && c.gen === this.gen) return c.op;
    let op;
    if (par.kind === 'bsr6') op = { kind: 'bsr6', units: par.nq, bp: this.put(par.bp), bc: this.put(par.bc), bv: this.put(par.bv), em6: this.put(par.em6) };
    else {
      const R = groupRefl(par, this.refl), mbo = par.kind === 'mbo';
      op = { kind: par.kind, units: par.n, n: par.n, ip: this.put(par.ip), ix: this.put(par.ix), dr: this.put(par.dr), di: this.put(par.di), em: this.put(par.em),
        rp: R ? this.put(R.rp) : 0, rj: R ? this.put(R.rj) : 0, rc: R ? this.put(R.rc) : 0, rs: R ? this.put(R.rs) : 0, m: mbo ? this.put(par.m) : 0, tau: mbo ? par.tau : 0, mbo: mbo ? 1 : 0 };
    }
    this.ops.set(par, { gen: this.gen, op }); return op;
  }
  /** the vectors of a solve with N unknowns (pointers), kept for the run */
  vectors(N) {
    let V = this.vecs.get(N); if (V) return V;
    const a = () => this.alloc(8 * N);
    V = { N, x: a(), r: a(), z: a(), p: a(), Ap: a(), dg: a(), pr1: a(), pr2: a(), pr3: a(), acc: this.alloc(16) };
    this.vecs.set(N, V); return V;
  }
  f64(p, len) { return new Float64Array(this.memory.buffer, p, len); }

  spinUntil(cond) {
    let t0 = 0;
    for (let s = 1; ; s++) {
      if (cond()) return;
      if ((s & 4095) === 0) {
        if (Atomics.load(this.ctl, ERR)) throw new Error('a solver thread failed');
        if (!t0) t0 = performance.now(); else if (performance.now() - t0 > 15000) throw new Error('the solver threads do not answer');
      }
    }
  }
  /** one solve: hand the helpers the operator and the vectors, wait until all of them listen */
  begin(op, V) {
    const { ctl } = this;
    Atomics.store(ctl, READY, 0);
    for (const w of this.workers) w.postMessage({ type: 'session', op, V });
    this.spinUntil(() => Atomics.load(ctl, READY) === this.K);
    this.op = op; this.V = V;
  }
  phase(kernel, nchunks) {
    const { ctl } = this;
    Atomics.store(ctl, KERNEL, kernel); Atomics.store(ctl, NCHUNK, nchunks); Atomics.store(ctl, NEXT, 0); Atomics.store(ctl, DONE, 0);
    const g = Atomics.add(ctl, GEN, 1) + 1;
    Atomics.notify(ctl, GEN);
    return g;
  }
  waitChunk(c, g) { if (Atomics.load(this.cdone, c) !== g) this.spinUntil(() => Atomics.load(this.cdone, c) === g); }
  waitDone() { if (Atomics.load(this.ctl, DONE) !== this.K) this.spinUntil(() => Atomics.load(this.ctl, DONE) === this.K); }
  end() { this.phase(K_END, 0); this.op = null; this.V = null; }

  // the three phases in the coordinator alone (the measurement, and a solve after the pool stepped aside): the same rows
  // and entries, and the same sums in the same order, as the phases on the pool below
  spmvDotAlone(op, V) {
    const ex = this.ex;
    if (op.kind === 'bsr6') ex.rows_bsr6(op.bp, op.bc, op.bv, op.em6, V.p, V.Ap, V.pr1, 0, op.units);
    else ex.rows_rosy(op.n, op.ip, op.ix, op.dr, op.di, op.em, op.rp, op.rj, op.rc, op.rs, op.m, op.tau, op.mbo, V.p, V.Ap, V.pr1, 0, op.units);
    return ex.sum_into(0, V.pr1, 0, V.N);
  }
  updateDotsAlone(V, alpha) {
    const ex = this.ex, acc = this.f64(V.acc, 2); acc[0] = 0; acc[1] = 0;
    ex.cg_update(V.x, V.p, V.r, V.Ap, V.dg, V.z, V.pr2, V.pr3, alpha, 0, V.N);
    ex.sum2_into(V.acc, V.pr2, V.pr3, 0, V.N); this.rz = acc[0]; this.rr = acc[1];
  }
  pupdateAlone(V, beta) { this.ex.cg_pupdate(V.p, V.z, beta, 0, V.N); }

  /** Ap = A p on the helpers; returns dot(p, Ap), summed here in index order while the chunks come in */
  spmvDot() {
    const op = this.op, V = this.V, ex = this.ex, nch = Math.ceil(op.units / ROWS), g = this.phase(K_SPMV, nch);
    let s = 0;
    if (op.kind === 'bsr6') {                             // row pair q holds entries 2q, 2q + 1: chunks are contiguous
      const N = 2 * op.units;
      for (let c = 0; c < nch; c++) { this.waitChunk(c, g); const lo = 2 * ROWS * c; s = ex.sum_into(s, V.pr1, lo, Math.min(N, lo + 2 * ROWS)); }
    } else {                                              // row i holds entries i and n + i: the re half first, then the im half
      const n = op.n;
      for (let c = 0; c < nch; c++) { this.waitChunk(c, g); const lo = ROWS * c; s = ex.sum_into(s, V.pr1, lo, Math.min(n, lo + ROWS)); }
      s = ex.sum_into(s, V.pr1, n, 2 * n);
    }
    this.waitDone();
    return s;
  }
  /** the CG update on the helpers; this.rz = dot(r, z), this.rr = dot(r, r), summed here in index order */
  updateDots(alpha) {
    const V = this.V, ex = this.ex, N = V.N, nch = Math.ceil(N / ELEMS); this.scal[0] = alpha;
    const acc = this.f64(V.acc, 2); acc[0] = 0; acc[1] = 0;
    const g = this.phase(K_UPDATE, nch);
    for (let c = 0; c < nch; c++) { this.waitChunk(c, g); const lo = ELEMS * c; ex.sum2_into(V.acc, V.pr2, V.pr3, lo, Math.min(N, lo + ELEMS)); }
    this.waitDone(); this.rz = acc[0]; this.rr = acc[1];
  }
  /** p = z + β p; nothing to sum, so the coordinator takes chunks too */
  pupdate(beta) {
    const V = this.V, N = V.N, nch = Math.ceil(N / ELEMS), { ctl, cdone, ex } = this; this.scal[1] = beta;
    const g = this.phase(K_PUPDATE, nch);
    for (;;) { const c = Atomics.add(ctl, NEXT, 1); if (c >= nch) break; runChunk(ex, K_PUPDATE, null, V, this.scal, c); Atomics.store(cdone, c, g); }
    this.waitDone();
  }
}

/** sparse.js · pcg on the pool (projectMean is not supported: those solves stay single-threaded). null = not here:
 *  sparse.js then runs the one-thread WebAssembly (also when the pool measured itself slower in this run, see MARGIN). */
export function pcgPool(pool, applyA, diag, b, x0, opts = {}) {
  const N = b.length;
  if (!pool || pool.broken || !pool.ex || N < MIN_N || !applyA.par) return null;
  let pk = pool.pick.get(N);
  if (!pk) pool.pick.set(N, pk = { use: pool.force === 'pool' || pool.force === 'alone' ? pool.force : null, tp: 0, t1: 0, wt: 0, wn: 0, debt: 0, late: 0 });
  if (pk.use === 'one' && wasmOn()) return null;
  return runPool(pool, applyA, diag, b, x0, opts, pk);
}

function runPool(pool, applyA, diag, b, x0, { maxIter = 500, tol = 1e-8 } = {}, pk) {
  const N = b.length;
  try {
    const op = pool.bind(applyA.par), V = pool.vectors(N);
    const x = pool.f64(V.x, N), r = pool.f64(V.r, N), z = pool.f64(V.z, N), p = pool.f64(V.p, N), Ap = pool.f64(V.Ap, N);
    if (x0) x.set(x0); else x.fill(0);
    const ex = pool.ex;                                    // A x0 here, once per solve
    if (op.kind === 'bsr6') ex.rows_bsr6(op.bp, op.bc, op.bv, op.em6, V.x, V.Ap, 0, 0, op.units);
    else ex.rows_rosy(op.n, op.ip, op.ix, op.dr, op.di, op.em, op.rp, op.rj, op.rc, op.rs, op.m, op.tau, op.mbo, V.x, V.Ap, 0, 0, op.units);
    for (let i = 0; i < N; i++) r[i] = b[i] - Ap[i];
    const bnorm = Math.sqrt(dot(b, b)) || 1.0;
    for (let i = 0; i < N; i++) { z[i] = r[i] / diag[i]; p[i] = z[i]; }
    let rz = dot(r, z), it = 0, rn = Math.sqrt(dot(r, r));
    if (it < maxIter && rn > tol * bnorm) {
      pool.f64(V.dg, N).set(diag);
      let helpers = pk.use === null || pk.use === 'pool', measure = pk.use === null, tA = 0, tB = 0, t1 = Infinity;
      const off = () => { pool.end(); helpers = false; pk.use = 'one'; };      // the rest of this solve alone
      if (helpers) { const tb = performance.now(); pool.begin(op, V); if (pk.t1) pk.wt += performance.now() - tb; }
      try {
        while (it < maxIter && rn > tol * bnorm) {
          const alone = !helpers || (measure && it >= WARM + PROBE), watch = helpers && !measure && pk.t1 > 0;
          const ti = measure || watch ? performance.now() : 0;
          if (measure && it === WARM) tA = ti;
          const pAp = alone ? pool.spmvDotAlone(op, V) : pool.spmvDot();
          if (!(pAp > 0)) break;
          const alpha = rz / pAp;
          if (alone) pool.updateDotsAlone(V, alpha); else pool.updateDots(alpha);
          const rzNew = pool.rz, beta = rzNew / rz; rz = rzNew;
          if (alone) pool.pupdateAlone(V, beta); else pool.pupdate(beta);
          rn = Math.sqrt(pool.rr); it++;
          if (measure) {                                     // the pool against the coordinator alone, per iteration
            const tj = performance.now();
            if (it === WARM + PROBE) tB = tj; else if (it > WARM + PROBE) t1 = Math.min(t1, tj - ti);
            if (it === WARM + 2 * PROBE) {
              measure = false; pk.tp = (tB - tA) / PROBE; pk.t1 = t1;
              pk.use = pool.force === 'late' || (pool.force !== 'one' && pk.tp < MARGIN * t1) ? 'pool' : 'one';
              if (pool.force === 'late') pk.t1 = 1e-9;
              if (pk.use === 'one') off();
            }
          } else if (watch) {                                // the pool keeps proving itself (WINDOW)
            pk.wt += performance.now() - ti;
            if (++pk.wn === WINDOW) {
              const q = pk.wt / pk.wn / pk.t1; pk.wt = 0; pk.wn = 0; pk.debt = Math.max(0, pk.debt + q - 1);
              if (pk.debt > DEBT) { pk.late = q; off(); }
            }
          }
        }
      } finally { if (helpers) pool.end(); }
    }
    return { x: Float64Array.from(pool.f64(V.x, N)), iters: it, relres: rn / bnorm };
  } catch (e) {
    pool.broken = true;                                    // from now on single-threaded; this solve restarts there
    console.warn('parallel CG off:', e && e.message || e);
    return null;
  }
}

/** The pool for this worker and the hook in sparse.js · pcg; null (single thread) without SharedArrayBuffer or when the
 *  page is not cross-origin isolated. threads: helpers (default HELPERS, fewer on small machines); spawn: Worker
 *  factory. Await pool.ready before the first solve (global_worker.js does). */
export function startPool({ threads = null, spawn = null } = {}) {
  if (typeof SharedArrayBuffer === 'undefined' || typeof WebAssembly === 'undefined' || globalThis.crossOriginIsolated === false) return null;
  const cores = (globalThis.navigator && navigator.hardwareConcurrency) || 4;
  const K = threads ?? Math.min(HELPERS, cores - 2);
  if (!(K >= (threads ? 1 : 2))) return null;              // one helper alone is slower than no pool
  if (!spawn) {
    // the helpers start from a blob module that imports pcg_worker.js: in a cross-origin isolated page, a worker
    // started from a worker with an http URL fails to load in Chromium 152 (the Claude app's browser), a blob one does
    const url = URL.createObjectURL(new Blob([`import ${JSON.stringify(new URL('./pcg_worker.js', import.meta.url).href)};`], { type: 'text/javascript' }));
    spawn = () => new Worker(url, { type: 'module' });
  }
  let pool;
  try { pool = new Pool(K, spawn); } catch (e) { console.warn('no solver threads:', e && e.message || e); return null; }
  pool.ready.catch((e) => { pool.broken = true; console.warn('solver threads off:', e && e.message || e); });
  pool.pcg = (applyA, diag, b, x0, opts) => pcgPool(pool, applyA, diag, b, x0, opts);
  setParallelPCG(pool.pcg);
  return pool;
}

// ------------------------------------------------------------------ the helper side (pcg_worker.js)
/** listen(handler): installs the message handler; reply(msg): posts to the coordinator. 'init' brings the control
 *  block, the compiled kernels and the shared memory; a 'session' message runs the phase loop until the coordinator
 *  ends the solve, then returns to the event loop for the next one. */
export function helperMain(listen, reply) {
  let ctl = null, cdone = null, scal = null, ex = null;
  listen((d) => {
    if (d.type === 'init') {
      ctl = d.ctl; cdone = d.cdone; scal = d.scal;
      ex = new WebAssembly.Instance(d.module, { env: { memory: d.memory } }).exports;
      reply({ type: 'ready' }); return;
    }
    if (d.type !== 'session') return;
    try {
      const { op, V } = d;
      let seen = Atomics.load(ctl, GEN);
      Atomics.add(ctl, READY, 1);
      for (;;) {
        for (let s = 0; Atomics.load(ctl, GEN) === seen; s++) if (s > SPIN) { Atomics.wait(ctl, GEN, seen, 100); s = 0; }
        seen = Atomics.load(ctl, GEN);
        const kernel = Atomics.load(ctl, KERNEL);
        if (kernel === K_END) return;
        const nch = Atomics.load(ctl, NCHUNK);
        for (;;) {
          const c = Atomics.add(ctl, NEXT, 1); if (c >= nch) break;
          runChunk(ex, kernel, op, V, scal, c);
          Atomics.store(cdone, c, seen);
        }
        Atomics.add(ctl, DONE, 1);
      }
    } catch (e) { Atomics.store(ctl, ERR, 1); throw e; }
  });
}
