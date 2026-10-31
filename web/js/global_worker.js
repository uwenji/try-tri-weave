// Worker for the no-seed pipeline (global_pipeline.js): smoothest 6-RoSy, cones and size field ĝ, stripe phases, u0.
// Message in: {id, key, base: {name, h, hBase, n, Fi, Li, K, B, V, N, F}, fine: {n, Fi, Li, B, V, N, F, parents} | null, q0, a, mode,
//              edits: [{vertex, on, valence, radius, angle, weight}] (the user's changes to the cones, by base vertex),
//              mbo: {core, iters} | null (unit-length field), flat: {W, deg, ramp} | null (ĝ flat on the developable parts)}.
// Messages out: {id, stage} progress, then {id, done: true, …globalPipeline result} or {id, error}.
// The base carrier (rebuilt from the arrays) and the pipeline cache are kept per `key` (carrier name + level + n), so a
// density change reruns only the cones (greedy / field) and the stripes, and the model-metric field is solved once.
// The CG solves of the rotation field and the stripes run on the fastest of (all with the same numbers, bit for bit):
//   the thread pool of parallel.js (WebAssembly kernels on shared memory; needs a cross-origin isolated page), else
//   the WebAssembly kernels in this thread (pcg_wasm.js; any host), else JS. The stage messages end with what ran
//   ("· 5 threads, wasm" / "· 1 thread, wasm" / "· 1 thread, js"; "(threads slower now)": the pool measured itself slower
//   than one thread in this run — a browser at low priority — and stepped aside, see parallel.js · MARGIN).
import { carrierFromPack, globalPipeline } from './global_pipeline.js';
import { startPool } from './parallel.js';
import { initWasm, resetWasm, wasmOn } from './pcg_wasm.js';

const pool = startPool();
const threads = () => (pool && !pool.broken ? pool.K + 1 : 1);
const wasmReady = initWasm();
// the helpers load through this worker's event loop, so they are awaited here, before a solve keeps this thread busy;
// a pool that does not start within 10 s is dropped (one thread)
const poolReady = !pool ? Promise.resolve() : Promise.race([pool.ready, new Promise((_, no) => setTimeout(() => no(new Error('timeout')), 10000))])
  .catch((err) => { pool.broken = true; console.warn('solver threads off:', err && err.message || err); });

let ent = { key: null, cb: null, cache: {} };
const asFine = (d) => ({ n: d.n, Fi: d.Fi, Li: d.Li, boundary: d.B, V: d.V, N: d.N, F: d.F });

self.onmessage = async (e) => {
  await wasmReady; await poolReady;
  const { id, key, base, fine, q0, a, mode, edits, mbo, flat } = e.data;
  try {
    resetWasm(); if (pool && !pool.broken) pool.reset();
    if (ent.key !== key || !ent.cb) ent = { key, cb: carrierFromPack(base), cache: {} };
    const used = () => (pool && pool.mode() === 'threads' ? threads() : 1);    // the pool steps aside when one thread measures faster
    const how = () => { const n = used(); return `${n} thread${n > 1 ? 's' : ''}, ${n > 1 || wasmOn() ? 'wasm' : 'js'}${n === 1 && pool && pool.mode() === 'one' ? ' (threads slower now)' : ''}`; };
    const out = globalPipeline(ent.cb, fine ? asFine(fine) : null, { q0, a, mode: mode || 'uniform', parents: fine ? fine.parents : null, edits: edits || null, mbo: mbo || null, flat: flat || null, cache: ent.cache, post: (stage) => self.postMessage({ id, stage: `${stage} · ${how()}` }) });
    out.t.threads = used(); out.t.wasm = wasmOn();
    // per solve size: what the pool chose and its time per iteration over the coordinator's alone (parallel.js · MARGIN)
    // ("threads 0.45 → one 2.10": the pool stayed, then averaged 2.1× one iteration alone and stepped aside)
    out.t.pick = pool ? [...pool.pick.values()].filter((p) => p.t1).map((p) => { const r = (p.tp / p.t1).toFixed(2);
      return p.late ? `threads ${r} → one ${p.late.toFixed(2)}` : `${p.use === 'pool' ? 'threads' : 'one'} ${r}`; }).join(', ') : '';
    self.postMessage({ id, done: true, ...out });
  } catch (err) { self.postMessage({ id, error: String(err && err.message || err) }); }
};
