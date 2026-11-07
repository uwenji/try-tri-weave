// The CG solves in WebAssembly, one thread (pcg_kernels.c, embedded in pcg_wasm_bin.js): needs no special headers, so
// it runs on any host. sparse.js · pcg hands each solve whose operator carries a kernel description (applyA.par, set
// by rosy.js and stripes.js) to pcgWasm: the same numbers as the JS loop, bit for bit, about twice as fast.
// Memory: one arena. The operators are copied in once per pipeline run (resetWasm() starts a run); a solve's vectors
// sit above them for that solve only.
import { ST } from './pcg_wasm_bin.js';
import { setWasmPCG } from './sparse.js';

export const wasmBytes = (b64) => Uint8Array.from(atob(b64), (c) => c.charCodeAt(0));
let W = null;                                         // {ex, memory, base, top, gen, cache, refl}

/** compile and register (async: browsers refuse a synchronous compile of > 4 KB on the main thread). true when on */
export async function initWasm() {
  if (W) return true;
  if (typeof WebAssembly === 'undefined') return false;
  try {
    const memory = new WebAssembly.Memory({ initial: 256 });                 // 16 MB, grows with the carrier
    const { instance } = await WebAssembly.instantiate(wasmBytes(ST), { env: { memory } });
    const base = (instance.exports.__heap_base.value + 15) & ~15;
    W = { ex: instance.exports, memory, base, top: base, gen: 1, cache: new WeakMap(), refl: new WeakMap() };
    setWasmPCG(pcgWasm);
    return true;
  } catch (e) { console.warn('WebAssembly kernels off:', e && e.message || e); return false; }
}
export const wasmOn = () => !!W;
/** a new pipeline run: the operators of the last one are dropped (the memory is reused) */
export function resetWasm() { if (W) { W.top = W.base; W.gen++; } }

function alloc(nbytes) {
  const p = (W.top + 15) & ~15, end = p + Math.max(8, nbytes), have = W.memory.buffer.byteLength;
  if (end > have) W.memory.grow(Math.ceil((end - have) / 65536));
  W.top = end; return p;
}
/** a typed array in the arena, copied once per run */
function put(a) {
  if (!a) return 0;
  const c = W.cache.get(a); if (c && c.gen === W.gen) return c.p;
  const p = alloc(a.byteLength);
  new Uint8Array(W.memory.buffer, p, a.byteLength).set(new Uint8Array(a.buffer, a.byteOffset, a.byteLength));
  W.cache.set(a, { gen: W.gen, p }); return p;
}
/** the seam reflections of a non-orientable carrier (rosy.js: flat [i, j, wc, ws, …], rows ascending) grouped by row,
 *  in their original order: {rp, rj, rc, rs}, or null */
export function groupRefl(par, cache = null) {
  const R = par.refl, nR = R ? R.length / 4 : 0; if (!nR) return null;
  let g = cache && cache.get(par); if (g) return g;
  const n = par.n, rp = new Int32Array(n + 1), rj = new Int32Array(nR), rc = new Float64Array(nR), rs = new Float64Array(nR);
  for (let k = 0; k < nR; k++) rp[R[4 * k] + 1]++;
  for (let i = 0; i < n; i++) rp[i + 1] += rp[i];
  const at = rp.slice(0, n);
  for (let k = 0; k < nR; k++) { const t = at[R[4 * k]]++; rj[t] = R[4 * k + 1]; rc[t] = R[4 * k + 2]; rs[t] = R[4 * k + 3]; }
  g = { rp, rj, rc, rs }; if (cache) cache.set(par, g); return g;
}

/** sparse.js · pcg in WebAssembly (projectMean is not supported: those solves stay in JS). null = not here. */
export function pcgWasm(applyA, diag, b, x0, { maxIter = 500, tol = 1e-8 } = {}) {
  const par = applyA.par, N = b.length;
  if (!W || !par) return null;
  let run;
  if (par.kind === 'bsr6') {
    const bp = put(par.bp), bc = put(par.bc), bv = put(par.bv), em6 = put(par.em6);
    run = (d, bb, x, w, res) => W.ex.pcg_bsr6(par.nq, bp, bc, bv, em6, d, bb, x, w, maxIter, tol, res);
  } else {
    const R = groupRefl(par, W.refl), mbo = par.kind === 'mbo';
    const p = [put(par.ip), put(par.ix), put(par.dr), put(par.di), put(par.em), R ? put(R.rp) : 0, R ? put(R.rj) : 0, R ? put(R.rc) : 0, R ? put(R.rs) : 0, mbo ? put(par.m) : 0];
    run = (d, bb, x, w, res) => W.ex.pcg_rosy(par.n, ...p, mbo ? par.tau : 0, mbo ? 1 : 0, d, bb, x, w, maxIter, tol, res);
  }
  const keep = W.top, pd = alloc(8 * N), pb = alloc(8 * N), px = alloc(8 * N), pw = alloc(32 * N), pr = alloc(16);
  const buf = W.memory.buffer;
  new Float64Array(buf, pd, N).set(diag); new Float64Array(buf, pb, N).set(b);
  if (x0) new Float64Array(buf, px, N).set(x0); else new Float64Array(buf, px, N).fill(0);
  run(pd, pb, px, pw, pr);
  const B = W.memory.buffer, res = new Float64Array(B, pr, 2), x = Float64Array.from(new Float64Array(B, px, N));
  W.top = keep;
  return { x, iters: res[0], relres: res[1] };
}
