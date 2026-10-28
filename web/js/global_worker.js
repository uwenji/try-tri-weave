// Worker for the no-seed pipeline (global_pipeline.js): smoothest 6-RoSy, cones and size field ĝ, stripe phases, u0.
// Message in: {id, key, base: {name, h, hBase, n, Fi, Li, K, B, V, N, F}, fine: {n, Fi, Li, B, V, N, F, parents} | null, q0, a, mode,
//              edits: [{vertex, on, valence, radius, angle, weight}] (the user's changes to the cones, by base vertex),
//              mbo: {core, iters} | null (unit-length field), flat: {W, deg, ramp} | null (ĝ flat on the developable parts)}.
// Messages out: {id, stage} progress, then {id, done: true, …globalPipeline result} or {id, error}.
// The base carrier (rebuilt from the arrays) and the pipeline cache are kept per `key` (carrier name + level + n), so a
// density change reruns only the cones (greedy / field) and the stripes, and the model-metric field is solved once.
import { carrierFromPack, globalPipeline } from './global_pipeline.js';

let ent = { key: null, cb: null, cache: {} };
const asFine = (d) => ({ n: d.n, Fi: d.Fi, Li: d.Li, boundary: d.B, V: d.V, N: d.N, F: d.F });

self.onmessage = (e) => {
  const { id, key, base, fine, q0, a, mode, edits, mbo, flat } = e.data;
  try {
    if (ent.key !== key || !ent.cb) ent = { key, cb: carrierFromPack(base), cache: {} };
    const out = globalPipeline(ent.cb, fine ? asFine(fine) : null, { q0, a, mode: mode || 'uniform', parents: fine ? fine.parents : null, edits: edits || null, mbo: mbo || null, flat: flat || null, cache: ent.cache, post: (stage) => self.postMessage({ id, stage }) });
    self.postMessage({ id, done: true, ...out });
  } catch (err) { self.postMessage({ id, error: String(err && err.message || err) }); }
};
