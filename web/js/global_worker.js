// Worker for the no-seed pipeline: smoothest 6-RoSy on the base carrier → lift to the refined carrier → stripe phases on
// the base → lift → u0 and branch points. Message in: {id, key, base: {n, Fi, Li, B, V, N}, fine: same | null, q0}.
// Messages out: {id, stage} progress, then {id, done: true, nBase, theta, nb0, reversed, singular, psi, u0, branch, t}.
// The base rotation field and its lift are cached by `key` (carrier name + level + n), so a density change only reruns
// the stripes.
import { smoothestRoSy, liftField, nb0Of } from './rosy.js';
import { stripes, liftStripes } from './stripes.js';

const cache = { key: null, rosy: null, lift: null };
const asCarrier = (d) => ({ n: d.n, Fi: d.Fi, Li: d.Li, boundary: d.B, V: d.V, N: d.N });

self.onmessage = (e) => {
  const { id, key, base, fine, q0 } = e.data, t = {};
  try {
    const cb = asCarrier(base), cf = fine ? asCarrier(fine) : null;
    let rosy, lift;
    if (cache.key === key && cache.rosy) { rosy = cache.rosy; lift = cache.lift; }
    else {
      self.postMessage({ id, stage: 'solving the rotation field on the base carrier…' });
      let t0 = performance.now(); rosy = smoothestRoSy(cb); t.rosy = performance.now() - t0;
      lift = null;
      if (cf) { self.postMessage({ id, stage: 'lifting the rotation field to the refined carrier…' }); t0 = performance.now(); lift = liftField(rosy.theta, rosy.frames, cb.n, cf, rosy.N); t.lift = performance.now() - t0; }
      cache.key = key; cache.rosy = rosy; cache.lift = lift;
    }
    self.postMessage({ id, stage: 'integrating the stripe phases on the base carrier…' });
    let t0 = performance.now(); const st = stripes(cb, rosy, q0); t.stripes = performance.now() - t0;
    let out;
    if (cf) {
      self.postMessage({ id, stage: 'lifting the stripe phases…' }); t0 = performance.now();
      const ls = liftStripes(st.psi, cb.n, cf, lift.theta, lift.frames, q0, rosy.N, { reversed: lift.reversed }); t.liftStripes = performance.now() - t0;
      out = { theta: lift.theta, nb0: nb0Of(lift.frames), reversed: lift.reversed, psi: ls.psi, u0: ls.u0, branch: ls.branch, branchOnBase: false };
    } else {
      out = { theta: rosy.theta, nb0: nb0Of(rosy.frames), reversed: rosy.reversed, psi: st.psi, u0: st.u0, branch: st.branch, branchOnBase: true };
    }
    self.postMessage({ id, done: true, nBase: cb.n, singular: rosy.singular, sumIndex: rosy.sumIndex, worstMatch: st.worstMatch, branchBase: st.branch.length, t, ...out });
  } catch (err) { self.postMessage({ id, error: String(err && err.message || err) }); }
};
