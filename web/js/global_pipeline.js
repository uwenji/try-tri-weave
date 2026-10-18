// The global (no-seed) pipeline, shared by global_worker.js (the app) and the Node tests:
//   smoothest 6-RoSy on the base carrier (model metric g)
//   → cones and the conformal size field ĝ = e^{2φ} g (cones.js; skipped in mode 'uniform'), with the user's edits of the
//     cones applied (off / valence / radius) and their alignment constraints (angle, weight) imposed on the ĝ field
//   → smoothest 6-RoSy in ĝ (its singularities sit at the cones) → lifted to the refined carrier
//   → stripe phases in ĝ (constant q0 there) → lifted → u0 and branch points → per-face phases for the stripes view.
// Everything heavy runs on the base carrier; the refined carrier only receives lifts. `cache` (one object per carrier)
// keeps the model-metric field and, per (mode, a), the automatic cones; per effective cone list the ĝ field and its lift.
import { Carrier } from './carrier.js';
import { smoothestRoSy, liftField, nb0Of } from './rosy.js';
import { stripes, liftStripes, stripeFaces, farthestVertex } from './stripes.js';
import { placeCones, conformalLengths, interpByParents, sizeFieldFor, applyConeEdits, discAround } from './cones.js';

/** a Carrier from the arrays app.js sends to the worker (same vertex order, already centred) */
export function carrierFromPack(d) {
  return new Carrier({ name: d.name, meta: { h: d.h, hBase: d.hBase || d.h }, V: d.V, F: d.F, Fi: d.Fi, Li: d.Li, K: d.K, B: d.B, sideChecked: true });
}

/**
 * globalPipeline(cb, cf, opts) → the seed and its diagnostics.
 *   cb: base Carrier · cf: {n, Fi, Li, boundary, V, N, F} of the refined carrier or null · opts.parents: midpoint parents
 *   per level (Carrier.refined) · opts.mode: 'uniform' | 'greedy' | 'field' · opts.edits: user edits of the cones
 *   [{vertex, on, valence, radius, angle, weight}] · opts.post(stage): progress.
 */
export function globalPipeline(cb, cf, { q0, a, mode = 'uniform', parents = null, edits = null, cache = {}, post = () => {}, log = null, coneOpts = { relax: 5 } } = {}) {
  const t = {}, now = () => performance.now();
  let t0 = now();
  if (!cache.rosyG) { post('solving the rotation field on the base carrier…'); cache.rosyG = smoothestRoSy(cb); t.rosy = now() - t0; }
  const rosyG = cache.rosyG;
  let placed = null, phiF = null, LenB = null, LenF = null, rosy = rosyG, lift = null, coneList = [], constraints = [];
  const key = `${mode}|${a.toPrecision(5)}`;
  if (mode !== 'uniform') {
    if (!cache.placed || cache.placed.key !== key) {
      post(`placing the cones of the size field ĝ (${mode})…`); t0 = now();
      const v = placeCones(cb, a, { mode, fieldSingular: rosyG.singular, log, ...coneOpts }); t.cones = now() - t0;
      cache.placed = { key, v };
    }
    const auto = cache.placed.v;
    // the user's edits: dropped cones, other valences and radii change φ; angles with a weight constrain the ĝ field
    const ed = applyConeEdits(auto.cones, edits, a); coneList = ed.all;
    const edKey = JSON.stringify(ed.active), same = edKey === JSON.stringify(auto.cones.map(q => ({ vertex: q.vertex, valence: q.valence })));
    if (!cache.phi || cache.phi.key !== edKey) { t0 = now(); cache.phi = { key: edKey, v: same ? auto : { ...auto, ...sizeFieldFor(cb, ed.active, a, { smear: 0.5 * a }) } }; t.cones = (t.cones || 0) + now() - t0; cache.rosyHat = null; cache.liftHat = null; }
    placed = cache.phi.v;
    LenB = conformalLengths(cb.Fi, cb.Li, placed.phi);
    if (cf) { phiF = interpByParents(placed.phi, parents || []); LenF = conformalLengths(cf.Fi, cf.Li, phiF); } else phiF = placed.phi;
    if (!cache.rosyFree || cache.rosyFree.key !== edKey) { post('solving the rotation field in ĝ…'); t0 = now(); cache.rosyFree = { key: edKey, v: smoothestRoSy(cb, { Len: LenB }) }; t.rosyHat = now() - t0; cache.rosyHat = null; cache.liftHat = null; }
    const free = cache.rosyFree.v;
    // alignment constraints: a ring of vertices 0.4–1.2 a around an edited cone is asked for the free direction turned by its angle
    for (const q of coneList) {
      if (!(q.weight > 0)) continue;
      const d = q.angle * Math.PI / 180;
      for (const [i, dist] of discAround(cb, q.vertex, 1.2 * a)) if (dist >= 0.4 * a && !cb.boundary[i]) constraints.push({ vertex: i, theta: free.theta[i] + d, weight: q.weight });
    }
    const cKey = edKey + '|' + JSON.stringify(coneList.filter(q => q.weight > 0).map(q => [q.vertex, q.angle, q.weight]));
    if (!cache.rosyHat || cache.rosyHat.key !== cKey) {
      if (constraints.length) { post(`aligning the ĝ field to ${coneList.filter(q => q.weight > 0).length} edited cones…`); t0 = now(); cache.rosyHat = { key: cKey, v: smoothestRoSy(cb, { Len: LenB, constraints, x0: free.z, log }) }; t.align = now() - t0; }
      else cache.rosyHat = { key: cKey, v: free };
      cache.liftHat = null;
    }
    rosy = cache.rosyHat.v;
    if (cf) {
      if (!cache.liftHat || cache.liftHat.n !== cf.n) { post('lifting the ĝ rotation field to the refined carrier…'); t0 = now(); cache.liftHat = { n: cf.n, v: liftField(rosy.theta, rosy.frames, cb.n, cf, rosy.N, { Len: LenF }) }; t.lift = now() - t0; }
      lift = cache.liftHat.v;
    }
  } else if (cf) {
    if (!cache.liftG || cache.liftG.n !== cf.n) { post('lifting the rotation field to the refined carrier…'); t0 = now(); cache.liftG = { n: cf.n, v: liftField(rosyG.theta, rosyG.frames, cb.n, cf, rosyG.N) }; t.lift = now() - t0; }
    lift = cache.liftG.v;
  }
  post(`integrating the stripe phases on the base carrier${LenB ? ' in ĝ' : ''}…`); t0 = now();
  const st = stripes(cb, rosy, q0, { Len: LenB }); t.stripes = now() - t0;
  let out, ls = null;
  if (cf) {
    post('lifting the stripe phases…'); t0 = now();
    ls = liftStripes(st.psi, cb.n, cf, lift.theta, lift.frames, q0, rosy.N, { reversed: lift.reversed, Len: LenF }); t.liftStripes = now() - t0;
    out = { theta: lift.theta, nb0: nb0Of(lift.frames), reversed: lift.reversed, psi: ls.psi, u0: ls.u0, branch: ls.branch, branchOnBase: false };
  } else {
    out = { theta: rosy.theta, nb0: nb0Of(rosy.frames), reversed: rosy.reversed, psi: st.psi, u0: st.u0, branch: st.branch, branchOnBase: true };
  }
  // phases per render face of the shown carrier, colour slots walked from the vertex farthest from the field singularities
  t0 = now();
  const cR = cf || cb, nbR = (cf ? lift.frames : rosy.frames).map(f => f.nb), src = [];
  for (const sg of rosy.singular) for (let k = 0; k < 3; k++) src.push(cb.Fi[3 * sg.face + k]);        // base ids are kept by the refinement
  out.stripeFaces = cR.F ? stripeFaces(cR.F, cR.n, nbR, cf ? ls.edge : st.edge, out.psi, farthestVertex(nbR, cR.n, src), cR.V) : null; t.faces = now() - t0;
  const coneInfo = placed ? { edges: placed.edges, spread: placed.stats.spread, residualDeg: placed.residualDeg, net: placed.net, need: placed.need, closed: placed.closed, history: placed.history, constrained: constraints.length } : null;
  return { nBase: cb.n, mode, singular: rosy.singular, sumIndex: rosy.sumIndex, singularG: rosyG.singular, sumIndexG: rosyG.sumIndex,
    worstMatch: st.worstMatch, branchBase: st.branch.length, t, ...out, phi: phiF, cones: coneList, coneInfo };
}
