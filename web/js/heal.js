// Healing bound 5–7 pairs in the SH field.
// A 5-7-5 whose two 5s sit in the 7's link with two regular sites between them is a pentagon of sites with its centre
// missing: the 7 holds both diagonals of the pentagon, its two neighbours on the pentagon have lost an edge. Putting a
// site at the pentagon's centre gives one 5 there and six corners of valence 6 (net charge +1 either way).
// The mirror, a 7-5-7 whose two 7s are separated by a single site B in the 5's link, is a heptagon holding two sites
// where one belongs: merging the 5 with B gives one 7 at the centre.
// Both are done inside the solver, not on the mesh: a bump (trough) is written into u at the centre and held for a few
// steps while the neighbourhood relaxes, then the field runs free and the lattice is read again; a motif that does not
// resolve is reverted. The remaining arrangements (a quad with the wrong diagonal, 5-5-7 / 7-7-5 triangles) are counted
// but left alone.
import { discAround } from './cones.js';
import { extractLattice } from './lattice.js';

/** cyclic link (ccw neighbour order) of every peak from the oriented dual triangles; null where the link is not one closed cycle */
export function peakLinks(lat) {
  const n = lat.peaks.length, T = lat.tris, nxt = Array.from({ length: n }, () => new Map());
  for (let t = 0; t < T.length / 3; t++) { const a = T[3 * t], b = T[3 * t + 1], c = T[3 * t + 2]; nxt[a].set(b, c); nxt[b].set(c, a); nxt[c].set(a, b); }
  const links = new Array(n).fill(null);
  for (let i = 0; i < n; i++) {
    const m = nxt[i]; if (!m.size) continue;
    const start = m.keys().next().value, cyc = [start]; let cur = m.get(start);
    while (cur !== undefined && cur !== start && cyc.length <= m.size) { cyc.push(cur); cur = m.get(cur); }
    if (cur === start && cyc.length === m.size) links[i] = cyc;
  }
  return links;
}

/**
 * findMotifs(lat) → [{kind, peaks, centreOf, ...}] — the bound pairs the lattice holds:
 *   insert: 5-7-5, the 5s separated by two regular sites in the 7's link → nucleate a site at the pentagon's centre
 *   merge:  7-5-7, the 7s separated by one regular site B in the 5's link → merge the 5 with B
 *   flip:   5-7-5 with one site between the 5s, or a 7-7-5 triangle (an edge flip would do; not healed here)
 *   tri:    5-5-7 triangle (no single local operation; not healed here)
 * `peaks` lists the motif's non-6 peaks; `corners` the polygon whose centre takes the new site (insert) or the pair to merge;
 * `pv`, `cv`: the same as carrier vertices (peak indices change between read-outs).
 */
export function findMotifs(lat, { links = peakLinks(lat) } = {}) {
  const val = lat.valence, inter = lat.interior, n = lat.peaks.length, out = [], used = new Uint8Array(n);
  const reg = (i) => inter[i] && val[i] === 6, sing = (i, v) => inter[i] && val[i] === v;
  const take = (m) => { if (m.peaks.some(i => used[i])) return; for (const i of m.peaks) used[i] = 1; out.push({ ...m, pv: m.peaks.map(i => lat.peaks[i]), cv: m.corners.map(i => lat.peaks[i]), av: m.alt ? m.alt.map(i => lat.peaks[i]) : null }); };
  for (let v = 0; v < n; v++) {
    const L = links[v]; if (!inter[v] || !L) continue;
    const k = L.length;
    if (val[v] === 7 && k === 7) {
      for (let i = 0; i < 7; i++) {
        const p = L[i]; if (!sing(p, 5)) continue;
        for (let s = 0; s <= 2; s++) {                                   // s regular sites between p and q going ccw
          const q = L[(i + s + 1) % 7]; if (!sing(q, 5)) continue;
          const between = []; for (let j = 1; j <= s; j++) between.push(L[(i + j) % 7]);
          if (!between.every(reg)) continue;
          if (s === 2) {
            const other = []; for (let j = 4; j <= 6; j++) other.push(L[(i + j) % 7]);         // the three link sites on the far side
            take({ kind: 'insert', peaks: [v, p, q], corners: [v, p, between[0], between[1], q], alt: [v, q, ...other, p], v7: v });
          }
          else if (s === 1) take({ kind: 'flip', peaks: [v, p, q], corners: [v, p, between[0], q], edge: [v, between[0]], opp: [p, q] });   // v7–x → p–q: x becomes the 5
          else take({ kind: 'tri', peaks: [v, p, q], corners: [v, p, q] });
        }
      }
    } else if (val[v] === 5 && k === 5) {
      for (let i = 0; i < 5; i++) {
        const p = L[i]; if (!sing(p, 7)) continue;
        for (let s = 0; s <= 1; s++) {
          const q = L[(i + s + 1) % 5]; if (!sing(q, 7)) continue;
          if (s === 1) { const B = L[(i + 1) % 5]; if (reg(B)) take({ kind: 'merge', peaks: [v, p, q], corners: [v, B], v5: v, B }); }
          else {
            // 7-7-5 triangle (v, p, q): flip p–q to v–w, w the far corner of the other triangle on p–q; w becomes the 7
            const Lp = links[p]; if (!Lp) continue;
            const j = Lp.indexOf(q), cand = [Lp[(j + 1) % Lp.length], Lp[(j - 1 + Lp.length) % Lp.length]].filter(w => w !== v);
            if (cand.length === 1 && reg(cand[0])) take({ kind: 'flip', peaks: [v, p, q], corners: [v, p, cand[0], q], edge: [p, q], opp: [v, cand[0]] });
          }
        }
      }
    }
  }
  return out;
}

/** carrier vertex nearest (in 3-space) to the centroid of the corner peaks `cs` (default: the motif's corners), among the vertices of their Voronoi cells */
export function motifCentre(lat, c, motif, cs = motif.corners) {
  const cx = [0, 0, 0];
  for (const p of cs) for (let d = 0; d < 3; d++) cx[d] += lat.P[3 * p + d] / cs.length;
  const inCells = new Set(cs); let best = -1, bd = Infinity;
  for (let i = 0; i < c.n; i++) {
    if (!inCells.has(lat.label[i])) continue;
    const dx = c.V[3 * i] - cx[0], dy = c.V[3 * i + 1] - cx[1], dz = c.V[3 * i + 2] - cx[2], d2 = dx * dx + dy * dy + dz * dz;
    if (d2 < bd) { bd = d2; best = i; }
  }
  return best;
}

/** where an insert motif can take its new site: the pentagon's centre (the new site is the 5) and the far side's hexagon
 *  centre (the new site is a 6 and the 7 becomes the 5), each with its distance to the nearest corner; roomiest first */
export function insertSites(lat, c, motif, a) {
  const out = [];
  for (const cs of [motif.corners, motif.alt].filter(Boolean)) {
    const v = motifCentre(lat, c, motif, cs); if (v < 0) continue;
    const disc = discAround(c, v, 2 * a); let dmin = Infinity;
    for (const p of cs) { const d = disc.get(lat.peaks[p]); if (d !== undefined && d < dmin) dmin = d; }
    out.push({ vertex: v, corners: cs, room: dmin / a });
  }
  return out.sort((x, y) => y.room - x.room);
}

/** median height of the peaks: the amplitude a nucleated site should have */
export function peakAmplitude(lat, u) {
  const h = Array.from(lat.peaks, p => u[p]).sort((x, y) => x - y);
  return h.length ? h[h.length >> 1] : 0;
}

/** a smooth bump profile amp·(1 − (d/ρ)²)² over the geodesic disc of radius ρ around v: Map vertex → value */
function bump(c, v, rho, amp) {
  const out = new Map();
  for (const [i, d] of discAround(c, v, rho)) { const t = 1 - (d / rho) * (d / rho); out.set(i, amp * t * t); }
  return out;
}

/**
 * kicks for a motif: [{vertex, value}] lists to write with max (bumps) and min (troughs).
 *   insert: a peak at the pentagon centre; merge: troughs on the two sites, a peak between them.
 */
export function motifKick(lat, c, motif, a, amp, { site = null, push = 0 } = {}) {
  if (motif.kind === 'insert') {
    const v = site !== null ? site : motifCentre(lat, c, motif);
    const kick = { centre: v, bumps: [bump(c, v, 0.42 * a, amp)], troughs: [] };
    if (push > 0) {
      // the 7 moves away from the new site: trough on its old place, peak at the point `push`·a farther out along c→7
      const v7 = lat.peaks[motif.v7], P = (i) => [c.V[3 * i], c.V[3 * i + 1], c.V[3 * i + 2]], C = P(v), S = P(v7);
      const dir = [S[0] - C[0], S[1] - C[1], S[2] - C[2]], len = Math.hypot(...dir) || 1;
      const T = dir.map((x, k) => S[k] + push * a * x / len);
      let best = -1, bd = Infinity;
      for (const [i] of discAround(c, v7, 1.2 * a)) { const Q = P(i), d2 = (Q[0] - T[0]) ** 2 + (Q[1] - T[1]) ** 2 + (Q[2] - T[2]) ** 2; if (d2 < bd) { bd = d2; best = i; } }
      if (best >= 0) { kick.troughs.push(bump(c, v7, 0.3 * a, -0.6 * amp)); kick.bumps.push(bump(c, best, 0.42 * a, amp)); kick.moved = best; }
    }
    return kick;
  }
  if (motif.kind === 'flip') {
    const P = (i) => [c.V[3 * i], c.V[3 * i + 1], c.V[3 * i + 2]], move = (from, towards, dist) => {
      const A = P(from), B = P(towards), d = [B[0] - A[0], B[1] - A[1], B[2] - A[2]], len = Math.hypot(...d) || 1, T = d.map((x, k) => A[k] + dist * x / len);
      let best = -1, bd = Infinity; for (const [i] of discAround(c, from, Math.abs(dist) + 0.3 * a)) { const Q = P(i), d2 = (Q[0] - T[0]) ** 2 + (Q[1] - T[1]) ** 2 + (Q[2] - T[2]) ** 2; if (d2 < bd) { bd = d2; best = i; } }
      return best;
    };
    const [i, j] = motif.edge.map(p => lat.peaks[p]), [k, l] = motif.opp.map(p => lat.peaks[p]), dlt = push > 0 ? push : 0.35;
    const ni = move(i, j, -dlt * a), nj = move(j, i, -dlt * a), nk = move(k, l, dlt * a), nl = move(l, k, dlt * a);
    if ([ni, nj, nk, nl].some(x => x < 0)) return null;
    const centre = motifCentre(lat, c, motif);
    return { centre, troughs: [i, j, k, l].map(v => bump(c, v, 0.3 * a, -0.6 * amp)), bumps: [ni, nj, nk, nl].map(v => bump(c, v, 0.42 * a, amp)) };
  }
  if (motif.kind === 'merge') {
    const v = motifCentre(lat, c, motif), [A, B] = motif.corners.map(p => lat.peaks[p]);
    return { centre: v, bumps: [bump(c, v, 0.42 * a, amp)], troughs: [bump(c, A, 0.3 * a, -0.6 * amp), bump(c, B, 0.3 * a, -0.6 * amp)] };
  }
  return null;
}

export function applyKick(u, kick) {
  for (const t of kick.troughs) for (const [i, val] of t) if (u[i] > val) u[i] = val;
  for (const b of kick.bumps) for (const [i, val] of b) if (u[i] < val) u[i] = val;
}

/** interior non-6 peaks within geodesic radius R of carrier vertex v, and their charge */
export function localDefects(lat, c, v, R) {
  const disc = discAround(c, v, R), pos = new Int32Array(c.n).fill(-1); lat.peaks.forEach((p, i) => { pos[p] = i; });
  let count = 0, charge = 0;
  for (const i of disc.keys()) { const p = pos[i]; if (p >= 0 && lat.interior[p] && lat.valence[p] !== 6) { count++; charge += 6 - lat.valence[p]; } }
  return { count, charge };
}

const nonSix = (lat) => { let k = 0; for (let i = 0; i < lat.valence.length; i++) if (lat.interior[i] && lat.valence[i] !== 6) k++; return k; };

/**
 * healPass(solver, carrier, a, opts) → report. Reads the lattice, finds the motifs, and heals the insert / merge ones one
 * after the other: kick → hold (the kick re-applied every `every` steps for `hold` steps) → `release` free steps → read
 * again; kept when the defects within 2.5 a of the centre went down and the total did not go up, else the field is put
 * back. opts.log(msg), opts.onStep(solver) (called after each stepping batch, for the UI), opts.max (motifs per pass).
 */
export async function healPass(solver, c, a, { hold = 60, every = 5, release = 400, max = 50, log = () => {}, onStep = null, kinds = ['insert', 'merge', 'flip'], readOpts = {}, bothSides = true, push = 0.4, pin = false, pins = null } = {}) {
  let lat = extractLattice(c, solver.u, a, readOpts);
  if (!lat) return { motifs: {}, healed: 0, tried: 0, results: [] };
  const counts = {}; for (const m of findMotifs(lat)) counts[m.kind] = (counts[m.kind] || 0) + 1;
  log(`motifs: ${JSON.stringify(counts)} · interior non-6 ${nonSix(lat)}`);
  const step = async (k) => { await solver.step(k); if (onStep) await onStep(solver); };
  const tried = [], results = []; let healed = 0, pinned = 0;
  const near = (v) => tried.some(t => { const dx = c.V[3 * t] - c.V[3 * v], dy = c.V[3 * t + 1] - c.V[3 * v + 1], dz = c.V[3 * t + 2] - c.V[3 * v + 2]; return dx * dx + dy * dy + dz * dz < a * a; });
  for (let round = 0; round < max; round++) {
    lat = extractLattice(c, solver.u, a, readOpts); if (!lat) break;
    // the next motif not yet attempted (motifs are re-read after every attempt: peaks move)
    let m = null, sites = null;
    for (const q of findMotifs(lat)) {
      if (!kinds.includes(q.kind)) continue;
      const ss = q.kind === 'insert' ? insertSites(lat, c, q, a).slice(0, bothSides ? 2 : 1) : [{ vertex: motifCentre(lat, c, q), corners: q.corners, room: 0 }];   // flip / merge: one attempt
      if (!ss.length || ss[0].vertex < 0 || near(ss[0].vertex)) continue;
      m = q; sites = ss; break;
    }
    if (!m) break;
    const amp = peakAmplitude(lat, solver.u);
    const u0 = Float64Array.from(solver.u), steps0 = solver.steps, snaps0 = solver.snapshots.length, totalBefore = nonSix(lat);
    let ok = false, pinKick = null;
    for (const site of sites) {
      const kick = motifKick(lat, c, m, a, amp, { site: site.vertex, push }); if (!kick) break;
      const before = localDefects(lat, c, kick.centre, 2.5 * a);
      const uk = Float64Array.from(solver.u); applyKick(uk, kick); writeU(solver, uk);
      for (let s = 0; s < hold; s += every) { await step(every); const uu = Float64Array.from(solver.u); applyKick(uu, kick); writeU(solver, uu); }
      // the held state: what a pin would keep — only an insert / merge whose held lattice is already cleaner qualifies
      const latH = extractLattice(c, solver.u, a, readOpts), held = latH ? localDefects(latH, c, kick.centre, 2.5 * a) : before;
      if (!pinKick && m.kind !== 'flip' && held.count < before.count && held.charge === before.charge) pinKick = { kick, held: held.count, before: before.count, u: Float64Array.from(solver.u), steps: solver.steps, snaps: solver.snapshots.length };
      await step(release);
      const lat2 = extractLattice(c, solver.u, a, readOpts), after = lat2 ? localDefects(lat2, c, kick.centre, 2.5 * a) : before, totalAfter = lat2 ? nonSix(lat2) : totalBefore;
      ok = after.count < before.count && totalAfter <= totalBefore && after.charge === before.charge;
      log(`${m.kind} at ${kick.centre} (${site.corners.length}-gon, room ${site.room.toFixed(2)} a): local non-6 ${before.count} → ${held.count} held → ${after.count} free, total ${totalBefore} → ${totalAfter} · ${ok ? 'kept' : 'reverted'}`);
      results.push({ kind: m.kind, centre: kick.centre, before: before.count, held: held.count, after: after.count, totalBefore, totalAfter, ok });
      if (ok) break;
      rewind(solver, u0, steps0, snaps0);
    }
    tried.push(sites[0].vertex);
    if (ok) healed++;
    else if (pin && pinKick && pins) {
      // held instead: back to the held state, and the kick stays applied after every batch of steps (app.js / the test driver re-apply `pins`)
      writeU(solver, pinKick.u); solver.steps = pinKick.steps; solver.snapshots.length = Math.min(solver.snapshots.length, pinKick.snaps);
      pins.push({ kind: m.kind, centre: pinKick.kick.centre, kick: pinKick.kick }); pinned++;
      log(`${m.kind} at ${pinKick.kick.centre}: pinned (local non-6 ${pinKick.before} → ${pinKick.held} while held)`);
    }
  }
  return { motifs: counts, tried: tried.length, healed, pinned, results };
}

/** re-apply the pins' kicks to the solver's field (after a batch of steps) */
export function applyPins(solver, pins) {
  if (!pins || !pins.length) return;
  const u = Float64Array.from(solver.u); for (const p of pins) applyKick(u, p.kick); writeU(solver, u);
}

/** write a field into the solver without resetting its step count */
export function writeU(solver, u) { if (solver.setU) solver.setU(u); else solver.u.set(u); }

/** put the solver back to a saved field and step count (dropping the snapshots taken since) */
export function rewind(solver, u0, steps0, snaps0) { writeU(solver, u0); solver.steps = steps0; solver.snapshots.length = snaps0; }
