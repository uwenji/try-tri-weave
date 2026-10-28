// Healing bound 5–7 pairs in the SH field.
// A 5-7-5 whose two 5s sit in the 7's link with two regular sites between them is a pentagon of sites with its centre
// missing: the 7 holds both diagonals of the pentagon, its two neighbours on the pentagon have lost an edge. Putting a
// site at the pentagon's centre gives one 5 there and six corners of valence 6 (net charge +1 either way).
// The mirror, a 7-5-7 whose two 7s are separated by a single site B in the 5's link, is a heptagon holding two sites
// where one belongs: merging the 5 with B gives one 7 at the centre.
// Both are done inside the solver, not on the mesh: a bump (trough) is written into u at the centre and held for a few
// steps while the neighbourhood relaxes, then the field runs free and the lattice is read again; a motif that does not
// resolve is reverted. A quad with the wrong diagonal (5-7-5 with one site between the 5s, 7-7-5 and 5-7-7-5) is
// flipped by moving its four sites. A 5-5-7 triangle is healed by merging its two 5s: the merged site has the union
// of the two links (6), the 7 loses two neighbours and gains the merged one (6), the other common neighbour of the two
// 5s loses two and gains one (5) — the 5 has moved one site and the dislocation is gone, net charge +1 as before.
// Free 5–7 pairs (dislocations) cannot be removed by any single local move: an edge flip conserves the Burgers vector
// and only glides the pair by one lattice vector. glidePass chains such flips, inside the solver, to carry a free pair
// to its anti-pair (annihilation), to a lone 5 or 7 (absorption through the motifs above) or to the rim.
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
            // 5-7-7-5: the far corner is the other pair's 5 — flipping the 7–7 edge makes four 6s (two anti-pairs annihilate)
            else if (cand.length === 1 && sing(cand[0], 5)) take({ kind: 'flip', peaks: [v, p, q, cand[0]], corners: [v, p, cand[0], q], edge: [p, q], opp: [v, cand[0]], annihilate: true });
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
  if (motif.kind === 'merge' || motif.kind === 'tri') {
    // merge: the 5 and the regular site between the 7s; tri: the two 5s of the 5-5-7 triangle (motif.peaks = [7, 5, 5])
    const cs = motif.kind === 'tri' ? [motif.peaks[1], motif.peaks[2]] : motif.corners;
    const v = site !== null ? site : motifCentre(lat, c, motif, cs), [A, B] = cs.map(p => lat.peaks[p]);
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

/** vertices where the surface is flat at the lattice scale: less than `maxDeg` of curvature per lattice cell, the angle
 *  defects smoothed over a/2 (cached on the carrier per a) — the K = 0 areas, where the lattice should be regular */
export function flatMask(c, a, maxDeg = 3) {
  const key = `${a.toPrecision(5)}|${maxDeg}`;
  if (c._flatMask && c._flatMask.key === key) return c._flatMask.v;
  const rho = 0.5 * a, u = c.heatSolve(Float64Array.from(c.K, (k, i) => c.boundary[i] ? 0 : k), 0.5 * rho * rho), cell = Math.sqrt(3) / 2 * a * a, out = new Uint8Array(c.n);
  for (let i = 0; i < c.n; i++) out[i] = !c.boundary[i] && Math.abs(u[i] * cell) * 180 / Math.PI < maxDeg ? 1 : 0;
  c._flatMask = { key, v: out }; return out;
}

/** lattice regularity around carrier vertex v within geodesic radius R: rms edge strain (chord / a − 1) over the dual
 *  edges and rms corner-angle deviation from 60° over the dual triangles whose peaks all lie inside; `irr` = strain +
 *  angle / 60°, over the whole disc, and `flat` = the same over the triangles whose three peaks stand on flat ground
 *  (flatMask: K ≈ 0 at the lattice scale) — the part where a weave must not twist or bend; 0 when the disc has none */
export function localRegularity(lat, c, v, R, a) {
  const disc = discAround(c, v, R), pos = new Int32Array(c.n).fill(-1); lat.peaks.forEach((p, i) => { pos[p] = i; });
  const inside = new Uint8Array(lat.peaks.length); for (const i of disc.keys()) { const p = pos[i]; if (p >= 0) inside[p] = 1; }
  const fm = flatMask(c, a), T = lat.tris, P = lat.P, seen = new Set(), len = (p, q) => Math.hypot(P[3 * p] - P[3 * q], P[3 * p + 1] - P[3 * q + 1], P[3 * p + 2] - P[3 * q + 2]);
  const acc = [[0, 0, 0, 0], [0, 0, 0, 0]];                                          // [all, flat] × [Σ strain², edges, Σ angle², angles]
  for (let t = 0; t < T.length / 3; t++) {
    const i = T[3 * t], j = T[3 * t + 1], k = T[3 * t + 2]; if (!inside[i] || !inside[j] || !inside[k]) continue;
    const flat = fm[lat.peaks[i]] && fm[lat.peaks[j]] && fm[lat.peaks[k]], lij = len(i, j), ljk = len(j, k), lki = len(k, i);
    for (const [p, q, l] of [[i, j, lij], [j, k, ljk], [k, i, lki]]) { const key = p < q ? p * 1e5 + q : q * 1e5 + p; if (seen.has(key)) continue; seen.add(key); const e2 = (l / a - 1) ** 2; acc[0][0] += e2; acc[0][1]++; if (flat) { acc[1][0] += e2; acc[1][1]++; } }
    for (const [x, y, z] of [[lij, lki, ljk], [lij, ljk, lki], [ljk, lki, lij]]) { const cs = Math.max(-1, Math.min(1, (x * x + y * y - z * z) / (2 * x * y))), d2 = (Math.acos(cs) - Math.PI / 3) ** 2; acc[0][2] += d2; acc[0][3]++; if (flat) { acc[1][2] += d2; acc[1][3]++; } }
  }
  const irrOf = (q) => q[1] ? Math.sqrt(q[0] / q[1]) + Math.sqrt(q[2] / Math.max(q[3], 1)) / (Math.PI / 3) : 0;
  return { strain: Math.sqrt(acc[0][0] / Math.max(acc[0][1], 1)), angle: Math.sqrt(acc[0][2] / Math.max(acc[0][3], 1)), irr: irrOf(acc[0]), flat: irrOf(acc[1]), edges: acc[0][1], flatEdges: acc[1][1] };
}
/** the regularity rule of every pass (tol = 0.02 by default): the flat part (K ≈ 0, where laths must not twist or bend)
 *  no worse than 2.5·tol — it is a small sample of a few triangles, so a little slack, ≈ 3 % strain or 2° angle rms — and
 *  the whole disc, crease and defect cores included, no worse than 4·tol */
const regOK = (r0, r1, tol) => r1.flat <= r0.flat + 2.5 * tol && r1.irr <= r0.irr + 4 * tol;
const regStr = (r0, r1) => `irregularity ${r0.irr.toFixed(3)} → ${r1.irr.toFixed(3)}` + (r0.flatEdges || r1.flatEdges ? ` (flat ${r0.flat.toFixed(3)} → ${r1.flat.toFixed(3)})` : '');

const nonSix = (lat) => { let k = 0; for (let i = 0; i < lat.valence.length; i++) if (lat.interior[i] && lat.valence[i] !== 6) k++; return k; };

/**
 * healPass(solver, carrier, a, opts) → report. Reads the lattice, finds the motifs, and heals the insert / merge ones one
 * after the other: kick → hold (the kick re-applied every `every` steps for `hold` steps) → `release` free steps → read
 * again; kept when the defects within 2.5 a of the centre went down and the total did not go up, else the field is put
 * back. opts.log(msg), opts.onStep(solver) (called after each stepping batch, for the UI), opts.max (motifs per pass).
 */
export async function healPass(solver, c, a, { hold = 60, every = 5, release = 400, verify = 500, max = 50, log = () => {}, onStep = null, kinds = ['insert', 'merge', 'flip', 'tri'], readOpts = {}, bothSides = true, push = 0.4, pin = false, pins = null, regTol = 0.02 } = {}) {
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
      const ss = q.kind === 'insert' ? insertSites(lat, c, q, a).slice(0, bothSides ? 2 : 1)
        : q.kind === 'tri' ? [{ vertex: motifCentre(lat, c, q, [q.peaks[1], q.peaks[2]]), corners: [q.peaks[1], q.peaks[2]], room: 0 }]      // merge the two 5s
        : [{ vertex: motifCentre(lat, c, q), corners: q.corners, room: 0 }];   // flip / merge: one attempt
      if (!ss.length || ss[0].vertex < 0 || near(ss[0].vertex)) continue;
      m = q; sites = ss; break;
    }
    if (!m) break;
    const amp = peakAmplitude(lat, solver.u);
    const u0 = Float64Array.from(solver.u), steps0 = solver.steps, snaps0 = solver.snapshots.length, totalBefore = nonSix(lat);
    let ok = false, pinKick = null;
    for (const site of sites) {
      const kick = motifKick(lat, c, m, a, amp, { site: site.vertex, push }); if (!kick) break;
      const before = localDefects(lat, c, kick.centre, 2.5 * a), reg0 = localRegularity(lat, c, kick.centre, 2.5 * a, a);
      const uk = Float64Array.from(solver.u); applyKick(uk, kick); writeU(solver, uk);
      for (let s = 0; s < hold; s += every) { await step(every); const uu = Float64Array.from(solver.u); applyKick(uu, kick); writeU(solver, uu); }
      // the held state: what a pin would keep — only an insert / merge whose held lattice is already cleaner qualifies
      const latH = extractLattice(c, solver.u, a, readOpts), held = latH ? localDefects(latH, c, kick.centre, 2.5 * a) : before;
      const regH = latH ? localRegularity(latH, c, kick.centre, 2.5 * a, a) : reg0;
      // a pin holds this state for good, so it must be cleaner *and* no less regular than what was there
      if (!pinKick && m.kind !== 'flip' && held.count < before.count && held.charge === before.charge && regOK(reg0, regH, regTol)) pinKick = { kick, held: held.count, before: before.count, u: Float64Array.from(solver.u), steps: solver.steps, snaps: solver.snapshots.length };
      await step(release);
      const lat2 = extractLattice(c, solver.u, a, readOpts), after = lat2 ? localDefects(lat2, c, kick.centre, 2.5 * a) : before, totalAfter = lat2 ? nonSix(lat2) : totalBefore;
      const reg1 = lat2 ? localRegularity(lat2, c, kick.centre, 2.5 * a, a) : reg0;
      // fewer defects around it, none anywhere else, the charge kept, and the lattice around it no less regular (edge lengths and angles)
      ok = after.count < before.count && totalAfter <= totalBefore && after.charge === before.charge && regOK(reg0, reg1, regTol);
      let verdict = ok ? 'kept' : 'reverted';
      if (ok && pin && pins && pinKick && verify > 0) {
        // a kick that holds through `release` may still undo itself later (the pair re-forms): check, and pin it instead
        await step(verify);
        const lat3 = extractLattice(c, solver.u, a, readOpts), later = lat3 ? localDefects(lat3, c, kick.centre, 2.5 * a) : after;
        if (later.count > after.count || (lat3 && nonSix(lat3) > totalAfter)) { ok = false; verdict = `undone within ${verify} more steps → pinned`; }
      }
      log(`${m.kind} at ${kick.centre} (${site.corners.length}-gon, room ${site.room.toFixed(2)} a): local non-6 ${before.count} → ${held.count} held → ${after.count} free, total ${totalBefore} → ${totalAfter}, ${regStr(reg0, reg1)} · ${verdict}`);
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

// ------------------------------------------------------------------ free 5–7 pairs: glide
/** free (unbound) 5–7 pairs: a 7 whose link holds exactly one 5 and no other non-6, that 5 holding no other non-6 either.
 *  Returns [{b: 5, c: 7}] as peak indices. */
export function freePairs(lat, links = peakLinks(lat)) {
  const out = [], val = lat.valence, inter = lat.interior;
  for (let v = 0; v < lat.peaks.length; v++) {
    if (!inter[v] || val[v] !== 7 || !links[v]) continue;
    const L = links[v], fives = L.filter(p => inter[p] && val[p] === 5), others = L.filter(p => inter[p] && val[p] !== 6 && val[p] !== 5);
    if (fives.length !== 1 || others.length) continue;
    const b = fives[0], Lb = links[b]; if (!Lb || Lb.some(p => p !== v && inter[p] && val[p] !== 6)) continue;
    out.push({ b, c: v });
  }
  return out;
}

/** Burgers vector of a free pair as a 3-D tangent vector: a 5–7 dislocation in a triangular lattice has |b| = a,
 *  perpendicular to its 5→7 bond — b = n × (r₇ − r₅), n the mean normal of the two sites. Two dislocations annihilate
 *  only when their Burgers vectors are opposite (a 5→7 bond pointing the other way); with b at 120° they merge into one
 *  dislocation, parallel ones never meet. Compared between nearby pairs by the 3-D dot product. */
export function burgers(lat, c, pr) {
  const V = c.V, N = c.N, B = lat.peaks[pr.b], C = lat.peaks[pr.c];
  const d = [V[3 * C] - V[3 * B], V[3 * C + 1] - V[3 * B + 1], V[3 * C + 2] - V[3 * B + 2]];
  const n = [N[3 * B] + N[3 * C], N[3 * B + 1] + N[3 * C + 1], N[3 * B + 2] + N[3 * C + 2]], nl = Math.hypot(...n) || 1;
  n[0] /= nl; n[1] /= nl; n[2] /= nl;
  return [n[1] * d[2] - n[2] * d[1], n[2] * d[0] - n[0] * d[2], n[0] * d[1] - n[1] * d[0]];
}
/** cosine of the angle between the Burgers vectors of two free pairs: −1 for an annihilating (anti-parallel) couple */
export function burgersCos(lat, c, p1, p2) {
  const b1 = burgers(lat, c, p1), b2 = burgers(lat, c, p2);
  return (b1[0] * b2[0] + b1[1] * b2[1] + b1[2] * b2[2]) / ((Math.hypot(...b1) * Math.hypot(...b2)) || 1);
}

/** lone non-6 sites of valence `val` (no other non-6 in the link): the absorbing targets */
export function loneSites(lat, links, val) {
  const out = [];
  for (let v = 0; v < lat.peaks.length; v++) if (lat.interior[v] && lat.valence[v] === val && links[v] && !links[v].some(p => lat.interior[p] && lat.valence[p] !== 6)) out.push(v);
  return out;
}

/** curvature the surface holds within radius R of carrier vertex v, in units of 60° (the lattice-scale demand for a 5 or a 7) */
export function discCharge(c, v, R) {
  let s = 0; for (const [i] of discAround(c, v, R)) if (!c.boundary[i]) s += c.K[i];
  return s / (Math.PI / 3);
}

/** the two glide moves of a free pair: flip 7–x → 5–y, with x a far corner of the 5–7 edge (next to the 5 in the 7's link)
 *  and y the vertex after x. After the flip the 7 and the 5 are 6s, x is the 5 and y the 7: the pair has moved by one
 *  lattice vector. Returns [{x, y, motif}] (a flip motif for motifKick); moves through non-regular sites are left out. */
export function glideMoves(lat, links, pr) {
  const L = links[pr.c], k = L.length, i = L.indexOf(pr.b), out = [];
  if (i < 0) return out;
  for (const dir of [1, -1]) {
    const x = L[(i + dir + k) % k], y = L[(i + 2 * dir + k) % k];
    if (x === pr.b || y === pr.b || lat.valence[x] !== 6 || lat.valence[y] !== 6 || !lat.interior[x] || !lat.interior[y]) continue;
    out.push({ x, y, motif: { kind: 'flip', edge: [pr.c, x], opp: [pr.b, y], corners: [pr.c, pr.b, x, y], peaks: [pr.c, pr.b] } });
  }
  return out;
}

const P3 = (c, i) => [c.V[3 * i], c.V[3 * i + 1], c.V[3 * i + 2]];
const dist3 = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
const pairMid = (lat, c, pr) => { const B = P3(c, lat.peaks[pr.b]), C = P3(c, lat.peaks[pr.c]); return [(B[0] + C[0]) / 2, (B[1] + C[1]) / 2, (B[2] + C[2]) / 2]; };

/**
 * glidePass(solver, c, a, opts) → report. Reads the lattice, lists the free 5–7 pairs and glides each one flip at a time
 * towards its target (with needDeg > 0 a pair under which the surface holds ≥ needDeg of curvature within 1.2 a is left
 * alone — off by default: on a crease the disc test is unreliable, and the move acceptance already refuses harm) —
 * target — the nearest other free pair (anti-pairs annihilate), else the nearest lone 5 or 7 (then a 5-5-7, 5-7-5, 7-5-7
 * or 5-7-7 forms, which healPass takes apart), else the nearest rim vertex on an open surface (the rim absorbs it).
 * Each move: kick → hold → `release` free steps → read; kept when nothing new appeared around it and the pair count, the
 * total or the distance to the target went down, else the field is put back. Once a pair sits by its target the field
 * runs `settle` steps and healPass runs on the motifs that formed. opts: hold, every, release, settle, maxMoves (per
 * pair), maxPairs, needDeg, log, onStep, pin, pins, readOpts.
 */
export async function glidePass(solver, c, a, { hold = 60, every = 5, release = 150, settle = 300, maxMoves = 12, maxPairs = 20, needDeg = 0, log = () => {}, onStep = null, pin = false, pins = null, readOpts = {}, regTol = 0.02 } = {}) {
  const step = async (k) => { await solver.step(k); if (onStep) await onStep(solver); };
  const read = () => { const lat = extractLattice(c, solver.u, a, readOpts); return lat ? { lat, links: peakLinks(lat), non6: nonSix(lat) } : null; };
  const open = c.boundary.some(b => b);
  const rimOf = (v) => { let best = null; for (const [i, d] of discAround(c, v, 6 * a)) if (c.boundary[i] && (!best || d < best.d)) best = { i, d }; return best; };
  const report = { pairs: 0, skipped: 0, annihilated: 0, absorbed: 0, toRim: 0, left: 0, moves: 0, log: [] };
  let r = read(); if (!r) return report;
  const need = (pr) => Math.abs(discCharge(c, r.lat.peaks[pr.c], 1.2 * a)) * 60;     // what the curvature under the pair asks for, in degrees
  const all = freePairs(r.lat, r.links); report.pairs = all.length;
  let anti = 0; for (let i = 0; i < all.length; i++) for (let j = i + 1; j < all.length; j++) if (burgersCos(r.lat, c, all[i], all[j]) < -0.5) anti++;
  report.antiCouples = anti;
  log(`free 5–7 pairs: ${all.length}` + (all.length ? ` (${anti} anti-parallel couple${anti === 1 ? '' : 's'} by their Burgers vectors) · curvature under them: ` + all.map(pr => `${need(pr).toFixed(0)}°`).join(', ') : ''));
  const done = new Set();                                                  // carrier vertices of pairs already handled
  for (let n = 0; n < maxPairs; n++) {
    r = read(); if (!r) break;
    const pairs = freePairs(r.lat, r.links), free = pairs.filter(pr => !done.has(r.lat.peaks[pr.c]));
    if (!free.length) break;
    const A = free[0], keyA = r.lat.peaks[A.c]; done.add(keyA);
    const demand = need(A);
    if (needDeg > 0 && demand >= needDeg) { report.skipped++; log(`pair at ${keyA}: the surface holds ${demand.toFixed(0)}° under it — left in place`); continue; }
    let outcome = 'left', pos = pairMid(r.lat, c, A);                      // the pair is followed by position: peaks move between read-outs
    const seq = { u: Float64Array.from(solver.u), steps: solver.steps, snaps: solver.snapshots.length, pins: pins ? pins.length : 0, reg: localRegularity(r.lat, c, keyA, 4 * a, a), non6: r.non6 };
    const findMe = (rr, ps) => { let best = null, bd = 0.9 * a; for (const pr of ps) { const d = dist3(pairMid(rr.lat, c, pr), pos); if (d < bd) { bd = d; best = pr; } } return best; };
    for (let mv = 0; mv < maxMoves; mv++) {
      r = read(); if (!r) break;
      const pairsNow = freePairs(r.lat, r.links), me = findMe(r, pairsNow);
      if (!me) {
        // the pair is no longer free: it met its target (a motif formed) or it is gone
        const before = r.non6, motifs = findMotifs(r.lat).length;
        if (motifs) { const rep = await healPass(solver, c, a, { hold, every, release: 400, max: 6, log: (m) => log('  heal: ' + m), onStep, pin, pins, readOpts }); r = read(); }
        else { await step(settle); r = read(); }
        const pairsAfter = r ? freePairs(r.lat, r.links) : [];
        outcome = (r && findMe(r, pairsAfter)) ? 'left' : (r && r.non6 <= before - 2 ? 'absorbed' : 'left');
        break;
      }
      pos = pairMid(r.lat, c, me); for (const e of [me.b, me.c]) done.add(r.lat.peaks[e]);
      // target: the nearest other free pair whose Burgers vector is opposite to this one's (they annihilate), else the
      // nearest lone 5 / 7 (absorption: the disclination moves one site), else the rim
      const mid = pairMid(r.lat, c, me); let target = null, bd = Infinity, tkind = '';
      for (const q of pairsNow) { if (q === me || burgersCos(r.lat, c, me, q) > -0.5) continue; const d = dist3(mid, pairMid(r.lat, c, q)); if (d < bd) { bd = d; target = pairMid(r.lat, c, q); tkind = 'anti-pair'; } }
      if (!target) for (const val of [5, 7]) for (const v of loneSites(r.lat, r.links, val)) { const d = dist3(mid, P3(c, r.lat.peaks[v])); if (d < bd) { bd = d; target = P3(c, r.lat.peaks[v]); tkind = `lone ${val}`; } }
      if (!target && open) { const rim = rimOf(r.lat.peaks[me.c]); if (rim) { target = P3(c, rim.i); bd = rim.d; tkind = 'rim'; } }
      if (!target) { log(`pair at ${keyA}: no target`); break; }
      if (tkind === 'anti-pair' && bd < 1.6 * a) {
        // adjacent anti-pairs annihilate on their own (the field's attraction beats the pinning at this range)
        const before = r.non6; await step(settle); r = read();
        const still = r && findMe(r, freePairs(r.lat, r.links));
        if (r && !still) { outcome = r.non6 <= before - 2 ? 'annihilated' : 'absorbed'; log(`pair at ${keyA}: met ${tkind} at ${(bd / a).toFixed(2)} a → ${outcome} (non-6 ${before} → ${r.non6})`); break; }
        if (!r) break;
      }
      const moves = glideMoves(r.lat, r.links, me); if (!moves.length) { log(`pair at ${keyA}: no regular quad to flip`); break; }
      const after = (m) => { const X = P3(c, r.lat.peaks[m.x]), Y = P3(c, r.lat.peaks[m.y]); return dist3([(X[0] + Y[0]) / 2, (X[1] + Y[1]) / 2, (X[2] + Y[2]) / 2], target); };
      moves.sort((m1, m2) => after(m1) - after(m2));
      const amp = peakAmplitude(r.lat, solver.u), lat0 = r.lat, links0 = r.links, non60 = r.non6;
      let moved = false;
      for (const [m, push] of moves.flatMap(m => [[m, 0.35], [m, 0.22]])) {              // the full quad move first, a gentler one if it misfires
        const kick = motifKick(lat0, c, m.motif, a, amp, { push }); if (!kick) continue;
        const u0 = Float64Array.from(solver.u), steps0 = solver.steps, snaps0 = solver.snapshots.length;
        const before = localDefects(lat0, c, kick.centre, 3 * a), reg0 = localRegularity(lat0, c, kick.centre, 3 * a, a);
        const uk = Float64Array.from(solver.u); applyKick(uk, kick); writeU(solver, uk);
        for (let k = 0; k < hold; k += every) { await step(every); const uu = Float64Array.from(solver.u); applyKick(uu, kick); writeU(solver, uu); }
        await step(release);
        const r2 = read(); if (!r2) { rewind(solver, u0, steps0, snaps0); continue; }
        const loc = localDefects(r2.lat, c, kick.centre, 3 * a), p2 = freePairs(r2.lat, r2.links);
        // where the pair went: the free pair nearest to the quad's centre (the 5 is now at x, the 7 at y)
        let me2 = null, md = 1.3 * a; for (const pr of p2) { const d = dist3(pairMid(r2.lat, c, pr), P3(c, kick.centre)); if (d < md) { md = d; me2 = pr; } }
        let nearer = p2.length < pairsNow.length || r2.non6 < non60;
        const d2 = me2 ? dist3(pairMid(r2.lat, c, me2), target) : NaN;
        if (!nearer && me2) nearer = d2 < bd - 0.1 * a;
        const reg1 = localRegularity(r2.lat, c, kick.centre, 3 * a, a);
        const ok = loc.count <= before.count && loc.charge === before.charge && r2.non6 <= non60 && nearer;       // regularity is judged on the whole sequence, below
        log(`  pair at ${keyA}: flip 7 ${lat0.peaks[me.c]}–${lat0.peaks[m.x]} → 5 ${lat0.peaks[me.b]}–${lat0.peaks[m.y]} (${push} a) towards ${tkind} (${(bd / a).toFixed(2)} → ${me2 ? (d2 / a).toFixed(2) : '—'} a): local ${before.count} → ${loc.count}, total ${non60} → ${r2.non6}, pairs ${pairsNow.length} → ${p2.length}, ${regStr(reg0, reg1)} · ${ok ? 'kept' : 'reverted'}`);
        if (ok) { moved = true; report.moves++; if (me2) pos = pairMid(r2.lat, c, me2); break; }
        rewind(solver, u0, steps0, snaps0);
      }
      if (!moved) break;
    }
    // the whole sequence (moves + the motif heal) must leave the lattice around the start no less regular than it found it
    r = read();
    if (r) {
      const reg1 = localRegularity(r.lat, c, keyA, 4 * a, a);
      if (!regOK(seq.reg, reg1, regTol)) {
        rewind(solver, seq.u, seq.steps, seq.snaps); if (pins) pins.length = seq.pins;
        log(`pair at ${keyA}: ${outcome}, but the lattice around it got less regular (${regStr(seq.reg, reg1)}) — put back`); outcome = 'left';
      } else log(`pair at ${keyA}: ${outcome} · ${regStr(seq.reg, reg1)}, non-6 ${seq.non6} → ${r.non6}`);
    }
    if (outcome === 'annihilated') report.annihilated++; else if (outcome === 'absorbed') report.absorbed++; else if (outcome === 'toRim') report.toRim++; else report.left++;
  }
  r = read(); report.freeAfter = r ? freePairs(r.lat, r.links).length : 0; report.non6 = r ? r.non6 : null;
  return report;
}

// ------------------------------------------------------------------ lone 5s and 7s: climb to the curvature
/**
 * climbPass(solver, c, a, opts) → report. A lone 5 (7) belongs where the surface holds its +60° (−60°): on a crease, a
 * tip, a dome (a saddle, a neck). One that sits on flat ground — less than `needDeg` (20°) of its sign within 1.2 a — is
 * walked one site at a time towards the nearest vertex whose disc of 1.2 a holds at least `fullDeg` (40°) of that sign:
 * a trough on its peak and a bump one lattice edge along the way (the kick that moves the 7 in the insert heal). A step
 * is kept when a lone site of the same valence now stands by the bump, the curvature under it went up, no defect
 * appeared within 3 a and the lattice there is no less regular; else the field is put back. A site that will not hold
 * its new place (undone within `verify` steps) is pinned when `pin` is set.
 */
export async function climbPass(solver, c, a, { hold = 60, every = 5, release = 400, verify = 500, maxSteps = 6, maxSites = 20, needDeg = 20, fullDeg = 40, log = () => {}, onStep = null, pin = false, pins = null, readOpts = {}, regTol = 0.02 } = {}) {
  const step = async (k) => { await solver.step(k); if (onStep) await onStep(solver); };
  const read = () => { const lat = extractLattice(c, solver.u, a, readOpts); return lat ? { lat, links: peakLinks(lat), non6: nonSix(lat) } : null; };
  const charge = (v) => discCharge(c, v, 1.2 * a) * 60;                                   // degrees of curvature under a site
  const nearLone = (rr, val, at, R) => { let best = -1, bd = R; for (const v of loneSites(rr.lat, rr.links, val)) { const d = dist3(P3(c, rr.lat.peaks[v]), at); if (d < bd) { bd = d; best = v; } } return best; };
  const report = { sites: 0, moved: 0, pinned: 0, left: 0, steps: 0 };
  let r = read(); if (!r) return report;
  const cand = [];
  for (const val of [5, 7]) for (const v of loneSites(r.lat, r.links, val)) { const q = charge(r.lat.peaks[v]) * (val === 5 ? 1 : -1); if (q < needDeg) cand.push({ v, val, q, at: P3(c, r.lat.peaks[v]) }); }
  report.sites = cand.length;
  log(`lone sites off their curvature: ${cand.length}` + (cand.length ? ' · ' + cand.map(x => `${x.val} at ${r.lat.peaks[x.v]} (${x.q.toFixed(0)}°)`).join(', ') : ''));
  for (let n = 0; n < Math.min(maxSites, cand.length); n++) {
    const { val } = cand[n], sign = val === 5 ? 1 : -1; let pos = cand[n].at, outcome = 'left', stepsHere = 0;
    for (let mv = 0; mv < maxSteps; mv++) {
      r = read(); if (!r) break;
      const me = nearLone(r, val, pos, 0.9 * a); if (me < 0) { outcome = 'lost'; break; }
      const vMe = r.lat.peaks[me], qMe = sign * charge(vMe); pos = P3(c, vMe);
      if (qMe >= fullDeg) { outcome = 'arrived'; break; }
      // target: the nearest vertex (within 6 a) whose disc holds fullDeg of the right sign; then one lattice edge towards it
      let target = -1, td = Infinity;
      for (const [i, d] of discAround(c, vMe, 6 * a)) { if (d < 0.8 * a || d >= td || c.boundary[i]) continue; if (sign * charge(i) >= fullDeg) { target = i; td = d; } }
      if (target < 0) { log(`${val} at ${vMe}: nothing of its sign within 6 a`); break; }
      const T = P3(c, target), dir = [T[0] - pos[0], T[1] - pos[1], T[2] - pos[2]], len = Math.hypot(...dir) || 1, stepLen = Math.min(a, len);
      const X = dir.map((x, k) => pos[k] + stepLen * x / len);
      let site = -1, sd = Infinity; for (const [i] of discAround(c, vMe, 1.5 * a)) { const Q = P3(c, i), d2 = (Q[0] - X[0]) ** 2 + (Q[1] - X[1]) ** 2 + (Q[2] - X[2]) ** 2; if (d2 < sd) { sd = d2; site = i; } }
      if (site < 0 || site === vMe) break;
      const amp = peakAmplitude(r.lat, solver.u), kick = { centre: site, troughs: [bump(c, vMe, 0.3 * a, -0.6 * amp)], bumps: [bump(c, site, 0.42 * a, amp)] };
      const u0 = Float64Array.from(solver.u), steps0 = solver.steps, snaps0 = solver.snapshots.length;
      const before = localDefects(r.lat, c, site, 3 * a), reg0 = localRegularity(r.lat, c, site, 3 * a, a), non60 = r.non6, S3 = P3(c, site);
      const uk = Float64Array.from(solver.u); applyKick(uk, kick); writeU(solver, uk);
      for (let k = 0; k < hold; k += every) { await step(every); const uu = Float64Array.from(solver.u); applyKick(uu, kick); writeU(solver, uu); }
      const rH = read(); let heldState = null;
      if (rH && nearLone(rH, val, S3, 0.9 * a) >= 0) { const ld = localDefects(rH.lat, c, site, 3 * a), rgH = localRegularity(rH.lat, c, site, 3 * a, a); if (ld.count <= before.count && ld.charge === before.charge && regOK(reg0, rgH, regTol)) heldState = { u: Float64Array.from(solver.u), steps: solver.steps, snaps: solver.snapshots.length }; }
      await step(release);
      const r2 = read(); let ok = false, me2 = -1, q2 = NaN, reg1 = reg0, loc = before;
      if (r2) {
        me2 = nearLone(r2, val, S3, 0.9 * a); loc = localDefects(r2.lat, c, site, 3 * a); reg1 = localRegularity(r2.lat, c, site, 3 * a, a);
        q2 = me2 >= 0 ? sign * charge(r2.lat.peaks[me2]) : NaN;
        ok = me2 >= 0 && q2 > qMe + 1 && loc.count <= before.count && loc.charge === before.charge && r2.non6 <= non60 && regOK(reg0, reg1, regTol);
      }
      let verdict = ok ? 'kept' : 'reverted';
      if (ok && pin && pins && heldState && verify > 0) {
        await step(verify); const r3 = read();
        if (!r3 || nearLone(r3, val, S3, 0.9 * a) < 0 || r3.non6 > non60) { ok = false; verdict = `undone within ${verify} more steps → pinned`; }
      }
      log(`${val} at ${vMe} (${qMe.toFixed(0)}° under it) → ${site} towards ${target}: ${me2 >= 0 ? `${val} now at ${r2.lat.peaks[me2]} (${q2.toFixed(0)}°)` : 'no lone ' + val + ' by the new site'}, local non-6 ${before.count} → ${loc.count}, ${regStr(reg0, reg1)} · ${verdict}`);
      if (ok) { report.steps++; stepsHere++; pos = P3(c, r2.lat.peaks[me2]); continue; }
      if (pin && pins && heldState && verdict.startsWith('undone')) { writeU(solver, heldState.u); solver.steps = heldState.steps; solver.snapshots.length = Math.min(solver.snapshots.length, heldState.snaps); pins.push({ kind: 'climb', centre: site, kick }); report.pinned++; report.steps++; stepsHere++; outcome = 'pinned'; break; }
      rewind(solver, u0, steps0, snaps0); break;
    }
    if (stepsHere > 0) report.moved++; else report.left++;
    log(`${val}: ${outcome}${stepsHere ? ` after ${stepsHere} step${stepsHere > 1 ? 's' : ''}` : ''}`);
  }
  r = read(); report.non6 = r ? r.non6 : null;
  return report;
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
