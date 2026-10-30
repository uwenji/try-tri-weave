// Mesh file → carrier data, entirely in the browser (the pipeline of export_carrier.py / carrier.py Carrier):
// read → weld → duplicate faces → T-junction cracks → orient → manifold corner table → loose parts → planar caps →
// small holes (< 3 h) → isotropic remesh at h (features kept) → centre → outward side → angle defects →
// intrinsic Delaunay (edge flips on lengths) → quality / meta. The result is what an .lshc file holds.
import { readMesh } from './meshio.js';
import { weld, dedupeFaces, repairCracks, orientFaces, manifoldTable, dropSmallParts, stripPlanarCaps, fillSmallHoles, bounds, faceAreas, cutPaths, zipCut } from './meshrepair.js';
import { remeshIsotropic, relaxRegion } from './remesh.js';
import { flipToDelaunay, angleDefects, orientSide, edgeTable } from './carrier.js';

export const IMPORT_DEFAULTS = { edges: 120, maxVerts: 26000, creaseDeg: 40, stripCaps: true, iters: 10, up: '+z' };

// Up axis of the file → +z (the page is z-up: Rhino / CAD). Proper rotations, so the faces keep their orientation;
// '+y' is the y-up export of glTF, Blender and most web models: (x, y, z) → (x, −z, y), the front stays the front.
export const UP_ROT = {
  '+z': [1, 0, 0, 0, 1, 0, 0, 0, 1], '-z': [1, 0, 0, 0, -1, 0, 0, 0, -1],
  '+y': [1, 0, 0, 0, 0, -1, 0, 1, 0], '-y': [1, 0, 0, 0, 0, 1, 0, -1, 0],
  '+x': [0, 0, -1, 0, 1, 0, 1, 0, 0], '-x': [0, 0, 1, 0, 1, 0, -1, 0, 0],
};
/** rotate positions in place by the row-major 3×3 matrix M */
export function rotatePoints(P, M) {
  for (let i = 0; i < P.length; i += 3) {
    const x = P[i], y = P[i + 1], z = P[i + 2];
    P[i] = M[0] * x + M[1] * y + M[2] * z; P[i + 1] = M[3] * x + M[4] * y + M[5] * z; P[i + 2] = M[6] * x + M[7] * y + M[8] * z;
  }
  return P;
}
/** the rotation taking a carrier made with up `from` to up `to`: R_to · R_fromᵀ */
export function upChange(from, to) {
  const A = UP_ROT[to] || UP_ROT['+z'], B = UP_ROT[from] || UP_ROT['+z'], M = new Array(9).fill(0);
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) M[3 * i + j] += A[3 * i + k] * B[3 * j + k];
  return M;
}

export function triangleQuality(V, F) {
  const m = F.length / 3, E = new Float64Array(3 * m); let amin = 180, skinny = 0;
  const ang = (x, y, z) => Math.acos(Math.max(-1, Math.min(1, (y * y + z * z - x * x) / Math.max(2 * y * z, 1e-300)))) * 180 / Math.PI;
  for (let f = 0; f < m; f++) {
    const p = 3 * F[3 * f], q = 3 * F[3 * f + 1], r = 3 * F[3 * f + 2];
    const a = Math.hypot(V[q] - V[p], V[q + 1] - V[p + 1], V[q + 2] - V[p + 2]);
    const b = Math.hypot(V[r] - V[q], V[r + 1] - V[q + 1], V[r + 2] - V[q + 2]);
    const c = Math.hypot(V[p] - V[r], V[p + 1] - V[r + 1], V[p + 2] - V[r + 2]);
    E[3 * f] = a; E[3 * f + 1] = b; E[3 * f + 2] = c;
    const t = Math.min(ang(a, b, c), ang(b, c, a), ang(c, a, b)); if (t < amin) amin = t; if (t < 15) skinny++;
  }
  const s = E.slice().sort();
  return { edge_median: s[s.length >> 1], min_angle: amin, skinny_pct: 100 * skinny / Math.max(m, 1) };
}

/** buf: file contents, name: file name. opts: {edges | h, maxVerts, creaseDeg, stripCaps, iters}. progress(frac, text).
 *  Returns carrier data {name, meta, V, F, Fi, Li, K, B, sideChecked} (V centred, F outward). */
export function prepareCarrier(buf, name, opts = {}, progress = () => {}) {
  const o = { ...IMPORT_DEFAULTS, ...opts }, notes = [], t0 = Date.now();
  const say = (f, t) => progress(f, t);
  say(0.01, 'reading');
  const raw = readMesh(buf, name);
  if (o.up !== '+z' && UP_ROT[o.up]) rotatePoints(raw.P, UP_ROT[o.up]);
  const bb0 = bounds(raw.P); if (!(bb0.diag > 0)) throw new Error('the mesh has no extent');
  let h = o.h || bb0.diag / o.edges;
  say(0.04, `welding ${raw.faces} triangles`);
  let { V, F } = weld(raw.P, raw.F);
  const dd = dedupeFaces(F, V.length / 3); F = dd.F; if (dd.removed) notes.push(`removed ${dd.removed} duplicate faces`);
  say(0.07, 'closing cracks');
  const cr = repairCracks(V, F); F = cr.F; if (cr.splits) notes.push(`closed cracks: ${cr.splits} T-junction splits`);
  say(0.09, 'orienting');
  // a non-orientable surface (a Möbius band) has no consistent winding: it is cut open along a short edge path for the
  // repair and the remesh (both need an oriented corner table) and glued back after the remesh
  const or = orientFaces(F, V.length / 3, V); F = or.F;
  let cut = null;
  if (or.conflicts) {
    const paths = cutPaths(V, or.cut);
    cut = { paths, edges: or.conflicts };
    notes.push(`non-orientable: cut along ${or.conflicts} edges (${or.cutLength.toPrecision(3)} long, shortened from ${or.cutLength0.toPrecision(3)}) for the remesh`);
  }
  let mt = manifoldTable(V, F); V = mt.V; F = mt.F; let O = mt.O;
  if (mt.nonManifoldEdges) notes.push(`cut ${mt.nonManifoldEdges} non-manifold edges`);
  if (mt.splitVerts) notes.push(`split ${mt.splitVerts} non-manifold vertices`);
  const sp = dropSmallParts(V, F, O, 4 * h * h);
  if (sp.dropped) { notes.push(`dropped ${sp.dropped} loose part${sp.dropped > 1 ? 's' : ''} smaller than 4 h²`); ({ V, F } = sp); ({ V, F, O } = manifoldTable(V, F)); }
  if (o.stripCaps) {
    const sc = stripPlanarCaps(V, F, O);
    if (sc.removed) { notes.push(`stripped ${sc.removed} planar cap faces`); ({ V, F } = sc); ({ V, F, O } = manifoldTable(V, F)); }
  }
  const base = { V, F, O }, inputVertices = V.length / 3;
  let area = 0; for (const a of faceAreas(V, F)) area += a;
  if (!(area > 1e-9 * bb0.diag * bb0.diag) || F.length < 3) throw new Error('the mesh has no area (no faces left after the repair)');
  // vertex cap: a regular mesh at edge h has ≈ 1.2 A / h² vertices (the zoo carriers: 1.15–1.38)
  const nEst = 1.2 * area / (h * h);
  if (nEst > o.maxVerts) h *= Math.sqrt(nEst / o.maxVerts) * 1.02;
  let res, filled = 0, remIn = null;
  for (let attempt = 0; attempt < 4; attempt++) {
    ({ V, F, O } = base);
    const fh = fillSmallHoles(V, F, O, 3 * h); filled = fh.filled;
    if (fh.filled) { ({ V, F, O } = manifoldTable(fh.V, fh.F)); }
    say(0.12, `remeshing at h = ${h.toPrecision(3)}`);
    // the cut's vertices slide along it (its zigzag along the input edges makes no corners); its ends on the rim stay
    let fixed = null, soft = null;
    if (cut && cut.paths) {
      const tag = new Map(), key = (P, i) => `${P[3 * i]},${P[3 * i + 1]},${P[3 * i + 2]}`;
      for (const p of cut.paths) { const k = p.P.length / 3; for (let i = 0; i < k; i++) tag.set(key(p.P, i), !p.closed && (i === 0 || i === k - 1) ? 2 : 1); }
      const nv = V.length / 3; fixed = new Uint8Array(nv); soft = new Uint8Array(nv);
      for (let v = 0; v < nv; v++) { const t = tag.get(key(V, v)); if (t === 2) fixed[v] = 1; else if (t === 1) soft[v] = 1; }
    }
    remIn = { V, F };
    res = remeshIsotropic(V, F, h, { O, iters: o.iters, creaseDeg: o.creaseDeg, fixed, soft, progress: (f, t) => say(0.12 + 0.78 * f, t) });
    const n = res.V.length / 3;
    if (n <= o.maxVerts) break;
    h *= Math.sqrt(n / o.maxVerts) * 1.05;
  }
  if (filled) notes.push(`filled ${filled} boundary loops shorter than 3 h`);
  if (h > (o.h || bb0.diag / o.edges) * 1.001) notes.push(`h coarsened to ${h.toPrecision(3)} for the ${o.maxVerts} vertex cap`);
  V = res.V; F = res.F;
  if (cut && !cut.paths) notes.push('the orientation cut branches: kept open');
  else if (cut) {
    say(0.91, 're-gluing the orientation cut');
    const open = cut.paths.filter(p => !p.closed);
    const z = open.length === cut.paths.length ? zipCut(V, F, open, h) : { error: 'a closed cut' };
    if (z.error) notes.push(`the orientation cut was kept open (${z.error})`);
    else {
      // smooth the seam and two rings around it (the merged sampling of its two sides), back onto the input surface
      const nz = z.V.length / 3, rank = new Uint8Array(nz), onB = new Uint8Array(nz), move = new Uint8Array(nz);
      for (let i = 0; i < z.idx.length; i++) if (z.idx[i] >= 0) rank[z.idx[i]] = res.rank[i];
      const ec = new Map(); for (let f = 0; f < z.F.length / 3; f++) for (let k = 0; k < 3; k++) { const a = z.F[3 * f + k], b = z.F[3 * f + (k + 1) % 3], kk = a < b ? a * 4294967296 + b : b * 4294967296 + a; ec.set(kk, (ec.get(kk) || 0) + 1); }
      for (const [kk, c] of ec) if (c === 1) { onB[Math.floor(kk / 4294967296)] = 1; onB[kk % 4294967296] = 1; }
      let ring = Uint8Array.from(z.seam);
      for (let r = 0; r < 2; r++) { const nx = Uint8Array.from(ring); for (let f = 0; f < z.F.length / 3; f++) { const t = [z.F[3 * f], z.F[3 * f + 1], z.F[3 * f + 2]]; if (t.some(v => ring[v])) for (const v of t) nx[v] = 1; } ring = nx; }
      for (let v = 0; v < nz; v++) move[v] = ring[v] && !onB[v] && (z.seam[v] || rank[v] === 0) ? 1 : 0;
      relaxRegion(z.V, z.F, move, remIn.V, remIn.F);
      V = z.V; F = z.F;
      notes.push(`re-glued the cut after the remesh (${z.nSeam} seam vertices): a non-orientable carrier`);
    }
  }
  const n = V.length / 3;
  // centre (bounding-box centre, as carrier.py), outward side
  const bb = bounds(V), cen = [0, 1, 2].map(k => 0.5 * (bb.lo[k] + bb.hi[k]));
  for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) V[3 * i + k] -= cen[k];
  say(0.92, 'intrinsic Delaunay');
  F = orientSide(V, Uint32Array.from(F));
  const { K, boundary } = angleDefects(V, F, n);
  const ne = edgeTable(F, n).ne, chi = n - ne + F.length / 3;
  let omega = 0; for (let i = 0; i < n; i++) if (!boundary[i]) omega += K[i];
  const { Fi, Li, flips } = flipToDelaunay(V, F);
  const q = triangleQuality(V, F);
  let areaOut = 0; for (const a of faceAreas(V, F)) areaOut += a;
  const st = res.stats, last = st.iters[st.iters.length - 1];
  const remesh = `remeshed in the browser in ${((Date.now() - t0) / 1000).toFixed(1)} s (${st.chains} feature line${st.chains === 1 ? '' : 's'}, ${st.fixed} fixed point${st.fixed === 1 ? '' : 's'})`;
  void last;
  const meta = {
    chi, omega_deg: omega * 180 / Math.PI, h, bbox: bb.diag, area: areaOut, n, nf: F.length / 3, flips,
    min_angle: q.min_angle, skinny_pct: q.skinny_pct, input_vertices: inputVertices, notes, a_default: bb0.diag / 13,
    source: name, remesh, build_ms: Date.now() - t0, input_faces: raw.faces, options: { edges: o.h ? null : o.edges, maxVerts: o.maxVerts, creaseDeg: o.creaseDeg, stripCaps: o.stripCaps, up: o.up },
  };
  say(1, `${n} vertices`);
  return { name: name.replace(/\.[^.]+$/, ''), meta, V, F, Fi, Li, K, B: boundary, sideChecked: true, stats: st };
}
