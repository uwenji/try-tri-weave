// The frozen lattice: a read-out taken out of the field. Its sites are points on the carrier, not carrier vertices.
//   freeze     the lattice of a read-out with every site at the peak of u between the vertices
//   relax      moves the sites on the surface until the laths run straight
//   lathBend   the kink of the two laths that cross at every interior lattice edge
//   lineBend   the kink of the three lattice lines through every regular site
//   shape      edge lengths, the largest move, triangles folded over
//   project    the point of the carrier nearest to a point
// A lath passes through the midpoints of the lattice edges. At the edge a–c between the triangles a b c and a c d it
// turns by the angle between c − b and d − a (the other lath there: a − b and d − c), with the second triangle
// unfolded about the edge into the plane of the first. Both turn by the same vector, (b + d − a − c) / 2: the laths
// are straight exactly where the two triangles form a parallelogram, and then the lattice lines are straight too.
// All angles are measured on the lattice itself (its triangles unfolded), so bending with the surface does not count.

const DEG = 180 / Math.PI;

/** vertex → incident faces of the drawn mesh, as CSR; kept on the carrier */
function vertexFaces(c) {
  if (c._vf) return c._vf;
  const F = c.F, n = c.n, ptr = new Int32Array(n + 1);
  for (let i = 0; i < F.length; i++) ptr[F[i] + 1]++;
  for (let i = 0; i < n; i++) ptr[i + 1] += ptr[i];
  const idx = new Int32Array(F.length), fill = ptr.slice(0, n);
  for (let f = 0; f < F.length / 3; f++) for (let k = 0; k < 3; k++) idx[fill[F[3 * f + k]]++] = f;
  return (c._vf = { ptr, idx });
}

const TP = new Float64Array(5);                       // closest point on a triangle: x y z and the weights of b and c
/** closest point of triangle (ia, ib, ic) of V to p (Ericson, Real-Time Collision Detection 5.1.5); returns the squared distance */
function closestOnTriangle(V, ia, ib, ic, px, py, pz) {
  const ax = V[3 * ia], ay = V[3 * ia + 1], az = V[3 * ia + 2], abx = V[3 * ib] - ax, aby = V[3 * ib + 1] - ay, abz = V[3 * ib + 2] - az, acx = V[3 * ic] - ax, acy = V[3 * ic + 1] - ay, acz = V[3 * ic + 2] - az;
  const apx = px - ax, apy = py - ay, apz = pz - az, d1 = abx * apx + aby * apy + abz * apz, d2 = acx * apx + acy * apy + acz * apz;
  let v, w;
  if (d1 <= 0 && d2 <= 0) { v = 0; w = 0; }
  else {
    const bpx = apx - abx, bpy = apy - aby, bpz = apz - abz, d3 = abx * bpx + aby * bpy + abz * bpz, d4 = acx * bpx + acy * bpy + acz * bpz;
    if (d3 >= 0 && d4 <= d3) { v = 1; w = 0; }
    else {
      const vc = d1 * d4 - d3 * d2;
      if (vc <= 0 && d1 >= 0 && d3 <= 0) { v = d1 / (d1 - d3); w = 0; }
      else {
        const cpx = apx - acx, cpy = apy - acy, cpz = apz - acz, d5 = abx * cpx + aby * cpy + abz * cpz, d6 = acx * cpx + acy * cpy + acz * cpz;
        if (d6 >= 0 && d5 <= d6) { v = 0; w = 1; }
        else {
          const vb = d5 * d2 - d1 * d6;
          if (vb <= 0 && d2 >= 0 && d6 <= 0) { v = 0; w = d2 / (d2 - d6); }
          else {
            const va = d3 * d6 - d5 * d4;
            if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) { w = (d4 - d3) / ((d4 - d3) + (d5 - d6)); v = 1 - w; }
            else { const den = 1 / (va + vb + vc); v = vb * den; w = vc * den; }
          }
        }
      }
    }
  }
  const x = ax + v * abx + w * acx, y = ay + v * aby + w * acy, z = az + v * abz + w * acz;
  TP[0] = x; TP[1] = y; TP[2] = z; TP[3] = v; TP[4] = w;
  return (x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2;
}

/** The point of the carrier nearest to (x, y, z): a walk down the vertex graph from `start` to the nearest vertex, then
 *  the closest point on the faces around it (and around its neighbours when that point lies on the far edge of a face).
 *  Writes the point, the surface normal there (on the side of the normal of that vertex) and the vertex into out[0..6]. */
export function project(c, start, x, y, z, out) {
  const V = c.V, N = c.N, F = c.F, G = c.graph, vf = vertexFaces(c);
  let v = start, bd = (V[3 * v] - x) ** 2 + (V[3 * v + 1] - y) ** 2 + (V[3 * v + 2] - z) ** 2;
  for (;;) {
    let best = -1;
    for (let k = G.indptr[v]; k < G.indptr[v + 1]; k++) { const q = G.indices[k], d = (V[3 * q] - x) ** 2 + (V[3 * q + 1] - y) ** 2 + (V[3 * q + 2] - z) ** 2; if (d < bd) { bd = d; best = q; } }
    if (best < 0) break; v = best;
  }
  let fd = Infinity, bf = -1, bx = V[3 * v], by = V[3 * v + 1], bz = V[3 * v + 2], bv = 0, bw = 0, far = false;
  const around = (w) => {
    for (let k = vf.ptr[w]; k < vf.ptr[w + 1]; k++) {
      const f = vf.idx[k], d = closestOnTriangle(V, F[3 * f], F[3 * f + 1], F[3 * f + 2], x, y, z);
      if (d < fd - 1e-18) { fd = d; bf = f; bx = TP[0]; by = TP[1]; bz = TP[2]; bv = TP[3]; bw = TP[4]; }
    }
  };
  around(v);
  if (bf >= 0) { const wv = F[3 * bf] === v ? 1 - bv - bw : F[3 * bf + 1] === v ? bv : bw; far = wv < 1e-9; }
  if (far || bf < 0) for (let k = G.indptr[v]; k < G.indptr[v + 1]; k++) around(G.indices[k]);
  let nx = N[3 * v], ny = N[3 * v + 1], nz = N[3 * v + 2];
  if (bf >= 0) {                                       // the normal between the corners, each taken on the side of N[v] (a one-sided carrier flips across its seam)
    const w = [1 - bv - bw, bv, bw]; let sx = 0, sy = 0, sz = 0;
    for (let k = 0; k < 3; k++) { const q = F[3 * bf + k], s = (N[3 * q] * nx + N[3 * q + 1] * ny + N[3 * q + 2] * nz) < 0 ? -w[k] : w[k]; sx += s * N[3 * q]; sy += s * N[3 * q + 1]; sz += s * N[3 * q + 2]; }
    const l = Math.hypot(sx, sy, sz); if (l > 1e-12) { nx = sx / l; ny = sy / l; nz = sz / l; }
  }
  out[0] = bx; out[1] = by; out[2] = bz; out[3] = nx; out[4] = ny; out[5] = nz; out[6] = v;
  return out;
}

/** The carrier continued past its boundary, for the sites a patch puts there (latedit.js): the quadric
 *  w = c0 + c1 u + c2 v + c3 u² + c4 u v + c5 v² over the tangent plane at the point s of the carrier (normal n), fitted
 *  to the carrier vertices within `rad` of s — a plane stays a plane, a cylinder or a cone keeps its curvature.
 *  Returns [s, e1, e2, n (3 each), c0 … c5]: what liftExt reads. */
export function fitExt(c, start, sx, sy, sz, nx, ny, nz, rad) {
  const V = c.V, G = c.graph, E = new Float64Array(18); let ex, ey, ez;
  if (Math.abs(nx) < 0.9) { ex = 1 - nx * nx; ey = -nx * ny; ez = -nx * nz; } else { ex = -ny * nx; ey = 1 - ny * ny; ez = -ny * nz; }
  const el = Math.hypot(ex, ey, ez) || 1; ex /= el; ey /= el; ez /= el;
  const fx = ny * ez - nz * ey, fy = nz * ex - nx * ez, fz = nx * ey - ny * ex;
  E[0] = sx; E[1] = sy; E[2] = sz; E[3] = ex; E[4] = ey; E[5] = ez; E[6] = fx; E[7] = fy; E[8] = fz; E[9] = nx; E[10] = ny; E[11] = nz;
  const A = Array.from({ length: 6 }, () => new Float64Array(7)), row = new Float64Array(6), seen = new Set([start]), r2 = rad * rad; let front = [start], m = 0;
  while (front.length) {
    const nf = [];
    for (const q of front) {
      const dx = V[3 * q] - sx, dy = V[3 * q + 1] - sy, dz = V[3 * q + 2] - sz, u = (dx * ex + dy * ey + dz * ez) / rad, v = (dx * fx + dy * fy + dz * fz) / rad, w = dx * nx + dy * ny + dz * nz;
      row[0] = 1; row[1] = u; row[2] = v; row[3] = u * u; row[4] = u * v; row[5] = v * v; m++;
      for (let r = 0; r < 6; r++) { for (let k = 0; k < 6; k++) A[r][k] += row[r] * row[k]; A[r][6] += row[r] * w; }
      for (let k = G.indptr[q]; k < G.indptr[q + 1]; k++) { const t = G.indices[k]; if (!seen.has(t) && (V[3 * t] - sx) ** 2 + (V[3 * t + 1] - sy) ** 2 + (V[3 * t + 2] - sz) ** 2 <= r2) { seen.add(t); nf.push(t); } }
    }
    front = nf;
  }
  if (m < 12) return E;                                  // too few vertices: the tangent plane
  for (let r = 0; r < 6; r++) {                          // Gauss–Jordan with partial pivoting
    let piv = r; for (let k = r + 1; k < 6; k++) if (Math.abs(A[k][r]) > Math.abs(A[piv][r])) piv = k;
    const t = A[r]; A[r] = A[piv]; A[piv] = t; const d = A[r][r]; if (Math.abs(d) < 1e-12) return E;
    for (let k = r; k < 7; k++) A[r][k] /= d;
    for (let q = 0; q < 6; q++) if (q !== r) { const f = A[q][r]; if (f) for (let k = r; k < 7; k++) A[q][k] -= f * A[r][k]; }
  }
  E[12] = A[0][6]; E[13] = A[1][6] / rad; E[14] = A[2][6] / rad; E[15] = A[3][6] / r2; E[16] = A[4][6] / r2; E[17] = A[5][6] / r2;
  return E;
}
/** (x, y, z) put on the continued carrier E (fitExt): writes the point and the normal there into out[0..5] */
export function liftExt(E, x, y, z, out) {
  const dx = x - E[0], dy = y - E[1], dz = z - E[2], u = dx * E[3] + dy * E[4] + dz * E[5], v = dx * E[6] + dy * E[7] + dz * E[8];
  const w = E[12] + E[13] * u + E[14] * v + E[15] * u * u + E[16] * u * v + E[17] * v * v, gu = E[13] + 2 * E[15] * u + E[16] * v, gv = E[14] + E[16] * u + 2 * E[17] * v;
  let nx = E[9] - gu * E[3] - gv * E[6], ny = E[10] - gu * E[4] - gv * E[7], nz = E[11] - gu * E[5] - gv * E[8]; const l = Math.hypot(nx, ny, nz) || 1;
  out[0] = E[0] + u * E[3] + v * E[6] + w * E[9]; out[1] = E[1] + u * E[4] + v * E[7] + w * E[10]; out[2] = E[2] + u * E[5] + v * E[8] + w * E[11];
  out[3] = nx / l; out[4] = ny / l; out[5] = nz / l;
  return out;
}

/** Sub-vertex peak positions: around every peak vertex p, u ≈ c0 + g·x + ½ xᵀH x is fitted by least squares over the
 *  vertices within `depth` rings (about 0.2 a), in the tangent plane at p; the site goes to the maximum, −H⁻¹g, at most
 *  `clamp`·h away, and back onto the carrier. Returns {P, nrm, anchor, shift (in h, per site)}. */
export function subVertexPeaks(c, u, lat, a, { clamp = 1.0 } = {}) {
  const V = c.V, N = c.N, G = c.graph, pk = lat.peaks, n = pk.length, P = new Float64Array(3 * n), nrm = new Float64Array(3 * n), anchor = new Int32Array(n), shift = new Float64Array(n);
  const depth = Math.max(2, Math.round(0.22 * a / c.h)), out = new Float64Array(7), A = Array.from({ length: 6 }, () => new Float64Array(7)), row = new Float64Array(6);
  for (let i = 0; i < n; i++) {
    const p = pk[i], nx = N[3 * p], ny = N[3 * p + 1], nz = N[3 * p + 2];
    let ex, ey, ez; if (Math.abs(nx) < 0.9) { ex = 1 - nx * nx; ey = -nx * ny; ez = -nx * nz; } else { ex = -ny * nx; ey = 1 - ny * ny; ez = -ny * nz; }
    const el = Math.hypot(ex, ey, ez) || 1; ex /= el; ey /= el; ez /= el;
    const fx = ny * ez - nz * ey, fy = nz * ex - nx * ez, fz = nx * ey - ny * ex;
    const ring = new Set([p]); let front = [p];
    for (let d = 0; d < depth; d++) { const nf = []; for (const q of front) for (let k = G.indptr[q]; k < G.indptr[q + 1]; k++) { const w = G.indices[k]; if (!ring.has(w)) { ring.add(w); nf.push(w); } } front = nf; }
    for (const r of A) r.fill(0);
    for (const q of ring) {
      const dx = V[3 * q] - V[3 * p], dy = V[3 * q + 1] - V[3 * p + 1], dz = V[3 * q + 2] - V[3 * p + 2], x = dx * ex + dy * ey + dz * ez, y = dx * fx + dy * fy + dz * fz;
      row[0] = 1; row[1] = x; row[2] = y; row[3] = x * x; row[4] = x * y; row[5] = y * y;
      for (let r = 0; r < 6; r++) { for (let s = 0; s < 6; s++) A[r][s] += row[r] * row[s]; A[r][6] += row[r] * u[q]; }
    }
    let ok = ring.size >= 7, sx = 0, sy = 0;
    for (let r = 0; r < 6 && ok; r++) {                                       // Gauss–Jordan with partial pivoting
      let piv = r; for (let s = r + 1; s < 6; s++) if (Math.abs(A[s][r]) > Math.abs(A[piv][r])) piv = s;
      const t = A[r]; A[r] = A[piv]; A[piv] = t; const d = A[r][r]; if (Math.abs(d) < 1e-14) { ok = false; break; }
      for (let s = r; s < 7; s++) A[r][s] /= d;
      for (let q = 0; q < 6; q++) if (q !== r) { const f = A[q][r]; if (f) for (let s = r; s < 7; s++) A[q][s] -= f * A[r][s]; }
    }
    if (ok) {
      const g1 = A[1][6], g2 = A[2][6], h11 = A[3][6], h12 = A[4][6], h22 = A[5][6], det = 4 * h11 * h22 - h12 * h12;
      if (det > 0 && h11 < 0) { sx = -(2 * h22 * g1 - h12 * g2) / det; sy = -(2 * h11 * g2 - h12 * g1) / det; const l = Math.hypot(sx, sy), lim = clamp * c.h; if (l > lim) { sx *= lim / l; sy *= lim / l; } }
    }
    project(c, p, V[3 * p] + sx * ex + sy * fx, V[3 * p + 1] + sx * ey + sy * fy, V[3 * p + 2] + sx * ez + sy * fz, out);
    P[3 * i] = out[0]; P[3 * i + 1] = out[1]; P[3 * i + 2] = out[2]; nrm[3 * i] = out[3]; nrm[3 * i + 1] = out[4]; nrm[3 * i + 2] = out[5]; anchor[i] = out[6];
    shift[i] = Math.hypot(P[3 * i] - V[3 * p], P[3 * i + 1] - V[3 * p + 1], P[3 * i + 2] - V[3 * p + 2]) / c.h;
  }
  return { P, nrm, anchor, shift };
}

/** every lattice edge shared by two triangles, as [a, c, b, d]: the edge a–c, b opposite in one triangle, d in the other */
export function rhombi(lat) {
  const T = lat.tris, n = lat.peaks.length, m = new Map(), out = [];
  for (let t = 0; t < T.length / 3; t++) for (let k = 0; k < 3; k++) {
    const a = T[3 * t + k], c = T[3 * t + (k + 1) % 3], o = T[3 * t + (k + 2) % 3], key = a < c ? a * n + c : c * n + a, e = m.get(key);
    if (e === undefined) m.set(key, [a, c, o, -1, 1]); else { e[3] = o; e[4]++; }
  }
  for (const e of m.values()) if (e[4] === 2) out.push(e[0], e[1], e[2], e[3]);
  return Int32Array.from(out);
}
/** where each rhombus lies: 0 among regular sites, 1 by a 5 or 7 (an interior site that is not 6-valent), 2 at the rim (a site with an open link) */
export function rhombusClass(lat, Q = rhombi(lat)) {
  const cls = new Uint8Array(Q.length / 4);
  for (let q = 0; q < cls.length; q++) {
    let rim = false, sing = false;
    for (let k = 0; k < 4; k++) { const i = Q[4 * q + k]; if (!lat.interior[i]) rim = true; else if (lat.valence[i] !== 6) sing = true; }
    cls[q] = rim ? 2 : sing ? 1 : 0;
  }
  return cls;
}

const RH = new Float64Array(13);
/** one rhombus a c | b d unfolded about a–c: writes the bend b + d′ − a − c as a 3-vector in the plane of a b c (0..2)
 *  and in the plane of a c d (3..5), its components along and across the edge (6, 7), the edge length (8) and b, d
 *  along / off the edge (9..12); returns false when a triangle is degenerate. */
function unfold(P, a, c, b, d) {
  const ax = P[3 * a], ay = P[3 * a + 1], az = P[3 * a + 2];
  let ex = P[3 * c] - ax, ey = P[3 * c + 1] - ay, ez = P[3 * c + 2] - az; const L = Math.hypot(ex, ey, ez); if (L < 1e-12) return false; ex /= L; ey /= L; ez /= L;
  let bx = P[3 * b] - ax, by = P[3 * b + 1] - ay, bz = P[3 * b + 2] - az; const tb = bx * ex + by * ey + bz * ez; bx -= tb * ex; by -= tb * ey; bz -= tb * ez; const hb = Math.hypot(bx, by, bz); if (hb < 1e-12) return false; bx /= hb; by /= hb; bz /= hb;
  let dx = P[3 * d] - ax, dy = P[3 * d + 1] - ay, dz = P[3 * d + 2] - az; const td = dx * ex + dy * ey + dz * ez; dx -= td * ex; dy -= td * ey; dz -= td * ez; const hd = Math.hypot(dx, dy, dz); if (hd < 1e-12) return false; dx /= hd; dy /= hd; dz /= hd;
  const rx = tb + td - L, ry = hb - hd;                                       // along the edge, and off it towards b
  RH[0] = rx * ex + ry * bx; RH[1] = rx * ey + ry * by; RH[2] = rx * ez + ry * bz;
  RH[3] = rx * ex - ry * dx; RH[4] = rx * ey - ry * dy; RH[5] = rx * ez - ry * dz;
  RH[6] = rx; RH[7] = ry; RH[8] = L; RH[9] = tb; RH[10] = hb; RH[11] = td; RH[12] = hd;
  return true;
}

/** rms, 95th percentile and maximum of a list of angles */
function stats(v) {
  const n = v.length; if (!n) return { rms: 0, p95: 0, max: 0, n: 0 };
  let s2 = 0; for (let i = 0; i < n; i++) s2 += v[i] * v[i];
  const s = Float64Array.from(v).sort();
  return { rms: Math.sqrt(s2 / n), p95: s[Math.min(n - 1, Math.floor(0.95 * n))], max: s[n - 1], n };
}

/** the kink of the two laths crossing at every interior lattice edge, in degrees: {kink (2 per rhombus), rms, p95,
 *  max, n} and, with the classes of the rhombi (rhombusClass), by: the same for the regular, by a 5 / 7 and rim ones */
export function lathBend(lat, Q = rhombi(lat), cls = null) {
  const P = lat.P, nq = Q.length / 4, kink = new Float64Array(2 * nq), all = [], by = cls ? [[], [], []] : null;
  const ang = (x1, y1, x2, y2) => Math.abs(Math.atan2(x1 * y2 - y1 * x2, x1 * x2 + y1 * y2)) * DEG;
  for (let q = 0; q < nq; q++) {
    if (!unfold(P, Q[4 * q], Q[4 * q + 1], Q[4 * q + 2], Q[4 * q + 3])) continue;
    const L = RH[8], tb = RH[9], hb = RH[10], td = RH[11], hd = RH[12];
    const k1 = ang(L - tb, -hb, td, -hd), k2 = ang(-tb, -hb, td - L, -hd);       // c − b → d′ − a   and   a − b → d′ − c
    kink[2 * q] = k1; kink[2 * q + 1] = k2; all.push(k1, k2); if (by) by[cls[q]].push(k1, k2);
  }
  return { kink, ...stats(all), by: by ? by.map(stats) : null };
}

/** the kink of the three lattice lines through every regular site (interior, six neighbours), in degrees: the six
 *  angles around the site, three on either side of a line, differ by twice its kink (0: the line runs straight on) */
export function lineBend(lat) {
  const P = lat.P, T = lat.tris, n = lat.peaks.length, kink = [];
  const ptr = new Int32Array(n + 1); for (let k = 0; k < T.length; k++) ptr[T[k] + 1]++;
  for (let i = 0; i < n; i++) ptr[i + 1] += ptr[i];
  const ref = new Int32Array(T.length), fill = ptr.slice(0, n); for (let k = 0; k < T.length; k++) ref[fill[T[k]]++] = k;
  const A = new Int32Array(6), B = new Int32Array(6), ring = new Int32Array(6), th = new Float64Array(6);
  const corner = (i, x, y) => {
    const ux = P[3 * x] - P[3 * i], uy = P[3 * x + 1] - P[3 * i + 1], uz = P[3 * x + 2] - P[3 * i + 2], vx = P[3 * y] - P[3 * i], vy = P[3 * y + 1] - P[3 * i + 1], vz = P[3 * y + 2] - P[3 * i + 2];
    return Math.atan2(Math.hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx), ux * vx + uy * vy + uz * vz) * DEG;
  };
  for (let i = 0; i < n; i++) {
    if (!lat.interior[i] || lat.valence[i] !== 6 || ptr[i + 1] - ptr[i] !== 6) continue;
    for (let k = 0; k < 6; k++) { const r = ref[ptr[i] + k], t = r - r % 3, s = r % 3; A[k] = T[t + (s + 1) % 3]; B[k] = T[t + (s + 2) % 3]; }
    // the six neighbours in order round the site: each triangle joins two of them
    ring[0] = A[0]; ring[1] = B[0]; let used = 1, ok = true;
    for (let m = 2; m < 6 && ok; m++) {
      const last = ring[m - 1]; let f = -1;
      for (let k = 1; k < 6; k++) if (!(used >> k & 1) && (A[k] === last || B[k] === last)) { f = k; break; }
      if (f < 0) ok = false; else { used |= 1 << f; ring[m] = A[f] === last ? B[f] : A[f]; }
    }
    if (!ok) continue;
    for (let k = 0; k < 6; k++) th[k] = corner(i, ring[k], ring[(k + 1) % 6]);
    for (let k = 0; k < 3; k++) kink.push(Math.abs(th[k] + th[k + 1] + th[k + 2] - th[(k + 3) % 6] - th[(k + 4) % 6] - th[(k + 5) % 6]) / 2);
  }
  return { kink: Float64Array.from(kink), ...stats(kink) };
}

/** A read-out taken out of the field: the same sites, triangles and valences, every site at the peak of u between the
 *  carrier vertices. P are the positions (relax moves them), P0 those as frozen; nrm the surface normal at each site,
 *  anchor the carrier vertex nearest to it; Q the rhombi and cls their classes; shift how far each site left the
 *  vertex of its peak (in h). `frozen: true` tells the view to draw the sites where they are. */
export function freeze(c, u, lat, a) {
  const sv = subVertexPeaks(c, u, lat, a), Q = rhombi(lat);
  return { ...lat, frozen: true, P: Float64Array.from(sv.P), P0: sv.P, nrm: Float64Array.from(sv.nrm), nrm0: sv.nrm, anchor: Int32Array.from(sv.anchor), anchor0: sv.anchor,
    shift: sv.shift, Q, cls: rhombusClass(lat, Q), sweeps: 0 };
}
/** the positions as frozen again */
export function reset(F) { F.P.set(F.P0); F.nrm.set(F.nrm0); F.anchor.set(F.anchor0); F.sweeps = 0; }

/** which sites relax may move (1) — all of them, unless the rim sites (open link) or the 5s and 7s are held */
export function movable(lat, { holdRim = false, holdSingular = false } = {}) {
  const n = lat.peaks.length, m = new Uint8Array(n);
  for (let i = 0; i < n; i++) m[i] = (holdRim && !lat.interior[i]) || (holdSingular && lat.interior[i] && lat.valence[i] !== 6) ? 0 : 1;
  return m;
}

/** site → the rhombi it belongs to (4·rhombus + its place in it: 0 a, 1 c, 2 b, 3 d), its neighbours and the edge
 *  lengths to them as frozen; kept on the lattice */
function tables(F) {
  if (F._tb) return F._tb;
  const Q = F.Q, n = F.peaks.length, T = F.tris, R = F.P0;
  const ptr = new Int32Array(n + 1); for (let k = 0; k < Q.length; k++) ptr[Q[k] + 1]++;
  for (let i = 0; i < n; i++) ptr[i + 1] += ptr[i];
  const ref = new Int32Array(Q.length), fill = ptr.slice(0, n); for (let k = 0; k < Q.length; k++) ref[fill[Q[k]]++] = k;
  const nb = Array.from({ length: n }, () => []);
  for (let t = 0; t < T.length / 3; t++) for (let k = 0; k < 3; k++) { const i = T[3 * t + k], j = T[3 * t + (k + 1) % 3]; if (!nb[i].includes(j)) nb[i].push(j); if (!nb[j].includes(i)) nb[j].push(i); }
  // an edge made by an edit (latedit.js) asks for the mean length around it (F.restFix, keyed i · 2^21 + j with i < j), the others for theirs as frozen
  const fix = F.restFix && F.restFix.size ? F.restFix : null, fixed = (i, j) => (fix ? fix.get(i < j ? i * 2097152 + j : j * 2097152 + i) : undefined);
  const rest = nb.map((l, i) => l.map(j => fixed(i, j) ?? Math.hypot(R[3 * j] - R[3 * i], R[3 * j + 1] - R[3 * i + 1], R[3 * j + 2] - R[3 * i + 2])));
  return (F._tb = { ptr, ref, nb, rest });
}

/**
 * relax(carrier, F, opts): Gauss–Seidel sweeps over the sites of a frozen lattice. Each site goes to where the bends
 * of the rhombi it belongs to cancel (the minimum of Σ |b + d′ − a − c|² in that one site, the others held, every
 * rhombus unfolded about its edge) and where its edges keep the lengths they were frozen with (springs of weight mu),
 * and is put back on the carrier. The connectivity does not change. opts: sweeps; from, the number of sweeps already
 * done (the direction alternates with it, so 20 sweeps are the same lattice in one call or in several); mu (0.3);
 * move (movable(F): Uint8Array, 0 = held). Returns the sweeps done in all.
 */
export function relax(c, F, { sweeps = 20, from = F.sweeps || 0, mu = 0.3, move = null } = {}) {
  const P = F.P, Q = F.Q, n = F.peaks.length, out = new Float64Array(7), { ptr, ref, nb, rest } = tables(F), ext = F.ext && F.ext.size ? F.ext : null;
  for (let s = from; s < from + sweeps; s++) {
    for (let step = 0; step < n; step++) {
      const i = s & 1 ? n - 1 - step : step; if (move && !move[i]) continue;   // alternate the direction: no drift along the numbering
      let dx = 0, dy = 0, dz = 0, w = 0;
      for (let k = ptr[i]; k < ptr[i + 1]; k++) {
        const q = ref[k] >> 2, role = ref[k] & 3;
        if (!unfold(P, Q[4 * q], Q[4 * q + 1], Q[4 * q + 2], Q[4 * q + 3])) continue;
        if (role === 2) { dx -= RH[0]; dy -= RH[1]; dz -= RH[2]; } else if (role === 3) { dx -= RH[3]; dy -= RH[4]; dz -= RH[5]; }
        else { dx += 0.5 * (RH[0] + RH[3]); dy += 0.5 * (RH[1] + RH[4]); dz += 0.5 * (RH[2] + RH[5]); }
        w++;
      }
      if (mu > 0) for (let k = 0; k < nb[i].length; k++) {
        const j = nb[i][k], ex = P[3 * j] - P[3 * i], ey = P[3 * j + 1] - P[3 * i + 1], ez = P[3 * j + 2] - P[3 * i + 2], l = Math.hypot(ex, ey, ez) || 1, f = mu * (1 - rest[i][k] / l);
        dx += f * ex; dy += f * ey; dz += f * ez; w += mu;
      }
      if (!w) continue;
      const E = ext && ext.get(i);                                             // back onto the carrier, or onto its continuation for a site patched beyond the boundary
      if (E) { liftExt(E, P[3 * i] + dx / w, P[3 * i + 1] + dy / w, P[3 * i + 2] + dz / w, out); out[6] = F.anchor[i]; }
      else project(c, F.anchor[i], P[3 * i] + dx / w, P[3 * i + 1] + dy / w, P[3 * i + 2] + dz / w, out);
      P[3 * i] = out[0]; P[3 * i + 1] = out[1]; P[3 * i + 2] = out[2]; F.nrm[3 * i] = out[3]; F.nrm[3 * i + 1] = out[4]; F.nrm[3 * i + 2] = out[5]; F.anchor[i] = out[6];
    }
  }
  return (F.sweeps = from + sweeps);
}

/** the shape of a frozen lattice against the one it was frozen with: mean edge length and its spread (cv), the largest
 *  move of a site (model units), and the triangles that turned over (normal reversed; after edits: facing away from the surface normal) */
export function shape(F) {
  const P = F.P, R = F.P0, T = F.tris, n = F.peaks.length, { nb } = tables(F);
  let s1 = 0, s2 = 0, ne = 0, move = 0, folded = 0;
  for (let i = 0; i < n; i++) {
    for (const j of nb[i]) if (j > i) { const l = Math.hypot(P[3 * j] - P[3 * i], P[3 * j + 1] - P[3 * i + 1], P[3 * j + 2] - P[3 * i + 2]); s1 += l; s2 += l * l; ne++; }
    const d = Math.hypot(P[3 * i] - R[3 * i], P[3 * i + 1] - R[3 * i + 1], P[3 * i + 2] - R[3 * i + 2]); if (d > move) move = d;
  }
  const nrmOf = (X, t, o) => {
    const i = T[3 * t], j = T[3 * t + 1], k = T[3 * t + 2], ux = X[3 * j] - X[3 * i], uy = X[3 * j + 1] - X[3 * i + 1], uz = X[3 * j + 2] - X[3 * i + 2], vx = X[3 * k] - X[3 * i], vy = X[3 * k + 1] - X[3 * i + 1], vz = X[3 * k + 2] - X[3 * i + 2];
    o[0] = uy * vz - uz * vy; o[1] = uz * vx - ux * vz; o[2] = ux * vy - uy * vx;
  };
  const n0 = [0, 0, 0], n1 = [0, 0, 0], N = F.nrm;
  for (let t = 0; t < T.length / 3; t++) {
    nrmOf(P, t, n1);
    if (F.edits) { for (let d = 0; d < 3; d++) n0[d] = N[3 * T[3 * t] + d] + N[3 * T[3 * t + 1] + d] + N[3 * T[3 * t + 2] + d]; }   // an edited lattice has triangles that were not there when it was frozen: against the surface normal
    else nrmOf(R, t, n0);
    if (n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2] < 0) folded++;
  }
  const mean = ne ? s1 / ne : 0;
  return { edge: mean, cv: ne && mean ? Math.sqrt(Math.max(s2 / ne - mean * mean, 0)) / mean : 0, move, folded, edges: ne };
}
