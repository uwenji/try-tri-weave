// Mesh file readers for the browser import (STL binary / ASCII, OBJ, PLY ascii / binary, OFF) and the .lshc writer
// (same layout as export_carrier.py write_lshc, read back by carrier.js parseLSHC).
// Readers return { P, F, faces }: P positions (Float64Array 3n), F triangle indices (Int32Array 3m) or null for a
// triangle soup (STL: every 3 consecutive points are one triangle), polygons fanned into triangles.

export const MESH_EXT = ['stl', 'obj', 'ply', 'off'];

export function readMesh(buf, name = '') {
  const ext = (name.split('.').pop() || '').toLowerCase();
  const u8 = new Uint8Array(buf);
  const head = latin1(u8.subarray(0, Math.min(u8.length, 512)));
  if (ext === 'ply' || head.startsWith('ply')) return readPLY(buf);
  if (ext === 'off' || /^\s*[CN]?OFF\b/.test(head)) return readOFF(buf);
  if (ext === 'obj') return readOBJ(buf);
  if (ext === 'stl' || isBinarySTL(buf)) return readSTL(buf);
  if (/^\s*(v|vn|vt|f|o|g|#|mtllib)\s/m.test(head)) return readOBJ(buf);
  throw new Error(`unknown mesh format: ${name || 'buffer'} (use .stl, .obj, .ply or .off)`);
}

function latin1(u8) { let s = ''; for (let i = 0; i < u8.length; i++) s += String.fromCharCode(u8[i]); return s; }
function isBinarySTL(buf) {
  if (buf.byteLength < 84) return false;
  const n = new DataView(buf).getUint32(80, true);
  return n > 0 && buf.byteLength === 84 + n * 50;
}

// ------------------------------------------------------------------ STL
export function readSTL(buf) {
  if (isBinarySTL(buf)) {
    const dv = new DataView(buf), n = dv.getUint32(80, true), P = new Float64Array(n * 9);
    for (let i = 0, o = 84; i < n; i++, o += 50) for (let k = 0; k < 9; k++) P[i * 9 + k] = dv.getFloat32(o + 12 + 4 * k, true);
    return { P, F: null, faces: n };
  }
  // ASCII: "vertex x y z" lines (anything else ignored)
  const txt = new TextDecoder().decode(new Uint8Array(buf)); const pts = [];
  const re = /vertex\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)/g; let m;
  while ((m = re.exec(txt))) pts.push(+m[1], +m[2], +m[3]);
  const nt = Math.floor(pts.length / 9);
  if (!nt) throw new Error('STL: no triangles found');
  return { P: Float64Array.from(pts.slice(0, nt * 9)), F: null, faces: nt };
}

// ------------------------------------------------------------------ OBJ
export function readOBJ(buf) {
  const txt = new TextDecoder().decode(new Uint8Array(buf));
  const P = [], F = []; let nv = 0;
  let i = 0; const L = txt.length;
  while (i < L) {
    let j = txt.indexOf('\n', i); if (j < 0) j = L;
    const c0 = txt.charCodeAt(i), c1 = txt.charCodeAt(i + 1);
    if (c0 === 118 && (c1 === 32 || c1 === 9)) {                         // 'v '
      const s = txt.slice(i + 2, j).trim().split(/\s+/);
      P.push(+s[0], +s[1], +s[2]); nv++;
    } else if (c0 === 102 && (c1 === 32 || c1 === 9)) {                  // 'f '
      const s = txt.slice(i + 2, j).trim().split(/\s+/); const idx = [];
      for (const t of s) {
        if (!t) continue; let v = parseInt(t, 10); if (Number.isNaN(v)) continue;
        idx.push(v < 0 ? nv + v : v - 1);
      }
      for (let k = 1; k + 1 < idx.length; k++) F.push(idx[0], idx[k], idx[k + 1]);
    }
    i = j + 1;
  }
  if (!F.length) throw new Error('OBJ: no faces found');
  return { P: Float64Array.from(P), F: Int32Array.from(F), faces: F.length / 3 };
}

// ------------------------------------------------------------------ OFF
export function readOFF(buf) {
  const txt = new TextDecoder().decode(new Uint8Array(buf));
  const tok = txt.replace(/#[^\n]*/g, ' ').trim().split(/\s+/); let p = 0;
  const h = tok[p++]; if (!/OFF$/.test(h)) throw new Error('OFF: bad header');
  const nv = +tok[p++], nf = +tok[p++]; p++;                              // edge count ignored
  const P = new Float64Array(3 * nv), F = [];
  // COFF / NOFF vertex lines carry a colour (4 values) or a normal (3) after x y z
  const extra = h === 'COFF' ? 4 : h === 'NOFF' ? 3 : h === 'CNOFF' ? 7 : 0;
  for (let i = 0; i < nv; i++) { P[3 * i] = +tok[p++]; P[3 * i + 1] = +tok[p++]; P[3 * i + 2] = +tok[p++]; p += extra; }
  for (let f = 0; f < nf; f++) {                                          // per-face colours are not supported
    const k = +tok[p++]; const idx = []; for (let s = 0; s < k; s++) idx.push(+tok[p++]);
    for (let s = 1; s + 1 < k; s++) F.push(idx[0], idx[s], idx[s + 1]);
  }
  return { P, F: Int32Array.from(F), faces: F.length / 3 };
}

// ------------------------------------------------------------------ PLY
const PLY_T = {
  char: ['getInt8', 1], int8: ['getInt8', 1], uchar: ['getUint8', 1], uint8: ['getUint8', 1],
  short: ['getInt16', 2], int16: ['getInt16', 2], ushort: ['getUint16', 2], uint16: ['getUint16', 2],
  int: ['getInt32', 4], int32: ['getInt32', 4], uint: ['getUint32', 4], uint32: ['getUint32', 4],
  float: ['getFloat32', 4], float32: ['getFloat32', 4], double: ['getFloat64', 8], float64: ['getFloat64', 8],
};

export function readPLY(buf) {
  const u8 = new Uint8Array(buf);
  // header ends at "end_header\n"
  let hend = -1; const probe = latin1(u8.subarray(0, Math.min(u8.length, 65536)));
  const m = /end_header\r?\n/.exec(probe); if (!m) throw new Error('PLY: no end_header');
  hend = m.index + m[0].length;
  const lines = probe.slice(0, m.index).split(/\r?\n/);
  let format = 'ascii'; const elems = [];
  for (const line of lines) {
    const s = line.trim().split(/\s+/);
    if (s[0] === 'format') format = s[1];
    else if (s[0] === 'element') elems.push({ name: s[1], count: +s[2], props: [] });
    else if (s[0] === 'property' && elems.length) {
      const e = elems[elems.length - 1];
      if (s[1] === 'list') e.props.push({ list: true, ct: s[2], it: s[3], name: s[4] });
      else e.props.push({ list: false, t: s[1], name: s[2] });
    }
  }
  const ev = elems.find(e => e.name === 'vertex'), ef = elems.find(e => e.name === 'face');
  if (!ev) throw new Error('PLY: no vertex element');
  const P = new Float64Array(3 * ev.count), F = [];
  const xi = ev.props.findIndex(p => p.name === 'x'), yi = ev.props.findIndex(p => p.name === 'y'), zi = ev.props.findIndex(p => p.name === 'z');
  const fi = ef ? ef.props.findIndex(p => p.list && (p.name === 'vertex_indices' || p.name === 'vertex_index')) : -1;
  if (format === 'ascii') {
    const tok = new TextDecoder().decode(u8.subarray(hend)).trim().split(/\s+/); let p = 0;
    for (const e of elems) for (let r = 0; r < e.count; r++) {
      for (let q = 0; q < e.props.length; q++) {
        const pr = e.props[q];
        if (pr.list) {
          const k = +tok[p++]; const idx = []; for (let s = 0; s < k; s++) idx.push(+tok[p++]);
          if (e === ef && q === fi) for (let s = 1; s + 1 < k; s++) F.push(idx[0], idx[s], idx[s + 1]);
        } else {
          const v = +tok[p++];
          if (e === ev) { if (q === xi) P[3 * r] = v; else if (q === yi) P[3 * r + 1] = v; else if (q === zi) P[3 * r + 2] = v; }
        }
      }
    }
  } else {
    const le = format === 'binary_little_endian'; if (!le && format !== 'binary_big_endian') throw new Error('PLY: format ' + format);
    const dv = new DataView(buf); let o = hend;
    const rd = (t) => { const T = PLY_T[t]; if (!T) throw new Error('PLY: type ' + t); const v = dv[T[0]](o, le); o += T[1]; return v; };
    for (const e of elems) {
      // fast path for all-scalar elements that are not the vertex element: skip by stride
      if (e !== ev && e !== ef && e.props.every(p => !p.list)) { o += e.count * e.props.reduce((s, p) => s + PLY_T[p.t][1], 0); continue; }
      for (let r = 0; r < e.count; r++) for (let q = 0; q < e.props.length; q++) {
        const pr = e.props[q];
        if (pr.list) {
          const k = rd(pr.ct);
          if (e === ef && q === fi) { const i0 = rd(pr.it); let prev = rd(pr.it); for (let s = 2; s < k; s++) { const cur = rd(pr.it); F.push(i0, prev, cur); prev = cur; } }
          else for (let s = 0; s < k; s++) rd(pr.it);
        } else {
          const v = rd(pr.t);
          if (e === ev) { if (q === xi) P[3 * r] = v; else if (q === yi) P[3 * r + 1] = v; else if (q === zi) P[3 * r + 2] = v; }
        }
      }
    }
  }
  if (!F.length) throw new Error('PLY: no faces (point clouds are not supported)');
  return { P, F: Int32Array.from(F), faces: F.length / 3 };
}

// ------------------------------------------------------------------ .lshc writer
/** data: {name, meta, V, F, Fi, Li, K, B} → ArrayBuffer in the export_carrier.py layout */
export function writeLSHC(data) {
  const arrays = [
    ['V', Float32Array.from(data.V), 'float32', [data.V.length / 3, 3]],
    ['F', Uint32Array.from(data.F), 'uint32', [data.F.length / 3, 3]],
    ['Fi', Uint32Array.from(data.Fi), 'uint32', [data.Fi.length / 3, 3]],
    ['Li', Float32Array.from(data.Li), 'float32', [data.Li.length / 3, 3]],
    ['K', Float32Array.from(data.K), 'float32', [data.K.length]],
    ['B', Uint8Array.from(data.B), 'uint8', [data.B.length]],
  ];
  const header = { name: data.name, meta: data.meta, arrays: {} }; let offset = 0;
  for (const [k, a, dtype, shape] of arrays) {
    header.arrays[k] = { dtype, shape, offset, bytes: a.byteLength };
    offset += a.byteLength + ((4 - (a.byteLength % 4)) % 4);
  }
  let hj = new TextEncoder().encode(JSON.stringify(header));
  const hpad = (4 - (hj.length % 4)) % 4;
  const out = new Uint8Array(8 + hj.length + hpad + offset);
  out.set([76, 83, 72, 67]); new DataView(out.buffer).setUint32(4, hj.length + hpad, true);
  out.set(hj, 8); for (let i = 0; i < hpad; i++) out[8 + hj.length + i] = 32;
  const base = 8 + hj.length + hpad;
  for (const [k, a] of arrays) out.set(new Uint8Array(a.buffer, a.byteOffset, a.byteLength), base + header.arrays[k].offset);
  return out.buffer;
}
