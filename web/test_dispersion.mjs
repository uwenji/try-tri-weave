// node test_dispersion.mjs [carrier] [edges-per-a] [steps] — realised lattice edge / a with and without the dispersion correction
import fs from 'node:fs'; import path from 'node:path';
import { Carrier, parseLSHC } from './js/carrier.js';
import { CpuSolver } from './js/solver_cpu.js';
import { extractLattice } from './js/lattice.js';
const name = process.argv[2] || 'Tube', epa = +(process.argv[3] || 8), steps = +(process.argv[4] || 4000);
const manifest = JSON.parse(fs.readFileSync(path.join(import.meta.dirname, 'manifest.json'), 'utf8'));
const entry = manifest.find(m => m.name === name) || manifest.find(m => m.file.includes(name));
const buf = fs.readFileSync(path.join(import.meta.dirname, entry.file));
const c = new Carrier(parseLSHC(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength)));
const a = epa * c.h; const det = c.detectSingularity(a); const valence = Math.max(3, Math.min(9, det.valence));
console.log(`${c.name}: ${c.n} v, h ${c.h.toFixed(4)}, a = ${a.toFixed(4)} (${epa} edges per a), seed v${valence} at ${det.centre}`);
for (const dispersion of [false, true]) {
  const s = new CpuSolver(c, { a, cgIters: 12, dispersion }); s.setInitial(c.seedField(det.centre, valence, a));
  const t0 = performance.now(); await s.step(steps);
  const lat = extractLattice(c, s.u, a);
  // lattice edge lengths from the dual triangles (peak-to-peak distances), interior only
  const P = lat.peaks, T = lat.tris, V = c.V, L = [];
  const pos = (k) => [V[3 * P[k]], V[3 * P[k] + 1], V[3 * P[k] + 2]];
  for (let t = 0; t < T.length; t += 3) for (let e = 0; e < 3; e++) { const p = pos(T[t + e]), q = pos(T[t + (e + 1) % 3]); L.push(Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2])); }
  L.sort((x, y) => x - y); const med = L[L.length >> 1], mean = L.reduce((x, y) => x + y, 0) / L.length;
  console.log(`dispersion ${dispersion ? 'on ' : 'off'}: iq2 factor ${(1 / s.disp).toFixed(4)}, ${P.length} peaks, median edge = ${(med / a).toFixed(4)} a, mean = ${(mean / a).toFixed(4)} a, hist ${JSON.stringify(Object.fromEntries(Object.entries(lat.hist || {})))}  [${((performance.now() - t0) / 1000).toFixed(1)} s]`);
}
