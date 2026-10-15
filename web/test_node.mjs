// node test_node.mjs [carrier-name] [steps] [valence]   — exercises carrier → detect → seed → CPU solver → lattice → strands
import fs from 'node:fs';
import path from 'node:path';
import { Carrier, parseLSHC } from './js/carrier.js';
import { CpuSolver } from './js/solver_cpu.js';
import { extractLattice, chargeSummary, traceStrands } from './js/lattice.js';

const name = process.argv[2] || 'v5Cone_sharp', steps = +(process.argv[3] || 3000), valArg = process.argv[4] ? +process.argv[4] : null;
const manifest = JSON.parse(fs.readFileSync(path.join(import.meta.dirname, 'manifest.json'), 'utf8'));
const entry = manifest.find(m => m.name === name) || manifest.find(m => m.file.includes(name));
const buf = fs.readFileSync(path.join(import.meta.dirname, entry.file));
let t0 = performance.now();
const c = new Carrier(parseLSHC(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength)));
console.log(`${c.name}: ${c.n} v, ${c.F.length / 3} f, chi ${c.chi}, Ω ${(c.omega * 180 / Math.PI).toFixed(1)}°, area ${c.area.toFixed(1)}, h ${c.h.toFixed(3)}  [${(performance.now() - t0).toFixed(0)} ms]`);
const a = process.argv[5] ? +process.argv[5] : entry.a_default;
t0 = performance.now();
const det = c.detectSingularity(a);
console.log(`detect: valence ${det.valence} at v${det.centre} (${[0, 1, 2].map(k => c.V[3 * det.centre + k].toFixed(2))}) ${det.info.mode}, ${det.info.withinA.toFixed(0)}% within a, defect ${det.info.defectAtCentre.toFixed(1)}°  [${(performance.now() - t0).toFixed(0)} ms]`);
const valence = valArg ?? Math.max(3, Math.min(9, det.valence));
t0 = performance.now();
const u0 = c.seedField(det.centre, valence, a);
let umax = 0; for (const v of u0) umax = Math.max(umax, v);
console.log(`seed v${valence}, a = ${a.toFixed(3)}: max ${umax.toFixed(3)}  [${(performance.now() - t0).toFixed(0)} ms]`);
const s = new CpuSolver(c, { a, cgIters: 12 });
s.setInitial(u0);
t0 = performance.now();
const every = 500;
for (let k = 0; k < steps; k += every) {
  await s.step(Math.min(every, steps - k));
  const lat = extractLattice(c, s.u, a); const cs = chargeSummary(lat);
  console.log(`step ${s.steps}: umax ${Math.max(...s.u).toFixed(3)} relres ${s.lastRelRes.toExponential(1)} peaks ${lat ? lat.peaks.length : 0} interior ${cs.nInterior} hist ${JSON.stringify(cs.hist)} Σ(6-n) ${cs.charge}  [${((performance.now() - t0) / 1000).toFixed(1)} s]`);
}
const lat = extractLattice(c, s.u, a);
const st = traceStrands(lat.peaks.length, lat.tris);
const lens = st.strands.map(x => x.knots.length);
console.log(`strands: ${st.strands.length} (closed ${st.strands.filter(x => x.closed).length}), knots per strand ${Math.min(...lens)}–${Math.max(...lens)}, families ${JSON.stringify([0, 1, 2].map(f => st.family.filter(x => x === f).length))}, colour conflicts ${st.colourConflicts}, parity conflicts ${st.parityConflicts}, lattice edges ${st.edges.length}, tris ${lat.tris.length / 3}`);
console.log(`steps/s: ${(steps / ((performance.now() - t0) / 1000)).toFixed(0)}`);
