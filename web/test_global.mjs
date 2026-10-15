// node test_global.mjs [carrier] [steps] [level] [a] — no-seed pipeline: smoothest 6-RoSy → stripe phases → u0 → SH23 → lattice
import fs from 'node:fs'; import path from 'node:path';
import { Carrier, parseLSHC } from './js/carrier.js';
import { CpuSolver, SQ3 } from './js/solver_cpu.js';
import { extractLattice, chargeSummary } from './js/lattice.js';
import { smoothestRoSy } from './js/rosy.js';
import { stripes } from './js/stripes.js';
const name = process.argv[2] || 'treeBig', steps = +(process.argv[3] || 2000), level = +(process.argv[4] || 0);
const manifest = JSON.parse(fs.readFileSync(path.join(import.meta.dirname, 'manifest.json'), 'utf8'));
const entry = manifest.find(m => m.name === name) || manifest.find(m => m.file.includes(name));
const buf = fs.readFileSync(path.join(import.meta.dirname, entry.file));
let c = new Carrier(parseLSHC(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength))); if (level) c = c.refined(level);
const a = process.argv[5] ? +process.argv[5] : Math.max(8 * c.h, entry.a_default), q0 = 4 * Math.PI / (SQ3 * a);
console.log(`${c.name}: ${c.n} v (×${2 ** level}), h ${c.h.toFixed(3)}, a = ${a.toFixed(3)} (${(a / c.h).toFixed(1)} edges per a, ${(c.area / (SQ3 / 2 * a * a)).toFixed(0)} sites), Ω/60° = ${(c.omega * 180 / Math.PI / 60).toFixed(1)}`);
let t0 = performance.now();
const r = smoothestRoSy(c); const hist = {}; for (const s of r.singular) hist[s.valence] = (hist[s.valence] || 0) + 1;
console.log(`rosy: ${r.singular.length} singular faces ${JSON.stringify(hist)}, Σ index ${r.sumIndex}/6  [${((performance.now() - t0) / 1000).toFixed(1)} s]`);
t0 = performance.now();
const st = stripes(c, r, q0, { log: (s) => process.stdout.write('  ' + s + '\n') });
const bh = {}; for (const b of st.branch) bh[b.winding] = (bh[b.winding] || 0) + 1;
console.log(`stripes: ${st.branch.length} branch points ${JSON.stringify(bh)}, worst branch match ${(st.worstMatch * 180 / Math.PI).toFixed(1)}°, u0 range ${Math.min(...st.u0).toFixed(2)}…${Math.max(...st.u0).toFixed(2)}  [${((performance.now() - t0) / 1000).toFixed(1)} s]`);
const s = new CpuSolver(c, { a, cgIters: 12 }); s.setInitial(st.u0);
t0 = performance.now();
for (let k = 0; k < steps; k += 500) {
  await s.step(Math.min(500, steps - k));
  const lat = extractLattice(c, s.u, a); const cs = chargeSummary(lat);
  console.log(`step ${s.steps}: umax ${Math.max(...s.u).toFixed(3)} peaks ${lat ? lat.peaks.length : 0} interior ${cs.nInterior} hist ${JSON.stringify(cs.hist)} Σ(6−n) ${cs.charge}  [${((performance.now() - t0) / 1000).toFixed(1)} s]`);
}
