// node test_lift.mjs [carrier] [steps] [level] [a] — solve rotation field + stripes on the BASE carrier, lift to the refined one, run SH
import fs from 'node:fs'; import path from 'node:path';
import { Carrier, parseLSHC } from './js/carrier.js';
import { CpuSolver, SQ3 } from './js/solver_cpu.js';
import { extractLattice, chargeSummary } from './js/lattice.js';
import { smoothestRoSy, liftField } from './js/rosy.js';
import { stripes, liftStripes } from './js/stripes.js';
const name = process.argv[2] || 'treeBig', steps = +(process.argv[3] || 2000), level = +(process.argv[4] || 1);
const manifest = JSON.parse(fs.readFileSync(path.join(import.meta.dirname, 'manifest.json'), 'utf8'));
const entry = manifest.find(m => m.name === name) || manifest.find(m => m.file.includes(name));
const buf = fs.readFileSync(path.join(import.meta.dirname, entry.file));
const base = new Carrier(parseLSHC(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength)));
const c = level ? base.refined(level) : base;
const a = process.argv[5] ? +process.argv[5] : Math.max(8 * c.h, entry.a_default), q0 = 4 * Math.PI / (SQ3 * a);
console.log(`${c.name}: base ${base.n} v → ×${2 ** level} ${c.n} v, a = ${a.toFixed(3)} (${(a / c.h).toFixed(1)} fine / ${(a / base.h).toFixed(1)} base edges per a)`);
let t0 = performance.now();
const r = smoothestRoSy(base); const hist = {}; for (const s of r.singular) hist[s.valence] = (hist[s.valence] || 0) + 1;
console.log(`rosy on base: ${r.singular.length} singular ${JSON.stringify(hist)}  [${((performance.now() - t0) / 1000).toFixed(1)} s]`);
t0 = performance.now(); const st = stripes(base, r, q0);
console.log(`stripes on base: ${st.branch.length} branch points  [${((performance.now() - t0) / 1000).toFixed(1)} s]`);
t0 = performance.now();
const lf = liftField(r.theta, r.frames, base.n, c, 6);
const ls = liftStripes(st.psi, base.n, c, lf.theta, lf.frames, q0, 6, { reversed: lf.reversed });
console.log(`lift to ×${2 ** level}: ${ls.branch.length} branch points on the fine mesh, u0 range ${Math.min(...ls.u0).toFixed(2)}…${Math.max(...ls.u0).toFixed(2)}  [${((performance.now() - t0) / 1000).toFixed(1)} s]`);
const s = new CpuSolver(c, { a, cgIters: 12 }); s.setInitial(ls.u0); t0 = performance.now();
for (let k = 0; k < steps; k += 500) {
  await s.step(Math.min(500, steps - k));
  const lat = extractLattice(c, s.u, a); const cs = chargeSummary(lat);
  console.log(`step ${s.steps}: peaks ${lat ? lat.peaks.length : 0} interior ${cs.nInterior} hist ${JSON.stringify(cs.hist)} Σ(6−n) ${cs.charge}  [${((performance.now() - t0) / 1000).toFixed(1)} s]`);
}
