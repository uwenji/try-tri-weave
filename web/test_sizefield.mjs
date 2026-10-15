// node test_sizefield.mjs [carrier] [steps] [valence] — uniform a vs conformal size field, CPU solver
import fs from 'node:fs';
import path from 'node:path';
import { Carrier, parseLSHC } from './js/carrier.js';
import { CpuSolver } from './js/solver_cpu.js';
import { extractLattice, chargeSummary, traceStrands, edgeStrain, orientationField } from './js/lattice.js';

const name = process.argv[2] || 'v4Cone_round', steps = +(process.argv[3] || 4000), valArg = process.argv[4] ? +process.argv[4] : null;
const manifest = JSON.parse(fs.readFileSync(path.join(import.meta.dirname, 'manifest.json'), 'utf8'));
const entry = manifest.find(m => m.name === name);
const buf = fs.readFileSync(path.join(import.meta.dirname, entry.file));
const c = new Carrier(parseLSHC(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength)));
const a = process.argv[5] ? +process.argv[5] : entry.a_default;
const det = c.detectSingularity(a); const valence = valArg ?? Math.max(3, Math.min(9, det.valence));
console.log(`${c.name}: ${c.n} v, Ω ${(c.omega * 180 / Math.PI).toFixed(1)}°, h ${c.h.toFixed(3)}, a = ${a.toFixed(3)} (${(a / c.h).toFixed(1)} edges per a); seed v${valence} at ${det.centre}`);
for (const conformal of [false, true]) {
  let ops = c, phi = null, t0 = performance.now();
  if (conformal) {
    const sf = c.sizeField([{ vertex: det.centre, valence }], a, det.centre); phi = sf.phi;
    ops = c.conformalOps(phi, a);
    console.log(`size field: φ ∈ [${sf.min.toFixed(2)}, ${sf.max.toFixed(2)}] (spacing a·e^-φ from ${(a * Math.exp(-sf.max)).toFixed(2)} to ${(a * Math.exp(-sf.min)).toFixed(2)}), residual ${sf.residualDeg.toFixed(1)}° (${sf.open ? 'rim' : 'spread'}), PCG ${sf.iters} it, ${ops.underResolvedPct.toFixed(1)}% edges > a/6 in ĝ, bad tris ${ops.badTriangles}, area ratio ${(ops.area / c.area).toFixed(2)}  [${(performance.now() - t0).toFixed(0)} ms]`);
  }
  const s = new CpuSolver(ops, { a, cgIters: 12 });
  s.setInitial(c.seedField(det.centre, valence, a, { phi }));
  t0 = performance.now();
  let lat, cs;
  for (let k = 0; k < steps; k += 1000) {
    await s.step(Math.min(1000, steps - k));
    lat = extractLattice(c, s.u, a); cs = chargeSummary(lat);
    console.log(`  ${conformal ? 'conformal' : 'uniform  '} step ${s.steps}: peaks ${lat ? lat.peaks.length : 0} interior ${cs.nInterior} ${JSON.stringify(cs.hist)} Σ(6-n) ${cs.charge}  [${((performance.now() - t0) / 1000).toFixed(1)} s]`);
  }
  if (lat && lat.tris.length) {
    const st = traceStrands(lat.peaks.length, lat.tris), es = edgeStrain(lat, st, a, phi), of = orientationField(lat, c);
    const byFam = [0, 1, 2].map(f => { let s = 0, s2 = 0, k = 0; for (let e = 0; e < es.fam.length; e++) if (es.fam[e] === f) { s += es.strain[e]; s2 += es.strain[e] ** 2; k++; } return k ? `${(100 * s / k).toFixed(1)}±${(100 * Math.sqrt(Math.max(s2 / k - (s / k) ** 2, 0))).toFixed(1)}% (${k})` : '—'; });
    const unk = Array.from(es.fam).filter(f => f < 0).length;
    console.log(`  strain per family A/B/C: ${byFam.join(' · ')} · unassigned edges ${unk} · orientation coherence ${of.coherence.toFixed(3)} over ${of.count} sites`);
  }
}
