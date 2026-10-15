// node test_predict.mjs <carrier> <steps> [r] [b] [dt] [anneal r2,b2,step] — tracked read-outs vs the growth prediction
import fs from 'node:fs';
import path from 'node:path';
import { Carrier, parseLSHC } from './js/carrier.js';
import { CpuSolver } from './js/solver_cpu.js';
import { extractLattice, chargeSummary } from './js/lattice.js';
import { predictGrowth, frontState } from './js/frontspeed.js';

const [name = 'v4Cone_sharp', stepsS = '3000', rS = '-0.08', bS = '1.2', dtS = '0.1', annS = ''] = process.argv.slice(2);
const steps = +stepsS, r = +rS, b = +bS, dt = +dtS;
const manifest = JSON.parse(fs.readFileSync(path.join(import.meta.dirname, 'manifest.json'), 'utf8'));
const entry = manifest.find(m => m.name === name);
const buf = fs.readFileSync(path.join(import.meta.dirname, entry.file));
const c = new Carrier(parseLSHC(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength)));
const a = +(process.env.A || entry.a_default), det = c.detectSingularity(a), valence = Math.max(3, Math.min(9, +(process.env.VAL || det.valence)));
const seedR = +(process.env.SEEDR || 2);
let anneal = { on: false };
if (annS) { const [r2, b2, st] = annS.split(',').map(Number); anneal = { on: true, r2, b2, when: 'step', step: st, startedAt: -1 }; }
const sol = new CpuSolver(c, { a, r, b, dt, cgIters: 12 });
sol.setInitial(c.seedField(det.centre, valence, a, { radius: seedR }));
const pred = predictGrowth({ dist: c.geodesicDistance(det.centre), m: c.m, K: c.K, boundary: c.boundary, rim: c.rimDistance() },
  { a, seedR, seedKind: 'wave', valence, r, b, dt, anneal, xmax: steps, every: 100 });
const fs0 = frontState(r, b, dt);
console.log(`${name}: a ${a.toFixed(3)} (${(a / c.h).toFixed(1)} edges/a), v${valence}, r ${r} b ${b} dt ${dt}: ${fs0.state}, ${fs0.vA100.toFixed(3)} a/100 steps, Rc ${fs0.RcA.toFixed(2)} a; Rmax ${(pred.Rmax / a).toFixed(2)} a, ${pred.total.toFixed(0)} sites; covered at ${pred.coveredAt}${anneal.on ? `, anneal at ${pred.switchAt}` : ''}`);
let e1 = 0, e2 = 0, k = 0;
for (let s = 0, i = 0; s <= steps; s += 100, i++) {
  if (s) { if (anneal.on && s - 100 < anneal.step && s >= anneal.step) { await sol.step(anneal.step - (s - 100)); sol.setParams({ r: anneal.r2, b: anneal.b2 }); await sol.step(s - anneal.step); } else await sol.step(100); }
  const lat = extractLattice(c, sol.u, a), cs = chargeSummary(lat), np = lat ? lat.peaks.length : 0;
  if (np > 5 && np < 0.95 * pred.total) { e1 += (pred.peaks[i] - np) ** 2; e2 += (pred.interior[i] - cs.nInterior) ** 2; k++; }
  if (s % 200 === 0 || s === steps) console.log(`${s}\tpeaks ${np}\tpred ${pred.peaks[i].toFixed(0)}\tinterior ${cs.nInterior}\tpred ${pred.interior[i].toFixed(0)}\tsingular ${cs.singular.length}\tfloor ${pred.floor[i]}\t${JSON.stringify(cs.hist)}`);
}
console.log(`rms(pred − tracked) during growth: peaks ${Math.sqrt(e1 / Math.max(k, 1)).toFixed(1)}, interior ${Math.sqrt(e2 / Math.max(k, 1)).toFixed(1)} over ${k} read-outs`);
