// node test_rosy.mjs [carrier] — smoothest 6-RoSy on the carrier: singularity count by valence, Σ index vs Ω/60°
import fs from 'node:fs'; import path from 'node:path';
import { Carrier, parseLSHC } from './js/carrier.js';
import { smoothestRoSy } from './js/rosy.js';
const name = process.argv[2] || 'treeBig';
const manifest = JSON.parse(fs.readFileSync(path.join(import.meta.dirname, 'manifest.json'), 'utf8'));
const entry = manifest.find(m => m.name === name) || manifest.find(m => m.file.includes(name));
const buf = fs.readFileSync(path.join(import.meta.dirname, entry.file));
let c = new Carrier(parseLSHC(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength)));
const level = +(process.argv[3] || 0); if (level) c = c.refined(level);
console.log(`${c.name}: ${c.n} v, χ ${c.chi}, Ω ${(c.omega * 180 / Math.PI).toFixed(0)}° → Ω/60° = ${(c.omega * 180 / Math.PI / 60).toFixed(1)}`);
const t0 = performance.now();
const r = smoothestRoSy(c, { log: (s) => process.stdout.write(s + '\n') });
const hist = {}; for (const s of r.singular) hist[s.valence] = (hist[s.valence] || 0) + 1;
console.log(`singular faces ${r.singular.length}, Σ index = ${r.sumIndex}/6, by valence ${JSON.stringify(hist)}, energy ${r.energy.toExponential(3)}, pcg iters ${r.iters}  [${((performance.now() - t0) / 1000).toFixed(1)} s]`);
const zs = c.V.filter((_, k) => k % 3 === 2);
for (const s of r.singular.slice(0, 40)) { const f = s.face, i = c.Fi[3 * f]; console.log(`  v${s.valence} at face ${f}, z = ${c.V[3 * i + 2].toFixed(0)}, |z| = ${r.mag[i].toFixed(3)}`); }
