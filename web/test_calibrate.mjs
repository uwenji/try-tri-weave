// node test_calibrate.mjs — equilibrium strain of the continuum SH23 ground state (ComFiT-style conf_strain_to_equilibrium):
// relaxes the hexagonal pattern on a periodic cell strained by ε, minimises the free-energy density over ε.
// Result: ε* ≈ +0.3 % at (r, b) = (−0.08, 1.2), +0.1 % at (−0.0267, 0.693) — the 2–3 % edge shift seen on the carriers is the
// mesh dispersion of the cotan Laplacian (corrected in the solvers, see test_dispersion.mjs), not a property of the equation.
import { equilibriumStrain, relaxedEnergy } from './calibrate.js';
for (const [r, b] of [[-0.08, 1.2], [-0.0267, 0.693], [-0.05, 0.9]]) {
  const t0 = Date.now(); const { eps, F, samples } = equilibriumStrain(r, b);
  console.log(`r = ${r}, b = ${b}: ε* = ${(100 * eps).toFixed(2)} %  (q*/q0 = ${(1 / (1 + eps)).toFixed(4)}),  F = ${F.toExponential(4)},  ${samples.length} evaluations, ${Date.now() - t0} ms`);
}
const row = []; for (let e = -0.04; e <= 0.0401; e += 0.01) row.push(`${(100 * e).toFixed(0)}%: ${relaxedEnergy(-0.08, 1.2, e).F.toExponential(3)}`);
console.log('F(ε) at (−0.08, 1.2):', row.join('  '));
