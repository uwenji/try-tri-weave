import { traceStrands } from './js/lattice.js';
// regular triangular lattice patch: points (i + j/2, j*sqrt3/2) for i,j in a rhombus, triangles (i,j),(i+1,j),(i,j+1) and (i+1,j),(i+1,j+1),(i,j+1)
const N = 7; const id = (i, j) => i * (N + 1) + j; const tris = [];
for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) { tris.push(id(i, j), id(i + 1, j), id(i, j + 1)); tris.push(id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)); }
const st = traceStrands((N + 1) * (N + 1), Int32Array.from(tris));
const lens = st.strands.map(s => s.knots.length);
console.log('strands', st.strands.length, 'lengths', lens.join(','), 'families', [0,1,2].map(f => st.family.filter(x => x === f).length), 'colour conflicts', st.colourConflicts, 'parity conflicts', st.parityConflicts);
// check straightness: in a regular lattice each strand's knots should be collinear
const P = []; for (let i = 0; i <= N; i++) for (let j = 0; j <= N; j++) P.push([i + j / 2, j * Math.sqrt(3) / 2]);
let maxdev = 0;
for (const s of st.strands) { const pts = s.knots.map(e => { const [a, b] = st.edges[e]; return [(P[a][0] + P[b][0]) / 2, (P[a][1] + P[b][1]) / 2]; });
  if (pts.length < 3) continue; const [x0, y0] = pts[0], [x1, y1] = pts[pts.length - 1]; const L = Math.hypot(x1 - x0, y1 - y0);
  for (const [x, y] of pts) maxdev = Math.max(maxdev, Math.abs((x1 - x0) * (y - y0) - (y1 - y0) * (x - x0)) / L); }
console.log('max deviation of knots from the strand chord (should be ~0):', maxdev.toFixed(4));
