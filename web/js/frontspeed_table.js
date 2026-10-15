// Front speed of the SH23 hexagon pattern into u = 0, measured on a flat periodic domain with the solver's own
// semi-implicit scheme (sim/front.py in localize/experiments, dt = 0.1, q0 = 1, planar front along a reciprocal vector;
// the other orientation agrees within ~10 %). v in q0-units of length per unit time; null = the seed dies out,
// 0 = pinned (|v| < 0.01: the patch survives but does not grow), < 0 = the front retreats.
// Circular fronts are slower: dR/dt = v − D/R (q0 units) with D ≈ 0.8 + 3.4 v, fitted to radial runs of the same
// seed at 8 (r, b) points (rms 0.03–0.1 a over R = 1.6 … 10 a). Critical radius R_c = D / v: smaller patches shrink.
export const FRONT = {
  r: [-0.5, -0.475, -0.45, -0.425, -0.4, -0.375, -0.35, -0.325, -0.3, -0.275, -0.25, -0.225, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0.0, 0.025, 0.05, 0.075, 0.1],
  b: [0.5, 0.6, 0.69, 0.75, 0.85, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
  v: [[null,null,null,null,null,null,null,null,null,null,0.0],[null,null,null,null,null,null,null,null,null,-0.029,0.0],[null,null,null,null,null,null,null,null,null,0.0,0.0],[null,null,null,null,null,null,null,null,null,0.0,0.0],[null,null,null,null,null,null,null,null,null,0.0,0.0],[null,null,null,null,null,null,null,null,null,0.0,0.038],[null,null,null,null,null,null,null,null,-0.092,0.0,0.21],[null,null,null,null,null,null,null,null,0.0,0.038,0.317],[null,null,null,null,null,null,null,null,0.0,0.18,0.41],[null,null,null,null,null,null,null,null,0.0,0.28,0.495],[null,null,null,null,null,null,null,null,0.133,0.368,0.577],[null,null,null,null,null,null,null,0.0,0.236,0.453,0.654],[null,null,null,null,null,null,null,0.072,0.329,0.535,0.729],[null,null,null,null,null,null,null,0.189,0.418,0.615,0.796],[null,null,null,null,null,null,0.0,0.293,0.504,0.693,0.872],[null,null,null,null,null,null,0.152,0.393,0.589,0.769,0.944],[null,null,null,null,null,-0.078,0.278,0.488,0.673,0.846,1.011],[null,null,null,null,null,0.156,0.392,0.582,0.753,0.917,1.078],[null,null,null,-0.098,0.129,0.311,0.501,0.67,0.834,0.993,1.146],[null,0.028,0.156,0.221,0.316,0.444,0.604,0.758,0.911,1.058,1.213],[0.25,0.309,0.364,0.401,0.467,0.567,0.705,0.846,0.99,1.141,1.28],[0.479,0.505,0.536,0.562,0.608,0.688,0.809,0.935,1.07,1.209,1.348],[0.786,0.757,0.736,0.732,0.754,0.816,0.916,1.031,1.152,1.285,1.414],[1.021,1.01,1.009,1.002,0.994,0.997,1.035,1.131,1.241,1.365,1.491],[1.199,1.198,1.201,1.208,1.198,1.21,1.234,1.257,1.34,1.443,1.558]],
  // relative front speed against dt at (r, b) = (−0.08, 1.2) and (−0.027, 0.69), same runs (dt → 0 limit = 1)
  dt: [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5],
  dtRel: [[0.9989, 0.9968, 0.9882, 0.9810, 0.9729, 0.9583, 0.9291], [0.9993, 0.9971, 0.9963, 0.9956, 0.9926, 0.9904, 0.9808]],
};
