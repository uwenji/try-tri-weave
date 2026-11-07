// Compact colour maps: viridis (field u) and a cool–warm diverging map (curvature, charge), as 256-entry LUTs.
const VIRIDIS = [[0.267, 0.005, 0.329], [0.283, 0.141, 0.458], [0.254, 0.265, 0.530], [0.207, 0.372, 0.553], [0.164, 0.471, 0.558], [0.128, 0.567, 0.551], [0.135, 0.659, 0.518], [0.267, 0.749, 0.441], [0.478, 0.821, 0.318], [0.741, 0.873, 0.150], [0.993, 0.906, 0.144]];
const COOLWARM = [[0.230, 0.299, 0.754], [0.406, 0.538, 0.934], [0.602, 0.731, 0.999], [0.790, 0.845, 0.939], [0.867, 0.867, 0.867], [0.941, 0.790, 0.700], [0.958, 0.603, 0.482], [0.867, 0.389, 0.282], [0.706, 0.016, 0.150]];

function buildLUT(anchors, n = 256) {
  const lut = new Float32Array(n * 3);
  for (let i = 0; i < n; i++) {
    const t = i / (n - 1) * (anchors.length - 1), k = Math.min(Math.floor(t), anchors.length - 2), f = t - k;
    for (let c = 0; c < 3; c++) lut[3 * i + c] = anchors[k][c] * (1 - f) + anchors[k + 1][c] * f;
  }
  return lut;
}
// cyclic (matplotlib twilight): lattice orientation θ mod 60° — the two ends are the same colour
const TWILIGHT = [[0.886, 0.85, 0.888], [0.768, 0.81, 0.833], [0.583, 0.71, 0.779], [0.448, 0.592, 0.756], [0.384, 0.461, 0.731], [0.37, 0.316, 0.677], [0.351, 0.166, 0.561], [0.27, 0.075, 0.361], [0.185, 0.079, 0.213], [0.291, 0.075, 0.258], [0.454, 0.116, 0.31], [0.597, 0.207, 0.313], [0.698, 0.338, 0.322], [0.763, 0.487, 0.387], [0.8, 0.641, 0.537], [0.849, 0.78, 0.747], [0.886, 0.85, 0.886]];
// misorientation to the neighbours, 0 → 15°: quiet grey inside a grain, amber, then wine red along seams
const MISOR = [[0.62, 0.64, 0.67], [0.70, 0.66, 0.52], [0.867, 0.667, 0.2], [0.80, 0.40, 0.28], [0.73, 0.2, 0.33]];
export const LUT = { viridis: buildLUT(VIRIDIS), coolwarm: buildLUT(COOLWARM), twilight: buildLUT(TWILIGHT), misor: buildLUT(MISOR) };
/** misorientation colour for the orientation legend, t = (mis − 2°) / 13° clipped to [0, 1]: the lattice-line colour of
 *  the theme (base, [r, g, b]) inside a grain (below 2°) → amber at 8.5° → wine red from 15° on (seams) */
export function misorRGB(t, base) {
  const amber = [0.867, 0.667, 0.2], red = [0.733, 0.2, 0.333], u = Math.max(0, Math.min(1, t));
  const [p, q, f] = u < 0.5 ? [base, amber, u / 0.5] : [amber, red, (u - 0.5) / 0.5];
  return [p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1]), p[2] + f * (q[2] - p[2])];
}
/** one colour of a LUT at t ∈ [0, 1] → [r, g, b] */
export function lutColor(name, t) { const lut = LUT[name]; let k = Math.round(Math.max(0, Math.min(1, t)) * 255); return [lut[3 * k], lut[3 * k + 1], lut[3 * k + 2]]; }

/** write colours of values into out (Float32Array 3n) using the LUT over [lo, hi] */
export function mapColors(values, out, name, lo, hi) {
  const lut = LUT[name], n = values.length, scale = 255 / Math.max(hi - lo, 1e-12);
  for (let i = 0; i < n; i++) {
    let k = Math.round((values[i] - lo) * scale); k = k < 0 ? 0 : k > 255 ? 255 : k;
    out[3 * i] = lut[3 * k]; out[3 * i + 1] = lut[3 * k + 1]; out[3 * i + 2] = lut[3 * k + 2];
  }
  return out;
}

export const CHARGE_COLORS = { 3: [0.62, 0.12, 0.12], 2: [0.80, 0.30, 0.18], 1: [0.84, 0.55, 0.22], 0: [0.70, 0.70, 0.70], '-1': [0.24, 0.48, 0.70], '-2': [0.14, 0.32, 0.56], '-3': [0.08, 0.20, 0.40] };
export function chargeColor(charge) { const c = Math.max(-3, Math.min(3, charge | 0)); return CHARGE_COLORS[c]; }

// Paul Tol's colour schemes (https://personal.sron.nl/~pault/)
const hex = (h) => [parseInt(h.slice(1, 3), 16) / 255, parseInt(h.slice(3, 5), 16) / 255, parseInt(h.slice(5, 7), 16) / 255];
/** high-contrast: deep blue, warm gold, muted crimson — the three lath families; the 4th entry is for unassigned strands */
export const FAMILY_HEX = ['#004488', '#DDAA33', '#BB5566', '#8C8C8C'];
export const FAMILY_COLORS = FAMILY_HEX.map(hex);
/** muted: one colour per strand when colouring "by strand" */
export const STRAND_HEX = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499'];
export const STRAND_COLORS = STRAND_HEX.map(hex);
/** light: fills that keep dark text readable — one per valence (singular vertex discs) */
export const VALENCE_HEX = { 3: '#EE8866', 4: '#EEDD88', 5: '#FFAABB', 6: '#DDDDDD', 7: '#77AADD', 8: '#99DDFF', 9: '#44BB99', 10: '#BBCC33', 11: '#AAAA00', 12: '#AAAA00' };
export function valenceHex(v) { return VALENCE_HEX[v] || (v < 6 ? '#EE8866' : '#44BB99'); }
export function valenceColor(v) { return hex(valenceHex(v)); }
export const PLAIN_LATH = hex('#8C8C8C');
