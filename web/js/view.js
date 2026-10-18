// three.js view: carrier mesh coloured by the field (with an opacity toward the plain mesh colour) or by the three
// stripe phases of the global seed, wireframe, lattice edges, peaks, singular vertices as labelled discs (one layer
// per valence), the centre marker, the selected vertex, and the lath ribbons (Kagome strands) with adjustable width,
// over/under offset, family colouring and a highlighted (selected) path. Z is up.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { mapColors, lutColor, FAMILY_COLORS, STRAND_COLORS, PLAIN_LATH, valenceHex } from './colormap.js';
import { NavCube } from './navcube.js';

const cssVar = (name, fallback) => getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback;
const cssColor = (name, fallback) => new THREE.Color(cssVar(name, fallback));
/** a CSS colour token as raw [0, 1] numbers: the convention of the carrier's vertex colours (no sRGB → linear step) */
const cssRaw = (name, fallback) => { const m = /^#?([0-9a-f]{6})$/i.exec(cssVar(name, fallback)) || /^#?([0-9a-f]{6})$/i.exec(fallback); const n = parseInt(m[1], 16); return new THREE.Vector3(((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255); };
/** a CSS colour token as linear RGB (shown as the token's own colour) */
const cssLinear = (name, fallback) => { const c = cssColor(name, fallback); return new THREE.Vector3(c.r, c.g, c.b); };
/** stripe colours, one per colour slot of the stripe phases: red, yellow, blue (the lath family tokens, true colour) */
const STRIPE_TOKENS = [['--fam-c', '#BB5566'], ['--fam-b', '#DDAA33'], ['--fam-a', '#004488']];
const Z_UP = new THREE.Vector3(0, 0, 1);
const HOME_DIR = new THREE.Vector3(-0.85, -1.15, 0.7).normalize();          // the default three-quarter view from the front-left, above
const POLE_EPS = 2e-3;                                                        // z-up: top/bottom views keep this tilt so the azimuth stays defined

export class View {
  constructor(canvas) {
    this.canvas = canvas;
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: false, powerPreference: 'high-performance' });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.scene = new THREE.Scene();
    this.persp = new THREE.PerspectiveCamera(38, 1, 0.01, 1000); this.persp.up.set(0, 0, 1);
    this.orthoH = 1;                                                          // ortho half-height at zoom 1; the zoom carries the scale
    this.ortho = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.01, 1000); this.ortho.up.set(0, 0, 1);
    this.camera = this.persp; this.projection = 'persp';
    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = true; this.controls.dampingFactor = 0.12; this.controls.rotateSpeed = 0.7;
    this.clock = new THREE.Clock(); this.anim = null; this.turntable = false; this.turnSpeed = 12;
    canvas.addEventListener('pointerdown', () => { this.anim = null; });      // a drag cancels a running camera move
    this.hemi = new THREE.HemisphereLight(0xffffff, 0x9a9ea8, 1.7); this.scene.add(this.hemi);
    this.key = new THREE.DirectionalLight(0xffffff, 1.6); this.key.position.set(-1.2, -1.6, 2.2); this.scene.add(this.key);
    this.fill = new THREE.DirectionalLight(0xffffff, 0.9); this.fill.position.set(1.5, 0.8, 0.6); this.scene.add(this.fill);
    this.group = new THREE.Group(); this.scene.add(this.group);
    this.layers = {};                                     // name → Object3D
    this.visible = { field: true, wire: false, probe: false, lattice: true, peaks: false, singular: true, centre: true, laths: false, strain: false, orient: false, rosy: false, cones: true };
    this.fieldAlpha = 1;                                  // carrier colours: 1 = the colour map, 0 = the plain mesh colour
    this.stripeAlpha = [1, 1, 1]; this.stripeHalfW = 0.65; // stripes: opacity per colour slot, band half-width in radians of phase
    this.carrierMode = 'vertex';                          // 'vertex' (vertex colours) | 'stripes' (per-face phases, shaded per pixel)
    this.showValence = {};                                // valence → bool (missing = true)
    this.sideSign = 1;                                    // +1: overlays on the side of the carrier normals, −1: the other side
    this.raycaster = new THREE.Raycaster();
    this.labelCache = new Map();
    this.applyTheme();
    this._ro = new ResizeObserver(() => this.resize()); this._ro.observe(canvas.parentElement);
    this.resize();
  }

  applyTheme() {
    this.dark = cssVar('--scheme', 'light') === 'dark';
    this.renderer.setClearColor(cssColor('--canvas', '#f6f6f5'), 1);
    this.colCanvas = cssColor('--canvas', '#f6f6f5');
    this.colMesh = cssColor('--mesh', '#e0e0de'); this.colEdge = cssColor('--mesh-edge', '#b8b8b5'); this.colKnot = cssColor('--knot', '#404045');
    this.colAccent = cssColor('--accent', '#2973c7'); this.colFg = cssColor('--fg', '#212324');
    if (this.layers.wire) { this.layers.wire.children[0].material.color.copy(this.colKnot); this.layers.wire.children[1].material.color.copy(this.colAccent); }
    if (this.layers.probe) this.layers.probe.traverse(o => { if (o.material && o.userData.role === 'fg') o.material.color.copy(this.colFg); if (o.material && o.userData.role === 'accent') o.material.color.copy(this.colAccent); });
    if (this.layers.lattice) this.layers.lattice.material.color.copy(this.colKnot);
    if (this.matField) this.matField.userData.u.uBase.value.copy(cssRaw('--mesh', '#e0e0de'));
    if (this.matStripes) { const u = this.matStripes.userData.u; u.uBase.value.copy(cssRaw('--mesh', '#e0e0de')); STRIPE_TOKENS.forEach(([n, f], k) => u.uFam.value[k].copy(cssLinear(n, f))); }
    this.labelCache.clear();
    if (this.nav) this.nav.applyTheme();
  }

  resize() {
    const el = this.canvas.parentElement, w = el.clientWidth || 800, h = el.clientHeight || 600, aspect = w / h;
    this.renderer.setSize(w, h, false);
    this.persp.aspect = aspect; this.persp.updateProjectionMatrix();
    this.ortho.left = -this.orthoH * aspect; this.ortho.right = this.orthoH * aspect; this.ortho.top = this.orthoH; this.ortho.bottom = -this.orthoH; this.ortho.updateProjectionMatrix();
    if (this.nav) this.nav.resize();
  }

  render() {
    const dt = Math.min(this.clock.getDelta(), 0.5);                   // CPU-solver batches block frames; keep the turntable rate honest
    this._stepAnim();
    this.controls.update(dt);
    this.renderer.render(this.scene, this.camera);
    if (this.nav) this.nav.render(this.camera, this.controls.target);
  }

  // --- navigation: view cube, canonical views, turntable, projection -------------------------------------------------
  attachNav(canvas) {
    this.nav = new NavCube(canvas, {
      onPick: (dir) => this.lookFrom(dir),
      onDrag: (dx, dy) => { const h = this.canvas.clientHeight || 600, k = 2 * Math.PI * this.controls.rotateSpeed / h; this.controls._rotateLeft(k * dx); this.controls._rotateUp(k * dy); },
    });
  }

  viewDir() { const d = this.camera.position.clone().sub(this.controls.target); return d.lengthSq() > 1e-18 ? d.normalize() : HOME_DIR.clone(); }
  viewDist() { return this.camera.position.distanceTo(this.controls.target); }
  centre() { const c = this.carrier; return c ? new THREE.Vector3((c.lo[0] + c.hi[0]) / 2, (c.lo[1] + c.hi[1]) / 2, (c.lo[2] + c.hi[2]) / 2) : new THREE.Vector3(); }
  /** camera distance that frames the carrier's bounding sphere in the perspective camera (also used as the ortho reference) */
  fitDistance() { const d = this.carrier ? this.carrier.bbox : 10, aspect = this.persp.aspect || 1; return Math.max(1.59 * d, 1.5 * d / aspect); }
  fitZoom() { const hh = this.fitDistance() * Math.tan(THREE.MathUtils.degToRad(this.persp.fov / 2)); return this.orthoH / hh; }

  /** turn the view to look from direction `dir` (array or Vector3, from the target toward the camera), animated.
   *  opts: {dist, target, zoom, ms, tiltFrom: Vector3} — top/bottom views get a tiny tilt toward −y (front at the
   *  bottom of the screen) or toward `tiltFrom`'s azimuth, because the up vector is always +z. */
  lookFrom(dir, opts = {}) {
    const d = (dir.isVector3 ? dir.clone() : new THREE.Vector3().fromArray(dir)).normalize();
    if (1 - Math.abs(d.z) < 1e-6) {
      const t = opts.tiltFrom && Math.hypot(opts.tiltFrom.x, opts.tiltFrom.y) > 1e-9 ? new THREE.Vector3(opts.tiltFrom.x, opts.tiltFrom.y, 0).normalize() : new THREE.Vector3(0, -1, 0);
      d.set(t.x * POLE_EPS, t.y * POLE_EPS, Math.sign(d.z)).normalize();
    }
    this.anim = {
      t0: performance.now(), ms: opts.ms ?? 450,
      dir0: this.viewDir(), dir1: d, dist0: this.viewDist(), dist1: opts.dist ?? this.viewDist(),
      tg0: this.controls.target.clone(), tg1: opts.target ? opts.target.clone() : this.controls.target.clone(),
      zoom0: this.ortho.zoom, zoom1: opts.zoom ?? this.ortho.zoom,
    };
    if (opts.ms === 0) this._stepAnim(performance.now() + 1);
  }

  _stepAnim(now = performance.now()) {
    const a = this.anim; if (!a) return;
    const t = Math.min(1, (now - a.t0) / Math.max(a.ms, 1)), e = t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
    const dir = slerp(a.dir0, a.dir1, e), dist = a.dist0 + (a.dist1 - a.dist0) * e;
    this.controls.target.copy(a.tg0).lerp(a.tg1, e);
    this.camera.position.copy(this.controls.target).addScaledVector(dir, dist);
    if (this.projection === 'ortho') { this.ortho.zoom = a.zoom0 + (a.zoom1 - a.zoom0) * e; this.ortho.updateProjectionMatrix(); }
    if (t >= 1) this.anim = null;
  }

  /** orbit by an azimuth / polar increment (radians) in the z-up frame; the polar angle stops at the poles */
  orbitBy(dAz, dPol) {
    const d0 = this.viewDir(), az = Math.atan2(d0.y, d0.x) + dAz;
    const pol = Math.min(Math.PI - POLE_EPS, Math.max(POLE_EPS, Math.acos(Math.max(-1, Math.min(1, d0.z))) + dPol));
    this.lookFrom(new THREE.Vector3(Math.sin(pol) * Math.cos(az), Math.sin(pol) * Math.sin(az), Math.cos(pol)), { tiltFrom: d0 });
  }

  home(ms = 500) { this.lookFrom(HOME_DIR, { dist: this.fitDistance(), target: this.centre(), zoom: this.fitZoom(), ms }); }
  fit(ms = 400) { this.lookFrom(this.viewDir(), { dist: this.fitDistance(), target: this.centre(), zoom: this.fitZoom(), ms, tiltFrom: this.viewDir() }); }

  /** continuous rotation about the vertical axis through the target (OrbitControls' autoRotate); speed in °/s */
  setTurntable(on, degPerSec = this.turnSpeed) {
    this.turntable = !!on; this.turnSpeed = degPerSec;
    this.controls.autoRotate = this.turntable; this.controls.autoRotateSpeed = degPerSec / 6;   // 2π/60·speed rad/s
  }

  setProjection(kind) {
    if (kind === this.projection) return;
    const dir = this.viewDir(), tg = this.controls.target.clone(), tanHalf = Math.tan(THREE.MathUtils.degToRad(this.persp.fov / 2));
    this.anim = null;
    if (kind === 'ortho') {
      const hh = this.viewDist() * tanHalf;                                  // visible half-height at the target: keep it
      this.ortho.zoom = this.orthoH / hh; this.ortho.position.copy(this.persp.position); this.ortho.updateProjectionMatrix();
      this.camera = this.ortho;
    } else {
      const hh = this.orthoH / this.ortho.zoom;
      this.persp.position.copy(tg).addScaledVector(dir, hh / tanHalf);
      this.camera = this.persp;
    }
    this.projection = kind; this.controls.object = this.camera; this.camera.lookAt(tg); this.controls.update();
  }

  _setLayer(name, obj) {
    if (this.layers[name]) { this.group.remove(this.layers[name]); disposeObject(this.layers[name]); }
    this.layers[name] = obj;
    if (obj) { obj.visible = !!this.visible[name]; this.group.add(obj); }
  }

  setVisible(name, on) { this.visible[name] = on; if (this.layers[name]) this.layers[name].visible = on; }

  setValenceVisible(v, on) {
    this.showValence[v] = on;
    const g = this.layers.singular; if (!g) return;
    for (const child of g.children) if (child.userData.valence === v) child.visible = on;
  }

  // --- carrier ----------------------------------------------------------------------------------------------
  /** opts.keepCamera: keep the current view (refinement of the same surface) */
  setCarrier(c, opts = {}) {
    this.carrier = c;
    for (const o of [this.geoField, this.stripeGeo, this.matField, this.matStripes]) o?.dispose?.();   // whichever the mesh is not showing
    this.stripeGeo = null; this.carrierMode = 'vertex';
    for (const k of Object.keys(this.layers)) this._setLayer(k, null);
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(Float32Array.from(c.V), 3));
    geo.setIndex(new THREE.BufferAttribute(Uint32Array.from(c.F), 1));
    geo.computeVertexNormals();
    this.colorAttr = new THREE.BufferAttribute(new Float32Array(c.n * 3).fill(0.85), 3);
    this.colorAttr.setUsage(THREE.DynamicDrawUsage);
    geo.setAttribute('color', this.colorAttr);
    this.geoField = geo; this.matField = this._fieldMaterial(); this.matStripes = this._stripeMaterial();
    const mesh = new THREE.Mesh(geo, this.matField); mesh.userData.kind = 'carrier'; this._setLayer('field', mesh); this.mesh = mesh;
    this._setLayer('wire', null); this._setLayer('probe', null);
    if (opts.keepCamera) { const d = c.bbox; for (const cam of [this.persp, this.ortho]) { cam.near = d * 0.01; cam.far = d * 20; cam.updateProjectionMatrix(); } }
    else this.fitCamera();
  }

  /** which side of the surface carries the lattice marks and the laths: +1 along the carrier normals, −1 opposite */
  setSide(sign) { this.sideSign = sign < 0 ? -1 : 1; }

  fitCamera() {
    const c = this.carrier; if (!c) return;
    const d = c.bbox;
    for (const cam of [this.persp, this.ortho]) { cam.near = d * 0.01; cam.far = d * 20; cam.updateProjectionMatrix(); }
    this.home(0);
  }

  setFieldColors(values, cmap, lo, hi) {
    if (!this.colorAttr) return;
    mapColors(values, this.colorAttr.array, cmap, lo, hi); this.colorAttr.needsUpdate = true;
  }

  setFlatColor(rgb) { if (!this.colorAttr) return; const a = this.colorAttr.array; for (let i = 0; i < a.length; i += 3) { a[i] = rgb[0]; a[i + 1] = rgb[1]; a[i + 2] = rgb[2]; } this.colorAttr.needsUpdate = true; }

  // polygon offset on both carrier materials: the mesh lines and the probe draw on the surface without z-fighting
  /** vertex colours mixed toward the plain mesh colour by uAlpha (the opacity slider): no recolouring per change */
  _fieldMaterial() {
    const u = { uAlpha: { value: this.fieldAlpha }, uBase: { value: cssRaw('--mesh', '#e0e0de') } };
    const mat = new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.9, metalness: 0.0, side: THREE.DoubleSide, flatShading: false, polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1 });
    mat.userData.u = u;
    mat.onBeforeCompile = (sh) => {
      Object.assign(sh.uniforms, mat.userData.u);
      sh.fragmentShader = sh.fragmentShader
        .replace('#include <color_pars_fragment>', '#include <color_pars_fragment>\nuniform float uAlpha;\nuniform vec3 uBase;')
        .replace('#include <color_fragment>', 'diffuseColor.rgb *= mix(uBase, vColor, uAlpha);');
    };
    return mat;
  }
  /** the three stripe families as bands where their laths would run (phase ≡ π), one colour and opacity per slot,
   *  over the plain mesh colour; anti-aliased in screen space. Attribute ph: the three phases at each face corner. */
  _stripeMaterial() {
    const u = { uBase: { value: cssRaw('--mesh', '#e0e0de') }, uFam: { value: STRIPE_TOKENS.map(([n, f]) => cssLinear(n, f)) }, uA: { value: new THREE.Vector3(...this.stripeAlpha) }, uHW: { value: this.stripeHalfW } };
    const mat = new THREE.MeshStandardMaterial({ roughness: 0.9, metalness: 0.0, side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1 });
    mat.userData.u = u;
    mat.onBeforeCompile = (sh) => {
      Object.assign(sh.uniforms, mat.userData.u);
      sh.vertexShader = sh.vertexShader
        .replace('#include <common>', '#include <common>\nattribute vec3 ph;\nvarying vec3 vPh;')
        .replace('#include <begin_vertex>', '#include <begin_vertex>\nvPh = ph;');
      sh.fragmentShader = sh.fragmentShader
        .replace('#include <common>', '#include <common>\nvarying vec3 vPh;\nuniform vec3 uBase;\nuniform vec3 uFam[3];\nuniform vec3 uA;\nuniform float uHW;')
        .replace('#include <color_fragment>', [
          'vec3 sc = uBase;',
          'for (int k = 0; k < 3; k++) {',
          '  float p = vPh[k], d = abs(mod(p, 6.2831853) - 3.1415927), aa = max(0.75 * fwidth(p), 1e-4);',
          '  sc = mix(sc, uFam[k], uA[k] * (1.0 - smoothstep(uHW - aa, uHW + aa, d)));',
          '}',
          'diffuseColor.rgb *= sc;'].join('\n'));
    };
    return mat;
  }
  /** opacity of the carrier colours (field u, curvature, field 5/7) over the plain mesh colour, 0 … 1 */
  setFieldAlpha(a) { this.fieldAlpha = a; if (this.matField) this.matField.userData.u.uAlpha.value = a; }
  /** stripes: alpha = [red, yellow, blue] opacities 0 … 1; halfWidth = band half-width in radians of phase */
  setStripeParams({ alpha, halfWidth } = {}) {
    if (alpha) this.stripeAlpha = alpha.slice(0, 3);
    if (halfWidth !== undefined) this.stripeHalfW = halfWidth;
    if (this.matStripes) { const u = this.matStripes.userData.u; u.uA.value.set(...this.stripeAlpha); u.uHW.value = this.stripeHalfW; }
  }
  /** stripe phases per render face (stripes.js · stripeFaces: 9 per face of carrier.F, unwrapped within the face):
   *  a non-indexed copy of the carrier with the phases at the face corners. Returns false if they do not fit. */
  setStripes(phase) {
    if (this.stripeGeo) { if (this.mesh && this.mesh.geometry === this.stripeGeo) this.mesh.geometry = this.geoField; this.stripeGeo.dispose(); this.stripeGeo = null; }
    const c = this.carrier;
    if (!c || !this.geoField || !phase || phase.length !== 3 * c.F.length) { if (this.carrierMode === 'stripes') this.setCarrierMode('vertex'); return false; }
    const F = c.F, nf = F.length / 3, V = c.V, Nv = this.geoField.attributes.normal.array, pos = new Float32Array(9 * nf), nrm = new Float32Array(9 * nf);
    for (let f = 0; f < nf; f++) for (let j = 0; j < 3; j++) {
      const v = F[3 * f + j], o = 9 * f + 3 * j;
      pos[o] = V[3 * v]; pos[o + 1] = V[3 * v + 1]; pos[o + 2] = V[3 * v + 2]; nrm[o] = Nv[3 * v]; nrm[o + 1] = Nv[3 * v + 1]; nrm[o + 2] = Nv[3 * v + 2];
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3)); g.setAttribute('normal', new THREE.BufferAttribute(nrm, 3));
    g.setAttribute('ph', new THREE.BufferAttribute(phase instanceof Float32Array ? phase : Float32Array.from(phase), 3));
    this.stripeGeo = g;
    if (this.carrierMode === 'stripes' && this.mesh) this.mesh.geometry = g;
    return true;
  }
  /** 'stripes' shows the stripe geometry (if set), anything else the vertex-coloured carrier; same mesh object, so
   *  picking and the layer switch are unchanged (face order is the same in both geometries) */
  setCarrierMode(mode) {
    if (!this.mesh) return;
    const stripes = mode === 'stripes' && !!this.stripeGeo;
    this.carrierMode = stripes ? 'stripes' : 'vertex';
    this.mesh.geometry = stripes ? this.stripeGeo : this.geoField;
    this.mesh.material = stripes ? this.matStripes : this.matField;
  }

  // --- labelled discs (sprites) ---------------------------------------------------------------------------------
  /** ring: a hollow disc with a thick border in fillHex (the cones of the size field, distinct from the lattice's own discs) */
  _labelTexture(text, fillHex, ring = false) {
    const key = text + '|' + fillHex + (ring ? '|ring' : '');
    if (this.labelCache.has(key)) return this.labelCache.get(key);
    const S = 128, cv = document.createElement('canvas'); cv.width = S; cv.height = S;
    const ctx = cv.getContext('2d');
    if (ring) {
      ctx.beginPath(); ctx.arc(S / 2, S / 2, S / 2 - 10, 0, 2 * Math.PI); ctx.fillStyle = this.dark ? 'rgba(30,32,34,0.82)' : 'rgba(255,255,255,0.88)'; ctx.fill();
      ctx.lineWidth = 14; ctx.strokeStyle = fillHex; ctx.stroke();
    } else {
      ctx.beginPath(); ctx.arc(S / 2, S / 2, S / 2 - 4, 0, 2 * Math.PI); ctx.fillStyle = fillHex; ctx.fill();
      ctx.lineWidth = 5; ctx.strokeStyle = this.dark ? 'rgba(255,255,255,0.85)' : 'rgba(30,32,34,0.75)'; ctx.stroke();
    }
    const px = text.length <= 1 ? 76 : text.length === 2 ? 62 : text.length === 3 ? 48 : 38;
    ctx.font = `700 ${px}px "IBM Plex Sans", "Helvetica Neue", Arial, sans-serif`; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.fillStyle = ring ? (this.dark ? '#e4e5e3' : '#1b1d1f') : '#1b1d1f'; ctx.fillText(text, S / 2, S / 2 + px * 0.06);
    const tex = new THREE.CanvasTexture(cv); tex.colorSpace = THREE.SRGBColorSpace; tex.anisotropy = 4;
    this.labelCache.set(key, tex); return tex;
  }

  _disc(text, fillHex, pos, size, userData, ring = false) {
    const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: this._labelTexture(text, fillHex, ring), transparent: true, depthWrite: false }));
    sp.position.set(pos[0], pos[1], pos[2]); sp.scale.set(size, size, 1); sp.renderOrder = 5; Object.assign(sp.userData, userData);
    return sp;
  }

  /** the cones of the size field ĝ as ring discs with their valence (grey when switched off, larger when selected);
   *  cones = [{vertex, valence, on}], picked as {kind: 'cone', idx} */
  setCones(cones, selected = -1) {
    const c = this.carrier;
    if (!c || !cones || !cones.length) { this._setLayer('cones', null); return; }
    const g = new THREE.Group(), da = this.a || 8 * c.h;
    cones.forEach((q, i) => {
      const sp = this._disc(String(q.valence), q.on === false ? '#8a8e95' : valenceHex(q.valence), this.lifted(q.vertex, 2.6 * c.h), (i === selected ? 1.4 : 1) * Math.max(2.2 * c.h, 0.34 * da), { kind: 'cone', idx: i }, true);
      sp.renderOrder = 7; g.add(sp);
    });
    this._setLayer('cones', g);
  }

  // --- lattice overlays ------------------------------------------------------------------------------------
  lifted(idx, k) { const c = this.carrier, s = k * this.sideSign; return [c.V[3 * idx] + s * c.N[3 * idx], c.V[3 * idx + 1] + s * c.N[3 * idx + 1], c.V[3 * idx + 2] + s * c.N[3 * idx + 2]]; }

  setLattice(lat, summary) {
    const c = this.carrier;
    this.lat = lat;
    if (!lat) { this._setLayer('lattice', null); this._setLayer('peaks', null); this._setLayer('singular', null); this._setLayer('selected', null); return; }
    const lift = 0.6 * c.h, npk = lat.peaks.length;
    const P = new Float32Array(npk * 3);
    for (let i = 0; i < npk; i++) { const q = this.lifted(lat.peaks[i], lift); P[3 * i] = q[0]; P[3 * i + 1] = q[1]; P[3 * i + 2] = q[2]; }
    this.latticeP = P;
    // edges
    const seen = new Set(), segs = [];
    const T = lat.tris;
    for (let t = 0; t < T.length / 3; t++) for (const [i, j] of [[T[3 * t], T[3 * t + 1]], [T[3 * t + 1], T[3 * t + 2]], [T[3 * t + 2], T[3 * t]]]) {
      const k = i < j ? i * 1e5 + j : j * 1e5 + i; if (seen.has(k)) continue; seen.add(k);
      segs.push(P[3 * i], P[3 * i + 1], P[3 * i + 2], P[3 * j], P[3 * j + 1], P[3 * j + 2]);
    }
    const lg = new THREE.BufferGeometry(); lg.setAttribute('position', new THREE.BufferAttribute(Float32Array.from(segs), 3));
    this._setLayer('lattice', new THREE.LineSegments(lg, new THREE.LineBasicMaterial({ color: this.colKnot, transparent: true, opacity: 0.9 })));
    // peaks (small spheres)
    const pk = spheres(P, npk, 0.28 * c.h, (i, col) => col.copy(this.colKnot)); pk.userData.kind = 'peaks';
    this._setLayer('peaks', pk);
    // singular vertices: labelled discs, one child per vertex, grouped by valence through userData
    const g = new THREE.Group();
    for (const i of summary.singular) {
      const da = this.a || 8 * c.h;                                        // disc size follows the lattice edge
      const v = lat.valence[i], q = this.lifted(lat.peaks[i], 1.6 * c.h);
      const sp = this._disc(String(v), valenceHex(v), q, Math.max(1.9 * c.h, 0.3 * da), { kind: 'vertex', peak: i, valence: v });
      sp.visible = this.showValence[v] !== false;
      g.add(sp);
    }
    this._setLayer('singular', g);
  }

  // --- global rotation field (no seed) -------------------------------------------------------------------------
  /** strokes of the three lattice directions at every carrier vertex (dirs: 3 × Float64Array(3n) unit vectors), length len */
  setRosy(dirs, len) {
    const c = this.carrier;
    if (!c || !dirs) { this._setLayer('rosy', null); return; }
    const n = c.n, lift = 0.4 * c.h, segs = new Float32Array(n * 3 * 6); let k = 0;
    for (let i = 0; i < n; i++) {
      const p = this.lifted(i, lift);
      for (const d of dirs) { const dx = d[3 * i] * len / 2, dy = d[3 * i + 1] * len / 2, dz = d[3 * i + 2] * len / 2;
        segs[k++] = p[0] - dx; segs[k++] = p[1] - dy; segs[k++] = p[2] - dz; segs[k++] = p[0] + dx; segs[k++] = p[1] + dy; segs[k++] = p[2] + dz; }
    }
    const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.BufferAttribute(segs, 3));
    this._setLayer('rosy', new THREE.LineSegments(g, new THREE.LineBasicMaterial({ color: this.colFg, transparent: true, opacity: 0.55 })));
  }
  /** larger disc with the lattice index on the selected peak (−1 clears) */
  setSelectedVertex(peak) {
    const c = this.carrier, lat = this.lat;
    if (!lat || peak < 0 || peak >= lat.peaks.length) { this._setLayer('selected', null); return; }
    const da = this.a || 8 * c.h;
    const v = lat.valence[peak], q = this.lifted(lat.peaks[peak], 1.9 * c.h);
    const sp = this._disc(String(peak), valenceHex(v), q, Math.max(2.9 * c.h, 0.46 * da), { kind: 'vertex', peak, valence: v });
    sp.renderOrder = 6;
    const g = new THREE.Group(); g.add(sp);
    this.visible.selected = true; this._setLayer('selected', g);
  }

  setCentre(idx) {
    const c = this.carrier;
    if (!c || idx < 0) { this._setLayer('centre', null); return; }
    const q = this.lifted(idx, 0.9 * c.h);                                  // below the valence disc of the same vertex
    this._setLayer('centre', spheres(Float32Array.from(q), 1, 0.42 * c.h, (i, col) => col.copy(this.colAccent)));
  }

  /** the carrier's intrinsic Delaunay triangulation: shared edges in the mesh colour, flipped (geodesic) edges in the
   *  accent colour. lines = carrier.intrinsicEdgeLines() */
  setMesh(lines) {
    if (!lines) { this._setLayer('wire', null); return; }
    const g = new THREE.Group();
    const mk = (arr, color, opacity) => { const geo = new THREE.BufferGeometry(); geo.setAttribute('position', new THREE.BufferAttribute(arr, 3));
      const ls = new THREE.LineSegments(geo, new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity })); ls.renderOrder = 1; return ls; };
    g.add(mk(lines.shared, this.colKnot, this.dark ? 0.32 : 0.28)); g.add(mk(lines.flipped, this.colAccent, 1));
    this._setLayer('wire', g);
  }

  /** scale probe at the centre: the seed's n neighbours at lattice distance a, spokes, the circle ρ = a */
  setProbe(probe) {
    const c = this.carrier;
    if (!c || !probe || !probe.sites.length) { this._setLayer('probe', null); return; }
    const g = new THREE.Group();
    const segs = []; for (const sp of probe.spokes) for (let k = 0; k + 1 < sp.length; k++) segs.push(...sp[k], ...sp[k + 1]);
    const lgeo = new THREE.BufferGeometry(); lgeo.setAttribute('position', new THREE.BufferAttribute(Float32Array.from(segs), 3));
    const spokes = new THREE.LineSegments(lgeo, new THREE.LineBasicMaterial({ color: this.colAccent })); spokes.renderOrder = 3; spokes.userData.role = 'accent';
    const rgeo = new THREE.BufferGeometry(); rgeo.setAttribute('position', new THREE.BufferAttribute(Float32Array.from(probe.ring.flat()), 3));
    const ring = new THREE.LineSegments(rgeo, new THREE.LineDashedMaterial({ color: this.colAccent, dashSize: 0.25 * c.h, gapSize: 0.2 * c.h })); ring.computeLineDistances(); ring.renderOrder = 3; ring.userData.role = 'accent';
    g.add(spokes, ring);
    const P = new Float32Array(probe.sites.length * 3); probe.sites.forEach((p, k) => { P[3 * k] = p[0]; P[3 * k + 1] = p[1]; P[3 * k + 2] = p[2]; });
    const dots = spheres(P, probe.sites.length, 0.45 * c.h, (i, col) => col.copy(this.colAccent)); dots.userData.role = 'accent-inst'; g.add(dots);
    // the label "a" beside the first spoke
    const sp0 = probe.spokes[0], mid = sp0[Math.floor(sp0.length / 2)];
    const lab = new THREE.Sprite(new THREE.SpriteMaterial({ map: this._textTexture('a'), transparent: true, depthWrite: false, depthTest: false }));
    lab.position.set(mid[0], mid[1], mid[2]); const s = 2.2 * c.h; lab.scale.set(s, s, 1); lab.renderOrder = 6; g.add(lab);
    this._setLayer('probe', g);
  }
  _textTexture(text) {
    const S = 96, cv = document.createElement('canvas'); cv.width = S; cv.height = S; const ctx = cv.getContext('2d');
    ctx.font = 'italic 600 64px "IBM Plex Sans", Arial, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.lineWidth = 10; ctx.strokeStyle = this.dark ? '#1d1f22' : '#ffffff'; ctx.strokeText(text, S / 2, S / 2);
    ctx.fillStyle = '#' + this.colAccent.getHexString(); ctx.fillText(text, S / 2, S / 2);
    const tex = new THREE.CanvasTexture(cv); tex.colorSpace = THREE.SRGBColorSpace; return tex;
  }

  /** lath ribbons. opts: {width, overUnder, colouring: 'family'|'strand'|'plain', selected: strand index or -1} */
  setLaths(lat, strandsData, opts) {
    const c = this.carrier;
    if (!lat || !strandsData || !strandsData.strands.length) { this._setLayer('laths', null); this._setLayer('path', null); this.lathFaceStrand = null; return; }
    const { edges, strands, family, over: overRaw } = strandsData, pk = lat.peaks, N = c.N, sg = this.sideSign;
    const over = sg > 0 ? overRaw : (si, i) => !overRaw(si, i);            // seen from the other side, over becomes under
    const w = opts.width, thk = opts.thickness ?? 0.35 * w, lift = 0.45 * c.h, sub = 2, sel = opts.selected ?? -1;
    const pos = [], nrm = [], col = [], idx = [], faceStrand = [];
    const dim = [this.colCanvas.r, this.colCanvas.g, this.colCanvas.b];
    const edgeMid = (e) => {
      const [a, b] = edges[e], va = pk[a], vb = pk[b];
      const m = [0.5 * (c.V[3 * va] + c.V[3 * vb]), 0.5 * (c.V[3 * va + 1] + c.V[3 * vb + 1]), 0.5 * (c.V[3 * va + 2] + c.V[3 * vb + 2])];
      let n = [N[3 * va] + N[3 * vb], N[3 * va + 1] + N[3 * vb + 1], N[3 * va + 2] + N[3 * vb + 2]]; const l = Math.hypot(...n) || 1; n = n.map(v => sg * v / l);
      return { m, n };
    };
    let selChain = null, selClosed = false;
    strands.forEach((s, si) => {
      const K = s.knots.length; if (K < 2) return;
      const nodes = s.knots.map(edgeMid);
      const chain = [];
      const nseg = s.closed ? K : K - 1;
      for (let i = 0; i < nseg; i++) {
        const a = nodes[i], b = nodes[(i + 1) % K];
        const offA = lift + (opts.overUnder ? (over(si, i) ? thk : 0) : thk / 2);
        chain.push({ m: a.m, n: a.n, off: offA });
        for (let k = 1; k < sub; k++) { const t = k / sub; chain.push({ m: lerp3(a.m, b.m, t), n: nrm3(lerp3(a.n, b.n, t)), off: lift + thk / 2 }); }
      }
      if (!s.closed) { const last = nodes[K - 1]; chain.push({ m: last.m, n: last.n, off: lift + (opts.overUnder ? (over(si, K - 1) ? thk : 0) : thk / 2) }); }
      const M = chain.length, base0 = pos.length / 3;
      let rgb = opts.colouring === 'plain' ? PLAIN_LATH : opts.colouring === 'strand' ? STRAND_COLORS[si % STRAND_COLORS.length] : FAMILY_COLORS[family[si] >= 0 ? family[si] % 3 : 3];
      if (this.dark) rgb = rgb.map(v => 0.82 * v + 0.18);                        // lift a little on a dark canvas
      if (sel >= 0 && si !== sel) rgb = rgb.map((v, k) => 0.42 * v + 0.58 * dim[k]); // dim everything but the selected path
      if (si === sel) { selChain = chain; selClosed = s.closed; }
      for (let i = 0; i < M; i++) {
        const prev = chain[(i - 1 + M) % M], next = chain[(i + 1) % M];
        const tp = s.closed ? next.m.map((v, k) => v - prev.m[k]) : (i === 0 ? next.m.map((v, k) => v - chain[i].m[k]) : i === M - 1 ? chain[i].m.map((v, k) => v - prev.m[k]) : next.m.map((v, k) => v - prev.m[k]));
        const t = nrm3(tp), n = chain[i].n, b = nrm3(cross3(n, t));
        const cm = chain[i].m.map((v, k) => v + chain[i].off * n[k]);
        for (const sgn of [-1, 1]) { pos.push(cm[0] + sgn * 0.5 * w * b[0], cm[1] + sgn * 0.5 * w * b[1], cm[2] + sgn * 0.5 * w * b[2]); nrm.push(n[0], n[1], n[2]); col.push(rgb[0], rgb[1], rgb[2]); }
      }
      const segsN = s.closed ? M : M - 1;
      for (let i = 0; i < segsN; i++) {
        const a0 = base0 + 2 * i, a1 = base0 + 2 * ((i + 1) % M);
        idx.push(a0, a0 + 1, a1, a0 + 1, a1 + 1, a1); faceStrand.push(si, si);
      }
    });
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(Float32Array.from(pos), 3));
    g.setAttribute('normal', new THREE.BufferAttribute(Float32Array.from(nrm), 3));
    g.setAttribute('color', new THREE.BufferAttribute(Float32Array.from(col), 3));
    g.setIndex(idx);
    const mat = new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.8, metalness: 0, side: THREE.DoubleSide });
    const mesh = new THREE.Mesh(g, mat); mesh.userData.kind = 'lath';
    this.lathFaceStrand = Int32Array.from(faceStrand);
    this._setLayer('laths', mesh);
    // highlighted path: a tube along the selected strand's centreline, in the foreground colour
    if (selChain) {
      const pts = selChain.map(q => new THREE.Vector3(q.m[0] + (q.off + 0.35 * thk) * q.n[0], q.m[1] + (q.off + 0.35 * thk) * q.n[1], q.m[2] + (q.off + 0.35 * thk) * q.n[2]));
      const curve = new THREE.CatmullRomCurve3(pts, selClosed, 'centripetal', 0.5);
      const tube = new THREE.Mesh(new THREE.TubeGeometry(curve, Math.max(16, pts.length * 4), 0.11 * Math.max(w, c.h), 6, selClosed),
        new THREE.MeshBasicMaterial({ color: this.colFg }));
      tube.renderOrder = 4; this.visible.path = true; this._setLayer('path', tube);
    } else this._setLayer('path', null);
  }

  // --- per-direction read-outs ---------------------------------------------------------------------------------
  /** lattice edges as thin bars coloured by strain ℓ/a_target − 1 (cool–warm over ±range), one family or all.
   *  es = edgeStrain(...) ; opts: {family: -1 (all) | 0 | 1 | 2, range: 0.15} */
  setStrain(lat, strandsData, es, opts = {}) {
    const c = this.carrier;
    if (!lat || !strandsData || !es) { this._setLayer('strain', null); return; }
    const fam = opts.family ?? -1, range = opts.range ?? 0.15, { edges } = strandsData, pk = lat.peaks;
    const keep = []; for (let e = 0; e < edges.length; e++) if (es.fam[e] >= 0 && (fam < 0 || es.fam[e] === fam)) keep.push(e);
    if (!keep.length) { this._setLayer('strain', null); return; }
    const da = this.a || 8 * c.h, radius = Math.max(0.035 * da, 0.22 * c.h), lift = 1.0 * c.h;
    const geo = new THREE.CylinderGeometry(radius, radius, 1, 6, 1), mat = new THREE.MeshStandardMaterial({ roughness: 0.7, metalness: 0 });
    const inst = new THREE.InstancedMesh(geo, mat, keep.length), m4 = new THREE.Matrix4(), q = new THREE.Quaternion(), col = new THREE.Color();
    const Y = new THREE.Vector3(0, 1, 0), A = new THREE.Vector3(), B = new THREE.Vector3(), D = new THREE.Vector3(), M = new THREE.Vector3(), SC = new THREE.Vector3();
    keep.forEach((e, k) => {
      const [i, j] = edges[e];
      A.fromArray(this.lifted(pk[i], lift)); B.fromArray(this.lifted(pk[j], lift));
      D.subVectors(B, A); const len = D.length() || 1e-9; D.divideScalar(len); M.addVectors(A, B).multiplyScalar(0.5);
      q.setFromUnitVectors(Y, D); SC.set(1, len * 0.92, 1);
      m4.compose(M, q, SC); inst.setMatrixAt(k, m4);
      const rgb = lutColor('coolwarm', 0.5 + 0.5 * es.strain[e] / range); col.setRGB(rgb[0], rgb[1], rgb[2]); inst.setColorAt(k, col);
    });
    inst.instanceMatrix.needsUpdate = true; if (inst.instanceColor) inst.instanceColor.needsUpdate = true;
    inst.userData.kind = 'strain';
    this._setLayer('strain', inst);
  }

  /** lattice orientation at the regular (6-valent, interior) sites as thin flat arrows in the tangent plane. The arrow
   *  is the local wave vector of the pattern — perpendicular to one family of lattice edges, i.e. 30° off the lattice
   *  directions, so it never hides under a lattice line — taken from the consistently transported lattice direction
   *  (orientationField().arrow). Colour: 'mis' = misorientation to the neighbours (the lattice-line colour inside a
   *  grain → amber → red along seams) or 'angle' = θ mod 60° (cyclic). */
  setOrientation(lat, orient, mode = 'mis') {
    const c = this.carrier;
    if (!lat || !orient || !orient.arrow) { this._setLayer('orient', null); return; }
    const { theta, arrow, mis } = orient, idx = []; for (let i = 0; i < theta.length; i++) if (Number.isFinite(theta[i])) idx.push(i);
    if (!idx.length) { this._setLayer('orient', null); return; }
    const da = this.a || 8 * c.h, back = 0.24 * da, front = 0.3 * da, headL = 0.13 * da, headW = 0.065 * da, shaftW = 0.016 * da;
    const lift = 1.0 * c.h, sg = this.sideSign, c30 = Math.cos(Math.PI / 6), s30 = Math.sin(Math.PI / 6);
    const pos = new Float32Array(idx.length * 27), col = new Float32Array(idx.length * 27);   // 3 triangles per arrow
    const amber = new THREE.Color('#DDAA33'), red = new THREE.Color('#BB3355'), cc = new THREE.Color();
    idx.forEach((i, k) => {
      const v = lat.peaks[i], M = this.lifted(v, lift), n = [sg * c.N[3 * v], sg * c.N[3 * v + 1], sg * c.N[3 * v + 2]];
      const l = [arrow[3 * i], arrow[3 * i + 1], arrow[3 * i + 2]], gl = cross3(n, l);
      const d = [0, 1, 2].map(q => c30 * l[q] + s30 * gl[q]), g = cross3(n, d);
      const P = (s, t) => [M[0] + s * d[0] + t * g[0], M[1] + s * d[1] + t * g[1], M[2] + s * d[2] + t * g[2]];
      const t0 = P(-back, -shaftW), t1 = P(-back, shaftW), h0 = P(front - headL, -shaftW), h1 = P(front - headL, shaftW);
      const w0 = P(front - headL, -headW), w1 = P(front - headL, headW), tip = P(front, 0);
      [t0, h0, h1, t0, h1, t1, w0, tip, w1].forEach((p, m) => pos.set(p, 27 * k + 3 * m));
      if (mode === 'angle') { const rgb = lutColor('twilight', (theta[i] + 30) / 60); cc.setRGB(rgb[0], rgb[1], rgb[2]); }
      else { const t = Math.max(0, Math.min(1, ((mis[i] || 0) - 2) / 13)); if (t < 0.5) cc.lerpColors(this.colKnot, amber, t / 0.5); else cc.lerpColors(amber, red, (t - 0.5) / 0.5); }
      for (let m = 0; m < 9; m++) col.set([cc.r, cc.g, cc.b], 27 * k + 3 * m);
    });
    const geo = new THREE.BufferGeometry(); geo.setAttribute('position', new THREE.BufferAttribute(pos, 3)); geo.setAttribute('color', new THREE.BufferAttribute(col, 3));
    const mesh = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide }));
    mesh.renderOrder = 3; mesh.userData.kind = 'orient'; mesh.userData.map = Int32Array.from(idx);   // 3 faces per arrow
    this._setLayer('orient', mesh);
  }

  // --- picking ------------------------------------------------------------------------------------------------
  /** returns {kind: 'lath', strand} | {kind: 'vertex', peak} | {kind: 'cone', idx} | {kind: 'carrier', vertex} | null */
  pick(clientX, clientY) {
    if (!this.mesh) return null;
    const r = this.canvas.getBoundingClientRect();
    const nd = new THREE.Vector2(((clientX - r.left) / r.width) * 2 - 1, -((clientY - r.top) / r.height) * 2 + 1);
    this.raycaster.setFromCamera(nd, this.camera);
    const targets = [];
    if (this.layers.selected && this.layers.selected.visible) targets.push(...this.layers.selected.children);
    if (this.layers.cones && this.layers.cones.visible) targets.push(...this.layers.cones.children);
    if (this.layers.singular && this.layers.singular.visible) targets.push(...this.layers.singular.children.filter(o => o.visible));
    if (this.layers.laths && this.layers.laths.visible) targets.push(this.layers.laths);
    if (this.layers.peaks && this.layers.peaks.visible) targets.push(this.layers.peaks);
    if (this.layers.orient && this.layers.orient.visible) targets.push(this.layers.orient);
    const overlayHits = this.raycaster.intersectObjects(targets, false), meshHits = this.raycaster.intersectObject(this.mesh, false);
    // the opaque field hides what is behind it; a hidden or wire-only carrier does not (what you see is what you click)
    const occluder = this.visible.field ? meshHits[0] : null;
    let h;
    if (overlayHits.length && (!occluder || overlayHits[0].distance <= occluder.distance + 1e-6)) h = overlayHits[0];
    else if (meshHits.length) h = meshHits[0];
    else return null;
    const o = h.object, kind = o.userData.kind;
    if (kind === 'cone') return { kind: 'cone', idx: o.userData.idx, point: h.point };
    if (kind === 'vertex') return { kind: 'vertex', peak: o.userData.peak, point: h.point };
    if (kind === 'lath') return { kind: 'lath', strand: this.lathFaceStrand ? this.lathFaceStrand[h.faceIndex] : -1, point: h.point };
    if (kind === 'peaks') return { kind: 'vertex', peak: h.instanceId, point: h.point };
    if (kind === 'orient') return { kind: 'vertex', peak: o.userData.map[Math.floor(h.faceIndex / 3)], point: h.point };
    const F = this.carrier.F, V = this.carrier.V, f = h.faceIndex;
    let best = -1, bd = Infinity;
    for (const vi of [F[3 * f], F[3 * f + 1], F[3 * f + 2]]) { const d = h.point.distanceTo(new THREE.Vector3(V[3 * vi], V[3 * vi + 1], V[3 * vi + 2])); if (d < bd) { bd = d; best = vi; } }
    return { kind: 'carrier', vertex: best, point: h.point };
  }

  screenshotDataURL() { this.render(); return this.canvas.toDataURL('image/png'); }
}

function spheres(P, count, radius, colorFn) {
  const geo = new THREE.SphereGeometry(radius, 12, 8), mat = new THREE.MeshStandardMaterial({ roughness: 0.6, metalness: 0 });
  const inst = new THREE.InstancedMesh(geo, mat, count), m = new THREE.Matrix4(), col = new THREE.Color();
  for (let i = 0; i < count; i++) { m.makeTranslation(P[3 * i], P[3 * i + 1], P[3 * i + 2]); inst.setMatrixAt(i, m); colorFn(i, col); inst.setColorAt(i, col); }
  inst.instanceMatrix.needsUpdate = true; if (inst.instanceColor) inst.instanceColor.needsUpdate = true;
  return inst;
}
/** spherical interpolation of unit vectors; antipodal pairs go the long way round through a vector near +z */
function slerp(a, b, t) {
  const dot = Math.max(-1, Math.min(1, a.dot(b))), om = Math.acos(dot);
  if (om < 1e-6) return b.clone();
  if (om > Math.PI - 1e-3) {
    let m = new THREE.Vector3().crossVectors(a, Math.abs(a.z) < 0.9 ? Z_UP : new THREE.Vector3(0, 1, 0)).normalize();
    m = new THREE.Vector3().crossVectors(m, a).normalize();                   // ⟂ a, in the plane of a and (roughly) z
    return t < 0.5 ? slerp(a, m, 2 * t) : slerp(m, b, 2 * t - 1);
  }
  const s = Math.sin(om);
  return a.clone().multiplyScalar(Math.sin((1 - t) * om) / s).addScaledVector(b, Math.sin(t * om) / s).normalize();
}
function disposeObject(o) { o.traverse?.(x => { x.geometry?.dispose?.(); if (x.material) (Array.isArray(x.material) ? x.material : [x.material]).forEach(m => m.dispose?.()); }); }
const lerp3 = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
const nrm3 = (a) => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
const cross3 = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
