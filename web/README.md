# LocalSH Weave Lab — the browser version

The localised Swift–Hohenberg solver of `../localize_sh_app.py`, in JavaScript: the time stepping runs on
**WebGPU** (plain-JS fallback), three.js draws the view, and a new **laths** layer shows the Kagome strands of the
extracted lattice as ribbons with adjustable width, over/under offset and family colouring. One seeded
singularity per surface, the same equation and the same three sliders (valence · steps · density) as the
Python app. `HOW_TO_START.md` is the short version of this file.

    u_t = r u − (1 + Δ/q0²)² u + b u² − u³        q0 = 4π/(√3 a),   r < 0 : the pattern only grows from the seed

## Folder

```
localize/
├─ export_carrier.py          Python: zoo .stl → web/carriers/*.lshc + web/manifest.json (remesh + intrinsic Delaunay)
└─ web/
   ├─ index.html              the page (complete document, for local use)
   ├─ artifact.html           the same page without the <html>/<head>/<body> wrapper — what gets published; made by build_artifact.py
   ├─ style.css               theme tokens (light/dark, WeaveSim idiom) and the layout
   ├─ manifest.json           list of carriers {name, file, n, omega_deg, chi, bbox, h, a_default}
   ├─ manifest_artifact.json  same list with .wasm file names (the artifact host only serves known types)
   ├─ carriers/*.lshc         prepared carriers (binary, ~0.6–1.5 MB each)
   ├─ js/
   │  ├─ app.js               state, UI bindings, frame loop, export
   │  ├─ carrier.js           .lshc loader, STL/weld/orient/Loop/Delaunay for dropped files, Laplacian, heat solves, detection, log map, seed
   │  ├─ sparse.js            CSR, SpMV, Jacobi-PCG, binary heap, Dijkstra
   │  ├─ solver_gpu.js        SH23 stepping in WGSL compute (Jacobi-PCG, scalars on the GPU)
   │  ├─ solver_cpu.js        the same stepping in typed arrays (fallback and reference)
   │  ├─ lattice.js           peaks → Voronoi dual → valences (closed-link test); strand tracing, parity, 3-colouring
   │  ├─ view.js              three.js scene: field, wire, lattice, peaks, singular, centre, laths, picking; camera moves, turntable, projection
   │  ├─ navcube.js           the CAD-style view cube (own small canvas): faces, edge strips, corners, drag-to-orbit, xyz triad
   │  ├─ frontspeed.js        front state at (r, b, dt), critical radius, growth prediction for the strip
   │  ├─ frontspeed_table.js  measured flat front speed v(r, b) and its dt dependence
   │  └─ colormap.js          viridis / cool–warm LUTs, charge and family colours
   ├─ vendor/                 three.js r170 (module build) + OrbitControls — no npm, no bundler
   ├─ serve.py                http / https server that also prints the LAN URL for other devices
   ├─ build_artifact.py       index.html → artifact.html
   ├─ test_node.mjs           carrier → detect → seed → CPU solver → lattice → strands, from the terminal
   ├─ test_predict.mjs        tracked read-outs against the growth prediction
   ├─ test_strands.mjs        strand tracing and colouring on a perfect flat lattice
   └─ test_gpu.html           CPU vs GPU solver on one seed (open it in the served page)
```

## Running

Local: `python serve.py` from `web/` (any static server works, e.g. `python -m http.server 8080`), then
<http://localhost:8080/>; `index.html#v4Cone_sharp` opens on a given surface. `serve.py` also prints the LAN URL
for other devices on the same Wi-Fi and has `--https` (self-signed certificate) because WebGPU needs a secure
context on anything but localhost. Online: the published artifact (private until shared).

Browsers: WebGPU is in Chrome/Edge 113+ (desktop, Android), Safari 26+ (macOS 26 / iOS 26) and recent
Firefox on Windows; elsewhere the page falls back to the CPU solver and says so in the status pill. WebGL 2 is
needed for the view in any case.

## The panel

The left bar starts with the **growth strip** (pinned, never folds: peaks · interior · singular against the step, solid =
tracked, dashed = the *ideal front* predicted from the flat front speed — see *Prediction* below; singular has its own right
axis; the x axis zooms to what is computed plus the predicted cover ▴; click or drag on it to scrub the shown step, the
shaded range is what has been computed; the line under it says how far the tracked front runs behind the ideal one). Below it the five sections fold
with a click on their header; a folded header still shows its key value (surface · a and sites · seed valence and centre
· steps · active layers), and what is folded is remembered in the browser.

| section | controls | what happens |
|---|---|---|
| 1 Surface | zoo list · Open STL… / drop · readouts | loads a prepared carrier (or builds one from an STL in the browser: weld, orient, Loop-subdivide to h, Delaunay flips — no remesher), shows vertices, h, min angle, flips, notes |
| 2 Density | `lattice edge a` (log) ⓘ · `≈ sites` · q0 · live diagram · `carrier ×1 ×2 ×4 ×8` ⓘ · on the model: `intrinsic mesh`, `scale a at the centre` | re-seeds at the new spacing; the diagram draws a, the lath pitch p = √3⁄2·a, the lath width and the carrier mesh at the true ratio (one tick per carrier edge along a), with sites = area / (√3⁄2·a²) worked out and a gauge of carrier edges per a (6–12 green). The slider spans 6…14 carrier edges per a, so a denser lattice needs a finer carrier: the × chips subdivide the exported carrier in the browser (see *Density* below). `intrinsic mesh` draws the intrinsic Delaunay triangulation the solver uses (flipped edges traced as geodesics across the faces, in the accent colour); `scale a` draws the seed's first ring at the centre (n sites at distance a, spokes, circle ρ = a) |
| 3 Singularity | number line (valence n · charge 6 − n · cone angle n·60°, detected Ω as a marker) · `valence at the centre` · `Use detected` · click on the surface · seed kind / radius (+ critical radius) / rotation | Ω = Σ K over interior vertices → charge round(Ω/60°) → valence 6 − charge; centre = maximum of the heat-smoothed curvature (flat surfaces: point farthest from the boundary); any change re-seeds. R_c = D / v at the current (r, b): a seed patch smaller than that shrinks |
| 4 Solve | `steps` (scrub, up to 20 000) · Run / +500 / Re-seed / Export lattice · (r, b) map · `r` `b` · front readout · `dt` ⓘ + dt graph · max steps · per frame · CG iterations | the field is kept when r, b, dt change (operator rebuilt); snapshots every ~20 steps make the slider scrub backwards; export = the dual triangulation as OBJ text (copy; download only works outside the artifact viewer). The (r, b) map is the measured front speed (hatched: the seed dies; tan: pinned; blue: advances, iso-lines 0.25/0.5/0.75/1 a per 100 steps; above r = 0 the zero state is unstable); drag the dot to set r and b, the diamond to set the anneal point. The dt graph: steps to cover the surface against how much the time step slows the front |
| 4 · anneal | `anneal after growth` · r₂ b₂ · switch when the front stops / at step N | two-stage schedule: grow at (r, b), then, once the peak count has stopped growing (3 read-outs) or at step N, continue with the same field at (r₂, b₂) nearer the hexagon saddle-node (default −0.027, 0.69 ≡ a k⁶-stiffened operator at the original parameters): dislocation dipoles annihilate 4–8× faster, the front is 3–4× slower. Removes isolated extra 5–7 pairs; holonomy seams (grain boundaries between fronts that met) stay. The strip tints the anneal stage; the r/b sliders drive the anneal stage while it runs; Re-seed restores (r, b) |
| 5 View & readouts | theme auto/light/dark · chips field · mesh (intrinsic) · scale a · lattice · peaks · singular · centre · laths · **orient** (a thin flat arrow at every regular site, in the tangent plane, along the local wave vector of the pattern — perpendicular to one family of lattice edges, 30° off the lattice directions, so it never hides under a lattice line; a direction of the lattice is only defined mod 60°, so the arrows are chosen consistently by transport outward from the site nearest the centre — parallel inside a grain, with one unavoidable 60° cut around a disclination; colour = misorientation to the neighbours after transport over one edge — the lattice-line colour in a grain, amber → red along seams — or `angle θ`, θ mod 60° against the projected global x in cyclic colours; coherence |⟨e^{6iθ}⟩| in the readout) · **strain** (lattice edges as bars coloured cool–warm by ℓ / a − 1, all or one lath family A/B/C, ± range; hover the chip or the ⓘ for the diagram: an edge of family A is twice the A-lath segment beside it, so ℓ > a means longer segments between crossings and larger hexagonal openings) · one chip per singular valence · `lath width [a]` · lath colouring · over/under · `weave side` (along the normals / flipped) · carrier shows (u / Gaussian curvature / smoothed / plain) · `turntable [°/s]` | singular vertices as discs with the valence digit (Paul Tol *light* fills: v3 orange, v4 yellow, v5 pink, v7 blue, v8 cyan, v9 mint); valence histogram of the *interior* peaks, Σ(6−n) vs Ω/60°, strands, knots, parity and colour conflicts |
| selection | click a lath → its path is highlighted (others dimmed, dark centreline), readout family / knots / over-under / open-closed; click a disc or peak → selected vertex with its lattice index in a larger disc; Esc clears; Shift-click or `Pick centre on surface…` moves the seed centre | selections are remembered by carrier vertex ids, so they survive the re-extraction that happens while the lattice grows |
| growth strip | click / drag scrubs the step | peaks · interior · singular against the step; the step where "singular" leaves 1 is where the surface stops being ideal for that lattice |
| view cube (top right of the view) | click a face / edge / corner · drag the cube · ▲▼◀▶ · ⌂ home · ⤢ fit · ⟳ turntable · ▱ perspective/orthographic · keys H F T P | the cube mirrors the camera (z up: FRONT faces −y, RIGHT +x, TOP +z); a click turns the view to look along that facet, animated; the arrows roll the model by 90° (▲ shows the face below, ◀ brings the right face to the front); the turntable spins the view about the vertical axis at the speed set in section 5 and pauses while you drag; projection switches keep the apparent size |

## Prediction (growth strip, 2026-10-15)

`js/frontspeed_table.js` holds the front speed v(r, b) of the SH23 hexagon pattern into u = 0, measured on a flat strip
with the solver's own semi-implicit step (dt 0.1, q0 = 1; r −0.5…0.1 × b 0.5…2; `../experiments/frontspeed/`). A circular
patch grows as dR/dt = v − D/R with D ≈ 0.8 + 3.4 v (fitted to radial runs; R_c = D/v is the critical radius). On the
carrier, R(t) is turned into covered area with the geodesic distance from the centre (`Carrier.geodesicDistance`, fast
marching on the intrinsic mesh — the heat distance compresses near the rim), then peaks = area / (√3⁄2·(0.975 a)²),
interior = the part farther than 0.85 a behind the front and from the rim, and the singular floor = 1 for the seed +
round(|Ω_enclosed/60° − (6 − n)|), the charges a strain-free lattice would need. The anneal stage switches to v(r₂, b₂).
Checked (time to half cover, tracked vs ideal): flat disk +10 %, ideal cone n=3 −12 %, v4Cone_sharp −21 %,
v3Cone_sharp −25 % — the gap grows with the mismatch between the cone angle and n·60°, which is the information.
`node test_predict.mjs <carrier> <steps> [r b dt] [r2,b2,step]` prints both side by side.

The conformal size field was taken out of the page on 2026-10-15 (sheets in `../shots/ui_2026-10-15/`);
`Carrier.sizeField / conformalOps` remain for `test_sizefield.mjs`. Every ⓘ shows a hover card (`initTips` in app.js:
the button's text, or a `<template>` named by `data-tip` — the strain diagram is one).

## How it works

**Carriers.** The zoo STLs are Rhino NURBS meshes (fans, slivers, no interior vertices). Remeshing has no
JavaScript implementation, so it stays in Python: `export_carrier.py` runs the pipeline of `../carrier.py`
(strip the flat cap of closed cones, isotropic remesh to h = bbox/120 with the boundary turned into a crease and
the apex snapped back, intrinsic Delaunay flips) and writes `.lshc`: vertices, extrinsic faces (for rendering and
picking), intrinsic faces and lengths (for the operators), angle defects, boundary flags. The page rebuilds the
cotan Laplacian and lumped mass from the intrinsic data; all interior weights are ≥ 0.

**Density.** `a` is the lattice edge in model units. *Sites* are the lattice vertices, ≈ area / (√3⁄2 · a²); a
Kagome weave has three strand crossings per site, its strand pitch (centre to centre of parallel laths) is √3⁄2 · a,
and the lath width is the slider × a — so for a material of width w woven at pitch/width ratio ρ, a = ρ · w / 0.866.
The carrier must resolve the lattice: at least 6 carrier edges per a for the fourth-order operator, and the
time-step system's condition number grows like (a/h)² (κ ≈ 13 at 8 edges per a, ≈ 46 at 14), which is why the
slider stops at 14 edges per a and why the CG iteration count follows (a/h)² automatically (12 at 9 edges per a;
type a value in *CG iterations* to override). Exported carriers have h = bbox/120, so their a ranges from bbox/20 to
bbox/8.6 whatever the model's scale — scaling the model in Rhino changes nothing. Two ways to go denser:
`carrier ×2 / ×4 / ×8` subdivides the carrier in the browser by **midpoint subdivision** (each triangle → 4 similar
ones: the surface is identical, angle defects unchanged, intrinsic Delaunay preserved, base vertex ids kept so the
centre and selections survive; 4× / 16× / 64× the vertices, a down to bbox/40 / 80 / 160); detection, curvature
smoothing and the seed chart stay on the base carrier and are interpolated. Or export finer:
`python ../export_carrier.py tree --edges 240 --max-verts 120000` remeshes at h = bbox/240. The GPU handles ×4
(100–250 k vertices) comfortably; ×8 (0.5–1 M) steps at tens of steps/s and the read-outs take seconds (they run
less often: every 100 × n/12000 steps). Loop or Butterfly subdivision would be wrong here: Loop shrinks and
rounds the sharp apex and moves the curvature, Butterfly ripples around irregular vertices — resolution is all
we need from the refinement, not smoothing (the remeshed carrier is already smooth at bbox/120).

**Detection and seed.** Same as the Python app. The chart around the centre is ρ = heat-method distance
(t = 4h², 0.2 % error) and θ from a discrete exponential map propagated upwind on ρ with the intrinsic corner
angles rescaled to 2π at each vertex. At a cone apex this is exactly the rescaled chart that the folding
ζ = ρ e^{iθ n/6} expects; a tangent-plane projection is not, because the normal of a remeshed sharp tip is tilted
by ~20°, and that alone produced dislocations on the ideal cone. Seed: three plane waves in ζ (or the Gaussian
spots of the base script) under exp(−(ρ/2a)⁴).

**Time step.** Semi-implicit weak form `[(1 − dt r) M + dt S M⁻¹ S] u⁺ = M (u + dt (b u² − u³))`,
`S = M − L/q0²`, solved by matrix-free Jacobi-PCG, 12 iterations warm-started from the previous step
(relative residual ~1e-6, κ(A) ≈ 13). On the GPU one step is three kernels plus four per CG iteration
(≈ 50 dispatches); α and β never leave the GPU — every workgroup re-reduces the small array of per-workgroup
partial sums. The CPU solver is the same algorithm in float64 and is the reference (`test_gpu.html` shows the
two agree to 1e-6). Batches are sized to ~120 ms of solver time and yield to the renderer between batches.

**Lattice.** Local maxima of u above 0.4·max → multi-source Dijkstra Voronoi on the intrinsic graph → dual
triangles where three cells meet, excluding void (smoothed rms of u below 0.15·max) and boundary vertices → a
peak is *interior* when its dual link is a closed cycle, so the growth front and the rim never contribute fake
valences.

**Side.** Open surfaces get their global orientation from a rule shared by `carrier.py` and `carrier.js` (applied
again on load, so older exports agree): closed → outward; open → the convex side when the total mean curvature is
clearly signed (tubes, branching tubes, domes face outward), else +z (saddles, sheets). A branching tube like the
*tree* used to come out facing inward (the old +z rule), which put the laths inside the trunk; `weave side` in
section 5 flips the side by hand, and over/under swaps with it. Picking follows what is on screen: an opaque
field hides what is behind it, a hidden or wire-only carrier does not.

**Camera.** OrbitControls with +z up, so the view never rolls. Canonical views are reached by animating the
camera position on the sphere around the target (slerp, 450 ms); top and bottom views keep a 0.1° tilt toward the
front (or toward the azimuth you came from with the arrows) because the azimuth is undefined at the pole. The
orthographic camera shares the target and direction; switching keeps the visible half-height at the target
(`zoom = 1 / (dist · tan(fov/2))`). The turntable is OrbitControls' `autoRotate` fed with the real frame time, so the
speed in °/s holds even when the CPU solver makes frames irregular.

**Colours.** Lath families use Paul Tol's high-contrast scheme (#004488 deep blue, #DDAA33 warm gold, #BB5566 muted
crimson; lifted 18 % toward white on the dark canvas), strands "by strand" use his muted scheme, singular discs
his light scheme so the dark digit stays readable; the growth chart uses the same three high-contrast colours.

**Laths.** Each lath passes through lattice edge midpoints; in every dual triangle it takes the mid-segment that
wraps the vertex it did not wrap in the previous triangle (a straight line on a regular lattice,
`test_strands.mjs` checks it). Width = slider × a; over/under alternates along each strand and a parity BFS makes
it consistent at the crossings (odd closed strands and frustrated cycles are counted as parity conflicts).
Families: strands crossing at a knot must differ; greedy + min-conflict sweeps. On a v3 cone the colouring is
perfect (rotation by 180° maps each family to itself); around a v4 or v5 the three families permute, so a cut
of conflicts is unavoidable — the panel shows the count rather than hiding it.

## Reference results (CPU solver, a = default, 3000–4000 steps, r = −0.08, b = 1.2, dt = 0.1)

| surface | seed | interior valences |
|---|---|---|
| ideal cone n=3 | v3 | {3:1, 6:163} at a = 2.0 · {3:1, 6:121} at 2.3 · {3:1, 6:97} at 2.55 |
| v4Cone_sharp | v4 | {4:1, 6:77} |
| v5Cone_sharp | v5 | {5:1, 6:70} |
| v4Cone_round (dome 192°) | v4, uniform a | {5:3, 6:88, 7:1} at a = 2.31 |
| v4Cone_round (dome 192°) | v4, **size field** | {4:1, 6:80} (node test_sizefield.mjs) |

## Tests

    node test_node.mjs v5Cone_sharp 3000 [valence] [a]   # whole pipeline from the terminal (Node ≥ 20)
    node test_sizefield.mjs v4Cone_round 4000 4          # uniform a vs conformal size field, strain per family, orientation coherence
    node test_predict.mjs v3Cone_sharp 2400              # tracked peaks / interior / singular against the growth prediction
    node test_strands.mjs                                # strands on a perfect flat lattice: straight, 0 conflicts
    http://localhost:8080/test_gpu.html?c=v5Cone_sharp&steps=600   # CPU vs GPU timing and max |difference|

## Troubleshooting

* Blank view / "load failed": the page was opened from `file://`. Serve the folder (section Running).
* Pill says `CPU (JS)`: no WebGPU in this browser; everything works, 100–300 steps/s instead of thousands.
  In Chrome check `chrome://gpu` → WebGPU; in Safari it is on by default from Safari 26.
* Slow with WebGPU: raise `per frame` (Solve) — batches are capped by that number.
* A surface is missing from the list: it is not in `manifest.json`; run `python ../export_carrier.py <name>`.
* `export_carrier.py` fails on import: run it inside the `weavesim` env with `gpytoolbox` and `potpourri3d`
  installed; it imports `../carrier.py` and `../localize_sh_app.py` (for the synthetic cones).
* Republishing: `python build_artifact.py`, then publish `artifact.html` with `style.css`, `js/`, `vendor/`,
  `manifest_artifact.json` as `manifest.json`, and each `carriers/X.lshc` as `carriers/X.wasm`.
