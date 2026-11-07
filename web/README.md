# Try-tri-weave — the remesh for weave mesh 

The  Swift–Hohenberg solver of `../localize_sh_app.py`, in JavaScript: the time stepping runs on
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
   ├─ index.html              the page (complete document)
   ├─ Tri-Weave-Lattice-Equations.html  the equations page, self-contained (left bar → 6 Explain; ← TryTriWeave goes back)
   ├─ artifact.html           the same page without the <html>/<head>/<body> wrapper — what gets published; made by build_artifact.py
   ├─ style.css               theme tokens (light/dark, WeaveSim idiom) and the layout
   ├─ manifest.json           list of carriers {name, file, n, omega_deg, chi, bbox, h, a_default}
   ├─ manifest_artifact.json  same list with .wasm file names (the artifact host only serves known types)
   ├─ carriers/*.lshc         prepared carriers (binary, ~0.6–1.5 MB each)
   ├─ js/
   │  ├─ app.js               state, UI bindings, frame loop, export
   │  ├─ carrier.js           .lshc loader, intrinsic Delaunay flips, Laplacian, heat solves, detection, log map, seed
   │  ├─ meshio.js            readers for STL (binary / ASCII), OBJ, PLY (ASCII / binary), OFF; the .lshc writer
   │  ├─ meshrepair.js        weld, duplicate faces, T-junction cracks, orientation, non-manifold edges / vertices, loose parts, planar caps, small holes
   │  ├─ remesh.js            isotropic remesher (Botsch–Kobbelt on a corner table, feature lines and fixed points, AABB-tree projection)
   │  ├─ import_pipeline.js   file → carrier data: the export_carrier.py pipeline in JavaScript
   │  ├─ import_worker.js     runs import_pipeline.js off the main thread
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
   ├─ test_import.mjs         the browser import on mesh files, beside the Python carriers
   ├─ test_predict.mjs        tracked read-outs against the growth prediction
   ├─ test_strands.mjs        strand tracing and colouring on a perfect flat lattice
   └─ test_gpu.html           CPU vs GPU solver on one seed (open it in the served page)
```

GitHub (`TryTriWeave/web/`, served by GitHub Pages) gets the page and what it loads: `index.html`, `style.css`,
`manifest.json`, `Tri-Weave-Lattice-Equations.html`, `carriers/`, `js/`, `vendor/`, plus `serve.py` and the two `.md`
files. `artifact.html`, `build_artifact.py`, `manifest_artifact.json`, `calibrate.js`, the `test_*` files and
`localsh-selfsigned.pem` (serve.py's private key) stay here.

## Running

Local: `python serve.py` from `web/` (any static server works, e.g. `python -m http.server 8080`), then
<http://localhost:8080/>; `index.html#v4Cone_sharp` opens on a given surface. `serve.py` also prints the LAN URL
for other devices on the same Wi-Fi and has `--https` (self-signed certificate) because WebGPU needs a secure
context on anything but localhost. Online: GitHub Pages, <https://uwenji.github.io/try-tri-weave/web/>
(repo `uwenji/try-tri-weave`, local clone `TryTriWeave`).

Browsers: WebGPU is in Chrome/Edge 113+ (desktop, Android), Safari 26+ (macOS 26 / iOS 26) and recent
Firefox on Windows; elsewhere the page falls back to the CPU solver and says so in the status pill. WebGL 2 is
needed for the view in any case.

## The panel

Layout (2026-10-18): status bar · left bar (solver strip + sections 1–4 + 5 Explain) · 3-D view · **Display** sidebar.
The left bar starts with the **solver strip** (pinned): the legend line (peaks · interior · singular, with a chevron that
folds the chart and a status mark), the **growth chart** (solid = tracked, dashed = the *ideal front* predicted from the
flat front speed — see *Prediction* below; singular has its own right axis; the x axis zooms to what is computed plus the
predicted cover ▴; click or drag on it to scrub the shown step), then `steps`, Run / +500 / Re-seed / Export, max steps and
per frame. Below it the sections fold with a click on their header; a folded header still shows its key value, and what is
folded is remembered in the browser. The **Display** sidebar (the old section 5, as in WeaveSim) sits beside the view: its
header folds it to a narrow vertical tab, ⇄ docks it on the other side of the view; both are remembered. On a phone it
stacks under the view.

**Status marks.** Read-outs no longer fill the panel: a coloured ⓘ carries them, and hovering (or tapping) it shows the
read-out in a card, live while it is open. Colour = how things stand: green fine, orange look at it, red a problem, grey
nothing yet. Growth (legend line): the singular count and the front against the ideal one (orange when it runs > 12 % off,
red when the pattern died out). File (zoo row): the carrier summary (orange when the import repaired cracks). Lattice edge
a (2 Density): the carrier resolution and the density diagram (green at 6–12 carrier edges per a). Heal: the last pass's report (green healed, orange pins held or nothing healed). Removed from the
panel: the computed-steps line, the carrier range line, the `on the model` chips (they are the `mesh` and `scale a` chips
of Display), the rotation-field / ĝ / stripes read-outs (the section summary and the status bar keep the essentials), the
valence number line, the cone editor, the view-cube hint and the lattice / selection read-outs (a click on a lath or a
vertex reports in the status bar; the valence chips in Display show the counts).

| section | controls | what happens |
|---|---|---|
| 1 File | zoo list · status mark · Open mesh… / drop anywhere · `up` axis of the file (z / y / −y / x / −x / −z) · Export lattice (the dual triangulation as OBJ text: copy, or download outside the artifact viewer) · import options (edges / bbox, max v, crease °, strip planar caps) · Remesh again · Download .lshc | loads a prepared carrier, or imports a mesh (.stl .obj .ply .off) in a worker: repair, isotropic remesh, intrinsic Delaunay — the same pipeline as `export_carrier.py` (see *Opening your own mesh*); a dropped `.lshc` loads as it is. The card shows vertices, h, min angle, flips, the repairs made |
| 2 Density | `lattice edge a` (log) with its status mark (the live diagram is in its card) · `≈ sites` · q0 · `correct q0` · `carrier ×1 ×2 ×4 ×8` ⓘ | re-seeds at the new spacing; the diagram draws a, the lath pitch p = √3⁄2·a, the lath width and the carrier mesh at the true ratio (one tick per carrier edge along a), with sites = area / (√3⁄2·a²) worked out and a gauge of carrier edges per a (6–12 green). The slider spans 6…14 carrier edges per a, so a denser lattice needs a finer carrier: the × chips subdivide the exported carrier in the browser (see *Density* below). `intrinsic mesh` draws the intrinsic Delaunay triangulation the solver uses (flipped edges traced as geodesics across the faces, in the accent colour); `scale a` draws the seed's first ring at the centre (n sites at distance a, spokes, circle ρ = a) |
| solver strip (pinned) | legend · fold chevron · status mark · growth chart · `steps` (scrub, up to 20 000) · Run / +500 / Re-seed · max steps · per frame | snapshots every ~20 steps make the slider scrub backwards |
| 3 Singularity | `seed` local / global rotation field · local: `valence at the centre` · `Use detected` · click on the surface · seed kind / radius (+ critical radius) / rotation · global: `size field ĝ` (cones from the curvature / uniform a) · Solve rotation field + seed | Ω = Σ K over interior vertices → charge round(Ω/60°) → valence 6 − charge; centre = maximum of the heat-smoothed curvature (flat surfaces: point farthest from the boundary); any change re-seeds. R_c = D / v at the current (r, b): a seed patch smaller than that shrinks. Global mode: see *Global mode* below |
| 4 Solve | (r, b) map · `r` `b` · front readout · `dt` ⓘ + dt graph · CG iterations | the field is kept when r, b, dt change (operator rebuilt). The (r, b) map is measured on a flat domain with the solver's own scheme, its edges by bisection (hatched: the seed dies; grey: pinned — the patch survives but does not grow, a band under 0.01 wide for b ≤ 1.1 that widens to 0.25 at b = 2; blue bands: the front advances, darker = faster, key in a per 100 steps; dots above r = 0: the zero state is unstable and patterns nucleate everywhere); drag the dot to set r and b, hover to read a point. The dt graph: steps to cover the surface against how much the time step slows the front |
| 4 · heal | `Heal 5-7-5` · `pin what will not hold` · `Clear pins` · status mark (the report) | takes the bound 5–7 pairs apart inside the solver (see *Healing bound pairs* below): 5-7-5 → one 5, 7-5-7 → one 7, wrong diagonals flipped; each kick is held for 60 steps and tested over 400 free steps, kept only if the defects around it went down; a kick that does not hold can be pinned (kept applied after every batch) |
| Display (sidebar) | theme auto/light/dark · chips field · mesh (intrinsic) · scale a · lattice · peaks · singular · centre · laths · **orient** (a thin flat arrow at every regular site, in the tangent plane, along the local wave vector of the pattern — perpendicular to one family of lattice edges, 30° off the lattice directions, so it never hides under a lattice line; a direction of the lattice is only defined mod 60°, so the arrows are chosen consistently by transport outward from the site nearest the centre — parallel inside a grain, with one unavoidable 60° cut around a disclination; colour = misorientation to the neighbours after transport over one edge — the lattice-line colour in a grain, amber → red along seams — or `angle θ`, θ mod 60° against the projected global x in cyclic colours; coherence |⟨e^{6iθ}⟩| in the readout) · **strain** (lattice edges as bars coloured cool–warm by ℓ / a − 1, all or one lath family A/B/C, ± range; hover the chip or the ⓘ for the diagram: an edge of family A is twice the A-lath segment beside it, so ℓ > a means longer segments between crossings and larger hexagonal openings) · one chip per singular valence · `lath width [a]` · lath colouring · over/under · `input normal` (along the normals / flipped) · carrier shows (u / Gaussian curvature / smoothed / field 5/7 / stripes, 3 phases / plain) with `colour opacity` toward the plain mesh colour (stripes: one slider each for red, yellow, blue) · `turntable [°/s]` | singular vertices as discs with the valence digit (Paul Tol *light* fills: v3 orange, v4 yellow, v5 pink, v7 blue, v8 cyan, v9 mint); read-out at the bottom: strands, knots, parity and colour conflicts |
| selection | click a lath → its path is highlighted (others dimmed, dark centreline), status bar: family / knots / over-under / open-closed; click a disc or peak → selected vertex with its lattice index in a larger disc (status bar: valence, interior or open link); Esc clears; Shift-click or `Pick centre on surface…` moves the seed centre | selections are remembered by carrier vertex ids, so they survive the re-extraction that happens while the lattice grows |
| growth strip | click / drag scrubs the step | peaks · interior · singular against the step; the step where "singular" leaves 1 is where the surface stops being ideal for that lattice |
| view cube (top right of the view) | click a face / edge / corner · drag the cube · ▲▼◀▶ · ⌂ home · ⤢ fit · ⟳ turntable · ▱ perspective/orthographic · keys H F T P | the cube mirrors the camera (z up: FRONT faces −y, RIGHT +x, TOP +z); a click turns the view to look along that facet, animated; the arrows roll the model by 90° (▲ shows the face below, ◀ brings the right face to the front); the turntable spins the view about the vertical axis at the speed set in Display and pauses while you drag; projection switches keep the apparent size |

## Prediction (growth strip, 2026-10-15)

`js/frontspeed_table.js` holds the front speed v(r, b) of the SH23 hexagon pattern into u = 0, measured on a flat strip
with the solver's own semi-implicit step (dt 0.1, q0 = 1; r −0.5…0.1 × b 0.5…2; `../experiments/frontspeed/`). A circular
patch grows as dR/dt = v − D/R with D ≈ 0.8 + 3.4 v (fitted to radial runs; R_c = D/v is the critical radius). On the
carrier, R(t) is turned into covered area with the geodesic distance from the centre (`Carrier.geodesicDistance`, fast
marching on the intrinsic mesh — the heat distance compresses near the rim), then peaks = area / (√3⁄2·(0.975 a)²),
interior = the part farther than 0.85 a behind the front and from the rim, and the singular floor = 1 for the seed +
round(|Ω_enclosed/60° − (6 − n)|), the charges a strain-free lattice would need.
Checked (time to half cover, tracked vs ideal): flat disk +10 %, ideal cone n=3 −12 %, v4Cone_sharp −21 %,
v3Cone_sharp −25 % — the gap grows with the mismatch between the cone angle and n·60°, which is the information.
`node test_predict.mjs <carrier> <steps> [r b dt] [r2,b2,step]` prints both side by side.

The conformal size field was taken out of the page on 2026-10-15 (sheets in `../shots/ui_2026-10-15/`);
`Carrier.sizeField / conformalOps` remain for `test_sizefield.mjs`. Every ⓘ shows a hover card (`initTips` in app.js:
the button's text, or a `<template>` named by `data-tip` — the strain diagram is one).

## How it works

**Carriers.** The zoo STLs are Rhino NURBS meshes (fans, slivers, no interior vertices). The zoo carriers are made in
Python (an opened mesh goes through the JavaScript port of the same pipeline, see *Opening your own mesh*): `export_carrier.py` runs the pipeline of `../carrier.py`
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
*tree* used to come out facing inward (the old +z rule), which put the laths inside the trunk; `input normal` in
Display flips the side by hand, and over/under swaps with it. Picking follows what is on screen: an opaque
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
    node test_import.mjs ../../model/zoo/goose.stl       # the browser import on a mesh file, beside carriers/goose.lshc
    node test_sizefield.mjs v4Cone_round 4000 4          # uniform a vs conformal size field, strain per family, orientation coherence
    node test_predict.mjs v3Cone_sharp 2400              # tracked peaks / interior / singular against the growth prediction
    node test_strands.mjs                                # strands on a perfect flat lattice: straight, 0 conflicts
    http://localhost:8080/test_gpu.html?c=v5Cone_sharp&steps=600   # CPU vs GPU timing and max |difference|

## Troubleshooting

* Blank view / "load failed": the page was opened from `file://`. Serve the folder (section Running).
* Pill says `CPU (JS)`: no WebGPU in this browser; everything works, 100–300 steps/s instead of thousands.
  In Chrome check `chrome://gpu` → WebGPU; in Safari it is on by default from Safari 26.
* Slow with WebGPU: raise `per frame` (Solve) — batches are capped by that number.
* A surface is missing from the list: open the mesh (`Open mesh…` or drop it), or add it to `manifest.json` with
  `python ../export_carrier.py <name>`.
* An import fails or looks wrong: the card of the `1 File` mark lists the repairs; a point cloud, a mesh without faces
  or a NURBS file (.3dm, .step) cannot be read — export a mesh (STL / OBJ / PLY) from the modeller.
* `export_carrier.py` fails on import: run it inside the `weavesim` env with `gpytoolbox` and `potpourri3d`
  installed; it imports `../carrier.py` and `../localize_sh_app.py` (for the synthetic cones).
* Publishing: copy the GitHub set (section Folder) into `TryTriWeave/web/`, commit and push; Pages serves it as is.
  The Claude artifact route still works: `python build_artifact.py`, then publish `artifact.html` with `style.css`,
  `js/`, `vendor/`, `manifest_artifact.json` as `manifest.json`, and each `carriers/X.lshc` as `carriers/X.wasm`.

## Global mode — no seed (2026-10-15)
Section 3 → `seed: global rotation field`. `js/rosy.js` solves the globally smoothest 6-fold direction field on the carrier
(Knöppel–Crane–Pinkall–Schröder 2013; intrinsic frames, connection Laplacian, inverse iteration) — its singularities are the
5s/7s the curvature demands (treeBig: exactly 30 sevens, six per fillet, twelve at the junction). `js/stripes.js` integrates the
three stripe phases along it as one rank-3 bundle (Knöppel et al. 2015; families permute around singularities, conjugated when
reversed); its branch points are the dislocations the geometry forces; u0 = ⅓ Σ Re ψ̂ is the hexagonal initial field everywhere.
Layer: `rotation field`. `carrier shows` → `field 5/7`: the field's singular faces as charges, red = 5, blue = 7, each a Gaussian of the geodesic distance (σ = 0.8 a, full colour at the face); → `stripes, 3 phases`: the seed's three stripe phases (before SH) as red / yellow / blue bands where each family's laths would run (band = lath width, phase ≡ π), one opacity per family; a band that stops is a branch point, two families stopping together make a 5–7; around a field 5/7 the colours swap across a seam (no consistent labelling exists there). Per-face phases come from `stripes.js · stripeFaces` in the worker. Headless: `node test_rosy.mjs treeBig`, `node test_global.mjs treeBig 2500 1 321.2`
(treeBig ×2 at 558 sites → {5:2, 6:469, 7:28, 8:2}, Σ(6−n) = −30). The two solves run in a Worker (`js/global_worker.js` →
`js/global_pipeline.js`, shared with `test_ghat.mjs`).

### Size field ĝ in global mode (2026-10-18)
Section 3 → `size field ĝ`. At ×2 and denser the stripe phases are forced to fork (branch points) wherever the field's
curl — the surface's curvature spread over the surface — leaves turning that no smooth phase can follow; SH keeps those forks
as 5–7 pairs and 5-5-7 clusters even where the lattice is nearly flat. The three steps: (1) count the forced forks
(`stripes.js · branchPoints`), (2) split the turning into what a size field can absorb and what topology fixes (6χ net
charge on a closed surface, the rim's turning on an open one), (3) put the turning into cones and let ĝ = e^{2φ} g be flat
between them: L φ = Σ_c Ω_c δ_c − K (`carrier.sizeField`, cones smeared over a/2 so φ has no log spikes). The rotation field,
the stripes and SH are then solved in ĝ (`conformalOps`: ℓ̂ = e^{(φi+φj)/2} ℓ, faces a steep φ would flatten are scaled
uniformly), q0 is constant in ĝ, and the lattice is regular there with edge a·e^{−φ} on the model.
`cones from the curvature` (`js/cones.js · placeCones`): the curvature is quantised to the lattice scale by error diffusion
along a tour of a-sized cells — the residual charge runs up cell by cell and every time it passes ±30° a cone of
round(charge/60°) (valence 5/4/3 or 7/8/9) is placed in the cell that contributed most; a closed surface ends with 6χ net;
five rounds of Lloyd relaxation then move each cone to the centroid of the curvature it stands for, and each cone is
smeared over a/2 along the curvature of its own sign (an apex keeps its cone, a dome spreads it), with φ normalised so the
site count equals that of uniform a. `uniform a`: no size field (the old behaviour).
`flat on developable parts` (2026-10-28, on by default): φ is solved with a Laplacian 30× stiffer ≥ 1 a from the curvature
(|K| ≥ 5° per lattice cell; `cones.js · developableWeights`, `Carrier.sizeField · weights`), so tubes, flat caps and strips
keep the edge a and the size change stays on the curved strips, where the 5s and 7s sit. Without it φ is harmonic there and
takes its level from the curved parts — along a tube it does not decay, so tubes shrink or taper (treeBig ×2: 0.95 a,
junctions up to 1.5 a) and a crease ring swings 0.7–1.6 a between its cones (cylinder). The field in ĝ is relaxed to unit
length like the model-metric one when `singularities on the curvature` is on.
Read-outs: the cones and their net charge, the model edge a·e^{−φ} (p5–p95), the share of ĝ edges longer than a/6 (refine when
it exceeds a few %), and the branch points of the ĝ stripes — the test of step 3: on the domes they go from 25 to 0.
`carrier shows` → `lattice edge a·e^{−φ}` draws the size field (red tighter, blue wider); `field 5/7` and `stripes` now also
work in local mode (computed for the display, the local seed stays). Headless: `node test_ghat.mjs v4Cone_round greedy 1 2500`
(mode greedy | field | uniform, level, steps, r, b) appends to `test_ghat_results.jsonl`.

**Cone editor (removed 2026-10-18).** An editor for the automatic cones (on/off, valence, radius, rotation with an
alignment weight) was tried and taken out again: too much hand work for what it gave. The worker and `global_pipeline.js`
still accept `edits` (`cones.js · applyConeEdits`, the constrained `smoothestRoSy`), the app sends none.

**Solver hand-over.** Switching to global mode (or changing the density) while the GPU solver is stepping used to freeze the
run: the old `GpuSolver` was disposed with a step in flight, its `mapAsync` rejected and the pump stopped the run at the
old step count. `GpuSolver.dispose` now waits for in-flight steps and the pump carries on with the replacement solver.
`cb-disp` (section 2): the solvers now correct q0 for the carrier's dispersion (edge = a; see `test_dispersion.mjs`, `test_calibrate.mjs`).

## Healing bound 5–7 pairs (2026-10-18)
At ×2 the 5s acquire bound dislocations — the beginning of grain-boundary scars (Bausch et al. 2003): a 5-7-5 where the
lattice could hold a single 5, a 7-5-7 where it could hold a single 7. Read as a mesh, each is one local operation away
from the single defect (`js/heal.js · findMotifs`, from the cyclic links of the dual triangulation):
- **insert** — 5-7-5 with two regular sites between the 5s in the 7's link: a pentagon (the 7, the two 5s and those two
  sites) missing its centre. The 7 holds both diagonals; a site at the centre makes it a 5 and the five corners 6s.
- **merge** — 7-5-7 with one regular site B between the 7s in the 5's link: a heptagon holding one site too many; merging
  the 5 with B leaves one 7.
- **flip** — 5-7-5 with one site between the 5s (quad with the diagonal on the 7), or a 7-7-5 triangle: an edge flip.
- **tri** — 5-5-7 triangles: no single flip (a dislocation beside a 5); since 2026-10-27 healed by merging the two 5s, and free
  5–7 pairs are glided (`glidePass`).
The operation is done in the field, not on the mesh (`healPass`): a bump of the peaks' median height (radius 0.42 a) is
written at the new site, the 7 is pushed 0.4 a outward (trough on its old place, bump on the new one) — without the push the
new peak and the 7 merge and slide back; merge: troughs on the two sites, bump between them; flip: the edge's ends moved
apart and the opposite corners together by 0.35 a. The kick is re-applied every 5 steps for 60 steps, then the field runs
400 steps free; the lattice is read again and the kick is kept when the non-6 peaks within 2.5 a went down with the same
charge and the total did not go up, else the field and step count are put back. Motifs are re-read after every attempt;
an insert tries the pentagon's centre first and the far side's hexagon centre (the 7 becomes the 5) second.
**Pins.** A 5-7-5 that re-forms whenever the field is free is one the physics wants — on Bunny the one case had a second
7 at 1.9 a (two dislocations, local charge 0). Its held state is still a clean single 5 (spokes 0.82–0.96 a), so with `pin
what will not hold` the kick is kept applied after every batch of steps (`applyPins` in the pump) and the 5 stays; only
inserts and merges whose held lattice is cleaner are pinned, never flips; `Clear pins` lets the field go back. Pins are
dropped with the field (re-seed, density, new surface).
Measured at ×2 with the uniform global seed, 3000 steps: Bunny (a 8.1, r −0.11, b 1.44) 67 → 59 non-6 (2 merges, 1 insert,
1 pinned; 8 s on the GPU); goose (a 7.08, r −0.24, b 1.70) 65 → 57 (1 flip, 3 inserts, nothing pinned; 6 s). CPU check:
`node test_heal.mjs Bunny 8.1 1 3000 -0.110 1.44` (the field after the growth is cached in `heal_<name>_<level>_<steps>.json`).
Also in this version: the pump yields through a MessageChannel in a hidden tab (Chrome throttles timers there to 1 Hz,
1/min after five minutes), and starts as soon as the target changes or a seed arrives instead of waiting for the tick.

## Opening your own mesh (2026-10-19)
`Open mesh…` (or a drop anywhere on the page) runs the whole `export_carrier.py` pipeline in a module worker, so the GitHub
Pages site takes a user's model without Python and without uploading it anywhere:
1. **read** (`meshio.js`) — STL binary / ASCII, OBJ (polygons fanned, negative indices), PLY ASCII / binary (any vertex
   properties, list faces), OFF. A `.lshc` is loaded directly. `up` turns the file's up axis onto +z with a proper rotation
   (y-up: (x, y, z) → (x, −z, y), the glTF / Blender convention, front kept); changing `up` while an opened model is shown
   rotates the carrier itself (`UP_ROT`, `upChange` in `import_pipeline.js`) — the remesh does not depend on the
   orientation — and the choice is stored in the carrier's meta (`options.up`) and in the browser for the next file.
2. **repair** (`meshrepair.js`) — weld (1e-6 of the bbox, hashed), drop degenerate and duplicate faces, close T-junction
   cracks (`repair_cracks`), orient consistently (BFS; a non-orientable mesh — a Möbius band — is cut along a short edge
   path and glued back after the remesh, see *Non-orientable carriers*), cut edges with more than two faces
   and split vertices whose corners form more than one fan, so the corner table is a manifold; drop loose parts smaller
   than 4 h²; strip planar caps (`strip_planar_caps`, not when the flat part is most of the model, e.g. a box); fill
   boundary loops shorter than 3 h (`fill_small_holes`).
3. **remesh** (`remesh.js`) — Botsch–Kobbelt at h = bbox/120 (`edges / bbox`), 10 iterations of split > 4⁄3 h, collapse
   < 4⁄5 h, valence flips (squared deviation from 6, 4 on the boundary, the angle-sum valence at fixed points),
   tangential relaxation and projection onto the input through an AABB tree (a vertex only takes a triangle facing its
   own side, so thin parts do not merge). The boundary and creases (dihedral > 40°) are feature lines: their vertices
   slide along the line and are projected onto it; line ends, junctions, corners (turning > 45° over h/2) and apices
   (|K| > 20°, the cone tips) are fixed. Two fixed points closer than h/2 merge. h is coarsened for the vertex cap
   (26 000) from the area (n ≈ 1.2 A/h²) before the remesh, and again if the result is still over. The orientation cut of
   a non-orientable mesh is a boundary here too, but its zigzag along the input edges makes no corners (its vertices
   slide along it) and only its two ends on the rim are fixed; after the remesh its two sides are re-glued.
4. **carrier** — centre, outward side (`orientSide`), angle defects, intrinsic Delaunay flips (`flipToDelaunay`), the
   same meta as the Python export (χ, ∫K, h, min angle, skinny %, notes). `Download .lshc` writes it in the carriers/
   format.

Against the Python carriers (same h, all 19 zoo STLs, `node test_import.mjs …`): χ identical everywhere, ∫K within 1°,
area within 0.03 %; 3–13 % fewer vertices (the mean edge lands on h); minimum angle 22–30° on the CAD meshes where
gpytoolbox leaves 1–12° slivers at the snapped apex and along the rim (v4Cone_sharp 22° vs 2.8°, saddle2 30.5° vs 0.9°,
developableSrf1 23.5° vs 1.5°), 11–14° on the scanned Bunny and the creased goose (Python 6.9° and 2.1°), slightly lower
on the smooth genus 2 / 3 (28° vs 33° / 29°). Time in Node: 0.5–1.6 s per zoo model, 3.3 s for tree (38 k vertices),
5.7 s for Armadillo (346 k triangles), 11 s for treeBig (290 MB ASCII STL).
The lattice does not tell the two carriers apart. Global SH on 10 zoo shapes, Python carrier vs browser carrier at the same a
(bench of `Tri-Weave-Lattice-Equations.html`, CPU, non-6 vertices): v4Cone_sharp 2 / 1, v3Cone_sharp 4 / 8, v3Cone_round 3 / 3,
v5Cone_sharp 1 / 1, v4Cone_round 9 / 9, saddle1 1 / 1, Tube 0 / 0, torus 24 / 18, genus 2 38 / 39, Bunny 45 / 40 — 127 / 120 in
all; the spread of the lattice edge (CV 0.05–0.11) is the same on both.

## Non-orientable carriers — a Möbius band (2026-10-30)
A Möbius band has no consistent winding. The import used to leave it cut open where the orientation BFS met itself, so
the carrier was a disc (χ = 1 instead of 0) with a free edge across the band, the two sides of that cut remeshed
independently (mobius_surface.stl: 106 of 374 cut vertices still coincided; min angle 3°, mobius.stl 0° with 4 % slivers).
SH saw a rim there: the pattern grew up to the cut from both sides without meeting — a seam line of 5s and 7s, and a 10–12
valent peak where the cut met the rim (the Voronoi cell of a rim peak wrapping round the cut end). SH itself needs no
orientation (cotan Laplacian and mass from lengths); what did:
* **import** (`meshrepair.js · orientFaces, cutPaths, zipCut`, `remesh.js · relaxRegion`) — the BFS conflicts are a zigzag
  through a row of faces or a staircase (mobius.stl: 117 long for a band 10 wide); single faces are turned over while that
  shortens the cut (→ 10.5), the disc is remeshed with the cut as a sliding feature line (ends on the rim fixed), then the
  two sides — both on the cut path, sampled independently — are merged into one sampling (samples of the two sides closer
  than h/2, mutual nearest, become one vertex at their mean arc length; every other sample is inserted into the opposite
  boundary edge, a fan in that edge's face) and welded, and the seam and two rings round it are smoothed tangentially back
  onto the input. The carrier is the real surface: χ = 0, one rim, faces wound alike except across the seam (min angle
  24–26°, 0 % slivers). The notes say `non-orientable: cut along … edges` and `re-glued the cut after the remesh`.
* **carrier** (`carrier.js`) — `seamEdges`: the edges whose two faces run the same way along them; normals per side at the
  seam; `flipToDelaunay` links such a pair and flips it like any edge (the winding of the flipped pair follows the first
  face, so the seam moves with the flip); `cutAlongSeam` makes the display copy (`carrier.cut`).
* **rings and frames** (`rosy.js · buildRings`, used by `intrinsicFrames` and `Carrier._rings`) — the walk round a vertex goes
  through the shared edges whatever the faces' winding, keeps the oriented rings unchanged where the faces agree, and at a
  seam vertex turns the way most faces run counter-clockwise about N; `sg = −1` marks the edges where two rings turn
  opposite ways, `eps` the faces a ring passes against their winding. Across such an edge the parallel transport is a
  reflection, t ↦ ρ − t: the smoothest RoSy couples z_j to e^{iNρ} z̄_i (a symmetric real 2×2 block, the rest of the
  connection Laplacian unchanged), its index is read per face in the face's own turning sense, `liftField` and the stripe
  coupling (`stripeCoupling`) transport by the reflection, and the local seed's log map flips its turning sense past the
  seam (`or`).
* **view** (`view.js`) — the carrier is drawn cut open along the seam (one copy of each seam vertex per side, colours mapped
  through `cut.map`), so both sides shade smoothly; the lattice, peaks, valence discs, cones, centre, RoSy strokes, strain
  bars and orientation arrows are drawn on both sides (a one-sided surface has no outside to choose). `heal.js · burgers`
  compares Burgers vectors with the normals on one side.
On every orientable carrier all of this reduces to the old code bit for bit (rings, frames, flips, normals, log map, RoSy,
stripes, lifts on cylinder, Tube, saddle1, genus2 and the reversed treeBig; imports of saddle1 and v3Cone_sharp identical).
Ji's files, a = bbox/13, SH (−0.08, 1.2), 3000+ steps, singular peaks (old cut carrier → re-glued):
mobius.stl global seed 1 (a 12 at the cut end, and 124 peaks for 72 sites) → 0 (no field singularity, no branch point;
×2: 0 of 115 interior sites); local seed (v7 at the curvature maximum) 2 → 1–2. mobius_surface.stl (a closed-looking bulb
with a throat, ∫K ≈ 0 but curved everywhere) global 20 → 17 (after Heal pairs 16 → 11, no 8), local 20 → 19 (14 → 15); the
remaining 5s and 7s sit on the bulb and the throat, not on the seam. Sheets: `localize/experiments/mobius/`.
Not done yet: the laths' over/under is still referred to the normals, which flip across the seam (a parity break there in
the display; whether a kagome weave closes up on a Möbius band is its own question).
