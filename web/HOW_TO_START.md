# How to start — TryTriWeave (web)

Three ways, from quickest to most flexible.

## 1 · Online, nothing to install

Open <https://uwenji.github.io/try-tri-weave/web/> (GitHub Pages) in Chrome, Edge or Safari. The pill at the top right should read `WebGPU · …`; if it reads `CPU (JS)` the browser has no
WebGPU and the plain-JavaScript solver runs instead (about 5–10× slower, same results).

## 2 · On your Mac, from this folder

```bash
conda activate weavesim                 # any Python 3 works, this one is already there
cd ".../TopOpiTriWeave/localize/web"
python serve.py                         # or: python -m http.server 8080
```

Then open <http://localhost:8080/> in Chrome (or Safari 26+). Opening `index.html` by double-click does
**not** work: the page fetches the carriers and its own modules, and browsers refuse that from `file://`.
`serve.py` also sends no-cache headers, so edited modules and re-exported carriers show up on a plain reload.

Stop the server with Ctrl-C. Another port: `python serve.py --port 8090`.

## 2b · Other devices on the same Wi-Fi (iPad, phone, a colleague's laptop)

`serve.py` prints a second URL of the form `http://192.168.x.x:8080/` — open that on the other device
(macOS may ask once to allow incoming connections for Python; allow). Over plain http the other device gets the
CPU solver, because browsers only enable WebGPU on "secure" pages (localhost or https). For WebGPU on the other
device either open the GitHub Pages site (it is https) or run

```bash
python serve.py --https             # https://192.168.x.x:8443/, self-signed certificate made once with openssl
```

and accept the certificate warning on the device the first time.

## 3 · First five minutes in the app

1. **File** — pick `v4Cone_sharp` from the zoo list. The carrier loads, the singularity is detected
   (Ω = ∫K → suggested valence and centre, blue marker), the seed is planted and 3000 steps start running.
2. **Density** — drag `lattice edge a` (model units; strand pitch = 0.866 a, lath width = slider × a). The pattern
   re-seeds at the new spacing. The slider stops at 6 carrier edges per a: for a denser weave click `carrier ×2 / ×4 /
   ×8` (subdivides the carrier in the browser; ×4 and ×8 want the GPU) or export a finer carrier
   (`python ../export_carrier.py name --edges 240 --max-verts 120000`).
3. **Singularity** — change the `valence at the centre` slider (3–9) to plant a different singularity;
   Shift-click the surface (or `Pick centre on surface…`, then click) to move the centre; `Use detected` goes back.
4. **Solve** — the solver controls sit at the top of the left bar, right under the growth strip, so they are
   always in reach: drag `steps` (up to 20 000) left to scrub back in time, right to compute further; `Run` keeps
   going, `+500` adds a batch, `Re-seed` (`Export lattice` is in `1 File`); `max steps` and `per frame` below them. Section 4
   holds the pattern parameters: `r`, `b`, `dt` change the equation without resetting the field.
5. **Display** (the sidebar beside the view) — `theme` auto / light / dark. Chips toggle layers: `field · scale a ·
   lattice · singular · centre · laths · orient · strain · bounding box`; the `singular` row has one chip per valence
   present (v3, v5, v7 …, with counts) so each kind of singularity can be shown or hidden on its own. Singular vertices
   are drawn as discs with the valence digit inside. Turn **laths** on and move `lath width [a]` and `thickness [a]`;
   set `carrier shows` to *plain* to see the weave alone.
   `input normal` flips the side of the surface that carries the laths and marks (a branching tube has no side that
   is outside everywhere). The **growth chart** at the top of the left bar plots peaks · interior · singular against
   the step; the chevron folds it, and clicking or dragging on it scrubs the shown step; the solver controls sit right
   under it. The five sections below fold; `6 Explain` opens the equations page.
6. **Selecting** — click a lath to follow its whole path (the rest dims, a dark line runs along it, the status bar
   gives family, knot count, over/under and open/closed); click a singular disc or a peak to select that vertex
   (a larger disc shows its lattice index, which is its line in the OBJ export). Esc clears. Moving the centre is
   now Shift-click, or `Pick centre on surface…` then click, so an ordinary click never re-seeds by accident.
7. **Navigating** — the **view cube** at the top right of the view turns with the model: click a face, an edge or
   a corner to look from there (FRONT is −y, RIGHT is +x, TOP is +z), drag the cube to orbit, ▲▼◀▶ roll the
   model by 90°, ⌂ is the home view, ⤢ fits the surface, ⟳ starts the **turntable** (speed in Display,
   pauses while you drag), the last button switches perspective / orthographic. Keys: `H` home, `F` fit,
   `T` turntable, `P` projection.

What "ideal" looks like: exactly one singular peak and Σ(6−n) = round(Ω/60°). Try the three
`ideal cone` entries first, then your own cones.

## Global seed with the size field ĝ (2026-10-18)

* Section 3 → `global rotation field` needs no centre. `size field ĝ` → `cones from the curvature` (default) puts the surface's
  curvature into a few 5s and 7s and solves the field, the stripes and SH in a metric that is flat between them, so the lattice
  comes out regular with those cones as its only defects where the surface is resolved (v4Cone_round ×2: 3 fives, no pairs,
  with the default r, b). The read-out lists the cones, the model edge a·e^{−φ} they imply and the branch points left (0 is the
  goal). `carrier shows → lattice edge a·e^{−φ}` draws the size field; `stripes` and `field 5/7` also work in local mode now.
* `flat on developable parts` (on by default, 2026-10-28): the size field leaves tubes, flat caps and strips at the edge a and
  changes the edge only on the curved strips; untick it for the free size field (tubes then shrink or taper).
* Features thinner than a lattice edge (a sharp beak, a leg of two sites) still produce clusters — refine (`×2`) or lower the
  density there.
* **Heal 5-7-5** (now `Heal pairs`, in section 5 since 2026-11-20) — after the surface is covered, click it: the bound pairs (5-7-5, 7-5-7, wrong diagonals) are
  taken apart one by one inside the solver and kept only where the lattice holds the single defect by itself; the read-out
  says what was found, healed and pinned. `pin what will not hold` keeps the single 5 where the field would re-form the pair
  (held after every batch of steps; `Clear pins` releases). Click again after more steps: new pairs may have
  formed. Bunny ×2: 67 → 59, goose ×2: 65 → 57 in a few seconds on the GPU.

## Edit lattice: heal, freeze, smooth, edit (2026-11-20, 2026-11-24)

Section `5 Edit lattice` holds what changes a finished lattice; `6 Explain` is the old section 5. Connectivity is edited
in the field (`Heal pairs`) or on a lattice taken out of the field (`Freeze lattice`, then the click tools); positions
on the frozen lattice (`smooth`). Since 24 Nov the panel carries controls and one line of numbers: no report, no log,
no help texts. What a click did is one line in the status bar at the top; hover a control for a sentence about it.

* **Heal pairs** — in the field, as before: the bound motifs, then `glide` (free 5–7 pairs) and `lone 5 · 7` (to their
  curvature, experimental) if ticked; `pin` holds what will not stay. `Motifs`, `Glide`, `Move 5 · 7` run one pass
  alone. The result is one line in the status bar (`heal: 3 of 4 motifs healed · …`); every kick is in the browser
  console (`[heal] …`).
* **Freeze lattice** — takes the lattice of the shown step out of the field. Every site goes from the carrier vertex of
  its peak to the peak itself, between the vertices, and stays a point on the surface. While a lattice is frozen the
  solver is locked (`Run`, `+500`, `Re-seed`, the steps, sections 2–4, the heal row; the headers of 2–4 say `locked`).
  `Unfreeze` goes back to the read-out of the field. Edits made on the frozen lattice are not in the field: with edits,
  `Unfreeze` asks twice (`Unfreeze: drop 3 edits?`). Opening another surface or re-seeding discards the frozen lattice.
* **smooth** (0–100 sweeps, 20) — moves the frozen sites on the surface until the laths run as straight as the lattice
  allows (rhombi unfolded, so bending with the surface does not count). `lengths` (on) holds every edge at its length
  by soft springs, `hold rim` and `hold 5 · 7` keep those sites where they are.
* **click: heal · to rim · flip · collapse · split · patch** — one of them on, a click on the canvas edits the frozen lattice
  (the same chip again, or Esc: off). After every edit the lattice is smoothed again with the sweeps set.
  * `heal` — click a 5 or a 7: its cluster (the 5s and 7s that are neighbours or have one site between them) goes down
    to its net charge: 5-7-5 → 5, 7-5-7 → 7, 5-5-7 → 5, two opposite 5–7 pairs side by side → nothing. The sequence is
    found by search (up to three flips, collapses and splits; among those that reach the same valences, the one with
    the best triangles). A lone 5–7 pair is a dislocation and stays; with another 5 or 7 within four edges it is
    walked to it where that takes only flips and collapses. Next to the rim the cluster is taken out altogether where
    a short sequence does it: rim sites have no valence to keep (`Heal all` leaves the rim alone).
  * `to rim` — click a 5 or a 7: it is walked to the nearest rim and taken out there, in one step of Undo. With a
    neighbour of the other kind it goes as a 5–7 pair: a flip glides the pair one site along its row (no site added
    or removed), a collapse or a split + flip takes it to the next row; each step leaves Σ (valence − 6)² as it is and
    brings the pair closer to the rim, a glide first if one does. At the rim the naked edge behind the 5 is split and
    the one behind the 7 collapsed. Alone, a 5 or 7 does not walk (moving it one site leaves a pair behind): it hops —
    one flip moves it a site outwards and leaves a pair, the pair is walked out, and so on — so a lone one far in costs
    many edits and changes the number of sites (nine rows in: about 85 edits). After every step the sites round it are
    relaxed, so long walks leave no folded triangle. It gives up after 1.5 s and keeps the furthest it got; on a closed
    surface there is no rim to walk to.
  * `flip` — click an edge: it turns to the other diagonal. Its two ends lose one neighbour, the two opposite sites gain one.
  * `collapse` — click an edge: its two sites become one (valence: the sum − 4); the two opposite sites lose one.
    On a naked edge of the rim: the site behind it loses one — a 7 there becomes a 6 — and the rim is one site shorter.
  * `split` — click an edge: a new site in its middle (valence 4); the two opposite sites gain one.
    On a naked edge: the site behind it gains one — a 5 there becomes a 6 — and the rim has one site more.
  * `patch` — click a site of the rim: its fan is completed to six triangles (one triangle, or up to four new sites).
* **Heal all** — every cluster that holds more than its net charge, one after the other. **Patch rim** — the rim is
  grown with regular sites until the lattice covers the surface: the notches of a staircase rim first, then rows, two
  rows per click. A new site lies on the surface; past the boundary it lies on the surface continued (a quadric fitted
  to the carrier within 2 a of the boundary: a plane stays a plane, a cylinder or cone keeps its curvature), at most
  1.2 a out. A fan is added only where the surface goes on more than 0.15 a outside a rim edge. **Undo** takes back
  one edit (40 are kept).
* **The line under them** — `kink 1.10° rms · 4° max`: how much the laths turn at their crossings. Hover it for the
  values as read out and between the vertices, the lattice lines, the edge lengths and the largest move.
* **Export** — a frozen lattice is written with its own positions and its edited connectivity, sites past the boundary
  included. The first line of the OBJ (a, scale) keeps its form, as WeaveTogether reads it; a last comment line says
  `# frozen at step …, n connectivity edits`, and the file name ends in `_frozen`.

Also on 24 Nov:

* **Picked centre on a refined carrier** — on ×2, ×4 or ×8 a centre picked on the surface is usually a vertex the ×1
  carrier does not have; the seed chart was asked for it there, came back NaN, and the field turned black and stood
  still. The centre now goes to the nearest ×1 vertex (at most about half a ×1 edge away; the status bar says how
  far), and a seed that is not finite is never handed to the solver.
* **Edge lengths in the growth legend** — after `predicted`: `edge avg · max · min` of the shown lattice, in model
  units (the units of `a` in `2 Density`), and below it in centimetres when the scale is not 1. `a` is in model units
  and does not change with the size set in `1 File → import options`; centimetres = model units × scale, for `a`
  (now also shown in `2 Density`), the edges and the bounding box alike. The mean edge is `a` within about 2 %.

* **Top bar** — the file name · the vertices of the lattice · the singular sites · the mean edge (model units, and
  centimetres when the scale is not 1), then what the last action did, and a light / dark switch at the right (as in
  WeaveTogether; it follows the system until it is used). The theme row left `Display`. The carrier's numbers, the
  step and the solver backend (WebGPU or CPU) are on hover: over the file name, the lattice numbers and `TryTriWeave`.

`node test_latedit.mjs` checks the edits on lattices with a known answer: the valence rules; each of flip, split and
collapse undone by `heal`; a lone 5 and a lone 7 restored from the clusters one flip makes of them; random scrambles
of one edit always healed back, of two or three always lowered (about two in three all the way: the rest end as a
5–7 pair three or four edges from the 5); sites taken out of a rim put back where they were; a lattice smaller than
its carrier grown until it covers it, the sites past the boundary on the plane or the cylinder continued.
Naked edges and `to rim` (test 6): a collapse and a split of a naked edge undone by `heal`; two 5–7 pairs four and
seven rows inside a regular lattice walked out in 5 and 7 edits, the lattice regular again; a lone 5 nine rows in
out in 85 edits, a lone 7 eight rows in in 46, no triangle folded. On Ji's read-out of the canopy (494 sites, a = 25.97,
picked centre: 4 × 5, 32 × 7, 1 × 8), one `to rim` click on each 5: all four out in 36 edits and 0.5 s (three as pairs,
one alone), 27 × 7 and the 8 left, no triangle folded, lath kink 4.78° → 4.04° rms at 20 sweeps — the same lattice
edited by hand: 109 edits, 1 × 5 and 33 × 7 left.
On `canopy_tree_round` ×4, a = 16.16, picked centre, 12 000 steps: `Heal all` 116 → 64 singular (12 clusters, 0.6 s,
Σ(6−n) unchanged), `Patch rim` 81 sites past the edge of the disc and the foot of the trunk (0.1 s), no triangle folded.

Smoothing in numbers (20 Nov; the panel showed this table then, now the line under the tools and its hover).
Lath kink, rms in degrees (global seed with uniform a, default r, b, dt, 2500 steps, carrier ×1; `node test_smooth.mjs <name>`):

| surface | sites · non-6 | as read out | sub-vertex | 20 sweeps | of that: regular · by a 5 or 7 · rim |
|---|---|---|---|---|---|
| HalfTube | 147 · 0 | 5.55 | 4.58 | 0.55 | 0.24 · – · 0.85 |
| developableSrf1 | 161 · 0 | 3.77 | 2.35 | 0.67 | 0.54 · – · 0.85 |
| v4Cone_sharp | 104 · 2 | 5.44 | 4.20 | 1.58 | 0.90 · 5.4 · 0.95 |
| saddle1 | 123 · 1 | 5.39 | 4.06 | 2.30 | 1.38 · 8.3 · 1.47 |
| v5Cone_round | 101 · 1 | 5.41 | 4.58 | 2.94 | 1.62 · 12.4 · 1.39 |
| cylinder (local seed, 7000 steps) | 259 · 14 | 5.66 | 4.55 | 4.00 | 1.67 · 8.7 · – |
| genus1_torus | 204 · 24 | 7.75 | 7.03 | 6.59 | 2.35 · 10.5 · – |
| Bunny ×2, a = 8.1 (3000 steps) | 346 · 67 | 8.46 | 7.83 | 6.83 | 3.18 · 9.4 · – |
| tree | 117 · 29 | 9.68 | 9.52 | 7.76 | 3.76 · 8.8 · 6.8 |

* Reading the table: snapping to carrier vertices is the zigzag between regular sites (3–5° at 9 carrier edges per a;
  the sub-vertex positions take it to 1–2° on the simpler surfaces, smoothing to under 1° on developable parts). On an
  open surface the field places the outermost row unevenly (4–8° in the rhombi at the rim); smoothing with the rim free
  takes that to about 1° (not on the tree). What is left sits at the 5s and 7s.
* What smoothing cannot remove: around a lone 5 or 7 every lath turns by up to a / (6 r) per crossing at distance r
  (4.8° two edges away, 1.9° five, 1° ten). `node test_smooth.mjs` (test 4) builds that lattice, checks the law and
  shows that the sweeps leave it. Only healing lowers it: fewer pairs, lone 5s and 7s on the curvature that pays for them.
* 100 sweeps gain little over 20 (HalfTube 0.55° → 0.45°) and move the sites further (Bunny ×2: largest move 0.36 a →
  0.56 a). A sweep takes about half a millisecond per 300 sites.
* Why this smoothing (20 sweeps, lath kink rms; sub-vertex is the start):

  | | HalfTube | cylinder | saddle1 | Bunny ×2 |
  |---|---|---|---|---|
  | sub-vertex, not smoothed | 4.58 | 4.55 | 4.06 | 7.83 |
  | to the centroid of the neighbours (Laplacian) | 14.2, edges 0.45 a | 6.06 | 14.1, edges 0.34 a | 19.2, 15 triangles folded |
  | the same, rim held | 4.27 | 6.06 | 3.88 | 19.2 |
  | parallelograms in space (not unfolded) | 0.41 | 5.55 | 2.28 | 9.75 |
  | parallelograms unfolded, no springs | 0.35 | 3.95 | 2.27 | 6.81 |
  | **parallelograms unfolded, springs 0.3 (built)** | 0.55 | 4.00 | 2.30 | 6.83 |
  | the same, rim held | 4.18 | – | 3.77 | – |
  | the same, 5s and 7s held | – | 4.10 | 2.30 | 7.30 |

  The Laplacian shortens the lattice lines instead of straightening them: with the rim free an open lattice collapses,
  with the rim held little changes, on a closed surface it gets worse. Parallelograms taken in space fight the bending
  of the surface (cylinder, Bunny). Unfolded they do not. The springs cost 0.2° on the HalfTube and nothing elsewhere;
  they keep the sites from wandering (Bunny ×2 after 400 sweeps: 0.74 a with springs, 1.14 a without, the same kink).
* Files: `js/smooth.js` (freeze, relax, lathBend, lineBend, shape, project), `test_smooth.mjs`; `app.js` (freezeNow,
  unfreeze, applySmooth, syncFrozenUI, healNow(which)); `view.js` (site, sitePos, siteNormal, sitesFlip). From the
  console: `LocalSH.frozen` is the frozen lattice ({P, nrm, tris, valence, …}), `LocalSH.api.setSmooth({sweeps, mu,
  holdRim, holdSing})`.

## Size, strips and the display chips (2026-11-08, after WeaveTogether)

* **Scale** — `1 File → import options → scale`: `× model units → cm` (× 1 reads the file's units as centimetres,
  × 0.1 as millimetres), or type one size of the bounding box (`size` W × D × H, cm) and the factor follows. The
  read-out gives the lattice edge a in centimetres. Only the sizes in cm follow: the carrier, the lattice and the
  export stay in model units. The factor is kept for the next file. `Export lattice` writes it on the first line of
  the OBJ (`scale = … cm/unit`); WeaveTogether reads it and opens the lattice at that size.
* **Bounding box** — a chip in Display: the carrier's box, dashed, with W / D / H in centimetres.
* **Laths as strips** — width × thickness swept along a smooth spline through the crossings (as WeaveTogether and
  WeaveSim draw them): under, a strip lies on the surface; over, it rides on the strip under it. `thickness [a]`
  sits under `lath width [a]` in Display; `over / under` off lays every strip on the surface.
* The `mesh`, `peaks` and `rotation field` chips are hidden (the layers are still in the code: remove `hidden` from
  their buttons in `index.html` to bring them back).

## Speed: WebAssembly and solver threads (2026-10-31)

* The CG solves of the global seed (rotation field, unit-length steps, stripes) run in WebAssembly (`js/pcg_wasm.js`;
  the kernels `pcg_kernels.c` are embedded in `js/pcg_wasm_bin.js`) and, where the page is cross-origin isolated, on
  5 threads (`js/parallel.js`, `js/pcg_worker.js`) — the same results bit for bit as the JS solves. The status line
  says what ran: `· 5 threads, wasm`, `· 1 thread, wasm` or `· 1 thread, js`, and so does the message at the end
  (`seeded … 1.3 s on 5 threads`).
* The threads check themselves: on the first solve of each size they time a few iterations against one thread, and
  they keep timing while they run. Where they are slower (the browser at low priority, e.g. Chrome behind another
  app) the rest runs in one thread, with the same results, and the message ends `on 1 thread (the threads measured
  slower)`. Bunny ×2 in a Chrome behind the Claude app: 12–13 s on the threads before this check, 5–7 s with it.
* Isolation: on GitHub Pages the service worker `coi-sw.js` provides it (the first visit reloads once by itself).
  `python serve.py` here sends no isolation headers (one thread); `python ../multi_threading/serve.py` (port 8081)
  does (threads). Locally, `?coi=1` in the address turns the service worker on, `?coi=0` removes it.
* Bunny ×2 (cones from the curvature): 6.3 s before 31 Oct → 1.8 s in one thread with WebAssembly → 1.3 s with the
  threads. Keep the tab in front while it solves: a background tab runs 3–6× slower.
* Rebuilding the kernels needs clang with the wasm32 target (Homebrew `llvm`): `python3 ../multi_threading/wasm/build.py`.

## Size field and the direction layers (2026-10-13)

* Section 3 → **conformal size field**: the lattice edge becomes a·e^{−φ} with φ from the surface curvature and the
  cones (the centre, plus any you add with `+ cone…`). On the domes this turns {5:10, 6:112, 7:8} into {4:1, 6:98}.
  Watch the readout: if it warns about resolution in ĝ, click `×2` in section 2.
* (The two-stage anneal of section 4 was removed on 2026-10-28: on the zoo it gave no lasting gain — use `Heal pairs`.)
* Section 5 → chips **strain** (edge length per direction family) and **orient** (lattice orientation: grains and
  seams at a glance).

## Add a surface

* In the page (online too): `Open mesh…` in `1 File`, or drop the file anywhere on the page. STL, OBJ, PLY and OFF are
  read, repaired and remeshed in the browser, in a background worker — nothing is uploaded. About 1 s for a zoo model,
  5 s for a 350 k-triangle scan. The repairs (cracks closed, caps stripped, holes filled, non-manifold parts split…) are
  listed in the card of the `1 File` mark. `up` (beside `Open mesh…`) says which axis of the file points up: `z` for Rhino
  and most CAD, `y` for most web models (glTF, Blender, Sketchfab); changing it turns the opened model at once (no
  remesh) and is kept for the next file. `import options`: `edges / bbox` (120 → h = bbox/120, as export_carrier.py),
  `max v` (26 000; h is coarsened to fit), `crease °` (40: sharper edges are kept as feature lines), `strip planar caps`;
  `Remesh again` reruns the last file with them. `Download .lshc` saves the carrier; dropping that file later loads it
  without the remesh. A non-orientable surface (a Möbius band) is cut open for the remesh and glued back: the card says
  `re-glued the cut after the remesh`, χ is 0, and the lattice runs across the seam (drawn on both sides).
* Into the zoo list: put the file in `../model/zoo/` and run `python ../export_carrier.py name` in the `weavesim` env
  (needs `pip install gpytoolbox potpourri3d`). This writes `carriers/name.lshc` and adds it to `manifest.json`.
  Reload the page. (A downloaded `.lshc` copied into `carriers/` and added to `manifest.json` works the same.)

## Update the online version

Copy `index.html`, `style.css`, `manifest.json`, `Tri-Weave-Lattice-Equations.html`, `coi-sw.js`, `carriers/`, `js/`,
`vendor/`, `serve.py` and the two `.md` files into `TryTriWeave/web/`, then commit and push; GitHub Pages serves the folder as is
(<https://uwenji.github.io/try-tri-weave/web/>). The test, calibration and artifact files stay in `localize/web`.

## Layout (2026-10-18)

* **Status marks** — the coloured ⓘ next to the growth legend, the zoo list, `lattice edge a` and
  `Heal 5-7-5` carry the read-outs: hover (or tap) for the card. (Since 2026-11-20 the heal read-out is in the panel, section 5.) Green fine · orange look at it · red a problem · grey
  nothing yet. The density diagram is in the card of `lattice edge a`.
* **Growth chart** — the chevron left of the legend folds it; the legend, the status mark and the solver controls stay.
* **Display** — the sidebar beside the view (theme, layers, singular chips, laths, input normal, carrier colours, turntable).
  Its header folds it to a narrow tab; ⇄ moves it to the other side of the view. Both are remembered.
* A click on a lath or a vertex reports in the status bar.
