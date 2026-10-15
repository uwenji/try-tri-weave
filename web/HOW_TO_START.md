# How to start — LocalSH Weave Lab (web)

Three ways, from quickest to most flexible.

## 1 · Online, nothing to install

Open the Claude artifact **LocalSH Weave Lab** (link in the chat, or in your artifact gallery) in Chrome, Edge
or Safari. The pill at the top right should read `WebGPU · …`; if it reads `CPU (JS)` the browser has no
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
device either open the Claude artifact (it is https) or run

```bash
python serve.py --https             # https://192.168.x.x:8443/, self-signed certificate made once with openssl
```

and accept the certificate warning on the device the first time.

## 3 · First five minutes in the app

1. **Surface** — pick `v4Cone_sharp` from the zoo list. The carrier loads, the singularity is detected
   (Ω = ∫K → suggested valence and centre, blue marker), the seed is planted and 3000 steps start running.
2. **Density** — drag `lattice edge a` (model units; strand pitch = 0.866 a, lath width = slider × a). The pattern
   re-seeds at the new spacing. The slider stops at 6 carrier edges per a: for a denser weave click `carrier ×2 / ×4 /
   ×8` (subdivides the carrier in the browser; ×4 and ×8 want the GPU) or export a finer carrier
   (`python ../export_carrier.py name --edges 240 --max-verts 120000`).
3. **Singularity** — change the `valence at the centre` slider (3–9) to plant a different singularity;
   Shift-click the surface (or `Pick centre on surface…`, then click) to move the centre; `Use detected` goes back.
4. **Solve** — drag `steps` (up to 20 000) left to scrub back in time, right to compute further; `Run` keeps
   going, `+500` adds a batch; `r`, `b`, `dt` change the equation without resetting the field.
5. **View & readouts** — `theme` auto / light / dark. Chips toggle layers: `field · wire · lattice · peaks ·
   singular · centre · laths`; the `singular` row has one chip per valence present (v3, v5, v7 …) so each kind
   of singularity can be shown or hidden on its own. Singular vertices are drawn as discs with the valence
   digit inside. Turn **laths** on and move `lath width [a]`; set `carrier shows` to *plain* to see the weave
   alone. The valence histogram, Σ(6−n) and the strand/knot counts update with the lattice.
   `weave side` flips the side of the surface that carries the laths and marks (a branching tube has no side that
   is outside everywhere). The **growth strip** at the top of the left bar plots peaks · interior · singular against
   the step; it never folds, and clicking or dragging on it scrubs the shown step. The five sections below fold.
6. **Selecting** — click a lath to follow its whole path (the rest dims, a dark line runs along it, the readout
   gives family, knot count, over/under and open/closed); click a singular disc or a peak to select that vertex
   (a larger disc shows its lattice index, which is its line in the OBJ export). Esc clears. Moving the centre is
   now Shift-click, or `Pick centre on surface…` then click, so an ordinary click never re-seeds by accident.
7. **Navigating** — the **view cube** at the top right of the view turns with the model: click a face, an edge or
   a corner to look from there (FRONT is −y, RIGHT is +x, TOP is +z), drag the cube to orbit, ▲▼◀▶ roll the
   model by 90°, ⌂ is the home view, ⤢ fits the surface, ⟳ starts the **turntable** (speed in section 5,
   pauses while you drag), the last button switches perspective / orthographic. Keys: `H` home, `F` fit,
   `T` turntable, `P` projection.

What "ideal" looks like: exactly one singular peak and Σ(6−n) = round(Ω/60°). Try the three
`ideal cone` entries first, then your own cones.

## Size field and the two-stage schedule (2026-10-13)

* Section 3 → **conformal size field**: the lattice edge becomes a·e^{−φ} with φ from the surface curvature and the
  cones (the centre, plus any you add with `+ cone…`). On the domes this turns {5:10, 6:112, 7:8} into {4:1, 6:98}.
  Watch the readout: if it warns about resolution in ĝ, click `×2` in section 2.
* Section 4 → **anneal after growth**: once the surface is covered, the run continues at (r₂, b₂) where dislocations
  are mobile. Use it for isolated extra 5–7 pairs; a seam between fronts that met is a holonomy mismatch and stays —
  move the centre instead (seed on the flat part, so the merge happens at the junctions).
* Section 5 → chips **strain** (edge length per direction family) and **orient** (lattice orientation: grains and
  seams at a glance).

## Add a surface

* Quick: drop an `.stl` on the view (welded, Loop-subdivided and Delaunay-flipped in the browser; not remeshed,
  the panel says so).
* Proper: put the file in `../model/zoo/` and run `python ../export_carrier.py name` in the `weavesim` env
  (needs `pip install gpytoolbox potpourri3d`). This remeshes it like the Python app, writes
  `carriers/name.lshc` and adds it to `manifest.json`. Reload the page.

## Update the online version

```bash
python build_artifact.py        # index.html → artifact.html (body only, as the artifact host wants it)
```

then publish `artifact.html` with the same supporting files (the carriers go up under the `.wasm` extension,
`manifest_artifact.json` is served as `manifest.json`) — ask Claude to republish, it has the file map.
