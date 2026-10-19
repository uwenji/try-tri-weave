# Try-Tri-Weave Tool
Open link: https://uwenji.github.io/try-tri-weave/web/

## Frist time in the app
the solver controls sit at the top of the left bar, right under the growth strip, so they are
   always in reach: drag `steps` (up to 20 000) left to scrub back in time, right to compute further; `Run` keeps
   going, `+500` adds a batch, `Re-seed`; `max steps` and `per frame` below them.
1. **File** — pick `v4Cone_sharp` from the zoo list or immport your model. `Export lattice` is in `File`. The carrier loads, the singularity is detected
   (Ω = ∫K → suggested valence and centre, blue marker), the seed is planted and 3000 steps start running.
2. **Density** — drag `lattice edge a` (model units; strand pitch = 0.866 a, lath width = slider × a). The pattern
   re-seeds at the new spacing. The slider stops at 6 carrier edges per a: for a denser weave click `carrier ×2 / ×4 /
   ×8` (subdivides the carrier in the browser; ×4 and ×8 want the GPU) or export a finer carrier
3. **Singularity** — change the `valence at the centre` slider (3–9) to plant a different singularity;
   Shift-click the surface (or `Pick centre on surface…`, then click) to move the centre; `Use detected` goes back.
4. **Solve** —  Section 4
   holds the pattern parameters: `r`, `b`, `dt` change the equation without resetting the field.
5. **Display** (the sidebar beside the view) — `theme` auto / light / dark. Chips toggle layers: `field · mesh · lattice ·
   peaks · singular · centre · laths`; the `singular` row has one chip per valence present (v3, v5, v7 …, with counts) so
   each kind of singularity can be shown or hidden on its own. Singular vertices are drawn as discs with the valence
   digit inside. Turn **laths** on and move `lath width [a]`; set `carrier shows` to *plain* to see the weave alone.
   `input normal` flips the side of the surface that carries the laths and marks (a branching tube has no side that
   is outside everywhere). The **growth chart** at the top of the left bar plots peaks · interior · singular against
   the step; the chevron folds it, and clicking or dragging on it scrubs the shown step; the solver controls sit right
   under it. The four sections below fold; `5 Explain` opens the equations page.
6. **Selecting** — click a lath to follow its whole path (the rest dims, a dark line runs along it, the status bar
   gives family, knot count, over/under and open/closed); click a singular disc or a peak to select that vertex
   (a larger disc shows its lattice index, which is its line in the OBJ export). Esc clears. Moving the centre is
   now Shift-click, or `Pick centre on surface…` then click, so an ordinary click never re-seeds by accident.
7. **Navigating** — the **view cube** at the top right of the view turns with the model: click a face, an edge or
   a corner to look from there (FRONT is −y, RIGHT is +x, TOP is +z), drag the cube to orbit, ▲▼◀▶ roll the
   model by 90°, ⌂ is the home view, ⤢ fits the surface, ⟳ starts the **turntable** (speed in Display,
   pauses while you drag), the last button switches perspective / orthographic. Keys: `H` home, `F` fit,
   `T` turntable, `P` projection.

more detail please read in [Tri-Weave-Lattice-Equations.html](https://uwenji.github.io/try-tri-weave/web/Tri-Weave-Lattice-Equations.html) (also `6 Explain` in the app's left bar; source `web/Tri-Weave-Lattice-Equations.html`)

More about our group: https://royaldanishacademy.com/en/case/kagome-architectures

Feedback is welcome - please open an issue.
