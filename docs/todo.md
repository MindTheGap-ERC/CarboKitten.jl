TODO
====

Project Summary
---------------

CarboKitten is models the production, accumulation and transport of sediment on carbonate platforms. It usually does this with time steps ranging from a decade to a few hundred years, over the span of a few million years. The model is written in Julia. The entire source code is contained in its documentation using **literate programming** with [Entangled](https://entangled.github.io), located at `docs/src`. Each model in CarboKitten is composed from a set of **components** identified in `docs/src/componenents` using a data composition approach from [ModuleMixins.jl](https://jhidding.github.io/ModuleMixins.jl). Components are combined into **models**, which are found in `docs/src/models`.

The principal model that first time users will run is called ALCAP, which combines the following components:

- The Bosscher & Schlager (1994) model for **sediment production**.
- A **cellular automaton** surrogate model to add spatial heterogeneity to the sediment production.
- An **active layer** modelling entrained sediment from production and disintegrated older sediments.
- An **advection** based finite-difference transport model for sediment dispersal within the active layer.

Additionally, there are interfaces to supply variable **insolation** and eustatic **sealevel** to the production model, and **subsidence rates** can be variable in both time and space.

CarboKitten supports multiple box topologies, with the goal of either modelling a strip of coast line, or a full atol.

Output can be written to HDF5 files or memory. There's a range of options for post-processing and visualization. The visualization is implemented as an optional extension in `docs/src/visualization` that is only loaded when `Makie` is available. Generic algorithms are in `docs/src/algorithms`.

The rest of the documentation is sometimes still a bit scattered. We're aiming to have a user targeted manual in `docs/src/manual`, although that is currently still scattered (see `docs/make.jl` for the index). A set of how-to, or case based examples in `docs/src/cases` is well under way.

Past Milestones
---------------

- 0.8.0 release summary:
    - major performance improvements
    - upgrade to improved component structure
    - variable subsidence feature
    - circular sediment buffer implementation
    - new data api with packs, bundles and data
    - improved denudation models, including `WithoutCADenutation`
    - many smaller bug fixes

Done since Latest Milestone
---------------------------

Tasks Remaining
---------------

The following is a list of tasks, not always with a issue attached. Before an issue is opened a task needs to be fleshed out. Add a reference to the issue here.

### Model Features

- [ ] Compactification (dynamic), and flexing (post-processing?). This could be achieved by adding a compactification component. This component adds a new field to the state, where for each cell in a 3-dimensional buffer we need to track its physical height (in contrast to the current sediment buffer, which treats a fixed height and mass as identical units of meters). In this case we need to check that the buffer doesn't overflow, or at the least that any sediment leaving the buffer is fully compactified. The `bathymetry` needs to be updated to model the additional subsidence caused by compactification.
- [ ] Highstand shedding. The idea is that production in shallow regions gets transported much further because sediment gets suspended in the entire water column, and is transported out by undercurrents. This transport mechanism is different from our current model, which assumes that the transport happens close to the sea floor (call this benthic transport). We'd need to add a second active layer that quantifies pelagic transport. We need to consult the literature on what is needed here, but I assume the finer granules get suspended here and are transported out by a large diffusion coefficient, and this transport being less dependent on slope steepness, so explicit diffusion is what's needed here.
- [ ] Denudation as an `AbstractDenudation` role, with one component per denudation model (`ChemicalDissolution`, `PhysicalErosion`, `EmpiricalDenudation`, `NoDenudation` as default). Once `NoDenudation` is the default, `WithDenudation` can be merged into ALCAP and `WithoutCADenudation` into `WithoutCA`. Coordinate with ongoing denudation work. Issue #293 (supersedes #225).
- [ ] Temperature and saturation-index dependent kinetics for chemical dissolution. Fits in the `ChemicalDissolution` component of #293. Issue #66.
- [ ] Carbonate factory delay: a factory density state variable that dies off under subaerial exposure and regrows after resubmersion. Issue #112.
- [ ] Use critical angle theory to determine the amount of entrained sediment in the active layer. Issue #101.
- [ ] Wave energy shielding by reef barriers; needs literature input on heuristics. Issue #268.
- [ ] Multi-scale workflow: a higher-level CA state that guides local producer establishment, for basin-scale trends. Issue #257.

### Post-processing and Visualization Features

- [ ] Subsidence plotted along sealevel: this is possible when subsidence is constant or only time-dependent, but more complex when there is a spatial dependence. In that case we could plot quartiles and extrema with some shading.
- [ ] Extract unconformities from data. This is already done for visualization, but needs to be accessible to users as well. Related: issue #235, where `plot_unconformities` still uses `water_depth < 0`, which is wrong in the presence of denudation.
- [ ] Identification of onlap and offlap regions.
- [ ] Issue #257 may have some overlap with the onlap/offlap item.
- [ ] Horizontal facies distribution plot at a given time step. PR #244 was closed so the time-based slicing (preserved and deposited) could be redone through the sediment buffer interface, which is now possible since #288. Issue #159.
- [ ] Post-deposition facies classification into user-defined facies (sediment proportions, palaeo-water depth, wave energy). Issue #260, PR #261 open.
- [ ] Rich `text/html` summary for Pluto/Jupyter; should target `Pack` rather than `Header`. Issue #130.

### Quality

- [ ] We need to be very explicit about naming and add a section in the manual about this.
- [ ] Complete the user manual, rearrange the files to reflect this arrangement.
- [ ] Remove dead code, issue #144. `Vectors` and the dead `Boxes` code are handled in PR #292 (open). Remaining: `Utility.select` and `src/Transport/ActiveLayer.jl`.
- [ ] Improve test coverage, issue #143. Remaining: `DataSets`, `Skeleton`, `CAProduction`, `ActiveLayer` header metadata, and a full `run_model` test for `WithDenudation`.
- [ ] Update the examples that still use the old facies API (10 scripts, plus `CustomProductionModel.jl` still on `@compose`), and add the denudation run scripts to the Entangled workflow. Bring the `brei` workflow (`uv run brei figures`) up to date and run it in CI before making a release; `.entangled/tasks.json` was missing after tangling with Entangled 2.5.0. Issue #216 (includes #142).
- [ ] Examples for initial topography: multi-linear ramp, conic (atoll), Gaussian random fields. Issue #177.
- [ ] Emphasise ALCAP as the main model in the docs. Wait for #293, which merges the denudation models into ALCAP. Issue #118.
- [ ] Grid edge effect at the landward edge, and artefacts in narrow boxes (10 cells wide). Reports predate v0.8.0; re-check on the current version first. Issue #91 (includes #227).

### Performance

- [ ] There could be a speedup in the transport model by transposing data, setting `facies` as the last dimension in arrays. See issue #289.
- [ ] Add PrecompileTools to improve latency, especially for the plotting extension. Issue #10.

### Stability

- [ ] Stability check in diagnostic mode. The advective CFL condition is handled by the adaptive transport substeps, but the model as a whole is diffusive (disintegration, transport and lithification on a forward Euler step), and that is not checked. A rough estimate is κ ≈ d·C, with stability limit κΔt/Δx² ≤ 1/4, which grows with Δt². Issue #221.


