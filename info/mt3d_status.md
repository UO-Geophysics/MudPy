# MudPy `mt3d` — volumetric (distributed point-source) MT inversion: current state

Summary of what exists in the repo today for the InSAR inversion with multiple
volumetrically distributed sources. Code dates to April 2017 (`4ceb568 started
mt3d module`), last meaningfully touched in the July 2019 Python-3 cleanup.
Nothing else in the package imports it and there are no driver scripts in
`run/`, `examples/`, or `scripts/` — it was always driven by ad-hoc scripts
calling `mudpy.mt3d` directly.

## Concept

- The "fault" file is a cloud of **point sources**, 4 columns: `num, lon, lat, z(km)`.
  No geometry, no slip vectors.
- Each source carries a **full 6-component moment tensor** (Mxx, Mxy, Mxz, Myy,
  Myz, Mzz) as free parameters → model vector is `6 × Nsources`.
- Alternate mode `forceMT=True`: the MT shape is fixed by the user and only one
  scalar amplitude per source is inverted.
- The pipeline mirrors the classic MudPy static slip inversion
  (GFs → synthetics → G assembly → regularized least squares) and reuses
  `inverse.py` for the data vector, weights, VR, and synthetics output.

## Files

| File | Role |
|---|---|
| `src/python/mudpy/mt3d.py` (909 lines) | Orchestration: GF/synthetics dispatch, G assembly, inversion, output writers |
| `src/python/mudpy/parallel.py` — `run_parallel_synthetics_mt3d` (line 838) | Per source–station synthetics via FK's `syn` executable |
| `src/python/mudpy/parallel.py` — `run_parallel_green` (line 6) | fk.pl static GFs; `insar=True` just changes the file suffix |
| `src/python/mudpy/inverse.py` | Reused as-is: `getdata`, `get_data_weights`, `get_VR`, `write_synthetics`, `mini_station_file` |

## Pipeline detail

1. **GFs** — `mt3d.inversionGFs` reads the `.gflist`; column 3 = static/GPS,
   column 7 = InSAR LOS. Calls `make_parallel_green` → MPI →
   `parallel.run_parallel_green` in static mode: `fk.pl -M<model>/<depth>/f -N1`
   per source depth. Output: `GFs/static/<model>.static.<depth>.sub<NNNN>.insar`
   (or `.gps`).

2. **Synthetics** — `make_parallel_synthetics` → MPI →
   `parallel.run_parallel_synthetics_mt3d`. For each source–station pair the
   static GFs are piped to FK's `syn` **six times**, once per unit MT component
   (unit moment Mw 5.0, `M0 = 10^(1.5·5+9.1)·1e7` dyne-cm). Radial/transverse
   rotated to N/E, cm → m, and for InSAR the (n, e, u) vector is dotted with a
   per-station LOS unit vector read from `data/statics/<sta>.los`. Per-process
   temp files are concatenated into `GFs/static/_Mxx.los … _Mzz.los`
   (or `_M??.neu` for GPS).

3. **G assembly** — `getG`/`makeG` in mt3d.py. InSAR: one row per station,
   6 columns per source, MT components **interleaved** (source k occupies
   columns 6k…6k+5). Optionally cached to `GFs/matrices/<G_name>.npy`.

4. **Inversion** — `run_inversion`: normal equations
   `K = GᵀG + λ²LᵀL` with **L = identity** (0th-order Tikhonov damping only —
   no spatial smoothing between sources), swept over a list of λ values.
   Solver: `numpy lstsq` or `scipy nnls`. Outputs per λ:
   - `.inv` file: per-source MT components (Nm), both nodal planes
     (via `analysis.MT`), per-source M0
   - `.log` file: L2, Lm, VR (static + InSAR), total M0, Mw
   - synthetics via `inverse.write_synthetics`

   Total moment: per-source `M0 = sqrt(0.5·Σλᵢ²)` from MT eigenvalues, summed,
   then `Mw = (2/3)(log10 M0 − 9.1)`.

## What actually works vs. doesn't

- ✅ **Full-MT (`forceMT=False`) + InSAR** is the only completed path — GF
  computation, synthetics, G assembly, inversion, and output all wired up.
- ❌ **Static/GPS G assembly is unfinished**: in `makeG` the block that fills G
  for the `.neu` case is commented out with a `##Needs work` marker
  (mt3d.py ~line 507). It returns a matrix of zeros, so a joint GPS+InSAR
  inversion would silently use zero GFs for the GPS rows.
- ❌ **`forceMT=True` + InSAR crashes**: references an undefined `kfault`
  (mt3d.py ~line 607) → `NameError`.
- ⚠️ `nnls` is applied to the normal-equations system `K·m = Gᵀd`, so
  positivity is enforced on the solution of the squared system — not the same
  problem as `min ‖Gm−d‖ s.t. m ≥ 0`.
- ⚠️ Regularization is identity-only; no Laplacian/smoothing operator across
  the source cloud was ever built for this module.
- ⚠️ `write_log` hardcodes `VR[0]` (static) and `VR[4]` (InSAR) — assumes the
  VR vector layout from `inverse.get_VR`.
- ⚠️ Model-vector indexing in `get_moment`/`write_model` uses
  `6*arange(len(sol)/6)` — float indices under Python 3; would need `//`.

## Inputs required (for the working path)

- `data/model_info/<source_name>`: `num lon lat z(km)` per point source
- `structure/<model_name>`: fk-format 1D velocity model
- `data/station_info/<GF_list>`: standard MudPy gflist, column 7 = 1 for InSAR
- `data/statics/<sta>.los`: per-station LOS data; the LOS unit vector is read
  from the same file (first value skipped, then n/e/u components)
- `$MUD` env var set; `fk.pl` and `syn` (FK package) on PATH; MPI via `mpiexec`
