# Geodetic Source Tomography (`geotomo`)

Current implementation status, November 10, 2026.

Geodetic imaging of distributed subsurface deformation sources, implemented
within MudPy. Imagery-derived displacement is primary; GNSS is complementary
and its response assembly remains unimplemented in this workflow.

## Entry points and files

```python
from mudpy import geotomo

geotomo.init_project('/absolute/path/new_project')
geotomo.run('/absolute/path/new_project/configs/experiment.toml')
```

| File | Role |
| --- | --- |
| `src/python/mudpy/geotomo.py` | GF/synthetic orchestration, matrix assembly, inversion, and model/log output |
| `src/python/mudpy/geotomo_config.py` | Schema-1 TOML, project initialization, input validation, provenance, and cache verification |
| `src/python/mudpy/geotomo_dc.py` | Pure-DC nonlinear least squares, analytic matrix-free Jacobian, multistart diagnostics and TOML adapter |
| `src/python/mudpy/geotomo_io.py` | HDF5 observations/responses, block projection, predictions, and checked result readers |
| `src/python/mudpy/view.py` | Reusable geotomo tradeoff, imagery fit and point-source viewers alongside existing slip viewers |
| `src/python/mudpy/analysis.py` | Checked ISO/DEV/DC/CLVD decomposition of saved geotomo tensors; existing analysis interfaces retained |
| `src/python/mudpy/parallel.py` | `run_parallel_green` and `run_parallel_synthetics_geotomo` workers |
| [Starter TOML](../run/geotomo/template.toml) | Project-relative inputs, forward settings, source model, damping, and run stages |
| [Example README](../examples/geotomo/README.md) | Usage, conventions, validation, limitations, and ownership |
| [Validation driver](../examples/geotomo/validate.py) | Homogeneous references and layered/off-grid experiments |

## Supported behavior

- Full six-component tensors or prescribed tensors with signed scalar amplitudes.
  Sources are integrated moments at specified points, not moment densities.
- Source-parallel FK static responses including both double-couple and explosion terms.
  Tensor coordinates are NED; displacement/projection vectors are NEU.
- Vectorized float64 moment-tensor synthesis evaluates all receivers and basis
  components per source using the `radiats.c` formulas, including explosion terms.
  Production runs require no `syn` executable. The slip path and C code are unchanged.
- Explicit TOML `execution.processes` caps concurrent source workers
  (default 1), with isolated source folders, bounded batches, and one response
  writer preserving source order. Solver execution and BLAS threading are unchanged.
  Serial/parallel response arrays and matrices agree exactly; all 79 regression
  methods pass, including overlap, ordering, single-writer and failure checks.
- TOML preflight rejects FK interface sources by default. Explicit
  `forward.interface_policy = "above"` or `"below"` and
  `interface_offset_km` (default 0.01 km) resolve them once, with adjacent-layer
  and FK precision checks. Requested geometry is preserved; effective geometry
  is archived, hashed, and used by forward calculations, caches and saved models.
  The small actual layered FK interface reference also passed at that milestone.
- Scalar imagery projections in metres, source-file ordering, and checked
  response completeness. Legacy slip-inversion interfaces retain their behavior.
- Weighted data-space least squares and NNLS, with identity damping through
  augmented rows. Spatial smoothing and other source penalties remain future work.
- Pure-DC nonlinear least squares for canonical HDF5 projects, initialized from
  one GCMT NDK event or explicit NED tensor. Positive moments/free rotations,
  scalar-moment L2 damping, deterministic sequential multistart and analytic
  TRF/LSMR Jacobian operators reuse the existing six-component matrix. The
  forward module and its cache hashes are unchanged; no cache migration occurs.
  Checked per-member solver diagnostics expose local convergence/failure and
  between-start tensor/prediction differences, plus per-start/member wall and CPU
  timings. Direct damped NNLS initializes amplitudes robustly for correlated
  volume responses. The Norcia DC sweep and comparison figures are complete.
  Five damping values and 15 starts took 6,205 s wall with four BLAS threads.
  The saved-model reader handles subnormal rounding on zero-moment sources.
- One project per site/event, descriptive run names, exact configuration
  snapshots, resolved input/code/executable identities, and verified cache reuse.
- HDF5 input from preprocessing packages, preserving mean look coefficients and
  pixel membership while modeling one representative location per sample.
  Source-block NEU response storage and direct unweighted NumPy matrix reuse.
- Optional HDF5 `execution.keep_fk_intermediates = false` removes known native
  files and source working folders only after whole-run success and closed-library
  validation. Default retention preserves existing behavior. Cleanup receipts and
  pruned/pending cache state preserve traceability; HDF5/matrix reuse remains valid.
- Exclusive `.geotomo.lock` ownership for top-level runs. Existing projects are
  never initialized again, and existing named results are protected.
- Read-only canonical model/log/prediction readers and quick Matplotlib viewers:
  `view.geotomo_tradeoff`, `view.geotomo_fit`, `view.geotomo_sources`. Figures/axes
  are returned without automatic saving/showing. Depth slices and summed horizontal
  projections distinguish point-source moments from summed-tensor moment; an
  optional external Mw reference does not select a preferred damping. Caller-owned
  axes, colour scales and normalization support project-specific figure layouts.
  Reader-only edits preserve numerical cache identity; complete module hashes
  remain in provenance. Existing Norcia response/matrix caches still validate.
- Source-component depth maps, full/DC beachballs, DC character or explicit-plane
  rake circles, individual-source mechanism sections and both DC planes versus
  depth. ISO/CLVD strengths use signed tensor norms; ISO balance is a defined
  0–1 index, not an admissibility criterion. Decomposition reconstructs each
  tensor and agrees with an independent Pyrocko reference for all 315 member-9
  tensors. No new inversion or legacy mechanism-interface changes were needed.

## Validation and limits

The synthetic forward profiling mode (`validate.py --profile-forward`) separates
FK, vectorized make-synthetics and HDF5 costs for three sources/1,507 receivers.
Make-synthetics took 0.64 s with either one or three workers, replacing the
229.24/123.49 s per-receiver executable path. Sampled C reference comparisons
pass a 6e-4 relative-L2 bound for legacy single precision/printed rounding;
independent homogeneous references and the 40-solution layered/off-grid study pass.

Independent homogeneous isotropic and Okada shear checks pass. Layered checks
cover the homogeneous limit, numerical convergence, superposition, and source
signs. Forty off-grid/noisy solutions expose source spreading and excess moment
despite small displacement residuals; individual-source recovery remains limited.

Run regression checks with:

```sh
python -B -m unittest discover -s src/python/tests -p 'test_geotomo*.py'
```

Independent layered benchmarks, finite sources, realistic correlated noise and
nuisance terms, production imagery ingestion, MPI, and numerical topography
remain open. Equivalent elastic moments alone do not determine volcanic pressure,
magma volume, groundwater extraction, or irreversible compaction.
