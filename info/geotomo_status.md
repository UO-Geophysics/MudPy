# Geodetic Source Tomography (`geotomo`)

Current implementation status, November 9, 2026.

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
| `src/python/mudpy/geotomo_io.py` | HDF5 observations, full NEU response libraries, block projection, and predictions |
| `src/python/mudpy/parallel.py` | `run_parallel_green` and `run_parallel_synthetics_geotomo` workers |
| [Starter TOML](../run/geotomo/template.toml) | Project-relative inputs, forward settings, source model, damping, and run stages |
| [Example README](../examples/geotomo/README.md) | Usage, conventions, validation, limitations, and ownership |
| [Validation driver](../examples/geotomo/validate.py) | Homogeneous references and layered/off-grid experiments |

## Supported behavior

- Full six-component tensors or prescribed tensors with signed scalar amplitudes.
  Sources are integrated moments at specified points, not moment densities.
- Source-parallel FK static responses including both double-couple and explosion terms.
  Tensor coordinates are NED; displacement/projection vectors are NEU.
- Explicit TOML `execution.processes` caps concurrent FK/syn subprocesses
  (default 1), with isolated source folders, bounded batches, and one response
  writer preserving source order. Solver execution and BLAS threading are unchanged.
  Serial/two-worker FK reference products agree exactly; all 46 regression
  methods pass, including overlap, ordering, single-writer and failure checks.
- Scalar imagery projections in metres, source-file ordering, and checked
  response completeness. Legacy slip-inversion interfaces retain their behavior.
- Weighted data-space least squares and NNLS, with identity damping through
  augmented rows. Spatial smoothing and other source penalties remain future work.
- One project per site/event, descriptive run names, exact configuration
  snapshots, resolved input/code/executable identities, and verified cache reuse.
- HDF5 input from preprocessing packages, preserving mean look coefficients and
  pixel membership while modeling one representative location per sample.
  Source-block NEU response storage and direct unweighted NumPy matrix reuse.
- Exclusive `.geotomo.lock` ownership for top-level runs. Existing projects are
  never initialized again, and existing named results are protected.

## Validation and limits

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
