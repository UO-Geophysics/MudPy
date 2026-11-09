# Geodetic Source Tomography (`geotomo`)

Geodetic imaging of distributed subsurface deformation sources. The current
implementation uses moment tensors and static imagery projections; the examples
below validate that source basis and its existing FK forward model.

## HDF5 observations and elastic responses

The starter TOML now selects HDF5 observations. `h5py` is required; the tested
environment has Python 3.8, NumPy 1.21.5, h5py 2.10.0, and HDF5 1.10.6. No new
dependency installation is performed by the workflow.

Import an existing geographic preprocessing package into the project:

```python
from mudpy import geotomo

geotomo.import_observations(
    '/absolute/path/preprocessing_run',
    '/absolute/path/Norcia2016/data/statics/T44.h5',
    id_prefix='T44_')
```

This reads `samples.txt`, `processed.npz`, and `processing.json`, verifies their
agreement, and preserves input hashes, processing metadata, and native-pixel
membership. Coordinates are the stored representative pixels, not recomputed
centroids. LOS remains in metres, positive toward the satellite. Mean ENU
coefficients are reordered to NEU without normalization. The importer requires
documented signs/units and geographic coordinates; it does not reproject or
invent an uncertainty estimate. Native rasters are not duplicated in this file.

The observation file stores `ids`, `coordinates` (longitude/latitude degrees),
`los_m`, `projection_neu`, optional calibrated `sigma_m`, and a `sampling` group.
Metadata explicitly retains unresolved source datum/reference information.
Imported tile scatter remains a diagnostic, separate from `sigma_m`.

Select the file and a descriptive response-library basename in TOML:

```toml
[inputs]
sources = "cloud.sources"
earth_model = "layered.mod"
observations = "T44.h5"          # In data/statics/
observation_format = "hdf5"

[responses]
name = "t44_layered_v1"         # GFs/static/t44_layered_v1.h5

[inversion]
source_model = "full_tensor"
solver = "lstsq"
weighting = "none"              # Until calibrated sigma_m is supplied
damping = [0.0, 0.1, 1.0]
```

HDF5 observation and response schemas are version 1, independent of TOML schema
version 1. The runner writes no per-observation `.los` files for HDF5 inputs.
It evaluates the elastic response once at each representative observation
location. A cell-mean observation is approximated by the prediction at that
point; retained membership is provenance, not an instruction to compute kernels
at all native pixels. Background correction remains the supplied preprocessing
choice, not an automatically fitted inversion nuisance model.

FK's complete native DC+EX files remain backend products. Make-synthetics
evaluates the `mt_radiat` formulas from `radiats.c` with NumPy across all receivers
and tensor components. Responses are rotated to NEU, converted to metres, and written directly into a float64 HDF5
array with axes `(source, receiver, NEU, tensor)`. Tensor components are
`xx,xy,xz,yy,yz,zz` in NED coordinates, with one coefficient scaling a component
by `10**16.6 N m`. Source IDs, coordinates, receiver IDs, conventions, and run
provenance accompany the array. Chunking uses one source, up to 1,024 receivers,
three displacement components, and all six tensor components; no compression
or reduced precision is introduced.

Production synthesis launches no `syn` processes. The C executable is retained
as a validation reference; it uses single precision and four-significant-digit
output, so float64 responses agree within its rounding tolerance rather than
bit-for-bit. Cache records identify the vectorized method and its numerical code;
old response/matrix caches must be regenerated. The slip path and C programs
retain their existing behavior.

One writer fills each complete source block and records completion. Failed
files stay marked incomplete and are rejected by matrix assembly. Named
libraries cannot be overwritten: choose a new `responses.name` to regenerate.
The library always stores all six basis responses, independent of look
coefficients, weighting, and full/fixed-tensor inversion choices.

Assembly reads source blocks, applies the supplied NEU coefficients, combines
components for fixed mechanisms if requested, and saves an **unweighted**
`GFs/matrices/<matrix.name>.npy`. Cached inversion loads that NumPy matrix once
for the whole damping sweep, bypassing response reads and assembly. Cache
dependency and product-hash validation still run. Changed LOS values, sigma,
solver, or damping reuse the same matrix. Changed projection coefficients or
mechanism require a new matrix but can reuse the physical library. Changed
receiver/source positions or Earth model require new physical responses.

Predictions are saved per damping value as
`output/inverse_models/statics/<run>.<number>.predictions.h5`, including IDs,
observations, coefficients, predictions, residuals, and run provenance. Tensor
models and logs retain their established outputs. The imagery variance-reduction
statistic retains MudPy's fractional reduction in the L1 residual norm.

The existing text-list path remains explicit with
`inputs.observation_format = "gflist"` (also the default for older TOMLs), with
the list in `data/station_info/`. It continues to support the existing numerical
reference examples. The slip-inversion workflow is unchanged.

Run the small HDF5 verification into a new destination:

```sh
PYTHONPATH="$PWD/src/python" python -B examples/geotomo/validate.py --hdf5 --project /absolute/existing-parent/new-project
```

Add `--processes 2` to also regenerate the small forward case in parallel and
require identical native GF files, text tables, HDF5 arrays, and projected matrix.
The validation directory must be new; this does not run a real-data inversion.

It runs the independent homogeneous reference, compares text/HDF5 matrices,
checks predictions and fixed mechanisms, and verifies NumPy-only reuse with
generation/assembly deliberately forbidden. `hdf5_validation.json` records
errors, code hashes, and small-case timings; those timings are not scaling claims.

## Forward profiling

For a forward-only performance check on synthetic data, use a new destination:

```sh
PYTHONPATH="$PWD/src/python" python -u -B examples/geotomo/validate.py --profile-forward --project /absolute/existing-parent/new-profile
```

This compares one and three workers for three sources and 1,507 synthetic
receivers. It times native FK generation and tensor conversion separately and
records parent/child CPU time, external call counts, and peak overlapping calls.
Production synthesis is required to launch zero `syn` processes. A separate
full-library write measures HDF5 cost. Serial/parallel response arrays and matrices
must agree exactly; sampled unit, isotropic, signed and mixed mechanisms are
checked against the original C `syn` at its output precision. `forward_profile.json`
records measurements and limitations; each run retains its TOML/provenance.
No inversion or Norcia data is used.

## Project initialization and named TOML runs

Each earthquake, volcano, or groundwater site remains one MudPy project, with
shared inputs and response libraries and multiple descriptive run names. Create
a project once from a short Python script or interactive session:

```python
from mudpy import geotomo

config = geotomo.init_project('/absolute/path/Norcia2016')
print(config)  # /absolute/path/Norcia2016/configs/template.toml
```

The destination must not exist. Initialization creates the standard MudPy tree,
including directories currently unused by static tomography, and adds `configs/`.
The four-column source template goes in `data/model_info/`, the HDF5 observation
example in `data/statics/`, and the velocity model in `structure/`. Text-list
examples are also retained for the supported text interface. All are explicitly
illustrative inputs.
The legacy `runslip.init(home, project_name)` retains its existing signature,
directory/template contents, and confirmation behavior; the shared directory
population helper avoids competing definitions of the tree.

Edit/replace project inputs, copy the starter to a descriptive experiment TOML,
and execute it from Python:

```python
from mudpy import geotomo

record = geotomo.run('/absolute/path/Norcia2016/configs/layered_fulltensor.toml')
print(record)  # .../output/inverse_models/models/<run_name>.resolved.json
```

Existing projects can use this entry point after adding a TOML; they do not need
initialization. The runner validates all required folders and primary inputs.
The launch interface is these two Python functions; there is no separate CLI.
Python before 3.11 requires the `tomli` reader (already present in the audited
Python 3.8 environment); Python 3.11+ uses `tomllib`. A missing reader raises an
import error rather than installing dependencies or selecting another format.

The commented [starter TOML](../../run/geotomo/template.toml) defines schema version
1. Required sections are `project`, `run`, `inputs`, `forward`, and `inversion`.
HDF5 inputs use `data/statics/`; optional `responses.name` identifies their
elastic library. Omitting `inputs.observation_format` selects `gflist` explicitly
as the established schema-1 default, rather than guessing from an extension.
`forward.earth_flattening`, `inversion.source_model`, and `inversion.damping` are
explicit. Optional values have recorded defaults: FK backend, `dk=0.05`,
`pmin=0`, `pmax=1`, `kmax=25`, `lstsq`, diagonal weighting, matrix assembly, and
one process. The default matrix basename follows `run.name`.

The TOML runner checks source depths against FK-disallowed layer interfaces
before launching workers or saving run products. The default
`forward.interface_policy = "error"` reports the source ID and interface depth.
To explicitly place interface sources on one side, use:

```toml
[forward]
earth_flattening = false
interface_policy = "below"  # "above" is also supported
interface_offset_km = 0.01  # 10 m from the interface, in physical depth
```

Sources away from interfaces are not shifted. FK depths are represented to
four decimal places (0.0001 km); half that increment is the interface tolerance.
Preflight uses FK's Vp-contrast test and accounts for its Earth-flattening
transformation when enabled. The rounded shifted point must lie strictly inside
the adjacent layer. Offsets below FK depth precision, shifts across another
boundary, and ambiguous hits on multiple interfaces fail explicitly. Above and
below select different material layers; this is a modeling choice, not a retry
strategy or an automatic choice of side.

The requested source file is preserved. When depths shift or require FK depth
rounding, the runner saves `data/model_info/_geotomo_<run>.sources`, reports each
change in the terminal, and records requested/effective coordinates and file
hashes in the resolved JSON and response-library provenance. Every calculation,
matrix check, cached solve, and saved model uses the effective geometry.
Physical caches include that geometry; changing sides/depths invalidates them.
Changing policy without changing effective geometry does not. Direct numerical
helpers do not apply this TOML preparation policy; the slip workflow is unchanged.

To parallelize forward calculations, set an explicit positive integer:

```toml
[execution]
processes = 4
```

This caps concurrent source workers for FK and vectorized make-synthetics;
it never selects all available cores automatically. Fewer sources use fewer
workers. Threads launch independent FK executables or evaluate NumPy source blocks,
so no MPI installation or multiprocessing entry-point guard is required.
Each FK source has a separate working folder. Batches contain at most this many
source blocks, and one calling thread writes HDF5 or text responses in source
order. Failed workers propagate their errors; active workers finish before the
project lock and temporary station input are released. Failed HDF5 libraries
remain marked incomplete and cannot be reused.

The requested count is saved in run provenance, but does not change the physical
response/matrix cache signatures. Solver damping sweeps remain sequential;
NumPy/BLAS thread settings are independent of `execution.processes`.

The project path resolves relative to the TOML file. Input filenames are
basenames in the fixed source/model/observation-list folders, respectively.
Observation-list data paths may be absolute or relative to the project root,
and, for text lists, must refer to `data/statics/<observation_id>.los`. The runner writes a
selected, absolute-path list at `data/station_info/_geotomo_<run_name>.gflist`;
the original observation list is not edited. Because that legacy format is
whitespace-delimited, project paths containing whitespace are rejected.

Configuration validation rejects unknown keys, unsafe names, incompatible
settings, malformed input arrays, invalid geometry/uncertainties, and unsupported
data modes. Source count is derived from the file. `full_tensor`
and `fixed_tensor` are supported; fixed mode requires six explicit components
in `xx,xy,xz,yy,yz,zz` order, with their supplied normalization preserved.
Only identity damping is currently supported. Damping may be an explicit list
or a logarithmic sweep:

```toml
[inversion]
source_model = "full_tensor"
solver = "lstsq"
weighting = "diagonal"
damping = { log10_min = -3.0, log10_max = -1.0, count = 15 }
```

`run.steps` is a nonempty, unique ordered subset of `greens`, `synthetics`, and
`invert`. A preparation-only run can request `greens` and `synthetics`; subsequent
named runs can request `invert`. Regenerating greens before inversion also
requires regenerating synthetics. These are fixed stages, not a general workflow
language.

For each run, files beside the existing numbered models/logs include:

| File | Contents |
| --- | --- |
| `<run>.config.toml` | Exact original bytes, including comments |
| `<run>.resolved.json` | Effective settings/defaults, expanded damping values, absolute paths, input hashes, source ordering, code identity, dependency versions, and executable identities |
| `<run>.0000.log`, etc. | The parameter combination and statistics for that solution, configuration-record paths/hashes, normalized observation-list path, and actual matrix hash |

Records are created after validation and before numerical execution. Changing
the editable TOML leaves archived records intact. A run name with existing
records/results is rejected. Computation errors propagate and the configuration
records remain available for the failed attempt; a normal return indicates that
the requested stages finished. Preserve the project inputs identified by the
hashes: configuration records do not duplicate large datasets. Primary inputs
must stay unchanged during a run; changes detected after computation raise an
error before stamping a cache as valid or returning success.

The runner holds the project lock across validation, archival, and all stages.
It also records cache identities alongside the existing GF tables and `.npy`
matrix. Synthetic reuse checks source/Earth-model/receiver identities, numerical
settings/code/executables. Text responses additionally depend on their baked-in
projection/mechanism; HDF5 physical responses do not. Matrix reuse checks
the same dependencies and the product hash. Observation values and uncertainty
changes are recorded as run inputs but do not invalidate a geometry-only matrix.
For an inversion using an existing verified matrix:

```toml
[run]
name = "different_damping"
steps = ["invert"]

[matrix]
name = "earlier_run"  # Existing GFs/matrices/earlier_run.npy
mode = "reuse"
```

An assembly target that exists is rejected. Missing, stale, or modified cache
records/products also raise errors; the runner never silently recomputes them.
Products without these records need an explicit rebuild through the runner.
Stage products are shared project resources and may
be regenerated by later named runs, while archived run configurations and
numbered results remain protected.

## Reproducing the reference

This example exercises the existing MudPy project formats, actual FK executables
and vectorized synthesis through an imagery-first workflow. It creates 24 receiver
locations with two representative LOS directions each and two fixed source
locations with intentionally reversed numeric IDs. All observations represent
one artificial deformation interval; these are not real satellite products.

Run from the MudPy checkout using the environment that imports this checkout:

```sh
python -B src/python/tests/test_geotomo.py
python -B examples/geotomo/validate.py --project /absolute/existing-parent/new-project
```

Exercise the same case through initialization and TOML, including preparation,
fixed/full-tensor runs, matrix reuse across different observation values, and
the single-receiver case:

```sh
python -B -m unittest discover -s src/python/tests -p 'test_geotomo*.py'
python -B examples/geotomo/validate.py --toml --project /absolute/existing-parent/new-project
```

The project path must not exist. Previous experiment files are never deleted.
FK's `fk.pl`, `st_fk`, and `trav` must be on PATH. `syn` is additionally required
for the C reference comparisons in `--hdf5` and `--profile-forward`, not for
production runs. The Python dependencies
are the existing NumPy/SciPy, ObsPy, and pyproj stack. No MPI or new package is
required. The committed Okada fixture removes the need to load an optional
reference solver during normal validation.

## Layered responses and off-grid recovery

Run the bounded layered/off-grid study in a new project:

```sh
python -B examples/geotomo/validate.py --layered --project /absolute/existing-parent/new-project
```

This mode always uses the TOML runner. Named experiments share one project,
receiver geometry, and verified response matrices. Core scientific algorithms
are the existing FK/full-tensor/identity-damped methods.
Run the driver as the sole owner of its new project: it updates shared inputs
between separately locked `geotomo.run` calls.

The experiment has 49 receiver locations on a 7-by-7 grid spanning north/east
coordinates from -6 to +6 km, with two representative LOS directions per point.
The homogeneous control uses `Vp=6 km/s`, `Vs=6/sqrt(3) km/s`, and density
`2.7 g/cm³`. Its identical-layer counterpart introduces a boundary at 2 km
without a material contrast. The contrasted model instead has a 2 km upper
layer with `Vp=4.3 km/s`, `Vs=2.4 km/s`, and density `2.4 g/cm³` over the
homogeneous half-space. Earth flattening is disabled throughout.

Forward checks include all six components at two sources, homogeneous-layer
equivalence, integration-spacing values `dk=[0.2,0.1,0.05,0.025]` at `kmax=35`,
and cutoff values `kmax=[15,25,35]` at `dk=0.025`. The finest computation is a
numerical reference, not an exact solution. The identical-layer limit and the
two finest spacings/cutoffs must agree within an initial `1e-3` relative L2
bound for every matrix column. Direct prescribed shear, positive/negative
isotropic, and mixed tensors additionally test superposition and sign reversal
in the genuinely layered medium. This is a pipeline/convergence check; an
independent layered-Earth physical benchmark remains outstanding.

Off-grid truth locations are `[-0.4,0.35,3.6]` and `[0.55,-0.25,4.4]` in
north/east/depth km. Inversion candidates lie at north/east coordinates ±1 km
and depths 3 and 5 km: eight sources and 48 full-tensor unknowns. IDs are
deliberately reverse ordered. None coincides with a true location.

Pure shear, inflation, and deflation activate only the first truth location;
the mixed case has shear at the first and isotropic expansion at the second.
Pure-source scalar moment is `10**16.6 N m`; the mixed isotropic contribution
is 0.6 times that value. Truth remains integrated point moment, not moment
density, pressure, or physical volume change.

Each case is run noiselessly and with one reproducible iid Gaussian noise
vector, sigma 0.3 mm and seed 20261106. The same realization is reused across
mechanisms for a controlled comparison. This is not a realistic correlated
InSAR uncertainty model. The bounded damping sweep is `[0,0.01,0.1,1,10]`, using
the existing raw six-component identity penalty and weights `1/sigma`.
No automated hyperparameter selection is performed. The example-map damping
value of 1 is fixed before computing recovery errors.

The known-grid control checks algebraic prediction consistency. For every
solution, the example verifies source IDs, recomposition of saved predictions
from saved tensors, stationarity of the stated quadratic objective, and the
expected increase in data residual/decrease in model norm with damping.
Noiseless inflation/deflation solutions must have opposite signs. These
acceptance checks do not require a good off-grid reconstruction: ambiguity,
excess inferred moment, and source spreading are scientific findings.

Recovery diagnostics include clean-displacement error, weighted data RMS,
scalar-moment-weighted centroid and spread, aggregate tensor error and angle,
and summed scalar moment relative to truth. The `cancellation_ratio` field is
`sum(M0_i)/M0(sum(M_i))`: it reflects differing tensor orientations as well as
opposing-source cancellation, and must be interpreted against the known truth.
Aggregate tensor/centroid agreement does not establish correct mechanisms at
individual source locations. Zero models have undefined centroid, spread,
angle, and scalar-to-net ratio, recorded as JSON nulls. The configured truth
models all have nonzero aggregate tensors.

Outputs include:

| Location | Result |
| --- | --- |
| `analysis/layered_offgrid.json` | Forward acceptance results, every recovery metric, condition number, noise specification, driver hash, timings, and limitations |
| `analysis/validate.executed.py` | Exact driver snapshot corresponding to `driver_sha256` |
| `analysis/<forward-run>.G.npy` | Preserved matrices for numerical comparisons after shared text tables change |
| `analysis/<inversion-run>.data.npy` | Exact observation vector used by that run |
| `analysis/synthetic_truth.npz` | Source/receiver coordinates, LOS directions, prescribed tensors, and noise realization |
| `analysis/grid_prepare.prof`, `analysis/mixed_noisy.prof` | Standard-library profiles of library construction and an inversion sweep |
| `plots/layered_convergence.{png,pdf}` | Integration-spacing/cutoff checks |
| `plots/layered_response.{png,pdf}` | Homogeneous/layered displacement and their difference |
| `plots/offgrid_tradeoffs.{png,pdf}` | Fit, centroid, aggregate tensor, and total-moment behavior |
| `plots/offgrid_example.{png,pdf}` | Mixed noisy observations, prediction/residual, and inferred source moments |

Working `.los` files contain the last dataset written. Earlier inputs are
preserved in the per-run observation arrays; the common receiver list and
`synthetic_truth.npz` provide the ordering/directions needed to reconstruct
their `.los` files with the driver's `%.12e` formatting. Re-running the whole
example into a fresh project reproduces the deterministic sequence without
overwriting earlier named results. These datasets use the same FK backend and
Earth model for truth and inversion, with different source discretizations;
they do not test Earth-model error, finite-source mismatch, real imagery,
correlated atmospheric errors, or topography.

The execution driver hash identifies the numerical run. If figures are later
rendered again after layout-only edits, retain that snapshot/hash and record
the plotting driver's hash separately as `plot_driver_sha256` in the summary.

## Homogeneous reference checks

- Generate both DC and explosion GFs, with Earth flattening disabled for the
  homogeneous flat half-space benchmark. In the API, `flat=True` retains the
  historical `/f` Earth-flattening option; use `flat=False` for this benchmark.
- Check isotropic displacement against a closed-form half-space solution and
  five double couples spanning the deviatoric tensor space against Okada.
- Check fixed-mechanism responses against combinations of six basis responses.
- Recover prescribed full tensors for mixed sources, inflation, and deflation;
  also recover signed fixed-mechanism amplitudes.
- Verify source IDs, physical scalar moment, cached-array equality, and saved
  predicted LOS values. Unit tests cover incomplete/duplicate responses, one
  observation, tensor norms, nodal planes, and the weighted/damped NNLS objective.

Forward comparisons allow a relative L2 error of `1e-3` over the observation
vector. This is an explicit initial acceptance bound, allowing for FK integration
error. Separate C `syn` comparisons allow `6e-4` relative L2 for its single
precision and four-significant-digit output. Same-matrix noiseless
recovery uses `1e-7`; that tight algebraic check does not measure real source
resolution. Agree and revisit scientific tolerances as models/data broaden.

`validation.json` records numerical results, limitations, Python/NumPy versions,
and hashes of the four modified numerical modules. The project retains inputs,
GF tables, matrices, model/log files, and predicted observations.

## Conventions and independent references

Tensor coordinates are north/east/down; stored component order is
`xx, xy, xz, yy, yz, zz`. Displacement and projection vectors use north/east/up.
Coordinates are longitude/latitude in degrees; source depth is positive in km.
Scalar observations are in meters, positive along the supplied unit vector.
The synthetic vectors are representative geometries rather than actual orbit
metadata. Legacy `.los` files store `[displacement, north, east, up]`.

Each coefficient scales a tensor component by `10**16.6 N m`, equivalent to
`10**23.6 dyne cm` in the C reference. This component scale is not a common scalar moment for
every unit basis tensor. For an isotropic moment `M=m I`, the independent surface
reference used here is

```text
u_NEU = (1 - 2*nu)*m / (2*pi*mu) * [north, east, depth] / R**3
R**2 = north**2 + east**2 + depth**2
```

It follows from the Mogi displacement with equivalent cavity volume change
`dV=m/(lambda+2*mu)` and all distances in meters. This equivalence does not turn
moment into magma volume, stress-free volume change, or a layered-Earth cavity
volume. The [USGS deformation-model reference](https://pubs.usgs.gov/tm/13/b1/)
provides the analytical modeling context.

The frozen `okada_reference.json` was generated from the locally available
Fortran DC3D0 implementation. Its library and source hashes are recorded. Before
generating the fixture, DC3D was checked against the published manual example,
yielding `[-37.8981438, 63.1789360, 14.9606857]`, consistent with the tabulated
answer. DC3D0 inputs use `alpha=2/3`, source depth 3000 m, and strike/dip potency
`M0/mu`. Local x follows strike, local y points to its left, and z is up.
The five recorded mechanisms check off-diagonal signs and diagonal differences.
See the [NIED DC3D/DC3D0 manual](https://www.bosai.go.jp/information/pdf/DC3Dmanual.pdf).

## Scope and ownership

Forward source calculations accept a positive `ncpus` count (or TOML
`execution.processes`), with serial execution as the default. Unsupported
GNSS/waveform requests fail explicitly. The standard slip-GF path has
not been changed to use the new tensor branch. No FK source code was modified.

`inversionGFs` and `run_inversion` acquire the exclusive `.geotomo.lock` project lock. A second
top-level run in that project fails before mutation. A killed process can leave
the lock behind; inspect it and the failed run before manually removing it.
Direct helper calls require exclusive caller ownership of the project, including
rank input files, response tables, and matrix caches. No process working-directory
changes are made in the repaired tensor path.

Direct legacy `getG` calls still lack cache provenance checks: after a source,
observation geometry, mechanism, or Earth-model change, rebuild with
`G_from_file=False`. The TOML runner checks its cache sidecars before reuse.
Old ten-column DC-only GF files are rejected; regenerate complete fourteen-column
DC+EX responses. No storage backend or general observation schema has been adopted.

This example checks two LOS views. Simple horizontal offset-direction projections
are also checked, but production radar/optical pixel-offset ingestion remains to
be designed. The layered mode adds numerical convergence, off-grid point sources,
and iid noise. Independent layered-Earth benchmarks, finite-source mismatch,
correlated noise, nuisance terms, real imagery, MPI, and numerical topography
remain later validation work. Pure isotropic/zero tensors have undefined DC nodal planes,
which are written as NaNs; other reported planes describe a best fitting DC only.
