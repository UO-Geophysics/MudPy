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

Native FK intermediates are retained by default. For HDF5 runs, opt in to
cleanup after successful synthesis with:

```toml
[execution]
processes = 4
keep_fk_intermediates = false
```

Cleanup runs at the end of a successful run while holding the project lock.
It verifies the closed HDF5 library's schema, completion flags, source/receiver
geometry, float64 shape and finite source blocks, plus product/cache hashes.
It then removes only the listed native `.insar` files and matching source
folders' copied model and `junk.p`/`junk.s` files. No recursive/glob deletion is
used. Unexpected files, symlinks or changed working models stop cleanup before
any deletion. Calculation/solver failures retain the intermediates, and a
greens-only run keeps them until a later successful HDF5 synthesis.

The response-library `.h5.json` sidecar records a cleanup receipt with file
hashes, paths, bytes, run name and completion time. `_geotomo_greens.json` marks
native availability as `pruned`. HDF5 libraries, NumPy matrices, project inputs
and run records are preserved; new matrix assembly and cached inversion need
no native files. Rebuilding the response library requires explicitly rerunning
`greens` after pruning. An interrupted filesystem cleanup remains marked
`cleanup_pending`, preventing native-cache reuse; its verified HDF5 library
remains available. Retention policy does not change physical cache signatures.
Text-list and legacy slip runs retain their intermediates.

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
It also exercises actual native-file cleanup, then checks fresh matrix assembly
and cached inversion against the saved reference after those files are gone.

## Quick result viewing

Reusable geotomo viewers live in `mudpy.view`. They read the canonical archived
models/logs/prediction HDF5 files, return `(figure, axes)` with a flat axes array,
and do not show, save or modify inversion results automatically:

```python
from mudpy import view

project = '/absolute/path/Norcia2016'
run_name = 't44_aci_volume_l2_v1'

fig, axes = view.geotomo_tradeoff(project, run_name, reference_mw=6.6)
fig.savefig('tradeoff.png', dpi=180)

fig, axes = view.geotomo_fit(project, run_name, 9, units='mm',
                             limits=(-250, 250), residual_limits=(-30, 30))

fig, axes = view.geotomo_sources(project, run_name, 9)
fig, axes = view.geotomo_sources(project, run_name, 9, kind='sections',
                                 center=(13.10, 42.83), strike=155)
```

These example reference values and member numbers are caller choices, not viewer
defaults or preferred solutions. The sweep viewer plots the coefficient/data L2
tradeoff, unweighted residual RMS, coefficient norm and moment diagnostics.
`highlight=[...]` labels particular members. RMS is computed directly from saved
residuals; the logged legacy "VR" statistic is not labelled as variance reduction.

The imagery viewer plots actual samples without interpolation. Observed and
predicted panels share limits; separate residual limits are explicit and colourbar
extensions indicate clipping. `units` is `m` or `mm`; `label` can describe another
scalar projection/offset. Source plots use point-source scalar moments in Nm.
Depth views share a normalization across depth levels. Sections sum scalar moments
at coincident projected nodes (1 mm numerical grouping) across the other horizontal
coordinate. They accept irregular source clouds and do not assume a fault mesh.
Pass a Matplotlib `norm` to compare solutions on a common linear/log scale.

All viewers accept `axes` for caller-owned layouts and preserve an
existing figure title. Standalone figures use local styling. Existing slip viewers
retain their interfaces and behavior. Project analysis scripts choose members,
reference values, map extents, shared colour scales, composites and export paths;
they can save PNG/PDF or call `plt.show()` as appropriate.

`geotomo_io.read_model(path)`, `read_predictions(path)`,
`read_solution(project, run_name, member)` and `read_run(project, run_name)` expose
arrays/dictionaries for custom analysis. Readers validate geometry, finite tensors,
scalar moments, observed-minus-predicted residuals, archived damping and log norms;
missing/inconsistent members fail. Undefined nodal planes may be NaN. Moments are
explicitly named `sum_scalar_moment_nm` and `net_scalar_moment_nm`: the first sums
point-source tensor norms, while the second is the scalar norm of their summed
tensor. These are different diagnostics; neither automatically selects a catalog-
consistent solution. `reference_mw` is optional and converted with MudPy's Mw
convention. The inversion's raw-coefficient L2 penalty is distinct from these
rotation-invariant scalar-moment norms.

Full source hashes remain in run provenance. Response/matrix cache identity uses
the numerical portion of `geotomo_io.py`, ending at its explicit result-reader
boundary. Reader-only changes preserve numerical caches; production storage/
assembly code belongs above the boundary. No cache records are migrated or
rewritten by plotting.

### Source components and mechanisms

```python
from mudpy.analysis import decompose_geotomo_tensors
from mudpy.geotomo_io import read_solution

solution = read_solution(project, run_name, 9)
parts = decompose_geotomo_tensors(solution['model']['tensors_nm'])

# Each quantity shares a colour normalization across its depth panels.
fig, axes = view.geotomo_sources(project, run_name, 9, quantity='iso')
fig, axes = view.geotomo_sources(project, run_name, 9, quantity='iso_balance')
fig, axes = view.geotomo_mechanisms(project, run_name, 9, depths=[4])
fig, axes = view.geotomo_mechanisms(project, run_name, 9,
                                  component='dc', style='circles')
fig, axes = view.geotomo_mechanisms(project, run_name, 9, component='dc',
                                  style='circles', color_by='rake', plane=1)
fig, axes = view.geotomo_mechanisms(project, run_name, 9, component='dc',
                                  kind='section', center=(13.10, 42.83), strike=155,
                                  section_offset_km=8, section_width_km=3)
fig, axes = view.geotomo_orientation(project, run_name, 9,
                                   center=(13.10, 42.83), strike=155,
                                   section_offset_km=8, section_width_km=3)
```

The decomposition accepts batches of six NED components in the saved model order:
`xx,xy,xz,yy,yz,zz`. It returns full, ISO, DEV, DC and CLVD tensors, their scalar
moments, signed ISO/CLVD strengths, ISO balance, DC planes/axes and mechanism
colours. ISO is `trace(M)/3 I`; DEV = DC + CLVD uses the conventional
[Jost–Herrmann decomposition](https://www.seismosoc.org/Publications/srl/SRL_83/srl_83-3_es_I/).
The existing `MT.compute_DC` is a best-fitting DC representation and retains its
behavior; saved best-fitting planes do not measure the standard decomposition's
DC amount. Pure CLVD has no DC mechanism under this decomposition.

Every component strength uses `sqrt(M:M/2)` in Nm, consistent with existing
geotomo moments. Component scalar norms are **not additive percentages**. ISO
sign follows trace (positive expansion, negative contraction); CLVD sign follows
its unique principal axis, positive for eigenvalues proportional to `[-1,-1,2]`.
CLVD is trace-free, so its sign does not mean inflation/deflation. Pure shear on
a curved fault remains locally DC; mixed ISO/DEV sources are also physically
possible, including opening/closing cracks.

Depth-map `quantity` accepts `moment` (default), `iso`, `dev`, `dc`, `clvd` or
`iso_balance`; `depths=[...]` selects existing model levels. ISO/CLVD defaults
are symmetric diverging scales. Balance is `M0(ISO)/(M0(ISO)+M0(DEV))`, fixed to
0–1, with symbol area proportional to full moment. A zero tensor has undefined
balance and is marked neutrally. Balance is a defined index, not a fraction of
the original norm or an admissibility test for physical sources.

Mechanism beachballs show the full tensor by default, or the decomposed DC.
They are lower-hemisphere, north-up focal-sphere symbols (black positive, white
negative radiation), drawn with ObsPy/MoPaD; they are not observed seismic data.
Symbol area scales with the selected component moment. `moment_reference_nm`
sets the moment represented by `diameter` in points; the default reference is
the maximum over the whole model, keeping sizes comparable across depths.
Zero/undefined mechanisms are neutral crosses. DC orientations below 64 machine
epsilons of full moment are numerical undefined values.

DC circles blend RGB = `(T_down**2, B_down**2, P_down**2)`: reverse red,
strike-slip green, normal blue, with intermediate colours representing oblique
axis orientations. They summarize principal-axis character, not a chosen fault
plane. Rake colouring uses a cyclic -180–180 degree scale and requires an explicit
`plane=1` or `2`; neither algebraic alternative identifies the physical fault.

Mechanism `kind='section'` selects an explicit along-strike strip and projects
its individual source positions across strike; it never sums tensors. Exact
projected overlaps (1 mm grouping) fail and require a narrower strip. Section
beachballs retain their **north-up plan-view orientation**, rather than rotating
the focal sphere into the section. `geotomo_orientation` shows both plane
alternatives' strike/dip/rake versus depth, alongside source position, with
moment-scaled symbols and an optional matching strip. No plane selection, depth
jitter or fitted flattening trend is imposed.

Tests cover pure ISO/DC/CLVD, signs, reconstruction, rotations, extreme scales,
zero/undefined sources, plane reconstruction, symbol sizing and section selection.
Independent comparison of all 315 Norcia member-9 tensors with Pyrocko's standard
decomposition agrees within 2e-15 relative tensor norm. The installed legacy MoPaD
has an inconsistent full/deviatoric eigensort for some mixed ISO tensors; regression
comparisons therefore use its trace-free DC/CLVD decomposition. Production
decomposition uses NumPy eigenanalysis; Pyrocko is not a runtime requirement.

## Pure double-couple inversion

The canonical HDF5 runner supports a pure-DC nonlinear least-squares mode:

```toml
[inversion]
source_model = "double_couple"
solver = "least_squares"
weighting = "none"
damping = [0.001]  # Illustrative; choose a fresh sweep for this penalty.

[inversion.double_couple]
initial_catalog = "event.ndk"
starts = 3
seed = 0
perturbation_deg = 30.0
max_nfev = 500
ftol = 1e-9
xtol = 1e-9
gtol = 1e-9
```

Place a **single-event GCMT NDK** file in `data/model_info`. The reader rejects
multiple events/mechanisms and partial parses. ObsPy converts NDK dyne-cm to Nm;
the runner explicitly converts USE coordinates into NED. The original catalog
tensor, event identity, normalized DC initialization and file hash are archived.
Alternatively supply `initial_tensor = [xx,xy,xz,yy,yz,zz]` in this table, using
six NED components with any common scale. Exactly one initialization is required.
Its standard DC component must be defined. Only its orientation initializes the
sources; the catalog event's moment is not assigned to every source or imposed
as a total-moment constraint.

Each source has a nonnegative amplitude `a = M0 / 10**16.6 Nm` and three rotation
coordinates. Rotating a normalized DC tensor preserves eigenvalues `[-a,0,a]`,
so ISO and CLVD remain zero throughout optimization. By default, the solver
fits each source mechanism independently. The first start uses the catalog DC;
additional starts independently perturb each source by seeded rotations whose
angles are uniform between zero and `perturbation_deg`.

The default objective is `||W(Gm(a,r)-d)||² + lambda² ||a||²`. This is nonlinear least
squares with **rotation-invariant scalar-moment damping**, unlike the existing
raw-component penalty. Retune damping rather than transferring a member's old
lambda. For each start, direct NNLS on the damped system initializes amplitudes;
a numerical floor of 1% of the initial amplitude scale frees rotations at zero
nodes. Final amplitudes have lower bound zero, without that initialization floor.
An exactly zero dataset has the explicit zero-moment global minimum.

Optional horizontal mechanism smoothing adds one weight, `horizontal_smoothing`,
under `[inversion.double_couple]`. It defaults to zero. When positive, set
`horizontal_neighbor_radius_km` to define the graph from source geometry. For
example, 4.5 km includes the side neighbors of a 4 km square grid and excludes the
diagonals. Pairs must share the same requested depth, before any FK adjustment,
and be within the WGS84 surface distance cutoff. The resolved record lists
the paired source IDs. A positive weight with no eligible pairs raises an error.

For each pair, the added residual is the smoothing weight times the difference
of the two unit-moment DC tensors, weighted by the harmonic mean
`2*a_i*a_j/(a_i+a_j)` of the normalized moments. The squared tensor norm is
`D:D/2`. Equal mechanisms incur zero penalty at any amplitudes, and a zero
moment disables that pair's penalty. Scalar-moment damping retains its own
weight. The initial NNLS solve still uses data and scalar-moment damping, and
the nonlinear solver then includes all three terms. Diagnostics record the
neighbor count, cutoff, weight, and each start's `horizontal_mechanism_l2` before
multiplication by the smoothing weight.

SciPy TRF/LSMR uses analytic rotation derivatives and a matrix-free Jacobian.
The data residual and regularization terms are divided by the data L2 norm for
convergence scaling, which preserves the objective minimizer. Saved residuals, costs
and moment norms use the original units. `scaled_optimality` refers to the scaled
objective. `max_nfev` applies to each start at each damping value. Starts execute
sequentially with local state; ordinary BLAS threading remains independent.

Only converged starts can be selected. Every member gets a `.solver.json` with
status, message, evaluations, cost, scaled optimality and moment/residual norms
for all starts, plus relative tensor/prediction differences between converged
starts and the selected model. These describe local-solution sensitivity, not a
global-optimality guarantee. If no start converges, failure diagnostics are saved,
the run raises, and that member receives no normal model/prediction output.
The lock is released and any requested FK cleanup does not run after failure.

Diagnostics also record wall/CPU time for each start, initialization and nonlinear
optimization wall time, and total solver wall/CPU time per damping value. CPU time
includes all BLAS threads. A full-volume Norcia attempt exposed the former bounded
TRF amplitude initializer's iteration limit; direct NNLS solves that same convex
subproblem and is checked for nonnegativity and KKT optimality on correlated
source responses. Nonlinear stopping criteria remain explicit TOML choices.

The existing 17-column model, log and prediction formats are retained. Result
readers verify the diagnostic checksum, selected-start norms and pure-DC tensor
eigenvalues. The tradeoff viewer labels the new moment-amplitude norm explicitly.
Full-tensor/fixed-tensor and legacy slip inversion paths retain their behavior.

The solver lives in `geotomo_dc.py`, separately fingerprinted in run provenance.
The forward/linear implementation is unchanged. Full-tensor and pure-DC modes
use the **same six-component matrix basis**, with identical physical cache
signatures; changing initialization or solver controls does not invalidate GFs
or that matrix. A fixed-tensor amplitude matrix remains a different basis.
Existing cache records are neither migrated nor rewritten on reuse.

Controlled validation (no real-data inversion):

```sh
PYTHONPATH="$PWD/src/python" python -B examples/geotomo/validate_dc.py --output /absolute/new-validation-directory
# Optionally pass --matrix /absolute/synthetic-six-component-matrix.npy.
```

The numerical `geotomo_dc.solve` API also accepts `moment_penalty="l1"`.
It minimizes `||Gm-d||² + damping*sum(a) + horizontal_smoothing²*S`
with nonnegative moments and analytic L-BFGS-B gradients. L1 uses the damping
weight directly. L2 uses its square. The amplitude coordinates are scaled
before optimization. Diagnostics identify the optimizer, moment L1 norm,
iterations, evaluations, exact active count and stopping condition.

Full-volume sparse synthetic validation uses an existing Norcia six-component
matrix and its checked project configuration. This experiment plants three
normal sources among all 315 candidates, adds 2 mm LOS noise, tests horizontal
smoothing, and repeats with an oblique starting mechanism. Fixed-mechanism
L-BFGS-B is checked against scikit-learn Lasso's objective and dual gap.
The numerical experiment uses two starts per case.

```sh
PATH="$PWD/src/fk:$PATH" PYTHONPATH="$PWD/src/python" MPLBACKEND=Agg \
python -u -B examples/geotomo/validate_dc.py \
  --sparse-config /absolute/Norcia2016/configs/t44_aci_volume_dc_gcmt_v4.toml \
  --output /absolute/existing-parent/new-sparse-validation
```

This command requires scikit-learn for the independent convex reference. It
preserves arrays, metrics, input/code fingerprints and PNG plots. Source counts
include both exact nonzero moments and the number carrying 95% of total moment.
Canonical TOML runs also accept `moment_penalty = "l1"` in
`[inversion.double_couple]` with `solver = "lbfgsb"` in `[inversion]`.
Existing DC TOMLs default to the L2/TRF solver. L1 settings, solver diagnostics,
model/prediction checksums and per-start timings are archived with each run.
The tradeoff viewer uses the recorded L1 moment norm for these runs. The log's
legacy `Lm` field retains its L2 meaning for existing readers.

The driver tests noiseless/noisy changing-dip sources and damping, records code/
input identities, and refuses an existing output directory. Numerical and TOML
integration tests cover derivative/adjoint checks, catalog conventions, weighted
sweeps, full-matrix reuse, diagnostics, failed convergence and read-only viewing.
Small known-location recovery does not establish resolution or performance for a
large source volume. Norcia pure-DC testing is a separate, deferred step.

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
