# Tensor imagery CPU reference

## Project initialization and named TOML runs

Each earthquake, volcano, or groundwater site remains one MudPy project, with
shared inputs and response libraries and multiple descriptive run names. Create
a project once from a short Python script or interactive session:

```python
from mudpy import mt3d

config = mt3d.init_project('/absolute/path/Norcia2016')
print(config)  # /absolute/path/Norcia2016/configs/template.toml
```

The destination must not exist. Initialization creates the standard MudPy tree,
including directories currently unused by static tomography, and adds `configs/`.
The four-column source template goes in `data/model_info/`, the imagery list in
`data/station_info/`, its example displacement/projection file in `data/statics/`,
and the velocity model in `structure/`. All are explicitly illustrative inputs.
The legacy `runslip.init(home, project_name)` retains its existing signature,
directory/template contents, and confirmation behavior; the shared directory
population helper avoids competing definitions of the tree.

Edit/replace project inputs, copy the starter to a descriptive experiment TOML,
and execute it from Python:

```python
from mudpy import mt3d

record = mt3d.run('/absolute/path/Norcia2016/configs/layered_fulltensor.toml')
print(record)  # .../output/inverse_models/models/<run_name>.resolved.json
```

Existing projects can use this entry point after adding a TOML; they do not need
initialization. The runner validates all required folders and primary inputs.
The launch interface is these two Python functions; there is no separate CLI.
Python before 3.11 requires the `tomli` reader (already present in the audited
Python 3.8 environment); Python 3.11+ uses `tomllib`. A missing reader raises an
import error rather than installing dependencies or selecting another format.

The commented [starter TOML](../../run/mt3d/template.toml) defines schema version
1. Required sections are `project`, `run`, `inputs`, `forward`, and `inversion`.
`forward.earth_flattening`, `inversion.source_model`, and `inversion.damping` are
explicit. Optional values have recorded defaults: FK backend, `dk=0.05`,
`pmin=0`, `pmax=1`, `kmax=25`, `lstsq`, diagonal weighting, matrix assembly, and
one process. The default matrix basename follows `run.name`.

The project path resolves relative to the TOML file. Input filenames are
basenames in the fixed source/model/observation-list folders, respectively.
Observation-list data paths may be absolute or relative to the project root,
and must refer to `data/statics/<observation_id>.los`. The runner writes a
selected, absolute-path list at `data/station_info/_mt3d_<run_name>.gflist`;
the original observation list is not edited. Because that legacy format is
whitespace-delimited, project paths containing whitespace are rejected.

Configuration validation rejects unknown keys, unsafe names, incompatible
settings, malformed input arrays, invalid geometry/uncertainties, and unsupported
parallel or data modes. Source count is derived from the file. `full_tensor`
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
settings/code/executables, and projection/mechanism identity. Matrix reuse checks
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
Legacy products without these records need an explicit rebuild through the
runner. These JSON sidecars document existing formats; they do not introduce a
new array storage backend. Stage products are shared project resources and may
be regenerated by later named runs, while archived run configurations and
numbered results remain protected.

## Reproducing the reference

This example exercises the existing MudPy project formats and actual FK/`syn`
executables through a serial, imagery-first workflow. It creates 24 receiver
locations with two representative LOS directions each and two fixed source
locations with intentionally reversed numeric IDs. All observations represent
one artificial deformation interval; these are not real satellite products.

Run from the MudPy checkout using the environment that imports this checkout:

```sh
python -B src/python/tests/test_mt3d.py
python -B examples/mt3d/validate.py --project /absolute/existing-parent/new-project
```

Exercise the same case through initialization and TOML, including preparation,
fixed/full-tensor runs, matrix reuse across different observation values, and
the single-receiver case:

```sh
python -B -m unittest discover -s src/python/tests -p 'test_mt3d*.py'
python -B examples/mt3d/validate.py --toml --project /absolute/existing-parent/new-project
```

The project path must not exist. Previous experiment files are never deleted.
FK's `fk.pl`, `st_fk`, `trav`, and `syn` must be on PATH. The Python dependencies
are the existing NumPy/SciPy, ObsPy, and pyproj stack. No MPI or new package is
required. The committed Okada fixture removes the need to load an optional
reference solver during normal validation.

## Checks

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
and the existing `syn` four-significant-digit text output. Same-matrix noiseless
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

Each coefficient scales a tensor component by `10**16.6 N m`, passed to `syn`
as `10**23.6 dyne cm`. This component scale is not a common scalar moment for
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

The first implementation requires `ncpus=1`. Unsupported GNSS/waveform requests
and parallel tensor execution fail explicitly. The standard slip-GF path has
not been changed to use the new tensor branch. No FK source code was modified.

`inversionGFs` and `run_inversion` acquire an exclusive project lock. A second
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
be designed. Layered-Earth benchmarks, finite-source/grid mismatch, correlated
noise, nuisance terms, real imagery, MPI, and numerical topography remain later
validation work. Pure isotropic/zero tensors have undefined DC nodal planes,
which are written as NaNs; other reported planes describe a best fitting DC only.
