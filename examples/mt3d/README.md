# Tensor imagery CPU reference

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

Legacy `.npy` caches still lack provenance checks: after a source, observation,
geometry, mechanism, or Earth-model change, rebuild with `G_from_file=False`.
Old ten-column DC-only GF files are rejected; regenerate complete fourteen-column
DC+EX responses. No storage backend or general observation schema has been adopted.

This example checks two LOS views. Simple horizontal offset-direction projections
are also checked, but production radar/optical pixel-offset ingestion remains to
be designed. Layered-Earth benchmarks, finite-source/grid mismatch, correlated
noise, nuisance terms, real imagery, MPI, and numerical topography remain later
validation work. Pure isotropic/zero tensors have undefined DC nodal planes,
which are written as NaNs; other reported planes describe a best fitting DC only.
