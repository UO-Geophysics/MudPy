# Native-grid InSAR preprocessing

This general MudPy workflow is independent of geodetic tomography. The functions
are in `mudpy.insartools`. Existing Matlab workflows remain unchanged. No Kite,
Shapely, raster-processing installation, or GUI is required by the array-based
core: it uses NumPy, with Matplotlib used only for diagnostics. The optional
recovered-Norcia reader uses Rasterio.

## Processing contract

`preprocess_insar(displacement, x, y, look_enu, ...)` takes a 2-D native raster,
x/y coordinate axes or matching coordinate arrays, and `(rows, columns, 3)` ENU
look vectors. NaN/masked pixels are missing; zeros remain measurements. A finite
missing-data sentinel must be supplied explicitly as `nodata`. An optional
boolean `quality_mask` excludes unsuitable observations and fitting pixels.

The caller supplies `units` (`m`, `cm`, `mm`), `los_positive` (`toward`, `away`,
`unknown`), `look_direction` (`ground_to_satellite`, `satellite_to_ground`,
`unknown`), and a nonempty `coordinate_system` description. Known conventions
are standardized to meters, positive toward satellite, with ground-to-satellite
vectors. Unknown directions stay marked unknown and are not guessed. Input
vectors used for observations must meet the stated unit-norm tolerance; none
are silently normalized. No original arrays or source files are modified.

**An output labeled unknown is suitable for diagnostic processing, not for
inversion.** Resolving sign, reference, acquisition history, quality and noise
assumptions remains the user's responsibility. This module does not run an
inversion or automatically export a MudPy project.

## Independent polygon roles

- `include_polygons`: the observation ROI. Only its pixels enter tile statistics
  and means, before downsampling. None retains the whole scene.
- `exclude_polygons`: deformation or other areas excluded from the background
  trend fit. Outside the union of these polygons is eligible for fitting.
- `background_polygons`: an optional explicit fitting region. None allows the
  whole scene outside the excluded polygons, including outside the observation ROI.

Quality/missing LOS masks apply to fitting. Look geometry is required for
sampled observations, but a valid LOS background pixel can still constrain the
surface even when its geometry is missing. Thus a small observation ROI or
partial geometry coverage does not silently eliminate the far field.

Polygons may be Nx2 vertex arrays, a sequence of such polygons (union), or
GeoJSON Polygon/MultiPolygon/Feature/FeatureCollection dictionaries. GeoJSON holes
are supported. Polygon boundaries, including hole boundaries, belong to the
polygon. Empty, self-intersecting and invalid rings raise. Coordinates and all
polygons must be in the same CRS; no implicit reprojection is performed.

## Optional background correction

`trend=None` (default) preserves the displacement reference. Other choices are:

- `offset`: constant.
- `plane`: constant, x and y.
- `quadratic`: constant, x, y, x squared, x*y and y squared.

`fit_surface` is also independently callable. It fits full-resolution background
LOS pixels before quadtree sampling. Coordinates are centered/scaled for numerical
conditioning, and their origin/scale and every polynomial coefficient are saved.
The basis uses the supplied coordinates; use local projected east/north if a
polynomial in physical distances is desired. No data are reprojected by the module.
Optional `fit_weights` are nonnegative least-squares weights; zero removes a
pixel from fitting. No robust fit or polynomial-order fallback is implicit.
Too few points, rank deficiency, or invalid weights raise an error.

The fitted surface is subtracted across the valid scene. A constant term changes
the displacement reference. A plane/quadratic can remove real broad deformation
even when fitted outside a polygon: the polygon states a scientific assumption,
not proof of non-deforming ground. Preserve and inspect the original, fit mask,
removed surface and corrected scene; compare reasonable background choices.
Surface-fit uncertainty and spatial noise covariance are not estimated here.

## Quadtree behavior and geometry

`quadtree` is independently callable. It uses native pixels with no griddata,
gap filling, smoothing or forced 2048x2048 interpolation. Rectangular and odd-sized
rasters use a conceptual power-of-two tree clipped at their actual edges.

- `threshold` is in meters in `preprocess_insar`, or in the supplied displacement
  units when calling `quadtree` directly.
- `criterion='range'` uses maximum minus minimum, matching the basic floating-point
  MATLAB `qtdecomp` criterion. `criterion='std'` is a different, explicitly selected
  population-standard-deviation criterion. Their thresholds are not interchangeable.
- `min_size`/`max_size` are powers of two in native pixels. Clipped edge tiles may
  be smaller. ROI boundary tiles split to minimum size.
- `min_valid_fraction` is valid/selected pixels within a tile; pixels outside the
  ROI do not count as missing. `min_pixels` gives a minimum contributing count.
  Low-coverage tiles split to recover usable subregions; unsuitable final tiles
  are dropped. Leaves still above threshold at minimum size are retained and
  explicitly flagged, never silently declared to meet the tolerance.
  The processing record explicitly counts valid selected pixels omitted by the
  coverage/count requirements.
- Means use only valid ROI members. Their coordinates, counts, coverage, bounds,
  range, standard deviation and contributing native indices are retained.
  `xy` is the actual member pixel nearest the valid-pixel centroid, so exported
  coordinates remain inside the ROI. `centroid_xy` is saved separately.

`sampling['look_enu']` is the arithmetic mean of the same pixels' projection
coefficients, **without renormalization**. It may be slightly shorter than one.
A tile average assigned to one position is a point-model approximation; a nearest
single look vector would also be an approximation. For scientifically consistent
tile predictions, use:

```python
tile_predictions = insartools.apply_quadtree(pixel_predictions, result['sampling'])
```

`pixel_predictions` starts with the original raster dimensions, optionally with
trailing source/component dimensions. The function uses the exact saved pixel
membership and equal averaging weights. This also works for vector fields.
Within-tile range/std describe variability, including deformation gradients;
neither is a calibrated measurement sigma or a justification for dividing by
the square root of pixel count under spatially correlated noise.

## Files and diagnostics

`write_preprocessing(result, new_directory)` refuses existing directories.
Exclusive directory creation establishes one writer; no shared caches or global
processing state are used. Calculations are serial and deterministic. Failed
writes surface as errors; partially written runs are not reused automatically.

- `processing.json`: units/signs, CRS description, source/array/code hashes,
  polygons, parameters, fit coefficients, masks/counts and limitations.
- `processed.npz`: original standardized LOS, corrected LOS, removed surface,
  coordinates, geometry, masks and sampling arrays. All arrays load with
  `allow_pickle=False`.
- `samples.txt`: x, y, LOS meters, mean E/N/U, with coordinate/sign headers.
  These six columns alone are not an exact representation of area averaging.
- `plot_preprocessing(result, filename)`: six-panel original, actual fitting
  mask, surface, corrected scene, tile means and tile residual diagnostic.
  Black ROI outlines are shown on original/corrected panels. Existing figures
  are refused; call Matplotlib plotting serially.

To reconstruct the averaging definition from the NPZ, use `source_shape`,
`sample_source_indices`, `sample_offsets` and `sample_counts`. Each interval of
indices between consecutive offsets is one tile's native flat-index membership.
`apply_quadtree` accepts those arrays as a sampling dictionary.

The legacy `quadtree2mudpy` helper remains unchanged for existing workflows. It
assumes geographic coordinates and positive-toward LOS. Do not feed it unknown
conventions, projected coordinates or unreviewed point approximations.

## Run the reference examples

From this MudPy worktree, with an existing NumPy/Matplotlib Python environment:

```sh
export PYTHONPATH="$PWD/src/python"
python examples/insar_preprocessing/preprocess.py --output /existing/parent/synthetic
```

The default synthetic reference has a known quadratic background, compact
deformation, missing pixels, separate ROI/exclusion masks and an odd rectangular
native grid. It checks recovery of the known deformation and saves
`validation.json` and `overview.png` alongside the processing outputs.

The optional recovered Norcia reader validates native cell identity before
placing each already supplied look vector at its native cell:

```sh
python examples/insar_preprocessing/preprocess.py \
  --norcia /Users/dmelgarm/Amatrice2016/InSAR/M6.6/T44_Italy_1027_1102 \
  --output /existing/parent/norcia --threshold 0.02 --min-size 4 --max-size 64
```

It explicitly uses the historical mm scale and defaults to **no trend removal**
and **unknown LOS/vector directions**. It does not settle historical conventions.
For a reviewed background correction, add `--trend plane` or `--trend quadratic`
and explicit `--exclude` and/or `--background` polygon files. `--include`,
`--exclude` and `--background` may be repeated and accept Nx2 text or GeoJSON.
Supply sign/direction flags only after independently establishing them.

Run the invariant tests without installing MudPy:

```sh
PYTHONPATH=src/python python -m unittest discover -s src/python/tests -p test_insartools.py -v
```

## Validation recorded November 9, 2026

All 14 invariant tests pass with the existing Python 3.8/NumPy 1.21 environment
and Python 3.12 environment. They cover polygon selection/holes, independent
background fitting, weighted fits, missing pixels, odd raster sizes, exact
averaging of forward predictions, conventions and reproducible exclusive outputs.

The default synthetic example recovers the known deformation after quadratic
removal to a maximum absolute error of `3.61e-16 m` before sampling. It produces
195 tiles; 35 reach minimum size above the range threshold, and four valid ROI
pixels are omitted by the coverage requirement.

The recovered Norcia example, with no trend removal, a `0.02 m` range threshold,
and nominal sizes 4–64 pixels, produces 13,697 tiles from 851,452 pixels. Of those
tiles, 2,820 still exceed threshold at minimum size; 729 of 852,181 valid pixels
are omitted by coverage. These are demonstration settings, not adopted inversion
settings. LOS and vector directions remain explicitly unknown. Scientific source
file hashes match the recovered inventory, and each run records its input and
processor hashes. No source data were modified or dependencies installed.
