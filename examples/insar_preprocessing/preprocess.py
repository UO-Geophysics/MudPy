#!/usr/bin/env python3
"""Run the synthetic reference or inspect the recovered Norcia pair.

Set PYTHONPATH to MudPy's src/python; no installation is needed. The output
directory must be new. Rasterio is used only by the optional Norcia reader.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mudpy import insartools


def read_polygons(paths):
    if paths is None:
        return None
    return [json.loads(p.read_text()) if p.suffix.lower() in ('.json', '.geojson')
            else np.loadtxt(p) for p in paths]


def synthetic():
    x, y = np.meshgrid(np.linspace(-20000, 20000, 97), np.linspace(-15000, 15000, 65))
    background = 0.08+2e-6*x-1e-6*y+1e-10*x*y+8e-11*x*x-5e-11*y*y
    radius = (x/7000)**2+(y/5000)**2
    signal = 0.25*np.maximum(0, 1-radius)**2
    data = background+signal
    data[5:12, 60:65] = np.nan
    east, north = -0.5+2e-6*x, -0.11+2e-7*y
    look = np.stack([east, north, np.sqrt(1-east**2-north**2)], axis=2)
    angle = np.linspace(0, 2*np.pi, 65)
    exclude = np.column_stack([8500*np.cos(angle), 6500*np.sin(angle)])
    include = [[-18000, -14000], [18000, -14000], [18000, 14000], [-18000, 14000]]
    return data, x, y, look, signal, include, exclude


def read_norcia(folder):
    import rasterio
    grid = folder/'los_ll.grd'
    table = folder/'T44_Italy_1027_1102.lltnde'
    with rasterio.Env(GDAL_PAM_ENABLED='NO'), rasterio.open(grid) as src:
        if src.transform.b != 0 or src.transform.d != 0:
            raise ValueError('This example expects the recovered unrotated Norcia grid')
        data = src.read(1, masked=True).filled(np.nan)
        x = src.transform.c+(np.arange(src.width)+0.5)*src.transform.a
        y = src.transform.f+(np.arange(src.height)+0.5)*src.transform.e
        points = np.loadtxt(table)
        if points.ndim != 2 or points.shape[1] != 8 or not np.isfinite(points).all():
            raise ValueError('Expected a finite eight-column recovered observation table')
        rr, cc = rasterio.transform.rowcol(src.transform, points[:, 0], points[:, 1])
        rr, cc = np.asarray(rr), np.asarray(cc)
        if ((rr < 0) | (rr >= src.height) | (cc < 0) | (cc >= src.width)).any():
            raise ValueError('Table coordinates fall outside the original raster')
        if len(np.unique(rr*src.width+cc)) != len(points) or len(points) != np.isfinite(data).sum():
            raise ValueError('Observation table does not cover native valid cells one-to-one')
        if not np.allclose(data[rr, cc], points[:, 6], rtol=0, atol=1e-6):
            raise ValueError('Displacements disagree between the grid and table')
    look = np.full(data.shape+(3,), np.nan)
    look[rr, cc] = points[:, 3:6]
    return data, x, y, look, [grid, table]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--norcia', type=Path, help='Recovered T44_Italy_1027_1102 folder')
    parser.add_argument('--include', type=Path, action='append', help='Observation polygon: vertex text or GeoJSON; repeat for union')
    parser.add_argument('--exclude', type=Path, action='append', help='Exclude deformation from surface fitting; repeat for union')
    parser.add_argument('--background', type=Path, action='append', help='Optional explicit fitting domain; repeat for union')
    parser.add_argument('--trend', choices=('offset', 'plane', 'quadratic'))
    parser.add_argument('--threshold', type=float, default=0.02, help='Tile LOS range threshold in meters')
    parser.add_argument('--min-size', type=int, default=4, help='Minimum nominal tile size in native pixels')
    parser.add_argument('--max-size', type=int, default=64, help='Maximum nominal tile size in native pixels')
    parser.add_argument('--los-positive', choices=('toward', 'away', 'unknown'), default='unknown')
    parser.add_argument('--look-direction', choices=('ground_to_satellite', 'satellite_to_ground', 'unknown'), default='unknown')
    args = parser.parse_args()
    include, exclude = read_polygons(args.include), read_polygons(args.exclude)
    background = read_polygons(args.background)
    if args.norcia:
        if args.trend is not None and exclude is None and background is None:
            parser.error('Norcia trend fitting requires an explicit --exclude or --background polygon')
        data, x, y, look, files = read_norcia(args.norcia)
        units, positive, direction = 'mm', args.los_positive, args.look_direction
        coordinate_system = 'Geographic longitude/latitude degrees; datum not encoded in source GMT grid'
        trend, truth = args.trend, None
        metadata = {'dataset': 'Recovered Norcia T44, nominal 2016-10-27 to 2016-11-02',
                    'unit_evidence': 'Historical project_GPS_M6.6.py divides table column 7 by 1000',
                    'geometry_evidence': 'Historical scripts identify table columns 4-6 as E/N/U',
                    'reference_status': 'Historical reference/correction history remains under audit'}
    else:
        data, x, y, look, truth, default_include, default_exclude = synthetic()
        include = default_include if include is None else include
        exclude = default_exclude if exclude is None else exclude
        units, positive, direction = 'm', 'toward', 'ground_to_satellite'
        coordinate_system, trend = 'Synthetic local east/north meters', args.trend or 'quadratic'
        files, metadata = [Path(__file__)], {'dataset': 'Known quadratic background plus compact deformation; noiseless'}
    files += (args.include or [])+(args.exclude or [])+(args.background or [])
    result = insartools.preprocess_insar(data, x, y, look, units=units,
                los_positive=positive, look_direction=direction, coordinate_system=coordinate_system,
                include_polygons=include, exclude_polygons=exclude, background_polygons=background,
                trend=trend, threshold=args.threshold, min_size=args.min_size, max_size=args.max_size,
                source_files=files, metadata=metadata)
    directory = insartools.write_preprocessing(result, args.output)
    insartools.plot_preprocessing(result, directory/'overview.png')
    validation = {'counts': result['provenance']['counts'],
                  'LOS_positive': result['provenance']['output_los_positive'],
                  'look_direction': result['provenance']['output_look_direction'],
                  'trend_model': trend, 'max_abs_tile_residual_m': float(np.max(np.abs(
                      result['corrected_m'].ravel()[result['sampling']['source_indices']]
                      -np.repeat(result['sampling']['los'], result['sampling']['counts']))))}
    if truth is not None:
        finite = np.isfinite(result['corrected_m'])
        error = float(np.max(np.abs(result['corrected_m'][finite]-truth[finite])))
        validation['maximum_known_signal_error_m'] = error
        if trend == 'quadratic' and error > 1e-12:
            raise AssertionError('Synthetic background removal did not recover known deformation')
    (directory/'validation.json').write_text(json.dumps(validation, indent=2, allow_nan=False)+'\n')
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    main()
