'''
InSAR tools go in here
Diego Melgar
University of Oregon
'''

import hashlib
import json
from pathlib import Path

import numpy as np


def quadtree2mudpy(home,project_name,quadtree_file,gflist_file,prefix):
    '''
    Convert quadtree Matlab generated file into a .sta file with station codes 
    and locations. Also generate individaul .neu files for each los displacement
    '''
    
    from numpy import genfromtxt,c_,savetxt
    from matplotlib import pyplot as plt
    from matplotlib import cm
    
    insar=genfromtxt(quadtree_file)
    plt.figure()
    vmin=insar[:,2].min()
    vmax=insar[:,2].max()
    #Write to file
    gflist=open(gflist_file,'w')
    #Write gflist header
    gflist.write('#station	lat	lon	static	disp	vel	tsun	strain	Static file	 displacement file	 velocity file	tsunami file	strain file	static sigmas(n	e	u)	displacement sigmas(n	e	u)	velocity sigmas(n	e	u)	tsunami sigma	strain sigmas(5 components?)\n')																						
    for k in range(len(insar)):
        sta=prefix+str(k).rjust(4,'0')
        #Now make .los file
        los_c=insar[k,2]
        N=insar[k,4]
        E=insar[k,3]
        U=insar[k,5]
        los=c_[los_c,N,E,U]
        savetxt(home+project_name+'/data/statics/'+sta+'.los',los,header='LOS(m),los unit vector (positive towards satellite) n,e,u')
        #Generate gflist file as well
        gflist.write('%s\t%10.6f\t%10.6f\t0\t0\t0\t0\t1\t/foo/bar\t/foo/bar\t/foo/bar\t/foo/bar\t%s\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\n'%(sta,insar[k,0],insar[k,1],home+project_name+'/data/statics/'+sta+'.los'))
    
    #Make plot  to verify
    plt.scatter(insar[:,0],insar[:,1],c=insar[:,2],cmap=cm.magma,s=80,lw=0,vmin=vmin,vmax=vmax)
    plt.legend(['LOS'])
    plt.colorbar()
    gflist.close() 
    plt.show()  
    
    
    
def un_nanify(infile,outfile):        
    from numpy import genfromtxt,savetxt,where,nan
    insar=genfromtxt(infile)
    i=where(infile[:,2]!=nan)[0] 
    insar=insar[i,:]
    savetxt(outfile,insar,fmt='%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f')


# Array-based preprocessing. These functions do not modify their inputs or use
# shared caches. Coordinates and polygon vertices must use the same CRS.
def _float_array(values):
    if np.iscomplexobj(values):
        raise ValueError('InSAR fields must be real-valued')
    return np.asarray(np.ma.filled(np.ma.asarray(values, dtype=float), np.nan))


def _coordinates(x, y, shape):
    x, y = _float_array(x), _float_array(y)
    if x.ndim == y.ndim == 1 and (len(y), len(x)) == shape:
        x, y = np.broadcast_arrays(x[None, :], y[:, None])
    elif x.shape != shape or y.shape != shape:
        raise ValueError('Coordinates must be x/y axes or arrays matching the raster')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Coordinates must be finite')
    return x, y


def _boolean_mask(mask, shape):
    if mask is None:
        return np.ones(shape, dtype=bool)
    mask = np.asarray(mask)
    if mask.shape != shape or mask.dtype != np.bool_:
        raise ValueError('Masks must be boolean arrays matching the raster')
    return mask.copy()


def _segments_intersect(a, b, c, d):
    def cross(p, q, r):
        u, v = q-p, r-p
        return float(u[0]*v[1]-u[1]*v[0])
    ab_c, ab_d = cross(a, b, c), cross(a, b, d)
    cd_a, cd_b = cross(c, d, a), cross(c, d, b)
    if ab_c*ab_d < 0 and cd_a*cd_b < 0:
        return True
    for p, q, r, orientation in ((a, b, c, ab_c), (a, b, d, ab_d),
                                  (c, d, a, cd_a), (c, d, b, cd_b)):
        if orientation == 0 and np.all(r >= np.minimum(p, q)) and np.all(r <= np.maximum(p, q)):
            return True
    return False


def _ring(vertices):
    ring = _float_array(vertices)
    if ring.ndim != 2 or ring.shape[1] != 2 or not np.isfinite(ring).all():
        raise ValueError('Polygon rings must contain finite x/y vertex pairs')
    if len(ring) > 1 and np.array_equal(ring[0], ring[-1]):
        ring = ring[:-1]
    if len(ring) < 3 or len(np.unique(ring, axis=0)) != len(ring):
        raise ValueError('A polygon ring needs at least three distinct vertices')
    shifted = ring-ring[0]
    area = np.sum(shifted[:, 0]*np.roll(shifted[:, 1], -1)
                  - shifted[:, 1]*np.roll(shifted[:, 0], -1))
    if area == 0:
        raise ValueError('Polygon ring has zero signed area')
    for i in range(len(ring)):
        for j in range(i+1, len(ring)):
            if j == i+1 or (i == 0 and j == len(ring)-1):
                continue
            if _segments_intersect(ring[i], ring[(i+1) % len(ring)],
                                   ring[j], ring[(j+1) % len(ring)]):
                raise ValueError('Self-intersecting polygon ring')
    return ring.copy()


def _inside_ring(x, y, ring):
    inside = np.zeros(x.shape, dtype=bool)
    boundary = np.zeros(x.shape, dtype=bool)
    # Explicit edge inclusion, independent of vertex orientation. A tolerance of
    # a few floating-point ulps handles coordinates on a supplied polygon edge.
    scale = max(1., float(np.max(np.abs(ring))), float(np.max(np.ptp(ring, axis=0))))
    tolerance = 16*np.finfo(float).eps*scale
    for a, b in zip(ring, np.roll(ring, -1, axis=0)):
        dx, dy = b-a
        cross = (x-a[0])*dy-(y-a[1])*dx
        boundary |= ((np.abs(cross) <= tolerance*np.hypot(dx, dy))
                     & (x >= min(a[0], b[0])-tolerance)
                     & (x <= max(a[0], b[0])+tolerance)
                     & (y >= min(a[1], b[1])-tolerance)
                     & (y <= max(a[1], b[1])+tolerance))
        if dy != 0:
            crossing = (a[1] > y) != (b[1] > y)
            inside ^= crossing & (x < a[0]+(y-a[1])*dx/dy)
    return inside | boundary, boundary


def _polygons(polygons):
    """Normalize simple rings/unions and GeoJSON polygons with holes."""
    if polygons is None:
        return []
    if isinstance(polygons, dict):
        kind = polygons.get('type')
        if kind == 'Feature':
            return _polygons(polygons['geometry'])
        if kind == 'FeatureCollection':
            return [p for feature in polygons['features'] for p in _polygons(feature)]
        if kind == 'Polygon':
            rings = polygons['coordinates']
            if not rings:
                raise ValueError('Empty Polygon')
            normalized = [_ring(r) for r in rings]
            outer = normalized[0]
            for i, hole in enumerate(normalized[1:]):
                inside, edge = _inside_ring(hole[:, 0], hole[:, 1], outer)
                if not inside.all() or edge.any():
                    raise ValueError('Polygon holes must lie strictly inside the exterior')
                for other in normalized[:i+1]:
                    if any(_segments_intersect(a, b, c, d)
                           for a, b in zip(hole, np.roll(hole, -1, axis=0))
                           for c, d in zip(other, np.roll(other, -1, axis=0))):
                        raise ValueError('Polygon rings intersect')
                    if other is not outer:
                        if (_inside_ring(hole[:1, 0], hole[:1, 1], other)[0].any()
                                or _inside_ring(other[:1, 0], other[:1, 1], hole)[0].any()):
                            raise ValueError('Polygon holes overlap')
            return [normalized]
        if kind == 'MultiPolygon':
            return [p for coordinates in polygons['coordinates']
                    for p in _polygons({'type': 'Polygon', 'coordinates': coordinates})]
        raise ValueError('Expected Polygon, MultiPolygon, Feature, or FeatureCollection')
    if not isinstance(polygons, (list, tuple, np.ndarray)):
        raise ValueError('Expected a vertex array, polygon sequence, or GeoJSON')
    if len(polygons) == 0:
        raise ValueError('Empty polygon collection; use None for no polygon restriction')
    first = np.asarray(polygons[0])
    if first.shape == (2,) and all(np.isscalar(v) for v in first):
        return [[_ring(polygons)]]
    return [p for polygon in polygons for p in _polygons(polygon)]


def polygon_mask(x, y, polygons):
    """Return membership in a polygon union, including polygon boundaries.

    x/y must be matching coordinate arrays (use meshgrid for 1-D axes). Accepts
    an Nx2 vertex array, a sequence of such polygons, or GeoJSON geometries/
    Features/FeatureCollections. GeoJSON holes are supported. None selects all
    points. No coordinate transformation or implicit CRS interpretation occurs.
    """
    x, y = _float_array(x), _float_array(y)
    if x.shape != y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Polygon test coordinates must be matching finite arrays')
    if polygons is None:
        return np.ones(x.shape, dtype=bool)
    normalized = _polygons(polygons)
    if not normalized:
        raise ValueError('Empty polygon collection')
    selected = np.zeros(x.shape, dtype=bool)
    for rings in normalized:
        member, _ = _inside_ring(x, y, rings[0])
        for hole in rings[1:]:
            in_hole, on_edge = _inside_ring(x, y, hole)
            member &= ~(in_hole & ~on_edge)
        selected |= member
    return selected


def fit_surface(displacement, x, y, fit_mask, model='plane', weights=None):
    """Fit an offset, plane, or full quadratic to explicitly selected pixels.

    The fit uses centered/scaled input coordinates and ordinary or weighted
    least squares. Weights are nonnegative least-squares weights (not their
    square roots). Rank-deficient fits raise; the polynomial order is never
    silently reduced. RMS/scatter are diagnostics, not calibrated uncertainties.
    """
    data = _float_array(displacement)
    if data.ndim != 2:
        raise ValueError('Displacement must be a 2-D raster')
    x, y = _coordinates(x, y, data.shape)
    mask = _boolean_mask(fit_mask, data.shape) & np.isfinite(data)
    powers = {'offset': [(0, 0)], 'plane': [(0, 0), (1, 0), (0, 1)],
              'quadratic': [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]}
    if model not in powers:
        raise ValueError('Surface model must be offset, plane, or quadratic')
    if weights is None:
        weight = np.ones(data.shape)
    else:
        weight = _float_array(weights)
        if weight.shape != data.shape or not np.isfinite(weight[mask]).all() or (weight[mask] < 0).any():
            raise ValueError('Fit weights must match the raster and be finite/nonnegative on fit pixels')
        mask &= weight > 0
    if mask.sum() < len(powers[model]):
        raise ValueError('Too few background pixels for the requested surface')
    origin = [float(x[mask].mean()), float(y[mask].mean())]
    scale = [float(np.ptp(x[mask])/2), float(np.ptp(y[mask])/2)]
    # A constant coordinate needs no scaling; rank checking below still rejects
    # a plane/quadratic whose geometry cannot determine all requested terms.
    scale = [s if s > 0 else 1. for s in scale]
    xn, yn = (x-origin[0])/scale[0], (y-origin[1])/scale[1]
    terms = powers[model]
    matrix = np.column_stack([xn[mask]**i*yn[mask]**j for i, j in terms])
    w = np.sqrt(weight[mask])
    coefficients, _, rank, singular = np.linalg.lstsq(matrix*w[:, None], data[mask]*w, rcond=None)
    if rank != len(terms):
        raise ValueError('Background surface is rank deficient; change the fit mask/model')
    surface = sum(c*xn**i*yn**j for c, (i, j) in zip(coefficients, terms))
    if not np.isfinite(surface).all():
        raise ValueError('Fitted surface is non-finite over the supplied scene')
    residual = data[mask]-surface[mask]
    return {'surface': surface, 'fit_mask': mask, 'model': model,
            'powers': [list(p) for p in terms], 'coefficients': coefficients.tolist(),
            'origin': origin, 'scale': scale, 'rank': int(rank),
            'condition_number': float(singular[0]/singular[-1]),
            'fit_pixels': int(mask.sum()), 'weighted': weights is not None,
            'fit_rms': float(np.sqrt(np.mean(residual**2))),
            'weighted_fit_rms': float(np.sqrt(np.average(residual**2, weights=weight[mask])))}


def apply_quadtree(values, sampling):
    """Apply the stored tile means to another native-grid field or predictions.

    Leading dimensions must match source_shape; trailing dimensions are kept.
    This applies exactly the observation's averaging to model predictions or
    geometry. Missing/non-finite values at contributing pixels raise.
    """
    values = _float_array(values)
    shape = tuple(sampling['source_shape'])
    if values.shape[:2] != shape:
        raise ValueError('Values must start with the original raster dimensions')
    flat = values.reshape((int(np.prod(shape)),) + values.shape[2:])
    selected = flat[sampling['source_indices']]
    if not np.isfinite(selected).all():
        raise ValueError('Non-finite values at quadtree member pixels')
    totals = np.add.reduceat(selected, sampling['offsets'][:-1], axis=0)
    divisor = sampling['counts'].reshape((-1,) + (1,)*(values.ndim-2))
    return totals/divisor


def quadtree(displacement, x, y, look_enu, *, selection_mask=None,
             valid_mask=None, threshold=0.03, min_size=4, max_size=128,
             criterion='range', min_valid_fraction=0.5, min_pixels=1):
    """Average native pixels in a deterministic quadtree; never interpolate.

    Tile sizes are powers of two in native pixels. The conceptual square tree
    is clipped at rectangular/odd image edges, with no synthetic pixels added.
    Polygon-selection boundaries split to min_size. Coverage is valid pixels /
    selected pixels, so pixels outside the ROI cannot bias averages or coverage.
    Leaves still exceeding threshold at min_size are returned and flagged.

    look_enu contains per-pixel ENU projection coefficients. Their tile means
    are not renormalized. xy is a real member pixel nearest the valid-pixel
    centroid; centroid_xy is also retained. A single-point forward calculation
    is an approximation: use apply_quadtree on pixel predictions for exact means.
    """
    data = _float_array(displacement)
    if data.ndim != 2 or min(data.shape) < 1:
        raise ValueError('Displacement must be a nonempty 2-D raster')
    x, y = _coordinates(x, y, data.shape)
    look = _float_array(look_enu)
    if look.shape != data.shape + (3,):
        raise ValueError('Look vectors must have shape (rows, columns, 3) in ENU order')
    for size in (min_size, max_size, min_pixels):
        if isinstance(size, (bool, np.bool_)) or not isinstance(size, (int, np.integer)) or size < 1:
            raise ValueError('Tile sizes and min_pixels must be positive integers')
    if (min_size & (min_size-1)) or (max_size & (max_size-1)) or min_size > max_size:
        raise ValueError('Tile sizes must be powers of two with min_size <= max_size')
    if not np.isfinite(threshold) or threshold < 0 or criterion not in ('range', 'std'):
        raise ValueError('Use a finite nonnegative threshold and range/std criterion')
    if not np.isfinite(min_valid_fraction) or not 0 <= min_valid_fraction <= 1:
        raise ValueError('min_valid_fraction must lie between zero and one')
    selected = _boolean_mask(selection_mask, data.shape)
    valid = (_boolean_mask(valid_mask, data.shape) & selected
             & np.isfinite(data) & np.isfinite(look).all(axis=2))
    height, width = data.shape
    flat_x, flat_y = x.ravel(), y.ravel()
    root_size = 1 << (max(height, width)-1).bit_length()
    stack = [(0, 0, root_size)]
    members, bounds, counts, fractions, spreads, stds, exceeded = [], [], [], [], [], [], []
    centroids, representatives = [], []
    while stack:
        r0, c0, size = stack.pop()
        r1, c1 = min(r0+size, height), min(c0+size, width)
        if r0 >= height or c0 >= width:
            continue
        eligible = int(selected[r0:r1, c0:c1].sum())
        mask = valid[r0:r1, c0:c1]
        count = int(mask.sum())
        if not eligible or count < min_pixels:
            continue
        values = data[r0:r1, c0:c1][mask]
        fraction = count/eligible
        spread, std = float(np.ptp(values)), float(np.std(values))
        metric = spread if criterion == 'range' else std
        boundary = eligible != (r1-r0)*(c1-c0)
        split = size > max_size or boundary or fraction < min_valid_fraction or metric > threshold
        if split and size > min_size:
            half = size//2
            stack.extend([(r0+half, c0+half, half), (r0+half, c0, half),
                          (r0, c0+half, half), (r0, c0, half)])
            continue
        if fraction < min_valid_fraction:
            continue
        rr, cc = np.nonzero(mask)
        indices = (rr+r0)*width+(cc+c0)
        xx, yy = flat_x[indices], flat_y[indices]
        center = [float(xx.mean()), float(yy.mean())]
        nearest = int(np.argmin((xx-center[0])**2+(yy-center[1])**2))
        representatives.append(int(indices[nearest]))
        centroids.append(center)
        members.append(indices)
        bounds.append([r0, r1, c0, c1])
        counts.append(count)
        fractions.append(fraction)
        spreads.append(spread)
        stds.append(std)
        exceeded.append(metric > threshold)
    if not members:
        raise ValueError('No quadtree leaves satisfy the masks and coverage requirements')
    result = {'source_shape': list(data.shape), 'source_indices': np.concatenate(members),
              'offsets': np.r_[0, np.cumsum(counts)], 'counts': np.asarray(counts),
              'bounds': np.asarray(bounds), 'valid_fraction': np.asarray(fractions),
              'range': np.asarray(spreads), 'std': np.asarray(stds),
              'threshold_exceeded': np.asarray(exceeded, dtype=bool),
              'centroid_xy': np.asarray(centroids),
              'representative_indices': np.asarray(representatives),
              'selected_pixels': int(selected.sum()), 'valid_selected_pixels': int(valid.sum()),
              'parameters': {'threshold': float(threshold), 'criterion': criterion,
                             'min_size': int(min_size), 'max_size': int(max_size),
                             'min_valid_fraction': float(min_valid_fraction), 'min_pixels': int(min_pixels)}}
    result['xy'] = np.column_stack([flat_x[representatives], flat_y[representatives]])
    result['los'] = apply_quadtree(data, result)
    result['look_enu'] = apply_quadtree(look, result)
    return result


def _array_identity(array):
    array = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(str((array.dtype.str, array.shape)).encode())
    digest.update(array.tobytes())
    return {'shape': list(array.shape), 'dtype': array.dtype.str, 'sha256': digest.hexdigest()}


def _file_identity(filename):
    path = Path(filename).resolve(strict=True)
    before = path.stat()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError('Source changed while hashing: '+str(path))
    return {'path': str(path), 'bytes': before.st_size, 'sha256': digest.hexdigest()}


def preprocess_insar(displacement, x, y, look_enu, *, units, los_positive,
                     look_direction, coordinate_system, include_polygons=None,
                     exclude_polygons=None, background_polygons=None,
                     quality_mask=None, nodata=None, trend=None, fit_weights=None,
                     threshold=0.03, min_size=4, max_size=128, criterion='range',
                     min_valid_fraction=0.5, min_pixels=1,
                     look_norm_tolerance=1e-3, source_files=(), metadata=None):
    """Standardize, optionally detrend, and downsample an InSAR raster.

    units is m/cm/mm. LOS sign is toward/away/unknown; look_direction is
    ground_to_satellite/satellite_to_ground/unknown. Known conventions are
    converted to positive-toward and ground-to-satellite. Unknown conventions
    remain explicitly unknown: no polarity is inferred from the displacement.

    include_polygons selects observations BEFORE tile statistics. The trend is
    fitted independently on valid pixels inside background_polygons (full scene
    if None) and OUTSIDE exclude_polygons. It is then subtracted over the scene.
    Fitting occurs before sampling. trend=None makes no reference adjustment.

    All results are new arrays. Original units/sign, masks, hashes, polynomial
    coefficients and exact tile membership are retained for reproducibility.
    No reprojection, interpolation, gap filling or smoothing is performed.
    """
    scales = {'m': 1., 'cm': 0.01, 'mm': 0.001}
    if units not in scales or los_positive not in ('toward', 'away', 'unknown'):
        raise ValueError('Explicit m/cm/mm units and toward/away/unknown LOS sign are required')
    if look_direction not in ('ground_to_satellite', 'satellite_to_ground', 'unknown'):
        raise ValueError('Explicit look-vector direction is required')
    if not isinstance(coordinate_system, str) or not coordinate_system.strip():
        raise ValueError('An explicit coordinate-system description is required')
    if not np.isfinite(look_norm_tolerance) or look_norm_tolerance < 0:
        raise ValueError('look_norm_tolerance must be finite and nonnegative')
    raw = _float_array(displacement).copy()
    if raw.ndim != 2 or min(raw.shape) < 1:
        raise ValueError('Displacement must be a nonempty 2-D raster')
    x, y = (a.copy() for a in _coordinates(x, y, raw.shape))
    look = _float_array(look_enu).copy()
    if look.shape != raw.shape+(3,):
        raise ValueError('Look vectors must have shape (rows, columns, 3) in ENU order')
    quality = _boolean_mask(quality_mask, raw.shape)
    missing = ~np.isfinite(raw)
    if nodata is not None:
        if not np.isscalar(nodata) or not np.isfinite(nodata):
            raise ValueError('nodata must be a finite numeric sentinel; NaNs are already masked')
        missing |= raw == nodata
    finite_look = np.isfinite(look).all(axis=2)
    valid_los = ~missing & quality
    valid = valid_los & finite_look
    keep = polygon_mask(x, y, include_polygons)
    if not (valid & keep).any():
        raise ValueError('No valid selected displacement/geometry pixels')
    norms = np.linalg.norm(look[valid & keep], axis=1)
    if np.max(np.abs(norms-1)) > look_norm_tolerance:
        raise ValueError('Input look vectors exceed the declared unit-length tolerance')
    raw_identity = _array_identity(raw)
    look_identity = _array_identity(look)
    original = raw*scales[units]*(-1 if los_positive == 'away' else 1)
    original[missing] = np.nan
    if look_direction == 'satellite_to_ground':
        look *= -1
    background = polygon_mask(x, y, background_polygons)
    if exclude_polygons is not None:
        background &= ~polygon_mask(x, y, exclude_polygons)
    background &= valid_los
    if trend is None:
        if fit_weights is not None:
            raise ValueError('fit_weights supplied while trend removal is disabled')
        surface, used = np.zeros(raw.shape), np.zeros(raw.shape, dtype=bool)
        fit = None
    else:
        fitted = fit_surface(original, x, y, background, model=trend, weights=fit_weights)
        surface, used = fitted['surface'], fitted['fit_mask']
        fit = {key: value for key, value in fitted.items() if key not in ('surface', 'fit_mask')}
    corrected = np.where(valid_los, original-surface, np.nan)
    sampling = quadtree(corrected, x, y, look, selection_mask=keep, valid_mask=valid,
                        threshold=threshold, min_size=min_size, max_size=max_size,
                        criterion=criterion, min_valid_fraction=min_valid_fraction,
                        min_pixels=min_pixels)
    def polygon_record(value):
        return None if value is None else {'type': 'MultiPolygon', 'coordinates': [
            [ring.tolist() for ring in p] for p in _polygons(value)]}
    provenance = {
        'coordinate_system': coordinate_system, 'input_units': units, 'output_units': 'm',
        'input_los_positive': los_positive,
        'output_los_positive': 'unknown' if los_positive == 'unknown' else 'toward',
        'input_look_direction': look_direction,
        'output_look_direction': 'unknown' if look_direction == 'unknown' else 'ground_to_satellite',
        'input_arrays': {'displacement': raw_identity, 'look_enu': look_identity,
                         'x': _array_identity(x), 'y': _array_identity(y),
                         'quality_mask': _array_identity(quality)},
        'source_files': [_file_identity(p) for p in
                         ([source_files] if isinstance(source_files, (str, Path)) else source_files)],
        'processor': _file_identity(__file__),
        'metadata': {} if metadata is None else json.loads(json.dumps(metadata, allow_nan=False)),
        'nodata': None if nodata is None else float(nodata),
        'look_norm_tolerance': float(look_norm_tolerance),
        'include_polygons': polygon_record(include_polygons),
        'exclude_polygons': polygon_record(exclude_polygons),
        'background_polygons': polygon_record(background_polygons),
        'fit_weights': None if fit_weights is None else _array_identity(_float_array(fit_weights)),
        'trend': fit, 'sampling': sampling['parameters'],
        'counts': {'source_pixels': raw.size, 'valid_pixels': int(valid.sum()),
                   'valid_los_pixels': int(valid_los.sum()),
                   'missing_displacement_pixels': int(missing.sum()),
                   'missing_geometry_pixels': int((~finite_look).sum()),
                   'quality_excluded_pixels': int((~quality).sum()),
                   'selected_pixels': int(keep.sum()), 'fit_pixels': int(used.sum()),
                   'valid_selected_pixels': sampling['valid_selected_pixels'],
                   'sampled_pixels': int(sampling['counts'].sum()),
                   'dropped_valid_selected_pixels': int(sampling['valid_selected_pixels']-sampling['counts'].sum()),
                   'leaves': len(sampling['counts']),
                   'leaves_exceeding_threshold_at_min_size': int(sampling['threshold_exceeded'].sum())},
        'limitations': ['Tile scatter is not measurement uncertainty.',
                       'A point-location forward model approximates a tile mean; membership is retained.',
                       'Surface-fit uncertainty and spatial noise covariance are not estimated.']}
    # Validate metadata now, rather than failing halfway through output writing.
    json.dumps(provenance, allow_nan=False)
    return {'x': x, 'y': y, 'original_m': original, 'corrected_m': corrected,
            'trend_m': surface, 'look_enu': look, 'valid_mask': valid,
            'los_valid_mask': valid_los,
            'keep_mask': keep, 'fit_mask': used, 'sampling': sampling,
            'fit_weights': None if fit_weights is None else _float_array(fit_weights).copy(),
            'provenance': provenance}


def write_preprocessing(result, output_directory):
    """Write an inspectable processing record and arrays to a NEW directory.

    Exclusive directory creation gives one writer ownership; reruns and other
    writers cannot overwrite an existing run. No inversion is launched. The
    six-column samples.txt is an approximation for point-based modeling, with
    unnormalized mean ENU coefficients; exact membership is in processed.npz.
    """
    text = json.dumps(result['provenance'], indent=2, allow_nan=False)+'\n'
    sampling = result['sampling']
    arrays = {key: result[key] for key in ('x', 'y', 'original_m', 'corrected_m',
               'trend_m', 'look_enu', 'valid_mask', 'los_valid_mask', 'keep_mask', 'fit_mask')}
    arrays.update({'sample_'+key: value for key, value in sampling.items()
                   if isinstance(value, np.ndarray)})
    arrays['source_shape'] = np.asarray(sampling['source_shape'])
    if result['fit_weights'] is not None:
        arrays['fit_weights'] = result['fit_weights']
    directory = Path(output_directory)
    directory.mkdir()  # Parents must already exist; existing runs are refused.
    (directory/'processing.json').write_text(text)
    np.savez_compressed(directory/'processed.npz', **arrays)
    header = ('x y LOS_m mean_E mean_N mean_U\n'
              'coordinate_system='+result['provenance']['coordinate_system']+'\n'
              'LOS_positive='+result['provenance']['output_los_positive']+'; look_direction='
              +result['provenance']['output_look_direction']+'\n'
              'Point representation approximates tile means; coefficients are not renormalized. '
              'Tile scatter is not sigma.')
    np.savetxt(directory/'samples.txt', np.column_stack([
        sampling['xy'], sampling['los'], sampling['look_enu']]), fmt='%.17g', header=header)
    return directory


def plot_preprocessing(result, filename):
    """Save six diagnostic panels; call serially (Matplotlib has global state)."""
    from matplotlib import pyplot as plt
    shape = result['original_m'].shape
    if min(shape) < 2:
        raise ValueError('Diagnostic maps require at least two rows and columns')
    sampling = result['sampling']
    tiled = np.full(shape, np.nan)
    tiled.ravel()[sampling['source_indices']] = np.repeat(sampling['los'], sampling['counts'])
    fields = [result['original_m'], result['fit_mask'].astype(float), result['trend_m'],
              result['corrected_m'], tiled, result['corrected_m']-tiled]
    limit = max(float(np.nanmax(np.abs(f))) for i, f in enumerate(fields) if i != 1)
    if limit == 0:
        limit = 1.  # A constant zero field still needs a defined display scale.
    titles = ['Original LOS (standardized units/sign)', 'Actual background-fit pixels',
              'Removed polynomial surface', 'Corrected LOS', 'Quadtree means',
              'Corrected LOS minus tile mean']
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    for i, (ax, field, title) in enumerate(zip(axes.flat, fields, titles)):
        palette = plt.get_cmap('gray' if i == 1 else 'RdBu_r').copy()
        palette.set_bad('#dedede')
        image = ax.pcolormesh(result['x'], result['y'], np.ma.masked_invalid(field),
                              shading='nearest', cmap=palette,
                              vmin=0 if i == 1 else -limit, vmax=1 if i == 1 else limit,
                              rasterized=True)
        if i in (0, 3) and result['keep_mask'].any() and not result['keep_mask'].all():
            ax.contour(result['x'], result['y'], result['keep_mask'], levels=[0.5],
                       colors='black', linewidths=0.6)
        if i == 4:
            for r0, r1, c0, c1 in sampling['bounds']:
                rr, cc = [r0, r0, r1-1, r1-1, r0], [c0, c1-1, c1-1, c0, c0]
                ax.plot(result['x'][rr, cc], result['y'][rr, cc], color='black', lw=0.25)
        ax.set(title=title, xlabel='x (input coordinates)', ylabel='y (input coordinates)')
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04).set_label('Selected' if i == 1 else 'm')
    counts = result['provenance']['counts']
    fig.suptitle('InSAR preprocessing | '+result['provenance']['coordinate_system']+'\n'
                 +str(counts['leaves'])+' tiles; '+str(counts['sampled_pixels'])
                 +' contributing pixels; LOS positive '+result['provenance']['output_los_positive'])
    path = Path(filename)
    format_name = path.suffix.lstrip('.')
    if format_name not in fig.canvas.get_supported_filetypes():
        plt.close(fig)
        raise ValueError('Use a supported figure filename suffix')
    try:
        with path.open('xb') as stream:
            fig.savefig(stream, format=format_name, dpi=160)
    finally:
        plt.close(fig)
