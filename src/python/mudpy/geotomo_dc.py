"""Pure double-couple inversion on the existing six-tensor matrix."""
import numpy as np

from mudpy.analysis import decompose_geotomo_tensors, double_couple_tensors
from mudpy.geotomo_io import REFERENCE_MOMENT_NM


_MOMENT_WEIGHTS = np.sqrt([.5, 1., 1., .5, 1., .5])


class DCConvergenceError(RuntimeError):
    def __init__(self, diagnostics):
        super().__init__('No DC solver start converged; increase budget or revise initialization')
        self.diagnostics = diagnostics


def dc_reference(tensor):
    """Normalize the standard DC part of an initial NED six-component tensor."""
    d = decompose_geotomo_tensors(tensor)
    if len(d['moment_nm']['full']) != 1:
        raise ValueError('Initialization requires exactly one reference tensor')
    moment = d['moment_nm']['dc'][0]
    if moment <= 64*np.finfo(float).eps*d['moment_nm']['full'][0]:
        raise ValueError('Initial tensor has no defined DC mechanism')
    return d['tensors_nm']['dc'][0]/moment


def _horizontal_neighbors(sources, radius_km):
    """Find unique same-depth pairs in source row order using WGS84 distance."""
    from pyproj import Geod, Transformer
    from scipy.spatial import cKDTree

    lon, lat, depth = sources[:, 1:].T
    radius = 1000*radius_km
    geod = Geod(ellps='WGS84')
    # Surface chords cheaply bound the geodesic distances used for selection.
    transform = Transformer.from_crs(4979, 4978, always_xy=True)
    xyz = np.column_stack(transform.transform(lon, lat, np.zeros_like(lon)))
    edges = []
    for z in np.unique(depth):
        rows = np.flatnonzero(depth == z)
        pairs = rows[cKDTree(xyz[rows]).query_pairs(radius, output_type='ndarray')]
        if len(pairs):
            i, j = pairs.T
            _, _, distance = geod.inv(lon[i], lat[i], lon[j], lat[j])
            if np.any(distance == 0):
                raise ValueError('Horizontal smoothing requires distinct positions within each depth layer')
            edges.extend(pairs[distance <= radius].tolist())
    if not edges:
        raise ValueError('Horizontal smoothing found no same-depth neighbors within horizontal_neighbor_radius_km')
    return np.array(sorted(edges), dtype=int)


def _horizontal_terms(amplitudes, derivative, pairs, *, jacobian=False):
    """Residual or endpoint Jacobian blocks for harmonic-moment weighted DC differences."""
    a, b = amplitudes[pairs].T
    # Normalize by the larger amplitude before taking ratios. A zero-zero
    # pair has zero squared cost/gradient; TRF evaluates at positive moments.
    scale = np.maximum(a, b)
    a_fraction = np.divide(a, scale, out=np.zeros_like(a), where=scale > 0)
    b_fraction = np.divide(b, scale, out=np.zeros_like(b), where=scale > 0)
    denominator = a_fraction+b_fraction
    a_fraction = np.divide(a_fraction, denominator, out=a_fraction, where=denominator > 0)
    b_fraction = np.divide(b_fraction, denominator, out=b_fraction, where=denominator > 0)
    i, j = pairs.T
    # d(a U)/da is U, including when a is zero.
    difference = (derivative[i, :, 0]-derivative[j, :, 0])*_MOMENT_WEIGHTS
    if not jacobian:
        weight = 2*(np.minimum(a, b)*np.maximum(a_fraction, b_fraction))
        return weight[:, None]*difference
    left, right = derivative[i], derivative[j]
    left[:, :, 0] = (2*b_fraction**2)[:, None]*difference
    right[:, :, 0] = (2*a_fraction**2)[:, None]*difference
    # Rotation derivatives already include the endpoint's own amplitude.
    left[:, :, 1:] *= (2*b_fraction)[:, None, None]*_MOMENT_WEIGHTS[None, :, None]
    right[:, :, 1:] *= (-2*a_fraction)[:, None, None]*_MOMENT_WEIGHTS[None, :, None]
    return left, right


def _jacobian_operator(matrix, derivative, damping, pairs=None, edge_blocks=None):
    """One immutable derivative per operator; no shared state between starts."""
    from scipy.sparse.linalg import LinearOperator
    ndata, nsource = len(matrix), len(derivative)
    end = ndata+nsource
    if pairs is not None:
        i, j = pairs.T
        left, right = edge_blocks

    def matvec(vector):
        h = np.asarray(vector).reshape(nsource, 4)
        data = matrix @ np.einsum('nck,nk->nc', derivative, h).ravel()
        if pairs is not None:
            edge = np.einsum('eck,ek->ec', left, h[i])+np.einsum('eck,ek->ec', right, h[j])
            return np.r_[data, damping*h[:, 0], edge.ravel()]
        return np.r_[data, damping*h[:, 0]]

    def rmatvec(vector):
        vector = np.asarray(vector).ravel()
        result = np.einsum('nck,nc->nk', derivative, (matrix.T @ vector[:ndata]).reshape(nsource, 6))
        result[:, 0] += damping*vector[ndata:end]
        if pairs is not None:
            edge = vector[end:].reshape(-1, 6)
            np.add.at(result, i, np.einsum('eck,ec->ek', left, edge))
            np.add.at(result, j, np.einsum('eck,ec->ek', right, edge))
        return result.ravel()

    nrows = end if pairs is None else end+6*len(pairs)
    return LinearOperator((nrows, 4*nsource), matvec=matvec, rmatvec=rmatvec, dtype=float)


def _initial_amplitudes(matrix, data, damping):
    """Solve the convex damped nonnegative amplitude subproblem directly."""
    from scipy.optimize import nnls
    count = matrix.shape[1]
    augmented = np.vstack((matrix, damping*np.eye(count)))
    target = np.r_[data, np.zeros(count)]
    return nnls(augmented, target, maxiter=10*count)[0]


def _l1_objective(parameters, matrix, data, bases, weight, smoothing, pairs):
    """Half the data/L1/smoothing objective and its analytic parameter gradient."""
    coefficients, derivative = double_couple_tensors(parameters, bases)
    residual = matrix @ coefficients.ravel()-data
    blocks = None
    rows = np.r_[residual, np.zeros(len(parameters))]
    if pairs is not None:
        edge = smoothing*_horizontal_terms(parameters[:, 0], derivative, pairs)
        blocks = tuple(smoothing*block for block in
                       _horizontal_terms(parameters[:, 0], derivative, pairs, jacobian=True))
        rows = np.r_[rows, edge.ravel()]
    operator = _jacobian_operator(matrix, derivative, 0., pairs, blocks)
    gradient = operator.rmatvec(rows).reshape(parameters.shape)
    gradient[:, 0] += .5*weight
    return .5*(rows @ rows+weight*parameters[:, 0].sum()), gradient.ravel()


def solve(matrix, data, damping, initial_tensor, *, starts=3, seed=0, perturbation_deg=30.,
          max_nfev=500, ftol=1e-9, xtol=1e-9, gtol=1e-9,
          horizontal_smoothing=0., horizontal_pairs=None, moment_penalty='l2', verbose=False):
    """Fit DC sources with scalar-moment damping and optional mechanism smoothing.

    G is in metres per 10**16.6 Nm coefficient; data/G may already be weighted.
    Amplitudes a are scalar moments divided by that reference, rotations radians.
    horizontal_pairs gives unique undirected pairs of zero-based source rows.
    The smoothing weight multiplies their harmonic-moment weighted differences
    of unit DC tensors, measured in the scalar-moment norm.
    moment_penalty='l2' uses damping**2 * sum(a**2) and TRF/LSMR.
    moment_penalty='l1' uses damping * sum(a) and L-BFGS-B, with analytic
    gradients in scaled amplitude coordinates. ftol/gtol apply to both methods;
    xtol applies to TRF. max_nfev also caps L-BFGS-B iterations.
    TRF/LSMR uses analytic matrix-free Jacobians. Starts run sequentially, each
    from the same catalog DC plus seeded independent source rotations, followed
    by damped NNLS amplitude initialization. No global-optimality claim.
    An unconverged start is recorded but cannot be selected. If all fail, return
    no model: RuntimeError carries the per-start diagnostics for archival.
    """
    from scipy.optimize import least_squares, minimize
    from scipy.spatial.transform import Rotation
    from time import perf_counter, process_time
    solve_wall, solve_cpu = perf_counter(), process_time()
    G, d = np.asarray(matrix, dtype=float), np.asarray(data, dtype=float).reshape(-1)
    if (G.ndim != 2 or not G.shape[0] or not G.shape[1] or G.shape[1] % 6 or
            len(d) != len(G) or not np.isfinite(G).all() or not np.isfinite(d).all() or
            not np.isfinite(damping) or damping < 0):
        raise ValueError('Expected finite matching data/six-component matrix and nonnegative damping')
    if (type(starts) is not int or starts < 1 or type(seed) is not int or seed < 0 or
            type(max_nfev) is not int or max_nfev < 1 or
            not np.isfinite(perturbation_deg) or not 0 <= perturbation_deg <= 180 or
            (starts > 1 and perturbation_deg == 0) or
            any(not np.isfinite(v) or v <= np.finfo(float).eps for v in (ftol, xtol, gtol))):
        raise ValueError('Invalid DC start count, seed, perturbation, evaluation budget or tolerances')
    base = dc_reference(initial_tensor)
    nsource = G.shape[1]//6
    if moment_penalty not in ('l1', 'l2'):
        raise ValueError('moment_penalty must be l1 or l2')
    if not np.isfinite(horizontal_smoothing) or horizontal_smoothing < 0:
        raise ValueError('horizontal_smoothing must be finite and nonnegative')
    pairs = None
    if horizontal_smoothing:
        pairs = np.asarray(horizontal_pairs)
        if (pairs.ndim != 2 or pairs.shape[1] != 2 or not len(pairs) or
                not np.issubdtype(pairs.dtype, np.integer) or
                (pairs < 0).any() or (pairs >= nsource).any() or
                (pairs[:, 0] == pairs[:, 1]).any()):
            raise ValueError('Horizontal smoothing requires nonempty pairs of distinct source rows')
        pairs = np.sort(pairs, axis=1)
        if len(np.unique(pairs, axis=0)) != len(pairs):
            raise ValueError('Horizontal smoothing pairs must be unique')
    if d.any() and not G.any():
        raise ValueError('Nonzero data have no responsive matrix columns')
    # Common residual scaling preserves the objective minimizer and makes
    # convergence tolerances meaningful for both metre-valued and unit tests.
    residual_scale = np.linalg.norm(d) or 1.
    scaled_G, scaled_d, scaled_damping = G/residual_scale, d/residual_scale, damping/residual_scale
    scaled_smoothing = horizontal_smoothing/residual_scale
    scaled_l1 = damping/residual_scale**2
    rng = np.random.default_rng(seed)
    diagnostics, candidates = [], []
    for start in range(starts):
        start_wall, start_cpu = perf_counter(), process_time()
        if verbose:
            print('DC %s start %d/%d: initializing' % (moment_penalty.upper(), start+1, starts), flush=True)
        bases = np.broadcast_to(base, (nsource, 3, 3)).copy()
        if start:
            axes = rng.normal(size=(nsource, 3))
            axes /= np.linalg.norm(axes, axis=1)[:, None]
            angle = np.deg2rad(perturbation_deg)*rng.uniform(size=nsource)
            rotation = Rotation.from_rotvec(axes*angle[:, None]).as_matrix()
            bases = rotation @ bases @ np.swapaxes(rotation, 1, 2)
        p = np.zeros((nsource, 4))
        p[:, 0] = 1
        unit, _ = double_couple_tensors(p, bases)
        projected = np.einsum('rnc,nc->rn', scaled_G.reshape(len(G), nsource, 6), unit)
        amplitude = _initial_amplitudes(projected, scaled_d,
                                        scaled_damping if moment_penalty == 'l2' else 0.)
        initialization_wall = perf_counter()-start_wall
        optimization_wall = perf_counter()
        if not d.any():
            # Known global minimum, including damping=0: no orientation is defined.
            p[:, 0] = 0
            coefficients, derivative = double_couple_tensors(p, bases)
            item = dict(start=start, success=True, status=1, message='Exact zero-data minimum',
                        nfev=0, njev=0, cost=0., scaled_optimality=0., residual_l2=0., moment_l2=0.)
        else:
            # A positive initialization floor frees rotations at zero NNLS nodes.
            # This is not a final lower bound or a prior on recovered moments.
            scale = max(np.linalg.norm(amplitude)/np.sqrt(nsource),
                        np.linalg.norm(scaled_d)/max(np.linalg.norm(projected), np.finfo(float).tiny))
            p[:, 0] = np.maximum(amplitude, .01*scale)
            def residual(x):
                parameters = x.reshape(nsource, 4)
                coefficients, derivative = double_couple_tensors(parameters, bases)
                data_residual = scaled_G @ coefficients.ravel()-scaled_d
                moment_residual = scaled_damping*parameters[:, 0]
                if pairs is not None:
                    edge = _horizontal_terms(parameters[:, 0], derivative, pairs)
                    return np.r_[data_residual, moment_residual, scaled_smoothing*edge.ravel()]
                return np.r_[data_residual, moment_residual]

            def jacobian(x):
                parameters = x.reshape(nsource, 4)
                _, derivative = double_couple_tensors(parameters, bases)
                blocks = None
                if pairs is not None:
                    blocks = tuple(scaled_smoothing*block for block in
                                   _horizontal_terms(parameters[:, 0], derivative, pairs, jacobian=True))
                return _jacobian_operator(scaled_G, derivative, scaled_damping, pairs, blocks)
            lower = np.tile([0, -np.inf, -np.inf, -np.inf], nsource)
            parameter_scale = np.tile([max(scale, 1e-12), 1, 1, 1], nsource)
            if moment_penalty == 'l1':
                def objective(q):
                    value, gradient = _l1_objective(
                        (q*parameter_scale).reshape(nsource, 4), scaled_G, scaled_d,
                        bases, scaled_l1, scaled_smoothing, pairs)
                    return value, gradient*parameter_scale

                result = minimize(objective, p.ravel()/parameter_scale, jac=True,
                    method='L-BFGS-B', bounds=list(zip(lower, np.full(len(lower), np.inf))),
                    options={'maxiter': max_nfev, 'maxfun': max_nfev, 'ftol': ftol, 'gtol': gtol})
                projected_gradient = result.jac.copy()
                at_bound = (result.x <= lower) & (projected_gradient > 0)
                projected_gradient[at_bound] = 0
                optimality = np.max(abs(projected_gradient))
                cost = result.fun
                p = (result.x*parameter_scale).reshape(nsource, 4)
            else:
                result = least_squares(residual, p.ravel(), jac=jacobian, bounds=(lower, np.inf),
                    method='trf', tr_solver='lsmr', tr_options={'atol': 1e-12, 'btol': 1e-12},
                    x_scale=parameter_scale, loss='linear',
                    max_nfev=max_nfev, ftol=ftol, xtol=xtol, gtol=gtol)
                p = result.x.reshape(nsource, 4)
                cost, optimality = result.cost, result.optimality
            optimization_wall = perf_counter()-optimization_wall
            coefficients, derivative = double_couple_tensors(p, bases)
            item = dict(start=start, success=bool(result.success), status=int(result.status),
                        message=result.message, nfev=int(result.nfev), njev=int(result.njev),
                        cost=float(cost*residual_scale**2), scaled_optimality=float(optimality),
                        residual_l2=float(np.linalg.norm(G @ coefficients.ravel()-d)),
                        moment_l2=float(np.linalg.norm(p[:, 0])))
        if moment_penalty == 'l1':
            item.update(moment_l1=float(p[:, 0].sum()), active_source_count=int(np.count_nonzero(p[:, 0])),
                        iterations=0 if not d.any() else int(result.nit))
        if pairs is not None:
            item['horizontal_mechanism_l2'] = float(np.linalg.norm(_horizontal_terms(p[:, 0], derivative, pairs)))
        item.update(wall_seconds=perf_counter()-start_wall, cpu_seconds=process_time()-start_cpu,
                    initialization_wall_seconds=initialization_wall,
                    optimization_wall_seconds=0. if not d.any() else optimization_wall)
        diagnostics.append(item)
        if verbose:
            print('DC start %d/%d: success=%s; %d evaluations; %.3f s wall; cost %.8g; %s' %
                  (start+1, starts, item['success'], item['nfev'], item['wall_seconds'],
                   item['cost'], item['message']), flush=True)
        if item['success']:
            candidates.append((item['cost'], start, coefficients.copy()))
    if not candidates:
        raise DCConvergenceError(diagnostics)
    _, selected, coefficients = min(candidates, key=lambda candidate: candidate[0])
    tensor_weights = np.sqrt([1, 2, 2, 1, 2, 1])
    model_norm = np.linalg.norm(coefficients*tensor_weights)
    prediction_norm = np.linalg.norm(G @ coefficients.ravel())
    for _, start, candidate in candidates:
        diagnostics[start]['tensor_difference_to_selected'] = (
            float(np.linalg.norm((candidate-coefficients)*tensor_weights)/model_norm) if model_norm else None)
        diagnostics[start]['prediction_difference_to_selected'] = (
            float(np.linalg.norm(G @ (candidate-coefficients).ravel())/prediction_norm) if prediction_norm else None)
    report = {'selected_start': selected, 'starts': diagnostics,
        'penalty': 'scalar_moment_'+moment_penalty, 'moment_reference_nm': REFERENCE_MOMENT_NM,
        'initial_dc_ned': base.tolist(), 'seed': seed, 'residual_scale': float(residual_scale),
        'initialization_method': 'damped_nnls' if moment_penalty == 'l2' else 'nnls',
        'optimizer': 'TRF/LSMR' if moment_penalty == 'l2' else 'L-BFGS-B',
        'wall_seconds': perf_counter()-solve_wall, 'cpu_seconds': process_time()-solve_cpu,
        'global_optimum_established': not d.any()}
    if pairs is not None:
        report.update(horizontal_smoothing=float(horizontal_smoothing), horizontal_pair_count=len(pairs))
    return coefficients.ravel()[:, None], report


def _run_inversion(config, sources, observations, matrix, initial_tensor, provenance, source_name,
                   horizontal_pairs):
    """TOML adapter; caller owns the project lock and immutable run snapshots."""
    from mudpy import geotomo
    from mudpy.geotomo_config import _hash, _json_bytes
    from mudpy.geotomo_io import write_predictions
    data = observations['los_m']
    if matrix.shape != (len(data), len(sources)*6):
        raise ValueError('Pure DC requires the full six-component matrix basis')
    weights = 1/observations['sigma_m'] if config.weighting == 'diagonal' else np.ones(len(data))
    WG, wd = matrix*weights[:, None], data*weights
    options = {key: value for key, value in config.double_couple.items()
               if key not in ('initial_tensor', 'initial_catalog', 'horizontal_neighbor_radius_km')}
    smoothing = {}
    if horizontal_pairs is not None:
        smoothing = {'horizontal_smoothing': options['horizontal_smoothing'],
                     'horizontal_neighbor_radius_km': config.double_couple['horizontal_neighbor_radius_km'],
                     'horizontal_pair_count': len(horizontal_pairs)}
    project = config.project
    home, name = str(project.parent)+'/', project.name
    for member, damping in enumerate(config.damping):
        stem = '%s.%04d' % (config.name, member)
        diagnostic_file = project/'output/inverse_models/models'/(stem+'.solver.json')
        try:
            coefficients, diagnostic = solve(WG, wd, damping, initial_tensor, verbose=True,
                                             horizontal_pairs=horizontal_pairs, **options)
        except DCConvergenceError as error:
            with diagnostic_file.open('xb') as stream:
                stream.write(_json_bytes({'success': False, 'run_name': config.name,
                    'member': member, 'damping': damping, 'starts': error.diagnostics, **smoothing}))
            raise
        selected = diagnostic['starts'][diagnostic['selected_start']]
        diagnostic.update(success=True, run_name=config.name, member=member, damping=damping, **smoothing)
        with diagnostic_file.open('xb') as stream:
            stream.write(_json_bytes(diagnostic))
        prediction = matrix @ coefficients[:, 0]
        total, _, magnitude = geotomo.get_moment(home, name, source_name, coefficients, False)
        vr = np.full(5, np.nan)
        vr[4] = 100*(1-abs(data-prediction).sum()/abs(data).sum()) if abs(data).sum() else np.nan
        geotomo.write_log(home, name, config.name, member, damping, selected['residual_l2'],
            selected['moment_l2'], vr, total, magnitude, config.earth_model, source_name,
            config.matrix_name, config.observations, provenance=dict(provenance,
                Solver='pure DC; '+diagnostic['optimizer'], Penalty=diagnostic['penalty'],
                Solver_diagnostics=diagnostic_file.name, **smoothing))
        geotomo.write_model(home, name, config.name, source_name, config.earth_model,
                            coefficients, member, False, ())
        write_predictions(project/'output/inverse_models/statics'/(stem+'.predictions.h5'), observations,
            prediction, {'run_name': config.name, 'damping': damping, 'provenance': provenance,
                         'source_model': 'double_couple', 'solver_sha256': _hash(diagnostic_file)})
        print('DC member %04d: %d/%d starts converged; selected %d; residual L2 %.6g' %
              (member, sum(item['success'] for item in diagnostic['starts']), len(diagnostic['starts']),
               diagnostic['selected_start'], selected['residual_l2']), flush=True)
