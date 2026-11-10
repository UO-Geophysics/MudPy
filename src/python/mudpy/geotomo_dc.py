"""Pure double-couple nonlinear least squares on the existing six-tensor matrix."""
import numpy as np

from mudpy.analysis import decompose_geotomo_tensors, double_couple_tensors
from mudpy.geotomo_io import REFERENCE_MOMENT_NM


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


def _jacobian_operator(matrix, derivative, damping):
    """One immutable derivative per operator; no shared state between starts."""
    from scipy.sparse.linalg import LinearOperator
    ndata, nsource = len(matrix), len(derivative)
    def matvec(vector):
        h = np.asarray(vector).reshape(nsource, 4)
        data = matrix @ np.einsum('nck,nk->nc', derivative, h).ravel()
        return np.r_[data, damping*h[:, 0]]
    def rmatvec(vector):
        vector = np.asarray(vector).ravel()
        result = np.einsum('nck,nc->nk', derivative, (matrix.T @ vector[:ndata]).reshape(nsource, 6))
        result[:, 0] += damping*vector[ndata:]
        return result.ravel()
    return LinearOperator((ndata+nsource, 4*nsource), matvec=matvec, rmatvec=rmatvec, dtype=float)


def solve(matrix, data, damping, initial_tensor, *, starts=3, seed=0, perturbation_deg=30.,
          max_nfev=500, ftol=1e-9, xtol=1e-9, gtol=1e-9):
    """Minimize ||G*m(a,r)-d||² + damping² ||a||², with a >= 0.

    G is in metres per 10**16.6 Nm coefficient; data/G may already be weighted.
    Amplitudes a are scalar moments divided by that reference, rotations radians.
    TRF/LSMR uses analytic matrix-free Jacobians. Starts run sequentially, each
    from the same catalog DC plus seeded independent source rotations, followed
    by bounded linear amplitude initialization. No global-optimality claim.
    An unconverged start is recorded but cannot be selected. If all fail, return
    no model: RuntimeError carries the per-start diagnostics for archival.
    """
    from scipy.optimize import least_squares, lsq_linear
    from scipy.spatial.transform import Rotation
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
    if d.any() and not G.any():
        raise ValueError('Nonzero data have no responsive matrix columns')
    # Common residual scaling preserves the objective minimizer and makes
    # convergence tolerances meaningful for both metre-valued and unit tests.
    residual_scale = np.linalg.norm(d) or 1.
    scaled_G, scaled_d, scaled_damping = G/residual_scale, d/residual_scale, damping/residual_scale
    rng = np.random.default_rng(seed)
    diagnostics, candidates = [], []
    for start in range(starts):
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
        augmented = np.vstack((projected, scaled_damping*np.eye(nsource)))
        target = np.r_[scaled_d, np.zeros(nsource)]
        amplitude = lsq_linear(augmented, target, bounds=(0, np.inf), tol=1e-12)
        if not amplitude.success:
            raise RuntimeError('DC amplitude initialization did not converge: '+amplitude.message)
        if not d.any():
            # Known global minimum, including damping=0: no orientation is defined.
            p[:, 0] = 0
            coefficients, _ = double_couple_tensors(p, bases)
            item = dict(start=start, success=True, status=1, message='Exact zero-data minimum',
                        nfev=0, njev=0, cost=0., scaled_optimality=0., residual_l2=0., moment_l2=0.)
        else:
            # A positive initialization floor frees rotations at zero NNLS nodes.
            # This is not a final lower bound or a prior on recovered moments.
            scale = max(np.linalg.norm(amplitude.x)/np.sqrt(nsource),
                        np.linalg.norm(scaled_d)/max(np.linalg.norm(projected), np.finfo(float).tiny))
            p[:, 0] = np.maximum(amplitude.x, .01*scale)
            def residual(x):
                coefficients, _ = double_couple_tensors(x.reshape(nsource, 4), bases)
                return np.r_[scaled_G @ coefficients.ravel()-scaled_d, scaled_damping*x.reshape(nsource, 4)[:, 0]]
            def jacobian(x):
                _, derivative = double_couple_tensors(x.reshape(nsource, 4), bases)
                return _jacobian_operator(scaled_G, derivative, scaled_damping)
            lower = np.tile([0, -np.inf, -np.inf, -np.inf], nsource)
            result = least_squares(residual, p.ravel(), jac=jacobian, bounds=(lower, np.inf),
                method='trf', tr_solver='lsmr', tr_options={'atol': 1e-12, 'btol': 1e-12},
                x_scale=np.tile([max(scale, 1e-12), 1, 1, 1], nsource), loss='linear',
                max_nfev=max_nfev, ftol=ftol, xtol=xtol, gtol=gtol)
            p = result.x.reshape(nsource, 4)
            coefficients, _ = double_couple_tensors(p, bases)
            item = dict(start=start, success=bool(result.success), status=int(result.status),
                        message=result.message, nfev=int(result.nfev), njev=int(result.njev),
                        cost=float(result.cost*residual_scale**2), scaled_optimality=float(result.optimality),
                        residual_l2=float(np.linalg.norm(G @ coefficients.ravel()-d)),
                        moment_l2=float(np.linalg.norm(p[:, 0])))
        diagnostics.append(item)
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
    return coefficients.ravel()[:, None], {'selected_start': selected, 'starts': diagnostics,
        'penalty': 'scalar_moment_l2', 'moment_reference_nm': REFERENCE_MOMENT_NM,
        'initial_dc_ned': base.tolist(), 'seed': seed, 'residual_scale': float(residual_scale),
        'global_optimum_established': not d.any()}


def _run_inversion(config, sources, observations, matrix, initial_tensor, provenance, source_name):
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
               if key not in ('initial_tensor', 'initial_catalog')}
    project = config.project
    home, name = str(project.parent)+'/', project.name
    for member, damping in enumerate(config.damping):
        stem = '%s.%04d' % (config.name, member)
        diagnostic_file = project/'output/inverse_models/models'/(stem+'.solver.json')
        try:
            coefficients, diagnostic = solve(WG, wd, damping, initial_tensor, **options)
        except DCConvergenceError as error:
            with diagnostic_file.open('xb') as stream:
                stream.write(_json_bytes({'success': False, 'run_name': config.name,
                    'member': member, 'damping': damping, 'starts': error.diagnostics}))
            raise
        selected = diagnostic['starts'][diagnostic['selected_start']]
        diagnostic.update(success=True, run_name=config.name, member=member, damping=damping)
        with diagnostic_file.open('xb') as stream:
            stream.write(_json_bytes(diagnostic))
        prediction = matrix @ coefficients[:, 0]
        total, _, magnitude = geotomo.get_moment(home, name, source_name, coefficients, False)
        vr = np.full(5, np.nan)
        vr[4] = 100*(1-abs(data-prediction).sum()/abs(data).sum()) if abs(data).sum() else np.nan
        geotomo.write_log(home, name, config.name, member, damping, selected['residual_l2'],
            selected['moment_l2'], vr, total, magnitude, config.earth_model, source_name,
            config.matrix_name, config.observations, provenance=dict(provenance,
                Solver='pure DC nonlinear least squares; TRF/LSMR', Penalty='scalar_moment_l2',
                Solver_diagnostics=diagnostic_file.name))
        geotomo.write_model(home, name, config.name, source_name, config.earth_model,
                            coefficients, member, False, ())
        write_predictions(project/'output/inverse_models/statics'/(stem+'.predictions.h5'), observations,
            prediction, {'run_name': config.name, 'damping': damping, 'provenance': provenance,
                         'source_model': 'double_couple', 'solver_sha256': _hash(diagnostic_file)})
        print('DC member %04d: %d/%d starts converged; selected %d; residual L2 %.6g' %
              (member, sum(item['success'] for item in diagnostic['starts']), len(diagnostic['starts']),
               diagnostic['selected_start'], selected['residual_l2']), flush=True)
