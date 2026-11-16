"""Controlled pure-DC recovery; optional pre-existing synthetic six-basis matrix."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from obspy.imaging.scripts.mopad import MomentTensor

from mudpy.analysis import decompose_geotomo_tensors
from mudpy.geotomo_dc import solve
from mudpy import analysis, geotomo_dc


def sparse_study(output, config_file):
    """Known sparse normal sources in the full Norcia observation/source geometry."""
    import shutil
    import os
    from time import perf_counter
    from scipy.optimize import minimize
    from sklearn.linear_model import Lasso
    from mudpy import geotomo_config as config
    from mudpy.geotomo_io import REFERENCE_MOMENT_NM

    output = Path(output).expanduser().resolve()
    output.mkdir()
    settings = config.RunConfig.load(config_file)
    sources, observations, geometry, hashes = config._inputs(settings)
    code = config._code_identity()
    executables = {name: {'path': path, 'sha256': config._hash(path)} if path else None
                   for name, path in ((name, shutil.which(name)) for name in ('fk.pl', 'st_fk', 'trav'))}
    specification = config._signatures(settings, sources, geometry, code, executables)[2]
    config._require_cache(settings.project, specification)
    matrix = np.load(specification[2][0], allow_pickle=False)
    if matrix.shape != (len(observations['ids']), 6*len(sources)):
        raise ValueError('Expected a full six-component matrix matching source/receiver geometry')
    pairs = geotomo_dc._horizontal_neighbors(sources, 4.5)
    planted_ids = [94, 95, 96]
    index = {int(row[0]): k for k, row in enumerate(sources)}
    planted = np.array([index[source_id] for source_id in planted_ids])
    moments = np.array([80., 100., 45.])
    reference = MomentTensor([155, 45, -90]).get_M(system='NED')[np.triu_indices(3)]
    reference /= np.linalg.norm(reference*np.sqrt([.5, 1, 1, .5, 1, .5]))
    truth = np.zeros((len(sources), 6))
    for row, amplitude, dip in zip(planted, moments, (50, 45, 40)):
        tensor = MomentTensor([155, dip, -90]).get_M(system='NED')[np.triu_indices(3)]
        truth[row] = amplitude*tensor/np.linalg.norm(tensor*np.sqrt([.5, 1, 1, .5, 1, .5]))
    clean = matrix @ truth.ravel()
    sigma = .002
    rng = np.random.default_rng(261112)
    noisy = clean+rng.normal(0, sigma, len(clean))
    np.savez(output/'synthetic_inputs.npz', sources=sources, receiver_coordinates=observations['coordinates'],
             truth=truth, clean=clean, noisy=noisy, horizontal_pairs=pairs)

    # Independent convex reference with a common fixed mechanism.
    amplitude_matrix = matrix.reshape(len(matrix), len(sources), 6) @ reference
    fixed_truth = np.zeros(len(sources))
    fixed_truth[planted] = moments
    fixed_data = amplitude_matrix @ fixed_truth
    weight = .01
    start = perf_counter()
    lasso = Lasso(alpha=weight/(2*len(matrix)), positive=True, fit_intercept=False,
                  max_iter=50000, tol=1e-9).fit(amplitude_matrix, fixed_data)
    lasso_seconds = perf_counter()-start
    if lasso.dual_gap_ > 1e-9*np.mean(fixed_data**2):
        raise RuntimeError('The fixed-mechanism Lasso reference did not reach its dual-gap tolerance')
    scale = max(np.linalg.norm(fixed_data)/np.linalg.norm(amplitude_matrix, axis=0).max(), 1.)
    residual_scale = np.linalg.norm(fixed_data)
    def amplitude_objective(q):
        amplitude = q*scale
        residual = amplitude_matrix @ amplitude-fixed_data
        value = .5*(residual @ residual+weight*amplitude.sum())/residual_scale**2
        gradient = (amplitude_matrix.T @ residual+.5*weight)*scale/residual_scale**2
        return value, gradient
    start = perf_counter()
    bounded = minimize(amplitude_objective, np.zeros(len(sources)), method='L-BFGS-B',
        jac=True, bounds=[(0, None)]*len(sources),
        options={'maxiter': 20000, 'maxfun': 25000, 'ftol': 1e-13, 'gtol': 1e-9})
    bounded_seconds = perf_counter()-start
    if not bounded.success:
        raise RuntimeError('Fixed-mechanism reference failed: '+bounded.message)
    reference_objective = amplitude_objective(lasso.coef_/scale)[0]
    if not np.isclose(reference_objective, bounded.fun, rtol=1e-6, atol=1e-10):
        raise RuntimeError('Fixed-mechanism L-BFGS-B disagrees with the Lasso objective')
    control = {'l1_weight': weight, 'lasso_wall_seconds': lasso_seconds, 'lasso_iterations': lasso.n_iter_,
               'lbfgsb_wall_seconds': bounded_seconds, 'lasso_dual_gap': float(lasso.dual_gap_),
               'lasso_objective': float(reference_objective*residual_scale**2),
               'lbfgsb_objective': float(bounded.fun*residual_scale**2),
               'prediction_relative_difference': float(np.linalg.norm(amplitude_matrix @
                   (bounded.x*scale-lasso.coef_))/np.linalg.norm(fixed_data))}
    print('Fixed-mechanism Lasso comparison:', json.dumps(control), flush=True)
    cases = []
    specifications = [
        ('clean_l1_weak', clean, .001, 0., [155, 45, -90]),
        ('noisy_l1_weak', noisy, .001, 0., [155, 45, -90]),
        ('noisy_l1_weak_hsmooth', noisy, .001, .003, [155, 45, -90]),
        ('noisy_l1', noisy, .01, 0., [155, 45, -90]),
        ('noisy_l1_hsmooth', noisy, .01, .003, [155, 45, -90]),
        ('noisy_l1_hsmooth_strong', noisy, .01, .01, [155, 45, -90]),
        ('noisy_l1_strong_hsmooth', noisy, .05, .003, [155, 45, -90]),
        ('noisy_l1_misaligned_start', noisy, .001, 0., [110, 65, -45]),
    ]
    for label, data, weight, smoothing, initial_sdr in specifications:
        print(f'Starting {label}: L1 {weight:g}, horizontal smoothing {smoothing:g}', flush=True)
        initial_tensor = MomentTensor(initial_sdr).get_M(system='NED')[np.triu_indices(3)]
        try:
            coefficients, diagnostic = solve(
                matrix, data, weight, initial_tensor, moment_penalty='l1', starts=2, seed=71,
                max_nfev=20000, ftol=1e-9, gtol=1e-7,
                horizontal_smoothing=smoothing, horizontal_pairs=pairs)
        except geotomo_dc.DCConvergenceError as error:
            (output/(label+'.failure.json')).write_text(json.dumps(error.diagnostics, indent=2)+'\n')
            raise
        recovered = coefficients.reshape(-1, 6)
        amplitude = np.linalg.norm(recovered*np.sqrt([.5, 1, 1, .5, 1, .5]), axis=1)
        prediction = matrix @ recovered.ravel()
        parts = decompose_geotomo_tensors(recovered)
        if np.max(parts['moment_nm']['iso']) > 1e-10 or np.max(parts['moment_nm']['clvd']) > 1e-10:
            raise RuntimeError('Recovered source violates the pure-DC constraint')
        n95 = int(np.searchsorted(np.cumsum(np.sort(amplitude)[::-1]), .95*amplitude.sum())+1)
        row = {'case': label, 'l1_weight': weight, 'horizontal_smoothing': smoothing,
               'initial_sdr_deg': initial_sdr,
               'residual_rms_mm': float(np.sqrt(np.mean((prediction-data)**2))*1000),
               'clean_prediction_rms_mm': float(np.sqrt(np.mean((prediction-clean)**2))*1000),
               'moment_sum_nm': float(amplitude.sum()*REFERENCE_MOMENT_NM),
               'moment_sum_ratio': float(amplitude.sum()/moments.sum()),
               'exact_nonzero_sources': int(np.count_nonzero(amplitude)), 'sources_carrying_95pct_moment': n95,
               'moment_on_planted_nodes_fraction': float(amplitude[planted].sum()/amplitude.sum()),
               'tensor_relative_error': float(np.linalg.norm(recovered-truth)/np.linalg.norm(truth)),
               'diagnostics': diagnostic}
        cases.append(row)
        np.savez(output/(label+'.npz'), coefficients=recovered, prediction=prediction)
        (output/(label+'.json')).write_text(json.dumps(row, indent=2, allow_nan=False)+'\n')
        print(f'Finished {label}: RMS {row["residual_rms_mm"]:.3f} mm, '
              f'{n95} sources carry 95% of moment, {diagnostic["wall_seconds"]:.2f} s', flush=True)
    config._verify_inputs(hashes)
    config._require_cache(settings.project, specification)
    report = {'matrix_shape': list(matrix.shape), 'source_ids': sources[:, 0].astype(int).tolist(),
              'planted_ids': planted_ids, 'planted_sdr_deg': [[155, dip, -90] for dip in (50, 45, 40)],
              'true_moment_nm': float(moments.sum()*REFERENCE_MOMENT_NM),
              'noise_sigma_mm': sigma*1000, 'noise_seed': 261112, 'initial_sdr_deg': [155, 45, -90],
              'starts_per_case': 2, 'horizontal_pair_count': len(pairs), 'neighbor_radius_km': 4.5,
              'fixed_mechanism_control': control, 'cases': cases, 'input_sha256': hashes,
              'matrix_sha256': config._hash(specification[2][0]),
              'config_path': str(settings.config_path), 'config_sha256': config._hash(settings.config_path),
              'code_sha256': {str(Path(path).resolve()): config._hash(path)
                              for path in (__file__, analysis.__file__, geotomo_dc.__file__)},
              'versions': {'numpy': np.__version__, 'scipy': scipy.__version__},
              'blas_thread_settings': {name: os.environ.get(name) for name in
                                      ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')},
              'interpretation': 'On-grid synthetic recovery using all candidate locations and one LOS geometry.'}
    (output/'sparse_validation.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    plot_sparse_study(output)
    return report


def plot_sparse_study(output):
    """Plot truth, support, DC mechanisms and LOS recovery for the sparse study."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.offsetbox import AnnotationBbox, DrawingArea
    from obspy.imaging.mopad_wrapper import beach
    from mudpy.geotomo_io import REFERENCE_MOMENT_NM

    output = Path(output)
    report = json.loads((output/'sparse_validation.json').read_text())
    inputs = np.load(output/'synthetic_inputs.npz', allow_pickle=False)
    sources, coordinates, truth = inputs['sources'], inputs['receiver_coordinates'], inputs['truth']
    cases = report['cases']
    selected = ['noisy_l1_weak', 'noisy_l1_weak_hsmooth', 'noisy_l1_misaligned_start']
    models = [truth]+[np.load(output/(name+'.npz'))['coefficients'] for name in selected]
    labels = ['Truth', 'L1 0.001', 'L1 0.001 + smoothing 0.003', 'L1 0.001 · oblique start']
    source_labels = ['Truth', 'λ = .001, μ = 0', 'λ = .001, μ = .003', 'Oblique start, λ = .001']
    moments = [np.linalg.norm(model*np.sqrt([.5, 1, 1, .5, 1, .5]), axis=1) for model in models]
    maximum = max(a.max() for a in moments)
    folder = output/'plots'
    folder.mkdir(exist_ok=True)
    files = []
    def save(figure, name):
        path = folder/(name+'.png')
        figure.savefig(path, dpi=160, bbox_inches='tight')
        files.append(str(path))
        plt.close(figure)
    def map_axis(axis):
        axis.set(xlim=(sources[:, 1].min()-.02, sources[:, 1].max()+.02),
                 ylim=(sources[:, 2].min()-.02, sources[:, 2].max()+.02),
                 xlabel='Longitude (degrees)', ylabel='Latitude (degrees)')
        axis.set_aspect(1/np.cos(np.deg2rad(sources[:, 2].mean())))
        axis.grid(alpha=.2)
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 9}):
        depths = np.unique(sources[:, 3])
        figure, axes = plt.subplots(len(models), len(depths), figsize=(17, 12), constrained_layout=True)
        for row, (label, amplitude) in enumerate(zip(source_labels, moments)):
            for col, depth in enumerate(depths):
                mask = sources[:, 3] == depth
                axis = axes[row, col]
                image = axis.scatter(sources[mask, 1], sources[mask, 2], c=amplitude[mask]*REFERENCE_MOMENT_NM,
                    s=28, cmap='magma', norm=Normalize(0, maximum*REFERENCE_MOMENT_NM))
                planted = mask & (moments[0] > 0)
                axis.scatter(sources[planted, 1], sources[planted, 2], marker='s',
                              s=75, facecolors='none', edgecolors='cyan', linewidths=1)
                map_axis(axis)
                axis.set_title(f'{label}\n{depth:g} km')
        figure.colorbar(image, ax=list(axes.flat), shrink=.65, label='Source moment (N m)')
        figure.suptitle('Sparse DC synthetic recovery · cyan squares mark the planted sources')
        save(figure, 'source_recovery_by_depth')

        figure, axes = plt.subplots(2, 2, figsize=(12, 11))
        figure.subplots_adjust(top=.90, bottom=.075, hspace=.35, wspace=.25)
        for axis, label, model, amplitude in zip(axes.flat, labels, models, moments):
            mask = (sources[:, 3] == 4) & (amplitude > 0)
            axis.scatter(sources[sources[:, 3] == 4, 1], sources[sources[:, 3] == 4, 2], s=6, color='.8')
            for k in np.flatnonzero(mask)[np.argsort(amplitude[mask])]:
                width = 34*np.sqrt(amplitude[k]/maximum)
                if width < .5:
                    continue
                tensor = model[k][[0, 3, 5, 1, 2, 4]]
                tensor = tensor/np.abs(tensor).max()
                drawing = DrawingArea(width, width, 0, 0)
                drawing.add_artist(beach(tensor, mopad_basis='NED', xy=(width/2, width/2),
                    width=width, linewidth=.4, facecolor='black', bgcolor='white', size=100))
                axis.add_artist(AnnotationBbox(drawing, sources[k, 1:3], frameon=False, pad=0,
                                               box_alignment=(.5, .5), annotation_clip=True))
            map_axis(axis)
            axis.set_title(label)
        figure.suptitle('DC mechanisms at 4 km · recovered from data with 2 mm noise\n'
                       'Common symbol area proportional to source moment')
        save(figure, 'mechanism_recovery_4km')

        figure, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)
        limit = np.max(abs(inputs['noisy']))*1000
        for row, name in enumerate(('noisy_l1_weak', 'noisy_l1_weak_hsmooth')):
            prediction = np.load(output/(name+'.npz'))['prediction']
            for col, (values, title) in enumerate(((inputs['noisy'], 'Synthetic observations'),
                                                   (prediction, 'Prediction'),
                                                   (inputs['noisy']-prediction, 'Residual'))):
                axis = axes[row, col]
                bound = 10 if col == 2 else limit
                image = axis.scatter(coordinates[:, 0], coordinates[:, 1], c=values*1000,
                    s=8, cmap='RdBu_r', vmin=-bound, vmax=bound)
                axis.set(xlabel='Longitude (degrees)', ylabel='Latitude (degrees)',
                          title=name.replace('_', ' ')+'\n'+title)
                figure.colorbar(image, ax=axis, label='mm', shrink=.8)
        figure.suptitle('2 mm noise · free DC inversion over all 315 candidates')
        save(figure, 'los_fit_and_residual')

        figure, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
        positions = np.arange(len(cases))
        short_labels = [('Clean' if case['case'].startswith('clean') else
                         'Oblique start' if case['case'].endswith('misaligned_start') else 'Noisy')+
                        '\nλ = %g\nμ = %g' % (case['l1_weight'], case['horizontal_smoothing'])
                        for case in cases]
        measures = [('residual_rms_mm', 'Fit RMS (mm)'), ('sources_carrying_95pct_moment', 'Sources carrying 95% of moment'),
                    ('moment_sum_ratio', 'Summed moment / truth')]
        for axis, (key, ylabel) in zip(axes.flat, measures):
            if key == 'sources_carrying_95pct_moment':
                axis.bar(positions-.18, [case['exact_nonzero_sources'] for case in cases],
                          width=.36, label='Nonzero sources')
                axis.bar(positions+.18, [case[key] for case in cases], width=.36, label='95% of moment')
                ylabel = 'Source count'
            else:
                axis.bar(positions, [case[key] for case in cases])
            axis.set(xticks=positions, xticklabels=short_labels, ylabel=ylabel)
            axis.grid(axis='y', alpha=.2)
        axes[0, 0].axhline(report['noise_sigma_mm'], color='black', linestyle='--', label='Noise sigma')
        axes[0, 0].legend()
        axes[0, 1].axhline(3, color='black', linestyle='--', label='Three true sources')
        axes[0, 1].legend()
        axes[1, 0].axhline(1, color='black', linestyle='--')
        axes[1, 1].bar(positions, [case['diagnostics']['wall_seconds'] for case in cases])
        axes[1, 1].set(xticks=positions, xticklabels=short_labels, ylabel='Solver wall time (s)')
        figure.suptitle('L1 synthetic validation · each case has two starts')
        save(figure, 'fit_sparsity_moment_and_timing')
    report['figures'] = files
    report['plot_driver_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (output/'sparse_validation.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print('Saved plots:', folder, flush=True)


def study(output, matrix_file=None):
    output = Path(output).expanduser().resolve()
    output.mkdir()  # Never overwrite a previous validation study.
    rng = np.random.default_rng(261110)
    if matrix_file is None:
        G = rng.normal(size=(98, 18))*.01
        provenance = {'matrix': 'deterministic dense full-rank control', 'seed': 261110}
    else:
        matrix_file = Path(matrix_file).expanduser().resolve()
        checksum = hashlib.sha256(matrix_file.read_bytes()).hexdigest()
        full = np.load(matrix_file, allow_pickle=False)
        if full.ndim != 2 or full.shape[1] % 6 or full.shape[1] < 18:
            raise ValueError('Validation matrix must have at least three six-component sources')
        indices = np.linspace(0, full.shape[1]//6-1, 3, dtype=int)
        G = np.concatenate([full[:, 6*k:6*k+6] for k in indices], axis=1)
        provenance = {'matrix': str(matrix_file), 'sha256': checksum,
                      'source_indices_zero_based': indices.tolist()}
    # Independent SDR construction: constant strike, progressively shallower dip.
    truth = np.array([MomentTensor([155, dip, -90]).get_M(system='NED')[np.triu_indices(3)]*a
                      for dip, a in zip((65, 45, 25), (2, 1.5, 1))])
    initial = MomentTensor([155, 45, -90]).get_M(system='NED')[np.triu_indices(3)]
    clean = G @ truth.ravel()
    sigma = np.sqrt(np.mean(clean**2))*.01
    noisy = clean+rng.normal(0, sigma, len(clean))
    damping_scale = np.linalg.svd(G, compute_uv=False)[0]
    rows = []
    for label, data, damping in [('noiseless', clean, 0), ('noisy_weak', noisy, damping_scale*.001),
                                 ('noisy_stronger', noisy, damping_scale*.1)]:
        coefficients, diagnostics = solve(G, data, damping, initial, starts=3, seed=71,
                                           max_nfev=1000, gtol=1e-11)
        parts = decompose_geotomo_tensors(coefficients.reshape(-1, 6))
        assert np.max(parts['moment_nm']['iso']) < 1e-12
        assert np.max(parts['moment_nm']['clvd']) < 1e-12
        tensor_error = np.linalg.norm(coefficients[:, 0]-truth.ravel())/np.linalg.norm(truth)
        if label == 'noiseless':
            assert tensor_error < 1e-4, tensor_error
        rows.append({'case': label, 'damping': float(damping), 'tensor_relative_error': float(tensor_error),
                     'residual_rms_m': float(np.sqrt(np.mean((G @ coefficients[:, 0]-data)**2))),
                     'moment_sum_ratio': float(parts['moment_nm']['full'].sum()/4.5),
                     'diagnostics': diagnostics})
        np.save(output/(label+'.coefficients.npy'), coefficients, allow_pickle=False)
    if matrix_file is not None:
        assert hashlib.sha256(matrix_file.read_bytes()).hexdigest() == checksum
    report = {'input': provenance, 'shape': list(G.shape), 'true_sdr_deg': [[155, z, -90] for z in (65, 45, 25)],
              'initial_sdr_deg': [155, 45, -90], 'true_moments_reference_units': [2, 1.5, 1],
              'noise_sigma_m': float(sigma), 'cases': rows,
              'versions': {'numpy': np.__version__, 'scipy': scipy.__version__},
              'code_sha256': {str(Path(path).resolve()): hashlib.sha256(Path(path).read_bytes()).hexdigest()
                              for path in (__file__, analysis.__file__, geotomo_dc.__file__)},
              'limits': ['Known source locations; three active sources; one controlled geometry.',
                         'Passing recovery does not establish resolution for a larger source volume.',
                         'Multistart explores local solutions, not a proof of global optimality.']}
    (output/'dc_validation.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(rows, indent=2, allow_nan=False), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--matrix', type=Path, help='Six-component matrix from a synthetic forward study')
    parser.add_argument('--sparse-config', type=Path, help='Norcia volume TOML for full-grid L1 synthetic validation')
    args = parser.parse_args()
    if args.sparse_config:
        if args.matrix:
            parser.error('--matrix and --sparse-config are mutually exclusive')
        sparse_study(args.output, args.sparse_config)
    else:
        study(args.output, args.matrix)
