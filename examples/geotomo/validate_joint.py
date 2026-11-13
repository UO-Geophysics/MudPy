"""Controlled joint DC/background recovery in a verified physical observation geometry."""
import argparse
import json
from pathlib import Path
import shutil
from time import perf_counter, process_time

import numpy as np
import matplotlib.pyplot as plt

from mudpy import analysis, geotomo_dc as dc, geotomo_io as io, geotomo_config as cfg


def main(config_file, output):
    wall, cpu = perf_counter(), process_time()
    output = Path(output).resolve()
    output.mkdir()
    config = cfg.RunConfig.load(config_file)
    if not config.joint:
        raise ValueError('Provide a joint TOML with an already verified unweighted matrix')
    sources, observed, geometry, hashes = cfg._inputs(config)
    code = cfg._code_identity()
    executables = {name: {'path': path, 'sha256': cfg._hash(path)} if path else None for name, path in
                   ((name, shutil.which(name)) for name in ('fk.pl', 'st_fk', 'trav'))}
    specification = cfg._signatures(config, sources, geometry, code, executables)[2]
    cache = cfg._require_cache(config.project, specification)
    physical = np.load(specification[2][0], mmap_mode='r', allow_pickle=False)
    layer = np.flatnonzero(np.isclose(sources[:, 3], 4.))
    distance = ((sources[layer, 1]-13.18)*np.cos(np.deg2rad(42.82)))**2+(sources[layer, 2]-42.82)**2
    selected = layer[np.argsort(distance)[:3]]
    columns = (6*selected[:, None]+np.arange(6)).ravel()
    G = physical[:, columns].copy()
    initial = io.read_initial_tensor(config.project/'data/model_info'/config.double_couple['initial_catalog'])['tensor_nm']
    base = dc.dc_reference(initial)
    parameters = np.array([[100., .10, -.08, .03], [80., -.07, .13, .1], [0., 0., 0., 0.]])
    truth, _ = analysis.double_couple_tensors(parameters, base)
    beta = np.array([.01, -.004, .006, .002, -.003, .001])
    count = len(observed['imagery']['ids'])
    basis = observed['background_basis']
    source = G @ truth.ravel()
    clean = source.copy()
    clean[:count] += basis @ beta
    rng = np.random.default_rng(20261112)
    noise_sigma = np.r_[np.full(count, .015), observed['gnss']['sigma_m'].ravel()]
    noisy = clean+rng.normal(size=len(clean))*noise_sigma
    np.savez(output/'inputs.npz', sources=sources[selected], matrix=G, truth_coefficients=truth,
              background_coefficients_m=beta, background_basis=basis, clean=clean, noisy=noisy,
              noise_sigma_m=noise_sigma, imagery_coordinates=observed['imagery']['coordinates'],
              gnss_coordinates=observed['gnss']['coordinates'])
    cases = []
    for label, data, damping in [('clean', clean, 0.), ('noisy', noisy, .0001)]:
        synthetic = dict(observed,
            imagery=dict(observed['imagery'], los_m=data[:count], metadata={'synthetic': True}),
            gnss=dict(observed['gnss'], displacement_m=data[count:].reshape(-1, 3),
                      metadata={'synthetic': True, 'sigma_interpretation': 'known independent Gaussian noise standard deviations'}))
        io.write_joint(output/(label+'.h5'), synthetic['imagery'], synthetic['gnss'], basis,
                        metadata={'synthetic': True, 'no_measured_displacements_used': True, 'LOS_noise_std_m': .015 if label == 'noisy' else 0.})
        for eta in (0., 1.):
            case_name = '%s_eta%d' % (label, eta)
            start = perf_counter()
            coefficients, diagnostic = dc.solve_joint(G, synthetic, damping, initial,
                gnss_weight=eta, gnss_scale_m=.015, moment_penalty='l1', starts=2, seed=71,
                perturbation_deg=30., max_nfev=20000, ftol=1e-13, gtol=1e-9, verbose=True)
            recovered = coefficients.reshape(-1, 6)
            recovered_beta = np.array(diagnostic['joint']['background_coefficients_m'])
            fitted = G @ coefficients[:, 0]
            fitted[:count] += basis @ recovered_beta
            error = float(np.linalg.norm(recovered-truth)/np.linalg.norm(truth))
            background_error = float(abs(recovered_beta-beta).max())
            parts = analysis.decompose_geotomo_tensors(recovered)
            moment = parts['moment_nm']['full']
            if np.any(parts['moment_nm']['iso']+parts['moment_nm']['clvd'] > 1e-8*np.maximum(moment, 1.)):
                raise RuntimeError('Synthetic solution violates the DC constraint')
            if label == 'clean' and (error > 1e-3 or background_error > 1e-4):
                raise RuntimeError('Clean joint/control recovery failed the tensor/background accuracy checks')
            normalized_gps = float(np.sqrt(np.mean(((fitted[count:]-data[count:])/noise_sigma[count:])**2)))
            if label == 'noisy' and eta and normalized_gps > 2:
                raise RuntimeError('Joint synthetic GNSS residuals exceed the known noise-scale sanity bound')
            row = {'case': case_name, 'gnss_weight': eta, 'lambda': damping,
                   'tensor_relative_error': error, 'background_max_error_mm': background_error*1000,
                   'source_moments_nm': (moment*io.REFERENCE_MOMENT_NM).tolist(),
                   'background_coefficients_m': recovered_beta.tolist(),
                   'imagery_residual_rms_mm': float(np.sqrt(np.mean((fitted[:count]-data[:count])**2))*1000),
                   'gnss_residual_rms_neu_mm': (np.sqrt(np.mean((fitted[count:]-data[count:]).reshape(-1, 3)**2, axis=0))*1000).tolist(),
                   'gnss_noise_scaled_rms': normalized_gps, 'wall_seconds': perf_counter()-start,
                   'diagnostics': diagnostic}
            cases.append(row)
            np.savez(output/(case_name+'.npz'), coefficients=recovered, background_coefficients_m=recovered_beta, prediction_m=fitted)
            print(case_name, 'tensor error', error, 'background error mm', background_error*1000, flush=True)
    figure, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for row, label in enumerate(('clean', 'noisy')):
        x = np.arange(3)
        axes[row, 0].bar(x-.25, parameters[:, 0]*io.REFERENCE_MOMENT_NM, width=.25, color='.5', label='Truth')
        for k, case in enumerate(cases[2*row:2*row+2]):
            axes[row, 0].bar(x+k*.25, case['source_moments_nm'], width=.25, label='η = %g' % case['gnss_weight'])
            axes[row, 1].plot(np.arange(6), np.array(case['background_coefficients_m'])*1000, 'o-', label='η = %g' % case['gnss_weight'])
        axes[row, 1].plot(np.arange(6), beta*1000, 'k--', label='Truth')
        axes[row, 0].set(xticks=x, xticklabels=sources[selected, 0].astype(int), xlabel='Source ID', ylabel='Scalar moment (N m)', title=label.title())
        axes[row, 1].set(xticks=np.arange(6), xticklabels=['1', 'x', 'y', 'x²', 'xy', 'y²'], ylabel='Background coefficient (mm)', title=label.title())
    for axis in axes.flat:
        axis.legend(fontsize=8); axis.grid(axis='y', alpha=.2)
    figure.suptitle('Synthetic DC sources and quadratic background · verified Norcia geometry')
    plot = output/'recovery_and_background.png'
    figure.savefig(plot, dpi=180, bbox_inches='tight'); plt.close(figure)
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    x = np.arange(3)
    for k, case in enumerate(cases[2:]):
        axes[0].bar(x+(k-.5)*.35, case['gnss_residual_rms_neu_mm'], width=.35, label='η = %g' % case['gnss_weight'])
        axes[1].bar(k, case['tensor_relative_error'], label='η = %g' % case['gnss_weight'])
    axes[0].set(xticks=x, xticklabels=['North', 'East', 'Up'], ylabel='GNSS residual RMS (mm)')
    axes[1].set(xticks=[0, 1], xticklabels=['Imagery-only control', 'Joint'], ylabel='Tensor relative error')
    for axis in axes: axis.grid(axis='y', alpha=.2)
    axes[0].legend()
    figure.suptitle('Known noisy synthetic data · fixed source locations, free DC mechanisms')
    plot2 = output/'noisy_comparison.png'
    figure.savefig(plot2, dpi=180, bbox_inches='tight'); plt.close(figure)
    cfg._verify_inputs(hashes)
    cfg._require_cache(config.project, specification)
    report = {'cases': cases, 'source_ids': sources[selected, 0].astype(int).tolist(),
              'truth_moment_nm': float(parameters[:, 0].sum()*io.REFERENCE_MOMENT_NM),
              'matrix_shape': list(G.shape), 'background_rank': int(np.linalg.matrix_rank(basis)),
              'reference_matrix_files': cache['files'], 'input_sha256': hashes, 'code': code,
              'script_sha256': cfg._hash(__file__), 'figures': [str(plot), str(plot2)],
              'wall_seconds': perf_counter()-wall, 'cpu_seconds': process_time()-cpu,
              'interpretation': 'Controlled numerical validation in verified elastic geometry; no real displacements inverted.'}
    (output/'validation.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print('Validated synthetic joint study:',output, 'wall',report['wall_seconds'],flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    main(args.config, args.output)
