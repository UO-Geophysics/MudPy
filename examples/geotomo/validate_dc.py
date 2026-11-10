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
    args = parser.parse_args()
    study(args.output, args.matrix)
