"""Reproducible CPU tensor/imagery checks using the existing MudPy formats.

Run with --project PATH. PATH must not exist: existing results are preserved.
No MPI, real imagery, GNSS, or optional reference solver is needed at runtime.
"""
import argparse
import hashlib
import sys
import json
from pathlib import Path

import numpy as np
from pyproj import Geod

from mudpy import mt3d


def run_with_toml(project, name, steps, mechanism=None, matrix_mode='assemble', matrix_name=None):
    """Write an explicit experiment configuration for this reference case."""
    source_model = 'fixed_tensor' if mechanism is not None else 'full_tensor'
    text = '''schema_version = 1
[project]
directory = ".."
[run]
name = %s
steps = %s
[inputs]
sources = "cloud.sources"
earth_model = "halfspace.mod"
observations = "imagery.gflist"
[forward]
backend = "fk"
earth_flattening = false
dk = 0.05
pmin = 0.0
pmax = 1.0
kmax = 25.0
[inversion]
source_model = %s
solver = "lstsq"
weighting = "diagonal"
damping = [0.0]
%s
[matrix]
name = %s
mode = %s
[execution]
processes = 1
''' % (json.dumps(name), json.dumps(steps), json.dumps(source_model),
       'fixed_tensor = '+json.dumps(list(mechanism)) if mechanism is not None else '',
       json.dumps(matrix_name or 'toml_'+name), json.dumps(matrix_mode))
    path = project/'configs'/(name+'.toml')
    with path.open('x') as stream:
        stream.write(text)
    return mt3d.run(path)


def validate_single_receiver(project, use_toml=False):
    """Exercise real one-source/one-observation files and project exclusion."""
    single = project / 'single_receiver'
    if use_toml:
        mt3d.init_project(single)
    else:
        single.mkdir()
        for directory in ('data/model_info', 'data/station_info', 'data/statics',
                          'structure', 'GFs/static', 'GFs/matrices',
                          'output/inverse_models/models', 'output/inverse_models/statics'):
            (single / directory).mkdir(parents=True)
    source = np.loadtxt(project / 'data/model_info/cloud.sources', ndmin=2)[:1]
    np.savetxt(single / 'data/model_info/cloud.sources', source,
               fmt=['%d', '%.9f', '%.9f', '%.9f'])
    (single / 'structure/halfspace.mod').write_bytes((project / 'structure/halfspace.mod').read_bytes())
    header, first, *_ = (project / 'data/station_info/imagery.gflist').read_text().splitlines()
    fields = first.split()
    observation = single / 'data/statics' / (fields[0]+'.los')
    look = np.loadtxt(fields[12])[1:]
    np.savetxt(observation, np.r_[0., look])
    fields[12] = str(observation)
    (single / 'data/station_info/imagery.gflist').write_text(header+'\n'+' '.join(fields)+'\n')
    home, name = str(single.parent)+'/', single.name
    mechanism = [1, 0, 0, 1, 0, 1]
    cwd = Path.cwd()
    if use_toml:
        run_with_toml(single, 'prepare', ['greens', 'synthetics'], mechanism)
    else:
        mt3d.inversionGFs(home, name, 'imagery.gflist', 'cloud.sources', 'halfspace.mod',
                         1, 1, .05, 0, 1, 25, 0, 1, True, mechanism, flat=False)
    G = mt3d.getG(home, name, 'cloud.sources', 'halfspace.mod', 'imagery.gflist', False, 'single', True)
    assert G.shape == (1, 1)
    np.savetxt(observation, np.r_[-.5*G[0, 0], look])
    arguments = (home, name, 'contraction', 'cloud.sources', 'halfspace.mod',
                 'imagery.gflist', False, 'single', [0], 1, 'lstsq', True, mechanism)
    if use_toml:
        run_with_toml(single, 'contraction', ['invert'], mechanism)
    else:
        mt3d.run_inversion(*arguments, weight=True)
    result = np.loadtxt(single / 'output/inverse_models/models/contraction.0000.inv', ndmin=2)
    np.testing.assert_allclose(result[:, 4:10]/10**16.6, [-.5*np.array(mechanism)], atol=1e-10)
    with mt3d._project_lock(home, name):
        try:
            mt3d.run_inversion(*arguments)
        except FileExistsError:
            pass  # Expected exclusion before matrix/output mutation.
        else:
            raise AssertionError('Concurrent project ownership was not rejected')
    assert Path.cwd() == cwd
    assert not (single / '.mt3d.lock').exists()
    return float(abs(result[0, 4]/10**16.6+.5)/.5)


def main(project, use_toml=False):
    reference = json.loads(Path(__file__).with_name('okada_reference.json').read_text())
    project = project.resolve()
    if use_toml:
        mt3d.init_project(project)
    else:
        project.mkdir()  # Refuse to overwrite a previous experiment.
        for directory in ('data/model_info', 'data/station_info', 'data/statics',
                          'structure', 'GFs/static', 'GFs/matrices',
                          'output/inverse_models/models', 'output/inverse_models/statics'):
            (project / directory).mkdir(parents=True)
    home, name = str(project.parent)+'/', project.name
    # Deliberately reverse numeric IDs to detect accidental sorting.
    sources = np.array([[2, 13., 43., 3.], [1, 13.02, 43.01, 5.]])
    np.savetxt(project / 'data/model_info/cloud.sources', sources,
               fmt=['%d', '%.9f', '%.9f', '%.9f'])
    vs = 6/np.sqrt(3)
    np.savetxt(project / 'structure/halfspace.mod', [[0, vs, 6, 2.7, 500, 1000]],
               fmt='%.12e')
    # Measurement directions point toward representative ascending/descending
    # sensors. Positive scalar displacement is along the specified NEU vector.
    looks = np.array([[.1, -.6, np.sqrt(.63)], [.1, .6, np.sqrt(.63)]])
    points = np.array(reference['offsets_ne_m'])
    geod = Geod(ellps='WGS84')
    observations, directions, receivers = [], [], []
    for k, (north, east) in enumerate(points):
        azimuth = np.degrees(np.arctan2(east, north))
        longitude, latitude, _ = geod.fwd(13., 43., azimuth, np.hypot(north, east))
        for view, look in zip(('asc', 'desc'), looks):
            station = '%s%02d' % (view, k)
            path = project / 'data/statics' / (station+'.los')
            np.savetxt(path, np.r_[0., look], fmt='%.12e')
            observations.append(path)
            directions.append(look)
            receivers.append((station, longitude, latitude))
    directions = np.array(directions)
    station_file = project / 'data/station_info/receivers.sta'
    station_file.write_text(''.join('%s %.9f %.9f\n' % row for row in receivers))
    gf_list = project / 'data/station_info/imagery.gflist'
    with gf_list.open('w') as stream:
        stream.write('# station lon lat static disp vel tsunami insar static_file disp_file '
                     'vel_file tsunami_file los_file sigma_n sigma_e sigma_u sigma_dn '
                     'sigma_de sigma_du sigma_vn sigma_ve sigma_vu sigma_tsun sigma_los '
                     'unused unused unused unused\n')
        for (station, lon, lat), path in zip(receivers, observations):
            stream.write(' '.join([station, str(lon), str(lat), '0', '0', '0', '0', '1']+
                                  ['unused']*4+[str(path)]+['1']*10+['0.001']+['1']*4)+'\n')
    if use_toml:
        run_with_toml(project, 'prepare', ['greens', 'synthetics'])
    else:
        mt3d.inversionGFs(home, name, gf_list.name, 'cloud.sources', 'halfspace.mod',
                         1, 1, .05, 0, 1, 25, 0, 1, False, None, flat=False)
    G = mt3d.getG(home, name, 'cloud.sources', 'halfspace.mod', gf_list.name,
                 False, 'imagery', False)
    assert G.shape == (48, 12)
    assert np.linalg.matrix_rank(G) == 12
    np.testing.assert_allclose(mt3d.getG(home, name, 'cloud.sources', 'halfspace.mod',
                               gf_list.name, True, 'imagery', False), G, rtol=0, atol=0)

    # Independent point-isotropic half-space reference. M = m I in N m:
    # u_NEU = (1-2nu)*m/(2*pi*mu) * [north,east,depth]/R^3.
    # Equivalent Mogi cavity dV = m/(lambda+2mu); this is not magma volume.
    offsets = np.c_[points, np.full(len(points), reference['depth_m'])]
    radii = np.linalg.norm(offsets, axis=1)
    isotropic_neu = ((1-2*reference['nu'])*10**16.6/(2*np.pi*reference['mu_Pa']) *
                     offsets/radii[:, None]**3)
    isotropic_los = np.sum(np.repeat(isotropic_neu, 2, axis=0)*directions, axis=1)
    fk_isotropic = G[:, 0]+G[:, 3]+G[:, 5]
    isotropic_error = np.linalg.norm(fk_isotropic-isotropic_los)/np.linalg.norm(isotropic_los)
    shear_neu = np.array(reference['displacement_neu_m'])
    shear_los = np.sum(np.repeat(shear_neu, 2, axis=0)*directions, axis=1)
    shear_error = np.linalg.norm(G[:, 1]-shear_los)/np.linalg.norm(shear_los)
    assert isotropic_error < 1e-3, isotropic_error
    assert shear_error < 1e-3, shear_error
    deviatoric_errors = {}
    for case in reference['deviatoric_references']:
        reference_neu = np.array(case['displacement_neu_m'])
        reference_los = np.sum(np.repeat(reference_neu, 2, axis=0)*directions, axis=1)
        prediction = G[:, :6].dot(case['tensor_components'])
        error = np.linalg.norm(prediction-reference_los)/np.linalg.norm(reference_los)
        assert error < 1e-3, (case['name'], error)
        deviatoric_errors[case['name']] = float(error)

    # Direct prescribed mechanisms take a separate syn call from the six basis
    # calls, checking superposition and the complete producer/consumer handoff.
    mechanism = np.array([.6, .2, -.1, -.4, .3, .25])
    mt3d.make_parallel_synthetics(home, name, station_file.name, 'cloud.sources',
                                'halfspace.mod', 0, 1, True, mechanism, insar=True)
    fixed = mt3d.getG(home, name, 'cloud.sources', 'halfspace.mod', gf_list.name,
                     False, 'fixed', True)
    expected = G.dot(np.tile(mechanism, 2))
    # Sum the two independently prescribed source responses.
    superposition_error = np.linalg.norm(fixed.sum(axis=1)-expected)/np.linalg.norm(expected)
    assert superposition_error < 1e-3, superposition_error
    amplitudes = np.array([1., -.4])
    for path, value, look in zip(observations, fixed.dot(amplitudes), directions):
        np.savetxt(path, np.r_[value, look], fmt='%.12e')
    if use_toml:
        run_with_toml(project, 'fixed', ['synthetics', 'invert'], mechanism)
    else:
        mt3d.run_inversion(home, name, 'fixed', 'cloud.sources', 'halfspace.mod',
                          gf_list.name, False, 'fixed', [0], len(sources),
                          'lstsq', True, mechanism, weight=True)
    fixed_output = np.loadtxt(project / 'output/inverse_models/models/fixed.0000.inv', ndmin=2)
    np.testing.assert_allclose(fixed_output[:, 4:10]/10**16.6,
                               amplitudes[:, None]*mechanism, rtol=1e-7, atol=1e-10)
    # Restore six-component tables after the fixed-mechanism check.
    mt3d.make_parallel_synthetics(home, name, station_file.name, 'cloud.sources',
                                'halfspace.mod', 0, 1, False, None, insar=True)

    models = {
        'mixed': np.array([[.6, .2, -.1, -.4, .3, .25], [1, 0, 0, 1, 0, 1]]),
        'inflation': np.array([[1, 0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0]]),
        'deflation': np.array([[-1, 0, 0, -1, 0, -1], [0, 0, 0, 0, 0, 0]])}
    errors = {}
    for run, model in models.items():
        data = G.dot(model.reshape(-1))
        for path, value, look in zip(observations, data, directions):
            np.savetxt(path, np.r_[value, look], fmt='%.12e')
        if use_toml:
            run_with_toml(project, run, ['synthetics', 'invert'] if run == 'mixed' else ['invert'],
                          matrix_mode='assemble' if run == 'mixed' else 'reuse',
                          matrix_name='toml_mixed')
        else:
            mt3d.run_inversion(home, name, run, 'cloud.sources', 'halfspace.mod',
                              gf_list.name, False, 'imagery', [0], len(sources),
                              'lstsq', False, None, weight=True)
        output = np.loadtxt(project / 'output/inverse_models/models' / (run+'.0000.inv'), ndmin=2)
        recovered = output[:, 4:10]/10**16.6
        errors[run] = float(np.linalg.norm(recovered-model)/np.linalg.norm(model))
        assert errors[run] < 1e-7, errors[run]
        np.testing.assert_array_equal(output[:, 0], sources[:, 0])
        saved = np.array([np.loadtxt(project / 'output/inverse_models/statics' /
                                    (run+'.0000.'+station+'.los'))[0]
                          for station, _, _ in receivers])
        np.testing.assert_allclose(saved, G.dot(recovered.reshape(-1)), rtol=1e-7, atol=1e-10)
        assert abs(mt3d.get_moment(home, name, 'cloud.sources', model, False)[0]-output[:, -1].sum()) < 1e7

    # Offset components use their own measurement directions. These simple
    # horizontal projections check axis/sign conventions without claiming a
    # production optical/radar offset file reader.
    offset_directions = np.array([[1, 0, 0], [0, 1, 0]])
    np.testing.assert_allclose(isotropic_neu.dot(offset_directions.T), isotropic_neu[:, :2])
    # Fixed scalar contraction has positive scalar moment and signed tensors.
    np.testing.assert_allclose(mt3d.get_moment(home, name, 'cloud.sources', [-1], True,
                                             [1, 0, 0, 1, 0, 1])[1],
                               [np.sqrt(1.5)*10**16.6])
    single_error = validate_single_receiver(project, use_toml)
    summary = {'observations': 48, 'source_ids': sources[:, 0].tolist(),
               'matrix_rank': int(np.linalg.matrix_rank(G)),
               'matrix_condition': float(np.linalg.cond(G)),
               'isotropic_reference_relative_error': float(isotropic_error),
               'okada_shear_reference_relative_error': float(shear_error),
               'deviatoric_reference_relative_errors': deviatoric_errors,
               'superposition_relative_error': float(superposition_error),
               'python': sys.version, 'numpy': np.__version__,
               'source_sha256': {module: hashlib.sha256(Path(mt3d.__file__).with_name(module+'.py').read_bytes()).hexdigest()
                                 for module in ('mt3d', 'parallel', 'inverse', 'analysis')},
               'recovery_relative_errors': errors,
               'single_receiver_recovery_relative_error': single_error,
               'toml_runner': use_toml,
               'limitations': ['same-grid noiseless recovery', 'one homogeneous material',
                               'serial execution', 'no real imagery or production pixel-offset ingestion']}
    (project / 'validation.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, required=True)
    parser.add_argument('--toml', action='store_true', help='Exercise project initialization and the TOML runner')
    args = parser.parse_args()
    main(args.project, args.toml)
