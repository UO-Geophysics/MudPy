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

from mudpy import geotomo


def run_with_toml(project, name, steps, mechanism=None, matrix_mode='assemble', matrix_name=None,
                  *, sources='cloud.sources', earth_model='halfspace.mod',
                  observations='imagery.gflist', damping=(0.,), dk=.05, kmax=25.,
                  observation_format='gflist', response_name='elastic', processes=1, keep_fk=True):
    """Write an explicit experiment configuration for this reference case."""
    source_model = 'fixed_tensor' if mechanism is not None else 'full_tensor'
    text = '''schema_version = 1
[project]
directory = ".."
[run]
name = %s
steps = %s
[inputs]
sources = %s
earth_model = %s
observations = %s
observation_format = %s
%s
[forward]
backend = "fk"
earth_flattening = false
dk = %s
pmin = 0.0
pmax = 1.0
kmax = %s
[inversion]
source_model = %s
solver = "lstsq"
weighting = "diagonal"
damping = %s
%s
[matrix]
name = %s
mode = %s
[execution]
processes = %s
keep_fk_intermediates = %s
''' % (json.dumps(name), json.dumps(steps), json.dumps(sources), json.dumps(earth_model),
       json.dumps(observations), json.dumps(observation_format),
       '[responses]\nname = '+json.dumps(response_name) if observation_format == 'hdf5' else '',
       dk, kmax, json.dumps(source_model), json.dumps(list(damping)),
       'fixed_tensor = '+json.dumps(list(mechanism)) if mechanism is not None else '',
       json.dumps(matrix_name or 'toml_'+name), json.dumps(matrix_mode), processes, str(keep_fk).lower())
    path = project/'configs'/(name+'.toml')
    with path.open('x') as stream:
        stream.write(text)
    return geotomo.run(path)


def hdf5_study(project, processes=1):
    """Small FK reference, HDF5 round trip, fixed basis, and NumPy-only reuse."""
    import time
    from unittest.mock import patch
    import h5py
    from mudpy import geotomo_io as io
    from mudpy.geotomo_config import _code_identity

    main(project, use_toml=True)
    project = project.resolve()
    syn_reference = check_syn_reference(project,
        np.loadtxt(project/'data/model_info/cloud.sources', ndmin=2),
        project/'data/station_info/receivers.sta', 'halfspace.mod')
    rows = [line.split() for line in (project/'data/station_info/imagery.gflist').read_text().splitlines()
            if line.strip() and not line.startswith('#')]
    values = np.array([np.loadtxt(row[12]) for row in rows])
    observation_file = project/'data/statics/imagery.h5'
    io.write_observations(observation_file, [r[0] for r in rows],
                          np.array([[float(r[1]), float(r[2])] for r in rows]),
                          values[:, 0], values[:, 1:], sigma_m=np.full(len(rows), .001),
                          metadata={'coordinate_system': 'WGS84 longitude/latitude degrees',
                                    'interval': 'synthetic', 'observation_kind': 'point LOS'})
    started = time.perf_counter()
    run_with_toml(project, 'hdf5_full', ['greens', 'synthetics', 'invert'],
                  observations=observation_file.name, observation_format='hdf5')
    full_time = time.perf_counter()-started
    reference = np.load(project/'GFs/matrices/imagery.npy')
    matrix_file = project/'GFs/matrices/toml_hdf5_full.npy'
    started = time.perf_counter()
    matrix = np.load(matrix_file, allow_pickle=False)
    load_time = time.perf_counter()-started
    relative = np.linalg.norm(matrix-reference)/np.linalg.norm(reference)
    assert relative < 1e-9, relative
    models = project/'output/inverse_models/models'
    truth = np.loadtxt(models/'deflation.0000.inv', ndmin=2)[:, 4:10]
    actual = np.loadtxt(models/'hdf5_full.0000.inv', ndmin=2)[:, 4:10]
    np.testing.assert_allclose(actual, truth, rtol=0, atol=1e-7*10**16.6)
    with h5py.File(str(project/'output/inverse_models/statics/hdf5_full.0000.predictions.h5'), 'r') as file:
        np.testing.assert_allclose(file['prediction_m'][()], matrix.dot((actual/10**16.6).reshape(-1)), atol=1e-11)
        np.testing.assert_allclose(file['residual_m'][()], file['los_m'][()]-file['prediction_m'][()], atol=0)

    with patch('mudpy.geotomo_io.assemble_matrix', side_effect=AssertionError('Unexpected matrix assembly')), \
            patch('mudpy.geotomo._inversion_gfs', side_effect=AssertionError('Unexpected GF generation')):
        started = time.perf_counter()
        run_with_toml(project, 'hdf5_reuse', ['invert'], matrix_mode='reuse',
                      matrix_name='toml_hdf5_full', observations=observation_file.name,
                      observation_format='hdf5', damping=(0., .01))
        reuse_time = time.perf_counter()-started
    np.testing.assert_array_equal(np.loadtxt(models/'hdf5_full.0000.inv'),
                                  np.loadtxt(models/'hdf5_reuse.0000.inv'))
    # The complete physical basis also supports prescribed mechanisms without syn.
    with patch('mudpy.geotomo._inversion_gfs', side_effect=AssertionError('Unexpected GF generation')):
        run_with_toml(project, 'hdf5_fixed', ['invert'], mechanism=[1, 0, 0, 1, 0, 1],
                      observations=observation_file.name, observation_format='hdf5')
    fixed = np.load(project/'GFs/matrices/toml_hdf5_fixed.npy')
    np.testing.assert_allclose(fixed, np.column_stack([matrix[:, 6*k:6*k+6][:, [0, 3, 5]].sum(axis=1)
                                                     for k in range(2)]), atol=1e-15)
    np.testing.assert_allclose(np.loadtxt(models/'hdf5_fixed.0000.inv', ndmin=2)[:, 4:10], truth,
                               rtol=0, atol=1e-7*10**16.6)
    assert not (project/'.geotomo.lock').exists()
    report = {'matrix_shape': list(matrix.shape), 'representative_points': len(rows),
              'source_count': 2, 'relative_matrix_difference_from_text': float(relative),
              'full_pipeline_s': full_time, 'matrix_load_s': load_time,
              'cached_two_damping_run_s': reuse_time,
              'reuse_bypasses_generation_and_assembly': True,
              'fixed_tensor_uses_existing_library': True, 'syn_reference': syn_reference,
              'source_sha256': _code_identity()['source_sha256'],
              'limitations': ['Tiny synthetic timing is not a large-library benchmark.',
                              'No real-data GF generation or inversion.']}
    if processes > 1:
        native = sorted((project/'GFs/static').glob('halfspace.mod.static.*.insar'))
        native_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in native}
        # Preserve the full-basis text outputs before regenerating in parallel.
        run_with_toml(project, 'text_serial', ['synthetics'])
        text_files = sorted((project/'GFs/static').glob('_M*.los'))
        text_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in text_files}
        run_with_toml(project, 'text_parallel', ['greens', 'synthetics'], processes=processes)
        assert text_hashes == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in text_files}
        started = time.perf_counter()
        run_with_toml(project, 'hdf5_parallel', ['greens', 'synthetics'],
                      observations=observation_file.name, observation_format='hdf5',
                      response_name='parallel_elastic', processes=processes)
        parallel_time = time.perf_counter()-started
        assert native_hashes == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in native}
        with h5py.File(str(project/'GFs/static/elastic.h5'), 'r') as serial, \
                h5py.File(str(project/'GFs/static/parallel_elastic.h5'), 'r') as parallel:
            assert serial.attrs['complete'] and parallel.attrs['complete']
            for key in ('sources', 'receiver_ids', 'receiver_coordinates', 'displacement', 'completed_sources'):
                np.testing.assert_array_equal(serial[key][()], parallel[key][()])
        parallel_matrix = io.assemble_matrix(project/'GFs/static/parallel_elastic.h5',
            io.read_observations(observation_file), np.loadtxt(project/'data/model_info/cloud.sources', ndmin=2))
        np.testing.assert_array_equal(parallel_matrix, matrix)
        report['parallel'] = {'requested_processes': processes, 'active_sources': 2,
                              'forward_s': parallel_time, 'native_files_identical': True,
                              'text_tables_identical': True,
                              'neu_arrays_and_matrix_identical': True}
    # Remove only this synthetic calculation's native products, then prove that
    # library assembly and direct matrix reuse remain independent of them.
    run_with_toml(project, 'hdf5_cleanup', ['greens', 'synthetics'],
                  observations=observation_file.name, observation_format='hdf5',
                  response_name='cleanup_elastic', processes=processes, keep_fk=False)
    sources = np.loadtxt(project/'data/model_info/cloud.sources', ndmin=2)
    for source in sources:
        assert not (project/('GFs/static/halfspace.mod.static.%.4f.sub%04d.insar' % (source[3], source[0]))).exists()
        assert not (project/('GFs/static/halfspace.mod_%.4f.sub%04d' % (source[3], source[0]))).exists()
    with patch('mudpy.geotomo._inversion_gfs', side_effect=AssertionError('Unexpected GF generation')):
        run_with_toml(project, 'cleanup_assemble', ['invert'],
                      observations=observation_file.name, observation_format='hdf5', response_name='cleanup_elastic')
        with patch('mudpy.geotomo_io.assemble_matrix', side_effect=AssertionError('Unexpected matrix assembly')):
            run_with_toml(project, 'cleanup_reuse', ['invert'], matrix_mode='reuse',
                          matrix_name='toml_cleanup_assemble', observations=observation_file.name,
                          observation_format='hdf5', response_name='cleanup_elastic')
    np.testing.assert_array_equal(matrix, np.load(project/'GFs/matrices/toml_cleanup_assemble.npy'))
    np.testing.assert_array_equal(np.loadtxt(models/'hdf5_full.0000.inv'),
                                  np.loadtxt(models/'cleanup_reuse.0000.inv'))
    native_record = json.loads((project/'GFs/static/_geotomo_greens.json').read_text())
    cleanup_record = json.loads((project/'GFs/static/cleanup_elastic.h5.json').read_text())
    assert native_record['availability'] == 'pruned'
    report['cleanup'] = {'receipt': cleanup_record['native_cleanup'],
                         'assembly_and_cached_reuse_verified': True,
                         'matrix_and_model_identical': True}
    (project/'hdf5_validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


def check_syn_reference(project, sources, station_file, model_name):
    """Compare float64 tensor conversion with rounded C syn on sampled receivers."""
    import subprocess
    from mudpy.green import src2sta, rt2ne

    mechanisms = np.vstack((np.eye(6), [1, 0, 0, 1, 0, 1], [-1, 0, 0, -1, 0, -1],
                            [1, .2, 0, -1, 0, 0], [-.3, .4, -.5, .8, -.9, .2], np.zeros(6)))
    errors, calls = [], 0
    for source in sources:
        _, azimuths = src2sta(str(station_file), source)
        native = project/('GFs/static/%s.static.%.4f.sub%04d.insar' % (model_name, source[3], source[0]))
        statics = np.loadtxt(native, ndmin=2)
        selected = np.linspace(0, len(statics)-1, min(12, len(statics)), dtype=int)
        vector = geotomo._tensor_static_neu(statics[selected], azimuths[selected], mechanisms)
        reference = np.empty_like(vector)
        for index, receiver in enumerate(selected):
            payload = ' '.join('%.10e' % value for value in statics[receiver])
            for component, mechanism in enumerate(mechanisms):
                command = ['syn', '-M'+str(10**23.6)+'/'+'/'.join(map(str, mechanism)),
                           '-A'+str(azimuths[receiver]), '-P']
                result = subprocess.run(command, input=payload, text=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                values = np.fromstring(result.stdout, sep=' ')
                if values.shape != (5,) or not np.isfinite(values).all():
                    raise ValueError('Invalid syn reference output: '+result.stdout+' '+result.stderr)
                north, east = rt2ne(values[3]/100, values[4]/100, azimuths[receiver])
                reference[index, :, component] = [north, east, values[2]/100]
                calls += 1
        for component in range(len(mechanisms)):
            scale = np.linalg.norm(reference[:, :, component])
            difference = np.linalg.norm(vector[:, :, component]-reference[:, :, component])
            if not scale:
                assert difference == 0
            else:
                errors.append(float(difference/scale))
    maximum = max(errors)
    assert maximum < 6e-4, maximum
    return {'calls': calls, 'max_mechanism_relative_l2': maximum, 'tolerance': 6e-4,
            'mechanisms': mechanisms.tolist(),
            'precision': 'C syn uses float32 and four-significant-digit output; NumPy retains float64'}


def profile_forward(project):
    """Time source-parallel FK and vectorized synthesis on a fixed synthetic case."""
    from collections import defaultdict
    from threading import Lock
    import resource
    import subprocess
    import time
    from unittest.mock import patch
    import h5py
    from mudpy import geotomo_io as io
    from mudpy.geotomo_config import _code_identity

    geotomo.init_project(project)
    project = project.resolve()
    sources = np.array([[9, 0, 0, 1.01], [2, .02, -.01, 7.], [5, -.01, .02, 13.]])
    np.savetxt(project/'data/model_info/profile.sources', sources,
               fmt=['%d', '%.17g', '%.17g', '%.4f'])
    (project/'structure/profile.mod').write_text(
        '2 2.6 4.5 2.4 300 600\n6 3.3 5.8 2.7 300 600\n0 3.6 6.3 2.8 300 600\n')
    count = 1507
    angles = np.arange(count)*np.pi*(3-np.sqrt(5))
    radii = .01+.3*np.sqrt((np.arange(count)+.5)/count)
    xy = np.column_stack((radii*np.cos(angles), radii*np.sin(angles)))
    ids = ['P%04d' % k for k in range(count)]
    io.write_observations(project/'data/statics/profile.h5', ids, xy,
        np.zeros(count), np.tile([.1, -.4, .8], (count, 1)), sigma_m=np.full(count, .001),
        metadata={'coordinate_system': 'Synthetic WGS84 longitude/latitude degrees',
                  'interval': 'synthetic profiling only'})
    station_file = project/'data/station_info/profile.sta'
    with station_file.open('x') as stream:
        for identifier, coordinate in zip(ids, xy):
            stream.write('%s %.17g %.17g\n' % (identifier, *coordinate))
    real_run, timings, lock = subprocess.run, {}, Lock()

    def measure(label, operation):
        calls, active, maximum = (defaultdict(int) for _ in range(3))
        def observed_run(command, *args, **kwargs):
            kind = Path(command[0]).name
            if kind not in ('syn', 'fk.pl'):
                return real_run(command, *args, **kwargs)
            if kind == 'syn':
                raise AssertionError('Production synthesis unexpectedly launched syn')
            with lock:
                calls[kind] += 1
                active[kind] += 1
                maximum[kind] = max(maximum[kind], active[kind])
            try:
                return real_run(command, *args, **kwargs)
            finally:
                with lock:
                    active[kind] -= 1
        def usage():
            own = resource.getrusage(resource.RUSAGE_SELF)
            child = resource.getrusage(resource.RUSAGE_CHILDREN)
            return np.array([own.ru_utime, own.ru_stime, child.ru_utime, child.ru_stime])
        print('Profiling '+label+' ...', flush=True)
        before, started = usage(), time.perf_counter()
        with patch('subprocess.run', new=observed_run):
            operation()
        elapsed, cpu = time.perf_counter()-started, usage()-before
        timings[label] = {'wall_s': elapsed, 'parent_user_s': cpu[0], 'parent_system_s': cpu[1],
                         'children_user_s': cpu[2], 'children_system_s': cpu[3],
                         'average_cpu_cores': float(cpu.sum()/elapsed), 'calls': dict(calls),
                         'maximum_concurrent_calls': dict(maximum)}
        print(json.dumps({label: timings[label]}, indent=2), flush=True)

    native_hashes = None
    for workers in (1, 3):
        def run_stage(stage):
            label = '%s_p%d' % (stage, workers)
            return run_with_toml(project, label, [stage], sources='profile.sources',
                earth_model='profile.mod', observations='profile.h5', observation_format='hdf5',
                response_name='responses_p%d' % workers, processes=workers)
        measure('fk_p%d' % workers, lambda: run_stage('greens'))
        current = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in (project/'GFs/static').glob('profile.mod.static.*.insar')}
        assert len(current) == 3
        if native_hashes is not None:
            assert native_hashes == current
        native_hashes = current
        measure('synthetics_p%d' % workers, lambda: run_stage('synthetics'))

    with h5py.File(str(project/'GFs/static/responses_p1.h5'), 'r') as serial, \
            h5py.File(str(project/'GFs/static/responses_p3.h5'), 'r') as parallel:
        assert serial.attrs['complete'] and parallel.attrs['complete']
        expected = serial['displacement'][()]
        assert np.isfinite(expected).all()
        np.testing.assert_array_equal(expected, parallel['displacement'][()])
        np.testing.assert_array_equal(parallel['sources'][()], sources)
        metadata = json.loads(serial.attrs['metadata_json'])
    observations = io.read_observations(project/'data/statics/profile.h5')
    serial_matrix = io.assemble_matrix(project/'GFs/static/responses_p1.h5', observations, sources)
    np.testing.assert_array_equal(serial_matrix,
        io.assemble_matrix(project/'GFs/static/responses_p3.h5', observations, sources))
    def write_library():
        with io.response_writer(project/'GFs/static/write_timing.h5', sources, ids, xy,
                                dict(metadata, profiling_only=True)) as (data, completed):
            for index in range(len(sources)):
                data[index] = expected[index]
                completed[index] = True
    measure('hdf5_write', write_library)
    reference = check_syn_reference(project, sources, station_file, 'profile.mod')
    report = {'sources': len(sources), 'receivers': count, 'syn_calls_per_conversion': 0,
              'synthesis': 'numpy_mt_radiat_float64', 'native_files_identical': True,
              'response_arrays_and_matrices_identical': True, 'syn_reference': reference,
              'timings': timings, 'source_sha256': _code_identity()['source_sha256'],
              'validation_driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'limitations': ['Single measurement per configuration; not a scaling law.',
                              'System load is uncontrolled.',
                              'Synthetic observations only; no Norcia changes or inversion.']}
    with (project/'forward_profile.json').open('x') as stream:
        stream.write(json.dumps(report, indent=2)+'\n')
    print('Profiling complete: '+str(project/'forward_profile.json'), flush=True)


def validate_single_receiver(project, use_toml=False):
    """Exercise real one-source/one-observation files and project exclusion."""
    single = project / 'single_receiver'
    if use_toml:
        geotomo.init_project(single)
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
        geotomo.inversionGFs(home, name, 'imagery.gflist', 'cloud.sources', 'halfspace.mod',
                         1, 1, .05, 0, 1, 25, 0, 1, True, mechanism, flat=False)
    G = geotomo.getG(home, name, 'cloud.sources', 'halfspace.mod', 'imagery.gflist', False, 'single', True)
    assert G.shape == (1, 1)
    np.savetxt(observation, np.r_[-.5*G[0, 0], look])
    arguments = (home, name, 'contraction', 'cloud.sources', 'halfspace.mod',
                 'imagery.gflist', False, 'single', [0], 1, 'lstsq', True, mechanism)
    if use_toml:
        run_with_toml(single, 'contraction', ['invert'], mechanism)
    else:
        geotomo.run_inversion(*arguments, weight=True)
    result = np.loadtxt(single / 'output/inverse_models/models/contraction.0000.inv', ndmin=2)
    np.testing.assert_allclose(result[:, 4:10]/10**16.6, [-.5*np.array(mechanism)], atol=1e-10)
    with geotomo._project_lock(home, name):
        try:
            geotomo.run_inversion(*arguments)
        except FileExistsError:
            pass  # Expected exclusion before matrix/output mutation.
        else:
            raise AssertionError('Concurrent project ownership was not rejected')
    assert Path.cwd() == cwd
    assert not (single / '.geotomo.lock').exists()
    return float(abs(result[0, 4]/10**16.6+.5)/.5)


def main(project, use_toml=False):
    reference = json.loads(Path(__file__).with_name('okada_reference.json').read_text())
    project = project.resolve()
    if use_toml:
        geotomo.init_project(project)
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
        geotomo.inversionGFs(home, name, gf_list.name, 'cloud.sources', 'halfspace.mod',
                         1, 1, .05, 0, 1, 25, 0, 1, False, None, flat=False)
    G = geotomo.getG(home, name, 'cloud.sources', 'halfspace.mod', gf_list.name,
                 False, 'imagery', False)
    assert G.shape == (48, 12)
    assert np.linalg.matrix_rank(G) == 12
    np.testing.assert_allclose(geotomo.getG(home, name, 'cloud.sources', 'halfspace.mod',
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

    # Direct prescribed mechanisms are evaluated separately from the six basis
    # calls, checking superposition and the complete producer/consumer handoff.
    mechanism = np.array([.6, .2, -.1, -.4, .3, .25])
    geotomo.make_parallel_synthetics(home, name, station_file.name, 'cloud.sources',
                                'halfspace.mod', 0, 1, True, mechanism, insar=True)
    fixed = geotomo.getG(home, name, 'cloud.sources', 'halfspace.mod', gf_list.name,
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
        geotomo.run_inversion(home, name, 'fixed', 'cloud.sources', 'halfspace.mod',
                          gf_list.name, False, 'fixed', [0], len(sources),
                          'lstsq', True, mechanism, weight=True)
    fixed_output = np.loadtxt(project / 'output/inverse_models/models/fixed.0000.inv', ndmin=2)
    np.testing.assert_allclose(fixed_output[:, 4:10]/10**16.6,
                               amplitudes[:, None]*mechanism, rtol=1e-7, atol=1e-10)
    # Restore six-component tables after the fixed-mechanism check.
    geotomo.make_parallel_synthetics(home, name, station_file.name, 'cloud.sources',
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
            geotomo.run_inversion(home, name, run, 'cloud.sources', 'halfspace.mod',
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
        assert abs(geotomo.get_moment(home, name, 'cloud.sources', model, False)[0]-output[:, -1].sum()) < 1e7

    # Offset components use their own measurement directions. These simple
    # horizontal projections check axis/sign conventions without claiming a
    # production optical/radar offset file reader.
    offset_directions = np.array([[1, 0, 0], [0, 1, 0]])
    np.testing.assert_allclose(isotropic_neu.dot(offset_directions.T), isotropic_neu[:, :2])
    # Fixed scalar contraction has positive scalar moment and signed tensors.
    np.testing.assert_allclose(geotomo.get_moment(home, name, 'cloud.sources', [-1], True,
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
               'source_sha256': {module: hashlib.sha256(Path(geotomo.__file__).with_name(module+'.py').read_bytes()).hexdigest()
                                 for module in ('geotomo', 'parallel', 'inverse', 'analysis')},
               'recovery_relative_errors': errors,
               'single_receiver_recovery_relative_error': single_error,
               'toml_runner': use_toml,
               'limitations': ['same-grid noiseless recovery', 'one homogeneous material',
                               'serial execution', 'no real imagery or production pixel-offset ingestion']}
    (project / 'validation.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


def tensor_metrics(locations, tensors, true_locations, true_tensors):
    """Aggregate tensor and scalar-moment-weighted spatial diagnostics.

    Locations are NE/depth in km; tensors scale 10**16.6 N m. Centroid/spread
    are undefined for a zero model. Aggregate agreement does not establish
    correct mechanisms at individual locations.
    """
    locations, tensors = np.asarray(locations), np.asarray(tensors)
    true_locations, true_tensors = np.asarray(true_locations), np.asarray(true_tensors)
    if (locations.shape != (len(tensors), 3) or tensors.ndim != 2 or tensors.shape[1] != 6 or
            true_locations.shape != (len(true_tensors), 3) or true_tensors.ndim != 2 or
            true_tensors.shape[1] != 6 or
            not all(np.isfinite(x).all() for x in (locations, tensors, true_locations, true_tensors))):
        raise ValueError('Expected finite location/tensor arrays of shapes (N,3)/(N,6)')

    def norm(values):
        return np.sqrt(np.sum(values*values*np.array([1, 2, 2, 1, 2, 1]), axis=-1))

    true_weights, weights = norm(true_tensors), norm(tensors)
    if not true_weights.sum() or not norm(true_tensors.sum(axis=0)):
        raise ValueError('Truth must have nonzero total scalar moment and aggregate tensor')
    true_centroid = np.average(true_locations, axis=0, weights=true_weights)
    true_spread = np.sqrt(np.average(np.sum((true_locations-true_centroid)**2, axis=1), weights=true_weights))
    summed, true_sum = tensors.sum(axis=0), true_tensors.sum(axis=0)
    summed_norm, true_norm = norm(summed), norm(true_sum)
    centroid = np.average(locations, axis=0, weights=weights) if weights.sum() else None
    spread = np.sqrt(np.average(np.sum((locations-centroid)**2, axis=1), weights=weights)) if centroid is not None else None
    cosine = np.sum(summed*true_sum*[1, 2, 2, 1, 2, 1])/(summed_norm*true_norm) if summed_norm else None
    return {
        'centroid_km': centroid.tolist() if centroid is not None else None,
        'true_centroid_km': true_centroid.tolist(),
        'centroid_error_km': float(np.linalg.norm(centroid-true_centroid)) if centroid is not None else None,
        'spread_km': float(spread) if spread is not None else None,
        'true_spread_km': float(true_spread),
        'aggregate_tensor_relative_error': float(norm(summed-true_sum)/true_norm),
        'aggregate_mechanism_angle_deg': float(np.degrees(np.arccos(np.clip(cosine, -1, 1)))) if cosine is not None else None,
        'total_scalar_moment_Nm': float(weights.sum()/np.sqrt(2)*10**16.6),
        'true_total_scalar_moment_Nm': float(true_weights.sum()/np.sqrt(2)*10**16.6),
        'scalar_moment_ratio': float(weights.sum()/true_weights.sum()),
        'cancellation_ratio': float(weights.sum()/summed_norm) if summed_norm else None,
        'aggregate_tensor_Nm': (summed*10**16.6).tolist()}


def plot_layered_study(project, summary, matrices, receiver_xy, grid_xyz, truth_xyz, decoded, predictions, mixed_data):
    """Export scientific PNG/PDF diagnostics; keep Matplotlib caches temporary."""
    import os
    import tempfile
    previous = os.environ.get('MPLCONFIGDIR')
    with tempfile.TemporaryDirectory(prefix='mudpy-matplotlib-') as cache:
        os.environ['MPLCONFIGDIR'] = cache
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D  # Registers the 3D projection.

            def save(fig, name):
                for extension in ('png', 'pdf'):
                    fig.savefig(project/'plots'/(name+'.'+extension), dpi=180, bbox_inches='tight')
                plt.close(fig)

            fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
            for ax, kind, values in zip(axes, ('dk', 'kmax'), ((.2, .1, .05, .025), (15, 25, 35))):
                names = ['layered_dk_'+str(v).replace('.', 'p') if kind == 'dk' else
                         ('layered_dk_0p025' if v == 35 else 'layered_kmax_'+str(v)) for v in values]
                errors = [summary['convergence'][name]['max_column_relative_l2'] for name in names]
                ax.plot(values, np.array(errors)*100, 'o-')
                ax.axhline(summary['forward_acceptance_bound']*100, color='0.5', linestyle='--', label='Initial acceptance bound')
                ax.set(xlabel='FK '+kind, ylabel='Largest column error (%)', title='Relative to dk=0.025, kmax=35')
                if kind == 'dk':
                    ax.set_xscale('log')
                    ax.set_xticks(values)
                    ax.set_xticklabels([str(v) for v in values])
                    ax.minorticks_off()
                ax.grid(alpha=.2)
            axes[0].legend(fontsize=8)
            fig.suptitle('Two-layer response convergence (float64 synthesis)')
            save(fig, 'layered_convergence')

            iso = np.array([1, 0, 0, 1, 0, 1])*np.sqrt(2/3)
            homogeneous = matrices['homogeneous'][:, :6].dot(iso)[::2]*1000
            layered = matrices['layered_dk_0p025'][:, :6].dot(iso)[::2]*1000
            fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
            maximum = max(abs(homogeneous).max(), abs(layered).max())
            for ax, values, title in zip(axes, (homogeneous, layered, layered-homogeneous),
                                         ('Homogeneous', 'Two layers', 'Layered minus homogeneous')):
                limit = maximum if ax is not axes[2] else abs(values).max()
                artist = ax.scatter(receiver_xy[:, 1], receiver_xy[:, 0], c=values, marker='s', s=65,
                                    cmap='RdBu_r', vmin=-limit, vmax=limit)
                ax.scatter(truth_xyz[0, 1], truth_xyz[0, 0], marker='*', c='black', s=90)
                ax.set(xlabel='East (km)', ylabel='North (km)', title=title, aspect='equal')
                fig.colorbar(artist, ax=ax, label='Ascending LOS (mm)', shrink=.8)
            fig.suptitle('Known off-grid isotropic source; same normalization and geometry')
            save(fig, 'layered_response')

            fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
            measures = (('clean_prediction_relative_error', 'Clean displacement error', True),
                        ('centroid_error_km', 'Moment-weighted centroid error (km)', False),
                        ('aggregate_tensor_relative_error', 'Aggregate tensor error', True),
                        ('scalar_moment_ratio', 'Total scalar moment / truth', True))
            colors = dict(zip(('shear', 'inflation', 'deflation', 'mixed'), ('C0', 'C1', 'C2', 'C3')))
            for ax, (measure, title, logarithmic) in zip(axes.flat, measures):
                for name, records in summary['recovery'].items():
                    if name == 'on_grid_control':
                        continue
                    mechanism = name.split('_')[0]
                    ax.plot([r['damping'] for r in records], [r[measure] for r in records],
                            color=colors[mechanism], linestyle='-' if name.endswith('_noisy') else '--',
                            marker='o', markersize=3, label=name.replace('_', ' '))
                ax.set_xscale('symlog', linthresh=.01)
                if logarithmic:
                    ax.set_yscale('log')
                ax.set(xlabel='Identity damping', ylabel=title)
                ax.grid(alpha=.2)
            axes[1, 1].axhline(1, color='black', linewidth=.8)
            handles, labels = axes[0, 0].get_legend_handles_labels()
            fig.legend(handles, labels, loc='lower center',
                       ncol=4, fontsize=8, bbox_to_anchor=(.5, .005))
            fig.suptitle('Off-grid imagery recovery: fit, location, integrated tensor and cancellation')
            fig.tight_layout(rect=(0, .10, 1, .95))
            save(fig, 'offgrid_tradeoffs')

            index = next(k for k, item in enumerate(summary['recovery']['mixed_noisy']) if item['damping'] == summary['display_damping'])
            prediction = predictions['mixed_noisy', index]
            fig = plt.figure(figsize=(10, 8), constrained_layout=True)
            maximum = max(abs(mixed_data).max(), abs(prediction).max())*1000
            for k, (values, title) in enumerate(((mixed_data, 'Synthetic observations'), (prediction, 'Predictions'),
                                               (mixed_data-prediction, 'Observed minus predicted'))):
                ax = fig.add_subplot(2, 2, k+1)
                limit = maximum if k < 2 else abs(values).max()*1000
                artist = ax.scatter(receiver_xy[:, 1], receiver_xy[:, 0], c=values[::2]*1000,
                                    cmap='RdBu_r', marker='s', s=65, vmin=-limit, vmax=limit)
                ax.scatter(truth_xyz[:, 1], truth_xyz[:, 0], c='black', marker='*', s=90)
                ax.set(xlabel='East (km)', ylabel='North (km)', title=title, aspect='equal')
                fig.colorbar(artist, ax=ax, label='Ascending LOS (mm)', shrink=.8)
            ax = fig.add_subplot(2, 2, 4, projection='3d')
            tensors = decoded['mixed_noisy', index]
            norms = np.sqrt(np.sum(tensors*tensors*[1, 2, 2, 1, 2, 1], axis=1))
            active = norms > 0
            fractions = tensors[active][:, [0, 3, 5]].sum(axis=1)/(np.sqrt(3)*norms[active])
            artist = ax.scatter(grid_xyz[active, 1], grid_xyz[active, 0], grid_xyz[active, 2],
                                s=300*norms[active]/norms.max(), c=fractions, cmap='RdBu_r', vmin=-1, vmax=1)
            ax.scatter(truth_xyz[:, 1], truth_xyz[:, 0], truth_xyz[:, 2], c='black', marker='*', s=110, label='True locations')
            ax.invert_zaxis()
            ax.set(xlabel='East (km)', ylabel='North (km)', zlabel='Depth (km)', title='Inferred source moments')
            ax.legend(fontsize=8)
            fig.colorbar(artist, ax=ax, label='Signed isotropic fraction', shrink=.65, pad=.15)
            fig.suptitle('Mixed sources, 0.3 mm iid noise, preselected damping = 1\nSource-marker area scales with inferred scalar moment')
            save(fig, 'offgrid_example')
        finally:
            if previous is None:
                os.environ.pop('MPLCONFIGDIR', None)
            else:
                os.environ['MPLCONFIGDIR'] = previous


def layered_study(project):
    """Layered convergence and off-grid/noisy recovery; one project, many runs."""
    import cProfile
    import time
    project = project.resolve()
    geotomo.init_project(project)
    analysis = project/'analysis'
    driver_bytes = Path(__file__).read_bytes()
    (analysis/'validate.executed.py').write_bytes(driver_bytes)
    geod = Geod(ellps='WGS84')
    sigma = .0003  # meters; synthetic independent noise, not an InSAR error model.
    damping = np.array([0., .01, .1, 1., 10.])
    forward_tolerance = 1e-3
    origin = (13., 43.)

    def source_rows(xyz, ids):
        azimuth = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
        lon, lat, _ = geod.fwd(np.full(len(xyz), origin[0]), np.full(len(xyz), origin[1]),
                             azimuth, np.hypot(xyz[:, 0], xyz[:, 1])*1000)
        return np.c_[ids, lon, lat, xyz[:, 2]]

    truth_xyz = np.array([[-.4, .35, 3.6], [.55, -.25, 4.4]])
    grid_xyz = np.array([[n, e, z] for z in (3., 5.) for n in (-1., 1.) for e in (-1., 1.)])
    truth = source_rows(truth_xyz, [23, 11])
    grid = source_rows(grid_xyz, np.arange(8, 0, -1))
    for filename, rows in (('truth.sources', truth), ('grid.sources', grid)):
        np.savetxt(project/'data/model_info'/filename, rows, fmt=['%d', '%.9f', '%.9f', '%.9f'])
    homogeneous = [0., 6/np.sqrt(3), 6., 2.7, 500., 1000.]
    for name, model in (
            ('homogeneous.mod', [homogeneous]),
            ('identical_layers.mod', [[2., *homogeneous[1:]], homogeneous]),
            ('layered.mod', [[2., 2.4, 4.3, 2.4, 500., 1000.], homogeneous])):
        np.savetxt(project/'structure'/name, model, fmt='%.12e')
    looks = np.array([[.1, -.6, np.sqrt(.63)], [.1, .6, np.sqrt(.63)]])
    receiver_xy = np.array([[n, e] for n in np.linspace(-6, 6, 7) for e in np.linspace(-6, 6, 7)])
    receivers, directions = [], []
    for k, (north, east) in enumerate(receiver_xy):
        lon, lat, _ = geod.fwd(*origin, np.degrees(np.arctan2(east, north)), np.hypot(north, east)*1000)
        for view, look in zip(('asc', 'desc'), looks):
            receivers.append(('%s%02d' % (view, k), lon, lat))
            directions.append(look)
    directions = np.array(directions)
    station_file = project/'data/station_info/receivers.sta'
    station_file.write_text(''.join('%s %.9f %.9f unused\n' % row for row in receivers))
    with (project/'data/station_info/imagery.gflist').open('x') as stream:
        stream.write('# Standard MudPy imagery list; synthetic 0.3 mm LOS uncertainties.\n')
        for station, lon, lat in receivers:
            stream.write(' '.join([station, str(lon), str(lat), '0', '0', '0', '0', '1']+
                                  ['unused']*4+['data/statics/'+station+'.los']+
                                  ['1']*10+[str(sigma)]+['1']*4)+'\n')

    def data_files(values, label=None):
        values = np.asarray(values)
        if values.shape != (len(receivers),) or not np.isfinite(values).all():
            raise ValueError('Invalid synthetic observation vector')
        if label is not None:
            with (analysis/(label+'.data.npy')).open('xb') as stream:
                np.save(stream, values)
        for (station, _, _), value, look in zip(receivers, values, directions):
            np.savetxt(project/'data/statics'/(station+'.los'), np.r_[value, look], fmt='%.12e')

    data_files(np.full(len(receivers), .001))  # Placeholder values during forward-only runs.
    timings = {}

    def forward(name, model='layered.mod', source='truth.sources', dk=.025, kmax=35., mechanism=None, profile=False):
        print('Forward:', name, flush=True)
        profiler = cProfile.Profile() if profile else None
        started = time.perf_counter()
        if profiler is not None:
            profiler.enable()
        run_with_toml(project, name, ['greens', 'synthetics'] if mechanism is None else ['synthetics'],
                      mechanism, sources=source, earth_model=model, dk=dk, kmax=kmax)
        if profiler is not None:
            profiler.disable()
            profiler.dump_stats(str(analysis/(name+'.prof')))
        compute_time = time.perf_counter()-started
        started = time.perf_counter()
        with geotomo._project_lock(str(project.parent)+'/', project.name):
            G = geotomo.makeG(str(project.parent)+'/', project.name, source, model,
                          station_file.name, 'insar', mechanism is not None)
        timings[name] = {'generation_and_io_s': compute_time, 'assembly_s': time.perf_counter()-started}
        if source == 'truth.sources':
            np.save(analysis/(name+'.G.npy'), G)
        return G

    matrices = {}
    matrices['homogeneous'] = forward('homogeneous', 'homogeneous.mod', dk=.05)
    matrices['identical_layers'] = forward('identical_layers', 'identical_layers.mod', dk=.05)
    for dk in (.2, .1, .05, .025):
        name = 'layered_dk_'+str(dk).replace('.', 'p')
        matrices[name] = forward(name, dk=dk)
    reference = matrices['layered_dk_0p025']
    for kmax in (15., 25.):
        name = 'layered_kmax_'+str(int(kmax))
        matrices[name] = forward(name, kmax=kmax)

    def difference(a, b):
        return {'relative_frobenius': float(np.linalg.norm(a-b)/np.linalg.norm(b)),
                'max_column_relative_l2': float(np.max(np.linalg.norm(a-b, axis=0)/np.linalg.norm(b, axis=0)))}

    equivalence = difference(matrices['identical_layers'], matrices['homogeneous'])
    convergence = {name: difference(value, reference) for name, value in matrices.items() if name.startswith('layered_')}
    assert equivalence['max_column_relative_l2'] < forward_tolerance, equivalence
    for name in ('layered_dk_0p05', 'layered_kmax_25'):
        assert convergence[name]['max_column_relative_l2'] < forward_tolerance, (name, convergence[name])
    contrast = difference(reference, matrices['homogeneous'])
    assert contrast['relative_frobenius'] > .01, contrast
    # Restore the finest GF cache before testing directly prescribed mechanisms.
    # The final cutoff run changed the shared GF products; use a named explicit run.
    restored = forward('layered_reference')
    np.testing.assert_array_equal(restored, reference)
    superposition = {}
    for name, mechanism in (
            ('layered_inflation', [1, 0, 0, 1, 0, 1]),
            ('layered_deflation', [-1, 0, 0, -1, 0, -1]),
            ('layered_shear', [0, 1, 0, 0, 0, 0]),
            ('layered_mixed', [.6, .2, -.1, -.4, .3, .25])):
        actual = forward(name, mechanism=mechanism)
        expected = np.column_stack([reference[:, 6*k:6*k+6].dot(mechanism) for k in range(2)])
        error = float(np.linalg.norm(actual-expected)/np.linalg.norm(expected))
        assert error < forward_tolerance, (name, error)
        superposition[name] = error
    np.testing.assert_allclose(np.load(analysis/'layered_inflation.G.npy'),
                               -np.load(analysis/'layered_deflation.G.npy'), rtol=0, atol=0)

    G = forward('grid_prepare', source='grid.sources', profile=True)
    singular = np.linalg.svd(G, compute_uv=False)
    nearest = np.min(np.linalg.norm(truth_xyz[:, None]-grid_xyz[None, :], axis=2), axis=1)
    assert (nearest > .1).all()
    iso = np.array([1, 0, 0, 1, 0, 1])*np.sqrt(2/3)  # Scalar moment = reference moment.
    shear = np.array([0, 1, 0, 0, 0, 0])
    zero = np.zeros(6)
    models = {'shear': np.array([shear, zero]), 'inflation': np.array([iso, zero]),
              'deflation': np.array([-iso, zero]), 'mixed': np.array([shear, .6*iso])}
    rng = np.random.default_rng(20261106)
    noise = rng.normal(0, sigma, len(receivers))
    np.savez(analysis/'synthetic_truth.npz', truth_locations_km=truth_xyz, grid_locations_km=grid_xyz,
             receiver_ne_km=receiver_xy, directions_neu=directions, noise_m=noise, **models)
    results, decoded, predictions = {}, {}, {}

    def invert(name, observed, clean, model, locations, values=damping, first=False, profile=False):
        print('Invert:', name, flush=True)
        data_files(observed, name)
        profiler = cProfile.Profile() if profile else None
        started = time.perf_counter()
        if profiler is not None:
            profiler.enable()
        run_with_toml(project, name, ['invert'], matrix_mode='assemble' if first else 'reuse',
                      matrix_name='recovery_grid', sources='grid.sources', earth_model='layered.mod',
                      damping=values, dk=.025, kmax=35.)
        if profiler is not None:
            profiler.disable()
            profiler.dump_stats(str(analysis/(name+'.prof')))
        timings[name] = {'inversion_and_io_s': time.perf_counter()-started}
        output = []
        for k, penalty in enumerate(values):
            saved = np.loadtxt(project/'output/inverse_models/models'/(name+'.%04d.inv' % k), ndmin=2)
            np.testing.assert_array_equal(saved[:, 0], grid[:, 0])
            coefficients = saved[:, 4:10]/10**16.6
            prediction = np.array([np.loadtxt(project/'output/inverse_models/statics'/
                                            (name+'.%04d.' % k+station+'.los'))[0] for station, _, _ in receivers])
            recomposed = G.dot(coefficients.reshape(-1))
            assert np.sqrt(np.mean((prediction-recomposed)**2))/sigma < 1e-4
            residual = (recomposed-observed)/sigma
            gradient = (G/sigma).T.dot(residual)+penalty**2*coefficients.reshape(-1)
            scale = np.linalg.norm((G/sigma).T.dot(observed/sigma))+penalty**2*np.linalg.norm(coefficients)
            stationarity = float(np.linalg.norm(gradient)/scale)
            assert stationarity < 1e-6, (name, penalty, stationarity)
            metrics = tensor_metrics(grid_xyz, coefficients, locations, model)
            metrics.update({'damping': float(penalty), 'weighted_rms': float(np.sqrt(np.mean(residual**2))),
                            'data_rms_mm': float(1000*np.sqrt(np.mean((recomposed-observed)**2))),
                            'clean_prediction_relative_error': float(np.linalg.norm(recomposed-clean)/np.linalg.norm(clean)),
                            'coefficient_l2': float(np.linalg.norm(coefficients)), 'stationarity': stationarity})
            output.append(metrics)
            decoded[name, k], predictions[name, k] = coefficients, prediction
        # Quadratic damping trades increasing data residual for decreasing model norm.
        assert (np.diff([item['data_rms_mm'] for item in output]) >= -1e-7).all()
        norms = np.array([item['coefficient_l2'] for item in output])
        assert (np.diff(norms) <= 1e-7*max(norms)).all()
        results[name] = output

    control = np.zeros((8, 6))
    control[0] = iso
    control_data = G.dot(control.reshape(-1))
    invert('on_grid_control', control_data, control_data, control, grid_xyz, values=(0.,), first=True)
    assert results['on_grid_control'][0]['clean_prediction_relative_error'] < 1e-7
    for name, model in models.items():
        clean = reference.dot(model.reshape(-1))
        for noisy in (False, True):
            label = name+('_noisy' if noisy else '_noiseless')
            invert(label, clean+noise if noisy else clean, clean, model, truth_xyz,
                   profile=label == 'mixed_noisy')
    for k in range(len(damping)):
        np.testing.assert_allclose(decoded['inflation_noiseless', k],
                                   -decoded['deflation_noiseless', k], rtol=1e-8, atol=1e-8)
    summary = {
        'forward_acceptance_bound': forward_tolerance, 'identical_layer_limit': equivalence,
        'layer_contrast_effect': contrast, 'convergence': convergence, 'superposition': superposition,
        'source_grid_shape': [8, 6], 'observation_count': len(receivers),
        'matrix_rank': int(np.linalg.matrix_rank(G)), 'matrix_condition': float(singular[0]/singular[-1]),
        'truth_to_nearest_candidate_km': nearest.tolist(), 'noise_sigma_m': sigma,
        'noise_seed': 20261106, 'noise_realization_rms_m': float(np.sqrt(np.mean(noise**2))),
        'penalty': 'identity on raw six-component integrated moment coefficients',
        'display_damping': 1., 'recovery': results, 'timings': timings,
        'driver_sha256': hashlib.sha256(driver_bytes).hexdigest(),
        'limitations': ['Layered checks establish homogeneous limit, convergence, and superposition, not an independent layered benchmark.',
                        'Off-grid truth and inversion share the same FK backend and Earth model.',
                        'Single iid noise realization; no atmospheric correlation or nuisance ramps.',
                        'Aggregate mechanism/centroid agreement does not establish individual source recovery.',
                        'No automatic choice of damping, real imagery, or numerical topography.']}
    (analysis/'layered_offgrid.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    plot_layered_study(project, summary, matrices, receiver_xy, grid_xyz, truth_xyz, decoded,
                       predictions, reference.dot(models['mixed'].reshape(-1))+noise)
    print(json.dumps({key: summary[key] for key in ('identical_layer_limit', 'layer_contrast_effect',
                     'matrix_rank', 'matrix_condition', 'truth_to_nearest_candidate_km')}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, required=True)
    parser.add_argument('--toml', action='store_true', help='Exercise project initialization and the TOML runner')
    parser.add_argument('--layered', action='store_true', help='Run layered convergence and off-grid/noisy imagery recovery through TOML')
    parser.add_argument('--hdf5', action='store_true', help='Check HDF5 responses and direct NumPy reuse on the small synthetic reference')
    parser.add_argument('--profile-forward', action='store_true',
                        help='Profile three synthetic sources and 1507 receivers with one/three workers')
    parser.add_argument('--processes', type=int, default=1,
                        help='With --hdf5, also compare serial and parallel forward calculations')
    args = parser.parse_args()
    if args.processes < 1 or (args.processes != 1 and not args.hdf5):
        parser.error('--processes must be positive and requires --hdf5 when not 1')
    if args.profile_forward:
        profile_forward(args.project)
    elif args.hdf5:
        hdf5_study(args.project, args.processes)
    elif args.layered:
        layered_study(args.project)
    else:
        main(args.project, args.toml)
