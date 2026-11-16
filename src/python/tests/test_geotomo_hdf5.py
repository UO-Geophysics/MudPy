"""Observation semantics, source-block assembly, and direct NumPy reuse."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from mudpy import geotomo, geotomo_io as io, insartools
from mudpy.geotomo_config import RunConfig, _inputs, _signatures, _code_identity


class ParallelForwardTests(unittest.TestCase):
    def test_source_workers_overlap_with_bounded_ordered_results(self):
        from threading import Barrier, Lock
        from mudpy.parallel import _geotomo_source_map
        barrier, lock = Barrier(2, timeout=5), Lock()
        active = 0
        maximum = 0
        def work(source):
            nonlocal active, maximum
            with lock:
                active += 1
                maximum = max(maximum, active)
            if source < 2:
                barrier.wait()
            with lock:
                active -= 1
            return source
        self.assertEqual(list(_geotomo_source_map(work, list(range(5)), 2)), list(range(5)))
        self.assertEqual(maximum, 2)
        # Requests larger than the source count and serial execution are valid.
        self.assertEqual(list(_geotomo_source_map(lambda x: x, [7], 8)), [7])
        self.assertEqual(list(_geotomo_source_map(lambda x: x, [7], 1)), [7])

    def test_worker_failure_waits_for_active_sources_and_stops_next_batch(self):
        from threading import Barrier, Event
        from mudpy.parallel import _geotomo_source_map
        barrier, finished = Barrier(2, timeout=5), Event()
        def work(source):
            if source == 0:
                barrier.wait()
                raise RuntimeError('failed source')
            if source == 1:
                barrier.wait()
                finished.set()
                return source
            self.fail('Scheduled another batch after failure')
        with self.assertRaisesRegex(RuntimeError, 'failed source'):
            list(_geotomo_source_map(work, [0, 1, 2], 2))
        self.assertTrue(finished.is_set())

    def test_parallel_synthetics_have_one_writer_and_leave_failures_incomplete(self):
        from threading import Barrier, get_ident
        from mudpy.parallel import run_parallel_synthetics_geotomo
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)/'site'
            for name in ('data/model_info', 'data/station_info', 'GFs/static'):
                (project/name).mkdir(parents=True)
            # Repeated depths and unsorted IDs must retain source-file order.
            sources = np.array([[19, 13, 43, 3], [3, 13.1, 43.1, 3], [8, 13.2, 43.2, 5]])
            np.savetxt(project/'data/model_info/mpi_source.0.fault', sources)
            stations = project/'data/station_info/one.sta'
            stations.write_text('A 13 43\n')
            for source in sources:
                values = np.zeros((1, 14))
                values[0, [0, 10]] = 1
                values[0, 1] = source[0]
                np.savetxt(project/('GFs/static/earth.mod.static.%.4f.sub%04d.insar' % (source[3], source[0])), values)
            barrier = Barrier(2, timeout=5)
            original_tensor = geotomo._tensor_static_neu
            def tensor(statics, azimuths, mechanisms):
                marker = statics[0, 1]
                if marker in (19, 3):
                    barrier.wait()
                return original_tensor(statics, azimuths, mechanisms)
            owner, writers = get_ident(), []
            original_setitem = h5py.Dataset.__setitem__
            def write(dataset, key, value):
                writers.append(get_ident())
                original_setitem(dataset, key, value)
            def run(name):
                run_parallel_synthetics_geotomo(temporary+'/', 'site', str(stations),
                    'earth.mod', False, (), True, 0, 1, response_file=project/name, processes=2)
            with patch('mudpy.green.src2sta', return_value=(np.array([1.]), np.array([0.]))), \
                    patch.object(geotomo, '_tensor_static_neu', side_effect=tensor), \
                    patch('subprocess.run', side_effect=AssertionError('Unexpected syn process')), \
                    patch.object(h5py.Dataset, '__setitem__', new=write):
                run('ok.h5')
            self.assertEqual(set(writers), {owner})
            with h5py.File(str(project/'ok.h5'), 'r') as file:
                self.assertTrue(file.attrs['complete'])
                np.testing.assert_array_equal(file['sources'][()], sources)
                np.testing.assert_allclose(file['displacement'][:, 0, 2, 0], -sources[:, 0]/6*(10**16.6*1e-15))
            with patch('mudpy.green.src2sta', return_value=(np.array([1.]), np.array([0.]))), \
                    patch.object(geotomo, '_tensor_static_neu', side_effect=ValueError('failed tensor calculation')):
                with self.assertRaisesRegex(ValueError, 'failed tensor calculation'):
                    run('failed.h5')
            with h5py.File(str(project/'failed.h5'), 'r') as file:
                self.assertFalse(file.attrs['complete'])
            observations = {'ids': np.array(['A']), 'coordinates': np.array([[13., 43.]]),
                            'projection_neu': np.array([[0., 0., 1.]])}
            with self.assertRaisesRegex(ValueError, 'Incomplete'):
                io.assemble_matrix(project/'failed.h5', observations, sources)


class StorageTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)

    def test_preprocessing_import_preserves_means_membership_and_enu_order(self):
        x, y = np.meshgrid(np.linspace(13, 13.01, 4), np.linspace(43, 43.01, 4))
        angles = np.linspace(.1, .3, 16).reshape(4, 4)
        look = np.stack((np.sin(angles), np.zeros_like(x), np.cos(angles)), axis=-1)
        result = insartools.preprocess_insar(np.arange(16.).reshape(4, 4), x, y, look,
            units='mm', los_positive='toward', look_direction='ground_to_satellite',
            coordinate_system='Geographic longitude/latitude degrees; datum unspecified',
            min_size=4, max_size=4, threshold=1.)
        package = insartools.write_preprocessing(result, self.folder/'package')
        filename = geotomo.import_observations(package, self.folder/'observations.h5')
        observations = io.read_observations(filename)
        np.testing.assert_array_equal(observations['los_m'], result['sampling']['los'])
        np.testing.assert_array_equal(observations['coordinates'], result['sampling']['xy'])
        np.testing.assert_array_equal(observations['projection_neu'], result['sampling']['look_enu'][:, [1, 0, 2]])
        self.assertLess(np.linalg.norm(observations['projection_neu'][0]), 1.)
        self.assertIsNone(observations['sigma_m'])
        with h5py.File(str(filename), 'r') as file:
            np.testing.assert_array_equal(file['sampling/source_indices'][()], result['sampling']['source_indices'])
            self.assertNotIn('original_m', file)
        with self.assertRaises(OSError):
            geotomo.import_observations(package, filename)
        table = package/'samples.txt'
        # Change one table value explicitly, independently of text formatting.
        samples = np.loadtxt(table, ndmin=2)
        samples[0, 2] += .001
        np.savetxt(table, samples)
        with self.assertRaisesRegex(ValueError, 'does not match'):
            geotomo.import_observations(package, self.folder/'mixed_package.h5')

    def test_observation_units_geometry_and_unknown_uncertainty_fail_explicitly(self):
        arguments = (self.folder/'obs.h5', ['A'], [[13, 43]], [.001])
        for projection in ([[0, 0, 0]], [[0, 0, 2]], [[0, np.nan, 1]]):
            with self.assertRaises(ValueError):
                io.write_observations(*arguments, projection, metadata={})
        io.write_observations(*arguments, [[.1, -.5, .8]], metadata={'datum': 'unspecified'})
        with h5py.File(str(arguments[0]), 'r+') as file:
            file.attrs['displacement_units'] = 'mm'
        with self.assertRaisesRegex(ValueError, 'units'):
            io.read_observations(arguments[0])

    def test_full_and_fixed_projection_and_ordering(self):
        sources = np.array([[9, 13, 43, 3], [2, 13.01, 43.01, 5.]])
        ids = ['B', 'A', 'C']
        coordinates = np.array([[13, 43], [13.02, 43.01], [13.03, 43.02]])
        look = np.array([[.1, -.5, .7], [.2, .3, .8], [0, 0, 1.]])
        path = self.folder/'obs.h5'
        io.write_observations(path, ids, coordinates, np.ones(3), look, metadata={})
        obs = io.read_observations(path)
        blocks = np.arange(108., dtype=float).reshape(2, 3, 3, 6)/1000
        library = self.folder/'elastic.h5'
        with io.response_writer(library, sources, ids, coordinates, {}) as (data, complete):
            for k in range(2):
                data[k] = blocks[k]
                complete[k] = True
        expected = np.array([[look[r].dot(blocks[s, r, :, t])
                              for s in range(2) for t in range(6)] for r in range(3)])
        np.testing.assert_allclose(io.assemble_matrix(library, obs, sources), expected, atol=1e-16)
        mechanism = np.array([-1, .2, 0, -1, 0, -1])
        np.testing.assert_allclose(io.assemble_matrix(library, obs, sources, mechanism),
                                  np.column_stack([expected[:, 6*k:6*k+6].dot(mechanism) for k in range(2)]))
        with self.assertRaisesRegex(ValueError, 'ordering'):
            io.assemble_matrix(library, obs, sources[::-1])
        with h5py.File(str(library), 'r+') as file:
            file['displacement'][0, 0, 0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, 'Nonfinite'):
            io.assemble_matrix(library, obs, sources)

    def test_failed_library_remains_incomplete_and_cannot_be_assembled(self):
        library = self.folder/'elastic.h5'
        with self.assertRaisesRegex(RuntimeError, 'worker failure'):
            with io.response_writer(library, [[1, 13, 43, 3]], ['A'], [[13, 43]], {}):
                raise RuntimeError('worker failure')
        with h5py.File(str(library), 'r') as file:
            self.assertFalse(file.attrs['complete'])
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            io.assemble_matrix(library, {'ids': ['A']}, [[1, 13, 43, 3]])


class HDF5RunTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.project = Path(self.temporary.name)/'site'
        self.template = geotomo.init_project(self.project)
        self.text = self.template.read_text()
        self.executable = self.project/'test_executable'
        self.executable.write_text('identity only; numerical worker is mocked')
        self.which = patch('mudpy.geotomo_config.shutil.which', return_value=str(self.executable))
        self.which.start()
        self.addCleanup(self.which.stop)

    def config(self, name, steps, mode='assemble', mechanism=None, keep_fk=True):
        text = self.text.replace('name = "example"', 'name = '+json.dumps(name))
        text = text.replace('steps = ["greens", "synthetics", "invert"]', 'steps = '+json.dumps(steps))
        text = text.replace('mode = "assemble"', 'name = "shared"\nmode = '+json.dumps(mode))
        text = text.replace('keep_fk_intermediates = true', 'keep_fk_intermediates = '+str(keep_fk).lower())
        if mechanism is not None:
            text = text.replace('source_model = "full_tensor"',
                                'source_model = "fixed_tensor"\nfixed_tensor = '+json.dumps(mechanism))
        path = self.project/'configs'/(name+'.toml')
        path.write_text(text)
        return path

    def worker(self, *args, **kwargs):
        sources = np.loadtxt(self.project/'data/model_info'/args[3], ndmin=2)
        obs = kwargs['observations']
        if args[5]:
            for source in sources:
                (self.project/'GFs/static'/('template.mod.static.%.4f.sub%04d.insar' %
                                          (source[3], source[0]))).write_text('native fixture')
                folder = self.project/('GFs/static/template.mod_%.4f.sub%04d' % (source[3], source[0]))
                folder.mkdir(exist_ok=True)
                (folder/'template.mod').write_bytes((self.project/'structure/template.mod').read_bytes())
                (folder/'junk.p').write_text('P travel-time fixture')
                (folder/'junk.s').write_text('S travel-time fixture')
        if not args[6]:
            return
        with io.response_writer(kwargs['response_file'], sources, obs['ids'], obs['coordinates'],
                                kwargs['response_metadata']) as (data, complete):
            data[:] = np.arange(1., 19.).reshape(1, 1, 3, 6)*.001
            complete[:] = True

    def test_native_retention_default_and_policy_independent_cache_signatures(self):
        path = self.config('retain', ['greens', 'synthetics'])
        kept = RunConfig.load(path)
        pruned_path = self.config('policy_only', ['greens', 'synthetics'], keep_fk=False)
        pruned = RunConfig.load(pruned_path)
        sources, _, geometry, _ = _inputs(kept)
        code = _code_identity()
        self.assertEqual(_signatures(kept, sources, geometry, code, {})[1][1],
                         _signatures(pruned, sources, geometry, code, {})[1][1])
        with patch.object(geotomo, '_inversion_gfs', side_effect=self.worker):
            geotomo.run(path)
        self.assertTrue((self.project/'GFs/static/template.mod.static.3.0000.sub0001.insar').exists())
        self.assertTrue((self.project/'GFs/static/template.mod_3.0000.sub0001/junk.p').exists())
        record = json.loads((self.project/'GFs/static/_geotomo_greens.json').read_text())
        self.assertNotIn('native_cleanup', record)

    def test_gnss_forward_uses_station_vectors_without_projection_or_inversion(self):
        io.write_gnss(self.project/'data/statics/gnss.h5', ['GPS1'], [[13., 43.]],
                      [[.01, -.02, .03]], [[.001, .002, .005]], [1000.], metadata={})
        path = self.config('gnss_forward', ['greens', 'synthetics'], keep_fk=False)
        text = path.read_text().replace('observations = "template.h5"', 'observations = "gnss.h5"')
        text = text.replace('observation_format = "hdf5"', 'observation_format = "gnss_hdf5"')
        path.write_text(text)
        config = RunConfig.load(path)
        sources, obs, geometry, _ = _inputs(config)
        self.assertNotIn('projection_neu', obs)
        self.assertEqual(geometry, [{'id': 'GPS1', 'lon': 13., 'lat': 43.}])
        with patch.object(geotomo, '_inversion_gfs', side_effect=self.worker), \
                patch.object(io, 'assemble_matrix', side_effect=AssertionError('Unexpected inversion')):
            record = geotomo.run(path)
        self.assertTrue(record.is_file())
        io.validate_responses(config.response_path, obs, sources)
        self.assertFalse((self.project/'GFs/static/template.mod.static.3.0000.sub0001.insar').exists())
        path.write_text(text.replace('steps = ["greens", "synthetics"]', 'steps = ["invert"]'))
        with self.assertRaisesRegex(ValueError, 'forward responses only'):
            RunConfig.load(path)

    def test_result_reader_changes_do_not_invalidate_numerical_caches(self):
        import hashlib
        from mudpy import geotomo_config
        config = RunConfig.load(self.template)
        sources, _, geometry, _ = _inputs(config)
        baseline = _code_identity()
        original_read = Path.read_bytes
        original_hash = geotomo_config._hash
        io_path = Path(io.__file__).resolve()
        def changed_reader(path):
            content = original_read(path)
            return content+b'\n# Reader-only change\n' if path.resolve() == io_path else content
        def changed_hash(path):
            return hashlib.sha256(changed_reader(Path(path))).hexdigest() if Path(path).resolve() == io_path else original_hash(path)
        with patch.object(Path, 'read_bytes', new=changed_reader), \
                patch.object(geotomo_config, '_hash', side_effect=changed_hash):
            changed = _code_identity()
        self.assertNotEqual(baseline['source_sha256']['geotomo_io'], changed['source_sha256']['geotomo_io'])
        self.assertEqual(baseline['io_cache_sha256'], changed['io_cache_sha256'])
        self.assertEqual(_signatures(config, sources, geometry, baseline, {})[1:],
                         _signatures(config, sources, geometry, changed, {})[1:])

    def test_pruning_preserves_library_assembly_fixed_projection_and_numpy_reuse(self):
        path = self.config('pruned', ['greens', 'synthetics'], keep_fk=False)
        unrelated = self.project/'GFs/static/user_notes.insar'
        unrelated.write_text('keep me')
        other_folder = self.project/'GFs/static/other_model_3.0000.sub0001'
        other_folder.mkdir()
        (other_folder/'notes.txt').write_text('also keep me')
        # Cleanup validates blocks without assembling a matrix.
        with patch.object(geotomo, '_inversion_gfs', side_effect=self.worker), \
                patch.object(io, 'assemble_matrix', side_effect=AssertionError('Unexpected matrix assembly')):
            resolved = geotomo.run(path)
        self.assertFalse((self.project/'GFs/static/template.mod.static.3.0000.sub0001.insar').exists())
        self.assertFalse((self.project/'GFs/static/template.mod_3.0000.sub0001').exists())
        self.assertEqual(unrelated.read_text(), 'keep me')
        self.assertEqual((other_folder/'notes.txt').read_text(), 'also keep me')
        native = json.loads((self.project/'GFs/static/_geotomo_greens.json').read_text())
        response = json.loads((self.project/'GFs/static/elastic.h5.json').read_text())
        self.assertEqual(native['availability'], 'pruned')
        self.assertEqual(response['native_cleanup']['state'], 'pruned')
        self.assertEqual(len(response['native_cleanup']['files']), 4)
        self.assertEqual(response['native_cleanup']['run'], 'pruned')
        self.assertFalse(json.loads(resolved.read_text())['settings']['execution']['keep_fk_intermediates'])
        assemble = self.config('assemble_pruned', ['invert'])
        reuse = self.config('reuse_pruned', ['invert'], 'reuse', keep_fk=False)
        with patch.object(geotomo, '_inversion_gfs', side_effect=AssertionError('Unexpected GF generation')):
            geotomo.run(assemble)
            with patch.object(io, 'assemble_matrix', side_effect=AssertionError('Unexpected assembly')):
                geotomo.run(reuse)
        models = self.project/'output/inverse_models/models'
        np.testing.assert_array_equal(np.loadtxt(models/'assemble_pruned.0000.inv'),
                                      np.loadtxt(models/'reuse_pruned.0000.inv'))
        obs = io.read_observations(self.project/'data/statics/template.h5')
        sources = np.loadtxt(self.project/'data/model_info/template.sources', ndmin=2)
        obs['projection_neu'][:] = [.2, -.3, .8]
        self.assertEqual(io.assemble_matrix(self.project/'GFs/static/elastic.h5', obs, sources,
                                            [1, 0, 0, 1, 0, 1]).shape, (1, 1))
        rebuild = self.config('needs_native', ['synthetics'])
        rebuild.write_text(rebuild.read_text().replace('name = "elastic"', 'name = "rebuilt_library"'))
        with self.assertRaisesRegex(ValueError, 'intentionally removed.*greens'):
            geotomo.run(rebuild)
        self.assertFalse((models/'needs_native.config.toml').exists())

    def test_greens_only_defers_cleanup_until_successful_hdf5_synthesis(self):
        greens = self.config('greens_only', ['greens'], keep_fk=False)
        synthesis = self.config('finish_synthesis', ['synthetics'], keep_fk=False)
        native = self.project/'GFs/static/template.mod.static.3.0000.sub0001.insar'
        with patch.object(geotomo, '_inversion_gfs', side_effect=self.worker):
            geotomo.run(greens)
            self.assertTrue(native.exists())
            self.assertFalse((self.project/'GFs/static/elastic.h5').exists())
            geotomo.run(synthesis)
        self.assertFalse(native.exists())
        self.assertTrue((self.project/'GFs/static/elastic.h5').exists())

    def test_forward_and_solver_failures_retain_native_intermediates(self):
        for stage in ('forward', 'solver'):
            with self.subTest(stage=stage):
                path = self.config(stage+'_failure', ['greens', 'synthetics'] if stage == 'forward'
                                   else ['greens', 'synthetics', 'invert'], keep_fk=False)
                path.write_text(path.read_text().replace('name = "elastic"', 'name = "'+stage+'_library"'))
                def forward(*args, **kwargs):
                    self.worker(*args, **kwargs)
                    if stage == 'forward':
                        raise RuntimeError('forward failed')
                with patch.object(geotomo, '_inversion_gfs', side_effect=forward), \
                        patch.object(geotomo, '_run_inversion', side_effect=RuntimeError('solver failed')):
                    with self.assertRaisesRegex(RuntimeError, stage+' failed'):
                        geotomo.run(path)
                self.assertTrue((self.project/'GFs/static/template.mod.static.3.0000.sub0001.insar').exists())
                self.assertTrue((self.project/'GFs/static/template.mod_3.0000.sub0001/junk.p').exists())
                self.assertFalse((self.project/'.geotomo.lock').exists())

    def test_invalid_library_prevents_cleanup_before_any_deletion(self):
        for corruption in ('incomplete', 'nan', 'geometry'):
            with self.subTest(corruption=corruption):
                path = self.config(corruption, ['greens', 'synthetics'], keep_fk=False)
                path.write_text(path.read_text().replace('name = "elastic"', 'name = "'+corruption+'_library"'))
                def invalid(*args, **kwargs):
                    self.worker(*args, **kwargs)
                    with h5py.File(str(kwargs['response_file']), 'r+') as file:
                        if corruption == 'incomplete':
                            file.attrs['complete'] = False
                        elif corruption == 'nan':
                            file['displacement'][0, 0, 0, 0] = np.nan
                        else:
                            file['sources'][0, 3] = 9
                with patch.object(geotomo, '_inversion_gfs', side_effect=invalid):
                    with self.assertRaises(ValueError):
                        geotomo.run(path)
                self.assertTrue((self.project/'GFs/static/template.mod.static.3.0000.sub0001.insar').exists())
                self.assertTrue((self.project/'GFs/static/template.mod_3.0000.sub0001/junk.p').exists())

    def test_unknown_files_and_symlinks_prevent_cleanup(self):
        outside = self.project/'personal.txt'
        outside.write_text('preserve this')
        folder = self.project/'GFs/static/template.mod_3.0000.sub0001'
        for kind in ('unknown', 'symlink'):
            path = self.config(kind, ['greens', 'synthetics'], keep_fk=False)
            path.write_text(path.read_text().replace('name = "elastic"', 'name = "'+kind+'_library"'))
            def modified(*args, **kwargs):
                self.worker(*args, **kwargs)
                if kind == 'unknown':
                    (folder/'personal.txt').write_text('keep this too')
                else:
                    (folder/'junk.s').unlink()
                    (folder/'junk.s').symlink_to(outside)
            with patch.object(geotomo, '_inversion_gfs', side_effect=modified):
                with self.assertRaisesRegex(ValueError, 'Unexpected FK working file'):
                    geotomo.run(path)
            self.assertTrue((self.project/'GFs/static/template.mod.static.3.0000.sub0001.insar').exists())
            self.assertTrue((folder/'template.mod').exists())
            self.assertEqual(outside.read_text(), 'preserve this')
            if kind == 'unknown':
                self.assertEqual((folder/'personal.txt').read_text(), 'keep this too')
                (folder/'personal.txt').unlink()  # Fixture only; permit the symlink case.

    def test_interrupted_cleanup_rejects_native_reuse_and_preserves_response_cache(self):
        path = self.config('interrupted', ['greens', 'synthetics'], keep_fk=False)
        model = (self.project/'GFs/static/template.mod_3.0000.sub0001/template.mod').resolve()
        original_unlink = Path.unlink
        def interrupted(target, *args, **kwargs):
            if target.resolve() == model:
                raise PermissionError('cleanup interrupted')
            return original_unlink(target, *args, **kwargs)
        with patch.object(geotomo, '_inversion_gfs', side_effect=self.worker), \
                patch.object(Path, 'unlink', new=interrupted):
            with self.assertRaisesRegex(PermissionError, 'cleanup interrupted'):
                geotomo.run(path)
        native_path = self.project/'GFs/static/_geotomo_greens.json'
        native = json.loads(native_path.read_text())
        response_path = self.project/'GFs/static/elastic.h5.json'
        response = json.loads(response_path.read_text())
        self.assertEqual(native['availability'], 'cleanup_pending')
        self.assertEqual(response['native_cleanup']['state'], 'cleanup_pending')
        with self.assertRaisesRegex(ValueError, 'cleanup was interrupted'):
            from mudpy.geotomo_config import _require_cache
            _require_cache(self.project, (native_path, native['inputs'],
                [self.project/name for name in native['files']]))
        _require_cache(self.project, (response_path, response['inputs'],
            [self.project/name for name in response['files']]))
        self.assertFalse((self.project/'.geotomo.lock').exists())

    def test_shifted_sources_follow_library_matrix_models_and_cached_reuse(self):
        self.text = self.text.replace('interface_policy = "error"', 'interface_policy = "below"')
        source_path = self.project/'data/model_info/template.sources'
        source_path.write_text('9 13 43 1\n')
        original_sources = source_path.read_bytes()
        (self.project/'structure/template.mod').write_text('1 3 5 2.5\n0 3.5 6 2.7\n')
        first = self.config('shifted', ['greens', 'synthetics', 'invert'])
        with patch.object(geotomo, '_inversion_gfs', side_effect=self.worker) as worker:
            resolved = geotomo.run(first)
        effective_file = self.project/'data/model_info/_geotomo_shifted.sources'
        self.assertEqual(worker.call_args[0][3], effective_file.name)
        self.assertEqual(source_path.read_bytes(), original_sources)
        self.assertEqual(np.loadtxt(effective_file)[3], 1.01)
        record = json.loads(resolved.read_text())
        self.assertEqual(record['source_geometry']['requested'][0][3], 1)
        self.assertNotIn('syn', record['executables'])
        self.assertEqual(record['synthesis']['method'], 'numpy_mt_radiat')
        self.assertEqual(record['source_geometry']['effective'][0][3], 1.01)
        self.assertIn(str(effective_file.resolve()), record['input_sha256'])
        with h5py.File(str(self.project/'GFs/static/elastic.h5'), 'r') as file:
            self.assertEqual(file['sources'][0, 3], 1.01)
            self.assertEqual(json.loads(file.attrs['metadata_json'])['source_geometry']['effective'][0][3], 1.01)
        model = self.project/'output/inverse_models/models/shifted.0000.inv'
        self.assertEqual(np.loadtxt(model)[3], 1.01)
        reuse = self.config('shifted_reuse', ['invert'], 'reuse')
        with patch.object(geotomo, '_inversion_gfs', side_effect=AssertionError('Unexpected GF generation')), \
                patch('mudpy.geotomo_io.assemble_matrix', side_effect=AssertionError('Unexpected assembly')):
            geotomo.run(reuse)
        reused = self.project/'output/inverse_models/models/shifted_reuse.0000.inv'
        np.testing.assert_array_equal(np.loadtxt(model), np.loadtxt(reused))
        # A different effective side cannot reuse the same matrix.
        wrong_side = self.config('wrong_side', ['invert'], 'reuse')
        wrong_side.write_text(wrong_side.read_text().replace('interface_policy = "below"',
                                                          'interface_policy = "above"'))
        with self.assertRaisesRegex(ValueError, 'do not match'):
            geotomo.run(wrong_side)
        self.assertFalse((self.project/'data/model_info/_geotomo_wrong_side.sources').exists())
        with self.assertRaisesRegex(FileExistsError, 'already exists'):
            geotomo.run(first)

    def test_complete_run_and_numpy_only_reuse_after_values_weights_solver_change(self):
        first = self.config('first', ['greens', 'synthetics', 'invert'])
        with patch.object(geotomo, '_inversion_gfs', side_effect=self.worker):
            geotomo.run(first)
        matrix = self.project/'GFs/matrices/shared.npy'
        original = matrix.read_bytes()
        # Even moving the elastic library cannot force a cached inversion to read it.
        (self.project/'GFs/static/elastic.h5').rename(self.project/'GFs/static/retained.h5')
        with h5py.File(str(self.project/'data/statics/template.h5'), 'r+') as file:
            file['los_m'][:] *= 2
            file['sigma_m'][:] *= 3
        reuse = self.config('reuse', ['invert'], 'reuse')
        reuse.write_text(reuse.read_text().replace('solver = "lstsq"', 'solver = "nnls"')
                         .replace('damping = [0.0]', 'damping = [0.01, 0.1]'))
        real_load = np.load
        with patch.object(geotomo, '_inversion_gfs', side_effect=AssertionError('Unexpected generation')), \
                patch('mudpy.geotomo_io.assemble_matrix', side_effect=AssertionError('Unexpected assembly')), \
                patch('mudpy.geotomo_config.np.load', wraps=real_load) as load:
            resolved = geotomo.run(reuse)
            self.assertEqual(load.call_count, 1)
            self.assertEqual(Path(load.call_args.args[0]).resolve(), matrix.resolve())
        self.assertEqual(matrix.read_bytes(), original)
        record = json.loads(resolved.read_text())
        self.assertIsNone(record['normalized_observations'])
        self.assertEqual(record['observation_count'], 1)
        for k in range(2):
            path = self.project/'output/inverse_models/statics'/('reuse.%04d.predictions.h5' % k)
            with h5py.File(str(path), 'r') as file:
                model = np.loadtxt(self.project/'output/inverse_models/models'/('reuse.%04d.inv' % k), ndmin=2)
                np.testing.assert_allclose(file['prediction_m'][()],
                    real_load(matrix).dot((model[:, 4:10]/10**16.6).reshape(-1)), rtol=1e-9)
        self.assertFalse((self.project/'.geotomo.lock').exists())

    def test_physical_library_signature_excludes_projection_and_fixed_mechanism(self):
        path = self.config('first', ['greens', 'synthetics'])
        with patch.object(geotomo, '_inversion_gfs', side_effect=self.worker):
            geotomo.run(path)
        first = RunConfig.load(path)
        def signatures(config):
            sources, _, geometry, _ = _inputs(config)
            return _signatures(config, sources, geometry, _code_identity(), {})
        baseline = signatures(first)
        with h5py.File(str(self.project/'data/statics/template.h5'), 'r+') as file:
            file['projection_neu'][:] = [.2, -.3, .8]
        fixed = self.config('fixed', ['invert'], mechanism=[-1, 0, 0, -1, 0, -1])
        changed = signatures(RunConfig.load(fixed))
        self.assertEqual(baseline[0][1], changed[0][1])
        self.assertEqual(baseline[1][1], changed[1][1])
        self.assertNotEqual(baseline[2][1], changed[2][1])
        with patch.object(geotomo, '_inversion_gfs', side_effect=AssertionError('Unexpected generation')):
            geotomo.run(fixed)
        self.assertEqual(np.load(self.project/'GFs/matrices/shared.npy').shape, (1, 1))

    def test_missing_sigma_and_duplicate_library_fail_before_run_records(self):
        with h5py.File(str(self.project/'data/statics/template.h5'), 'r+') as file:
            del file['sigma_m']
        config = self.config('unknown_sigma', ['greens'])
        with self.assertRaisesRegex(ValueError, 'sigma_m'):
            geotomo.run(config)
        self.assertFalse((self.project/'output/inverse_models/models/unknown_sigma.config.toml').exists())
        config.write_text(config.read_text().replace('weighting = "diagonal"', 'weighting = "none"')
                          .replace('steps = ["greens"]', 'steps = ["greens", "synthetics"]'))
        (self.project/'GFs/static/elastic.h5').write_text('existing product')
        with self.assertRaisesRegex(FileExistsError, 'already exists'):
            geotomo.run(config)

    def test_inputs_changed_during_assembly_do_not_certify_the_matrix(self):
        prepare = self.config('prepare', ['greens', 'synthetics'])
        with patch.object(geotomo, '_inversion_gfs', side_effect=self.worker):
            geotomo.run(prepare)
        assemble = self.config('assemble', ['invert'])
        real_assemble = io.assemble_matrix
        def changed_input(*args):
            matrix = real_assemble(*args)
            with h5py.File(str(self.project/'data/statics/template.h5'), 'r+') as file:
                file['los_m'][:] *= 2
            return matrix
        with patch('mudpy.geotomo_io.assemble_matrix', side_effect=changed_input):
            with self.assertRaisesRegex(RuntimeError, 'Input changed'):
                geotomo.run(assemble)
        self.assertFalse((self.project/'GFs/matrices/shared.npy.json').exists())
        self.assertFalse((self.project/'.geotomo.lock').exists())


if __name__ == '__main__':
    unittest.main()
