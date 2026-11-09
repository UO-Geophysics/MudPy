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

    def config(self, name, steps, mode='assemble', mechanism=None):
        text = self.text.replace('name = "example"', 'name = '+json.dumps(name))
        text = text.replace('steps = ["greens", "synthetics", "invert"]', 'steps = '+json.dumps(steps))
        text = text.replace('mode = "assemble"', 'name = "shared"\nmode = '+json.dumps(mode))
        if mechanism is not None:
            text = text.replace('source_model = "full_tensor"',
                                'source_model = "fixed_tensor"\nfixed_tensor = '+json.dumps(mechanism))
        path = self.project/'configs'/(name+'.toml')
        path.write_text(text)
        return path

    def worker(self, *args, **kwargs):
        sources = np.loadtxt(self.project/'data/model_info/template.sources', ndmin=2)
        obs = kwargs['observations']
        for source in sources:
            (self.project/'GFs/static'/('template.mod.static.%.4f.sub%04d.insar' %
                                      (source[3], source[0]))).write_text('native fixture')
        with io.response_writer(kwargs['response_file'], sources, obs['ids'], obs['coordinates'],
                                kwargs['response_metadata']) as (data, complete):
            data[:] = np.arange(1., 19.).reshape(1, 1, 3, 6)*.001
            complete[:] = True

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
