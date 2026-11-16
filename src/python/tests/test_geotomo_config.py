"""Configuration, ownership, provenance, and legacy-layout compatibility."""
from dataclasses import FrozenInstanceError
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from mudpy import geotomo, runslip
from mudpy.geotomo_config import (
    RunConfig, _inputs, _signatures, _record_cache, _require_cache,
    _code_identity, _effective_sources,
)


class ConfigTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.parent = Path(self.temporary.name).resolve()
        self.project = self.parent/'site'
        self.config_file = geotomo.init_project(self.project)
        # Continue exercising the supported text-list interface independently
        # of the canonical HDF5 starter, which has its own integration tests.
        self.original = self.config_file.read_text().replace('observations = "template.h5"',
            'observations = "template.gflist"').replace('observation_format = "hdf5"',
            'observation_format = "gflist"')
        self.original = self.original.replace('[responses]\n# Basename without .h5, in GFs/static. A new library requires a new name.\nname = "elastic"\n', '')
        self.config_file.write_text(self.original)

    def configure(self, *changes):
        text = self.original
        for before, after in changes:
            self.assertIn(before, text)
            text = text.replace(before, after)
        self.config_file.write_text(text)
        return RunConfig.load(self.config_file)

    def products(self, config):
        source, _, geometry, _ = _inputs(config)
        code = {'source_sha256': {name: 'code' for name in ('geotomo', 'parallel', 'green')}}
        return _signatures(config, source, geometry, code, {})

    def test_project_lock_excludes_competing_runs_and_preserves_existing_owner(self):
        lock = self.project/'.geotomo.lock'
        with geotomo._project_lock(str(self.parent)+'/', self.project.name):
            self.assertTrue(lock.exists())
            with self.assertRaises(FileExistsError):
                geotomo.run(self.config_file)
            self.assertTrue(lock.exists())
        self.assertFalse(lock.exists())
        lock.write_text('existing owner')
        with self.assertRaises(FileExistsError):
            geotomo.run(self.config_file)
        self.assertEqual(lock.read_text(), 'existing owner')

    def test_minimal_version_one_toml_records_canonical_provenance(self):
        # Independent minimal schema-1 configuration, with no workflow-name key.
        original = b'''schema_version = 1
[project]
directory = ".."
[run]
name = "minimal_config"
steps = ["greens"]
[inputs]
sources = "template.sources"
earth_model = "template.mod"
observations = "template.gflist"
[forward]
earth_flattening = false
[inversion]
source_model = "full_tensor"
damping = [0.0]
'''
        self.config_file.write_bytes(original)
        config = RunConfig.load(self.config_file)

        def worker(*args):
            for product in self.products(config)[0][2]:
                product.write_text('verified-test-product\n')

        with patch.object(geotomo, '_inversion_gfs', side_effect=worker):
            resolved = geotomo.run(self.config_file)
        self.assertEqual(resolved.with_name('minimal_config.config.toml').read_bytes(), original)
        record = json.loads(resolved.read_text())
        self.assertEqual(Path(record['normalized_observations']).name,
                         '_geotomo_minimal_config.gflist')
        self.assertIn('geotomo', record['code']['source_sha256'])
        self.assertTrue((self.project/'GFs/static/_geotomo_greens.json').is_file())

    def test_provenance_hashes_numerical_and_configuration_implementations(self):
        from mudpy import geotomo_config
        code = _code_identity()['source_sha256']
        for name, module in (('geotomo', geotomo), ('geotomo_config', geotomo_config)):
            self.assertEqual(code[name], hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest())

    def test_initializer_refuses_existing_project_and_templates_are_valid(self):
        config = RunConfig.load(self.config_file)
        self.assertEqual(config.project, self.project.resolve())
        self.assertEqual(config.steps, ('greens', 'synthetics', 'invert'))
        self.assertEqual(config.processes, 1)
        self.assertTrue(config.keep_fk_intermediates)
        self.assertEqual(_inputs(config)[0].shape, (1, 4))
        with self.assertRaises(FrozenInstanceError):
            config.name = 'different'
        self.assertEqual(np.loadtxt(self.project/'data/model_info/template.sources').shape, (4,))
        with self.assertRaises(FileExistsError):
            geotomo.init_project(self.project)
        self.assertEqual(self.config_file.read_text(), self.original)
        self.assertTrue((self.project/'GFs/dynamic').is_dir())

    def test_worker_count_changes_provenance_without_invalidating_physical_caches(self):
        serial = self.configure()
        serial_signatures = self.products(serial)
        parallel = self.configure(('processes = 1', 'processes = 8'))
        self.assertEqual(parallel.settings()['execution']['processes'], 8)
        self.assertEqual(self.products(parallel), serial_signatures)
        self.assertTrue((self.project/'analysis/frequency').is_dir())

    def test_strict_schema_values_and_cross_field_validation(self):
        invalid = (
            ('schema_version = 1', 'schema_version = 2'),
            ('schema_version = 1', 'schema_version = true'),
            ('backend = "fk"', 'bakend = "fk"'),
            ('source_model = "full_tensor"', 'source_model = "fixed_tensor"'),
            ('source_model = "full_tensor"', 'source_model = "full_tensor"\nfixed_tensor = [1,0,0,1,0,1]'),
            ('processes = 1', 'processes = true'),
            ('processes = 1', 'processes = 0'),
            ('processes = 1', 'processes = -2'),
            ('processes = 1', 'processes = 1.5'),
            ('keep_fk_intermediates = true', 'keep_fk_intermediates = "false"'),
            ('keep_fk_intermediates = true', 'keep_fk_intermediates = 0'),
            ('keep_fk_intermediates = true', 'keep_fk_intermediates = false'),
            ('earth_flattening = false', 'earth_flattening = "false"'),
            ('interface_policy = "error"', 'interface_policy = "guess"'),
            ('interface_offset_km = 0.01', 'interface_offset_km = 0.0'),
            ('interface_offset_km = 0.01', 'interface_offset_km = 0.00001'),
            ('interface_offset_km = 0.01', 'interface_offset_km = true'),
            ('interface_offset_km = 0.01', 'interface_offset_km = nan'),
            ('dk = 0.05', 'dk = -1'),
            ('pmax = 1.0', 'pmax = 0.0'),
            ('damping = [0.0]', 'damping = [true]'),
            ('damping = [0.0]', 'damping = [-1.0]'),
            ('damping = [0.0]', 'damping = [inf]'),
            ('name = "example"', 'name = "../escape"'),
            ('sources = "template.sources"', 'sources = "../template.sources"'),
            ('steps = ["greens", "synthetics", "invert"]', 'steps = ["invert", "greens"]'),
            ('steps = ["greens", "synthetics", "invert"]', 'steps = ["greens", "invert"]'),
            ('mode = "assemble"', 'mode = "reuse"'),
        )
        for before, after in invalid:
            with self.subTest(change=after), self.assertRaises(ValueError):
                self.configure((before, after))
        self.config_file.write_text(self.original+'\n[unknown]\nvalue = 1\n')
        with self.assertRaisesRegex(ValueError, 'Unknown root'):
            RunConfig.load(self.config_file)
        self.assertEqual(list((self.project/'output/inverse_models/models').iterdir()), [])
        for exponent in (-1000, 1000):
            with self.assertRaises(FloatingPointError):
                self.configure(('damping = [0.0]', 'damping = {log10_min=%d, log10_max=%d, count=2}' %
                                (exponent, exponent)))

    def test_logspace_defaults_fixed_tensor_and_path_independence(self):
        config = self.configure(('damping = [0.0]', 'damping = {log10_min=-3, log10_max=-1, count=3}'),
                                ('source_model = "full_tensor"',
                                 'source_model = "fixed_tensor"\nfixed_tensor = [1,0,0,1,0,1]'))
        np.testing.assert_allclose(config.damping, [.001, .01, .1])
        self.assertEqual(config.fixed_tensor, (1, 0, 0, 1, 0, 1))
        cwd = Path.cwd()
        os.chdir(self.parent)
        try:
            self.assertEqual(RunConfig.load(self.config_file).project, self.project.resolve())
            _, rows, _, _ = _inputs(config)
            self.assertEqual(Path(rows[0][12]), self.project/'data/statics/template.los')
        finally:
            os.chdir(cwd)
        config = self.configure(('name = "example"', 'name = "new_experiment"'))
        self.assertEqual(config.matrix_name, 'new_experiment')

    def test_interface_default_fails_before_workers_or_run_records(self):
        (self.project/'structure/template.mod').write_text('1 3 5 2.5\n0 3.5 6 2.7\n')
        (self.project/'data/model_info/template.sources').write_text('9 13 43 1\n')
        with patch.object(geotomo, '_inversion_gfs') as worker:
            with self.assertRaisesRegex(ValueError, 'Source 9.*interface 1 km'):
                geotomo.run(self.config_file)
            worker.assert_not_called()
        self.assertEqual(list((self.project/'output/inverse_models/models').iterdir()), [])
        self.assertFalse((self.project/'.geotomo.lock').exists())

    def test_interface_sides_quantization_and_adjacent_layer_limits(self):
        model = np.array([[.5, 2, 4], [.5, 2.1, 4.1], [.5, 2.2, 4.2], [0, 3, 5]])
        requested = np.array([[9., 13., 43., 1.], [2., 13., 43., 1.24996]])
        for side, expected in [('above', .99), ('below', 1.01)]:
            config = self.configure(('interface_policy = "error"', 'interface_policy = "'+side+'"'))
            effective = _effective_sources(config, requested, model)
            np.testing.assert_array_equal(effective[:, 0], [9, 2])
            np.testing.assert_array_equal(effective[:, 3], [expected, 1.25])
            np.testing.assert_array_equal(requested[:, 3], [1, 1.24996])
            rounded_hit = requested[:1].copy()
            rounded_hit[0, 3] = .99996
            self.assertEqual(_effective_sources(config, rounded_hit, model)[0, 3], expected)
            for offset in (.5, .75):
                bad = self.configure(('interface_policy = "error"', 'interface_policy = "'+side+'"'),
                                     ('interface_offset_km = 0.01', 'interface_offset_km = '+str(offset)))
                with self.assertRaisesRegex(ValueError, 'adjacent layer'):
                    _effective_sources(bad, requested, model)
        config = self.configure(('interface_policy = "error"', 'interface_policy = "below"'))
        # Bottom interface may shift into the unbounded half-space.
        requested[0, 3] = 1.5
        self.assertEqual(_effective_sources(config, requested, model)[0, 3], 1.51)
        # Identical layers are not real FK interfaces.
        identical = np.array([[1., 3., 5.], [0., 3., 5.]])
        requested[0, 3] = 1.
        np.testing.assert_array_equal(_effective_sources(self.configure(), requested, identical),
                                      np.array([[9, 13, 43, 1.], [2, 13, 43, 1.25]]))
        requested[0, 3] = .00001
        with self.assertRaisesRegex(ValueError, 'nonpositive'):
            _effective_sources(config, requested, model)

    def test_interface_earth_flattening_and_physical_cache_identity(self):
        config = self.configure(('earth_flattening = false', 'earth_flattening = true'),
                                ('interface_policy = "error"', 'interface_policy = "below"'))
        model = np.array([[10., 3., 5.], [0., 3.5, 6.]])
        boundary = -6371*np.expm1(-(10*6371/(6371-5))/6371)
        sources = np.array([[1., 13., 43., float('%.4f' % boundary)]])
        result = _effective_sources(config, sources, model)
        self.assertEqual(result[0, 3], float('%.4f' % (boundary+.01)))
        strict = self.configure(('earth_flattening = false', 'earth_flattening = true'))
        with self.assertRaisesRegex(ValueError, 'FK interface'):
            _effective_sources(strict, sources, model)
        # Policy alone does not invalidate responses when effective geometry is unchanged.
        serial = self.configure()
        unshifted = self.products(serial)
        below = self.configure(('interface_policy = "error"', 'interface_policy = "below"'))
        self.assertEqual(self.products(below), unshifted)
        (self.project/'structure/template.mod').write_text('3 3 5\n0 3.5 6\n')
        above = self.configure(('interface_policy = "error"', 'interface_policy = "above"'))
        self.assertNotEqual(self.products(above)[0][1], self.products(below)[0][1])

    def test_invalid_input_and_stale_cache_fail_before_run_artifacts(self):
        config = self.configure(('steps = ["greens", "synthetics", "invert"]', 'steps = ["invert"]'))
        with self.assertRaisesRegex(ValueError, 'Missing cache record'):
            geotomo.run(self.config_file)
        self.assertEqual(list((self.project/'output/inverse_models/models').iterdir()), [])
        self.assertFalse((self.project/'.geotomo.lock').exists())
        data = self.project/'data/statics/template.los'
        data.write_text('0\n0\n0\n0\n')
        with self.assertRaisesRegex(ValueError, 'unit NEU'):
            geotomo.run(self.config_file)
        self.assertEqual(list((self.project/'output/inverse_models/models').iterdir()), [])

    def test_cache_signatures_separate_values_weights_geometry_and_mechanisms(self):
        config = RunConfig.load(self.config_file)
        greens, responses, matrix = self.products(config)
        product = matrix[2][0]
        np.save(product, np.array([[1.]]))
        _record_cache(self.project, matrix)
        _require_cache(self.project, matrix)
        observation = self.project/'data/statics/template.los'
        observation.write_text(observation.read_text().replace('0.001', '0.002'))
        self.assertEqual(self.products(config)[2][1], matrix[1])
        station_list = self.project/'data/station_info/template.gflist'
        station_list.write_text(station_list.read_text().replace(' 0.001 ', ' 0.002 '))
        self.assertEqual(self.products(config)[2][1], matrix[1])
        observation.write_text('0.002\n0\n1\n0\n')
        changed = self.products(config)
        self.assertEqual(changed[0][1], greens[1])
        self.assertNotEqual(changed[1][1], responses[1])
        with self.assertRaisesRegex(ValueError, 'do not match'):
            _require_cache(self.project, changed[2])
        np.save(product, np.array([[2.]]))
        with self.assertRaisesRegex(ValueError, 'changed cached product'):
            _require_cache(self.project, matrix)

    def test_snapshot_normalization_and_exclusion_survive_working_file_edits(self):
        config = self.configure(('steps = ["greens", "synthetics", "invert"]', 'steps = ["greens"]'),
                                ('processes = 1', 'processes = 3'))
        exact = self.config_file.read_bytes()

        def worker(*args):
            folder = self.project/'output/inverse_models/models'
            self.assertEqual((folder/'example.config.toml').read_bytes(), exact)
            self.assertTrue((self.project/'.geotomo.lock').exists())
            with self.assertRaises(FileExistsError):
                geotomo.run(self.config_file)
            for product in self.products(config)[0][2]:
                product.write_text('verified-test-product\n')

        with patch.object(geotomo, '_inversion_gfs', side_effect=worker) as execute:
            resolved = geotomo.run(self.config_file)
            self.assertEqual(execute.call_count, 1)
            self.assertEqual(execute.call_args[0][12], 3)
        record = json.loads(resolved.read_text())
        self.assertEqual(record['configuration_sha256'], hashlib.sha256(exact).hexdigest())
        self.assertEqual(record['source_count'], 1)
        self.assertEqual(record['settings']['execution']['processes'], 3)
        self.assertEqual(len(record['input_sha256']), 4)
        normalized = Path(record['normalized_observations'])
        self.assertIn(str(self.project/'data/statics/template.los'), normalized.read_text())
        self.config_file.write_text(self.config_file.read_text()+'\n# Later edit\n')
        self.assertEqual(resolved.with_name('example.config.toml').read_bytes(), exact)
        self.assertEqual(json.loads(resolved.read_text()), record)
        with self.assertRaisesRegex(FileExistsError, 'already exists'):
            geotomo.run(self.config_file)
        self.assertFalse((self.project/'.geotomo.lock').exists())

    def test_computation_failure_retains_snapshot_and_propagates(self):
        self.configure(('steps = ["greens", "synthetics", "invert"]', 'steps = ["greens"]'))
        with patch.object(geotomo, '_inversion_gfs', side_effect=RuntimeError('FK failed')):
            with self.assertRaisesRegex(RuntimeError, 'FK failed'):
                geotomo.run(self.config_file)
        self.assertTrue((self.project/'output/inverse_models/models/example.config.toml').is_file())
        self.assertFalse((self.project/'.geotomo.lock').exists())

    def test_input_changes_during_computation_do_not_get_a_cache_record(self):
        config = self.configure(('steps = ["greens", "synthetics", "invert"]', 'steps = ["greens"]'))

        def worker(*args):
            for product in self.products(config)[0][2]:
                product.write_text('test product\n')
            (self.project/'data/statics/template.los').write_text('0.002\n0\n1\n0\n')

        with patch.object(geotomo, '_inversion_gfs', side_effect=worker):
            with self.assertRaisesRegex(RuntimeError, 'Input changed'):
                geotomo.run(self.config_file)
        self.assertFalse((self.project/'GFs/static/_geotomo_greens.json').exists())
        self.assertFalse((self.project/'.geotomo.lock').exists())


class LegacyInitializationTests(unittest.TestCase):
    def test_exact_legacy_directories_templates_and_declined_reinitialization(self):
        # This fixture specifies the pre-refactor tree independently of its helper.
        expected = {'GFs', 'GFs/static', 'GFs/dynamic', 'GFs/matrices', 'GFs/tsunami', 'GFs/STFs',
                    'data', 'data/waveforms', 'data/statics', 'data/station_info', 'data/model_info',
                    'structure', 'plots', 'scripts', 'forward_models', 'output', 'output/inverse_models',
                    'output/inverse_models/statics', 'output/inverse_models/waveforms',
                    'output/inverse_models/models', 'output/forward_models', 'logs', 'analysis', 'analysis/frequency'}
        repo = Path(runslip.__file__).resolve().parents[3]
        with tempfile.TemporaryDirectory() as folder, patch.dict(os.environ, {'MUD': str(repo)}):
            runslip.init(folder+'/', 'legacy')
            project = Path(folder)/'legacy'
            self.assertEqual({str(path.relative_to(project)) for path in project.rglob('*') if path.is_dir()}, expected)
            for filename, destination in (('template.fault', 'data/model_info'), ('template.gflist', 'data/station_info'),
                                          ('template.sta', 'data/station_info'), ('template.mod', 'structure')):
                self.assertEqual((project/destination/filename).read_bytes(), (repo/'run'/filename).read_bytes())
            sentinel = project/'keep.txt'
            sentinel.write_text('keep')
            with patch('builtins.input', return_value='n') as prompt:
                runslip.init(folder+'/', 'legacy')
                self.assertEqual(prompt.call_count, 1)
            self.assertEqual(sentinel.read_text(), 'keep')
            with patch('builtins.input', side_effect=['y', 'n']) as prompt:
                runslip.init(folder+'/', 'legacy')
                self.assertEqual(prompt.call_count, 2)
            self.assertEqual(sentinel.read_text(), 'keep')
            with patch('builtins.input', side_effect=['y', 'y']) as prompt:
                runslip.init(folder+'/', 'legacy')
                self.assertEqual(prompt.call_count, 2)
            self.assertFalse(sentinel.exists())
            self.assertEqual({str(path.relative_to(project)) for path in project.rglob('*') if path.is_dir()}, expected)


if __name__ == '__main__':
    unittest.main()
