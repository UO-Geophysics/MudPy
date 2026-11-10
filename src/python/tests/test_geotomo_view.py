"""Canonical result readers, moment cancellation, and read-only quick viewers."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

with matplotlib.rc_context():
    from mudpy import view
from mudpy import geotomo_io as io


class ResultViewTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.addCleanup(plt.close, 'all')
        self.project = Path(self.temporary.name)
        self.models = self.project/'output/inverse_models/models'
        self.statics = self.project/'output/inverse_models/statics'
        self.models.mkdir(parents=True)
        self.statics.mkdir()
        self.run = 'example'
        self.settings = {'run': {'name': self.run},
                         'inversion': {'damping': [.1, 1.], 'weighting': 'none'}}
        (self.models/'example.resolved.json').write_text(json.dumps({'settings': self.settings}))
        self.observations = {'ids': np.array(['A', 'B', 'C']),
            'coordinates': np.array([[0, .01], [.02, .03], [-.03, .02]]),
            'los_m': np.array([.01, -.02, .03]), 'projection_neu': np.tile([0., 0., 1.], (3, 1)),
            'sigma_m': None, 'metadata': {}}
        for member, damping in enumerate([.1, 1.]):
            model = np.zeros((2, 17))
            model[:, :4] = [[1, 0, 0, 3], [2, .01, 0, 3]]
            model[:, 4] = np.array([1e19, -1e19])/(member+1)
            model[:, 10:16] = np.nan
            model[:, 16] = abs(model[:, 4])/np.sqrt(2)
            np.savetxt(self.models/('example.%04d.inv' % member), model)
            predicted = self.observations['los_m']/(2+member)
            l2 = np.linalg.norm(self.observations['los_m']-predicted)
            norm = np.linalg.norm(model[:, 4:10]/io.REFERENCE_MOMENT_NM)
            (self.models/('example.%04d.log' % member)).write_text(
                'lambda_spatial = %s\nL2 = %s\nLm = %s\n' % (damping, l2, norm))
            io.write_predictions(self.statics/('example.%04d.predictions.h5' % member),
                self.observations, predicted, {'run_name': self.run, 'damping': damping})

    def test_readers_distinguish_source_sum_from_tensor_sum(self):
        model = io.read_model(self.models/'example.0000.inv')
        self.assertAlmostEqual(model['sum_scalar_moment_nm']/1e19, np.sqrt(2))
        self.assertEqual(model['net_scalar_moment_nm'], 0)
        self.assertTrue(np.isnan(model['nodal_planes_deg']).all())
        solution = io.read_solution(self.project, self.run, '0000')
        self.assertEqual(solution['member'], 0)
        np.testing.assert_array_equal(solution['predictions']['residual_m'], self.observations['los_m']/2)
        rows = io.read_run(self.project, self.run)['solutions']
        self.assertEqual([r['member'] for r in rows], [0, 1])

    def test_inconsistent_predictions_models_and_missing_members_fail(self):
        filename = self.statics/'example.0001.predictions.h5'
        with h5py.File(str(filename), 'r+') as file:
            file['residual_m'][0] += 1
        with self.assertRaisesRegex(ValueError, 'residuals'):
            io.read_predictions(filename)
        filename.unlink()
        with self.assertRaises(OSError):
            io.read_run(self.project, self.run)
        model_path = self.models/'example.0000.inv'
        table = np.loadtxt(model_path)
        table[0, 16] *= 2
        np.savetxt(model_path, table)
        with self.assertRaisesRegex(ValueError, 'scalar moments'):
            io.read_model(model_path)
        with self.assertRaises(ValueError):
            io.read_solution(self.project, '../unsafe', 0)

    def test_fit_units_scales_and_no_file_or_global_style_side_effects(self):
        before = sorted(str(p) for p in self.project.rglob('*'))
        fonts = matplotlib.rcParams['font.family'].copy()
        with patch.object(plt, 'show', side_effect=AssertionError('Unexpected show')):
            figure, axes = view.geotomo_fit(self.project, self.run, 0, residual_limits=(-20, 20))
        self.assertEqual(len(axes), 3)
        np.testing.assert_allclose(axes[0].collections[0].get_array(), self.observations['los_m']*1000)
        np.testing.assert_allclose(axes[2].collections[0].get_array(), self.observations['los_m']*500)
        self.assertEqual(axes[0].collections[0].get_clim(), axes[1].collections[0].get_clim())
        self.assertEqual(axes[2].collections[0].get_clim(), (-20, 20))
        self.assertEqual(matplotlib.rcParams['font.family'], fonts)
        self.assertEqual(before, sorted(str(p) for p in self.project.rglob('*')))
        self.assertIs(figure, axes[0].figure)

    def test_tradeoff_and_source_projection_preserve_moment_budget(self):
        figure, axes = view.geotomo_tradeoff(self.project, self.run, reference_mw=6.6, highlight=[0])
        self.assertEqual(len(axes), 4)
        reference = axes[3].lines[-1]
        np.testing.assert_allclose(reference.get_ydata(), [1e19, 1e19], rtol=1e-14)
        with self.assertRaises(ValueError):
            view.geotomo_tradeoff(self.project, self.run, highlight=[3])
        figure, axes = view.geotomo_sources(self.project, self.run, 0)
        self.assertEqual(len(axes), 1)
        self.assertEqual(len(axes[0].collections[0].get_array()), 2)
        figure, axes = view.geotomo_sources(self.project, self.run, 0, kind='sections', center=(0, 0), strike=0)
        self.assertEqual(len(axes), 2)
        # Both sources project onto along=0, depth=3; none is hidden by overplotting.
        self.assertEqual(len(axes[0].collections[0].get_array()), 1)
        self.assertAlmostEqual(float(axes[0].collections[0].get_array()[0])/1e19, np.sqrt(2))
        self.assertAlmostEqual(float(axes[1].collections[0].get_array().sum())/1e19, np.sqrt(2))
        self.assertTrue(axes[0].yaxis_inverted())
        supplied, ax = plt.subplots()
        supplied.suptitle('Custom layout title')
        result, returned = view.geotomo_sources(self.project, self.run, 0, axes=[ax])
        self.assertIs(result, supplied)
        self.assertIs(returned[0], ax)
        self.assertEqual(supplied._suptitle.get_text(), 'Custom layout title')

    def test_zero_damping_zero_moments_and_numpy_member_indices(self):
        self.settings['inversion']['damping'][0] = 0.
        (self.models/'example.resolved.json').write_text(json.dumps({'settings': self.settings}))
        for member, damping in enumerate([0., 1.]):
            filename = self.models/('example.%04d.inv' % member)
            model = np.loadtxt(filename)
            model[:, 4:10] = 0
            model[:, 16] = 0
            np.savetxt(filename, model)
            (self.models/('example.%04d.log' % member)).write_text(
                'lambda_spatial = %s\nL2 = %s\nLm = 0\n' % (damping, np.linalg.norm(self.observations['los_m'])))
            with h5py.File(str(self.statics/('example.%04d.predictions.h5' % member)), 'r+') as file:
                file['prediction_m'][:] = 0
                file['residual_m'][:] = self.observations['los_m']
                metadata = json.loads(file.attrs['metadata_json'])
                metadata['run']['damping'] = damping
                file.attrs['metadata_json'] = json.dumps(metadata)
        figure, axes = view.geotomo_tradeoff(self.project, self.run)
        self.assertEqual(axes[1].get_xscale(), 'symlog')
        self.assertEqual(axes[3].get_yscale(), 'linear')
        figure, axes = view.geotomo_sources(self.project, self.run, np.int64(0))
        self.assertTrue(np.all(axes[0].collections[0].get_array() == 0))
        with self.assertRaises(ValueError):
            io.read_solution(self.project, self.run, True)


if __name__ == '__main__':
    unittest.main()
