"""Pure-DC geometry, analytic Jacobians and controlled nonlinear recovery."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from mudpy.analysis import double_couple_tensors, decompose_geotomo_tensors
from mudpy.geotomo_dc import solve, dc_reference, _jacobian_operator
from mudpy import geotomo, geotomo_io as io, geotomo_config as config


NDK = '''PDEW 2006/04/09 20:50:46.0 -20.45  -70.24  34.6 5.5 5.8 NEAR COAST OF NORTHERN C
C200604092050A   B: 88  166  40 S: 96  189  50 M: 41   52 125 CMT: 1 TRIHD:  1.8
CENTROID:      5.3 0.1 -20.46 0.01  -70.73 0.01  39.0  0.4 FREE S-20060726112355
24  4.180 0.069 -1.700 0.046 -2.480 0.060 -1.050 0.052 -2.410 0.075 -2.280 0.038
V10   4.975 73 100   0.120  8 216  -5.095 15 308   5.035  49 30  106 211 61   81
'''


class DoubleCoupleTests(unittest.TestCase):
    def test_geometry_derivatives_and_adjoint(self):
        base = np.diag([-1., 0, 1])
        p = np.array([[2., .3, -.5, .1], [.7, 0, 0, 0], [1, np.pi, -.2, .4]])
        coefficients, derivative = double_couple_tensors(p, base)
        for k in range(4):
            plus, minus = p.copy(), p.copy()
            plus[:, k] += 1e-6
            minus[:, k] -= 1e-6
            numeric = (double_couple_tensors(plus, base)[0]-double_couple_tensors(minus, base)[0])/2e-6
            np.testing.assert_allclose(derivative[:, :, k], numeric, atol=2e-9, rtol=2e-8)
        parts = decompose_geotomo_tensors(coefficients)
        np.testing.assert_allclose(parts['moment_nm']['full'], p[:, 0])
        np.testing.assert_allclose(parts['moment_nm']['iso'], 0, atol=1e-15)
        np.testing.assert_allclose(parts['moment_nm']['clvd'], 0, atol=2e-15)
        rng = np.random.default_rng(25)
        G = rng.normal(size=(30, 18))
        J = _jacobian_operator(G, derivative, .2)
        h, v = rng.normal(size=12), rng.normal(size=33)
        self.assertAlmostEqual(float(v @ J.matvec(h)), float(h @ J.rmatvec(v)), places=12)

    def test_noiseless_recovery_damping_units_and_reproducible_starts(self):
        rng = np.random.default_rng(8)
        G = rng.normal(size=(60, 18))
        p = np.array([[2., .2, -.1, .1], [.7, -.2, .1, 0], [1.3, .4, 0, -.15]])
        truth, _ = double_couple_tensors(p, np.diag([-1., 0, 1]))
        data = G @ truth.ravel()
        reference = [-1, 0, 0, 0, 0, 1]
        solution, report = solve(G, data, 0, reference, starts=2, seed=71)
        np.testing.assert_allclose(solution[:, 0], truth.ravel(), atol=1e-7)
        again, repeated = solve(G*1e-4, data*1e-4, 0, reference, starts=2, seed=71)
        np.testing.assert_allclose(again, solution, atol=1e-7)
        self.assertEqual(report['selected_start'], repeated['selected_start'])
        damped, info = solve(G, data, 2., reference, starts=2, seed=71)
        self.assertGreater(np.linalg.norm(G @ damped[:, 0]-data), 0.01)
        self.assertLess(info['starts'][info['selected_start']]['moment_l2'], np.linalg.norm(p[:, 0]))
        selected = info['starts'][info['selected_start']]
        self.assertAlmostEqual(2*selected['cost'], selected['residual_l2']**2+4*selected['moment_l2']**2)

    def test_failed_convergence_is_loud_and_zero_data_are_exact(self):
        G = np.random.default_rng(12).normal(size=(30, 6))
        p = np.array([[2., .7, .3, -.1]])
        truth, _ = double_couple_tensors(p, np.diag([-1., 0, 1]))
        with self.assertRaisesRegex(RuntimeError, 'No DC solver start') as failed:
            solve(G, G @ truth.ravel(), 0, [-1, 0, 0, 0, 0, 1], starts=2, max_nfev=1)
        self.assertEqual(len(failed.exception.diagnostics), 2)
        self.assertFalse(any(item['success'] for item in failed.exception.diagnostics))
        solution, report = solve(G, np.zeros(30), 0, [-1, 0, 0, 0, 0, 1], starts=1)
        np.testing.assert_array_equal(solution, 0)
        self.assertEqual(report['starts'][0]['cost'], 0)
        self.assertTrue(report['global_optimum_established'])
        for bad in ([0]*6, [1, 0, 0, 1, 0, 1], [-1, 0, 0, -1, 0, 2]):
            with self.assertRaisesRegex(ValueError, 'no defined DC'):
                dc_reference(bad)


class DCProjectTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.project = Path(self.temporary.name).resolve()/'site'
        template = geotomo.init_project(self.project)
        text = template.read_text().replace('steps = ["greens", "synthetics", "invert"]', 'steps = ["invert"]')
        text = text.replace('mode = "assemble"', 'name = "shared"\nmode = "reuse"')
        text = text.replace('weighting = "diagonal"', 'weighting = "none"')
        self.full = self.project/'configs/full.toml'
        self.full.write_text(text)
        np.savetxt(self.project/'data/model_info/template.sources', [[19, 13, 43, 3], [4, 13.01, 43, 5]])
        self.matrix = np.random.default_rng(93).normal(size=(30, 12))*.01
        truth, _ = double_couple_tensors(np.array([[2, .2, -.1, 0], [1, -.2, 0, .1]]), np.diag([-1., 0, 1]))
        io.write_observations(self.project/'data/statics/synthetic.h5', ['R%d' % k for k in range(30)],
            np.column_stack((np.linspace(12.9, 13.1, 30), np.full(30, 43.))), self.matrix @ truth.ravel(),
            np.tile([0, 0, 1.], (30, 1)), sigma_m=np.full(30, .001), metadata={'synthetic': True})
        text = text.replace('observations = "template.h5"', 'observations = "synthetic.h5"')
        self.full.write_text(text)
        self.matrix_file = self.project/'GFs/matrices/shared.npy'
        np.save(self.matrix_file, self.matrix, allow_pickle=False)
        baseline = config.RunConfig.load(self.full)
        sources, _, geometry, _ = config._inputs(baseline)
        specification = config._signatures(baseline, sources, geometry, config._code_identity(), {})[2]
        config._record_cache(self.project, specification)
        self.cache_file = specification[0]
        self.cache_before = self.cache_file.read_bytes()
        self.dc_text = text.replace('name = "example"', 'name = "dc_synthetic"')
        self.dc_text = self.dc_text.replace('source_model = "full_tensor"', 'source_model = "double_couple"')
        self.dc_text = self.dc_text.replace('solver = "lstsq"', 'solver = "least_squares"')
        self.dc_text += '\n[inversion.double_couple]\ninitial_tensor = [-1,0,0,0,0,1]\nstarts = 2\nseed = 71\n'
        self.dc = self.project/'configs/dc.toml'
        self.dc.write_text(self.dc_text)

    def test_toml_reuses_full_matrix_and_archives_checked_dc_results(self):
        with patch.object(config.shutil, 'which', return_value=None), \
             patch.object(geotomo, '_inversion_gfs', side_effect=AssertionError('Forward calculation requested')), \
             patch.object(geotomo, '_run_inversion', side_effect=AssertionError('Linear solver requested')):
            record = json.loads(geotomo.run(self.dc).read_text())
        self.assertEqual(record['settings']['inversion']['penalty'], 'scalar_moment_l2')
        self.assertIn('geotomo_dc', record['code']['source_sha256'])
        result = io.read_solution(self.project, 'dc_synthetic', 0)
        np.testing.assert_array_equal(result['model']['ids'], [19, 4])
        self.assertLess(result['residual_rms_mm'], 1e-4)
        self.assertTrue(result['solver_diagnostics']['success'])
        self.assertEqual(self.cache_file.read_bytes(), self.cache_before)
        np.testing.assert_array_equal(np.load(self.matrix_file), self.matrix)
        self.assertFalse((self.project/'.geotomo.lock').exists())
        path = self.project/'output/inverse_models/models/dc_synthetic.0000.solver.json'
        path.write_text(path.read_text()+' ')
        with self.assertRaisesRegex(ValueError, 'diagnostics changed'):
            io.read_solution(self.project, 'dc_synthetic', 0)

    def test_exhausted_solver_writes_diagnostics_without_a_model(self):
        self.dc.write_text(self.dc_text+'max_nfev = 1\n')
        with patch.object(config.shutil, 'which', return_value=None):
            with self.assertRaisesRegex(RuntimeError, 'No DC solver start'):
                geotomo.run(self.dc)
        folder = self.project/'output/inverse_models/models'
        diagnostic = json.loads((folder/'dc_synthetic.0000.solver.json').read_text())
        self.assertFalse(diagnostic['success'])
        self.assertFalse((folder/'dc_synthetic.0000.inv').exists())
        self.assertFalse((self.project/'.geotomo.lock').exists())

    def test_weighted_sweep_readers_and_penalty_label(self):
        self.dc.write_text(self.dc_text.replace('weighting = "none"', 'weighting = "diagonal"')
                           .replace('damping = [0.0]', 'damping = [0.0, 10.0]'))
        with patch.object(config.shutil, 'which', return_value=None):
            geotomo.run(self.dc)
        rows = io.read_run(self.project, 'dc_synthetic')['solutions']
        self.assertGreater(rows[1]['residual_l2'], rows[0]['residual_l2'])
        self.assertLess(rows[1]['model_l2'], rows[0]['model_l2'])
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        with matplotlib.rc_context():
            from mudpy import view
        figure, axes = view.geotomo_tradeoff(self.project, 'dc_synthetic')
        self.assertEqual(axes[0].get_xlabel(), 'L2 norm of M0 / Mref')
        plt.close(figure)

    def test_catalog_units_coordinates_preflight_and_hashes(self):
        path = self.project/'data/model_info/event.ndk'
        path.write_text(NDK)
        initial = io.read_initial_tensor(path)
        np.testing.assert_allclose(initial['tensor_nm'], np.array([-1.7, 2.28, -1.05, -2.48, 2.41, 4.18])*1e17)
        self.dc.write_text(self.dc_text.replace('initial_tensor = [-1,0,0,0,0,1]', 'initial_catalog = "event.ndk"'))
        with patch.object(config.shutil, 'which', return_value=None):
            record = json.loads(geotomo.run(self.dc).read_text())
        self.assertEqual(record['input_sha256'][str(path)], hashlib.sha256(path.read_bytes()).hexdigest())
        self.assertIn('C200604092050A', record['double_couple_initialization']['event_id'])
        path.write_text(NDK+NDK)
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            io.read_initial_tensor(path)

    def test_configuration_and_solver_identity_are_separate_from_matrix_basis(self):
        for old, new in [('starts = 2', 'starts = true'), ('seed = 71', 'seed = -1'),
                         ('initial_tensor = [-1,0,0,0,0,1]', 'initial_tensor = [1,0,0,1,0,1]'),
                         ('solver = "least_squares"', 'solver = "lstsq"')]:
            self.dc.write_text(self.dc_text.replace(old, new))
            with self.assertRaises(ValueError):
                config.RunConfig.load(self.dc)
        self.dc.write_text(self.dc_text)
        dc = config.RunConfig.load(self.dc)
        sources, _, geometry, _ = config._inputs(dc)
        code = config._code_identity()
        full_spec = config._signatures(config.RunConfig.load(self.full), sources, geometry, code, {})
        dc_spec = config._signatures(dc, sources, geometry, code, {})
        self.assertEqual(full_spec, dc_spec)
        code['source_sha256']['geotomo_dc'] = 'different solver implementation'
        self.assertEqual(dc_spec, config._signatures(dc, sources, geometry, code, {}))
        code['source_sha256']['geotomo'] = 'different forward implementation'
        self.assertNotEqual(dc_spec[2], config._signatures(dc, sources, geometry, code, {})[2])


if __name__ == '__main__':
    unittest.main()
