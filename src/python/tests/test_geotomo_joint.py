"""Typed joint rows, exact nuisance elimination, DC recovery and canonical TOML runs."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from mudpy import geotomo, geotomo_io as io, geotomo_dc as dc, geotomo_config as cfg
from mudpy.analysis import double_couple_tensors


class JointTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.project = Path(self.temporary.name).resolve()/'site'
        template = geotomo.init_project(self.project)
        self.template = template.read_text()
        self.sources = np.array([[19, 13, 43, 3], [4, 13.01, 43, 5.]])
        np.savetxt(self.project/'data/model_info/template.sources', self.sources)
        x, y = np.meshgrid(np.linspace(-1, 1, 6), np.linspace(-1, 1, 4))
        x, y = x.ravel(), y.ravel()
        self.basis = np.column_stack((np.ones(len(x)), x, y, x*x, x*y, y*y))
        self.image = dict(ids=np.array(['I%d' % k for k in range(len(x))]),
            coordinates=np.c_[13+x*.01, 43+y*.01], projection_neu=np.tile([.2, -.4, .8], (len(x), 1)),
            los_m=np.zeros(len(x)), sigma_m=None, metadata={})
        self.gps = dict(ids=np.array(['G1', 'G2', 'G3']), coordinates=np.array([[13, 43], [13.1, 43], [13, 43.1]]),
            displacement_m=np.zeros((3, 3)), sigma_m=np.tile([.001, .002, .005], (3, 1)), height_m=np.zeros(3), metadata={})
        rng = np.random.default_rng(71)
        self.blocks = rng.normal(size=(2, len(x), 3, 6))*.001
        self.gps_blocks = rng.normal(size=(2, 3, 3, 6))*.001
        self.G = np.vstack((np.column_stack([np.einsum('rc,rct->rt', self.image['projection_neu'], b) for b in self.blocks]),
                            np.array([[self.gps_blocks[s, station, component, tensor]
                                       for s in range(2) for tensor in range(6)]
                                      for station in range(3) for component in range(3)])))
        self.initial = dc.dc_reference([-1, 0, 0, 0, 0, 1])
        p = np.array([[2., .1, -.15, .04], [.8, -.2, .1, .05]])
        self.truth, _ = double_couple_tensors(p, self.initial)
        self.beta = np.array([.01, -.005, .003, .002, -.001, .004])
        self.image['los_m'] = self.G[:len(x)] @ self.truth.ravel()+self.basis @ self.beta
        self.gps['displacement_m'] = (self.G[len(x):] @ self.truth.ravel()).reshape(-1, 3)
        self.data_file = self.project/'data/statics/joint.h5'
        io.write_joint(self.data_file, self.image, self.gps, self.basis, metadata={})
        self.obs = io.read_joint(self.data_file)

    def solve(self, obs=None, eta=1):
        return dc.solve_joint(self.G, self.obs if obs is None else obs, 1e-10, self.initial[np.triu_indices(3)],
            gnss_weight=eta, gnss_scale_m=.015, moment_penalty='l1', starts=2, seed=71,
            perturbation_deg=30, max_nfev=3000, ftol=1e-13, gtol=1e-9)

    def test_component_rows_and_background_projection_match_independent_dense_fit(self):
        for path, receivers, blocks in [(self.project/'GFs/static/image.h5', self.image, self.blocks),
                                        (self.project/'GFs/static/gps.h5', self.gps, self.gps_blocks)]:
            with io.response_writer(path, self.sources, receivers['ids'], receivers['coordinates'], {}) as (data, done):
                data[:] = blocks
                done[:] = True
        assembled = io.assemble_joint_matrix(self.project/'GFs/static/image.h5', self.project/'GFs/static/gps.h5',
                                             self.obs, self.sources)
        np.testing.assert_allclose(assembled, self.G, rtol=0, atol=1e-18)
        before = self.G.copy()
        G, d, w, Q, R = dc.joint_system(self.G, self.obs, gnss_weight=.3, gnss_scale_m=.015)
        np.testing.assert_allclose(w[-9:].reshape(-1, 3), np.sqrt(.3)*.015/self.gps['sigma_m'])
        B = np.vstack((self.basis, np.zeros((9, 6))))
        target = np.r_[self.image['los_m'], self.gps['displacement_m'].ravel()]
        trial = self.truth.ravel()*.8
        beta = np.linalg.lstsq(w[:, None]*B, w*(target-self.G @ trial), rcond=None)[0]
        np.testing.assert_allclose(G @ trial-d, w*(self.G @ trial+B @ beta-target), atol=1e-15)
        np.testing.assert_allclose(Q.T @ G, 0, atol=1e-16)
        np.testing.assert_array_equal(self.G, before)

    def test_synthetic_recovers_dc_sources_and_six_background_coefficients(self):
        coefficients, diagnostic = self.solve()
        np.testing.assert_allclose(coefficients.ravel(), self.truth.ravel(), rtol=2e-4, atol=2e-5)
        np.testing.assert_allclose(diagnostic['joint']['background_coefficients_m'], self.beta, atol=2e-7)
        self.assertTrue(diagnostic['starts'][diagnostic['selected_start']]['success'])
        self.assertEqual(diagnostic['joint']['background_rank'], 6)
        self.assertEqual(diagnostic['joint']['gnss_weight'], 1)

    def test_profiled_gradient_matches_independent_background_refits(self):
        G, d, weights, _, _ = dc.joint_system(self.G, self.obs, gnss_weight=3, gnss_scale_m=.015)
        p = np.array([[1.7, .2, -.1, .07], [1.1, -.1, .3, -.1]])
        bases = np.broadcast_to(self.initial, (2, 3, 3)).copy()
        value, gradient = dc._l1_objective(p, G, d, bases, 1e-5, 0, None)
        B = np.vstack((self.basis, np.zeros((9, 6))))*weights[:, None]
        target = np.r_[self.image['los_m'], self.gps['displacement_m'].ravel()]
        def dense_cost(parameters):
            coefficients, _ = double_couple_tensors(parameters.reshape(2, 4), bases)
            residual = weights*(self.G @ coefficients.ravel()-target)
            residual -= B @ np.linalg.lstsq(B, residual, rcond=None)[0]
            return .5*(residual @ residual+1e-5*parameters.reshape(2, 4)[:, 0].sum())
        self.assertAlmostEqual(value, dense_cost(p), places=14)
        numerical = []
        for k in range(p.size):
            step = np.zeros(p.size); step[k] = 1e-6
            numerical.append((dense_cost(p.ravel()+step)-dense_cost(p.ravel()-step))/2e-6)
        np.testing.assert_allclose(gradient, numerical, rtol=1e-6, atol=1e-11)

    def test_control_is_independent_of_gnss_values_and_bad_weights_rank_fail(self):
        c1, d1 = self.solve(eta=0)
        changed = dict(self.obs, gnss=dict(self.obs['gnss'], displacement_m=self.gps['displacement_m']+100))
        c2, d2 = self.solve(changed, eta=0)
        np.testing.assert_array_equal(c1, c2)
        self.assertEqual(d2['joint']['gnss_scaled_residual_l2_m'], 0)
        np.testing.assert_array_equal(d1['joint']['background_coefficients_m'], d2['joint']['background_coefficients_m'])
        with self.assertRaisesRegex(ValueError, 'Rank-deficient'):
            dc.joint_system(self.G, dict(self.obs, background_basis=np.ones((24, 6))), gnss_weight=1, gnss_scale_m=.015)
        with self.assertRaises(ValueError):
            dc.joint_system(self.G, self.obs, gnss_weight=-1, gnss_scale_m=.015)

    def config(self, name, eta, mode):
        text = self.template.replace('name = "example"', 'name = '+json.dumps(name))
        text = text.replace('steps = ["greens", "synthetics", "invert"]', 'steps = ["invert"]')
        text = text.replace('observations = "template.h5"', 'observations = "joint.h5"')
        text = text.replace('observation_format = "hdf5"', 'observation_format = "joint_hdf5"')
        text = text.replace('name = "elastic"', 'name = "image"\ngnss_name = "gps"')
        text = text.replace('source_model = "full_tensor"', 'source_model = "double_couple"')
        text = text.replace('solver = "lstsq"', 'solver = "lbfgsb"')
        text = text.replace('weighting = "diagonal"', 'weighting = "none"')
        text = text[:text.index('damping =')]+('damping = [1e-10]\n\n'+
            '[inversion.double_couple]\ninitial_tensor = [-1, 0, 0, 0, 0, 1]\nmoment_penalty = "l1"\n'+
            'starts = 2\nseed = 71\nmax_nfev = 3000\nftol = 1e-13\ngtol = 1e-9\n'+
            '[inversion.joint]\ngnss_weight = '+str(eta)+'\ngnss_scale_m = 0.015\n\n'+
            '[matrix]\nname = "joint"\nmode = "'+mode+'"\n')
        path = self.project/'configs'/(name+'.toml')
        path.write_text(text)
        return path

    def test_canonical_joint_run_matrix_reuse_and_reader_integrity(self):
        first = self.config('control', 0, 'assemble')
        second = self.config('pilot', 1, 'reuse')
        c = cfg.RunConfig.load(first)
        sources, obs, geometry, _ = cfg._inputs(c)
        specs = cfg._joint_response_specs(c, sources, geometry, cfg._code_identity(), {})
        for path, receivers, blocks, specification in [(c.response_path, self.image, self.blocks, specs[0][1]),
                (c.gnss_response_path, self.gps, self.gps_blocks, specs[1][1])]:
            with io.response_writer(path, sources, receivers['ids'], receivers['coordinates'], {}) as (data, done):
                data[:] = blocks
                done[:] = True
            cfg._record_cache(self.project, specification)
        with patch.object(cfg.shutil, 'which', return_value=None), \
                patch.object(geotomo, '_inversion_gfs', side_effect=AssertionError('Unexpected GF calculation')):
            geotomo.run(first)
            cache_before = (self.project/'GFs/matrices/joint.npy.json').read_bytes()
            with patch.object(io, 'assemble_joint_matrix', side_effect=AssertionError('Unexpected assembly')):
                geotomo.run(second)
        self.assertEqual(cache_before, (self.project/'GFs/matrices/joint.npy.json').read_bytes())
        result = io.read_solution(self.project, 'pilot', 0)
        self.assertLess(result['residual_rms_mm'], .001)
        self.assertEqual(result['joint_settings']['gnss_weight'], 1)
        np.testing.assert_allclose(result['joint_predictions']['background_coefficients_m'], self.beta, atol=2e-7)
        path = self.project/'output/inverse_models/statics/pilot.0000.predictions.h5'
        np.testing.assert_array_equal(io.read_gnss_predictions(path)['ids'], self.gps['ids'])
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        with plt.rc_context():
            from mudpy import view
        before = sorted(self.project.rglob('*'))
        figure, axes = view.geotomo_joint_background(self.project, 'pilot', 0)
        np.testing.assert_allclose(axes[1].collections[0].get_array(),
                                   result['predictions']['background_m']*1000)
        figure.canvas.draw()
        plt.close(figure)
        self.assertEqual(before, sorted(self.project.rglob('*')))
        with h5py.File(str(path), 'r+') as file:
            file['background_coefficients_m'][0] += .001
        with self.assertRaisesRegex(ValueError, 'background'):
            io.read_solution(self.project, 'pilot', 0)
