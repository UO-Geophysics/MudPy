"""Scientific invariants and regression cases for the tensor imagery path."""
import io
from contextlib import nullcontext
import unittest
from unittest.mock import patch

import numpy as np

from mudpy import geotomo, inverse
from mudpy.analysis import MT


class TensorTests(unittest.TestCase):
    def test_vectorized_tensor_isotropic_shear_signs_units_and_linearity(self):
        azimuths = np.array([0., 90., 180., 270.])
        statics = np.arange(56., dtype=float).reshape(4, 14)/10
        statics[:, 10] = statics[:, 0]
        scale = 10**16.6*1e-15
        isotropic = geotomo._tensor_static_neu(statics, azimuths, [[1, 0, 0, 1, 0, 1]])[:, :, 0]
        np.testing.assert_allclose(isotropic[:, 2], statics[:, 11]*scale, rtol=1e-15)
        angles = np.deg2rad(azimuths)
        np.testing.assert_allclose(isotropic[:, 0], statics[:, 12]*np.cos(angles)*scale)
        np.testing.assert_allclose(isotropic[:, 1], statics[:, 12]*np.sin(angles)*scale)
        np.testing.assert_array_equal(geotomo._tensor_static_neu(statics, azimuths, [[0]*6]), np.zeros((4, 3, 1)))

        # Independent cardinal-direction checks for the off-diagonal entries.
        terms = np.zeros((2, 14))
        terms[:, [4, 5, 6, 9]] = [2, 3, 5, 7]
        basis = geotomo._tensor_static_neu(terms, [0, 90], np.eye(6))/scale
        np.testing.assert_allclose(basis[0, :, 1], [0, 7, 0], atol=1e-14)  # Mxy
        np.testing.assert_allclose(basis[1, :, 1], [7, 0, 0], atol=1e-14)
        np.testing.assert_allclose(basis[0, :, 2], [-3, 0, -2], atol=1e-14)  # Mxz
        np.testing.assert_allclose(basis[1, :, 2], [5, 0, 0], atol=1e-14)
        np.testing.assert_allclose(basis[0, :, 4], [0, 5, 0], atol=1e-14)  # Myz
        np.testing.assert_allclose(basis[1, :, 4], [0, -3, -2], atol=1e-14)

        rng = np.random.default_rng(571)
        mechanisms = rng.normal(size=(8, 6))
        all_basis = geotomo._tensor_static_neu(statics, azimuths, np.eye(6))
        combined = geotomo._tensor_static_neu(statics, azimuths, mechanisms)
        np.testing.assert_allclose(combined, np.einsum('nck,pk->ncp', all_basis, mechanisms), atol=2e-13, rtol=1e-14)
        self.assertEqual(combined.dtype, np.dtype('float64'))

    def test_vectorized_tensor_rejects_invalid_native_rows_and_shapes(self):
        statics = np.zeros((1, 14))
        for table, azimuth, mechanism in (
                (statics[:, :10], [0], np.eye(6)), (statics, [np.nan], np.eye(6)),
                (statics, [0, 1], np.eye(6)), (statics, [0], [[1]*5]),
                (np.full((1, 14), np.nan), [0], np.eye(6))):
            with self.assertRaises(ValueError):
                geotomo._tensor_static_neu(table, azimuth, mechanism)
        statics[0, 10] = 1
        with self.assertRaises(ValueError):
            geotomo._tensor_static_neu(statics, [0], np.eye(6))

    def test_moment_norm_and_signed_fixed_tensor(self):
        unit = 10**16.6
        _, moment, _ = geotomo.get_moment('', '', '', [0, 1, 0, 0, 0, 0], False)
        np.testing.assert_allclose(moment, [unit])
        _, moment, _ = geotomo.get_moment('', '', '', [-2], True, [1, 0, 0, 1, 0, 1])
        np.testing.assert_allclose(moment, [2*np.sqrt(1.5)*unit])
        self.assertEqual(geotomo.get_moment('', '', '', np.zeros(6), False)[2], -np.inf)

    def test_double_couple_and_undefined_isotropic_planes(self):
        def tensor_from_plane(plane):
            strike, dip, rake = np.radians(plane)
            n = np.array([-np.sin(dip)*np.sin(strike),
                          np.sin(dip)*np.cos(strike), -np.cos(dip)])
            slip = np.array([np.cos(rake)*np.cos(strike)+np.sin(rake)*np.cos(dip)*np.sin(strike),
                             np.cos(rake)*np.sin(strike)-np.sin(rake)*np.cos(dip)*np.cos(strike),
                             -np.sin(rake)*np.sin(dip)])
            return np.outer(n, slip)+np.outer(slip, n)

        for angles in ((0,90,0), (0,90,90), (90,90,90), (45,90,0),
                       (0,45,90), (20,30,80), (240,60,-90)):
            expected = tensor_from_plane(angles)
            components = expected[np.triu_indices(3)]
            tensor = MT(*components, 0, 0, 3, mt_style='xyz')
            tensor.get_nodal_planes()
            np.testing.assert_allclose(tensor.tensor_DC, expected, atol=1e-14)
            # Inverse trig at horizontal/vertical limits loses about sqrt(eps).
            for plane in (tensor.nodal_plane1, tensor.nodal_plane2):
                np.testing.assert_allclose(tensor_from_plane(plane), expected, atol=1e-7)
        for components in ([1, 0, 0, 1, 0, 1], np.zeros(6)):
            isotropic = MT(*components, 0, 0, 3, mt_style='xyz')
            isotropic.get_nodal_planes()
            self.assertTrue(np.isnan(isotropic.nodal_plane1).all())
            self.assertTrue(np.isnan(isotropic.nodal_plane2).all())

    def test_source_order_single_observation_and_incomplete_tables(self):
        root = '/audit/project/'
        files = {root+'data/model_info/cloud': '2 13 43 3\n1 14 44 5\n',
                 root+'data/station_info/one.sta': 'A 15 45 unused\n'}
        components = ('Mxx', 'Mxy', 'Mxz', 'Myy', 'Myz', 'Mzz')
        for k, component in enumerate(components):
            files[root+'GFs/static/_'+component+'.los'] = (
                'A 1 15 45 14 44 5 %d\nA 2 15 45 13 43 3 %d\n' % (100+k, 200+k))
        real_load, real_gen = np.loadtxt, np.genfromtxt

        def reader(real):
            def read(path, *args, **kwargs):
                return real(io.StringIO(files[str(path)]), *args, **kwargs)
            return read

        with patch('numpy.loadtxt', side_effect=reader(real_load)), \
                patch('numpy.genfromtxt', side_effect=reader(real_gen)):
            matrix = geotomo.makeG('/audit/', 'project', 'cloud', '', 'one.sta', 'insar', False)
            np.testing.assert_array_equal(matrix, [[200, 201, 202, 203, 204, 205,
                                                   100, 101, 102, 103, 104, 105]])
            np.testing.assert_array_equal(
                geotomo.makeG('/audit/', 'project', 'cloud', '', 'one.sta', 'insar', True), [[200, 100]])
            files[root+'GFs/static/_Mxx.los'] = 'A 1 15 45 14 44 5 100\n'
            with self.assertRaisesRegex(ValueError, 'Missing, duplicate'):
                geotomo.makeG('/audit/', 'project', 'cloud', '', 'one.sta', 'insar', False)
            files[root+'GFs/static/_Mxx.los'] *= 2
            with self.assertRaisesRegex(ValueError, 'Missing, duplicate'):
                geotomo.makeG('/audit/', 'project', 'cloud', '', 'one.sta', 'insar', False)
        with self.assertRaises(NotImplementedError):
            geotomo.makeG('', '', '', '', '', 'static', False)

    def test_nnls_uses_weighted_data_objective_and_augmented_damping(self):
        G = np.array([[1., 2.], [0., 1.]])
        d = np.array([[0.], [1.]])
        for weighted, damping, denominator in ((False, 0., 5.), (True, .5, 17.25)):
            with patch.object(geotomo, '_project_lock', return_value=nullcontext()), \
                    patch.object(geotomo, 'getG', return_value=G), \
                    patch('numpy.loadtxt', return_value=np.zeros((2, 4))), \
                    patch.object(inverse, 'getdata', return_value=d), \
                    patch.object(inverse, 'get_data_weights', return_value=np.array([2., 1.])), \
                    patch.object(inverse, 'get_VR', return_value=(np.zeros(5), np.zeros(5))), \
                    patch.object(inverse, 'write_synthetics'), \
                    patch.object(geotomo, 'write_log') as log, \
                    patch.object(geotomo, 'write_model') as output:
                geotomo.run_inversion('', '', '', '', '', '', False, '', [damping],
                                  2, 'nnls', True, [1, 0, 0, 1, 0, 1], weight=weighted)
                solution = output.call_args.args[5]
                np.testing.assert_allclose(solution[:, 0], [0, 1/denominator], atol=1e-14)
                w = np.array([2., 1.]) if weighted else np.ones(2)
                self.assertAlmostEqual(log.call_args.args[5], np.linalg.norm(w[:, None]*(G@solution-d)))


if __name__ == '__main__':
    unittest.main()
