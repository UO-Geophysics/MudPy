"""Check the physical interpretation of the example recovery diagnostics."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np


path = Path(__file__).resolve().parents[3]/'examples/mt3d/validate.py'
spec = importlib.util.spec_from_file_location('mt3d_validation_example', path)
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


class RecoveryMetricTests(unittest.TestCase):
    def test_exact_source_and_signed_contraction(self):
        location = np.array([[.2, -.3, 4.]])
        tensor = np.array([[1., 0, 0, 1, 0, 1]])
        for sign in (-1, 1):
            result = example.tensor_metrics(location, sign*tensor, location, sign*tensor)
            self.assertEqual(result['centroid_error_km'], 0.)
            self.assertEqual(result['spread_km'], 0.)
            self.assertEqual(result['aggregate_tensor_relative_error'], 0.)
            self.assertAlmostEqual(result['scalar_moment_ratio'], 1.)
            self.assertLess(result['aggregate_mechanism_angle_deg'], 1e-5)
            self.assertAlmostEqual(result['total_scalar_moment_Nm']/10**16.6, np.sqrt(1.5))
        opposite = example.tensor_metrics(location, -tensor, location, tensor)
        self.assertAlmostEqual(opposite['aggregate_mechanism_angle_deg'], 180., places=5)
        self.assertEqual(opposite['aggregate_tensor_relative_error'], 2.)

    def test_neighbor_cancellation_is_visible_despite_correct_aggregate(self):
        locations = np.array([[-1., 0, 3], [1., 0, 3]])
        tensors = np.array([[2., 0, 0, 2, 0, 2], [-1., 0, 0, -1, 0, -1]])
        true_location = np.array([[-1/3, 0, 3]])
        true_tensor = np.array([[1., 0, 0, 1, 0, 1]])
        result = example.tensor_metrics(locations, tensors, true_location, true_tensor)
        self.assertEqual(result['aggregate_tensor_relative_error'], 0.)
        self.assertLess(result['centroid_error_km'], 1e-15)
        self.assertAlmostEqual(result['scalar_moment_ratio'], 3.)
        self.assertAlmostEqual(result['cancellation_ratio'], 3.)
        self.assertAlmostEqual(result['spread_km'], np.sqrt(8/9))
        self.assertEqual(result, example.tensor_metrics(locations[::-1], tensors[::-1], true_location, true_tensor))

    def test_zero_prediction_has_undefined_spatial_and_mechanism_metrics(self):
        location = np.array([[0., 0, 3]])
        truth = np.array([[0., 1, 0, 0, 0, 0]])
        result = example.tensor_metrics(location, np.zeros((1, 6)), location, truth)
        self.assertIsNone(result['centroid_km'])
        self.assertIsNone(result['spread_km'])
        self.assertIsNone(result['aggregate_mechanism_angle_deg'])
        self.assertIsNone(result['cancellation_ratio'])
        self.assertEqual(result['scalar_moment_ratio'], 0.)
        self.assertEqual(result['aggregate_tensor_relative_error'], 1.)

    def test_invalid_arrays_and_undefined_truth_fail(self):
        location = np.array([[0., 0, 3]])
        tensor = np.array([[0., 1, 0, 0, 0, 0]])
        for bad in (np.zeros((1, 5)), np.full((1, 6), np.nan)):
            with self.assertRaises(ValueError):
                example.tensor_metrics(location, bad, location, tensor)
        with self.assertRaisesRegex(ValueError, 'Truth must have nonzero'):
            example.tensor_metrics(location, tensor, location, np.zeros((1, 6)))


if __name__ == '__main__':
    unittest.main()
