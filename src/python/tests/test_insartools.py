"""Scientific invariants for native-grid InSAR preprocessing."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from mudpy import insartools as insar


def geometry(shape):
    y, x = np.indices(shape, dtype=float)
    east = -0.4-0.2*x/max(1, shape[1]-1)
    north = np.full(shape, -0.1)
    return np.stack([east, north, np.sqrt(1-east**2-north**2)], axis=2)


def process(data, **options):
    defaults = dict(units='m', los_positive='toward', look_direction='ground_to_satellite',
                    coordinate_system='local EN meters', min_size=1, max_size=8,
                    threshold=0.01)
    defaults.update(options)
    return insar.preprocess_insar(data, np.arange(data.shape[1]), np.arange(data.shape[0]),
                                  geometry(data.shape), **defaults)


class PolygonTests(unittest.TestCase):
    def test_union_orientation_boundaries_and_holes(self):
        x, y = np.meshgrid(np.arange(7.), np.arange(6.))
        outer = [[0, 0], [5, 0], [5, 5], [0, 5], [0, 0]]
        hole = [[1, 1], [4, 1], [4, 4], [1, 4], [1, 1]]
        polygon = {'type': 'Polygon', 'coordinates': [outer, hole]}
        mask = insar.polygon_mask(x, y, polygon)
        self.assertFalse(mask[2, 2])
        self.assertTrue(mask[1, 2])  # Hole boundary belongs to the polygon.
        self.assertTrue(mask[0, 0])
        self.assertFalse(mask[0, 6])
        np.testing.assert_array_equal(insar.polygon_mask(x, y, outer),
                                      insar.polygon_mask(x, y, outer[::-1]))
        union = {'type': 'FeatureCollection', 'features': [
            {'type': 'Feature', 'geometry': polygon},
            {'type': 'Feature', 'geometry': {'type': 'MultiPolygon', 'coordinates': [
                [[[5, 0], [6, 0], [6, 1], [5, 1], [5, 0]]]]}}]}
        self.assertTrue(insar.polygon_mask(x, y, union)[0, 6])
        for invalid in ([], 'polygon', [[0, 0], [1, 1], [2, 2]],
                        [[0, 0], [2, 2], [0, 2], [2, 0]],
                        {'type': 'Polygon', 'coordinates': [outer,
                         [[4, 4], [7, 4], [7, 7], [4, 7]]]}):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                insar.polygon_mask(x, y, invalid)

    def test_outside_pixels_cannot_influence_kept_observations(self):
        data = np.full((8, 12), 1e6)
        data[2:6, 2:6] = 0.125
        polygon = [[1.5, 1.5], [5.5, 1.5], [5.5, 5.5], [1.5, 5.5]]
        result = process(data, include_polygons=polygon, min_size=2)
        sampling = result['sampling']
        np.testing.assert_allclose(sampling['los'], 0.125)
        self.assertEqual(sampling['counts'].sum(), 16)
        self.assertTrue(insar.polygon_mask(sampling['xy'][:, 0], sampling['xy'][:, 1], polygon).all())
        self.assertTrue(result['keep_mask'].ravel()[sampling['source_indices']].all())


class SurfaceTests(unittest.TestCase):
    def test_background_quadratic_recovery_preserves_excluded_signal(self):
        y, x = np.indices((33, 47), dtype=float)
        surface = 0.2+0.004*x-0.003*y+0.0002*x*y+0.0001*x*x-0.0003*y*y
        signal = np.zeros_like(x)
        signal[10:23, 15:32] = 0.4+0.01*x[10:23, 15:32]
        exclusion = [[14.5, 9.5], [31.5, 9.5], [31.5, 22.5], [14.5, 22.5]]
        # Observation ROI lies wholly inside the excluded deformation area.
        # The independent fit must still use the surrounding full scene.
        result = process(surface+signal, include_polygons=exclusion,
                         exclude_polygons=exclusion, trend='quadratic')
        self.assertFalse((result['fit_mask'] & result['keep_mask']).any())
        self.assertEqual(result['fit_mask'].sum(), surface.size-13*17)
        np.testing.assert_allclose(result['trend_m'], surface, atol=2e-14)
        np.testing.assert_allclose(result['corrected_m'], signal, atol=2e-14)

    def test_weighted_plane_uses_explicit_background_and_zero_weights(self):
        y, x = np.indices((9, 12), dtype=float)
        plane = 0.3+0.02*x-0.01*y
        data = plane.copy()
        data[0, 0] = 999
        weights = np.ones_like(x)
        weights[0, 0] = 0
        fit = insar.fit_surface(data, x+1e7, y-2e7, np.ones(x.shape, bool), weights=weights)
        np.testing.assert_allclose(fit['surface'], plane, atol=1e-14)
        self.assertFalse(fit['fit_mask'][0, 0])
        self.assertTrue(fit['weighted'])
        small_region = [[0, 0], [5, 0], [5, 8], [0, 8]]
        result = process(plane, background_polygons=small_region, trend='plane')
        self.assertTrue((result['fit_mask'][:, :6]).all())
        self.assertFalse((result['fit_mask'][:, 6:]).any())
        np.testing.assert_allclose(result['corrected_m'], 0, atol=1e-14)

    def test_bad_background_fails_without_order_reduction(self):
        data = np.zeros((5, 7))
        mask = np.zeros(data.shape, dtype=bool)
        mask[0, :] = True
        with self.assertRaisesRegex(ValueError, 'rank deficient'):
            insar.fit_surface(data, np.arange(7), np.arange(5), mask, model='quadratic')
        with self.assertRaisesRegex(ValueError, 'Too few'):
            insar.fit_surface(data, np.arange(7), np.arange(5), np.zeros_like(mask))
        for weights in (-np.ones(data.shape), np.full(data.shape, np.nan)):
            with self.assertRaises(ValueError):
                insar.fit_surface(data, np.arange(7), np.arange(5), np.ones_like(mask), weights=weights)


class QuadtreeTests(unittest.TestCase):
    def test_odd_rectangle_partition_and_missing_pixels(self):
        data = np.zeros((5, 7))
        data[2, 3] = np.nan
        result = process(data)
        sampling = result['sampling']
        indices = sampling['source_indices']
        self.assertEqual(len(indices), 34)
        self.assertEqual(len(np.unique(indices)), 34)
        np.testing.assert_array_equal(np.sort(indices), np.flatnonzero(np.isfinite(data)))
        np.testing.assert_array_equal(result['corrected_m'], data)
        self.assertTrue((sampling['bounds'][:, 1] <= 5).all())
        self.assertTrue((sampling['bounds'][:, 3] <= 7).all())

    def test_tile_means_preserve_linear_forward_projection(self):
        y, x = np.indices((8, 8), dtype=float)
        look = geometry(x.shape)
        displacement = np.stack([np.sin(x), y*y/100, np.cos(y)], axis=2)
        projected = np.sum(look*displacement, axis=2)
        sampling = insar.quadtree(np.zeros(x.shape), x, y, look, min_size=1,
                                  max_size=8, threshold=1)
        self.assertEqual(len(sampling['counts']), 1)
        np.testing.assert_allclose(insar.apply_quadtree(projected, sampling), projected.mean())
        np.testing.assert_allclose(sampling['look_enu'][0], look.mean(axis=(0, 1)))
        self.assertLess(np.linalg.norm(sampling['look_enu'][0]), 1)
        np.testing.assert_allclose(insar.apply_quadtree(displacement, sampling)[0],
                                   displacement.mean(axis=(0, 1)))
        broken = projected.copy()
        broken[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, 'Non-finite'):
            insar.apply_quadtree(broken, sampling)

    def test_range_std_minimum_size_and_coverage_are_explicit(self):
        data = np.zeros((4, 4))
        data[0, 0] = 0.9
        options = dict(min_size=1, max_size=4, threshold=0.5)
        ranged = process(data, criterion='range', **options)['sampling']
        std = process(data, criterion='std', **options)['sampling']
        self.assertGreater(len(ranged['counts']), len(std['counts']))
        stopped = process(data, min_size=4, max_size=4, threshold=0.01)['sampling']
        self.assertTrue(stopped['threshold_exceeded'][0])
        sparse = np.full((4, 4), np.nan)
        sparse[0, 0] = 0
        sampling = process(sparse, min_size=1, min_valid_fraction=0.8)['sampling']
        self.assertEqual(sampling['counts'].sum(), 1)
        with self.assertRaisesRegex(ValueError, 'No quadtree leaves'):
            process(sparse, min_size=4, min_valid_fraction=0.8)

    def test_polygon_holes_are_not_filled_and_locations_are_members(self):
        polygon = {'type': 'Polygon', 'coordinates': [
            [[0, 0], [7, 0], [7, 7], [0, 7]], [[1, 1], [6, 1], [6, 6], [1, 6]]]}
        result = process(np.zeros((8, 8)), include_polygons=polygon, min_size=2)
        sample = result['sampling']
        self.assertTrue(insar.polygon_mask(sample['xy'][:, 0], sample['xy'][:, 1], polygon).all())
        self.assertTrue(result['keep_mask'].ravel()[sample['source_indices']].all())
        self.assertEqual(sample['counts'].sum(), result['keep_mask'].sum())
        # The saved polygon is directly reusable, including hole semantics.
        record = json.loads(json.dumps(result['provenance']))
        for polygon in record['include_polygons']['coordinates']:
            for ring in polygon:
                self.assertEqual(ring[0], ring[-1])
        replay = process(np.zeros((8, 8)), include_polygons=record['include_polygons'], min_size=2)
        np.testing.assert_array_equal(replay['keep_mask'], result['keep_mask'])
        np.testing.assert_array_equal(replay['sampling']['source_indices'], sample['source_indices'])


class PipelineTests(unittest.TestCase):
    def test_units_polarity_and_inputs_are_preserved(self):
        y, x = np.indices((6, 9), dtype=float)
        data = 10*x+y
        look = -geometry(data.shape)
        data_before, look_before = data.copy(), look.copy()
        result = insar.preprocess_insar(data, x, y, look, units='mm', los_positive='away',
                   look_direction='satellite_to_ground', coordinate_system='local EN meters')
        np.testing.assert_array_equal(data, data_before)
        np.testing.assert_array_equal(look, look_before)
        np.testing.assert_allclose(result['original_m'], -data/1000)
        np.testing.assert_allclose(result['look_enu'], -look)
        self.assertEqual(result['provenance']['output_los_positive'], 'toward')
        unknown = process(data, los_positive='unknown', look_direction='unknown')
        self.assertEqual(unknown['provenance']['output_los_positive'], 'unknown')
        np.testing.assert_array_equal(unknown['original_m'], data)

    def test_explicit_sentinel_quality_and_geometry_masks(self):
        data = np.zeros((8, 8))
        data[0, 0] = -999
        data[0, 1] = np.nan
        quality = np.ones(data.shape, dtype=bool)
        quality[1, 0] = False
        result = process(data, nodata=-999, quality_mask=quality)
        sample = result['sampling']
        self.assertEqual(sample['counts'].sum(), 61)
        self.assertEqual(result['provenance']['counts']['missing_displacement_pixels'], 2)
        # Geometry outside the observation ROI is not required to fit a LOS ramp.
        look = geometry(data.shape)
        look[:, :4] = np.nan
        include = [[4, 0], [7, 0], [7, 7], [4, 7]]
        result = insar.preprocess_insar(np.zeros(data.shape), np.arange(8), np.arange(8), look,
                    units='m', los_positive='toward', look_direction='ground_to_satellite',
                    coordinate_system='local EN meters', include_polygons=include, trend='plane')
        self.assertEqual(result['fit_mask'].sum(), 64)
        self.assertEqual(result['sampling']['counts'].sum(), 32)

    def test_invalid_parameters_fail_loudly(self):
        data = np.zeros((4, 4))
        for option in (dict(units='radians'), dict(los_positive='guess'),
                       dict(look_direction='guess'), dict(coordinate_system=''),
                       dict(min_size=3), dict(max_size=2, min_size=4),
                       dict(threshold=-1), dict(min_valid_fraction=1.1),
                       dict(quality_mask=np.ones(data.shape)), dict(trend='cubic'),
                       dict(fit_weights=np.ones(data.shape))):
            with self.subTest(option=option), self.assertRaises(ValueError):
                process(data, **option)
        with self.assertRaisesRegex(ValueError, 'real-valued'):
            process(data.astype(complex)+1j)
        look = geometry(data.shape)*2
        with self.assertRaisesRegex(ValueError, 'unit-length'):
            insar.preprocess_insar(data, np.arange(4), np.arange(4), look, units='m',
                los_positive='toward', look_direction='ground_to_satellite', coordinate_system='local')

    def test_artifacts_are_reconstructable_and_exclusive(self):
        data = np.arange(35.).reshape(5, 7)/100
        metadata = {'acquisitions': ['pre', 'post']}
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary)/'source.txt'
            np.savetxt(source, data)
            result = process(data, source_files=[source], metadata=metadata)
            metadata['acquisitions'][0] = 'changed after processing'
            self.assertEqual(result['provenance']['metadata']['acquisitions'][0], 'pre')
            destination = Path(temporary)/'run'
            insar.write_preprocessing(result, destination)
            with np.load(destination/'processed.npz', allow_pickle=False) as archive:
                for name in archive.files:
                    self.assertNotEqual(archive[name].dtype.kind, 'O')
                sampling = {'source_shape': archive['source_shape'],
                            'source_indices': archive['sample_source_indices'],
                            'offsets': archive['sample_offsets'], 'counts': archive['sample_counts']}
                np.testing.assert_allclose(insar.apply_quadtree(data, sampling), archive['sample_los'])
            record = json.loads((destination/'processing.json').read_text())
            self.assertEqual(record['source_files'][0]['path'], str(source.resolve()))
            np.testing.assert_allclose(np.loadtxt(destination/'samples.txt')[:, 2], result['sampling']['los'])
            before = (destination/'processing.json').read_bytes()
            with self.assertRaises(FileExistsError):
                insar.write_preprocessing(result, destination)
            self.assertEqual(before, (destination/'processing.json').read_bytes())
            insar.plot_preprocessing(result, destination/'overview.png')
            self.assertGreater((destination/'overview.png').stat().st_size, 1000)
            with self.assertRaises(FileExistsError):
                insar.plot_preprocessing(result, destination/'overview.png')
            from matplotlib import pyplot as plt
            self.assertEqual(plt.get_fignums(), [])

    def test_weighted_fit_is_reproducible_from_saved_arrays(self):
        y, x = np.indices((6, 9), dtype=float)
        data = 0.2+0.01*x-0.02*y
        data[0, 0] = 100
        weight = np.ones(data.shape)
        weight[0, 0] = 0
        background = [[0, 0], [8, 0], [8, 5], [0, 5]]
        result = process(data, trend='plane', fit_weights=weight, background_polygons=background)
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)/'weighted'
            insar.write_preprocessing(result, destination)
            with np.load(destination/'processed.npz', allow_pickle=False) as archive:
                fit = insar.fit_surface(archive['original_m'], archive['x'], archive['y'],
                           archive['fit_mask'], weights=archive['fit_weights'])
                np.testing.assert_allclose(fit['surface'], archive['trend_m'], atol=1e-14)


if __name__ == '__main__':
    unittest.main()
