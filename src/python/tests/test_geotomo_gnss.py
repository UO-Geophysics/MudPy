"""Explicit GNSS units/order, fixed-model NEU prediction and read-only viewers."""
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mudpy import geotomo_io as io
with matplotlib.rc_context():
    from mudpy import view


class GNSSTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.addCleanup(plt.close, 'all')
        self.folder = Path(self.temporary.name)
        self.table = self.folder/'offsets.dat'
        self.table.write_text('A 13 43 1200 10 2 -20 3 30 8\nB 13.1 43.1 300 -4 1 5 2 -6 7\n')

    def observations(self):
        path = io.import_gnss_table(self.table, self.folder/'gps.h5',
                                    metadata={'uncertainty_confidence': 'unknown'})
        return io.read_gnss(path)

    def test_import_preserves_units_neu_order_heights_and_explicit_selection(self):
        path = io.import_gnss_table(self.table, self.folder/'gps.h5', metadata={}, exclude={'B': 'duplicate under review'})
        data = io.read_gnss(path)
        np.testing.assert_array_equal(data['ids'], ['A'])
        np.testing.assert_allclose(data['displacement_m'], [[-.02, .01, .03]])
        np.testing.assert_allclose(data['sigma_m'], [[.003, .002, .008]])
        np.testing.assert_array_equal(data['height_m'], [1200])
        self.assertEqual(data['metadata']['excluded_stations'], {'B': 'duplicate under review'})
        self.assertEqual(len(data['metadata']['input_sha256'][str(self.table.resolve())]), 64)
        with self.assertRaises(OSError):
            io.import_gnss_table(self.table, path, metadata={})
        with self.assertRaisesRegex(ValueError, 'existing station IDs'):
            io.import_gnss_table(self.table, self.folder/'bad.h5', metadata={}, exclude={'C': 'unknown'})
        self.assertFalse((self.folder/'bad.h5').exists())

    def test_invalid_station_table_and_changed_hdf5_conventions_fail(self):
        original = self.table.read_text()
        for text in (original.replace('B ', 'A '), original.replace('10 2 ', '10 0 '),
                     original.replace('10 2 ', 'nan 2 '), original.replace('A 13 ', 'A 190 ')):
            self.table.write_text(text)
            with self.assertRaises(ValueError):
                io.import_gnss_table(self.table, self.folder/'invalid.h5', metadata={})
        self.table.write_text(original)
        self.observations()
        with h5py.File(str(self.folder/'gps.h5'), 'r+') as file:
            file.attrs['components'] = 'ENU'
        with self.assertRaisesRegex(ValueError, 'units, components'):
            io.read_gnss(self.folder/'gps.h5')

    def test_fixed_model_prediction_checks_source_and_station_order_and_reference_moment(self):
        obs = self.observations()
        sources = np.array([[9, 13, 43, 4], [2, 13.1, 43.1, 7.]])
        blocks = np.arange(72.).reshape(2, 2, 3, 6)*1e-4
        response = self.folder/'responses.h5'
        with io.response_writer(response, sources, obs['ids'], obs['coordinates'], {}) as (data, complete):
            data[:] = blocks
            complete[:] = True
        coefficient = np.array([[1, -.2, .3, -1, .5, .1], [-.4, 1, 0, .2, -.5, 2]])
        model = dict(ids=sources[:, 0], coordinates=sources[:, 1:3], depth_km=sources[:, 3],
                     tensors_nm=coefficient*io.REFERENCE_MOMENT_NM)
        expected = sum(np.array([[sum(blocks[s, r, c, t]*coefficient[s, t] for t in range(6))
                                  for c in range(3)] for r in range(2)]) for s in range(2))
        np.testing.assert_allclose(io.predict_gnss(response, obs, model), expected, atol=1e-16)
        reversed_model = {key: value[::-1] for key, value in model.items()}
        with self.assertRaisesRegex(ValueError, 'ordering'):
            io.predict_gnss(response, obs, reversed_model)
        with self.assertRaisesRegex(ValueError, 'ordering'):
            io.predict_gnss(response, dict(obs, ids=obs['ids'][::-1]), model)
        with h5py.File(str(response), 'r+') as file:
            file.attrs['complete'] = False
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            io.predict_gnss(response, obs, model)

    def test_prediction_archive_and_vector_views_keep_components_and_residual_sign(self):
        obs = self.observations()
        predicted = obs['displacement_m']*.5
        path = io.write_gnss_predictions(self.folder/'prediction.h5', obs, predicted,
                                         provenance={'model_sha256': 'test'})
        data = io.read_gnss_predictions(path)
        np.testing.assert_allclose(data['residual_m'], obs['displacement_m']*.5)
        before = sorted(self.folder.iterdir())
        figure, axes = view.geotomo_gnss_fit(path, labels=['A'], vector_scale=100.)
        self.assertEqual(axes.shape, (2, 3))
        np.testing.assert_allclose(axes[0, 0].collections[1].U, obs['displacement_m'][:, 1]*1000)
        np.testing.assert_allclose(axes[0, 0].collections[1].V, obs['displacement_m'][:, 0]*1000)
        np.testing.assert_allclose(axes[1, 2].collections[0].get_array(), obs['displacement_m'][:, 2]*500)
        figure.canvas.draw()
        figure, axes = view.geotomo_gnss_components(path)
        self.assertEqual(len(axes), 3)
        figure.canvas.draw()
        self.assertEqual(before, sorted(self.folder.iterdir()))
        with h5py.File(str(path), 'r+') as file:
            file['residual_m'][0, 0] *= -1
        with self.assertRaisesRegex(ValueError, 'residuals'):
            io.read_gnss_predictions(path)


if __name__ == '__main__':
    unittest.main()
