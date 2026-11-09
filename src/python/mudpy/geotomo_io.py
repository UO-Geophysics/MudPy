"""HDF5 observations/responses for geotomo; NumPy matrices remain separate."""
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import re

import h5py
import numpy as np


TENSOR_COMPONENTS = ('xx', 'xy', 'xz', 'yy', 'yz', 'zz')
REFERENCE_MOMENT_NM = 10**16.6


def _json(value):
    return json.dumps(value, sort_keys=True, allow_nan=False)


def _hash(filename):
    digest = hashlib.sha256()
    with Path(filename).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _header(file, kind):
    if file.attrs.get('format') != 'mudpy.geotomo.'+kind or file.attrs.get('schema_version') != 1:
        raise ValueError('Unsupported geotomo '+kind+' format/schema')


def _ids(values):
    values = np.asarray(values, dtype=str)
    if (values.ndim != 1 or not len(values) or len(set(values)) != len(values) or
            any(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value) is None for value in values)):
        raise ValueError('Expected unique nonempty observation IDs without whitespace')
    return values


def _coordinates(values, count):
    values = np.asarray(values, dtype=float)
    if (values.shape != (count, 2) or not np.isfinite(values).all() or
            (abs(values[:, 0]) > 180).any() or (abs(values[:, 1]) > 90).any()):
        raise ValueError('Expected finite longitude/latitude coordinates in degrees')
    return values


def _observations(ids, coordinates, los, projection, sigma, metadata):
    ids = _ids(ids)
    coordinates = _coordinates(coordinates, len(ids))
    los, projection = np.asarray(los, dtype=float), np.asarray(projection, dtype=float)
    if los.shape != (len(ids),) or not np.isfinite(los).all():
        raise ValueError('Expected one finite LOS displacement in metres per observation')
    if projection.shape != (len(ids), 3) or not np.isfinite(projection).all():
        raise ValueError('Expected finite NEU projection coefficients')
    norms = np.linalg.norm(projection, axis=1)
    if (norms <= 0).any() or (norms > 1+1e-5).any():
        raise ValueError('Projection coefficients must be nonzero unit vectors or unnormalized means of unit vectors')
    if sigma is not None:
        sigma = np.asarray(sigma, dtype=float)
        if sigma.shape != los.shape or not np.isfinite(sigma).all() or (sigma <= 0).any():
            raise ValueError('Expected positive finite observation standard deviations in metres')
    if not isinstance(metadata, dict):
        raise ValueError('Observation metadata must be a JSON object')
    _json(metadata)
    return {'ids': ids, 'coordinates': coordinates, 'los_m': los,
            'projection_neu': projection, 'sigma_m': sigma, 'metadata': metadata}


def write_observations(filename, ids, coordinates, los_m, projection_neu, *,
                       metadata, sigma_m=None, sampling=None):
    """Write a new scalar-imagery package. Coefficients are never normalized."""
    data = _observations(ids, coordinates, los_m, projection_neu, sigma_m, metadata)
    with h5py.File(str(filename), 'w-') as file:
        file.attrs.update(format='mudpy.geotomo.observations', schema_version=1,
                          coordinate_order='longitude,latitude', coordinate_units='degrees',
                          displacement_units='m', projection_order='NEU',
                          positive='along supplied projection',
                          forward_sampling='one representative point per observation',
                          metadata_json=_json(metadata))
        file.create_dataset('ids', data=data['ids'].astype('S'))
        for key in ('coordinates', 'los_m', 'projection_neu'):
            file.create_dataset(key, data=data[key])
        if sigma_m is not None:
            file.create_dataset('sigma_m', data=data['sigma_m'])
        if sampling is not None:
            group = file.create_group('sampling')
            for key, value in sampling.items():
                group.create_dataset(key, data=value)
    return Path(filename)


def read_observations(filename):
    """Read the downsampled observations only; native membership stays on disk."""
    with h5py.File(str(filename), 'r') as file:
        _header(file, 'observations')
        expected = {'coordinate_order': 'longitude,latitude', 'coordinate_units': 'degrees',
                    'displacement_units': 'm', 'projection_order': 'NEU',
                    'positive': 'along supplied projection',
                    'forward_sampling': 'one representative point per observation'}
        if any(file.attrs.get(key) != value for key, value in expected.items()):
            raise ValueError('Unsupported observation units, ordering, sign, or sampling')
        return _observations(file['ids'][()].astype(str), file['coordinates'][()],
                             file['los_m'][()], file['projection_neu'][()],
                             file['sigma_m'][()] if 'sigma_m' in file else None,
                             json.loads(file.attrs['metadata_json']))


def import_preprocessing(directory, filename, *, id_prefix='pixel'):
    """Import the existing preprocessing package, preserving representative pixels.

    The archive carries geographic coordinates, metre LOS, mean ENU coefficients,
    and membership. No native-pixel elastic responses are computed or stored.
    """
    directory = Path(directory).resolve()
    paths = [directory/name for name in ('samples.txt', 'processed.npz', 'processing.json')]
    hashes = {str(path): _hash(path) for path in paths}
    metadata = json.loads(paths[2].read_text())
    if (metadata['output_units'] != 'm' or metadata['output_los_positive'] != 'toward' or
            metadata['output_look_direction'] != 'ground_to_satellite'):
        raise ValueError('Import requires documented metre LOS, positive toward, and ground-to-satellite coefficients')
    system = metadata['coordinate_system'].lower()
    if not (('longitude' in system and 'latitude' in system) or system == 'epsg:4326'):
        raise ValueError('Import requires explicitly geographic longitude/latitude coordinates')
    table = np.loadtxt(paths[0], ndmin=2)
    with np.load(paths[1], allow_pickle=False) as archive:
        xy, los, look = (archive[key] for key in ('sample_xy', 'sample_los', 'sample_look_enu'))
        expected = np.column_stack((xy, los, look))
        if table.shape != expected.shape or not np.allclose(table, expected, rtol=0, atol=1e-12):
            raise ValueError('samples.txt does not match the archived observations/geometry')
        sampling = {key[len('sample_'):]: archive[key] for key in archive.files if key.startswith('sample_')}
        sampling['source_shape'] = archive['source_shape']
    offsets, counts, indices = (sampling[key] for key in ('offsets', 'counts', 'source_indices'))
    shape = sampling['source_shape']
    if (shape.shape != (2,) or (shape <= 0).any() or
            offsets.shape != (len(los)+1,) or counts.shape != (len(los),) or
            any(not np.issubdtype(a.dtype, np.integer) for a in (shape, offsets, counts, indices)) or
            indices.ndim != 1 or offsets[0] != 0 or offsets[-1] != len(indices) or
            not np.array_equal(np.diff(offsets), counts) or (counts <= 0).any() or
            (indices < 0).any() or (indices >= np.prod(shape)).any() or
            len(np.unique(indices)) != len(indices)):
        raise ValueError('Invalid archived native-pixel membership')
    if metadata['counts']['leaves'] != len(los):
        raise ValueError('Observation count disagrees with preprocessing metadata')
    if any(_hash(path) != hashes[str(path)] for path in paths):
        raise RuntimeError('Preprocessing inputs changed during import')
    return write_observations(filename, [id_prefix+str(k).zfill(6) for k in range(len(los))],
                              xy, los, look[:, [1, 0, 2]], sampling=sampling,
                              metadata={'preprocessing': metadata, 'input_sha256': hashes,
                                        'coordinate_system': metadata['coordinate_system'],
                                        'approximation': 'tile mean observed at its stored representative pixel',
                                        'uncertainty': 'not calibrated; tile scatter is not sigma'})


@contextmanager
def response_writer(filename, sources, ids, coordinates, metadata):
    """One writer installs complete source blocks; failed files remain incomplete."""
    sources = np.asarray(sources, dtype=float)
    ids = _ids(ids)
    coordinates = _coordinates(coordinates, len(ids))
    if (sources.ndim != 2 or sources.shape[1] != 4 or not len(sources) or
            not np.isfinite(sources).all() or (sources[:, 3] <= 0).any() or
            (sources[:, 0] <= 0).any() or (sources[:, 0] != sources[:, 0].astype(int)).any() or
            len(set(sources[:, 0])) != len(sources)):
        raise ValueError('Expected unique positive source IDs and finite buried source coordinates')
    _coordinates(sources[:, 1:3], len(sources))
    with h5py.File(str(filename), 'w-') as file:
        file.attrs.update(format='mudpy.geotomo.responses', schema_version=1, complete=False,
                          axes='source,receiver,NEU,tensor', tensor_coordinates='NED',
                          tensor_components=','.join(TENSOR_COMPONENTS),
                          reference_moment_Nm=REFERENCE_MOMENT_NM,
                          displacement_units='m per reference tensor coefficient',
                          source_depth_units='km', metadata_json=_json(metadata))
        file.create_dataset('sources', data=sources)
        file.create_dataset('receiver_ids', data=ids.astype('S'))
        file.create_dataset('receiver_coordinates', data=coordinates)
        responses = file.create_dataset('displacement', (len(sources), len(ids), 3, 6), dtype='f8',
                                        chunks=(1, min(len(ids), 1024), 3, 6), fillvalue=np.nan)
        complete = file.create_dataset('completed_sources', data=np.zeros(len(sources), dtype=bool))
        yield responses, complete
        if not complete[()].all():
            raise ValueError('Incomplete elastic response library')
        file.attrs['complete'] = True


def assemble_matrix(filename, observations, sources, mechanism=None):
    """Read one source block at a time and project into an unweighted matrix."""
    sources = np.asarray(sources)
    fixed = None if mechanism is None else np.asarray(mechanism, dtype=float)
    if fixed is not None and (fixed.shape != (6,) or not np.isfinite(fixed).all() or not fixed.any()):
        raise ValueError('Expected six finite nonzero-in-norm tensor components')
    count = 6 if fixed is None else 1
    matrix = np.empty((len(observations['ids']), len(sources)*count))
    with h5py.File(str(filename), 'r') as file:
        _header(file, 'responses')
        expected = {'axes': 'source,receiver,NEU,tensor', 'tensor_coordinates': 'NED',
                    'tensor_components': ','.join(TENSOR_COMPONENTS),
                    'reference_moment_Nm': REFERENCE_MOMENT_NM,
                    'displacement_units': 'm per reference tensor coefficient', 'source_depth_units': 'km'}
        if (not file.attrs.get('complete') or not file['completed_sources'][()].all() or
                any(file.attrs.get(key) != value for key, value in expected.items())):
            raise ValueError('Incomplete response library or incompatible conventions')
        if (not np.array_equal(file['receiver_ids'][()].astype(str), observations['ids']) or
                file['sources'].shape != sources.shape or
                not np.allclose(file['sources'][()], sources, rtol=0, atol=1e-9) or
                not np.allclose(file['receiver_coordinates'][()], observations['coordinates'], rtol=0, atol=1e-9)):
            raise ValueError('Response source/receiver IDs, ordering, or coordinates do not match')
        data = file['displacement']
        if data.shape != (len(sources), len(observations['ids']), 3, 6):
            raise ValueError('Invalid response-library shape')
        for k in range(len(sources)):
            block = data[k]
            if not np.isfinite(block).all():
                raise ValueError('Nonfinite elastic response block')
            projected = np.einsum('rc,rct->rt', observations['projection_neu'], block)
            matrix[:, k*count:(k+1)*count] = projected if fixed is None else projected.dot(fixed)[:, None]
    return matrix


def write_predictions(filename, observations, prediction, provenance):
    prediction = np.asarray(prediction).reshape(-1)
    if prediction.shape != observations['los_m'].shape or not np.isfinite(prediction).all():
        raise ValueError('Expected one finite predicted displacement per observation')
    write_observations(filename, observations['ids'], observations['coordinates'], observations['los_m'],
                       observations['projection_neu'], sigma_m=observations['sigma_m'],
                       metadata={'observation_metadata': observations['metadata'], 'run': provenance})
    with h5py.File(str(filename), 'r+') as file:
        file.create_dataset('prediction_m', data=prediction)
        file.create_dataset('residual_m', data=observations['los_m']-prediction)
