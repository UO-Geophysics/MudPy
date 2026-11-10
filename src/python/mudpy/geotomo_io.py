"""HDF5 observations/responses for geotomo; NumPy matrices remain separate."""
from contextlib import contextmanager, closing
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


def _response_blocks(filename, observations, sources):
    """Read finite blocks with common schema/geometry/completeness checks."""
    sources = np.asarray(sources)
    with h5py.File(str(filename), 'r') as file:
        _header(file, 'responses')
        expected = {'axes': 'source,receiver,NEU,tensor', 'tensor_coordinates': 'NED',
                    'tensor_components': ','.join(TENSOR_COMPONENTS),
                    'reference_moment_Nm': REFERENCE_MOMENT_NM,
                    'displacement_units': 'm per reference tensor coefficient', 'source_depth_units': 'km'}
        if (not file.attrs.get('complete') or file['completed_sources'].shape != (len(sources),) or
                not file['completed_sources'][()].all() or
                any(file.attrs.get(key) != value for key, value in expected.items())):
            raise ValueError('Incomplete response library or incompatible conventions')
        if (not np.array_equal(file['receiver_ids'][()].astype(str), observations['ids']) or
                file['sources'].shape != sources.shape or
                not np.allclose(file['sources'][()], sources, rtol=0, atol=1e-9) or
                not np.allclose(file['receiver_coordinates'][()], observations['coordinates'], rtol=0, atol=1e-9)):
            raise ValueError('Response source/receiver IDs, ordering, or coordinates do not match')
        data = file['displacement']
        if data.shape != (len(sources), len(observations['ids']), 3, 6) or data.dtype != np.dtype('float64'):
            raise ValueError('Invalid response-library shape')
        for k in range(len(sources)):
            block = data[k]
            if not np.isfinite(block).all():
                raise ValueError('Nonfinite elastic response block')
            yield block


def validate_responses(filename, observations, sources):
    """Verify a closed library without assembling an inversion matrix."""
    with closing(_response_blocks(filename, observations, sources)) as blocks:
        for _ in blocks:
            pass


def assemble_matrix(filename, observations, sources, mechanism=None):
    """Read one source block at a time and project into an unweighted matrix."""
    fixed = None if mechanism is None else np.asarray(mechanism, dtype=float)
    if fixed is not None and (fixed.shape != (6,) or not np.isfinite(fixed).all() or not fixed.any()):
        raise ValueError('Expected six finite nonzero-in-norm tensor components')
    count = 6 if fixed is None else 1
    matrix = np.empty((len(observations['ids']), len(sources)*count))
    with closing(_response_blocks(filename, observations, sources)) as blocks:
        for k, block in enumerate(blocks):
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


# Read-only result readers (excluded from response/matrix cache identity).


def read_initial_tensor(filename):
    """Read exactly one GCMT NDK event; ObsPy converts dyne-cm to Nm.

    Reject warnings/partial parses and multiple events or focal mechanisms.
    Return the original tensor in saved geotomo NED order, plus catalog identity.
    No event selection, coordinate guess, network access or silent fallback.
    """
    import warnings
    from obspy import read_events
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        catalog = read_events(str(filename), format='NDK')
    if len(catalog) != 1 or len(catalog[0].focal_mechanisms) != 1:
        raise ValueError('Initial catalog must contain exactly one event and focal mechanism')
    event = catalog[0]
    moment = event.focal_mechanisms[0].moment_tensor
    if moment is None or moment.tensor is None:
        raise ValueError('Initial catalog event has no moment tensor')
    t = moment.tensor
    values = np.array([t.m_tt, -t.m_tp, t.m_rt, t.m_pp, -t.m_rp, t.m_rr], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('Catalog initial tensor contains nonfinite components')
    return {'tensor_nm': values.tolist(), 'event_id': str(event.resource_id),
            'format': 'GCMT NDK', 'input_coordinates': 'USE', 'output_coordinates': 'NED',
            'units': 'Nm; NDK dyne-cm conversion performed by ObsPy'}


def read_model(filename):
    """Read a geotomo .inv model; tensors/moments are in Nm and depths in km.

    Undefined nodal-plane angles may be NaN. Scalar moment of the summed tensor
    and sum of point-source scalar moments are distinct, explicitly named values.
    """
    from mudpy.geotomo import get_moment
    table = np.loadtxt(filename, ndmin=2)
    if (table.shape[1] != 17 or not len(table) or
            not np.isfinite(table[:, :10]).all() or not np.isfinite(table[:, 16]).all() or
            np.isinf(table[:, 10:16]).any() or (table[:, 0] <= 0).any() or
            (table[:, 0] != table[:, 0].astype(int)).any() or
            len(np.unique(table[:, 0])) != len(table) or (table[:, 3] <= 0).any()):
        raise ValueError('Expected finite geotomo model rows with 17 columns and unique positive source IDs')
    _coordinates(table[:, 1:3], len(table))
    tensors = table[:, 4:10]
    total, moments, _ = get_moment('', '', '', tensors/REFERENCE_MOMENT_NM, False)
    if not np.allclose(table[:, 16], moments, rtol=1e-8, atol=0):
        raise ValueError('Stored scalar moments disagree with the model tensors')
    summed_tensor = tensors.sum(axis=0)
    net, _, _ = get_moment('', '', '', summed_tensor/REFERENCE_MOMENT_NM, False)
    return {'ids': table[:, 0].astype(int), 'coordinates': table[:, 1:3],
            'depth_km': table[:, 3], 'tensors_nm': tensors, 'nodal_planes_deg': table[:, 10:16],
            'scalar_moment_nm': moments, 'sum_scalar_moment_nm': float(total),
            'summed_tensor_nm': summed_tensor, 'net_scalar_moment_nm': float(net)}


def read_predictions(filename):
    """Read/validate observed, predicted and observed-minus-predicted samples."""
    result = read_observations(filename)
    with h5py.File(str(filename), 'r') as file:
        prediction, residual = file['prediction_m'][()], file['residual_m'][()]
    if (prediction.shape != result['los_m'].shape or residual.shape != prediction.shape or
            not np.isfinite(prediction).all() or not np.isfinite(residual).all() or
            not np.allclose(residual, result['los_m']-prediction, rtol=1e-12, atol=1e-14)):
        raise ValueError('Invalid predictions or observed-minus-predicted residuals')
    return dict(result, prediction_m=prediction, residual_m=residual)


def _result_paths(project, run_name, member):
    _ids([run_name])
    if isinstance(member, str) and member.isdecimal():
        member = int(member)
    if isinstance(member, (bool, np.bool_)) or not isinstance(member, (int, np.integer)) or member < 0:
        raise ValueError('Member must be a nonnegative integer or digit string')
    member = int(member)
    project = Path(project).expanduser().resolve()
    stem = '%s.%04d' % (run_name, member)
    return project, member, project/'output/inverse_models/models'/stem, project/'output/inverse_models/statics'/stem


def read_solution(project, run_name, member):
    """Read one canonical model/log/prediction set without modifying results."""
    project, member, model_path, prediction_path = _result_paths(project, run_name, member)
    model = read_model(str(model_path)+'.inv')
    prediction = read_predictions(str(prediction_path)+'.predictions.h5')
    fields = {}
    for line in Path(str(model_path)+'.log').read_text().splitlines():
        if ' = ' in line:
            key, value = line.split(' = ', 1)
            if key in ('lambda_spatial', 'L2', 'Lm'):
                if key in fields:
                    raise ValueError('Duplicate inversion metric: '+key)
                fields[key] = float(value)
    if (set(fields) != {'lambda_spatial', 'L2', 'Lm'} or
            any(not np.isfinite(v) or v < 0 for v in fields.values())):
        raise ValueError('Missing or invalid damping, residual norm or coefficient norm in log')
    archived = json.loads((project/'output/inverse_models/models'/(run_name+'.resolved.json')).read_text())
    settings = archived['settings']
    damping = settings['inversion']['damping']
    if settings['run']['name'] != run_name or member >= len(damping) or fields['lambda_spatial'] != damping[member]:
        raise ValueError('Log member/damping disagrees with archived run settings')
    metadata = prediction['metadata']['run']
    if metadata['run_name'] != run_name or metadata['damping'] != fields['lambda_spatial']:
        raise ValueError('Prediction identity/damping disagrees with the log')
    weighting = settings['inversion']['weighting']
    if weighting == 'diagonal':
        if prediction['sigma_m'] is None:
            raise ValueError('Weighted run predictions lack sigma_m')
        residual = prediction['residual_m']/prediction['sigma_m']
    elif weighting == 'none':
        residual = prediction['residual_m']
    else:
        raise ValueError('Unsupported archived data weighting')
    if not np.isclose(np.linalg.norm(residual), fields['L2'], rtol=1e-8, atol=1e-12):
        raise ValueError('Prediction residual norm disagrees with the log')
    solver = None
    if settings['inversion'].get('source_model') == 'double_couple':
        path = Path(str(model_path)+'.solver.json')
        if metadata['solver_sha256'] != _hash(path):
            raise ValueError('DC solver diagnostics changed after writing predictions')
        solver = json.loads(path.read_text())
        if (not solver['success'] or solver['run_name'] != run_name or solver['member'] != member or
                solver['damping'] != fields['lambda_spatial'] or solver['penalty'] != 'scalar_moment_l2' or
                solver['moment_reference_nm'] != REFERENCE_MOMENT_NM):
            raise ValueError('Invalid DC solver identity, convergence or moment normalization')
        selected = solver['starts'][solver['selected_start']]
        if (not selected['success'] or not np.isclose(selected['residual_l2'], fields['L2'], rtol=1e-8, atol=1e-12) or
                not np.isclose(selected['moment_l2'], fields['Lm'], rtol=1e-8, atol=1e-12) or
                not np.isclose(np.linalg.norm(model['scalar_moment_nm']/REFERENCE_MOMENT_NM), fields['Lm'], rtol=1e-8, atol=1e-12)):
            raise ValueError('DC selected solution disagrees with logged/model norms')
        tensor = model['tensors_nm'][:, [[0, 1, 2], [1, 3, 4], [2, 4, 5]]]
        # Check tensor shape at unit scale. The stored moment can underflow
        # to zero for a nonzero tensor when its components are squared.
        scale = np.abs(tensor).max(axis=(1, 2), keepdims=True)
        normalized = np.divide(tensor, scale, out=np.zeros_like(tensor), where=scale > 0)
        eigenvalues = np.linalg.eigvalsh(normalized)
        moment = np.linalg.norm(normalized, axis=(1, 2)) / np.sqrt(2)
        expected = moment[:, None]*np.array([-1., 0, 1])
        # Coefficients are formed in reference-moment units before conversion
        # to Nm. Subnormal coefficients can lose precision at that step.
        quantum = np.nextafter(0., 1.)*REFERENCE_MOMENT_NM
        rounding = np.divide(8*quantum, scale.ravel(), out=np.zeros_like(moment),
                             where=scale.ravel() > 0)
        if not np.all(abs(eigenvalues-expected) <= (1e-8*moment+rounding)[:, None]):
            raise ValueError('Saved DC solution contains ISO or CLVD')
    return {'member': member, 'run_name': run_name, 'model': model, 'predictions': prediction,
            'damping': fields['lambda_spatial'], 'residual_l2': fields['L2'], 'model_l2': fields['Lm'],
            'weighting': weighting, 'solver_diagnostics': solver,
            'residual_rms_mm': float(np.sqrt(np.mean(prediction['residual_m']**2))*1000)}


def read_run(project, run_name):
    """Tabulate every archived sweep member; missing/inconsistent results fail."""
    project, _, _, _ = _result_paths(project, run_name, 0)
    archived = json.loads((project/'output/inverse_models/models'/(run_name+'.resolved.json')).read_text())
    count = len(archived['settings']['inversion']['damping'])
    if not count:
        raise ValueError('Archived run has no inversion members')
    rows, first_model, first_observations = [], None, None
    for member in range(count):
        solution = read_solution(project, run_name, member)
        model, observations = solution['model'], solution['predictions']
        if first_model is None:
            first_model, first_observations = model, observations
        elif (not np.array_equal(model['ids'], first_model['ids']) or
              not np.array_equal(model['coordinates'], first_model['coordinates']) or
              not np.array_equal(model['depth_km'], first_model['depth_km']) or
              any(not np.array_equal(observations[key], first_observations[key])
                  for key in ('ids', 'coordinates', 'los_m', 'projection_neu', 'sigma_m'))):
            raise ValueError('Sweep member source or observation geometry/data changed')
        rows.append({key: solution[key] for key in ('member', 'damping', 'residual_l2', 'model_l2', 'residual_rms_mm')})
        rows[-1].update(sum_scalar_moment_nm=model['sum_scalar_moment_nm'],
                        net_scalar_moment_nm=model['net_scalar_moment_nm'])
    return {'project': project, 'run_name': run_name, 'weighting': solution['weighting'],
            'solutions': rows, 'settings': archived['settings']}
