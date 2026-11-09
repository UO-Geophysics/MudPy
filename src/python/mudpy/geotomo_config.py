"""Project setup and version-1 TOML runs for Geodetic Source Tomography."""
from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np


def _keys(table, allowed, required, label):
    if not isinstance(table, dict):
        raise ValueError(label+' must be a TOML table')
    unknown, missing = set(table)-set(allowed), set(required)-set(table)
    if unknown:
        raise ValueError('Unknown '+label+' keys: '+', '.join(sorted(unknown)))
    if missing:
        raise ValueError('Missing '+label+' keys: '+', '.join(sorted(missing)))


def _name(value, label):
    if not isinstance(value, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value):
        raise ValueError(label+' must start with a letter/digit and contain only letters, digits, _, - or .')
    return value


def _choice(value, choices, label):
    if not isinstance(value, str) or value not in choices:
        raise ValueError(label+' must be one of '+', '.join(choices))
    return value


def _real(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(label+' must be a finite number')
    return float(value)


def _damping(value):
    if isinstance(value, dict):
        _keys(value, ('log10_min', 'log10_max', 'count'),
              ('log10_min', 'log10_max', 'count'), 'inversion.damping')
        start = _real(value['log10_min'], 'inversion.damping.log10_min')
        stop = _real(value['log10_max'], 'inversion.damping.log10_max')
        count = value['count']
        if type(count) is not int or count < 1:
            raise ValueError('inversion.damping.count must be a positive integer')
        with np.errstate(over='raise', under='raise'):
            value = np.logspace(start, stop, count).tolist()
    if not isinstance(value, list) or not value:
        raise ValueError('inversion.damping must be a nonempty array or a log10_min/log10_max/count table')
    values = tuple(_real(item, 'inversion.damping') for item in value)
    if min(values) < 0:
        raise ValueError('Damping must be nonnegative')
    return values


@dataclass(frozen=True)
class RunConfig:
    config_path: Path
    original: bytes = field(repr=False)
    project: Path
    name: str
    steps: tuple
    sources: str
    earth_model: str
    observations: str
    earth_flattening: bool
    dk: float
    pmin: float
    pmax: float
    kmax: float
    source_model: str
    fixed_tensor: tuple
    solver: str
    weighting: str
    damping: tuple
    matrix_name: str
    matrix_mode: str
    observation_format: str
    response_name: str
    processes: int
    interface_policy: str
    interface_offset_km: float
    keep_fk_intermediates: bool

    @property
    def observation_path(self):
        folder = 'data/statics' if self.observation_format == 'hdf5' else 'data/station_info'
        return self.project/folder/self.observations

    @property
    def response_path(self):
        return self.project/'GFs/static'/(self.response_name+'.h5')

    @property
    def effective_source_path(self):
        return self.project/'data/model_info'/('_geotomo_'+self.name+'.sources')

    @classmethod
    def load(cls, filename):
        # The reader is selected by Python version, never by catching parse errors.
        if sys.version_info >= (3, 11):
            import tomllib as toml
        else:
            import tomli as toml
        path = Path(filename).expanduser().resolve()
        original = path.read_bytes()
        data = toml.loads(original.decode('utf-8'))
        _keys(data, ('schema_version', 'project', 'run', 'inputs', 'forward',
                     'inversion', 'matrix', 'execution', 'responses'),
              ('schema_version', 'project', 'run', 'inputs', 'forward', 'inversion'), 'root')
        if type(data['schema_version']) is not int or data['schema_version'] != 1:
            raise ValueError('Only schema_version = 1 is supported')
        project, run, inputs = data['project'], data['run'], data['inputs']
        forward, inversion = data['forward'], data['inversion']
        matrix, execution = data.get('matrix', {}), data.get('execution', {})
        responses = data.get('responses', {})
        _keys(project, ('directory',), ('directory',), 'project')
        _keys(run, ('name', 'steps'), ('name', 'steps'), 'run')
        _keys(inputs, ('sources', 'earth_model', 'observations', 'observation_format'),
              ('sources', 'earth_model', 'observations'), 'inputs')
        _keys(forward, ('backend', 'earth_flattening', 'dk', 'pmin', 'pmax', 'kmax',
                        'interface_policy', 'interface_offset_km'),
              ('earth_flattening',), 'forward')
        _keys(inversion, ('source_model', 'fixed_tensor', 'solver', 'weighting', 'damping'),
              ('source_model', 'damping'), 'inversion')
        _keys(matrix, ('name', 'mode'), (), 'matrix')
        _keys(execution, ('processes', 'keep_fk_intermediates'), (), 'execution')
        _keys(responses, ('name',), (), 'responses')
        directory = project['directory']
        if not isinstance(directory, str) or not directory.strip():
            raise ValueError('project.directory must be a nonempty path')
        directory = Path(directory).expanduser()
        directory = (directory if directory.is_absolute() else path.parent/directory).resolve()
        if any(char.isspace() for char in str(directory)):
            raise ValueError('The legacy whitespace-delimited observation format requires a project path without whitespace')
        name = _name(run['name'], 'run.name')
        stages = ('greens', 'synthetics', 'invert')
        steps = run['steps']
        if (not isinstance(steps, list) or not steps or
                any(not isinstance(step, str) or step not in stages for step in steps) or
                steps != [step for step in stages if step in steps]):
            raise ValueError('run.steps must be a nonempty, unique ordered subset of greens, synthetics, invert')
        filenames = {key: _name(inputs[key], 'inputs.'+key) for key in ('sources', 'earth_model', 'observations')}
        observation_format = _choice(inputs.get('observation_format', 'gflist'),
                                     ('gflist', 'hdf5'), 'inputs.observation_format')
        response_name = _name(responses.get('name', 'elastic'), 'responses.name')
        if response_name.endswith('.h5'):
            raise ValueError('responses.name is a basename without the .h5 extension')
        if responses and observation_format != 'hdf5':
            raise ValueError('responses settings require inputs.observation_format = "hdf5"')
        _choice(forward.get('backend', 'fk'), ('fk',), 'forward.backend')
        flattening = forward['earth_flattening']
        if type(flattening) is not bool:
            raise ValueError('forward.earth_flattening must be a boolean')
        interface_policy = _choice(forward.get('interface_policy', 'error'),
                                   ('error', 'above', 'below'), 'forward.interface_policy')
        interface_offset = _real(forward.get('interface_offset_km', .01), 'forward.interface_offset_km')
        if interface_offset < .0001:
            raise ValueError('forward.interface_offset_km must be at least 0.0001 km (FK depth precision)')
        dk = _real(forward.get('dk', .05), 'forward.dk')
        pmin = _real(forward.get('pmin', 0.), 'forward.pmin')
        pmax = _real(forward.get('pmax', 1.), 'forward.pmax')
        kmax = _real(forward.get('kmax', 25.), 'forward.kmax')
        if dk <= 0 or kmax <= 0 or not 0 <= pmin < pmax:
            raise ValueError('Require dk > 0, kmax > 0, and 0 <= pmin < pmax')
        source_model = _choice(inversion['source_model'], ('full_tensor', 'fixed_tensor'), 'inversion.source_model')
        mechanism = inversion.get('fixed_tensor', [])
        if not isinstance(mechanism, list):
            raise ValueError('inversion.fixed_tensor must be an array')
        mechanism = tuple(_real(item, 'inversion.fixed_tensor') for item in mechanism)
        if source_model == 'fixed_tensor' and (len(mechanism) != 6 or not any(mechanism)):
            raise ValueError('Fixed-tensor mode requires six finite, nonzero-in-norm components in xx,xy,xz,yy,yz,zz order')
        if source_model == 'full_tensor' and 'fixed_tensor' in inversion:
            raise ValueError('inversion.fixed_tensor is only valid in fixed_tensor mode')
        solver = _choice(inversion.get('solver', 'lstsq'), ('lstsq', 'nnls'), 'inversion.solver')
        weighting = _choice(inversion.get('weighting', 'diagonal'), ('diagonal', 'none'), 'inversion.weighting')
        matrix_name = _name(matrix.get('name', name), 'matrix.name')
        if matrix_name.endswith('.npy'):
            raise ValueError('matrix.name is a basename without the .npy extension')
        matrix_mode = _choice(matrix.get('mode', 'assemble'), ('assemble', 'reuse'), 'matrix.mode')
        processes = execution.get('processes', 1)
        if type(processes) is not int or processes < 1:
            raise ValueError('execution.processes must be a positive integer')
        keep_fk = execution.get('keep_fk_intermediates', True)
        if type(keep_fk) is not bool:
            raise ValueError('execution.keep_fk_intermediates must be a boolean')
        if not keep_fk and observation_format != 'hdf5':
            raise ValueError('FK intermediate cleanup requires HDF5 responses')
        if 'greens' in steps and 'invert' in steps and 'synthetics' not in steps:
            raise ValueError('New greens must be followed by synthetics before inversion')
        if matrix_mode == 'reuse' and steps != ['invert']:
            raise ValueError('matrix.mode = reuse requires run.steps = ["invert"]')
        return cls(path, original, directory, name, tuple(steps), filenames['sources'],
                   filenames['earth_model'], filenames['observations'], flattening,
                   dk, pmin, pmax, kmax, source_model, mechanism, solver, weighting,
                   _damping(inversion['damping']), matrix_name, matrix_mode,
                   observation_format, response_name, processes, interface_policy, interface_offset, keep_fk)

    def settings(self):
        settings = {
            'schema_version': 1, 'project': {'directory': str(self.project)},
            'run': {'name': self.name, 'steps': self.steps},
            'inputs': {'sources': str(self.project/'data/model_info'/self.sources),
                       'earth_model': str(self.project/'structure'/self.earth_model),
                       'observations': str(self.observation_path),
                       'observation_format': self.observation_format},
            'forward': {'backend': 'fk', 'earth_flattening': self.earth_flattening,
                        'dk': self.dk, 'pmin': self.pmin, 'pmax': self.pmax, 'kmax': self.kmax,
                        'interface_policy': self.interface_policy,
                        'interface_offset_km': self.interface_offset_km},
            'inversion': {'source_model': self.source_model, 'fixed_tensor': self.fixed_tensor,
                          'solver': self.solver, 'weighting': self.weighting, 'damping': self.damping},
            'matrix': {'name': self.matrix_name, 'mode': self.matrix_mode,
                       'path': str(self.project/'GFs/matrices'/(self.matrix_name+'.npy'))},
            'execution': {'processes': self.processes, 'keep_fk_intermediates': self.keep_fk_intermediates}}
        if self.observation_format == 'hdf5':
            settings['responses'] = {'name': self.response_name, 'path': str(self.response_path),
                                     'format': 'hdf5', 'basis': 'six reference moment tensors'}
        return settings


def init_project(project_directory):
    from mudpy.runslip import _populate_project
    project = Path(project_directory).expanduser().resolve()
    if any(char.isspace() for char in str(project)):
        raise ValueError('The legacy observation format requires a project path without whitespace')
    templates = Path(__file__).resolve().parents[3]/'run/geotomo'
    files = {
        'data/model_info/template.sources': templates/'template.sources',
        'data/station_info/template.gflist': templates/'template.gflist',
        'data/statics/template.los': templates/'template.los',
        'structure/template.mod': templates/'template.mod',
        'configs/template.toml': templates/'template.toml'}
    for source in files.values():
        if not source.is_file():
            raise FileNotFoundError(source)
    project.mkdir(parents=True)  # Existing projects are never initialized again.
    (project/'configs').mkdir()
    _populate_project(project, files)
    from mudpy.geotomo_io import write_observations
    write_observations(project/'data/statics/template.h5', ['template'], [[13., 43.]],
                       [.001], [[0., 0., 1.]], sigma_m=[.001],
                       metadata={'coordinate_system': 'Geographic longitude/latitude degrees; illustrative only',
                                 'interval': 'replace with acquisition dates', 'example_only': True})
    return project/'configs/template.toml'


def _hash(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n').encode('utf-8')


def _verify_inputs(input_hashes):
    for path, digest in input_hashes.items():
        if _hash(path) != digest:
            raise RuntimeError('Input changed during execution: '+path)


def _digest(value):
    return hashlib.sha256(_json_bytes(value)).hexdigest()


def _effective_sources(config, sources, model):
    """Resolve FK interface hits once, in physical depth coordinates.

    Match fk.pl's Vp boundary test and midpoint Earth flattening. Depths are
    quantized to the four decimals actually passed to FK. Half a depth increment
    is the interface tolerance; the chosen offset must stay in the adjacent layer.
    """
    effective = sources.copy()
    effective[:, 3] = [float('%.4f' % depth) for depth in sources[:, 3]]
    thickness, vp = model[:, 0].copy(), model[:, 2].copy()
    if config.earth_flattening:
        if thickness.sum() >= 6371 or (sources[:, 3] >= 6371).any():
            raise ValueError('FK Earth-flattening depths must be less than 6371 km')
        scale = 6371 / (6371 - np.cumsum(thickness) + .5*thickness)
        thickness *= scale
        vp *= scale
    boundaries = np.cumsum(thickness)[:-1]
    if config.earth_flattening:
        boundaries = -6371*np.expm1(-boundaries/6371)
    real = vp[:-1] != vp[1:]
    tolerance = .00005 + 1e-12
    for source in effective:
        depth = source[3]
        if depth <= 0:
            raise ValueError('Source %d rounds to a nonpositive FK depth' % source[0])
        hits = np.flatnonzero(real & (abs(boundaries-depth) <= tolerance))
        if not len(hits):
            continue
        if len(hits) != 1:
            raise ValueError('Source %d is too close to multiple FK interfaces' % source[0])
        index = hits[0]
        boundary = boundaries[index]
        label = 'Source %d at %.4f km meets FK interface %.10g km' % (source[0], depth, boundary)
        if config.interface_policy == 'error':
            raise ValueError(label+'; choose forward.interface_policy = above or below explicitly')
        above = config.interface_policy == 'above'
        source[3] = float('%.4f' % (boundary + (-1 if above else 1)*config.interface_offset_km))
        lower = (boundaries[index-1] if index else 0) if above else boundary
        upper = boundary if above else (boundaries[index+1] if index+1 < len(boundaries) else np.inf)
        if not lower+tolerance < source[3] < upper-tolerance:
            raise ValueError(label+'; interface offset must remain strictly inside the adjacent layer')
    return effective


def _inputs(config):
    """Validate all primary inputs and normalize only selected imagery rows."""
    project = config.project
    paths = [project/'data/model_info'/config.sources, project/'structure'/config.earth_model,
             config.observation_path]
    sources = np.loadtxt(paths[0], ndmin=2)
    if (sources.shape[1] != 4 or not np.isfinite(sources).all() or
            (sources[:, 0] <= 0).any() or (sources[:, 0] != sources[:, 0].astype(int)).any() or
            len(np.unique(sources[:, 0])) != len(sources) or (sources[:, 3] <= 0).any()):
        raise ValueError('Source file must have finite ID/lon/lat/depth_km rows with unique positive integer IDs and positive depth')
    if (abs(sources[:, 1]) > 180).any() or (abs(sources[:, 2]) > 90).any():
        raise ValueError('Source longitude/latitude are outside geographic bounds')
    model_text = paths[1].read_text()
    if '#' in model_text:
        raise ValueError('The FK velocity-model reader requires numeric rows without comments')
    model = np.loadtxt(paths[1], ndmin=2)
    if (not 3 <= model.shape[1] <= 6 or not np.isfinite(model).all() or
            (model[:, 0] < 0).any() or model[-1, 0] != 0 or
            (model[:, 1] <= 0).any() or (model[:, 2] <= model[:, 1]).any()):
        raise ValueError('Invalid FK model: require thickness/Vs/Vp [rho/Qs/Qp], positive speeds, Vp > Vs and a final zero-thickness half-space')
    sources = _effective_sources(config, sources, model)
    if config.observation_format == 'hdf5':
        from mudpy.geotomo_io import read_observations
        observations = read_observations(paths[2])
        if config.weighting == 'diagonal' and observations['sigma_m'] is None:
            raise ValueError('Diagonal weighting requires sigma_m; use weighting = "none" for uncalibrated data')
        geometry = [{'id': name, 'lon': xy[0], 'lat': xy[1], 'direction_neu': direction.tolist()}
                    for name, xy, direction in zip(observations['ids'], observations['coordinates'],
                                                   observations['projection_neu'])]
        return sources, observations, geometry, {str(path.resolve()): _hash(path) for path in paths}
    rows, geometry, seen = [], [], set()
    for line in paths[2].read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if not line:
            continue
        row = line.split()
        if len(row) < 28:
            raise ValueError('Observation-list rows must contain the standard 28 MudPy columns')
        station = _name(row[0], 'Observation identifier')
        if station in seen:
            raise ValueError('Duplicate observation identifier: '+station)
        seen.add(station)
        numbers = np.array(row[1:8], dtype=float)
        if (not np.isfinite(numbers).all() or abs(numbers[0]) > 180 or abs(numbers[1]) > 90 or
                not np.isin(numbers[2:], (0, 1)).all()):
            raise ValueError('Invalid observation coordinates or binary flags for '+station)
        if numbers[2:6].any():
            raise ValueError('Only imagery scalar projections are supported by the current geotomo runner')
        if numbers[6] != 1:
            continue
        observation = Path(row[12]).expanduser()
        observation = (observation if observation.is_absolute() else project/observation).resolve()
        expected = (project/'data/statics'/(station+'.los')).resolve()
        if observation != expected:
            raise ValueError('Observation '+station+' must use data/statics/'+station+'.los')
        values = np.loadtxt(observation)
        if (values.shape != (4,) or not np.isfinite(values).all() or
                not np.isclose(values[1:].dot(values[1:]), 1, rtol=0, atol=1e-6)):
            raise ValueError('Expected finite displacement_m and a unit NEU direction for '+station)
        if config.weighting == 'diagonal':
            sigma = float(row[23])
            if not math.isfinite(sigma) or sigma <= 0:
                raise ValueError('LOS uncertainty must be positive and finite for '+station)
        row[12] = str(observation)
        rows.append(row)
        geometry.append({'id': station, 'lon': numbers[0], 'lat': numbers[1],
                         'direction_neu': values[1:].tolist()})
        paths.append(observation)
    if not rows:
        raise ValueError('No imagery observations requested')
    return sources, rows, geometry, {str(path.resolve()): _hash(path) for path in paths}


def _code_identity():
    folder = Path(__file__).resolve().parent
    modules = ('geotomo', 'geotomo_config', 'geotomo_io', 'parallel', 'green', 'inverse', 'analysis', 'runslip')
    identity = {'source_sha256': {name: _hash(folder/(name+'.py')) for name in modules},
                'git_revision': None, 'git_status': None}
    repo = folder.parents[2]
    if (repo/'.git').exists():
        identity['git_revision'] = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=repo,
                                                capture_output=True, text=True, check=True).stdout.strip()
        identity['git_status'] = subprocess.run(['git', 'status', '--short'], cwd=repo,
                                              capture_output=True, text=True, check=True).stdout.splitlines()
    return identity


def _signatures(config, sources, geometry, code, executables):
    project = config.project
    forward = config.settings()['forward'].copy()
    del forward['interface_policy'], forward['interface_offset_km']
    base = {'sources': _hash(project/'data/model_info'/config.sources),
            'effective_sources': sources.tolist(),
            'earth_model': _hash(project/'structure'/config.earth_model),
            'receivers': [{key: row[key] for key in ('id', 'lon', 'lat')} for row in geometry],
            'forward': forward,
            'code': {name: code['source_sha256'][name] for name in ('geotomo', 'parallel', 'green')},
            'executables': {name: executables.get(name) for name in ('fk.pl', 'st_fk', 'trav')}}
    responses = {'greens': _digest(base), 'geometry': geometry,
                 'source_model': config.source_model, 'fixed_tensor': config.fixed_tensor,
                 'synthesis': 'numpy_mt_radiat_float64'}
    matrix = {'responses': _digest(responses), 'assembly_code': code['source_sha256']['geotomo']}
    greens_files = [project/'GFs/static'/('%s.static.%.4f.sub%04d.insar' %
                    (config.earth_model, source[3], source[0])) for source in sources]
    components = ('Mxx',) if config.source_model == 'fixed_tensor' else ('Mxx', 'Mxy', 'Mxz', 'Myy', 'Myz', 'Mzz')
    response_files = [project/'GFs/static'/('_'+component+'.los') for component in components]
    if config.observation_format == 'hdf5':
        # Physical NEU responses depend on positions, not look coefficients,
        # data, uncertainty, fixed mechanisms, or inversion choices.
        responses = {'greens': _digest(base), 'synthesis': 'numpy_mt_radiat_float64',
                     'io_code': code['source_sha256']['geotomo_io'],
                     'basis': 'xx,xy,xz,yy,yz,zz; NED; 10**16.6 Nm', 'displacement': 'NEU'}
        matrix = {'responses': _digest(responses), 'geometry': geometry,
                  'source_model': config.source_model, 'fixed_tensor': config.fixed_tensor,
                  'assembly_code': code['source_sha256']['geotomo_io']}
        return ((project/'GFs/static/_geotomo_greens.json', base, greens_files),
                (config.response_path.with_suffix('.h5.json'), responses, [config.response_path]),
                (project/'GFs/matrices'/(config.matrix_name+'.npy.json'), matrix,
                 [project/'GFs/matrices'/(config.matrix_name+'.npy')]))
    return ((project/'GFs/static/_geotomo_greens.json', base, greens_files),
            (project/'GFs/static/_geotomo_responses.json', responses, response_files),
            (project/'GFs/matrices'/(config.matrix_name+'.npy.json'), matrix,
             [project/'GFs/matrices'/(config.matrix_name+'.npy')]))


def _record_cache(project, specification):
    path, inputs, products = specification
    record = {'signature': _digest(inputs), 'inputs': inputs,
              'files': {str(product.relative_to(project)): _hash(product) for product in products}}
    path.write_bytes(_json_bytes(record))
    return record


def _require_cache(project, specification):
    path, inputs, products = specification
    if not path.is_file():
        raise ValueError('Missing cache record '+str(path)+'; explicitly recompute this stage')
    record = json.loads(path.read_text())
    if not isinstance(record, dict) or record.get('signature') != _digest(inputs):
        raise ValueError('Cached inputs/settings do not match '+str(path)+'; explicitly recompute this stage')
    availability = record.get('availability', 'available')
    if availability in ('cleanup_pending', 'pruned'):
        raise ValueError('Native FK intermediates were intentionally removed or cleanup was interrupted; '+
                         'explicitly rerun the greens stage to rebuild them: '+str(path))
    if availability != 'available':
        raise ValueError('Invalid cache availability: '+str(path))
    expected = {str(product.relative_to(project)) for product in products}
    if not isinstance(record.get('files'), dict) or set(record['files']) != expected:
        raise ValueError('Incomplete cache product record: '+str(path))
    for product in products:
        if not product.is_file() or _hash(product) != record['files'][str(product.relative_to(project))]:
            raise ValueError('Missing or changed cached product: '+str(product))
    return record


def _cleanup_fk_intermediates(config, sources, observations, greens, responses):
    """Verify the closed library and all deletion targets before pruning anything.

    The caller holds the project lock. No recursive/glob deletion is used. Cache
    state is marked pending before deletion so interruption cannot certify native
    reuse; the response library itself remains independently reusable.
    """
    from mudpy.geotomo_io import validate_responses
    project = config.project
    static = project/'GFs/static'
    static.resolve().relative_to(project)
    if static.is_symlink():
        raise ValueError('Refusing FK cleanup through a symlink: '+str(static))
    native_record = _require_cache(project, greens)
    response_record = _require_cache(project, responses)
    validate_responses(config.response_path, observations, sources)
    files = []
    for path in greens[2]:
        if path.is_symlink() or not path.is_file():
            raise ValueError('Unsafe native FK cleanup target: '+str(path))
        files.append({'path': str(path.relative_to(project)), 'size_bytes': path.stat().st_size,
                      'sha256': native_record['files'][str(path.relative_to(project))]})
    directories = []
    model_hash = _hash(project/'structure'/config.earth_model)
    for source in sources:
        folder = static/('%s_%.4f.sub%04d' % (config.earth_model, source[3], source[0]))
        if folder.is_symlink():
            raise ValueError('Unsafe FK working-directory symlink: '+str(folder))
        if not folder.exists():
            continue
        if not folder.is_dir():
            raise ValueError('Invalid FK working directory: '+str(folder))
        for path in sorted(folder.iterdir()):
            if path.name not in (config.earth_model, 'junk.p', 'junk.s') or path.is_symlink() or not path.is_file():
                raise ValueError('Unexpected FK working file; cleanup stopped before deletion: '+str(path))
            digest = _hash(path)
            if path.name == config.earth_model and digest != model_hash:
                raise ValueError('FK working model differs from project model: '+str(path))
            files.append({'path': str(path.relative_to(project)), 'size_bytes': path.stat().st_size,
                          'sha256': digest})
        directories.append(str(folder.relative_to(project)))
    receipt = {'state': 'cleanup_pending', 'run': config.name,
               'verified_response': str(config.response_path.relative_to(project)),
               'files': files, 'directories': directories,
               'size_bytes': sum(item['size_bytes'] for item in files)}
    native_record.update(availability='cleanup_pending', native_cleanup=receipt)
    response_record['native_cleanup'] = receipt
    # Both receipts exist before the first unlink; partial cleanup stays explicit.
    greens[0].write_bytes(_json_bytes(native_record))
    responses[0].write_bytes(_json_bytes(response_record))
    for item in files:
        (project/item['path']).unlink()
    for directory in directories:
        (project/directory).rmdir()
    receipt.update(state='pruned', completed_utc=datetime.now(timezone.utc).isoformat())
    native_record['availability'] = 'pruned'
    greens[0].write_bytes(_json_bytes(native_record))
    responses[0].write_bytes(_json_bytes(response_record))
    print('Removed %d native FK files and %d working folders (%d bytes); verified library: %s' %
          (len(files), len(directories), receipt['size_bytes'], config.response_path.name), flush=True)


def _collisions(config, normalized, matrix_specification):
    folder = config.project/'output/inverse_models/models'
    artifacts = [folder/(config.name+'.config.toml'), folder/(config.name+'.resolved.json'),
                 config.effective_source_path]
    if normalized is not None:
        artifacts.append(normalized)
    if config.observation_format == 'hdf5' and 'synthetics' in config.steps:
        artifacts += [config.response_path, config.response_path.with_suffix('.h5.json')]
    numbered = re.compile(re.escape(config.name)+r'\.\d{4,}\.')
    for directory in ('models', 'statics', 'waveforms'):
        path = config.project/'output/inverse_models'/directory
        if path.is_dir():
            artifacts.extend(item for item in path.iterdir() if numbered.match(item.name))
    if 'invert' in config.steps and config.matrix_mode == 'assemble':
        artifacts += [matrix_specification[0], *matrix_specification[2]]
    for path in artifacts:
        if path.exists():
            raise FileExistsError('Run or matrix output already exists: '+str(path)+'; choose a new run/matrix name or explicitly reuse a verified matrix')


def run(config_file):
    from mudpy import geotomo
    import scipy
    import obspy
    import pyproj
    config = RunConfig.load(config_file)
    hdf5 = config.observation_format == 'hdf5'
    project = config.project
    required = ('data/model_info', 'data/station_info', 'data/statics', 'structure',
                'GFs/static', 'GFs/matrices', 'output/inverse_models/models',
                'output/inverse_models/statics')
    for directory in required:
        if not (project/directory).is_dir():
            raise FileNotFoundError(project/directory)
    home, name = str(project.parent)+'/', project.name
    # One lock covers validation, snapshots, every compute stage, and logs.
    with geotomo._project_lock(home, name):
        sources, rows, geometry, input_hashes = _inputs(config)
        requested_sources = np.loadtxt(project/'data/model_info'/config.sources, ndmin=2)
        changed = np.flatnonzero(np.any(requested_sources != sources, axis=1))
        source_name = config.effective_source_path.name if len(changed) else config.sources
        code = _code_identity()
        executables = {}
        needed = (('fk.pl', 'st_fk', 'trav') if 'greens' in config.steps else ())
        for executable in ('fk.pl', 'st_fk', 'trav'):
            path = shutil.which(executable)
            if path is None and executable in needed:
                raise FileNotFoundError('Required executable is not on PATH: '+executable)
            executables[executable] = {'path': path, 'sha256': _hash(path)} if path else None
        greens, responses, matrix = _signatures(config, sources, geometry, code, executables)
        normalized = None if hdf5 else project/'data/station_info'/('_geotomo_'+config.name+'.gflist')
        _collisions(config, normalized, matrix)
        if 'synthetics' in config.steps and 'greens' not in config.steps:
            _require_cache(project, greens)
        if 'invert' in config.steps:
            if config.matrix_mode == 'reuse':
                matrix_record = _require_cache(project, matrix)
            elif 'synthetics' not in config.steps:
                _require_cache(project, responses)
        original = project/'output/inverse_models/models'/(config.name+'.config.toml')
        resolved = original.with_name(config.name+'.resolved.json')
        normalized_bytes = None if hdf5 else (
            '# Normalized selected imagery observations; paths are absolute.\n'+
            ''.join(' '.join(row)+'\n' for row in rows)).encode('utf-8')
        record = {'settings': config.settings(), 'configuration_source': str(config.config_path),
                  'configuration_sha256': hashlib.sha256(config.original).hexdigest(),
                  'normalized_observations': None if hdf5 else str(normalized),
                  'normalized_observations_sha256': None if hdf5 else hashlib.sha256(normalized_bytes).hexdigest(),
                  'source_count': len(sources), 'observation_count': len(rows['ids']) if hdf5 else len(rows),
                  'source_ids': sources[:, 0].astype(int).tolist(), 'input_sha256': input_hashes,
                  'code': code, 'executables': executables,
                  'versions': {'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__,
                               'obspy': obspy.__version__, 'pyproj': pyproj.__version__},
                  'created_utc': datetime.now(timezone.utc).isoformat(),
                  'synthesis': {'method': 'numpy_mt_radiat', 'dtype': 'float64',
                                'reference_moment_Nm': 10**16.6}}
        record['source_geometry'] = {
            'columns': ['id', 'longitude_deg', 'latitude_deg', 'depth_km'],
            'requested': requested_sources.tolist(), 'effective': sources.tolist(),
            'effective_file': str(project/'data/model_info'/source_name),
            'adjustments': [{'source_id': int(sources[k, 0]),
                             'requested_depth_km': float(requested_sources[k, 3]),
                             'effective_depth_km': float(sources[k, 3])} for k in changed]}
        if len(changed):
            with config.effective_source_path.open('x') as stream:
                np.savetxt(stream, sources, fmt=['%d', '%.17g', '%.17g', '%.4f'],
                           header='Effective FK geometry: ID longitude_deg latitude_deg depth_km')
            input_hashes[str(config.effective_source_path)] = _hash(config.effective_source_path)
            for adjustment in record['source_geometry']['adjustments']:
                print('Source %d: requested depth %.10g km -> effective FK depth %.4f km (interface policy: %s)' %
                      (adjustment['source_id'], adjustment['requested_depth_km'],
                       adjustment['effective_depth_km'], config.interface_policy), flush=True)
        if hdf5:
            import h5py
            record['versions'].update(h5py=h5py.__version__, hdf5=h5py.version.hdf5_version)
            record['observation_metadata'] = rows['metadata']
        # Immutable records are installed before the first numerical operation.
        with original.open('xb') as stream:
            stream.write(config.original)
        with resolved.open('xb') as stream:
            stream.write(_json_bytes(record))
        if normalized is not None:
            with normalized.open('xb') as stream:
                stream.write(normalized_bytes)
        fixed = config.source_model == 'fixed_tensor'
        if 'greens' in config.steps or 'synthetics' in config.steps:
            options = ({'observations': rows, 'response_file': config.response_path,
                        'response_metadata': record} if hdf5 else {})
            geotomo._inversion_gfs(home, name, config.observations if hdf5 else normalized.name,
                source_name, config.earth_model,
                int('greens' in config.steps), int('synthetics' in config.steps), config.dk,
                config.pmin, config.pmax, config.kmax, 0, config.processes, False if hdf5 else fixed,
                () if hdf5 else config.fixed_tensor, config.earth_flattening, **options)
            _verify_inputs(input_hashes)
            if 'greens' in config.steps:
                _record_cache(project, greens)
            if 'synthetics' in config.steps:
                _record_cache(project, responses)
        if 'invert' in config.steps:
            G = None
            if config.matrix_mode == 'assemble':
                if hdf5:
                    from mudpy.geotomo_io import assemble_matrix
                    G = assemble_matrix(config.response_path, rows, sources,
                                        config.fixed_tensor if fixed else None)
                    with matrix[2][0].open('xb') as stream:
                        np.save(stream, G, allow_pickle=False)
                else:
                    geotomo.getG(home, name, source_name, config.earth_model, normalized.name,
                               False, config.matrix_name, fixed)
                _verify_inputs(input_hashes)
                matrix_record = _record_cache(project, matrix)
            elif hdf5:
                # Only the .npy matrix is loaded: no response-library read or assembly.
                G = np.load(matrix[2][0], allow_pickle=False)
            provenance = {'Configuration': original.name, 'Configuration SHA256': _hash(original),
                          'Resolved settings': resolved.name, 'Resolved settings SHA256': _hash(resolved),
                          'Observations': str(config.observation_path) if hdf5 else str(normalized),
                          'Matrix SHA256': matrix_record['files'][str(matrix[2][0].relative_to(project))]}
            options = {'observations': rows, 'matrix': G} if hdf5 else {}
            geotomo._run_inversion(home, name, config.name, source_name, config.earth_model,
                config.observations if hdf5 else normalized.name,
                True, config.matrix_name, config.damping, len(sources),
                config.solver, fixed, config.fixed_tensor, weight=config.weighting == 'diagonal',
                provenance=provenance, **options)
        _verify_inputs(input_hashes)
        if hdf5 and 'synthetics' in config.steps and not config.keep_fk_intermediates:
            _cleanup_fk_intermediates(config, sources, rows, greens, responses)
        return resolved
