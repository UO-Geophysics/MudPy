"""Project setup and version-1 TOML runs for the existing mt3d workflow."""
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
                     'inversion', 'matrix', 'execution'),
              ('schema_version', 'project', 'run', 'inputs', 'forward', 'inversion'), 'root')
        if type(data['schema_version']) is not int or data['schema_version'] != 1:
            raise ValueError('Only schema_version = 1 is supported')
        project, run, inputs = data['project'], data['run'], data['inputs']
        forward, inversion = data['forward'], data['inversion']
        matrix, execution = data.get('matrix', {}), data.get('execution', {})
        _keys(project, ('directory',), ('directory',), 'project')
        _keys(run, ('name', 'steps'), ('name', 'steps'), 'run')
        _keys(inputs, ('sources', 'earth_model', 'observations'),
              ('sources', 'earth_model', 'observations'), 'inputs')
        _keys(forward, ('backend', 'earth_flattening', 'dk', 'pmin', 'pmax', 'kmax'),
              ('earth_flattening',), 'forward')
        _keys(inversion, ('source_model', 'fixed_tensor', 'solver', 'weighting', 'damping'),
              ('source_model', 'damping'), 'inversion')
        _keys(matrix, ('name', 'mode'), (), 'matrix')
        _keys(execution, ('processes',), (), 'execution')
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
        filenames = {key: _name(inputs[key], 'inputs.'+key) for key in inputs}
        _choice(forward.get('backend', 'fk'), ('fk',), 'forward.backend')
        flattening = forward['earth_flattening']
        if type(flattening) is not bool:
            raise ValueError('forward.earth_flattening must be a boolean')
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
        if type(execution.get('processes', 1)) is not int or execution.get('processes', 1) != 1:
            raise ValueError('execution.processes must be 1; parallel tensor execution is not validated')
        if 'greens' in steps and 'invert' in steps and 'synthetics' not in steps:
            raise ValueError('New greens must be followed by synthetics before inversion')
        if matrix_mode == 'reuse' and steps != ['invert']:
            raise ValueError('matrix.mode = reuse requires run.steps = ["invert"]')
        return cls(path, original, directory, name, tuple(steps), filenames['sources'],
                   filenames['earth_model'], filenames['observations'], flattening,
                   dk, pmin, pmax, kmax, source_model, mechanism, solver, weighting,
                   _damping(inversion['damping']), matrix_name, matrix_mode)

    def settings(self):
        return {
            'schema_version': 1, 'project': {'directory': str(self.project)},
            'run': {'name': self.name, 'steps': self.steps},
            'inputs': {'sources': str(self.project/'data/model_info'/self.sources),
                       'earth_model': str(self.project/'structure'/self.earth_model),
                       'observations': str(self.project/'data/station_info'/self.observations)},
            'forward': {'backend': 'fk', 'earth_flattening': self.earth_flattening,
                        'dk': self.dk, 'pmin': self.pmin, 'pmax': self.pmax, 'kmax': self.kmax},
            'inversion': {'source_model': self.source_model, 'fixed_tensor': self.fixed_tensor,
                          'solver': self.solver, 'weighting': self.weighting, 'damping': self.damping},
            'matrix': {'name': self.matrix_name, 'mode': self.matrix_mode,
                       'path': str(self.project/'GFs/matrices'/(self.matrix_name+'.npy'))},
            'execution': {'processes': 1}}


def init_project(project_directory):
    from mudpy.runslip import _populate_project
    project = Path(project_directory).expanduser().resolve()
    if any(char.isspace() for char in str(project)):
        raise ValueError('The legacy observation format requires a project path without whitespace')
    templates = Path(__file__).resolve().parents[3]/'run/mt3d'
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


def _inputs(config):
    """Validate all primary inputs and normalize only selected imagery rows."""
    project = config.project
    paths = [project/'data/model_info'/config.sources, project/'structure'/config.earth_model,
             project/'data/station_info'/config.observations]
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
            raise ValueError('Only imagery scalar projections are supported by the current mt3d runner')
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
    modules = ('mt3d', 'mt3d_config', 'parallel', 'green', 'inverse', 'analysis', 'runslip')
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
    base = {'sources': _hash(project/'data/model_info'/config.sources),
            'earth_model': _hash(project/'structure'/config.earth_model),
            'receivers': [{key: row[key] for key in ('id', 'lon', 'lat')} for row in geometry],
            'forward': config.settings()['forward'],
            'code': {name: code['source_sha256'][name] for name in ('mt3d', 'parallel', 'green')},
            'executables': {name: executables.get(name) for name in ('fk.pl', 'st_fk', 'trav')}}
    responses = {'greens': _digest(base), 'geometry': geometry,
                 'source_model': config.source_model, 'fixed_tensor': config.fixed_tensor,
                 'syn': executables.get('syn')}
    matrix = {'responses': _digest(responses), 'assembly_code': code['source_sha256']['mt3d']}
    greens_files = [project/'GFs/static'/('%s.static.%.4f.sub%04d.insar' %
                    (config.earth_model, source[3], source[0])) for source in sources]
    components = ('Mxx',) if config.source_model == 'fixed_tensor' else ('Mxx', 'Mxy', 'Mxz', 'Myy', 'Myz', 'Mzz')
    response_files = [project/'GFs/static'/('_'+component+'.los') for component in components]
    return ((project/'GFs/static/_mt3d_greens.json', base, greens_files),
            (project/'GFs/static/_mt3d_responses.json', responses, response_files),
            (project/'GFs/matrices'/(config.matrix_name+'.npy.json'), matrix,
             [project/'GFs/matrices'/(config.matrix_name+'.npy')]))


def _record_cache(project, specification):
    path, inputs, products = specification
    record = {'signature': _digest(inputs), 'inputs': inputs,
              'files': {str(product.relative_to(project)): _hash(product) for product in products}}
    path.write_bytes(_json_bytes(record))


def _require_cache(project, specification):
    path, inputs, products = specification
    if not path.is_file():
        raise ValueError('Missing cache record '+str(path)+'; explicitly recompute this stage')
    record = json.loads(path.read_text())
    if not isinstance(record, dict) or record.get('signature') != _digest(inputs):
        raise ValueError('Cached inputs/settings do not match '+str(path)+'; explicitly recompute this stage')
    expected = {str(product.relative_to(project)) for product in products}
    if not isinstance(record.get('files'), dict) or set(record['files']) != expected:
        raise ValueError('Incomplete cache product record: '+str(path))
    for product in products:
        if not product.is_file() or _hash(product) != record['files'][str(product.relative_to(project))]:
            raise ValueError('Missing or changed cached product: '+str(product))


def _collisions(config, normalized, matrix_specification):
    folder = config.project/'output/inverse_models/models'
    artifacts = [folder/(config.name+'.config.toml'), folder/(config.name+'.resolved.json'), normalized]
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
    from mudpy import mt3d
    import scipy
    import obspy
    import pyproj
    config = RunConfig.load(config_file)
    project = config.project
    required = ('data/model_info', 'data/station_info', 'data/statics', 'structure',
                'GFs/static', 'GFs/matrices', 'output/inverse_models/models',
                'output/inverse_models/statics')
    for directory in required:
        if not (project/directory).is_dir():
            raise FileNotFoundError(project/directory)
    home, name = str(project.parent)+'/', project.name
    # One lock covers validation, snapshots, every compute stage, and logs.
    with mt3d._project_lock(home, name):
        sources, rows, geometry, input_hashes = _inputs(config)
        code = _code_identity()
        executables = {}
        needed = (('fk.pl', 'st_fk', 'trav') if 'greens' in config.steps else ())
        needed += ('syn',) if 'synthetics' in config.steps else ()
        for executable in ('fk.pl', 'st_fk', 'trav', 'syn'):
            path = shutil.which(executable)
            if path is None and executable in needed:
                raise FileNotFoundError('Required executable is not on PATH: '+executable)
            executables[executable] = {'path': path, 'sha256': _hash(path)} if path else None
        greens, responses, matrix = _signatures(config, sources, geometry, code, executables)
        normalized = project/'data/station_info'/('_mt3d_'+config.name+'.gflist')
        _collisions(config, normalized, matrix)
        if 'synthetics' in config.steps and 'greens' not in config.steps:
            _require_cache(project, greens)
        if 'invert' in config.steps:
            if config.matrix_mode == 'reuse':
                _require_cache(project, matrix)
            elif 'synthetics' not in config.steps:
                _require_cache(project, responses)
        original = project/'output/inverse_models/models'/(config.name+'.config.toml')
        resolved = original.with_name(config.name+'.resolved.json')
        normalized_bytes = ('# Normalized selected imagery observations; paths are absolute.\n'+
                            ''.join(' '.join(row)+'\n' for row in rows)).encode('utf-8')
        record = {'settings': config.settings(), 'configuration_source': str(config.config_path),
                  'configuration_sha256': hashlib.sha256(config.original).hexdigest(),
                  'normalized_observations': str(normalized),
                  'normalized_observations_sha256': hashlib.sha256(normalized_bytes).hexdigest(),
                  'source_count': len(sources), 'observation_count': len(rows),
                  'source_ids': sources[:, 0].astype(int).tolist(), 'input_sha256': input_hashes,
                  'code': code, 'executables': executables,
                  'versions': {'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__,
                               'obspy': obspy.__version__, 'pyproj': pyproj.__version__},
                  'created_utc': datetime.now(timezone.utc).isoformat()}
        # Immutable records are installed before the first numerical operation.
        with original.open('xb') as stream:
            stream.write(config.original)
        with resolved.open('xb') as stream:
            stream.write(_json_bytes(record))
        with normalized.open('xb') as stream:
            stream.write(normalized_bytes)
        fixed = config.source_model == 'fixed_tensor'
        if 'greens' in config.steps or 'synthetics' in config.steps:
            mt3d._inversion_gfs(home, name, normalized.name, config.sources, config.earth_model,
                int('greens' in config.steps), int('synthetics' in config.steps), config.dk,
                config.pmin, config.pmax, config.kmax, 0, 1, fixed, config.fixed_tensor,
                config.earth_flattening)
            _verify_inputs(input_hashes)
            if 'greens' in config.steps:
                _record_cache(project, greens)
            if 'synthetics' in config.steps:
                _record_cache(project, responses)
        if 'invert' in config.steps:
            if config.matrix_mode == 'assemble':
                mt3d.getG(home, name, config.sources, config.earth_model, normalized.name,
                           False, config.matrix_name, fixed)
                _record_cache(project, matrix)
            provenance = {'Configuration': original.name, 'Configuration SHA256': _hash(original),
                          'Resolved settings': resolved.name, 'Resolved settings SHA256': _hash(resolved),
                          'Normalized observations': str(normalized),
                          'Matrix SHA256': _hash(matrix[2][0])}
            mt3d._run_inversion(home, name, config.name, config.sources, config.earth_model,
                normalized.name, True, config.matrix_name, config.damping, len(sources),
                config.solver, fixed, config.fixed_tensor, weight=config.weighting == 'diagonal',
                provenance=provenance)
        _verify_inputs(input_hashes)
        return resolved
