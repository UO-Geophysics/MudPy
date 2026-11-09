'''
Geodetic Source Tomography: distributed deformation-source inversion.

The current source basis uses full or prescribed moment tensors with static
imagery observations.

D. Melgar
Apr 4th 2017
'''



from contextlib import contextmanager


@contextmanager
def _project_lock(home, project_name):
    """Fail before mutation if another top-level geotomo run owns the project."""
    from pathlib import Path
    lock = Path(home) / project_name / '.geotomo.lock'
    with lock.open('x'):
        try:
            yield
        finally:
            lock.unlink()


def inversionGFs(home,project_name,GF_list,source_name,model_name,
        green_flag,synth_flag,dk,pmin,pmax,kmax,hot_start,ncpus,forceMT,mt,flat=True):
    """Generate tensor responses with exclusive project ownership.

    flat=True retains historical Earth flattening; use False for a homogeneous
    flat half-space benchmark. ncpus caps concurrent FK/syn subprocesses.
    """
    with _project_lock(home,project_name):
        return _inversion_gfs(home,project_name,GF_list,source_name,model_name,
            green_flag,synth_flag,dk,pmin,pmax,kmax,hot_start,ncpus,forceMT,mt,flat)


def _inversion_gfs(home,project_name,GF_list,source_name,model_name,
        green_flag,synth_flag,dk,pmin,pmax,kmax,hot_start,ncpus,forceMT,mt,flat=True,
        *, observations=None, response_file=None, response_metadata=None):
    """Compute source responses concurrently; one run owns the project at a time.

    ``flat`` retains the historical FK Earth-flattening option. Set False for
    a flat homogeneous half-space benchmark. Static GNSS is not implemented.
    """
    from numpy import atleast_1d, atleast_2d
    from pathlib import Path
    from numpy import genfromtxt, where, isfinite

    if type(ncpus) is not int or ncpus < 1:
        raise ValueError('ncpus must be a positive integer')
    project = Path(home) / project_name
    if observations is None:
        gf_file = project / 'data/station_info' / GF_list
        stations = atleast_1d(genfromtxt(gf_file, usecols=0, dtype='U'))
        flags = atleast_2d(genfromtxt(gf_file, usecols=range(1,8)))
        if not isfinite(flags).all() or not ((flags[:,2:] == 0) | (flags[:,2:] == 1)).all():
            raise ValueError('Observation flags must be finite and binary')
        if flags[:,2:6].any():
            raise NotImplementedError('geotomo currently supports scalar imagery projections only, not GNSS or waveforms')
        selected = where(flags[:,6] == 1)[0]
        coordinates = flags[:, :2]
    else:
        stations = observations['ids']
        coordinates = observations['coordinates']
        selected = range(len(stations))
    if not len(selected):
        raise ValueError('No imagery observations requested')
    if len(set(stations)) != len(stations):
        raise ValueError('Observation identifiers must be unique')
    station_file = project / 'data/station_info/temp.sta'
    # Exclusive creation makes concurrent top-level runs fail before writing inputs.
    with station_file.open('x') as stream:
        for k in selected:
            stream.write('%s %.17g %.17g\n' % (stations[k], *coordinates[k]))
    try:
        if green_flag:
            make_parallel_green(home,project_name,station_file.name,source_name,
                model_name,1.0,64,hot_start,dk,pmin,pmax,kmax,ncpus,
                insar=True,flat=flat)
        if synth_flag:
            make_parallel_synthetics(home,project_name,station_file.name,source_name,
                model_name,hot_start,ncpus,forceMT,mt,insar=True,
                response_file=response_file, response_metadata=response_metadata)
    finally:
        station_file.unlink()


def make_parallel_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,
        hot_start,dk,pmin,pmax,kmax,ncpus,insar=False,flat=True):
    """Generate complete DC+EX static GFs with independent source workers.

    Direct callers must own the project exclusively; inversionGFs enforces this.
    """
    from pathlib import Path
    from numpy import loadtxt, savetxt, isfinite, unique
    from mudpy.parallel import run_parallel_green

    if type(ncpus) is not int or ncpus < 1:
        raise ValueError('ncpus must be a positive integer')
    project = Path(home) / project_name
    source = loadtxt(project / 'data/model_info' / fault_name, ndmin=2)
    if source.shape[1] != 4 or not isfinite(source).all():
        raise ValueError('Expected finite source rows: ID, longitude, latitude, depth_km')
    ids = source[:,0]
    if (ids <= 0).any() or (ids != ids.astype(int)).any() or len(unique(ids)) != len(ids):
        raise ValueError('Source IDs must be positive unique integers')
    if (source[:,3] <= 0).any() or not 0 <= hot_start < len(source):
        raise ValueError('Sources must be buried and hot_start must select a source')
    savetxt(project / 'data/model_info/mpi_source.0.fault', source[hot_start:],
            fmt=['%d','%.17g','%.17g','%.17g'])
    run_parallel_green(home,project_name,str(project / 'data/station_info' / station_file),
        model_name,dt,NFFT,1,dk,pmin,pmax,kmax,False,insar,0,1,
        full_mt=True,flat=flat,processes=ncpus)


def make_parallel_synthetics(home,project_name,station_file,source_name,model_name,
        hot_start,ncpus,forceMT,mt,insar=False,*,response_file=None,response_metadata=None):
    """Generate source responses concurrently and assemble with one writer.

    Direct callers must own the project exclusively; inversionGFs enforces this.
    """
    from pathlib import Path
    from numpy import loadtxt, savetxt
    from mudpy.parallel import run_parallel_synthetics_geotomo

    if type(ncpus) is not int or ncpus < 1:
        raise ValueError('ncpus must be a positive integer')
    project = Path(home) / project_name
    source = loadtxt(project / 'data/model_info' / source_name, ndmin=2)
    if not 0 <= hot_start < len(source):
        raise ValueError('hot_start must select a source')
    savetxt(project / 'data/model_info/mpi_source.0.fault', source[hot_start:],
            fmt=['%d','%.17g','%.17g','%.17g'])
    run_parallel_synthetics_geotomo(home,project_name,
        str(project / 'data/station_info' / station_file),model_name,forceMT,mt,insar,0,1,
        response_file=response_file, response_metadata=response_metadata,processes=ncpus)
    if response_file is not None:
        return
    suffix = 'los' if insar else 'neu'
    components = ('Mxx',) if forceMT else ('Mxx','Mxy','Mxz','Myy','Myz','Mzz')
    for component in components:
        folder = project / 'GFs/static'
        with (folder / ('_'+component+'.'+suffix)).open('w') as output:
            output.write('# station source_id lon_sta lat_sta lon_src lat_src depth_km '+
                         ('projection_m' if insar else 'north_m east_m up_m')+'\n')
            output.write((folder / ('tmp_'+component+'_process0')).read_text())


def getG(home,project_name,source_name,model_name,GF_list,G_from_file,G_name,forceMT):
    """Assemble scalar imagery responses in observation and source-file order.

    The legacy .npy cache has no provenance validation. Use G_from_file=False
    after any source, geometry, Earth-model, or response change. Direct callers
    must own the project exclusively; run_inversion enforces this.
    """
    from numpy import atleast_1d, atleast_2d
    from pathlib import Path
    from numpy import genfromtxt, loadtxt, load, save, isfinite

    project = Path(home) / project_name
    gf_file = project / 'data/station_info' / GF_list
    stations = atleast_1d(genfromtxt(gf_file,usecols=0,dtype='U'))
    flags = atleast_2d(genfromtxt(gf_file,usecols=range(1,8)))
    if not isfinite(flags).all() or not ((flags[:,2:] == 0) | (flags[:,2:] == 1)).all():
        raise ValueError('Observation flags must be finite and binary')
    if flags[:,2:6].any():
        raise NotImplementedError('geotomo supports scalar imagery projections only')
    selected = flags[:,6] == 1
    if not selected.any() or len(set(stations)) != len(stations):
        raise ValueError('Expected imagery observations with unique identifiers')
    sources = loadtxt(project / 'data/model_info' / source_name,ndmin=2)
    expected = (int(selected.sum()),len(sources)*(1 if forceMT else 6))
    matrix_file = project / 'GFs/matrices' / (G_name+'.npy')
    if G_from_file:
        G = load(matrix_file)
    else:
        mini = project / 'data/station_info/tempG.sta'
        with mini.open('x') as stream:
            for sta, row in zip(stations[selected],flags[selected]):
                stream.write('%s %.9f %.9f unused\n' % (sta,row[0],row[1]))
        try:
            G = makeG(home,project_name,source_name,model_name,mini.name,'insar',forceMT)
        finally:
            mini.unlink()
        save(matrix_file,G)
    if G.shape != expected or not isfinite(G).all():
        raise ValueError('Response matrix has wrong shape or nonfinite entries')
    return G


def makeG(home,project_name,source_name,model_name,station_file,gftype,forceMT):
    """Assemble complete response rows, preserving explicit source-file order."""
    from numpy import atleast_1d, atleast_2d
    from pathlib import Path
    from numpy import genfromtxt, loadtxt, zeros, isfinite, unique

    if gftype.lower() != 'insar':
        raise NotImplementedError('GNSS tensor assembly is unfinished; zero response rows are not permitted')
    project = Path(home) / project_name
    source = loadtxt(project / 'data/model_info' / source_name,ndmin=2)
    if source.shape[1] != 4 or not isfinite(source).all():
        raise ValueError('Expected finite four-column source rows')
    source_ids = source[:,0]
    if (source_ids != source_ids.astype(int)).any() or len(unique(source_ids)) != len(source_ids):
        raise ValueError('Source IDs must be unique integers')
    stations = atleast_1d(genfromtxt(project / 'data/station_info' / station_file,
                         dtype='U',usecols=0))
    if len(set(stations)) != len(stations):
        raise ValueError('Observation identifiers must be unique')
    coordinates = loadtxt(project / 'data/station_info' / station_file,usecols=(1,2),ndmin=2)
    components = ('Mxx',) if forceMT else ('Mxx','Mxy','Mxz','Myy','Myz','Mzz')
    G = zeros((len(stations),len(source)*len(components)))
    for component_index, component in enumerate(components):
        path = project / 'GFs/static' / ('_'+component+'.los')
        names = atleast_1d(genfromtxt(path,usecols=0,dtype='U'))
        values = loadtxt(path,usecols=range(1,8),ndmin=2)
        if not isfinite(values).all():
            raise ValueError('Nonfinite response in '+str(path))
        for row, station in enumerate(stations):
            matches = values[names == station]
            if len(matches) != len(source) or set(matches[:,0]) != set(source_ids):
                raise ValueError('Missing, duplicate, or unknown source responses for '+station+' in '+str(path))
            if not (abs(matches[:,1:3]-coordinates[row]) < 1e-6).all():
                raise ValueError('Response receiver coordinates do not match observation ID')
            by_id = {item[0]:item for item in matches}
            ordered = [by_id[source_id] for source_id in source_ids]
            # Responses contain source longitude, latitude, depth after station coordinates.
            for source_row, response in zip(source,ordered):
                if not all(abs(response[3:6]-source_row[1:4]) < 1e-6):
                    raise ValueError('Response source coordinates do not match source ID')
            G[row,component_index::len(components)] = [item[6] for item in ordered]
    return G


def get_moment(home,project_name,source_name,sol,forceMT,mt=None):
    """Return nonnegative tensor moments in N m and aggregate Mw.

    Coefficients scale the Mw-5 reference component amplitude. A zero total
    moment has Mw=-inf. A fixed tensor's norm is required for normalization.
    """
    from numpy import asarray, sqrt, log10, isfinite

    coefficients = asarray(sol,dtype=float).reshape(-1)
    if forceMT:
        mechanism = asarray(mt,dtype=float)
        if mechanism.shape != (6,) or not isfinite(mechanism).all():
            raise ValueError('Fixed tensor requires six finite components')
        tensors = coefficients[:,None]*mechanism
    else:
        if len(coefficients) % 6:
            raise ValueError('Full tensor model must have six components per source')
        tensors = coefficients.reshape(-1,6)
    if not isfinite(tensors).all():
        raise ValueError('Model contains nonfinite tensor components')
    M0 = sqrt((tensors*tensors*[1,2,2,1,2,1]).sum(axis=1)/2)*10**16.6
    total = M0.sum()
    Mw = (2./3)*(log10(total)-9.1) if total > 0 else -float('inf')
    return total,M0,Mw


def write_log(home,project_name,run_name,k,lambda_spatial,
        L2,Lm,VR,Mo,Mw,velmod,source,g_name,gflist,provenance=None):
    '''
    Write inversion sumamry to .log file
    
    IN:
        home: Home directory location
        project_name: Name of the problem
        run_name: Name of inversion run
        k: Inversion run number
        rupture_speed: Fastest rupture speed allowed
        num_windows: Number of temporal rupture windows
        lambda_spatial: Spatial regularization parameter
        lambda_temporal: Temporal regularization parameter
        beta: Angular offset applied to rake
        L2: L2 norm of ivnersion L2=||Gm-d||
        Lm: Model norm Lm=||L*m||
        VR: Variance reduction
        ABIC: Value of Akaike's Bayesian ifnormation criterion
        Mo: Moment in N-m
        Mw: Moment magnitude
        velmod: Earth structure model used
        fault: Fault model used
        g_name: GF matrix used
        gflist: GF control file sued
        solver: Type of solver used
    OUT:
        Nothing
    '''
    
    num=str(k).rjust(4,'0')
    f=open(home+project_name+'/output/inverse_models/models/'+run_name+'.'+num+'.log','w')
    f.write('Project: '+project_name+'\n')
    f.write('Run name: '+run_name+'\n')
    f.write('Run number: '+num+'\n')
    f.write('Source model: '+source+'\n')
    f.write('G name: '+g_name+'\n')
    f.write('GF list: '+gflist+'\n')
    f.write('lambda_spatial = '+repr(lambda_spatial)+'\n')
    f.write('L2 = '+repr(L2)+'\n')
    f.write('VR static(%) = '+repr(VR[0])+'\n')
    f.write('VR InSAR LOS(%) = '+repr(VR[4])+'\n')
    f.write('Lm = '+repr(Lm)+'\n')
    f.write('M0(N-m) = '+repr(Mo)+'\n')
    f.write('Mw = '+repr(Mw)+'\n')
    if provenance is not None:
        for key,value in provenance.items():
            f.write(str(key)+': '+str(value)+'\n')
    f.close()



def write_model(home,project_name,run_name,source_name,model_name,sol,num,forceMT,mt):
    """Write tensors in N m; nodal planes describe the best fitting DC only.

    Isotropic and zero tensors have undefined nodal planes, written as NaNs.
    """
    from pathlib import Path
    from numpy import loadtxt, asarray, full, c_, savetxt
    from mudpy.analysis import MT

    project = Path(home) / project_name
    source = loadtxt(project / 'data/model_info' / source_name,ndmin=2)
    coefficients = asarray(sol).reshape(-1)
    tensors = (coefficients[:,None]*asarray(mt) if forceMT else coefficients.reshape(-1,6))
    if len(tensors) != len(source):
        raise ValueError('Model and source counts differ')
    _, moments, _ = get_moment(home,project_name,source_name,sol,forceMT,mt)
    tensors = tensors*10**16.6
    planes = full((len(source),6),float('nan'))
    for k, tensor in enumerate(tensors):
        mechanism = MT(*tensor,*source[k,1:4],mt_style='xyz')
        mechanism.get_nodal_planes()
        planes[k,:3] = mechanism.nodal_plane1
        planes[k,3:] = mechanism.nodal_plane2
    output = project / 'output/inverse_models/models' / (run_name+'.'+str(num).zfill(4)+'.inv')
    savetxt(output,c_[source,tensors,planes,moments],
            fmt=['%d']+['%.10e']*16,
            header='No lon lat depth_km Mxx_Nm Mxy_Nm Mxz_Nm Myy_Nm Myz_Nm Mzz_Nm '+
                   'DC_strike1 DC_dip1 DC_rake1 DC_strike2 DC_dip2 DC_rake2 M0_Nm')


def run_inversion(home,project_name,run_name,source_name,model_name,GF_list,G_from_file,G_name,
        reg_spatial,nsources,solver,forceMT,mt,weight=False,Ltype=0):
    """Solve the weighted, identity-damped objective with exclusive ownership."""
    with _project_lock(home,project_name):
        return _run_inversion(home,project_name,run_name,source_name,model_name,
            GF_list,G_from_file,G_name,reg_spatial,nsources,solver,forceMT,mt,weight,Ltype)


def _run_inversion(home,project_name,run_name,source_name,model_name,GF_list,G_from_file,G_name,
        reg_spatial,nsources,solver,forceMT,mt,weight=False,Ltype=0,*,provenance=None,
        observations=None,matrix=None):
    """Minimize ||W(Gm-d)||^2 + lambda^2 ||m||^2 in data space.

    Identity damping is the only supported penalty. NNLS imposes componentwise
    nonnegativity only when explicitly requested; signed tensors use lstsq.
    """
    from mudpy import inverse as inv
    from numpy import asarray, isfinite, eye, vstack, zeros, loadtxt, full
    from numpy.linalg import lstsq, norm
    from scipy.optimize import nnls

    if Ltype != 0:
        raise NotImplementedError('Only identity damping is implemented')
    if solver not in ('lstsq','nnls'):
        raise ValueError('Unknown solver: '+str(solver))
    regularization = asarray(reg_spatial,dtype=float).reshape(-1)
    if not len(regularization) or not isfinite(regularization).all() or (regularization < 0).any():
        raise ValueError('Damping values must be finite and nonnegative')
    G = (getG(home,project_name,source_name,model_name,GF_list,G_from_file,G_name,forceMT)
         if matrix is None else matrix)
    if G.shape != (G.shape[0], nsources*(1 if forceMT else 6)) or not isfinite(G).all():
        raise ValueError('Matrix columns must match the finite source basis')
    source = loadtxt(home+project_name+'/data/model_info/'+source_name,ndmin=2)
    if len(source) != nsources:
        raise ValueError('nsources does not match source file')
    d = (inv.getdata(home,project_name,GF_list,decimate=None,bandpass=[None,None,None],quiet=True)
         if observations is None else observations['los_m'][:, None])
    if d.shape != (G.shape[0],1) or not isfinite(d).all():
        raise ValueError('Data and matrix rows differ or observations are nonfinite')
    if weight:
        if observations is not None and observations['sigma_m'] is None:
            raise ValueError('Diagonal weighting requires calibrated observation sigma_m')
        w = (inv.get_data_weights(home,project_name,GF_list,d,None)
             if observations is None else 1/observations['sigma_m'])
        if w.shape != (len(d),) or not isfinite(w).all() or (w <= 0).any():
            raise ValueError('Observation standard deviations must be positive and finite')
        WG,wd = G*w[:,None],d*w[:,None]
    else:
        WG,wd = G,d
    for k, damping in enumerate(regularization):
        A,b = WG,wd
        if damping:
            A = vstack((WG,damping*eye(G.shape[1])))
            b = vstack((wd,zeros((G.shape[1],1))))
        sol = (lstsq(A,b,rcond=None)[0] if solver == 'lstsq' else nnls(A,b[:,0])[0][:,None])
        ds = G.dot(sol)
        L2,Lmodel = norm(WG.dot(sol)-wd),norm(sol)
        if observations is None:
            VR,_ = inv.get_VR(home,project_name,GF_list,sol,d,ds,None,WG,wd)
        else:
            VR = full(5, float('nan'))
            # Retain MudPy's imagery statistic: fractional reduction in L1 residual.
            VR[4] = 100*(1-abs(d-ds).sum()/abs(d).sum()) if abs(d).sum() else float('nan')
        total,_,magnitude = get_moment(home,project_name,source_name,sol,forceMT,mt)
        write_log(home,project_name,run_name,k,damping,L2,Lmodel,VR,total,magnitude,
                  model_name,source_name,G_name,GF_list,provenance=provenance)
        if observations is None:
            inv.write_synthetics(home,project_name,run_name,GF_list,G,sol,ds,k,None)
        else:
            from pathlib import Path
            from mudpy.geotomo_io import write_predictions
            filename = Path(home)/project_name/'output/inverse_models/statics'/(
                run_name+'.'+str(k).zfill(4)+'.predictions.h5')
            write_predictions(filename, observations, ds, {'run_name': run_name,
                              'damping': float(damping), 'provenance': provenance})
        write_model(home,project_name,run_name,source_name,model_name,sol,k,forceMT,mt)


def init_project(project_directory):
    """Create a new standard MudPy project with tensor templates and TOML."""
    from mudpy.geotomo_config import init_project as initialize
    return initialize(project_directory)


def run(config_file):
    """Validate and execute a named tensor experiment from TOML."""
    from mudpy.geotomo_config import run as execute
    return execute(config_file)


def import_observations(preprocessing_directory, filename, *, id_prefix='pixel'):
    """Import a geographic preprocessing package into a new HDF5 observation file."""
    from mudpy.geotomo_io import import_preprocessing
    return import_preprocessing(preprocessing_directory, filename, id_prefix=id_prefix)
