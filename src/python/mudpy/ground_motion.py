


def ground_motion_fq(ruptfile,Mw,lon_sta,lat_sta,vs30file=None,do_pgv=True,distance_saturation=6.0):
    '''
    # Script to compute ground motions based on perimeter of fakequake rupture as it evolves in time (based on rise times).
    # Christine J. Ruhl, May 11, 2017
    # Edited August 4, 2017
    '''
    
    import numpy as np
    from mudpy import ruptfunctions as rf
    from pyproj import Geod
    import os
        
    # set projection
    g=Geod(ellps='WGS84') 
    
    # set stations        
    slon=lon_sta
    slat=lat_sta
        
    # set vs30      
    print('Calculating Vs30 for '+str(len(slon))+' stations...')
    vs30_calc=get_vs30(slon,slat,vs30file)

    
    #load fault
    fault=np.genfromtxt(ruptfile)
    
    #disregard faults with slip==0
    rise_time=fault[:,7]
    i=np.where(rise_time>0)[0]
    fault=fault[i,:]
      
    Rjb=rf.calc_rjb(slon,slat,fault)
    
    if distance_saturation!=None:
        i=np.where(Rjb<distance_saturation)[0]
        Rjb[i]=0
        
    pga=np.nan*slon.copy()
    pga_sigma=np.nan*slon.copy()
 
    print('Estimating ground motions for '+str(len(slon))+' stations...')
    for i in range(len(slon)):
        print(i)
        pga[i],pga_sigma[i] = rf.bssa14_scalar(Mw,Rjb[i],vs30_calc[i],U=0,RS=1,NS=0,SS=0,Z1=None,intensity_measure='PGA') # pga in g   
    

    
    return slon,slat,pga
    
    
    
    
def get_vs30(slon,slat,vs30file):
    '''
    Function to get vs30 value from a gmt .grd grid file and station lat and long.
    Inputs:
        slon, slat = station lon and lat
        vs30file = vs30.grd file
    Outputs:
        vs30 (scalar or array depending on slon and slat)
    '''

    from netCDF4 import Dataset
    from numpy import meshgrid,genfromtxt,unravel_index,ones
    import matplotlib.pyplot as plt

    grd = Dataset(vs30file, 'r', format='NETCDF4')

    x=grd.variables['lon'][:]
    y=grd.variables['lat'][:]
    vs30=grd.variables['z'][:]
    X,Y=meshgrid(x,y)
    
    vs30_out=760*ones(len(slon))
    
    for k in range(len(slon)):
        d=((slon[k]-X)**2+(slat[k]-Y)**2)**0.5
        r,c=unravel_index(d.argmin(), d.shape)
        vs30_out[k]=vs30[r,c]
        
    return vs30_out

        
        
        
