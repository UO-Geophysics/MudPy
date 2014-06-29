import numpy as np
from scipy.interpolate import griddata
from mudpy.inverse import d2epi
from string import rjust

epicenter=np.array([142.68,38.19,21.0]) 
f=np.genfromtxt(u'/Users/dmelgarm/Research/Slip_Inv/tohoku_10s/data/model_info/tohoku.fault')
out='/Users/dmelgarm/scripts/GMT/RTOkada/'
lon=f[:,1]
lat=f[:,2]
z=f[:,3]

loni,lati=np.mgrid[min(lon):max(lon):600j,min(lat):max(lat):600j]
zi=griddata((lon,lat),z,(loni,lati))
i,j=np.where(np.isnan(zi)!=True)

lon=loni[i,j]
lat=lati[i,j]
z=zi[i,j]

source=np.c_[lon,lat,z]
d=d2epi(epicenter,source)
t=np.arange(1,205,1)
v=1.5
for k in range(len(t)):
    i1=np.where(d<v*t[k]+1)[0]
    i2=np.where(d>v*t[k]-1)[0]
    i=np.intersect1d(i1,i2)
    lonout=lon[i]
    latout=lat[i]
    outdata=np.c_[lonout,latout]
    np.savetxt(out+'v'+str(v)+'.'+rjust(str(k),4,'0')+'.velfile',outdata,fmt='%.4f\t%.4f')
    
    
    
    
