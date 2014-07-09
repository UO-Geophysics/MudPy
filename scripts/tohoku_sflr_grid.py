import numpy as np
from string import rjust

df=0.1
dc=0.25

#Make fine boxes
lon=np.arange(141,143+df,df)
lat=np.arange(35.5,35.8+df,df)
k=0
lonout=np.zeros(len(lon)*len(lat))
latout=lonout.copy()
for klon in range(len(lon)):
    for klat in range(len(lat)):
        lonout[k]=lon[klon]
        latout[k]=lat[klat]
        k+=1
lon1=lonout ; lat1=latout

lon=np.arange(141.5,143.4+df,df)
lat=np.arange(35.8+df,36.3+df,df)
k=0
lonout=np.zeros(len(lon)*len(lat))
latout=lonout.copy()
for klon in range(len(lon)):
    for klat in range(len(lat)):
        lonout[k]=lon[klon]
        latout[k]=lat[klat]
        k+=1    
lon2=lonout ; lat2=latout


lon=np.arange(142,143.6+df,df)
lat=np.arange(36.3+df,36.8+df,df)
k=0
lonout=np.zeros(len(lon)*len(lat))
latout=lonout.copy()
for klon in range(len(lon)):
    for klat in range(len(lat)):
        lonout[k]=lon[klon]
        latout[k]=lat[klat]
        k+=1
lon3=lonout ; lat3=latout

lon=np.arange(142.4,144+df,df)
lat=np.arange(36.8+df,37.2+df,df)
k=0
lonout=np.zeros(len(lon)*len(lat))
latout=lonout.copy()
for klon in range(len(lon)):
    for klat in range(len(lat)):
        lonout[k]=lon[klon]
        latout[k]=lat[klat]
        k+=1
lon4=lonout ; lat4=latout

lon=np.arange(142.8,144.3+df,df)
lat=np.arange(36.8+df,38.1+df,df)
k=0
lonout=np.zeros(len(lon)*len(lat))
latout=lonout.copy()
for klon in range(len(lon)):
    for klat in range(len(lat)):
        lonout[k]=lon[klon]
        latout[k]=lat[klat]
        k+=1
lon5=lonout ; lat5=latout

lon=np.arange(143,144.5+df,df)
lat=np.arange(38.1+df,39.5+df,df)
k=0
lonout=np.zeros(len(lon)*len(lat))
latout=lonout.copy()
for klon in range(len(lon)):
    for klat in range(len(lat)):
        lonout[k]=lon[klon]
        latout[k]=lat[klat]
        k+=1
lon6=lonout ; lat6=latout

lon=np.arange(143.5,144.5+df,df)
lat=np.arange(39.5+df,40.4+df,df)
k=0
lonout=np.zeros(len(lon)*len(lat))
latout=lonout.copy()
for klon in range(len(lon)):
    for klat in range(len(lat)):
        lonout[k]=lon[klon]
        latout[k]=lat[klat]
        k+=1
lon7=lonout ; lat7=latout

#This is the coarse grid
lon=np.arange(140,145+dc,dc)
lat=np.arange(35,41+dc,dc)
k=0
lonout=np.zeros(len(lon)*len(lat))
latout=lonout.copy()
for klon in range(len(lon)):
    for klat in range(len(lat)):
        lonout[k]=lon[klon]
        latout[k]=lat[klat]
        k+=1
lon8=lonout ; lat8=latout

lon=np.arange(140.5,142+df,df)
lat=np.arange(37,38.5+df,df)
k=0
lonout=np.zeros(len(lon)*len(lat))
latout=lonout.copy()
for klon in range(len(lon)):
    for klat in range(len(lat)):
        lonout[k]=lon[klon]
        latout[k]=lat[klat]
        k+=1
lon9=lonout ; lat9=latout

lon=np.arange(141.0,143+df,df)
lat=np.arange(38.5,40.5+df,df)
k=0
lonout=np.zeros(len(lon)*len(lat))
latout=lonout.copy()
for klon in range(len(lon)):
    for klat in range(len(lat)):
        lonout[k]=lon[klon]
        latout[k]=lat[klat]
        k+=1
lon10=lonout ; lat10=latout

lon=np.arange(142.1,142.7+df,df)
lat=np.arange(37.4,38.4+df,df)
k=0
lonout=np.zeros(len(lon)*len(lat))
latout=lonout.cospy()
for klon in range(len(lon)):
    for klat in range(len(lat)):
        lonout[k]=lon[klon]
        latout[k]=lat[klat]
        k+=1
lon11=lonout ; lat11=latout

lon=np.r_[lon1,lon2,lon3,lon4,lon5,lon6,lon7,lon8,lon9,lon10,lon11]
lat=np.r_[lat1,lat2,lat3,lat4,lat5,lat6,lat7,lat8,lat9,lat10,lat11]


f=open('/Volumes/Kanagawa/Slip_Inv/tohoku_tsunami_norefine/data/station_info/tsunami.sta','w')
for k in range(len(lon)):
    f.write('%s\t%.4f\t%.4f\n' %(rjust(str(k),4,'0'),lon[k],lat[k]))
f.close()

