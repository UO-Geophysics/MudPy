from scipy import arange,zeros,r_
from string import rjust


lat=arange(35,41.1,0.25)
lon=arange(140,142.1,0.25)
lonout1=zeros(len(lat)*len(lon))
latout1=lonout1.copy()
k=0
for klat in range(len(lat)):
    for klon in range(len(lon)):
        lonout1[k]=lon[klon]
        latout1[k]=lat[klat]
        k+=1

lat=arange(35,41.1,0.2)
lon=arange(142.2,145.1,0.2)
lonout2=zeros(len(lat)*len(lon))
latout2=lonout2.copy()
k=0
for klat in range(len(lat)):
    for klon in range(len(lon)):
        lonout2[k]=lon[klon]
        latout2[k]=lat[klat]
        k+=1
        
lonout=r_[lonout1,lonout2]
latout=r_[latout1,latout2]

f=open('/Volumes/Kanagawa/Slip_Inv/tohoku_10s/data/station_info/tohoku.tgf','w')
for k in range(len(lonout)):
     f.write('%s\t%.4f\t%.4f\n' %('t'+rjust(str(int(k)),3,'0'),lonout[k],latout[k]))
f.close()