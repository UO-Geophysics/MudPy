import numpy as np
import matplotlib.pyplot as pl

stafile='/Users/dmelgarm/Research/Slip_Inv/dsp0/data/station_info/dsp0.sta'
syndir='/Users/dmelgarm/Research/Slip_Inv/dsp0/output/forward_models/'
edfile='/Users/dmelgarm/Research/Slip_Inv/ed/dsp0.disp'
    
    
#Read stations
stations=np.genfromtxt(stafile,dtype="S6",usecols=0)
#Read edcmp
ed=np.loadtxt(edfile)
#Loop, read data and plot
ecoseis=np.empty(20)
ncoseis=np.empty(20)
ucoseis=np.empty(20)
efkcoseis=np.empty(20)
nfkcoseis=np.empty(20)
ufkcoseis=np.empty(20)
e=np.empty(20)
n=np.empty(20)
for j in range(len(stations)):#Read coseismic results from EDCMP
    sta=stations[j]
    print j
    n[j]=ed[j,0]
    e[j]=ed[j,1]
    ecoseis[j]=ed[j,3]
    ncoseis[j]=ed[j,2]
    ucoseis[j]=-ed[j,4]
    #Read coseismics from fk
    temp=np.loadtxt(syndir+sta+'.static.enu')
    efkcoseis[j]=temp[0]
    nfkcoseis[j]=temp[1]
    ufkcoseis[j]=temp[2]
    #Compute difference
    de=100*(ecoseis[j]-efkcoseis[j])/ecoseis[j]
    dn=100*(ncoseis[j]-nfkcoseis[j])/ncoseis[j]
    du=100*(ucoseis[j]-ufkcoseis[j])/ucoseis[j]
    print 'Difference at station '+stations[j]+' is:'
    print '  east, '+str(de)+'%'
    print '  north, '+str(dn)+'%'
    print '  up, '+str(du)+'%'
pl.close("all")
pl.quiver(e,n,efkcoseis,nfkcoseis)
pl.quiver(e,n,ecoseis,ncoseis,color='r')
pl.grid()
pl.show()
    