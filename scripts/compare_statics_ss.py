import numpy as np

stafile='/Users/dmelgarm/Research/Slip_Inv/ssp0/data/station_info/ssp0.sta'
syndir='/Users/dmelgarm/Research/Slip_Inv/ssp0/GFs/BJ97.mod_10.0/'
edfile='/Users/dmelgarm/Research/Slip_Inv/ed/ssp0.disp'
    
    
#Read stations
stations=np.genfromtxt(stafile,dtype="S6",usecols=0)
#Read edcmp
ed=np.loadtxt(edfile)
#Loop, read data and plot
for j in range(len(stations)):#Read coseismic results from EDCMP
    sta=stations[j]
    ecoseis=ed[j,2]
    ncoseis=ed[j,3]
    ucoseis=-ed[j,4]
    #Read coseismics from fk
    temp=np.loadtxt(syndir+sta+'.static.enu')
    efkcoseis=temp[0]
    nfkcoseis=temp[1]
    ufkcoseis=temp[2]
    #Compute difference
    de=100*(ecoseis-efkcoseis)/ecoseis
    dn=100*(ncoseis-nfkcoseis)/ncoseis
    du=100*(ucoseis-ufkcoseis)/ucoseis
    print 'Difference at station '+stations[j]+' is:'
    print '  east, '+str(de)+'%'
    print '  north, '+str(dn)+'%'
    print '  up, '+str(du)+'%'