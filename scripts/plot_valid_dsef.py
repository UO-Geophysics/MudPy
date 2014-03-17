'''
Plot computed synthetics versus known solutions from the inversion validation workshop results
'''

import numpy as np
import matplotlib.pyplot as pl
from obspy import read
from obspy.signal.filter import lowpass
from scipy.integrate import cumtrapz

stafile='/Users/dmelgarm/Research/Slip_Inv/dsef/data/station_info/dsef.sta'
syndir='/Users/dmelgarm/Research/Slip_Inv/dsef/GFs/BJ97.mod_10.0/'
fwddir='/Users/dmelgarm/Research/Slip_Inv/dsef/output/forward_models/'
edfile='/Users/dmelgarm/Research/Slip_Inv/ed/ssp0.disp'
val1='/Users/dmelgarm/Research/Slip_Inv/fk_validation/dsef0_AX/'
val2='/Users/dmelgarm/Research/Slip_Inv/fk_validation/dsef0_CS1/'
val3='/Users/dmelgarm/Research/Slip_Inv/fk_validation/dsef0_CS2/'
dt1=0.0225
dt2=0.005
dt3=0.01
tmin=0
tmax=20
#Filter stuff
freq=5
#Integrate?
integrate=1


pl.close("all")
#Read stations
stations=np.genfromtxt(stafile,dtype="S6",usecols=0)
x=np.genfromtxt(stafile,dtype="f8",usecols=1)
y=np.genfromtxt(stafile,dtype="f8",usecols=2)
efkcoseis=np.zeros(stations.shape)
nfkcoseis=np.zeros(stations.shape)
ufkcoseis=np.zeros(stations.shape)
#Loop, read data and plot
for j in range(len(stations)):
    k=j
    sta=stations[k]
    print j
    #Read validation 1
    e_v1=lowpass(np.genfromtxt(val1+sta.lower()+'.syn',dtype="f8",usecols=0,skip_header=4),freq,df=1/dt1,zerophase=True)
    n_v1=lowpass(np.genfromtxt(val1+sta.lower()+'.syn',dtype="f8",usecols=1,skip_header=4),freq,df=1/dt1,zerophase=True)
    u_v1=lowpass(np.genfromtxt(val1+sta.lower()+'.syn',dtype="f8",usecols=2,skip_header=4),freq,df=1/dt1,zerophase=True)
    t_v1=np.arange(0,len(e_v1)*dt1,dt1)#-0.8
    #Read validation 2
    e_v2=lowpass(np.genfromtxt(val2+sta.lower()+'.syn',dtype="f8",usecols=0,skip_header=4),freq,df=1/dt2,zerophase=True)
    n_v2=lowpass(np.genfromtxt(val2+sta.lower()+'.syn',dtype="f8",usecols=1,skip_header=4),freq,df=1/dt2,zerophase=True)
    u_v2=lowpass(np.genfromtxt(val2+sta.lower()+'.syn',dtype="f8",usecols=2,skip_header=4),freq,df=1/dt2,zerophase=True)
    t_v2=np.arange(0,len(e_v2)*dt2,dt2)
    #Read validation 3
    e_v3=lowpass(np.genfromtxt(val3+sta.lower()+'.syn',dtype="f8",usecols=0,skip_header=4),freq,df=1/dt3,zerophase=True)
    n_v3=lowpass(np.genfromtxt(val3+sta.lower()+'.syn',dtype="f8",usecols=1,skip_header=4),freq,df=1/dt3,zerophase=True)
    u_v3=lowpass(np.genfromtxt(val3+sta.lower()+'.syn',dtype="f8",usecols=2,skip_header=4),freq,df=1/dt3,zerophase=True)
    t_v3=np.arange(0,len(e_v3)*dt3,dt3)
    #Read my computations
    if integrate==1: #Go to displacement-land
        temp=read(fwddir+sta+'.disp.e')
        tbegin=temp[0].stats.starttime.second+0.52387
        e=lowpass(temp[0].data,freq,df=temp[0].stats.sampling_rate,zerophase=True)
        temp=read(fwddir+sta+'.disp.n')
        n=lowpass(temp[0].data,freq,df=temp[0].stats.sampling_rate,zerophase=True)
        temp=read(fwddir+sta+'.disp.z')
        u=lowpass(temp[0].data,freq,df=temp[0].stats.sampling_rate,zerophase=True)
        t=temp[0].times()+tbegin
        e_v1=cumtrapz(e_v1,t_v1)
        n_v1=cumtrapz(n_v1,t_v1)
        u_v1=cumtrapz(u_v1,t_v1)
        t_v1=t_v1[1:]
        e_v2=cumtrapz(e_v2,t_v2)
        n_v2=cumtrapz(n_v2,t_v2)
        u_v2=cumtrapz(u_v2,t_v2)
        t_v2=t_v2[1:]
        e_v3=cumtrapz(e_v3,t_v3)
        n_v3=cumtrapz(n_v3,t_v3)
        u_v3=cumtrapz(u_v3,t_v3)
        t_v3=t_v3[1:]
    else:
        temp=read(fwddir+sta+'.vel.e')
        tbegin=temp[0].stats.sac.b
        e=lowpass(temp[0].data,freq,df=temp[0].stats.sampling_rate,zerophase=True)
        temp=read(fwddir+sta+'.vel.n')
        n=lowpass(temp[0].data,freq,df=temp[0].stats.sampling_rate,zerophase=True)
        temp=read(fwddir+sta+'.vel.z')
        u=lowpass(temp[0].data,freq,df=temp[0].stats.sampling_rate,zerophase=True)
        t=temp[0].times()+tbegin
    #Read coseismic results from EDCMP
    ed=np.loadtxt(edfile)
    ecoseis=ed[:,2]
    ncoseis=ed[:,3]
    ucoseis=-ed[:,4]
    #Read coseismics from fk
    temp=np.loadtxt(fwddir+sta+'.static.enu')
    efkcoseis[j]=temp[0]
    nfkcoseis[j]=temp[1]
    ufkcoseis[j]=temp[2]
    #Plot
    pl.figure()
    #EAST
    pl.subplot(311)
    pl.plot(t,e,'k')
    pl.plot(t_v1,e_v1,t_v2,e_v2,t_v3,e_v3,t,np.ones(t.shape)*efkcoseis[j])
    pl.ylabel('East m/s')
    pl.grid()
    pl.legend(['My Syntetics','Axitra','CompSyn','CompSyn2','FK-Coseis'])
    pl.xlim((tmin,tmax))
    #NORTH
    pl.subplot(312)
    pl.plot(t,n,'k')
    pl.plot(t_v1,n_v1,t_v2,n_v2,t_v3,n_v3,t,np.ones(t.shape)*nfkcoseis[j])
    pl.ylabel('North m/s')
    pl.grid()
    pl.xlim((tmin,tmax))
    #UP
    pl.subplot(313)
    pl.plot(t,u,'k')
    pl.plot(t_v1,u_v1,t_v2,u_v2,t_v3,u_v3,t,np.ones(t.shape)*ufkcoseis[j])
    pl.ylabel('Up m/s')
    pl.xlabel('Time (s)')
    pl.grid()
    pl.xlim((tmin,tmax))
    #Done plotting
    pl.show()
    pl.subplot(311)
    pl.title('Station '+sta)
pl.figure()
pl.quiver(x,y,efkcoseis,nfkcoseis)
pl.grid()
pl.show()