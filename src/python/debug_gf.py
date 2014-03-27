from runslip import make_green
from numpy import loadtxt
import matplotlib.pyplot as pl
from obspy import read
from glob import glob
from os.path import split
from obspy.core.utcdatetime import UTCDateTime
t0=UTCDateTime(1970, 1, 1, 0, 0, 0, 0)



home='/Users/dmelgarm/Research/Slip_Inv/'
project_name='_debug'
model_name='gil7.mod'   #Velocity model
fault_name='tohoku_0001.fault'    #Fault geometry
station_file='tohoku_one.sta'    #Station distribution
# Displacement and velocity waveform parameters
NFFT=512
dt=1
static=0
coord_type=1
hot_start=0
compute=0
plot=1
xl=[0,512]

dir1='/Users/dmelgarm/Research/Slip_Inv/_debug/structure/tohoku_dk0.1dt0.25N2048kmax20z60/'
dir2='/Users/dmelgarm/Research/Slip_Inv/_debug/structure/tohoku_dk0.1dt0.25N2048kmax20z60/'
dist1='100'
dist2='100'

if compute==1:
    #Make GFs
    #station_file=home+project_name+'/data/station_info/'+station_file
    make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,hot_start,coord_type)
if plot==1:
    pl.close("all")
    pl.subplot(331)
    fname=glob(dir1+dist1+'*.0')[0]
    st1=read(fname) ; st1[0].integrate()
    fname=glob(dir2+dist2+'*.0')[0]
    st2=read(fname) ; st2[0].integrate()
    pl.plot(st1[0].times(),st1[0].data,st2[0].times(),st2[0].data)
    pl.grid()
    pl.xlim(xl)
    pl.legend(['dk=0.1,dt=1','dk=0.001,dt=1'])
    pl.ylabel('Vel m/s')
    pl.title('grn.0')
    pl.subplot(332)
    fname=glob(dir1+dist1+'*.1')[0]
    st1=read(fname) ; st1[0].integrate()
    fname=glob(dir2+dist2+'*.1')[0]
    st2=read(fname) ; st2[0].integrate()
    pl.plot(st1[0].times(),st1[0].data,st2[0].times(),st2[0].data)
    pl.grid()
    pl.xlim(xl)
    pl.title('grn.1')
    pl.subplot(333)
    fname=glob(dir1+dist1+'*.2')[0]
    st1=read(fname) ; st1[0].integrate()
    fname=glob(dir2+dist2+'*.2')[0]
    st2=read(fname) ; st2[0].integrate()
    pl.plot(st1[0].times(),st1[0].data,st2[0].times(),st2[0].data)
    pl.grid()
    pl.xlim(xl)
    pl.title('grn.2')
    
    pl.subplot(334)
    fname=glob(dir1+dist1+'*.3')[0]
    st1=read(fname) ; st1[0].integrate()
    fname=glob(dir2+dist2+'*.3')[0]
    st2=read(fname) ; st2[0].integrate()
    pl.plot(st1[0].times(),st1[0].data,st2[0].times(),st2[0].data)
    pl.grid()
    pl.xlim(xl)
    pl.ylabel('Vel m/s')
    pl.title('grn.3')
    pl.subplot(335)
    fname=glob(dir1+dist1+'*.4')[0]
    st1=read(fname) ; st1[0].integrate()
    fname=glob(dir2+dist2+'*.4')[0]
    st2=read(fname) ; st2[0].integrate()
    pl.plot(st1[0].times(),st1[0].data,st2[0].times(),st2[0].data)
    pl.grid()
    pl.xlim(xl)
    pl.title('grn.4')
    pl.subplot(336)
    fname=glob(dir1+dist1+'*.5')[0]
    st1=read(fname) ; st1[0].integrate()
    fname=glob(dir2+dist2+'*.5')[0]
    st2=read(fname) ; st2[0].integrate()
    pl.plot(st1[0].times(),st1[0].data,st2[0].times(),st2[0].data)
    pl.grid()
    pl.xlim(xl)
    pl.title('grn.5')

    
    pl.subplot(337)
    fname=glob(dir1+dist1+'*.6')[0]
    st1=read(fname) ; st1[0].integrate()
    fname=glob(dir2+dist2+'*.6')[0]
    st2=read(fname) ; st2[0].integrate()
    pl.plot(st1[0].times(),st1[0].data,st2[0].times(),st2[0].data)
    pl.grid()
    pl.xlim(xl)
    pl.xlabel('Seconds')
    pl.ylabel('Vel m/s')
    pl.title('grn.6')
    pl.subplot(338)
    fname=glob(dir1+dist1+'*.7')[0]
    st1=read(fname) ; st1[0].integrate()
    fname=glob(dir2+dist2+'*.7')[0]
    st2=read(fname) ; st2[0].integrate()
    pl.plot(st1[0].times(),st1[0].data,st2[0].times(),st2[0].data)
    pl.grid()
    pl.xlim(xl)
    pl.xlabel('Seconds')
    pl.title('grn.7')
    pl.subplot(339)
    fname=glob(dir1+dist1+'*.8')[0]
    st1=read(fname) ; st1[0].integrate()
    fname=glob(dir2+dist2+'*.8')[0]
    st2=read(fname) ; st2[0].integrate()
    pl.plot(st1[0].times(),st1[0].data,st2[0].times(),st2[0].data)
    pl.grid()
    pl.xlim(xl)
    pl.show()
    pl.xlabel('Seconds')
    pl.title('grn.8')
    