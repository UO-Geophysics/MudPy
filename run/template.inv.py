'''
Diego Melgar, 02/2014

Parameter file for inverse problem
Project: Tohoku
Comment: Inversion of Tohoku data with 10s rise times, epicenter location from Chu et al. (EPS,2011)
'''

from mudpy import runslip
import numpy as np
from obspy.core import UTCDateTime
########                            GLOBALS                             ########
home='/Volumes/Kanagawa/Slip_Inv/'
project_name='tohoku_10s'
run_name='20win_42_fine2'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=1 #Initalize project
make_green=0 #Compute GFs
make_synthetics=0 #Compute synthetics for a given model at given stations
G_from_file=0 # =0 read GFs and create a new G, =1 load G from file
invert=0  # =1 runs inversion, =0 does nothing
###############################################################################

###############            Green function parameters               #############
coord_type=1 # =0 for cartesian, =1 for lat/lon (will use Earth flattening transform)
hot_start=0  #Start at a certain subfault number
model_name='fnet_cmt.mod'   #Velocity model
fault_name='tohoku.fault'    #Fault geometry
station_file='tohoku.sta'    #Station distribution
GF_list='tohoku.dvt.gflist'#What GFs are to be computed for each station
tgf_file='tohoku.tgf'
G_name='fnet_20win_vr3.5_300s_dvt_mul1disp_mul40vel_mul1.1_OB2tsun.g.npy' #Either name of GF matrix to load or name to save GF matrix with
# Displacement and velocity waveform parameters
NFFT=2048 ; dt=0.25
#Tsunami deformation parameters
tsunNFFT=128 ; tsun_dt=4.0
#fk-parameters
dk=0.1 ; pmin=0 ; pmax=1 ; kmax=20
################################################################################

#############               Inversion Parameters               ################# 
epicenter=np.array([142.68,38.19,21.0])    #lon,lat,depth(positive in km)
time_epi=UTCDateTime('2011-03-11T05:46:23')
rupture_speed=3.5 #Fastest rupture allowed in km/s
num_windows=20
reg_spatial=np.logspace(0.25,0.35,num=2) #Set to False if you don't want to use it
reg_temporal=np.logspace(0,0,num=1) #Set to False if don't want to use it
nstrike=21 ; ndip=9 ; nfaults=(nstrike,ndip)
beta=45 #Rotational offset (in degrees) applied to rake (0 for normal)
solver='nnls' # 'lstsq','nnls'
top='free' ; bottom='locked' ; left='locked' ; right='locked' #'locked' or 'free'
bounds=(top,bottom,left,right) 
################################################################################

########      Run-time modifications to the time series             ############
decimate=4  #Decimate by constant (=None for NO decimation)
lowpass=None #Low pass corner frequency in Hz =None if no filter is desired
################################################################################

#Initalize project folders
if init==1:
    runslip.init(home,project_name)
    
# Run green functions          
if make_green==1 or make_synthetics==1:
    runslip.inversionGFs(home,project_name,GF_list,tgf_file,fault_name,model_name,
        dt,tsun_dt,NFFT,tsunNFFT,coord_type,make_green,make_synthetics,dk,pmin,
        pmax,kmax,beta,time_epi,hot_start)   

#Run inversion
if invert==1:
    runslip.run_inversion(home,project_name,run_name,fault_name,model_name,GF_list,G_from_file,
            G_name,epicenter,rupture_speed,num_windows,coord_type,reg_spatial,reg_temporal,
            nfaults,beta,decimate,lowpass,solver,bounds)
