'''
Diego Melgar, 03/2017

Parameter file for inverse problem
Project: 2015 Nepal earthquake

Version 2 with updated parameter file options
'''

from mudpy import runslip
import numpy as np
from obspy.core import UTCDateTime

########                            GLOBALS                             ########
home='/Users/dmelgar/Slip_inv/'
project_name='Nepal_example'
run_name='example'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####close   

init=0 #Initalize project
make_green=1 #Compute GFs
make_synthetics=1 #Compute synthetics for a given model at given stations
G_from_file=0# =0 read GFs and create a new G, =1 load G from file
invert=1  # =1 runs inversion, =0 does nothing
###############################################################################

###############          view  Green function parameters               #############
ncpus=8
hot_start=0 #Start at a certain subfault number
model_name='avouac.mod'   #Velocity model
fault_name='nepal_10.fault'    #Fault geometry
GF_list='nepal_simple_v2.gflist'
tgf_file=None
G_name='nepal' #Either name of GF matrix to load or name to save GF matrix with
# Displacement and velocity wcloseaveform parameters
NFFT=512; dt=0.2
#Tsunami deformation parameters
tsunNFFT=64 ; tsun_dt=2.0
#fk-parameters
dk=0.1 ; pmin=0 ; pmax=1 ; kmax=20
custom_stf=None
################################################################################

#############               Inversion Parameters               ################# 
time_epi=UTCDateTime('2015-04-25T06:11:26')
epicenter=np.array([84.708,28.147,15]) 
rupture_speed=3.2 #Fastest rupture allowed in km/s
num_windows=1
reg_spatial=np.logspace(-2,1,num=20)#np.logspace(-0.3,0,num=5) #Set to False if you don't want to use it np.array([0.0372759])
reg_temporal=None#np.logspace(-4,0,num=20)#Set to False if don't want to use it
nstrike=20 ; ndip=15 ; nfaults=(nstrike,ndip) #set nstrike to total no. of faults and ndip to 1 if using Tikh
beta=45 #Rotational offset (in degrees) applied to rake (0 for normal)
Ltype=2 # 0 for Tikhonov and 2 for Laplacian
solver='nnls' # 'lstsq','nnls'
top='locked' ; bottom='locked' ; left='locked' ; right='locked' #'locked' or 'free'
bounds=(top,bottom,left,right)
################################################################################e=

########      Run-time modifications to the time series             ############
weight=None
decimate=None  #Decimate by constant (=None for NO decimation)
# #Corner frequencies in Hz =None if no filter is desired
 # [0.5] is a low pass filter
 # [0.02,0.5] is a bdan pass filter
 # [0.02,np.inf] is a high pass filter
displacement_bandpass=np.array([0.5]) 
velocity_bandpass=None
tsunami_bandpass=None
bandpass=[displacement_bandpass,velocity_bandpass,tsunami_bandpass]
################################################################################

#Initalize project folders
if init==1:
    runslip.init(home,project_name)
    
# Run green functions          
if make_green==1 or make_synthetics==1:
    runslip.inversionGFs(home,project_name,GF_list,tgf_file,fault_name,model_name,
        dt,tsun_dt,NFFT,tsunNFFT,make_green,make_synthetics,dk,pmin,
        pmax,kmax,beta,time_epi,hot_start,ncpus,custom_stf)   

#Run inversion
if invert==1:
    runslip.run_inversion(home,project_name,run_name,fault_name,model_name,GF_list,G_from_file,
            G_name,epicenter,rupture_speed,num_windows,reg_spatial,reg_temporal,
            nfaults,beta,decimate,bandpass,solver,bounds,weight,Ltype)

