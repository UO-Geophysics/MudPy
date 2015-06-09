'''
Diego Melgar, 05/2015

Parameter file for inverse problem
Project: 2015 Nepal earthquake 
Comment: 
'''

from mudpy import runslip
import numpy as np
from obspy.core import UTCDateTime

########                            GLOBALS                             ########
home='/Users/dmelgar/Slip_inv/'
project_name='Nepal_tutorial'
run_name='GPS_only'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####close   

init=0 #Initalize project
make_green=1#Compute GFs
make_synthetics=0 #Compute synthetics for a given model at given stations
G_from_file=0 # =0 read GFs and create a new G, =1 load G from file
invert=0  # =1 runs inversion, =0 does nothing
###############################################################################

###############           Green function parameters               #############
hot_start=0  #Start at a certain subfault number
model_name='avouac.mod'   #Velocity model
fault_name='nepal_10.fault'    #Fault geometry
GF_list='nepal_simple.gflist'#What GFs are to be computed for each station
tgf_file=None#'seafloor.sta'
G_name='nepal_3' #Either name of GF matrix to load or name to save GF matrix with
# Displacement and velocity waveform parameters
NFFT=512; dt=0.2
#Tsunami deformation parameters
tsunNFFT=128 ; tsun_dt=2.0
#fk-parameters
dk=0.1 ; pmin=0 ; pmax=1 ; kmax=20
################################################################################

#############               Inversion Parameters               ################# 
time_epi=UTCDateTime('2015-04-25T06:11:26')
epicenter=np.array([84.708,28.147,15]) 
rupture_speed=2.72 #Fastest rupture allowed in km/s
num_windows=1
reg_spatial=np.logspace(-4,0,num=10) #Set to False if you don't want to use it
reg_temporal=None #np.logspace(-1.5,-1.3,num=2) #Set to False if don't want to use it
nstrike=20 ; ndip=15 ; nfaults=(nstrike,ndip) #set nstrike to total no. of faults and ndip to 1 if using Tikh
beta=45 #Rotational offset (in degrees) applied to rake (0 for normal)
Ltype=0 # 0 for Tikhonov and 2 for Laplacian
solver='nnls' # 'lstsq','nnls'
top='free' ; bottom='locked' ; left='locked' ; right='locked' #'locked' or 'free'
bounds=(top,bottom,left,right)
################################################################################e=
########      Run-time modifications to the time series             ############
weight=True
decimate=None  #Decimate by constant (=None for NO decimation)
bandpass=None#np.array([0.5]) #Corner frequencies in Hz =None if no filter is desired
################################################################################

#Initalize project folders
if init==1:
    runslip.init(home,project_name)
    
# Run green functions          
if make_green==1 or make_synthetics==1:
    runslip.inversionGFs(home,project_name,GF_list,tgf_file,fault_name,model_name,
        dt,tsun_dt,NFFT,tsunNFFT,make_green,make_synthetics,dk,pmin,
        pmax,kmax,beta,time_epi,hot_start)   

#Run inversion
if invert==1:
    runslip.run_inversion(home,project_name,run_name,fault_name,model_name,GF_list,G_from_file,
            G_name,epicenter,rupture_speed,num_windows,reg_spatial,reg_temporal,
            nfaults,beta,decimate,bandpass,solver,bounds,weight,Ltype)

