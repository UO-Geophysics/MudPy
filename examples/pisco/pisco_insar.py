'''
Diego Melgar, 08/2015

Parameter file for inverse problem
Project: 2007 Pisco earthquake
Comment: Inversion using only InSAR data
'''

from mudpy import runslip
import numpy as np
from obspy.core import UTCDateTime
########                            GLOBALS                             ########
home='/Users/dmelgar/Slip_inv/'
project_name='Pisco_insar'
run_name='insar'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####close   

init=1 #Initalize project
make_green=0 #Compute GFs
make_synthetics=0 #Compute synthetics for a given model at given stations
G_from_file=0# =0 read GFs and create a new G, =1 load G from file
invert=0  # =1 runs inversion, =0 does nothing
###############################################################################

###############          view  Green function parameters               #############
ncpus=1
hot_start=0  #Start at a certain subfault number
model_name='pisco.mod'   #Velocity model
fault_name='pisco.fault'    #Fault geometry
GF_list='psico.gflist'#What GFs are to be computed for each station
tgf_file=None
G_name='insar_only' #Either name of GF matrix to load or name to save GF matrix with
# Displacement and velocity waveform parameters
NFFT=1024; dt=0.2
#Tsunami deformation parameters
tsunNFFT=128 ; tsun_dt=2.0
#fk-parameters
dk=0.1 ; pmin=0 ; pmax=1 ; kmax=20
custom_stf=None
################################################################################

#############               Inversion Parameters               ################# 
time_hypo=UTCDateTime('2007-08-15T23:40:57')
hypocenter=np.array([-76.509,-13.354,39]) 
rupture_speed=3.3 #Fastest rupture allowed in km/s
num_windows=1
reg_spatial=np.logspace(-8,2,num=30) #Set to None if you don't want to use it
reg_temporal=None #Set to None if don't want to use it
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
bandpass=np.array([0.5]) #Corner frequencies in Hz =None if no filter is desired
################################################################################

#Initalize project folders
if init==1:
    runslip.init(home,project_name)
    
# Run green functions          
if make_green==1 or make_synthetics==1:
    runslip.inversionGFs(home,project_name,GF_list,tgf_file,fault_name,model_name,
        dt,tsun_dt,NFFT,tsunNFFT,make_green,make_synthetics,dk,pmin,
        pmax,kmax,beta,time_hypo,hot_start,ncpus,custom_stf)   

#Run inversion
if invert==1:
    runslip.run_inversion(home,project_name,run_name,fault_name,model_name,GF_list,G_from_file,
            G_name,hypocenter,rupture_speed,num_windows,reg_spatial,reg_temporal,
            nfaults,beta,decimate,bandpass,solver,bounds,weight,Ltype)

