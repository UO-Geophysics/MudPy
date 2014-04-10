'''
Diego Melgar, 02/2014

Parameter file for inverse problem
Project: Tohoku
Comment: Inversion of Tohoku data with 10s rise times in a halfspace. Using a rectangle instead of a triangle
'''

import runslip
import numpy as np
from obspy.core import UTCDateTime

########                            GLOBALS                             ########

home='/Users/dmelgarm/Research/Slip_Inv/'
project_name='tohoku_halfspace_rect'
run_name='halfspace'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=0 #Initalize project
make_green=0 #Compute GFs
make_synthetics=1 #Compute synthetics for a given model at given stations
G_from_file=0 # =0 read GFs and create a new G, =1 load G from file
invert=0  # =1 runs inversion, =0 does nothing
###############################################################################

###############            Green function parameters               #############
coord_type=1 # =0 for cartesian, =1 for lat/lon (will use Earth flattening transform)
hot_start=0   #Start at a certain subfault number
model_name='halfspace.mod'   #Velocity model
fault_name='tohoku.fault'    #Fault geometry
station_file='tohoku.sta'    #Station distribution
GF_list='tohoku.kaldisp_small.gflist' #What GFs are to be computed for each station
G_name='halfspace_1vr.g' #Either name of GF matrix to load or name to save GF matrix with
# Displacement and velocity waveform parameters
NFFT=2048 ; dt=0.25
#fk-parameters
dk=0.1 ; pmin=0 ; pmax=1 ; kmax=20
################################################################################

#############               Inversion Parameters               ################# 
epicenter=np.array([142.435,38.305,26.3305])    #lon,lat,depth(positive in km)
time_epi=UTCDateTime('2011-03-11T05:46:23')
rupture_speeds=np.array([2.0])#np.array([1.6,1.8,2.0,2.2,2.4,2.6,2.8])   #In km/s
regularization_parameter=np.logspace(-2,-2,num=1)
regularization_type='laplace'
top='free' ; bottom='locked' ; left='locked' ; right='locked' ; bounds=(top,bottom,left,right) #'locked' or 'free'
nstrike=21 ; ndip=9 ; nfaults=(nstrike,ndip)
################################################################################


#Initalize project folders
if init==1:
    runslip.init(home,project_name)

# Run green functions          
if make_green==1 or make_synthetics==1:
    runslip.inversionGFs(home,project_name,GF_list,fault_name,model_name,dt,NFFT,
                        coord_type,make_green,make_synthetics,dk,pmin,pmax,kmax,time_epi,hot_start)
    
#Run inversion
if invert==1:
    runslip.run_inversion(home,project_name,run_name,fault_name,model_name,GF_list,G_from_file,
                           G_name,epicenter,rupture_speeds,coord_type,bounds,regularization_type,
                           regularization_parameter,nfaults)
    


    
    
    