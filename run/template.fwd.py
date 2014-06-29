'''
Diego Melgar, 06/2014

Parameter file
Project: Tohoku forward tsunami modeling
Comment: 
'''

from mudpy import runslip,forward
from obspy.core import UTCDateTime

########                            GLOBALS                             ########

home='/Volumes/Kanagawa/Slip_Inv/'
project_name='tohoku_tsunami'
run_name='dv'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=0 #Initalize project
make_green=0 #Compute GFs
make_synthetics=0 #Compute synthetics for a given model at given stations
solve=1  # =1 solves forward problem or runs inverse calculation, =0 does nothing
###############################################################################

###############            Green function parameters               #############
coord_type=1 #(=0 for cartesian, =1 for lat/lon (will use Earth flattening transform)
hot_start=75  #Start at a certain subfault number
static=0  #=1 computes static GFs only, =0 computes the complete waveform
tsunami=1
model_name='fnet_cmt.mod'   #Velocity model
rupture_name='tohoku_dv.rupt'   #Rupture model, not needed for inversion
fault_name='tohoku_tsunami.fault'    #Fault geometry
station_file='tsunami.sta'   #Station distribution
NFFT=128 ; dt=4  #Time parameters
dk=0.2 ; pmin=0 ; pmax=1 ; kmax=10   #fk integration parameters
################################################################################

############                 Synthetics parameters               ###############
time_epi=UTCDateTime('2011-03-11T05:46:23')
resample=2 #Resample synthetics to this rate (in Hz)
integrate=1 #=0 produces velocities, =1 makes displacements
beta=0 #Rake offset, usually a good idea to keep at zero
################################################################################


#Initalize project folders
if init==1:
    runslip.init(home,project_name)

# Run green functions          
if make_green==1: 
    runslip.make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,
                tsunami,beta,hot_start,coord_type,dk,pmin,pmax,kmax)  

#Now make synthetics for source/station pairs
if make_synthetics==1:
    #runslip.rupt2fault(home,project_name,rupture_name) 
    runslip.make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,
                static,tsunami,beta,hot_start,coord_type,time_epi)
    
#Run forward comptuation or solve for inverse problem
if solve==1:
    if static==0: #Forward problem (full waveforms)
        forward.waveforms(home,project_name,rupture_name,station_file,model_name,run_name,integrate,tsunami,hot_start,resample)
    if static==1: #Forward problem (coseismics)
        forward.coseismics(home,project_name,rupture_name,run_name,station_file)
    
    
    