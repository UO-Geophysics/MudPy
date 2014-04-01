'''
Diego Melgar, 02/2014

Parameter file
Project: Tohoku tests
Comment: 
'''

import runslip,forward

########                            GLOBALS                             ########

home='/Users/dmelgarm/Research/Slip_Inv/'
project_name='tohoku'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=0 #Initalize project
make_green=0 #Compute GFs
make_synthetics=0 #Compute synthetics for a given model at given stations
solve=1  # =1 solves forward problem or runs inverse calculation, =0 does nothing
###############################################################################

###############            Green function parameters               #############
coord_type=1 #(=0 for cartesian, =1 for lat/lon (will use Earth flattening transform)
hot_start=0  #Start at a certain subfault number
static=1  #=1 computes static GFs only, =0 computes the complete waveform
model_name='gil7.mod'   #Velocity model
rupture_name='tohoku.allgps.rupt'   #Rupture model, not needed for inversion
fault_name='tohoku.allgps.fault'    #Fault geometry
station_file='tohoku.sta'    #Station distribution
NFFT=1024 ; dt=0.8  #Time parameters
dk=0.2 ; pmin=0 ; pmax=1 ; kmax=10   #fk integration parameters
################################################################################

############                 Synthetics parameters               ###############
resample=1 #Resample synthetics to this rate (in Hz)
integrate=1 #=0 produces velocities, =1 makes displacements
################################################################################


#Initalize project folders
if init==1:
    runslip.init(home,project_name)

# Run green functions          
if make_green==1: 
    runslip.forward_setup(home,project_name,rupture_name) 
    runslip.make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,hot_start,coord_type,dk,pmin,pmax,kmax)  

#Now make synthetics for source/station pairs
if make_synthetics==1:
    hot_start=0
    runslip.make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,hot_start,coord_type)
    
#Run forward comptuation or solve for inverse problem
if solve==1:
    if static==0: #Forward problem (full waveforms)
        forward.waveforms(home,project_name,rupture_name,station_file,model_name,integrate,hot_start,resample)
    if static==1: #Forward problem (coseismics)
        forward.coseismics(home,project_name,rupture_name,station_file,model_name)
    
    
    