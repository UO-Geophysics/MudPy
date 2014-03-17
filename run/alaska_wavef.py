'''
Diego Melgar, 02/2014

Parameter file
Project: Alaska Shield
Comment: Synthetics for the 1964 simulation
'''

import runslip,forward

########                            GLOBALS                             ########

home='/Users/dmelgarm/Research/Slip_Inv/'
project_name='alaska'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=0 #Initalize project
make_green=1 #Compute GFs
make_synthetics=1 #Compute synthetics for a given model at given stations
direction=1  #=1 for forward modeling, =0 for inversion
solve=1  # =1 solves forward problem or runs inverse calculation, =0 does nothing
hot_start=0
###############################################################################

###############            Green function parameters               #############
F=5
static=0  #=0 computes static GFs only, =1 computes the completer waveform
model_name='PS9.mod'   #Velocity model
rupture_name='alaska_small.rupt'   #Rupture model, not needed for inversion
fault_name='alaska_small.fault'    #Fault geometry
station_file='alaska.sta'    #Station distribution
NFFT=2048
dt=0.2
################################################################################

############                 Synthetics parameters               ###############

integrate=1 #=0 produces velocities, =1 makes dispalcements
################################################################################


#Initalize project folders
if init==1:
    runslip.init(home,project_name)

#Forward modelling?
if direction==1:
    runslip.forward_setup(home,project_name,rupture_name)

# Run green functions          
if make_green==1:  
    runslip.make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,hot_start)  

#Now make synthetics for source/station pairs
if make_synthetics==1:
    runslip.make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,hot_start)
    
#Run forward comptuation or solve for inverse problem
if solve==1:
    if direction==1 and static==0: #Forward problem (full waveforms)
        forward.waveforms(home,project_name,rupture_name,station_file,model_name,integrate,F)
    if direction==1 and static==1: #Forward problem (coseismics)
        forward.coseismics(home,project_name,rupture_name,station_file,model_name)
    
    
    