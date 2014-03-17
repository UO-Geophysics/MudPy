'''
Diego Melgar, 02/2014

Parameter file
Project: DSEF
Comment: Source inversion validation dip-slip finite-fault
'''


import runslip,forward

########                            GLOBALS                             ########

home='/Users/dmelgarm/Research/Slip_Inv/'
project_name='dsef'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=1 #Initalize project
make_green=0 #Compute GFs
make_synthetics=0 #Compute synthetics for a given model at given stations
direction=0  #=1 for forward modeling, =0 for inversion
solve=0 # =1 solves forward problem or runs inverse calculation, =0 does nothing
###############################################################################

###############            Green function parameters               #############
freq=5  #Lowpass corner
static=0  #=0 computes static GFs only, =1 computes the completer waveform
model_name='BJ97.mod'   #Velocity model
rupture_name='dsef.rupt'   #Rupture model, not needed for inversion
fault_name='sdef.fault'    #Fault geometry
station_file='dsef.sta'    #Station distribution
NFFT=1024
dt=0.02
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
    runslip.make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static)  

#Now make synthetics for source/station pairs
if make_synthetics==1:
    runslip.make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static)
    
#Run forward comptuation or solve for inverse problem
if solve==1:
    if direction==1 and static==0: #Forward problem (full waveforms)
        forward.waveforms(home,project_name,rupture_name,station_file,model_name,integrate,freq)
    if direction==1 and static==1: #Forward problem (coseismics)
        forward.coseismics(home,project_name,rupture_name,station_file,model_name)