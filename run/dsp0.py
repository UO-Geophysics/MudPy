'''
Diego Melgar, 02/2014

Parameter file
Project: DSP0
Comment: Source inversion validation dip-slip point source excercise
'''

import runslip

########                            GLOBALS                             ########

home='/Users/dmelgarm/Research/Slip_Inv/'
project_name='dsp0'
################################################################################

#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=0 #Initalize project
make_green=1 #Compute GFs
make_synthetics=0 #Compute synthetics for a given model at given stations
###############################################################################

###############            Green function parameters               #############

model_name='BJ97.mod'   #Velocity model
NFFT=1024
dt=0.02

min_depth=10
max_depth=10
delta_depth=0 # =0 if you only want one depth
################################################################################

############                 Synthetics parameters               ###############

fault_name='dsp0_multi.fault'    #Source model
station_file='dsp0.sta'    #Station distribution
integrate=1 #=0 produces velocities, =1 makes dispalcements
################################################################################


#Initalize project folders
if init==1:
    runslip.init(home,project_name)

# Run green functions          
if make_green==1:  
    runslip.make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT)  

#Now make synthetics for source/station pairs
if make_synthetics==1:
    runslip.make_synthetics(home,project_name,station_file,fault_name,model_name,delta_depth,delta_distance,integrate)
    
    
    