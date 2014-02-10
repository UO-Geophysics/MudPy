'''
Diego Melgar, 02/2014

Parameter file
Project: SSP0
Comment: Source inversion validation strike-slip point source excercise
'''

import runslip

########                            GLOBALS                             ########

home='/Users/dmelgarm/Research/Slip_Inv/'
project_name='ssp0'
################################################################################

#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=0 #Initalize project
make_green=1 #Compute GFs
make_synthetics=1 #Compute synthetics for a given model at given stations
###############################################################################

###############            Green function parameters               #############

model_name='BJ97.mod'   #Velocity model
NFFT=1024
dt=0.02

min_distance=0.8
max_distance=20
delta_distance=0.05

min_depth=10
max_depth=10
delta_depth=0 # =0 if you only want one depth
################################################################################

############                 Synthetics parameters               ###############

fault_name='ssp0.fault'    #Source model
station_file='ssp0.sta'    #Station distribution
integrate=1 #=0 produces velocities, =1 makes dispalcements
################################################################################


#Initalize project folders
if init==1:
    runslip.init(home,project_name)

# Run green functions          
if make_green==1:  
    runslip.make_green(home,project_name,model_name,min_depth,max_depth,delta_depth,min_distance,max_distance,delta_distance,dt,NFFT)  

#Now make synthetics for source/station pairs
if make_synthetics==1:
    runslip.make_synthetics(home,project_name,station_file,fault_name,model_name,delta_depth,delta_distance,integrate)
    
    
    