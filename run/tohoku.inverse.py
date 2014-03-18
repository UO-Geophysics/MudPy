'''
Diego Melgar, 02/2014

Parameter file for inverse problem
Project: Tohoku
Comment: Inversion fo Tohoku data
'''

import runslip,inverse

########                            GLOBALS                             ########

home='/Users/dmelgarm/Research/Slip_Inv/'
project_name='tohoku'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=0 #Initalize project
make_green=0 #Compute GFs
make_synthetics=1 #Compute synthetics for a given model at given stations
invert=0  # =1 runs inversion, =0 does nothing
###############################################################################

###############            Green function parameters               #############
coord_type=1 # =0 for cartesian, =1 for lat/lon (will use Earth flattening transform)
hot_start=0   #Start at a certain subfault number
static=2  # =2 computes statics and dynamics, =1 computes static GFs only, =0 computes only dynamics
model_name='gil7.mod'   #Velocity model
fault_name='tohoku.fault'    #Fault geometry
station_file='tohoku.sta'    #Station distribution
GF_list='tohoku.gflist'    #What GFs are to be computed for each station
# Displacement and velocity waveform parameters
NFFT=2048
dt=0.2
################################################################################



#Initalize project folders
if init==1:
    runslip.init(home,project_name)

# Run green functions          
if make_green==1 or make_synthetics==1:
    runslip.inversionGFs(home,project_name,GF_list,fault_name,model_name,dt,NFFT,coord_type,make_green,make_synthetics)
    
#Run inverse computation
#if solve==1:
#    if static==0: #Forward problem (full waveforms)
#        forward.waveforms(home,project_name,rupture_name,station_file,model_name,integrate)
#    if static==1: #Forward problem (coseismics)
#        forward.coseismics(home,project_name,rupture_name,station_file,model_name)
    
    
    