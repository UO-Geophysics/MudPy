'''
Diego Melgar, 04/2016

Parameter file
Project: Cascadia tsunami scenarios
Comment: First attempt
'''

from mudpy import runslip,forward
from obspy.core import UTCDateTime
from numpy import array

########                            GLOBALS                             ########

home='/Users/dmelgar/FakeQuakes/Forward_runs/'
project_name='Cascadia'
run_name='static'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=0 #Initalize project
make_green=1 #Compute GFs
make_synthetics=1 #Compute synthetics for a given model at given stations
solve=0  # =1 solves forward problem or runs inverse calculation, =0 does nothing
###############################################################################

###############            Green function parameters               #############
hot_start=0  #Start at a certain subfault number
static=1  #=1 computes static GFs only, =0 computes the complete waveform
tsunami=False
model_name='cascadia30.mod'   #Velocity model
rupture_name='usgs.rupt'   #Rupture model, not needed for inversion
fault_name='cascadia30.fault'    #Fault geometry
station_file='seafloor.sta'   #Station distribution
GF_list='nepal.gflist'#What GFs are to be computed for each station
NFFT=256 ; dt=1.0  #Time parameters
dk=0.2 ; pmin=0 ; pmax=1 ; kmax=10   #fk integration parameters
################################################################################

############                 Synthetics parameters               ###############
time_epi=UTCDateTime('2015-04-25T06:11:26')
epicenter=array([84.708,28.147,15]) 
resample=1 #Resample synthetics to this rate (in Hz)
integrate=1 #=0 produces velocities, =1 makes displacements
beta=0 #Rake offset, usually a good idea to keep at zero
rupture_speed=3.0 #Fastest rupture allowed in km/s
num_windows=1
################################################################################


#Initalize project folders
if init==1:
    runslip.init(home,project_name)

# Run green functions          
if make_green==1: 
    runslip.make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,
                tsunami,hot_start,dk,pmin,pmax,kmax)  

#Now make synthetics for source/station pairs
if make_synthetics==1:
    #runslip.rupt2fault(home,project_name,rupture_name) 
    runslip.make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,
                static,tsunami,beta,hot_start,time_epi)
    
#Run forward comptuation or solve for inverse problem
if solve==1:
    if static==0: #Forward problem (full waveforms)
        forward.waveforms_matrix(home,project_name,fault_name,rupture_name,station_file,GF_list,model_name,
                run_name,epicenter,time_epi,integrate,tsunami,hot_start,resample,beta,rupture_speed,
                num_windows,dt,NFFT)
    if static==1: #Forward problem (coseismics)
        forward.coseismics(home,project_name,rupture_name,station_file)
    
    
    