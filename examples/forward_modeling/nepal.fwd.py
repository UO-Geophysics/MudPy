'''
Diego Melgar, 09/2015

Parameter file
Project: Coquimbo
Comment: 
'''

from mudpy import runslip,forward
from obspy.core import UTCDateTime
from numpy import array

########                            GLOBALS                             ########

home='/Users/dmelgar/Slip_inv/'
project_name='Nepal_forward_test'
run_name='fwd'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=0 #Initalize project
make_green=1 #Compute GFs
make_synthetics=1 #Compute synthetics for a given model at given stations
solve=0 # =1 solves forward problem or runs inverse calculation, =0 does nothing
###############################################################################

###############            Green function parameters               #############
ncpus=4
hot_start=0  #Start at a certain subfault number
static=0  #=1 computes static GFs only, =0 computes the complete waveform
tsunami=None
model_name='avouac.mod'   #Velocity model
rupture_name='nepal.rupt'   #Rupture model, not needed for inversion
fault_name='nepal.fault'    #Fault geometry
station_file='nepal.sta'   #Station distribution
GF_list='gps.gflist'#What GFs are to be computed for each station
G_from_file=True
G_name='nepal'
NFFT=512 ; dt=0.2  #Time parameters
dk=0.2 ; pmin=0 ; pmax=1 ; kmax=10   #fk integration parameters
custom_stf=None
################################################################################

############                 Synthetics parameters               ###############
time_epi=UTCDateTime('2015-04-25T06:11:26')
epicenter=array([84.708,28.147,15]) 
resample=None #Resample synthetics to this rate (in Hz)
integrate=1 #=0 produces velocities, =1 makes displacements
beta=0 #Rake offset, usually a good idea to keep at zero
num_windows=1
rupture_speed=3.0  #Only necessary if onset times are not identified in rupt file
stf_type='dreger'
################################################################################


#Initalize project folders
if init==1:
    runslip.init(home,project_name)

# Run green functions          
if make_green==1: 
    if ncpus<2:
        runslip.make_green(home,project_name,station_file,fault_name,model_name,
            dt,NFFT,static,tsunami,hot_start,dk,pmin,pmax,kmax)
    else:
        runslip.make_parallel_green(home,project_name,station_file,fault_name,
            model_name,dt,NFFT,static,tsunami,hot_start,dk,pmin,pmax,kmax,ncpus)


#Now make synthetics for source/station pairs
if make_synthetics==1:
    if ncpus<2: 
        runslip.make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,
                static,tsunami,beta,hot_start,time_epi)
    else:
        runslip.make_parallel_synthetics(home,project_name,station_file,fault_name,model_name,integrate,
                static,tsunami,beta,hot_start,time_epi,ncpus,custom_stf=None)
    
#Run forward comptuation or solve for inverse problem
if solve==1:
    if static==0: #Forward problem (full waveforms)
        forward.waveforms_fakequakes(home,project_name,fault_name,None,GF_list,
                model_name,run_name,dt,NFFT,G_from_file,G_name,source_time_function=stf_type,
                stf_falloff_rate=4.0,rupture_name=rupture_name,epicenter=epicenter,time_epi=time_epi)
    if static==1: #Forward problem (coseismics)
        forward.coseismics_matrix(home,project_name,rupture_name,station_file,G_from_file,G_name)
    
    
    