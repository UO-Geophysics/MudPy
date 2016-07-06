'''
Parameter file for fakequakes run, stdev is set at 0.9mu as opposed to
0.6mu as in all previous runs. Also rejection percentage is 13% as opposed to 10%.
I'm using larger rise timescaling from Sommerville '99

USING 50 modes!!!
'''


from mudpy import fakequakes,runslip,forward,viewFQ
import numpy as np
from obspy.core import UTCDateTime



########                            GLOBALS                             ########
home='/Users/dmelgar/fakequakes/'
project_name='Cascadia_final1'
run_name='cascadia'
################################################################################


##############             What do you want to do??           ##################
init=0
make_ruptures=0
make_GFs=0
make_synthetics=0
make_waveforms=1
# Things that only need to be done once
load_distances=1
G_from_file=1
###############################################################################


#############                 Run-time parameters            ##################
ncpus=8

model_name='cascadia.mod'   # Velocity model
fault_name='cascadia30.fault'    # Fault geometry
slab_name='cascadia_slab.xyz'    # Slab 1.0 Ascii file
mesh_name='cascadia30.mshout'    # GMSH output file
distances_name='cascadia30_dist' # Name of dist matrix
rupture_list='onerup.list'
UTM_zone='10T'

#Station information
GF_list='cascadia_small.gflist'
G_name='small'

Nrealizations=100 # Number of fake ruptures to generate per magnitude bin
target_Mw=np.linspace(8.0,9.2,13) # Of what approximate magnitudes
max_slip=60 #Maximum slip (m) allowed in the model

# Displacement and velocity waveform parameters
NFFT=512 ; dt=1.0
#fk-parameters
dk=0.1 ; pmin=0 ; pmax=1 ; kmax=20
custom_stf=None

# Correlation function parameters
hurst=0.75
Ldip='auto'
Lstrike='auto'
lognormal=True
slip_standard_deviation=0.9

# Rupture parameters
time_epi=UTCDateTime('2016-09-07T14:42:26')
source_time_function='cosine' # options are 'triangle' or 'cosine'
num_modes=50
rake=90.0
rise_time_depths=[10,15] #Transition depths for rise time scaling
buffer_factor=0.5
###############################################################################



#Initalize project folders
if init==1:
    fakequakes.init(home,project_name)
    
#Generate rupture models
if make_ruptures==1:
    fakequakes.generate_ruptures(home,project_name,run_name,fault_name,slab_name,
            mesh_name,load_distances,distances_name,UTM_zone,target_Mw,model_name,
            hurst,Ldip,Lstrike,num_modes,Nrealizations,rake,buffer_factor,
            rise_time_depths,time_epi,max_slip,source_time_function,lognormal,
            slip_standard_deviation)
    print '-- Now making plots --\n'
    viewFQ.plot_KLslip(home,project_name,run_name,mesh_name)
                
# Prepare waveforms and synthetics       
if make_GFs==1 or make_synthetics==1:
    runslip.inversionGFs(home,project_name,GF_list,None,fault_name,model_name,
        dt,None,NFFT,None,make_GFs,make_synthetics,dk,pmin,
        pmax,kmax,0,time_epi,0,ncpus,custom_stf,impulse=True) 
        
if make_waveforms==1:
    forward.waveforms_fakequakes(home,project_name,fault_name,rupture_list,GF_list,
                model_name,run_name,dt,NFFT,G_from_file,G_name,source_time_function)