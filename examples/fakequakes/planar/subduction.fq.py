'''
Fakequakes runs with a simple planar geometry
'''


from mudpy import fakequakes,runslip,forward
import numpy as np
from obspy.core import UTCDateTime


########                            GLOBALS                             ########
home='/Users/dmelgarm/FakeQuakes/'
project_name='subduction_test'
run_name='subduction'
################################################################################


##############             What do you want to do??           ##################
init=1
make_ruptures=0
make_GFs=0
make_synthetics=0
make_waveforms=0
# Things that only need to be done once
load_distances=0
G_from_file=0
###############################################################################


#############                 Run-time parameters            ##################
ncpus=8

model_name='constants.mod'   # Velocity model
fault_name='synth_Mw9_2km_shifted.fault'    # Fault geometry
slab_name=None   # Slab 1.0 Ascii file
mesh_name=None    # GMSH output file
distances_name='test' # Name of dist matrix
rupture_list='ruptures.list'
UTM_zone='31N'
scaling_law='T' # T for thrust, S for strike-slip, N for normal

#Station information
GF_list='kermadec.gflist'
G_name='GFs'

Nrealizations=50 # Number of fake ruptures to generate per magnitude bin
target_Mw=np.arange(7.5,9.1,0.1) # Of what approximate magnitudes
#target_Mw=np.array([9.0])
max_slip=90 #Maximum slip (m) allowed in the model

# Displacement and velocity waveform parameters
NFFT=1024 ; dt=1.0
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
force_magnitude=True
time_epi=UTCDateTime('2016-09-07T14:42:26')
hypocenter=[0.8301,0.01,27.67]
force_hypocenter=True

source_time_function='dreger' # options are 'triangle' or 'cosine' or 'dreger'
num_modes=500
rake=90.0
rise_time_depths=[10,15] #Transition depths for rise time scaling
buffer_factor=0.5

#Force subfault to be center of model
no_random=True
subfault_hypo=3524
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
            slip_standard_deviation,scaling_law,force_magnitude=force_magnitude,
            hypocenter=hypocenter,force_hypocenter=force_hypocenter,no_random=no_random,
            shypo=subfault_hypo)
            
                
# Prepare waveforms and synthetics       
if make_GFs==1 or make_synthetics==1:
    runslip.inversionGFs(home,project_name,GF_list,None,fault_name,model_name,
        dt,None,NFFT,None,make_GFs,make_synthetics,dk,pmin,
        pmax,kmax,0,time_epi,0,ncpus,custom_stf,impulse=True) 
        
if make_waveforms==1:
    forward.waveforms_fakequakes(home,project_name,fault_name,rupture_list,GF_list,
                model_name,run_name,dt,NFFT,G_from_file,G_name,source_time_function)