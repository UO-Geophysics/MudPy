'''
Parameter file for fakequakes run, with Christine's Hayward discretization

For a new project, call fakequakes.init(home, project_name) once. Copy the
bundled hayward.fault to data/model_info/, bbp_norcal.mod to structure/,
and four_stations.gflist to data/station_info/ before running this driver.
'''


from os.path import expanduser
from mudpy import fakequakes,runslip,forward
import numpy as np
from obspy.core import UTCDateTime



########                            GLOBALS                             ########
home=expanduser('~/fakequakes/')
project_name='hayward'
run_name='hayward'
################################################################################


##############             What do you want to do??           ##################
init=0
make_ruptures=1
make_GFs=0
make_synthetics=0
make_waveforms=0
# Things that only need to be done once
load_distances=0
G_from_file=0
###############################################################################


#############                 Run-time parameters            ##################
ncpus=1  # Increase after the first run and set load_distances=1 to reuse distances.

model_name='bbp_norcal.mod'   # Velocity model
fault_name='hayward.fault'    # Fault geometry
slab_name=None    # Slab 1.0 Ascii file, set to None for simple geometry
mesh_name=None    # GMSH output file, set to None for simple geometry
distances_name='heyward_dist' # Name of dist matrix
rupture_list='ruptures.list'
UTM_zone='10'
scaling_law='S' # T for thrust, S for strike-slip, N for normal

#Station information
GF_list='four_stations.gflist'
G_name='bard'

Nrealizations=10 # Number of fake ruptures to generate per magnitude bin
target_Mw=np.linspace(6.0,7.4,8) # Of what approximate magnitudes
max_slip=8 #Maximum slip (m) allowed in the model

# Displacement and velocity waveform parameters
NFFT=256 ; dt=1.0
#fk-parameters
dk=0.1 ; pmin=0 ; pmax=1 ; kmax=20
custom_stf=None

# Correlation function parameters
hurst=0.75
Ldip='MB2002'  # Mai & Beroza correlation scaling used by the original Hayward example
Lstrike='MB2002'
lognormal=True
slip_standard_deviation=0.9

# Rupture parameters
time_epi=UTCDateTime('2016-09-07T14:42:26')
source_time_function='dreger' # options are 'triangle' or 'cosine'
stf_falloff_rate=6 #Only affects Dreger STF, choose 4-8 are reasonable values
num_modes=500
rake=180
rise_time_depths=[5,8] #Transition depths for rise time scaling
rise_time='S1999'  # Somerville (1999) rise-time scaling
shear_wave_fraction_shallow=0.56
shear_wave_fraction_deep=0.8
###############################################################################



#Initalize project folders
if init==1:
    fakequakes.init(home,project_name)
    
#Generate rupture models
if make_ruptures==1:
    fakequakes.generate_ruptures(home,project_name,run_name,fault_name,slab_name,
            mesh_name,load_distances,distances_name,UTM_zone,target_Mw,model_name,
            hurst,Ldip,Lstrike,num_modes,Nrealizations,rake,rise_time=rise_time,
            rise_time_depths=rise_time_depths,time_epi=time_epi,max_slip=max_slip,
            source_time_function=source_time_function,lognormal=lognormal,
            slip_standard_deviation=slip_standard_deviation,scaling_law=scaling_law,ncpus=ncpus,
            shear_wave_fraction_shallow=shear_wave_fraction_shallow,
            shear_wave_fraction_deep=shear_wave_fraction_deep)
                
# Prepare waveforms and synthetics       
if make_GFs==1 or make_synthetics==1:
    runslip.inversionGFs(home,project_name,GF_list,None,fault_name,model_name,
        dt,None,NFFT,None,make_GFs,make_synthetics,dk,pmin,
        pmax,kmax,0,time_epi,0,ncpus,custom_stf,impulse=True) 
        
if make_waveforms==1:
    forward.waveforms_fakequakes(home,project_name,fault_name,rupture_list,GF_list,
                model_name,run_name,dt,NFFT,G_from_file,G_name,source_time_function,
                stf_falloff_rate=stf_falloff_rate)
