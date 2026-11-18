'''
Fakequakes runs with a simple planar geometry

For a new project, call fakequakes.init(home, project_name) once. Copy the
bundled subduction_lowres.fault to data/model_info/, cascadia.mod to structure/,
and coastal_gps.gflist to data/station_info/ before running this driver.
'''


from os.path import expanduser
from mudpy import fakequakes,runslip,forward
import numpy as np
from obspy.core import UTCDateTime


########                            GLOBALS                             ########
home=expanduser('~/fakequakes/')
project_name='subduction_test'
run_name='subduction'
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

model_name='cascadia.mod'   # Velocity model
fault_name='subduction_lowres.fault'    # Fault geometry
slab_name=None   # Slab 1.0 Ascii file (only used for 3D fault)
mesh_name=None    # GMSH output file (only used for 3D fault)
distances_name='planar_subduction' # Name of distance matrix
rupture_list='ruptures.list'  # Don't change this (unless you know waht you're doing!)
UTM_zone='31'
scaling_law='T' # T for thrust, S for strike-slip, N for normal


#slip parameters
Nrealizations=5 # Number of fake ruptures to generate per magnitude bin
target_Mw=np.arange(7.5,8.4,0.1) # Of what approximate magnitudes
max_slip=60 #Maximum slip (m) allowed in the model

# Correlation function parameters
hurst=0.4
Ldip='auto'      # Melgar & Hayes (2019) correlation length scaling
Lstrike='auto'
lognormal=True
slip_standard_deviation=0.9
num_modes=500
rake=90.0

# Rupture parameters
force_magnitude=True    #Make the magnitudes EXACTLY the value in target_Mw
force_area=False   #Forces using the entire fault area defined by the .fault file as opposed to the scaling laws
no_random=False   #If true uses median length/width if false draws from prob. distribution
time_epi=UTCDateTime('2016-09-07T14:42:26')  #Defines the hypocentral time
hypocenter=None    # For a prescribed hypocenter set [lon, lat, depth_km] and force_hypocenter=True.
force_hypocenter=False     # Forces hypocenter to occur at specified lcoationa s opposed to random
mean_slip = None     #Provide path to file name of .rupt to be used as mean slip pattern
use_hypo_fraction = False  #If true use hypocenter PDF positions from Melgar & Hayes 2019, if false then selects at random

# Kinematic parameters
source_time_function='dreger' # options are 'triangle' or 'cosine' or 'dreger'
rise_time_depths=[10,15] #Transition depths for rise time scaling
rise_time='MH2017'
shear_wave_fraction_shallow=0.49
shear_wave_fraction_deep=0.8

#Station information (only used when syntehsizing waveforms)
GF_list='coastal_gps.gflist'
G_name='GFs'

# Displacement and velocity waveform parameters
NFFT=128 ; dt=1.0

#fk-parameters
dk=0.1 ; pmin=0 ; pmax=1 ; kmax=20
custom_stf=None

###############################################################################



#Initalize project folders
if init==1:
    fakequakes.init(home,project_name)
    
#Generate rupture models
if make_ruptures==1:
    fakequakes.generate_ruptures(home,project_name,run_name,fault_name,slab_name,
        mesh_name,load_distances,distances_name,UTM_zone,target_Mw,model_name,hurst,Ldip,
        Lstrike,num_modes,Nrealizations,rake,rise_time=rise_time,rise_time_depths=rise_time_depths,
        time_epi=time_epi,max_slip=max_slip,source_time_function=source_time_function,lognormal=lognormal,
        slip_standard_deviation=slip_standard_deviation,scaling_law=scaling_law,ncpus=ncpus,
        force_magnitude=force_magnitude,force_area=force_area,mean_slip_name=mean_slip,
        hypocenter=hypocenter,force_hypocenter=force_hypocenter,no_random=no_random,
        use_hypo_fraction=use_hypo_fraction,shear_wave_fraction_shallow=shear_wave_fraction_shallow,
        shear_wave_fraction_deep=shear_wave_fraction_deep)
            
                
# Prepare waveforms and synthetics       
if make_GFs==1 or make_synthetics==1:
    runslip.inversionGFs(home,project_name,GF_list,None,fault_name,model_name,
        dt,None,NFFT,None,make_GFs,make_synthetics,dk,pmin,
        pmax,kmax,0,time_epi,0,ncpus,custom_stf,impulse=True) 
       
# Synthesize the waveforms
if make_waveforms==1:
    forward.waveforms_fakequakes(home,project_name,fault_name,rupture_list,GF_list,
                model_name,run_name,dt,NFFT,G_from_file,G_name,source_time_function)
