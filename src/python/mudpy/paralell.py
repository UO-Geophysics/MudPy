'''
Module for routines that use paralell computing
'''


def run_paralell_green(home,project_name,station_file,model_name,dt,NFFT,static,dk,pmin,pmax,kmax,tsunami,rank,size):
    '''
    Compute GFs using Zhu & Rivera code for a given velocity model, source depth
    and station file. This function will make an external system call to fk.pl
    
    IN:
        source: 1-row numpy array containig informaiton aboutt he source, lat, lon, depth, etc...
        station_file: File name with the station coordinates
        dt: Desired sampling interval for waveforms
        NFFT: No. of samples requested in waveform (must be power of 2)
        static: =0 if computing full waveforms, =1 if computing only the static field
        coord_type: =0 if problem is in cartesian coordinates, =1 if problem is in lat/lon

    OUT:
        log: Sysytem standard output and standard error for log
    '''
    import subprocess
    from os import chdir
    from shutil import copy,rmtree
    from numpy import genfromtxt
    from string import rjust
    from shlex import split
    from shutil import copy
    from glob import glob
    from mudpy.green import src2sta

    #What parameters are we using?
    if rank==0:
        out='''Running all processes with:
        home = %s
        project_name = %s
        station_file = %s
        model_name = %s
        static = %s
        tsunami = %s
        dt = %.3f
        NFFT = %d
        dk = %.3f
        pmin = %.3f
        pmax = %.3f
        kmax = %.3f
        ''' %(home,project_name,station_file,model_name,str(static),str(tsunami),dt,NFFT,dk,pmin,pmax,kmax)
        print out
    #read your corresponding source file
    source=genfromtxt(home+project_name+'/data/model_info/mpi_source.'+str(rank)+'.fault')
    for ksource in range(len(source)):
        #Where should I be working boss?
        depth='%.4f' % source[ksource,3]
        subfault=rjust(str(int(source[ksource,0])),4,'0')
        if tsunami==False and static==0:
            subfault_folder=home+project_name+'/GFs/dynamic/'+model_name+'_'+depth+'.sub'+subfault
        elif tsunami==True and static==0:
            subfault_folder=home+project_name+'/GFs/tsunami/'+model_name+'_'+depth+'.sub'+subfault
        elif static==1:
            subfault_folder=home+project_name+'/GFs/static/'
        #Copy velocity model file
        copy(home+project_name+'/structure/'+model_name,subfault_folder+'/'+model_name)
        #Move to work folder
        chdir(subfault_folder)
        #Get station distances to source
        d,az=src2sta(station_file,source[ksource,:])
        #Make distance string for system call
        diststr=''
        for k in range(len(d)):
            diststr=diststr+' %.3f' % d[k] #Truncate distance to 3 decimal palces (meters)
        # Keep the user informed, lest he get nervous
        print 'MPI: processor #',rank,'is now working on subfault',int(source[ksource,0]),'(',ksource+1,'/',len(source),')'
        #Make the calculation
        if static==0: #Compute full waveform
            command=split("fk.pl -M"+model_name+"/"+depth+"/f -N"+str(NFFT)+"/"+str(dt)+'/1/'+repr(dk)+' -P'+repr(pmin)+'/'+repr(pmax)+'/'+repr(kmax)+diststr)
            p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            p.communicate() 
            # Move files up one level and delete folder created by fk
            files_list=glob(subfault_folder+'/'+model_name+'_'+depth+'/*.grn*')
            for f in files_list:
                newf=subfault_folder+'/'+f.split('/')[-1]
                copy(f,newf)
            rmtree(subfault_folder+'/'+model_name+'_'+depth)
        else: #Compute only statics
            write_file=subfault_folder+model_name+'.static.'+depth+'.sub'+subfault
            command=split("fk.pl -M"+model_name+"/"+depth+"/f -N1 "+diststr)
            p=subprocess.Popen(command,stdout=open(write_file,'w'),stderr=subprocess.PIPE)
            p.communicate()
            
            
            
#If main entry point
if __name__ == '__main__':
    import sys
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # Map command line arguments to function arguments.
    if sys.argv[1]=='run_paralell_green':
        #Parse command line arguments
        home=sys.argv[2]
        project_name=sys.argv[3]
        station_file=sys.argv[4]
        model_name=sys.argv[5]
        dt=float(sys.argv[6])
        NFFT=int(sys.argv[7])
        static=int(sys.argv[8])
        dk=float(sys.argv[9])
        pmin=int(sys.argv[10])
        pmax=int(sys.argv[11])
        kmax=int(sys.argv[12])
        tsunami=sys.argv[13]=='True'
        run_paralell_green(home,project_name,station_file,model_name,dt,NFFT,static,dk,pmin,pmax,kmax,tsunami,rank,size)
    else:
        print 'ERROR: You''re not allowed to run '+sys.argv[1]+' from the shell or it does not exist'