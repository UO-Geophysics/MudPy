'''
Module for routines that use paralell computing
'''


def run_parallel_green(home,project_name,station_file,model_name,dt,NFFT,static,dk,pmin,pmax,kmax,tsunami,rank,size):
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
    import os

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
            file_is_empty=True
            while file_is_empty:
                p=subprocess.Popen(command,stdout=open(write_file,'w'),stderr=subprocess.PIPE)
                p.communicate()
                if os.stat(write_file).st_size!=0: #File is NOT empty
                    file_is_empty=False
                else:
                    print 'Warning: I just had a mini-seizure and made an empty GF file on first try, re-running'
            #If file is empty run again
            

def run_parallel_synthetics(home,project_name,station_file,model_name,integrate,static,tsunami,time_epi,
                beta,custom_stf,impulse,rank,size,okada,mu_okada):
    '''
    Use green functions and compute synthetics at stations for a single source
    and multiple stations. This code makes an external system call to syn.c first it
    will make the external call for the strike-slip component then a second externall
    call will be made for the dip-slip component. The unit amount of moment is 1e15
    which corresponds to Mw=3.9333...
    
    IN:
        source: 1-row numpy array containig informaiton aboutt he source, lat, lon, depth, etc...
        station_file: File name with the station coordinates
        green_path: Directopry where GFs are stored
        model_file: File containing the Earth velocity structure
        integrate: =0 if youw ant velocity waveforms, =1 if you want displacements
        static: =0 if computing full waveforms, =1 if computing only the static field
        subfault: String indicating the subfault being worked on
        coord_type: =0 if problem is in cartesian coordinates, =1 if problem is in lat/lon
        
    OUT:
        log: Sysytem standard output and standard error for log
    '''

    import os
    import subprocess
    from mudpy.forward import get_mu
    from string import rjust
    from numpy import array,genfromtxt,loadtxt,savetxt,log10
    from obspy import read
    from shlex import split
    from mudpy.green import src2sta,rt2ne,origin_time,okada_synthetics
    from glob import glob
    from mudpy.green import silentremove
    from os import remove
    
    #What parameters are we using?
    if rank==0:
        out='''Running all processes with:
        home = %s
        project_name = %s
        station_file = %s
        model_name = %s
        integrate = %s
        static = %s
        tsunami = %s
        time_epi = %s
        beta = %d
        custom_stf = %s
        impulse = %s
        okada = %s
        mu = %.2e
        ''' %(home,project_name,station_file,model_name,str(integrate),str(static),str(tsunami),str(time_epi),beta,custom_stf,impulse,okada,mu_okada)
        print out
        
    #Read your corresponding source file
    mpi_source=genfromtxt(home+project_name+'/data/model_info/mpi_source.'+str(rank)+'.fault')
    
    #Constant parameters
    rakeDS=90+beta #90 is thrust, -90 is normal
    rakeSS=0+beta #0 is left lateral, 180 is right lateral
    tb=50 #Number of samples before first arrival (should be 50, NEVER CHANGE, if you do then adjust in fk.pl)
    
    #Figure out custom STF
    if custom_stf.lower()!='none':
        custom_stf=home+project_name+'/GFs/STFs/'+custom_stf
    else:
        custom_stf=None
    
    #Load structure
    model_file=home+project_name+'/structure/'+model_name
    structure=loadtxt(model_file,ndmin=2)
    
    for ksource in range(len(mpi_source)):
        
        source=mpi_source[ksource,:]
        
        #Parse the soruce information
        num=rjust(str(int(source[0])),4,'0')
        xs=source[1]
        ys=source[2]
        zs=source[3]
        strike=source[4]
        dip=source[5]
        rise=source[6]
        if impulse==True:
            duration=0
        else:
            duration=source[7]
        ss_length=source[8]
        ds_length=source[9]
        ss_length_in_km=ss_length/1000.
        ds_length_in_km=ds_length/1000.
        strdepth='%.4f' % zs
        subfault=rjust(str(int(source[0])),4,'0')
        if static==0 and tsunami==0:  #Where to save dynamic waveforms
            green_path=home+project_name+'/GFs/dynamic/'+model_name+"_"+strdepth+".sub"+subfault+"/"
        if static==0 and tsunami==1:  #Where to save dynamic waveforms
            green_path=home+project_name+'/GFs/tsunami/'+model_name+"_"+strdepth+".sub"+subfault+"/"
        if static==1:  #Where to save statics
            green_path=home+project_name+'/GFs/static/'
        staname=genfromtxt(station_file,dtype="S6",usecols=0)
        if staname.shape==(): #Single staiton file
            staname=array([staname])
        
        #Compute distances and azimuths
        d,az,lon_sta,lat_sta=src2sta(station_file,source,output_coordinates=True)
        
        #Get moment corresponding to 1 meter of slip on subfault
        mu=get_mu(structure,zs)
        Mo=mu*ss_length*ds_length*1.0
        Mw=(2./3)*(log10(Mo)-9.1)
        
        #Move to output folder
        os.chdir(green_path)
        print 'Processor '+str(rank)+' is working on subfault '+str(int(source[0]))+' and '+str(len(d))+' stations '
        for k in range(len(d)):
            if static==0: #Compute full waveforms
                diststr='%.3f' % d[k] #Need current distance in string form for external call
                #Form the strings to be used for the system calls according to user desired options
                if integrate==1: #Make displ.
                    #First Stike-Slip GFs
                    if custom_stf==None:
                        commandSS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -D"+str(duration)+ \
                            "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".SS.disp.x -G"+green_path+diststr+".grn.0"
                        commandSS=split(commandSS) #Split string into lexical components for system call
                        #Now dip slip
                        commandDS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -D"+str(duration)+ \
                            "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".DS.disp.x -G"+green_path+diststr+".grn.0"
                        commandDS=split(commandDS)
                    else:
                        commandSS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -S"+custom_stf+ \
                            " -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".SS.disp.x -G"+green_path+diststr+".grn.0"
                        commandSS=split(commandSS) #Split string into lexical components for system call
                        #Now dip slip
                        commandDS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -S"+custom_stf+ \
                            " -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".DS.disp.x -G"+green_path+diststr+".grn.0"
                        commandDS=split(commandDS)
                else: #Make vel.
                    #First Stike-Slip GFs
                    if custom_stf==None:
                        commandSS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -D"+str(duration)+ \
                            "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".SS.vel.x -G"+green_path+diststr+".grn.0"
                        commandSS=split(commandSS)
                        #Now dip slip
                        commandDS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -D"+str(duration)+ \
                            "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".DS.vel.x -G"+green_path+diststr+".grn.0"
                        commandDS=split(commandDS)
                    else:
                        commandSS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -S"+custom_stf+ \
                            " -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".SS.vel.x -G"+green_path+diststr+".grn.0"
                        commandSS=split(commandSS)
                        #Now dip slip
                        commandDS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -S"+custom_stf+ \
                            " -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".DS.vel.x -G"+green_path+diststr+".grn.0"
                        commandDS=split(commandDS)
                #Run the strike- and dip-slip commands (make system calls)
                p=subprocess.Popen(commandSS)
                p.communicate() 
                p=subprocess.Popen(commandDS)
                p.communicate()
                #Result is in RTZ system (+Z is down) rotate to NEZ with +Z up and scale to m or m/s
                if integrate==1: #'tis displacememnt
                    #Strike slip
                    if duration>0: #Is there a source time fucntion? Yes!
                        r=read(staname[k]+".subfault"+num+'.SS.disp.r')
                        t=read(staname[k]+".subfault"+num+'.SS.disp.t')
                        z=read(staname[k]+".subfault"+num+'.SS.disp.z')
                    else: #No! This is the impulse response!
                        r=read(staname[k]+".subfault"+num+'.SS.disp.ri')
                        t=read(staname[k]+".subfault"+num+'.SS.disp.ti')
                        z=read(staname[k]+".subfault"+num+'.SS.disp.zi')
                    ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                    #Scale to m and overwrite with rotated waveforms
                    n=r.copy()
                    n[0].data=ntemp/100
                    e=t.copy()
                    e[0].data=etemp/100
                    z[0].data=z[0].data/100
                    n=origin_time(n,time_epi,tb)
                    e=origin_time(e,time_epi,tb)
                    z=origin_time(z,time_epi,tb)
                    n.write(staname[k]+".subfault"+num+'.SS.disp.n',format='SAC')
                    e.write(staname[k]+".subfault"+num+'.SS.disp.e',format='SAC')
                    z.write(staname[k]+".subfault"+num+'.SS.disp.z',format='SAC')
                    silentremove(staname[k]+".subfault"+num+'.SS.disp.r')
                    silentremove(staname[k]+".subfault"+num+'.SS.disp.t')
                    if impulse==True:
                        silentremove(staname[k]+".subfault"+num+'.SS.disp.ri')
                        silentremove(staname[k]+".subfault"+num+'.SS.disp.ti')
                        silentremove(staname[k]+".subfault"+num+'.SS.disp.zi')
                    #Dip Slip
                    if duration>0: #Is there a source time fucntion? Yes!
                        r=read(staname[k]+".subfault"+num+'.DS.disp.r')
                        t=read(staname[k]+".subfault"+num+'.DS.disp.t')
                        z=read(staname[k]+".subfault"+num+'.DS.disp.z')
                    else: #No! This is the impulse response!
                        r=read(staname[k]+".subfault"+num+'.DS.disp.ri')
                        t=read(staname[k]+".subfault"+num+'.DS.disp.ti')
                        z=read(staname[k]+".subfault"+num+'.DS.disp.zi')
                    ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                    n=r.copy()
                    n[0].data=ntemp/100
                    e=t.copy()
                    e[0].data=etemp/100
                    z[0].data=z[0].data/100
                    n=origin_time(n,time_epi,tb)
                    e=origin_time(e,time_epi,tb)
                    z=origin_time(z,time_epi,tb)
                    n.write(staname[k]+".subfault"+num+'.DS.disp.n',format='SAC')
                    e.write(staname[k]+".subfault"+num+'.DS.disp.e',format='SAC')
                    z.write(staname[k]+".subfault"+num+'.DS.disp.z',format='SAC')
                    silentremove(staname[k]+".subfault"+num+'.DS.disp.r')
                    silentremove(staname[k]+".subfault"+num+'.DS.disp.t')
                    if impulse==True:
                        silentremove(staname[k]+".subfault"+num+'.DS.disp.ri')
                        silentremove(staname[k]+".subfault"+num+'.DS.disp.ti')
                        silentremove(staname[k]+".subfault"+num+'.DS.disp.zi')
                else: #Waveforms are velocity, as before, rotate from RT-Z to NE+Z and scale to m/s
                    #Strike slip
                    if duration>0: #Is there a source time fucntion? Yes!
                        r=read(staname[k]+".subfault"+num+'.SS.vel.r')
                        t=read(staname[k]+".subfault"+num+'.SS.vel.t')
                        z=read(staname[k]+".subfault"+num+'.SS.vel.z')
                    else: #No! This is the impulse response!
                        r=read(staname[k]+".subfault"+num+'.SS.vel.ri')
                        t=read(staname[k]+".subfault"+num+'.SS.vel.ti')
                        z=read(staname[k]+".subfault"+num+'.SS.vel.zi')
                    ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                    n=r.copy()
                    n[0].data=ntemp/100
                    e=t.copy()
                    e[0].data=etemp/100
                    z[0].data=z[0].data/100
                    n=origin_time(n,time_epi,tb)
                    e=origin_time(e,time_epi,tb)
                    z=origin_time(z,time_epi,tb)
                    n.write(staname[k]+".subfault"+num+'.SS.vel.n',format='SAC')
                    e.write(staname[k]+".subfault"+num+'.SS.vel.e',format='SAC')
                    z.write(staname[k]+".subfault"+num+'.SS.vel.z',format='SAC')
                    silentremove(staname[k]+".subfault"+num+'.SS.vel.r')
                    silentremove(staname[k]+".subfault"+num+'.SS.vel.t')
                    if impulse==True:
                        silentremove(staname[k]+".subfault"+num+'.SS.vel.ri')
                        silentremove(staname[k]+".subfault"+num+'.SS.vel.ti')
                        silentremove(staname[k]+".subfault"+num+'.SS.vel.zi')
                    #Dip Slip
                    if duration>0: #Is there a source time fucntion? Yes!
                        r=read(staname[k]+".subfault"+num+'.DS.vel.r')
                        t=read(staname[k]+".subfault"+num+'.DS.vel.t')
                        z=read(staname[k]+".subfault"+num+'.DS.vel.z')
                    else: #No! This is the impulse response!
                        r=read(staname[k]+".subfault"+num+'.DS.vel.ri')
                        t=read(staname[k]+".subfault"+num+'.DS.vel.ti')
                        z=read(staname[k]+".subfault"+num+'.DS.vel.zi')
                    ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                    n=r.copy()
                    n[0].data=ntemp/100
                    e=t.copy()
                    e[0].data=etemp/100
                    z[0].data=z[0].data/100
                    n=origin_time(n,time_epi,tb)
                    e=origin_time(e,time_epi,tb)
                    z=origin_time(z,time_epi,tb)
                    n.write(staname[k]+".subfault"+num+'.DS.vel.n',format='SAC')
                    e.write(staname[k]+".subfault"+num+'.DS.vel.e',format='SAC')
                    z.write(staname[k]+".subfault"+num+'.DS.vel.z',format='SAC')
                    silentremove(staname[k]+".subfault"+num+'.DS.vel.r')
                    silentremove(staname[k]+".subfault"+num+'.DS.vel.t')
                    if impulse==True:
                        silentremove(staname[k]+".subfault"+num+'.DS.vel.ri')
                        silentremove(staname[k]+".subfault"+num+'.DS.vel.ti')
                        silentremove(staname[k]+".subfault"+num+'.DS.vel.zi')
            else: #Compute static synthetics
                if okada==False:  #Layered earth model
                    temp_pipe=[]
                    diststr='%.1f' % d[k] #Need current distance in string form for external call
                    green_file=model_name+".static."+strdepth+".sub"+subfault #Output dir
                    statics=loadtxt(green_file) #Load GFs
                    if len(statics)<1:
                        print 'ERROR: Empty GF file'
                        break
                    #Print static GFs into a pipe and pass into synthetics command
                    try:
                        temp_pipe=statics[k,:]
                    except:
                        temp_pipe=statics
                    inpipe=''
                    for j in range(len(temp_pipe)):
                        inpipe=inpipe+' %.6e' % temp_pipe[j]
                    #Form command for external call
                    commandDS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+\
                            " -A"+str(az[k])+" -P"
                    commandSS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+\
                            " -A"+str(az[k])+" -P"
                    commandSS=split(commandSS) #Lexical split
                    commandDS=split(commandDS)
                    #Make system calls, one for DS, one for SS
                    ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
                    p=subprocess.Popen(commandSS,stdin=ps.stdout,stdout=open(staname[k]+'.subfault'+num+'.SS.static.rtz','w'))     
                    p.communicate()  
                    ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
                    p=subprocess.Popen(commandDS,stdin=ps.stdout,stdout=open(staname[k]+'.subfault'+num+'.DS.static.rtz','w'))     
                    p.communicate()       
                    #Rotate radial/transverse to East/North, correct vertical and scale to m
                    statics=loadtxt(staname[k]+'.subfault'+num+'.SS.static.rtz')
                    u=statics[2]/100
                    r=statics[3]/100
                    t=statics[4]/100
                    ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
                    n=ntemp[0]
                    e=etemp[0]
                    savetxt(staname[k]+'.subfault'+num+'.SS.static.neu',(n,e,u,beta))
                    statics=loadtxt(staname[k]+'.subfault'+num+'.DS.static.rtz')
                    u=statics[2]/100
                    r=statics[3]/100
                    t=statics[4]/100
                    ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
                    n=ntemp[0]
                    e=etemp[0]
                    savetxt(staname[k]+'.subfault'+num+'.DS.static.neu',(n,e,u,beta),header='north(m),east(m),up(m),beta(degs)')     
                else: #Okada half space solutions
                #SS
                    n,e,u=okada_synthetics(strike,dip,rakeSS,ss_length_in_km,ds_length_in_km,xs,ys,
                        zs,lon_sta[k],lat_sta[k],mu_okada)
                    savetxt(staname[k]+'.subfault'+num+'.SS.static.neu',(n,e,u,beta),header='north(m),east(m),up(m),beta(degs)')
                    #DS
                    n,e,u=okada_synthetics(strike,dip,rakeDS,ss_length_in_km,ds_length_in_km,xs,ys,
                        zs,lon_sta[k],lat_sta[k],mu_okada)
                    savetxt(staname[k]+'.subfault'+num+'.DS.static.neu',(n,e,u,beta),header='north(m),east(m),up(m),beta(degs)')
     

                      
                                            
                                                                  
def run_parallel_synthetics_mt3d(home,project_name,station_file,model_name,rank,size):
    '''
    Use green functions and compute synthetics at stations for a single source
    and multiple stations. This code makes an external system call to syn.c first it
    will make the external call for the strike-slip component then a second externall
    call will be made for the dip-slip component. The unit amount of moment is 1e15
    which corresponds to Mw=3.9333...
    
    IN:
        source: 1-row numpy array containig informaiton aboutt he source, lat, lon, depth, etc...
        station_file: File name with the station coordinates
        green_path: Directopry where GFs are stored
        model_file: File containing the Earth velocity structure
        integrate: =0 if youw ant velocity waveforms, =1 if you want displacements
        static: =0 if computing full waveforms, =1 if computing only the static field
        subfault: String indicating the subfault being worked on
        coord_type: =0 if problem is in cartesian coordinates, =1 if problem is in lat/lon
        
    OUT:
        log: Sysytem standard output and standard error for log
    '''

    import os
    import subprocess
    from mudpy.forward import get_mu
    from string import rjust
    from numpy import array,genfromtxt,loadtxt,savetxt,log10
    from obspy import read
    from shlex import split
    from mudpy.green import src2sta,rt2ne,origin_time,okada_synthetics
    from glob import glob
    from mudpy.green import silentremove
    from os import remove
    
    #What parameters are we using?
    if rank==0:
        out='''Running all processes with:
        home = %s
        project_name = %s
        station_file = %s
        model_name = %s
        ''' %(home,project_name,station_file,model_name)
        print out
        
    #temporary outoput files to be merged later, these will hold every soruce this process runs
    tmp_Mxx='tmp_Mxx_process'+str(rank)
    tmp_Mxy='tmp_Mxy_process'+str(rank)
    tmp_Mxz='tmp_Mxz_process'+str(rank)
    tmp_Myy='tmp_Myy_process'+str(rank)
    tmp_Myz='tmp_Myz_process'+str(rank)
    tmp_Mzz='tmp_Mzz_process'+str(rank)
    
    #temproary throw away files that will contain only one source
    tmp_small_Mxx='tMxx_proc'+str(rank)
    tmp_small_Mxy='tMxy_proc'+str(rank)
    tmp_small_Mxz='tMxz_proc'+str(rank)
    tmp_small_Myy='tMyy_proc'+str(rank)
    tmp_small_Myz='tMyz_proc'+str(rank)
    tmp_small_Mzz='tMzz_proc'+str(rank)    
        
    #Read your corresponding source file
    mpi_source=genfromtxt(home+project_name+'/data/model_info/mpi_source.'+str(rank)+'.fault')
    
    #Constant parameters
    tb=50 #Number of samples before first arrival (should be 50, NEVER CHANGE, if you do then adjust in fk.pl)
    
    #Load structure
    model_file=home+project_name+'/structure/'+model_name
    structure=loadtxt(model_file,ndmin=2)
    
    #Where the data
    green_path=home+project_name+'/GFs/static/'
    
    #delete files from rpevious runs
    try:
        remove(green_path+tmp_Mxx)
    except:
        pass
    try:
        remove(green_path+tmp_Mxy)
    except:
        pass
    try:
        remove(green_path+tmp_Mxz)
    except:
        pass
    try:
        remove(green_path+tmp_Myy)
    except:
        pass
    try:
        remove(green_path+tmp_Myz)
    except:
        pass
    try:
        remove(green_path+tmp_Mzz)
    except:
        pass
    
    
    #Create output files
    f_Mxx=open(green_path+tmp_Mxx,'a+')
    f_Mxy=open(green_path+tmp_Mxy,'a+')
    f_Mxz=open(green_path+tmp_Mxz,'a+')
    f_Myy=open(green_path+tmp_Myy,'a+')
    f_Myz=open(green_path+tmp_Myz,'a+')
    f_Mzz=open(green_path+tmp_Mzz,'a+')
    
    #Off we go now
    for ksource in range(len(mpi_source)):
        
        source=mpi_source[ksource,:]
        
        #Parse the soruce information
        num=rjust(str(int(source[0])),4,'0')
        xs=source[1]
        ys=source[2]
        zs=source[3]
 
        strdepth='%.4f' % zs
        subfault=rjust(str(int(source[0])),4,'0')
        staname=genfromtxt(station_file,dtype="S6",usecols=0)
        
        if staname.shape==(): #Single staiton file
            staname=array([staname])
        
        #Compute distances and azimuths
        d,az,lon_sta,lat_sta=src2sta(station_file,source,output_coordinates=True)
        
        #Get moment corresponding to 1 meter of slip on subfault
        Mw=5.0
        M0=10**(5.0*1.5+9.1)*1e7 #to dyne-cm
        
        #Move to output folder
        os.chdir(green_path)
        print 'Processor '+str(rank)+' is working on subfault '+str(int(source[0]))+' and '+str(len(d))+' stations '
        
        #Go one station at a time for that subfault
        for k in range(len(d)):
            
            
            temp_pipe=[]
            diststr='%.1f' % d[k] #Need current distance in string form for external call
            green_file=model_name+".static."+strdepth+".sub"+subfault #Output dir
            statics=loadtxt(green_file) #Load GFs
            if len(statics)<1:
                print 'ERROR: Empty GF file'
                break
            #Print static GFs into a pipe and pass into synthetics command
            try:
                temp_pipe=statics[k,:]
            except:
                temp_pipe=statics
            inpipe=''
            for j in range(len(temp_pipe)):
                inpipe=inpipe+' %.6e' % temp_pipe[j]
            
            #Form command for external call (remember syn.c and mudpy coordiante systems are not the same)
            # order of elelments in syn.c is M0/Mxx/Mxy/Mxz/Myy/Myz/Mzz
            Mxx=0 ; Myy=1 ; Mzz=0 ;Mxy=0; Mxz=0; Myz=0
            command_Mxx="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                    " -A"+str(az[k])+" -P"
            
            Mxx=0 ; Myy=0 ; Mzz=0 ;Mxy=1; Mxz=0; Myz=0
            command_Mxy="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                    " -A"+str(az[k])+" -P"
                    
            Mxx=0 ; Myy=0 ; Mzz=0 ;Mxy=0; Mxz=0; Myz=-1
            command_Mxz="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                    " -A"+str(az[k])+" -P"
                    
            Mxx=1 ; Myy=0 ; Mzz=0 ;Mxy=0; Mxz=0; Myz=0
            command_Myy="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                    " -A"+str(az[k])+" -P"
                    
            Mxx=0 ; Myy=0 ; Mzz=0 ;Mxy=0; Mxz=-1; Myz=0
            command_Myz="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                    " -A"+str(az[k])+" -P"
                    
            Mxx=0 ; Myy=0 ; Mzz=-1 ;Mxy=0; Mxz=0; Myz=0
            command_Mzz="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                    " -A"+str(az[k])+" -P"
                    
            command_Mxx=split(command_Mxx) #Lexical split
            command_Mxy=split(command_Mxy)
            command_Mxz=split(command_Mxz)
            command_Myy=split(command_Myy)
            command_Myz=split(command_Myz)
            command_Mzz=split(command_Mzz)
            
            #Make system calls, one for each MT component (rememebr to delete file when youa re done
            ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
            p=subprocess.Popen(command_Mxx,stdin=ps.stdout,stdout=open(tmp_small_Mxx,'w'))     
            p.communicate()  
            p.wait()

            ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
            p=subprocess.Popen(command_Mxy,stdin=ps.stdout,stdout=open(tmp_small_Mxy,'w'))     
            p.communicate()  
            p.wait()

            ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
            p=subprocess.Popen(command_Mxz,stdin=ps.stdout,stdout=open(tmp_small_Mxz,'w'))     
            p.communicate()  
            p.wait()

            ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
            p=subprocess.Popen(command_Myy,stdin=ps.stdout,stdout=open(tmp_small_Myy,'w'))     
            p.communicate()  
            p.wait()

            ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
            p=subprocess.Popen(command_Myz,stdin=ps.stdout,stdout=open(tmp_small_Myz,'w'))     
            p.communicate()  
            p.wait()

            ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
            p=subprocess.Popen(command_Mzz,stdin=ps.stdout,stdout=open(tmp_small_Mzz,'w'))     
            p.communicate()  
            p.wait()
                            
            
            #Rotate radial/transverse to East/North, correct vertical and scale to m
            statics=loadtxt(tmp_small_Mxx)
            u=statics[2]/100
            r=statics[3]/100
            t=statics[4]/100
            ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
            n=ntemp[0]
            e=etemp[0]
            line='%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (rjust(str(subfault),5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
            f_Mxx.write(line)
           
            statics=loadtxt(tmp_small_Mxy)
            u=statics[2]/100
            r=statics[3]/100
            t=statics[4]/100
            ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
            n=ntemp[0]
            e=etemp[0]
            line='%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (rjust(str(subfault),5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
            f_Mxy.write(line)
            
            statics=loadtxt(tmp_small_Mxz)
            u=statics[2]/100
            r=statics[3]/100
            t=statics[4]/100
            ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
            n=ntemp[0]
            e=etemp[0]
            line='%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (rjust(str(subfault),5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
            f_Mxz.write(line)
           
            statics=loadtxt(tmp_small_Myy)
            u=statics[2]/100
            r=statics[3]/100
            t=statics[4]/100
            ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
            n=ntemp[0]
            e=etemp[0]
            line='%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (rjust(str(subfault),5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
            f_Myy.write(line)
           
            statics=loadtxt(tmp_small_Myz)
            u=statics[2]/100
            r=statics[3]/100
            t=statics[4]/100
            ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
            n=ntemp[0]
            e=etemp[0]
            line='%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (rjust(str(subfault),5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
            f_Myz.write(line)
            
            statics=loadtxt(tmp_small_Mzz)
            u=statics[2]/100
            r=statics[3]/100
            t=statics[4]/100
            ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
            n=ntemp[0]
            e=etemp[0]
            line='%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (rjust(str(subfault),5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
            f_Mzz.write(line)
    
                  
    f_Mxx.close()
    f_Mxy.close()
    f_Mxz.close()
    f_Myy.close()
    f_Myz.close()
    f_Mzz.close()        
            
#If main entry point
if __name__ == '__main__':
    import sys
    from mpi4py import MPI
    from obspy import UTCDateTime

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # Map command line arguments to function arguments.
    if sys.argv[1]=='run_parallel_green':
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
        run_parallel_green(home,project_name,station_file,model_name,dt,NFFT,static,dk,pmin,pmax,kmax,tsunami,rank,size)
    elif sys.argv[1]=='run_parallel_synthetics':
        #Parse command line arguments
        home=sys.argv[2]
        project_name=sys.argv[3]
        station_file=sys.argv[4]
        model_name=sys.argv[5]
        integrate=int(sys.argv[6])
        static=int(sys.argv[7])
        tsunami=sys.argv[8]=='True'
        time_epi=UTCDateTime(sys.argv[9])
        beta=float(sys.argv[10])
        custom_stf=sys.argv[11]
        impulse=sys.argv[12]
        if impulse=='True':
            impulse=True
        elif impulse=='False':
            impulse=False
        okada=sys.argv[13]
        if okada=='True':
            okada=True
        elif okada=='False':
            okada=False
        mu_okada=float(sys.argv[14])
        run_parallel_synthetics(home,project_name,station_file,model_name,integrate,static,tsunami,time_epi,beta,custom_stf,impulse,rank,size,okada,mu_okada)
    elif sys.argv[1]=='run_parallel_synthetics_mt3d':
        #Parse command line arguments
        home=sys.argv[2]
        project_name=sys.argv[3]
        station_file=sys.argv[4]
        model_name=sys.argv[5]
        run_parallel_synthetics_mt3d(home,project_name,station_file,model_name,rank,size)
    else:
        print 'ERROR: You''re not allowed to run '+sys.argv[1]+' from the shell or it does not exist'
        
