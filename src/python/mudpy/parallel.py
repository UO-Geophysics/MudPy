'''
Module for routines that use paralell computing
'''


def run_parallel_green(home,project_name,station_file,model_name,dt,NFFT,static,dk,pmin,pmax,kmax,tsunami,insar,rank,size):
    '''
    Compute GFs using Zhu & Rivera code for a given velocity model, source depth
    and station file. This function will make an external system call to fk.pl
    
    IN:
        source: 1-row numpy array containig information about the source, lat, lon, depth, etc...
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
    from numpy import genfromtxt,zeros
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
        insar = %s
        ''' %(home,project_name,station_file,model_name,str(static),str(tsunami),dt,NFFT,dk,pmin,pmax,kmax,str(insar))
        print(out)
    #read your corresponding source file
    source=genfromtxt(home+project_name+'/data/model_info/mpi_source.'+str(rank)+'.fault')
    for ksource in range(len(source)):
        #Where should I be working boss?
        depth='%.4f' % source[ksource,3]
        subfault=str(int(source[ksource,0])).rjust(4,'0')
        if tsunami==False and static==0:
            subfault_folder=home+project_name+'/GFs/dynamic/'+model_name+'_'+depth+'.sub'+subfault
        elif tsunami==True and static==1:
            subfault_folder=home+project_name+'/GFs/tsunami/'+model_name+'_'+depth+'.sub'+subfault
        elif static==1:
            subfault_folder=home+project_name+'/GFs/static/'+model_name+'_'+depth+'.sub'+subfault
        
        #Check if subfault folder exists, if not create it
        if os.path.exists(subfault_folder+'/')==False:
            os.makedirs(subfault_folder+'/')
        
        #Copy velocity model file
        copy(home+project_name+'/structure/'+model_name,subfault_folder+'/'+model_name)
        #Move to work folder
        chdir(subfault_folder)
        #Get station distances to source
        d,az=src2sta(station_file,source[ksource,:])
        #Make distance string for system call
        diststr=''
        for k in range(len(d)):
            diststr=diststr+' %.6f' % d[k] #Truncate distance to 6 decimal palces (meters)
        # Keep the user informed, lest they get nervous
        print('MPI: processor #',rank,'is now working on subfault',int(source[ksource,0]),'(',ksource+1,'/',len(source),')')
        #Make the calculation
        if static==0: #Compute full waveform
            command = "fk.pl -M"+model_name+"/"+depth+"/f -N"+str(NFFT)+"/"+str(dt)+'/1/'+repr(dk)+' -P'+repr(pmin)+'/'+repr(pmax)+'/'+repr(kmax)+diststr
            command=split(command)
            p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            p.communicate() 
            # Move files up one level and delete folder created by fk
            files_list=glob(subfault_folder+'/'+model_name+'_'+depth+'/*.grn*')
            for f in files_list:
                newf=subfault_folder+'/'+f.split('/')[-1]
                copy(f,newf)
            rmtree(subfault_folder+'/'+model_name+'_'+depth)
        else: #Compute only statics
            if insar==True:
                suffix='insar'
            else:
                suffix='gps'
            write_file=subfault_folder+'/'+model_name+'.static.'+depth+'.sub'+subfault+'.'+suffix
            command=split("fk.pl -M"+model_name+"/"+depth+"/f -N1 "+diststr)
            file_is_empty=True
            while file_is_empty:
                p=subprocess.Popen(command,stdout=open(write_file,'w'),stderr=subprocess.PIPE)
                p.communicate()
                if os.stat(write_file).st_size!=0: #File is NOT empty
                    file_is_empty=False
                else:
                    print('Warning: I just had a mini-seizure and made an empty GF file on first try, re-running')
            #If file is empty run again
            

def run_parallel_synthetics(home,project_name,station_file,model_name,integrate,static,quasistatic2dynamic,tsunami,
                            time_epi,beta,custom_stf,impulse,NFFT,dt,rank,size,insar=False,okada=False,mu_okada=45e9,):
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
    from pandas import DataFrame as df
    from mudpy.forward import get_mu
    from numpy import array,genfromtxt,loadtxt,savetxt,log10,zeros,sin,cos,ones,deg2rad
    from obspy import read,Stream,Trace
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
        quasi2dynamic = %s
        time_epi = %s
        beta = %d
        custom_stf = %s
        impulse = %s
        insar = %s
        okada = %s
        mu = %.2e
        ''' %(home,project_name,station_file,model_name,str(integrate),str(static),str(tsunami),str(quasistatic2dynamic),str(time_epi),beta,custom_stf,impulse,insar,okada,mu_okada)
        print(out)
        
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
    
    #this keeps track of statics dataframe
    write_df=False
    
    for ksource in range(len(mpi_source)):
        
        source=mpi_source[ksource,:]
        
        #Parse the soruce information
        num=str(int(source[0])).rjust(4,'0')
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
        subfault=str(int(source[0])).rjust(4,'0')
        
        # print('static is ' +str(static))
        # print('tsunami is ' +str(tsunami))
        # print('quasi is ' +str(quasistatic2dynamic))
        
        if static==0 and tsunami==0 and quasistatic2dynamic==0:  #Where to save dynamic waveforms
            green_path=home+project_name+'/GFs/dynamic/'+model_name+"_"+strdepth+".sub"+subfault+"/"
        if static==1 and tsunami==1 and quasistatic2dynamic==0:  #Where to save dynamic waveforms
            green_path=home+project_name+'/GFs/tsunami/'+model_name+"_"+strdepth+".sub"+subfault+"/"
        if static==1 and tsunami==0 and quasistatic2dynamic==0:  #Where to save statics
            green_path=home+project_name+'/GFs/static/'+model_name+"_"+strdepth+".sub"+subfault+"/"
        if static==0 and quasistatic2dynamic==1:
            #Where to save quasistatic2dynamic "waveforms"
            green_path=home+project_name+'/GFs/dynamic/'+model_name+"_"+strdepth+".sub"+subfault+"/"
            #Where to read the statics from to spoof the waveforms
            read_statics_path=home+project_name+'/GFs/static/'+model_name+"_"+strdepth+".sub"+subfault+"/"
            

            
        staname=genfromtxt(station_file,dtype="U",usecols=0)
        if staname.shape==(): #Single staiton file
            staname=array([staname])
        #Compute distances and azimuths
        d,az,lon_sta,lat_sta=src2sta(station_file,source,output_coordinates=True)
        
        #Get moment corresponding to 1 meter of slip on subfault
        mu=get_mu(structure,zs)
        Mo=mu*ss_length*ds_length*1.0
        Mw=(2./3)*(log10(Mo)-9.1)
        
        #Move to output folder if it doesn't exist create it
        #Check fist if folder exists
        dir_exists = os.path.exists(green_path)
        if dir_exists: #no need to do anything
            pass
        else: #This only happens in quasistatic2dyanmic case
            os.makedirs(green_path)
        
        #ok now move there
        os.chdir(green_path)
        print('Processor '+str(rank)+' is working on subfault '+str(int(source[0]))+' and '+str(len(d))+' stations ')
        
        
        #If we are doibng quasistatic2dynamic we only need to read the syntehtics file once per source
        if static==0 and quasistatic2dynamic==1 : #Convert statics to ramp waveforms
                
                #read alls tations and statics
                fileSS = 'subfault'+subfault +'.SS.static.neu'
                fileDS = 'subfault'+subfault +'.DS.static.neu'
  
                station_names = genfromtxt(read_statics_path+fileSS,usecols=0,dtype='U')
                SSstatics = genfromtxt(read_statics_path+fileSS)
                DSstatics = genfromtxt(read_statics_path+fileDS)
        
        #This is looping over "sites"
        for k in range(len(d)):
            
            if static==0 and quasistatic2dynamic==0: #Compute full waveforms
                diststr='%.6f' % d[k] #Need current distance in string form for external call
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
                    # get rid of numerical "noise" in the first tb samples
                    n[0].data[0:tb]=0
                    e[0].data[0:tb]=0
                    z[0].data[0:tb]=0
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
            
            elif static==0 and quasistatic2dynamic==1 : #Convert statics to ramp waveforms
            
                #site name
                sta_name = station_names[k]
                
                #offsets
                nSS = SSstatics[k,1] * ones(NFFT)
                eSS = SSstatics[k,2] * ones(NFFT)
                zSS = SSstatics[k,3] * ones(NFFT)
                
                nDS = DSstatics[k,1] * ones(NFFT)
                eDS = DSstatics[k,2] * ones(NFFT)
                zDS = DSstatics[k,3] * ones(NFFT)
                
                #initalize streams
                st_nSS = Stream(Trace()) ; st_nDS = Stream(Trace())
                st_eSS = Stream(Trace()) ; st_eDS = Stream(Trace())
                st_zSS = Stream(Trace()) ; st_zDS = Stream(Trace())
                
                #assign data and metadata
                st_nSS[0].data = nSS
                st_nSS[0].stats.delta = dt
                st_nSS[0].stats.starttime = time_epi
                
                st_eSS[0].data = eSS
                st_eSS[0].stats.delta = dt
                st_eSS[0].stats.starttime = time_epi
                
                st_zSS[0].data = zSS
                st_zSS[0].stats.delta = dt
                st_zSS[0].stats.starttime = time_epi
                
                st_nDS[0].data = nDS
                st_nDS[0].stats.delta = dt
                st_nDS[0].stats.starttime = time_epi
                
                st_eDS[0].data = eDS
                st_eDS[0].stats.delta = dt
                st_eDS[0].stats.starttime = time_epi
                
                st_zDS[0].data = zDS
                st_zDS[0].stats.delta = dt
                st_zDS[0].stats.starttime = time_epi
                
                #Write files
                ss_file_out = sta_name+'.subfault'+subfault+'.SS.disp.'
                ds_file_out = sta_name+'.subfault'+subfault+'.DS.disp.'
                
                st_nSS.write(ss_file_out+'n',format='SAC')
                st_eSS.write(ss_file_out+'e',format='SAC')
                st_zSS.write(ss_file_out+'z',format='SAC')
                
                st_nDS.write(ds_file_out+'n',format='SAC')
                st_eDS.write(ds_file_out+'e',format='SAC')
                st_zDS.write(ds_file_out+'z',format='SAC')
            
            
            else: #Compute static synthetics
                if okada==False:  #Layered earth model
                    
                    #this is because when I first wrote this code it processed each
                    #source/station pair independently but now that it's vectorized
                    #it's does ALL stations in one fell swoop, given the logic it's 
                    #easier to just keep this inside the for loop and use the if to
                    #run it jsut the first time for all sites
                    if k==0:   
                        
                        #initalize output variables
                        staticsSS = zeros((len(d),4))
                        staticsDS = zeros((len(d),4))
                        write_df=True
                        
                        #read GFs file
                        if insar==True:
                            green_file=green_path+model_name+".static."+strdepth+".sub"+subfault+'.insar' #Output dir
                        else: #GPS
                            green_file=green_path+model_name+".static."+strdepth+".sub"+subfault+'.gps' #Output 

                        statics=loadtxt(green_file) #Load GFs
                        Nsites=len(statics) 
                        
                        if len(statics)<1:
                            print('ERROR: Empty GF file')
                            break
                        
                        #Now get radiation pattern terms, there will be 3 terms
                        #for each direction so 9 terms total. THis comes from funtion
                        #dc_radiat() in radiats.c from fk
                        radiation_pattern_ss = zeros((Nsites,9))
                        radiation_pattern_ds = zeros((Nsites,9))
                        
                        rakeSSrad = deg2rad(rakeSS)
                        rakeDSrad = deg2rad(rakeDS)
                        dip_rad = deg2rad(dip)
                        pseudo_strike = deg2rad(az-strike)
                        
                        #Let's do SS first
                        r = rakeSSrad
                        
                        #trigonometric terms following nomenclature used in radiats.c
                        sstk=sin(pseudo_strike) ; cstk=cos(pseudo_strike)
                        sdip=sin(dip_rad) ; cdip=cos(dip_rad)
                        srak=sin(r) ; crak=cos(r)
                        sstk2=2*sstk*cstk ; cstk2=cstk*cstk-sstk*sstk
                        sdip2=2*sdip*cdip ; cdip2=cdip*cdip-sdip*sdip
                        
                        # terms for up component
                        u_dd = 0.5*srak*sdip2*ones(Nsites)
                        u_ds = -sstk*srak*cdip2+cstk*crak*cdip
                        u_ss = -sstk2*crak*sdip-0.5*cstk2*srak*sdip2
                        
                        #terms for r component
                        r_dd = u_dd.copy()
                        r_ds = u_ds.copy()
                        r_ss = u_ss.copy()
                        
                        #terms for t component
                        t_dd = zeros(Nsites)
                        t_ds = cstk*srak*cdip2+sstk*crak*cdip
                        t_ss = cstk2*crak*sdip-0.5*sstk2*srak*sdip2

                        #assemble in one variable
                        radiation_pattern_ss[:,0] = u_dd
                        radiation_pattern_ss[:,1] = u_ds
                        radiation_pattern_ss[:,2] = u_ss
                        
                        radiation_pattern_ss[:,3] = r_dd
                        radiation_pattern_ss[:,4] = r_ds
                        radiation_pattern_ss[:,5] = r_ss
                        
                        radiation_pattern_ss[:,6] = t_dd
                        radiation_pattern_ss[:,7] = t_ds
                        radiation_pattern_ss[:,8] = t_ss
                        
                        
                        #Now radiation pattern for DS
                        r = rakeDSrad
                        
                        #trigonometric terms following nomenclature used in radiats.c
                        sstk=sin(pseudo_strike) ; cstk=cos(pseudo_strike)
                        sdip=sin(dip_rad) ; cdip=cos(dip_rad)
                        srak=sin(r) ; crak=cos(r)
                        sstk2=2*sstk*cstk ; cstk2=cstk*cstk-sstk*sstk
                        sdip2=2*sdip*cdip ; cdip2=cdip*cdip-sdip*sdip
                        
                        # terms for up component
                        u_dd = 0.5*srak*sdip2*ones(Nsites)
                        u_ds = -sstk*srak*cdip2+cstk*crak*cdip
                        u_ss = -sstk2*crak*sdip-0.5*cstk2*srak*sdip2
                        
                        #terms for r component
                        r_dd = u_dd.copy()
                        r_ds = u_ds.copy()
                        r_ss = u_ss.copy()
                        
                        #terms for t component
                        t_dd = zeros(Nsites)
                        t_ds = cstk*srak*cdip2+sstk*crak*cdip
                        t_ss = cstk2*crak*sdip-0.5*sstk2*srak*sdip2
                        
                        #assemble in one variable
                        radiation_pattern_ds[:,0] = u_dd
                        radiation_pattern_ds[:,1] = u_ds
                        radiation_pattern_ds[:,2] = u_ss
                        
                        radiation_pattern_ds[:,3] = r_dd
                        radiation_pattern_ds[:,4] = r_ds
                        radiation_pattern_ds[:,5] = r_ss
                        
                        radiation_pattern_ds[:,6] = t_dd
                        radiation_pattern_ds[:,7] = t_ds
                        radiation_pattern_ds[:,8] = t_ss
                        
                        
                        #Now define the scalng based on magnitude this is variable
                        #"coef" in the syn.c original source code
                        scale = 10**(1.5*Mw+16.1-20) #definition used in syn.c
                        
                        #Scale radiation patterns accordingly
                        radiation_pattern_ss *= scale
                        radiation_pattern_ds *= scale
                        
                        #Now multiply each GF component by the appropriate SCALED
                        #radiation pattern term and add em up to get the displacements
                        # also /100 to convert  to meters
                        up_ss = radiation_pattern_ss[:,0:3]*statics[:,[1,4,7]] 
                        up_ss = up_ss.sum(axis=1) / 100
                        up_ds = radiation_pattern_ds[:,0:3]*statics[:,[1,4,7]] 
                        up_ds = up_ds.sum(axis=1)    / 100                     

                        radial_ss = radiation_pattern_ss[:,3:6]*statics[:,[2,5,8]] 
                        radial_ss = radial_ss.sum(axis=1) / 100
                        radial_ds = radiation_pattern_ds[:,3:6]*statics[:,[2,5,8]] 
                        radial_ds = radial_ds.sum(axis=1) / 100 
                        
                        tangential_ss = radiation_pattern_ss[:,6:9]*statics[:,[3,6,9]] 
                        tangential_ss = tangential_ss.sum(axis=1) / 100
                        tangential_ds = radiation_pattern_ds[:,6:9]*statics[:,[3,6,9]] 
                        tangential_ds = tangential_ds.sum(axis=1) / 100
                        
                        #rotate to neu
                        n_ss,e_ss=rt2ne(radial_ss,tangential_ss,az)
                        n_ds,e_ds=rt2ne(radial_ds,tangential_ds,az)

                        #put in output variables
                        staticsSS[:,0]=n_ss
                        staticsSS[:,1]=e_ss
                        staticsSS[:,2]=up_ss
                        staticsSS[:,3]=beta*ones(Nsites)
                        
                        staticsDS[:,0]=n_ds
                        staticsDS[:,1]=e_ds
                        staticsDS[:,2]=up_ds
                        staticsDS[:,3]=beta*ones(Nsites)
                        
                    else:
                        pass
            

                else: #Okada half space solutions
                #SS
                    n,e,u=okada_synthetics(strike,dip,rakeSS,ss_length_in_km,ds_length_in_km,xs,ys,
                        zs,lon_sta[k],lat_sta[k],mu_okada)
                    savetxt(staname[k]+'.subfault'+num+'.SS.static.neu',(n,e,u,beta),header='north(m),east(m),up(m),beta(degs)')
                    #DS
                    n,e,u=okada_synthetics(strike,dip,rakeDS,ss_length_in_km,ds_length_in_km,xs,ys,
                        zs,lon_sta[k],lat_sta[k],mu_okada)
                    savetxt(staname[k]+'.subfault'+num+'.DS.static.neu',(n,e,u,beta),header='north(m),east(m),up(m),beta(degs)')

        
        if write_df==True and static ==1: #Note to self: stop using 0,1 and swithc to True/False
            
            #Strike slip
            SSdf = df(data=None, index=None, columns=['staname','n','e','u','beta'])
            SSdf.staname=staname
            SSdf.n=staticsSS[:,0]
            SSdf.e=staticsSS[:,1]
            SSdf.u=staticsSS[:,2]
            SSdf.beta=staticsSS[:,3]
            if insar == False: # GNSS statics file
                SSdf.to_csv(green_path+'subfault'+num+'.SS.static.neu',sep='\t',index=False,header=False)
            else:  # InSAR file
                SSdf.to_csv(green_path+'subfault'+num+'.SS.insar.neu',sep='\t',index=False,header=False)
            
            DSdf = df(data=None, index=None, columns=['staname','n','e','u','beta'])     
            DSdf.staname=staname
            DSdf.n=staticsDS[:,0]
            DSdf.e=staticsDS[:,1]
            DSdf.u=staticsDS[:,2]
            DSdf.beta=staticsDS[:,3]
            if insar == False: # GNSS statics file
                DSdf.to_csv(green_path+'subfault'+num+'.DS.static.neu',sep='\t',index=False,header=False)
            else:
                DSdf.to_csv(green_path+'subfault'+num+'.DS.insar.neu',sep='\t',index=False,header=False)

                      




def run_parallel_teleseismics_green(home,project_name,time_epi,station_file,model_name,teleseismic_vel_mod,endtime,rank,size):
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
    import requests
    from mudpy.forward import get_mu
    from numpy import array,genfromtxt
    from obspy import read
    from mudpy.green import src2sta
    
    #What parameters are we using?
    if rank==0:
        out='''Running all processes with:
        home = %s
        project_name = %s
        station_file = %s
        model_name = %s
        time_epi = %s
        end_time = %s
        ''' %(home,project_name,station_file,model_name,str(time_epi),str(endtime))
        print(out)
        
    #url for web request
    url='http://service.iris.edu/irisws/syngine/1/query'
        
    #Read your corresponding source file
    mpi_source=genfromtxt(home+project_name+'/data/model_info/mpi_source.'+str(rank)+'.fault')
    
    #Constant parameters
    rakeDS=90 #90 is thrust, -90 is normal
    rakeSS=0 #0 is left lateral, 180 is right lateral
    
    #Load structure
    model_file=home+project_name+'/structure/'+model_name
    structure=genfromtxt(model_file)
    
    for ksource in range(len(mpi_source)):
        
        source=mpi_source[ksource,:]
        
        #Parse the soruce information
        xs=source[1]
        ys=source[2]
        zs=source[3]
        zs_in_meters=int(zs*1000)
        strike=source[4]
        dip=source[5]
            
        #check longitude because iris only allows +/-180
        if xs > 180:
            xs -= 360
        
        ss_length=source[8]
        ds_length=source[9]
        
        strdepth='%.4f' % zs
        subfault=str(int(source[0])).rjust(4,'0')
        
        #Where to save dynamic waveforms
        green_path=home+project_name+'/GFs/dynamic/'+model_name+"_"+strdepth+".sub"+subfault+"/"
       
        #check if folder exists if not make it
        if os.path.isdir(green_path) == False:
            os.mkdir(green_path)
       
        staname=genfromtxt(station_file,dtype="U",usecols=0)
        if staname.shape==(): #Single staiton file
            staname=array([staname])
        
        #Compute distances and azimuths
        d,az,lon_sta,lat_sta=src2sta(station_file,source,output_coordinates=True)
        
        #Get moment corresponding to 1 meter of slip on subfault
        mu=get_mu(structure,zs)
        Mo=mu*ss_length*ds_length*1.0
        
        #Move to output folder
        os.chdir(green_path)
        print('Processor '+str(rank)+' is working on subfault '+str(int(source[0]))+' and '+str(len(d))+' stations ')
        
        #This is looping over "sites"
        for ksta in range(len(d)):
            
            #cehck longitude for valid range
            if lon_sta[ksta] > 180:
                lon_sta[ksta] -= 360
            
            #Form web request for SS syntehtic
            parameters = {'model': teleseismic_vel_mod,
                         'sourcelatitude':str(ys),
                         'sourcelongitude':str(xs),
                         'sourcedepthinmeters':str(zs_in_meters),
                         'sourcedoublecouple':str(strike)+','+str(dip)+','+str(rakeSS)+','+str(Mo),
                         'receiverlatitude':str(lat_sta[ksta]),
                         'receiverlongitude':str(lon_sta[ksta]),
                         'format':'miniseed',
                         'origintime':str(time_epi),
                         'endtime':str(endtime)}
                
            web_request = requests.get(url, params = parameters)
            
            #Filename for temporary write
            out=green_path+staname[ksta]+'.subfault'+subfault+'.SS.mseed'
            f=open(out,'wb')
            f.write(web_request.content)
            
            
            #Form web request for DS syntehtic
            parameters = {'model': teleseismic_vel_mod,
                         'sourcelatitude':str(ys),
                         'sourcelongitude':str(xs),
                         'sourcedepthinmeters':str(zs_in_meters),
                         'sourcedoublecouple':str(strike)+','+str(dip)+','+str(rakeDS)+','+str(Mo),
                         'receiverlatitude':str(lat_sta[ksta]),
                         'receiverlongitude':str(lon_sta[ksta]),
                         'format':'miniseed',
                         'origintime':str(time_epi),
                         'endtime':str(endtime)}
                
            web_request = requests.get(url, params = parameters)
            
            #Filename for temporary write
            out=green_path+staname[ksta]+'.subfault'+subfault+'.DS.mseed'
            f=open(out,'wb')
            f.write(web_request.content)
        






                                            
                                                                  
def run_parallel_synthetics_mt3d(home,project_name,station_file,model_name,forceMT,mt,insar,rank,size):
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
        forceMT = %s
        mt = %s
        insar = %s
        ''' %(home,project_name,station_file,model_name,str(forceMT),str(mt),str(insar))
        print(out)
        
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
    
    #Make moment tensor components
    if forceMT==True:
        Mxx=mt[0] ; Mxy=mt[1] ; Mxz=mt[2] ;Myy=mt[3]; Myz=mt[4]; Mzz=mt[5]
    
    #Off we go now
    for ksource in range(len(mpi_source)):
        
        source=mpi_source[ksource,:]
        
        #Parse the soruce information
        num=str(int(source[0])).rjust(4,'0')
        xs=source[1]
        ys=source[2]
        zs=source[3]
 
        strdepth='%.4f' % zs
        subfault=str(int(source[0])).rjust(4,'0')
        staname=genfromtxt(station_file,dtype="U",usecols=0)
        
        if staname.shape==(): #Single staiton file
            staname=array([staname])
        
        #Compute distances and azimuths
        d,az,lon_sta,lat_sta=src2sta(station_file,source,output_coordinates=True)
        
        #Get moment corresponding to 1 meter of slip on subfault
        Mw=5.0
        M0=10**(5.0*1.5+9.1)*1e7 #to dyne-cm
        
        #Load LOS vector for projection
        los_path=home+project_name+'/data/statics/'
        
        #Move to output folder
        os.chdir(green_path)
        print('Processor '+str(rank)+' is working on subfault '+str(int(source[0]))+' and '+str(len(d))+' stations ')
        
        #Go one station at a time for that subfault
        for k in range(len(d)):
            
            #Read los vector for this subfault
            if insar==True:
                los=genfromtxt(los_path+staname[k]+'.los')
                los=los[1:]
            
            # Load the GFs
            if insar==False:
                green_file=model_name+".static."+strdepth+".sub"+subfault+'.gps' #Output dir
            else:
                green_file=model_name+".static."+strdepth+".sub"+subfault+'.insar' #Output dir
            statics=loadtxt(green_file) #Load GFs
            if len(statics)<1:
                print('ERROR: Empty GF file')
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
            
            if forceMT==True: #Only run one thing
                command_Mxx="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                    " -A"+str(az[k])+" -P"  
                    
                command_Mxx=split(command_Mxx) #Lexical split
            
                #Make system calls, one for each MT component (rememebr to delete file when youa re done
                ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
                p=subprocess.Popen(command_Mxx,stdin=ps.stdout,stdout=open(tmp_small_Mxx,'w'))     
                p.communicate()  
                p.wait()        
                
                #Rotate radial/transverse to East/North, correct vertical and scale to m
                statics=loadtxt(tmp_small_Mxx)
                u=statics[2]/100
                r=statics[3]/100
                t=statics[4]/100
                ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
                #now project onto LOS if doing insar
                if insar==True:
                    los_mt=los.dot(array([ntemp[0],etemp[0],u]))
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,los_mt)
                else:
                    n=ntemp[0]
                    e=etemp[0]
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
                #write to file
                f_Mxx.write(line)      
            
            else: #Stuff to dow he computing entire MT
                Mxx=1 ; Myy=0 ; Mzz=0 ;Mxy=0; Mxz=0; Myz=0
                command_Mxx="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                        " -A"+str(az[k])+" -P"
                
                Mxx=0 ; Myy=0 ; Mzz=0 ;Mxy=1; Mxz=0; Myz=0
                command_Mxy="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                        " -A"+str(az[k])+" -P"
                        
                Mxx=0 ; Myy=0 ; Mzz=0 ;Mxy=0; Mxz=1; Myz=0
                command_Mxz="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                        " -A"+str(az[k])+" -P"
                       
                Mxx=0 ; Myy=1 ; Mzz=0 ;Mxy=0; Mxz=0; Myz=0
                command_Myy="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                        " -A"+str(az[k])+" -P"
                        
                Mxx=0 ; Myy=0 ; Mzz=0 ;Mxy=0; Mxz=0; Myz=1
                command_Myz="syn -M"+str(M0)+"/"+str(Mxx)+"/"+str(Mxy)+"/"+str(Mxz)+"/"+str(Myy)+"/"+str(Myz)+"/"+str(Mzz)+\
                        " -A"+str(az[k])+" -P"
                        
                Mxx=0 ; Myy=0 ; Mzz=1 ;Mxy=0; Mxz=0; Myz=0
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
                #now project onto LOS if doing insar
                if insar==True:
                    los_mt=los.dot(array([ntemp[0],etemp[0],u]))
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,los_mt)
                else:
                    n=ntemp[0]
                    e=etemp[0]
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
                #write to file
                f_Mxx.write(line)
            
                statics=loadtxt(tmp_small_Mxy)
                u=statics[2]/100
                r=statics[3]/100
                t=statics[4]/100
                ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
                n=ntemp[0]
                e=etemp[0]
                #now project onto LOS if doing insar
                if insar==True:
                    los_mt=los.dot(array([ntemp[0],etemp[0],u]))
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,los_mt)
                else:
                    n=ntemp[0]
                    e=etemp[0]
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
                #write to file
                f_Mxy.write(line)
                
                statics=loadtxt(tmp_small_Mxz)
                u=statics[2]/100
                r=statics[3]/100
                t=statics[4]/100
                ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
                n=ntemp[0]
                e=etemp[0]
                #now project onto LOS if doing insar
                if insar==True:
                    los_mt=los.dot(array([ntemp[0],etemp[0],u]))
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,los_mt)
                else:
                    n=ntemp[0]
                    e=etemp[0]
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
                #write to file
                f_Mxz.write(line)
            
                statics=loadtxt(tmp_small_Myy)
                u=statics[2]/100
                r=statics[3]/100
                t=statics[4]/100
                ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
                n=ntemp[0]
                e=etemp[0]
                #now project onto LOS if doing insar
                if insar==True:
                    los_mt=los.dot(array([ntemp[0],etemp[0],u]))
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,los_mt)
                else:
                    n=ntemp[0]
                    e=etemp[0]
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
                #write to file
                f_Myy.write(line)
            
                statics=loadtxt(tmp_small_Myz)
                u=statics[2]/100
                r=statics[3]/100
                t=statics[4]/100
                ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
                n=ntemp[0]
                e=etemp[0]
                #now project onto LOS if doing insar
                if insar==True:
                    los_mt=los.dot(array([ntemp[0],etemp[0],u]))
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,los_mt)
                else:
                    n=ntemp[0]
                    e=etemp[0]
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
                #write to file
                f_Myz.write(line)
                
                statics=loadtxt(tmp_small_Mzz)
                u=statics[2]/100
                r=statics[3]/100
                t=statics[4]/100
                ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
                n=ntemp[0]
                e=etemp[0]
                #now project onto LOS if doing insar
                if insar==True:
                    los_mt=los.dot(array([ntemp[0],etemp[0],u]))
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,los_mt)
                else:
                    n=ntemp[0]
                    e=etemp[0]
                    line='%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (staname[k],str(subfault).rjust(5,'0'),lon_sta[k],lat_sta[k],xs,ys,zs,n,e,u)
                #write to file
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
        insar=sys.argv[14]
        if insar=='True':
            insar=True
        elif insar=='False':
            insar=False
        run_parallel_green(home,project_name,station_file,model_name,dt,NFFT,static,dk,pmin,pmax,kmax,tsunami,insar,rank,size)

    
    elif sys.argv[1]=='run_parallel_synthetics':
        #Parse command line arguments
        home=sys.argv[2]
        project_name=sys.argv[3]
        station_file=sys.argv[4]
        model_name=sys.argv[5]
        integrate=int(sys.argv[6])
        static=int(sys.argv[7])
        quasistatic2dynamic=int(sys.argv[8])
        tsunami=sys.argv[9]=='True'
        time_epi=UTCDateTime(sys.argv[10])
        beta=float(sys.argv[11])
        custom_stf=sys.argv[12]
        impulse=sys.argv[13]
        if impulse=='True':
            impulse=True
        elif impulse=='False':
            impulse=False
        insar=sys.argv[14]
        if insar=='True':
            insar=True
        elif insar=='False':
            insar=False
        okada=sys.argv[15]
        if okada=='True':
            okada=True
        elif okada=='False':
            okada=False
        mu_okada=float(sys.argv[16])
        NFFT = int(sys.argv[17])
        dt = float(sys.argv[18])

        run_parallel_synthetics(home,project_name,station_file,model_name,integrate,static,quasistatic2dynamic,tsunami,time_epi,beta,custom_stf,impulse,NFFT,dt,rank,size,insar,okada,mu_okada)
    
    elif sys.argv[1]=='run_parallel_teleseismics_green':
        home=sys.argv[2]
        project_name=sys.argv[3]
        time_epi=UTCDateTime(sys.argv[4])
        station_file=sys.argv[5]
        model_name=sys.argv[6]
        teleseismic_vel_mod=sys.argv[7]
        endtime=sys.argv[8]
    
        run_parallel_teleseismics_green(home,project_name,time_epi,station_file,model_name,teleseismic_vel_mod,endtime,rank,size)
    
    elif sys.argv[1]=='run_parallel_synthetics_mt3d':
        #Parse command line arguments
        home=sys.argv[2]
        project_name=sys.argv[3]
        station_file=sys.argv[4]
        model_name=sys.argv[5]
        forceMT=sys.argv[6]
        if forceMT=='True':
            forceMT=True
        elif forceMT=='False':
            forceMT=False
        Mxx=float(sys.argv[7])
        Mxy=float(sys.argv[8])
        Mxz=float(sys.argv[9])
        Myy=float(sys.argv[10])
        Myz=float(sys.argv[11])
        Mzz=float(sys.argv[12])
        mt=[Mxx,Mxy,Mxz,Myy,Myz,Mzz]
        
        insar=sys.argv[13]
        if insar=='True':
            insar=True
        elif insar=='False':
            insar=False
        run_parallel_synthetics_mt3d(home,project_name,station_file,model_name,forceMT,mt,insar,rank,size)
    else:
        print('ERROR: You''re not allowed to run '+sys.argv[1]+' from the shell or it does not exist')
        