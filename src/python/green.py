'''
Diego Melgar 01/2014
Green functions routines for kinematic source models
'''

def run_green(source,station_file,model_name,dt,NFFT,static):
    '''
    Compute GFs using Zhu & Rivera code for a given velocity model, source depth
    and station file
    
    params:
        
    returns:
        nothing
    '''
    import subprocess
    from shlex import split
    
    depth=source[3]
    print("--> Computing GFs for source depth "+str(depth)+"km")
    #Get station distances to source
    d=src2sta(station_file,source)
    #Make distance string for system call
    diststr=''
    for k in range(len(d)):
        diststr=diststr+' %.3f' % d[k] #Truncate distance to 3 decimal palces (meters)
    if static==0: #Compute full waveform
        command=split("fk.pl -M"+model_name+"/"+str(depth)+" -N"+str(NFFT)+"/"+str(dt)+diststr)
        print "fk.pl -M"+model_name+"/"+str(depth)+" -N"+str(NFFT)+"/"+str(dt)+diststr
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err=p.communicate() 
    else: #Compute only statics
        command=split("fk.pl -M"+model_name+"/"+str(depth)+" -N1 "+diststr)
        print "fk.pl -M"+model_name+"/"+str(depth)+" -N1 "+diststr
        p=subprocess.Popen(command,stdout=open('staticgf','w'),stderr=subprocess.PIPE)
        out,err=p.communicate() 
    print out
    print err
    
    
def run_syn(source,station_file,green_dir,integrate,static):
    '''
    Use green functions and compute synthetics at stations for a given source
    distribution and station configuration.
    
    params:
        
    returns:
        a sense of well being
    '''

    import os
    import subprocess
    from numpy import array,genfromtxt,rad2deg,loadtxt,savetxt
    from obspy import read
    from obspy.signal.rotate import rotate_RT_NE
    from shlex import split
    
    xs=source[1]
    ys=source[2]
    zs=source[3]
    strike=source[4]
    dip=source[5]
    rake=source[6]
    Mw=source[7]
    rise=source[8]
    duration=source[9]
    strdepth='%.1f' % zs 
    green_path=green_dir+"_"+strdepth+"/"
    print("--> Computing synthetics at stations for the source at ("+str(xs)+" , "+str(ys)+")")
    staname=genfromtxt(station_file,dtype="S6",usecols=0)
    x=genfromtxt(station_file,dtype="f8",usecols=1)
    y=genfromtxt(station_file,dtype="f8",usecols=2)
    #Compute distances and azimuths
    d=src2sta(station_file,source)
    az=cartesian_azimuth(x,y,xs,ys)
    #Move to output folder
    os.chdir(green_path)
    for k in range(len(d)):
        if static==0: #Compute full waveforms
            diststr='%.3f' % d[k] #Need current distance in string form for external call
            if integrate==1: #Make displ.
                command="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rake)+" -D"+str(duration)+ \
                    "/"+str(rise)+" -A"+str(rad2deg(az[k]))+" -O"+staname[k]+".disp.x -G"+green_path+diststr+".grn.0"
                print command
                command=split(command)
            else: #Make vel.
                command="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rake)+" -D"+str(duration)+ \
                    "/"+str(rise)+" -A"+str(rad2deg(az[k]))+" -O"+staname[k]+".vel.x -G"+green_path+diststr+".grn.0"
                print command
                command=split(command)
            p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err=p.communicate() 
            #Rotate to NE
            if integrate==1: #'tis displacememnt
                r=read(staname[k]+'.disp.r')
                t=read(staname[k]+'.disp.t')
                ntemp,etemp=rotate_RT_NE(r[0].data,t[0].data,rad2deg(az[k]))
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                n.write(staname[k]+'.disp.n',format='SAC')
                e.write(staname[k]+'.disp.e',format='SAC')
                #Correct polarity in vertical, ZR code has down=positive, we don't like that precious do we?
                down=read(staname[k]+'.disp.z')
                up=down.copy()
                up[0].data=down[0].data/-100
                up.write(staname[k]+'.disp.z',format='SAC')
            else: #'tis velocity sire
                r=read(staname[k]+'.vel.r')
                t=read(staname[k]+'.vel.t')
                ntemp,etemp=rotate_RT_NE(r[0].data,t[0].data,rad2deg(az[k]))
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                n.write(staname[k]+'.vel.n',format='SAC')
                e.write(staname[k]+'.vel.e',format='SAC')
                #Correct polarity in vertical, ZR code has down=positive, we don't like that precious do we?
                down=read(staname[k]+'.vel.z')
                up=down.copy()
                up[0].data=down[0].data/-100
                up.write(staname[k]+'.vel.z',format='SAC')
        else: #Compute statics
            diststr='%.1f' % d[k] #Need current distance in string form for external call
            green_file=green_dir+".static."+strdepth
            statics=loadtxt(green_file)
            #Print static GFs into a pipe and pass into synthetics command
            temp_pipe=statics[k,:]
            inpipe=''
            for j in range(len(temp_pipe)):
                inpipe=inpipe+' %.6e' % temp_pipe[j]
            command="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rake)+\
                    " -A"+str(rad2deg(az[k]))+" -P"
            print staname[k]
            print command
            command=split(command)
            ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE,stderr=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
            p=subprocess.Popen(command,stdin=ps.stdout,stdout=open(staname[k]+'.static.rtz','w'),stderr=subprocess.PIPE)     
            out,err=p.communicate()       
            #Rotate radial/transverse to East/North
            statics=loadtxt(staname[k]+'.static.rtz')
            u=-statics[2]/100
            r=statics[3]/100
            t=statics[4]/100
            ntemp,etemp=rotate_RT_NE(array([r,r]),array([t,t]),rad2deg(az[k]))
            n=ntemp[0]
            e=etemp[0]
            savetxt(staname[k]+'.static.enu',(e,n,u))


##########                   Utilities and stuff                      ##########          
        
def cartesian_azimuth(x,y,xs,ys):
    '''
    Compute source to station azimuths when sources given in cartesianc oordinates
    '''
    
    from numpy import arctan,nonzero,pi
    
    az=arctan((x-xs)/(y-ys))
    i=nonzero(y<0)
    az[i]=pi+az[i]
    return az
    
def src2sta(station_file,source):
    '''
    Compute source to station distances for all station/source pairs.
    params:
        station_file - Path to station file
        source - numpy 1d array with source info read from file
    returns:
        d - sorted distances vector
    '''
    
    from numpy import genfromtxt,array
    
    #Read station file
    #staname=genfromtxt(home+station_file,dtype="S6",usecols=0)
    x=genfromtxt(station_file,dtype="f8",usecols=1)
    y=genfromtxt(station_file,dtype="f8",usecols=2)
    d=array([])
    xs=source[1]
    ys=source[2]
    d=((x-xs)**2+(y-ys)**2)**0.5
    return d