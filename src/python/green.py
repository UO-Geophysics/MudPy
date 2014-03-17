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
    
    depth='%.4f' % source[3]
    print("--> Computing GFs for source depth "+str(depth)+"km")
    #Get station distances to source
    d=src2sta(station_file,source)
    #Make distance string for system call
    diststr=''
    for k in range(len(d)):
        diststr=diststr+' %.3f' % d[k] #Truncate distance to 3 decimal palces (meters)
    if static==0: #Compute full waveform
        command=split("fk.pl -M"+model_name+"/"+depth+" -N"+str(NFFT)+"/"+str(dt)+diststr)
        print "fk.pl -M"+model_name+"/"+depth+" -N"+str(NFFT)+"/"+str(dt)+diststr
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err=p.communicate() 
    else: #Compute only statics
        command=split("fk.pl -M"+model_name+"/"+depth+" -N1 "+diststr)
        print "fk.pl -M"+model_name+"/"+depth+" -N1 "+diststr
        p=subprocess.Popen(command,stdout=open('staticgf','w'),stderr=subprocess.PIPE)
        out,err=p.communicate() 
    print out
    print err
    log=str(out)+str(err)
    return log
    
    
    
def run_syn(source,station_file,green_path,model_name,integrate,static,subfault):
    '''
    Use green functions and compute synthetics at stations for a given source
    distribution and station configuration.
    
    params:
        
    returns:
        a sense of well being
    '''

    import os
    import subprocess
    from string import rjust
    from numpy import array,genfromtxt,rad2deg,loadtxt,savetxt
    from obspy import read
    from obspy.signal.rotate import rotate_RT_NE
    from shlex import split
    
    #Constant parameters
    rakeDS=-90 #-90 is thrust, 90 is normal
    rakeSS=180
    Mw=3.933333333 #This is used as unit magnitude, corresponds to 1e15 N-m
    
    num=rjust(str(int(source[0])),4,'0')
    xs=source[1]
    ys=source[2]
    zs=source[3]
    strike=source[4]
    dip=source[5]
    rise=source[6]
    duration=source[7]

    strdepth='%.4f' % zs
    if static==0: 
        green_path=green_path+'dynamic/'+model_name+"_"+strdepth+".sub"+subfault+"/"
    print("--> Computing synthetics at stations for the source at ("+str(xs)+" , "+str(ys)+")")
    staname=genfromtxt(station_file,dtype="S6",usecols=0)
    x=genfromtxt(station_file,dtype="f8",usecols=1)
    y=genfromtxt(station_file,dtype="f8",usecols=2)
    #Compute distances and azimuths
    d=src2sta(station_file,source)
    az=cartesian_azimuth(x,y,xs,ys)
    #Move to output folder
    log=''
    os.chdir(green_path)
    for k in range(len(d)):
        if static==0: #Compute full waveforms
            diststr='%.3f' % d[k] #Need current distance in string form for external call
            if integrate==1: #Make displ.
                #First Stike-Slip GFs
                commandSS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -D"+str(duration)+ \
                    "/"+str(rise)+" -A"+str(rad2deg(az[k]))+" -O"+staname[k]+".subfault"+num+".SS.disp.x -G"+green_path+diststr+".grn.0"
                print commandSS
                log=log+commandSS+'\n'
                commandSS=split(commandSS)
                #Now dip slip
                commandDS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -D"+str(duration)+ \
                    "/"+str(rise)+" -A"+str(rad2deg(az[k]))+" -O"+staname[k]+".subfault"+num+".DS.disp.x -G"+green_path+diststr+".grn.0"
                print commandDS
                log=log+commandDS+'\n'
                commandDS=split(commandDS)
            else: #Make vel.
                #First Stike-Slip GFs
                commandSS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -D"+str(duration)+ \
                    "/"+str(rise)+" -A"+str(rad2deg(az[k]))+" -O"+staname[k]+".subfault"+num+".SS.vel.x -G"+green_path+diststr+".grn.0"
                print commandSS
                log=log+commandSS+'\n'
                commandSS=split(commandSS)
                #Now dip slip
                commandDS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -D"+str(duration)+ \
                    "/"+str(rise)+" -A"+str(rad2deg(az[k]))+" -O"+staname[k]+".subfault"+num+".DS.vel.x -G"+green_path+diststr+".grn.0"
                print commandDS
                log=log+commandDS+'\n'
                commandDS=split(commandDS)
            #Run them yo
            p=subprocess.Popen(commandSS,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err=p.communicate() 
            p=subprocess.Popen(commandDS,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err=p.communicate()
            #Rotate to NE
            if integrate==1: #'tis displacememnt
                #Strike slip
                r=read(staname[k]+".subfault"+num+'.SS.disp.r')
                t=read(staname[k]+".subfault"+num+'.SS.disp.t')
                ntemp,etemp=rotate_RT_NE(r[0].data,t[0].data,rad2deg(az[k]))
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                n.write(staname[k]+".subfault"+num+'.SS.disp.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.SS.disp.e',format='SAC')
                #Correct polarity in vertical, ZR code has down=positive, we don't like that precious do we?
                down=read(staname[k]+".subfault"+num+'.SS.disp.z')
                up=down.copy()
                up[0].data=down[0].data/-100
                up.write(staname[k]+".subfault"+num+'.SS.disp.z',format='SAC')
                #Dip Slip
                r=read(staname[k]+".subfault"+num+'.DS.disp.r')
                t=read(staname[k]+".subfault"+num+'.DS.disp.t')
                ntemp,etemp=rotate_RT_NE(r[0].data,t[0].data,rad2deg(az[k]))
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                n.write(staname[k]+".subfault"+num+'.DS.disp.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.DS.disp.e',format='SAC')
                #Correct polarity in vertical, ZR code has down=positive, we don't like that precious do we?
                down=read(staname[k]+".subfault"+num+'.DS.disp.z')
                up=down.copy()
                up[0].data=down[0].data/-100
                up.write(staname[k]+".subfault"+num+'.DS.disp.z',format='SAC')
            else: #'tis velocity sire
                #Strike slip
                r=read(staname[k]+".subfault"+num+'.SS.vel.r')
                t=read(staname[k]+".subfault"+num+'.SS.vel.t')
                ntemp,etemp=rotate_RT_NE(r[0].data,t[0].data,rad2deg(az[k]))
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                n.write(staname[k]+".subfault"+num+'.SS.vel.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.SS.vel.e',format='SAC')
                #Correct polarity in vertical, ZR code has down=positive, we don't like that precious do we?
                down=read(staname[k]+".subfault"+num+'.SS.vel.z')
                up=down.copy()
                up[0].data=down[0].data/-100
                up.write(staname[k]+".subfault"+num+'.SS.vel.z',format='SAC')
                #Dip Slip
                r=read(staname[k]+".subfault"+num+'.DS.vel.r')
                t=read(staname[k]+".subfault"+num+'.DS.vel.t')
                ntemp,etemp=rotate_RT_NE(r[0].data,t[0].data,rad2deg(az[k]))
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                n.write(staname[k]+".subfault"+num+'.DS.vel.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.DS.vel.e',format='SAC')
                #Correct polarity in vertical, ZR code has down=positive, we don't like that precious do we?
                down=read(staname[k]+".subfault"+num+'.DS.vel.z')
                up=down.copy()
                up[0].data=down[0].data/-100
                up.write(staname[k]+".subfault"+num+'.DS.vel.z',format='SAC')
        else: #Compute statics
            os.chdir(green_path+'static/')
            diststr='%.1f' % d[k] #Need current distance in string form for external call
            green_file=model_name+".static."+strdepth+".sub"+subfault
            print green_file
            log=log+green_file+'\n'
            statics=loadtxt(green_file)
            #Print static GFs into a pipe and pass into synthetics command
            temp_pipe=statics[k,:]
            inpipe=''
            for j in range(len(temp_pipe)):
                inpipe=inpipe+' %.6e' % temp_pipe[j]
            commandDS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+\
                    " -A"+str(rad2deg(az[k]))+" -P"
            commandSS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+\
                    " -A"+str(rad2deg(az[k]))+" -P"
            print staname[k]
            print commandSS
            print commandDS
            log=log+staname[k]+'\n'+commandSS+'\n'+commandDS+'\n'
            commandSS=split(commandSS)
            commandDS=split(commandDS)
            ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE,stderr=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
            p=subprocess.Popen(commandSS,stdin=ps.stdout,stdout=open(staname[k]+'.subfault'+num+'.SS.static.rtz','w'),stderr=subprocess.PIPE)     
            out,err=p.communicate()  
            log=log+str(out)+str(err)
            ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE,stderr=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
            p=subprocess.Popen(commandDS,stdin=ps.stdout,stdout=open(staname[k]+'.subfault'+num+'.DS.static.rtz','w'),stderr=subprocess.PIPE)     
            out,err=p.communicate() 
            log=log+str(out)+str(err)       
            #Rotate radial/transverse to East/North
            statics=loadtxt(staname[k]+'.subfault'+num+'.SS.static.rtz')
            u=-statics[2]/100
            r=statics[3]/100
            t=statics[4]/100
            ntemp,etemp=rotate_RT_NE(array([r,r]),array([t,t]),rad2deg(az[k]))
            n=ntemp[0]
            e=etemp[0]
            savetxt(staname[k]+'.subfault'+num+'.SS.static.enu',(e,n,u))
            statics=loadtxt(staname[k]+'.subfault'+num+'.DS.static.rtz')
            u=-statics[2]/100
            r=statics[3]/100
            t=statics[4]/100
            ntemp,etemp=rotate_RT_NE(array([r,r]),array([t,t]),rad2deg(az[k]))
            n=ntemp[0]
            e=etemp[0]
            savetxt(staname[k]+'.subfault'+num+'.DS.static.enu',(e,n,u))
    return log


##########                   Utilities and stuff                      ##########          
        
def cartesian_azimuth(x,y,xs,ys):
    '''
    Compute source to station azimuths when sources given in cartesianc oordinates
    '''
    
    from numpy import arctan,nonzero,pi,intersect1d
    
    az=arctan((x-xs)/(y-ys))
    #Adjsut elements in 2nd quadrant
    ix=nonzero(x-xs<0)[0]
    iy=nonzero(y-ys>0)[0]
    i=intersect1d(ix,iy)
    az[i]=2*pi+az[i]
    #Adjust elements in 3rd and 4th quadrants
    i=nonzero(y-ys<0) 
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