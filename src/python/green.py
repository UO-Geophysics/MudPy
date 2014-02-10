'''
Diego Melgar 01/2014
Green functions routines for kinematic source models
'''

deltad=0.1

def run_green(source_depths,station_distances,model_name,dt,NFFT):
    '''
    Compute GFs using Zhu & Rivera code for a given velocity model, source depth
    and receiver distances
    '''
    import subprocess
    from shlex import split
    
    print("--> Computing GFs for source depth "+str(source_depths)+"km")
    diststr=''
    for k in range(len(station_distances)):
        diststr=diststr+' %.3f' % station_distances[k]
    command=split("fk.pl -M"+model_name+"/"+str(source_depths)+" -N"+str(NFFT)+"/"+str(dt)+diststr)
    #print "fk.pl -M"+model_name+"/"+str(source_depth)+" -N"+str(NFFT)+"/"+str(dt)+diststr
    p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err=p.communicate() 
    print out
    print err
    
    
def run_syn(source,station_file,delta_distance,green_dir,integrate):
    '''
    Make synthetics at stations
    '''

    import os
    import subprocess
    from numpy import genfromtxt,round,rad2deg
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
    green_dir=green_dir+"_"+str(int(zs))+"/"
    print("--> Computing synthetics at stations for the source at ("+str(xs)+" , "+str(ys)+")")
    staname=genfromtxt(station_file,dtype="S6",usecols=0)
    x=genfromtxt(station_file,dtype="f8",usecols=1)
    y=genfromtxt(station_file,dtype="f8",usecols=2)
    #Compute distances and azimuths
    d=((x-xs)**2+(y-ys)**2)**0.5
    az=cartesian_azimuth(x,y,xs,ys)
    #round dsitances to GF precision
    d=round(d/delta_distance)*delta_distance
    #Move to output folder
    os.chdir(green_dir)
    for k in range(len(d)):
        diststr='%.3f' % d[k] #Need current distance in string form for external call
        if integrate==1: #Make displ.
            command=split("syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rake)+" -D"+str(duration)+ \
                "/"+str(rise)+" -A"+str(rad2deg(az[k]))+" -O"+staname[k]+".disp.x -G"+green_dir+diststr+".grn.0")
        else: #Make vel.
            command=split("syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rake)+" -D"+str(duration)+ \
                "/"+str(rise)+" -A"+str(rad2deg(az[k]))+" -O"+staname[k]+".vel.x -G"+green_dir+diststr+".grn.0")
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err=p.communicate() 
        #Rotate to NE
        if integrate==1: #'tis displacememnt
            r=read(staname[k]+'.disp.r')
            t=read(staname[k]+'.disp.t')
            ntemp,etemp=rotate_RT_NE(r[0].data,t[0].data,rad2deg(az[k]))
            n=r.copy()
            n[0].data=ntemp
            e=t.copy()
            e[0].data=etemp
            n.write(staname[k]+'.disp.n',format='SAC')
            e.write(staname[k]+'.disp.e',format='SAC')
            #Correct polarity in vertical, ZR code has down=positive, we don't like that precious do we?
            down=read(staname[k]+'.disp.z')
            up=down.copy()
            up[0].data=down[0].data*-1
            up.write(staname[k]+'.disp.z',format='SAC')
        else: #'tis velocity sire
            r=read(staname[k]+'.vel.r')
            t=read(staname[k]+'.vel.t')
            ntemp,etemp=rotate_RT_NE(r[0].data,t[0].data,rad2deg(az[k]))
            n=r.copy()
            n[0].data=ntemp
            e=t.copy()
            e[0].data=etemp
            n.write(staname[k]+'.vel.n',format='SAC')
            e.write(staname[k]+'.vel.e',format='SAC')
            #Correct polarity in vertical, ZR code has down=positive, we don't like that precious do we?
            down=read(staname[k]+'.vel.z')
            up=down.copy()
            up[0].data=down[0].data*-1
            up.write(staname[k]+'.vel.z',format='SAC')


##########                   Utilities and stuff                      ##########          
        
def cartesian_azimuth(x,y,xs,ys):
    from numpy import arctan,nonzero,pi
    az=arctan((x-xs)/(y-ys))
    i=nonzero(y<0)
    az[i]=pi+az[i]
    return az