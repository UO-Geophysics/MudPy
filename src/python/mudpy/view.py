# -*- coding: utf-8 -*-
'''
D.Melgar
04/2014

Some routines to make quick and not so quick plots of the forward modeling and 
inversion results
'''

import matplotlib
cdict = {'red': ((0., 1, 1),
                 (0.05, 1, 1),
                 (0.11, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.05, 1, 1),
                   (0.11, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.05, 1, 1),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}

whitejet = matplotlib.colors.LinearSegmentedColormap('whitejet',cdict,256)


def quick_model(rupt):
    '''
    Quick and dirty plot of a .rupt file
    '''
    
    from numpy import genfromtxt,unique,where,zeros
    import matplotlib.pyplot as plt
    
    f=genfromtxt(rupt)
    num=f[:,0]
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    ss=zeros(len(unum))
    ds=zeros(len(unum))
    for k in range(len(unum)):
        i=where(unum[k]==num)
        ss[k]=all_ss[i].sum()
        ds[k]=all_ds[i].sum()
    #Sum them
    slip=(ss**2+ds**2)**0.5
    #Get other parameters
    lon=f[0:len(unum),1]
    lat=f[0:len(unum),2]
    strike=f[0:len(unum),4]
    #Get projection of rake vector
    x,y=slip2geo(ss,ds,strike)
    #Plot
    plt.figure()
    plt.scatter(lon,lat,marker='o',c=slip,s=250,cmap=whitejet)
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    cb=plt.colorbar()
    cb.set_label('Slip (m)')
    plt.quiver(lon,lat,x,y,color='green',width=0.0013)
    plt.grid()
    plt.title(rupt)
    plt.show()
    
def quick_static(gflist,datapath,run_name,run_num,c):
    '''
    Make quick quiver plot of static fields
    
    IN:
        gflist: Tod ecide which stations to plot
        datapath: Where are the data files
    '''
    import matplotlib.pyplot as plt
    from numpy import genfromtxt,where,zeros,meshgrid,linspace
    from matplotlib.mlab import griddata
    from matplotlib.colors import Colormap

    
    GF=genfromtxt(gflist,usecols=3)
    sta=genfromtxt(gflist,usecols=0,dtype='S')
    lon=genfromtxt(gflist,usecols=1,dtype='f')
    lat=genfromtxt(gflist,usecols=2,dtype='f')
    #Read coseismcis
    i=where(GF!=0)[0]
    lon=lon[i]
    lat=lat[i]
    n=zeros(len(i))
    e=zeros(len(i))
    u=zeros(len(i))
    #Get data
    if run_name!='' or run_num!='':
        run_name=run_name+'.'
        run_num=run_num+'.'
    for k in range(len(i)):
        neu=genfromtxt(datapath+run_name+run_num+sta[i[k]]+'.static.neu')
        n[k]=neu[0]#/((neu[0]**2+neu[1]**2)**0.5)
        e[k]=neu[1]#/((neu[0]**2+neu[1]**2)**0.5)
        u[k]=neu[2]#/(2*abs(neu[2]))

            
    #Plot
    plt.figure()
    xi = linspace(min(lon), max(lon), 500)
    yi = linspace(min(lat), max(lat), 500)
    Z = griddata(lon, lat, u, xi, yi)
    X, Y = meshgrid(xi, yi)
    #c=Colormap('bwr')
    #plt.contourf(X,Y,Z,100)
    #plt.colorbar()
    Q=plt.quiver(lon,lat,e,n,width=0.001,color=c)
    plt.scatter(lon,lat,color='b')
    plt.grid()
    plt.title(datapath+run_name+run_num)
    plt.show()
    qscale_en=1
    plt.quiverkey(Q,X=0.1,Y=0.9,U=qscale_en,label=str(qscale_en)+'m')
    
def tile_slip(rupt,nstrike,ndip):
    '''
    Quick and dirty plot of a .rupt file
    '''
    
    from numpy import genfromtxt,unique,where,zeros
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    f=genfromtxt(rupt)
    num=f[:,0]
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    ss=zeros(len(unum))
    ds=zeros(len(unum))
    for k in range(len(unum)):
        i=where(unum[k]==num)
        ss[k]=all_ss[i].sum()
        ds[k]=all_ds[i].sum()
    #Sum them
    slip=(ss**2+ds**2)**0.5
    #Get unit rake vector
    rakess=ss/slip
    rakeds=ds/slip
    #Get indices for plot
    istrike=zeros(nstrike*ndip)
    idip=zeros(nstrike*ndip)
    k=0
    for i in range(ndip):
         for j in range(nstrike):
             istrike[k]=nstrike-j
             idip[k]=ndip-i
             k+=1           
    #Plot
    plt.figure()
    plt.scatter(istrike,idip,marker='o',c=slip,s=250,cmap=whitejet)
    plt.ylabel('Along-dip index')
    plt.xlabel('Along-strike index')
    cb=plt.colorbar()
    cb.set_label('Slip (m)')
    plt.axis('equal')
    plt.xlim(istrike.min()-1,istrike.max()+1)
    plt.ylim(idip.min()-1,idip.max()+1)
    plt.quiver(istrike,idip,rakess,rakeds,color='green',width=0.002)
    plt.grid()
    plt.title(rupt)
    plt.show()

        
def tile_moment(rupt,epicenter,nstrike,ndip):
    '''
    Tile plot of subfault source-time functions
    '''
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from numpy import genfromtxt,unique,zeros,where,meshgrid,linspace
    from mudpy.forward import get_source_time_function,add2stf
    from mudpy.inverse import d2epi
    
    f=genfromtxt(rupt)
    num=f[:,0]
    nfault=nstrike*ndip
    #Get slips
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    #Count number of windows
    nwin=len(where(num==unum[0])[0])
    #Get rigidities
    mu=f[0:len(unum),13]
    #Get rise times
    rise_time=f[0:len(unum),7]
    #Get areas
    area=f[0:len(unum),10]*f[0:len(unum),11]
    #Get coordinates and compute distances
    source=f[0:len(unum),1:4]
    d=d2epi(epicenter,source)
    #Define velocity limits
    vfast=3.6
    vslow=1.0
    
    #Get indices for plot
    istrike=zeros(nstrike*ndip)
    idip=zeros(nstrike*ndip)
    k=0
    for i in range(ndip):
         for j in range(nstrike):
             istrike[k]=nstrike-j-1
             idip[k]=i
             k+=1  
    #Define canvas
    fig, axarr = plt.subplots(ndip, nstrike)
    #Loop over subfaults
    Mmax=0
    for kfault in range(nfault):
        if kfault%10==0:
            print '... working on subfault '+str(kfault)+' of '+str(nfault)
        #Get rupture times for subfault windows
        i=where(num==unum[kfault])[0]
        trup=f[i,12]
        #Get slips on windows
        ss=all_ss[i]
        ds=all_ds[i]
        #Add it up
        slip=(ss**2+ds**2)**0.5
        #Get first source time function
        t1,M1=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip[0])
        #Loop over windows
        for kwin in range(nwin-1):
            #Get next source time function
            t2,M2=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip[kwin+1])
            #Add the soruce time functions
            t1,M1=add2stf(t1,M1,t2,M2)
        #Track maximum moment
        Mmax=max(Mmax,M1.max())
        #Done now plot them
        #get current axis
        ax=axarr[idip[kfault], istrike[kfault]]
        #Make contourf
        Mc=linspace(0,0.98*max(M1),100)
        T,M=meshgrid(t1,Mc)
        i=where(T==0)[0]
        T[i]=0.01
        V=d[kfault]/T
        im=ax.contourf(T,M,V,100,vmin=vslow,vmax=vfast,cmap=cm.spectral)
        #Cover upper part
        ax.fill_between(t1,y1=M1,y2=1.01*M1.max(),color='white')
        #Plot curve
        ax.plot(t1, M1,color='k')
        ax.grid()
        ax.set_xlim([t1[0],t1[-1]])
        ax.xaxis.set_ticks(linspace(t1[0],t1[-1],5))
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
    #Go back and rescale all subplots by maximum moment
    for k in range(ndip):
        for k2 in range(nstrike):
            ax=axarr[k,k2]
            ax.set_ylim([0,Mmax])
    #Fix subplot arrangement
    plt.subplots_adjust(left=0.02, bottom=0.02, right=0.9, top=0.98, wspace=0, hspace=0)
    #Add colorbar
    cbar_ax = fig.add_axes([0.91, 0.15, 0.01, 0.7])
    cb=fig.colorbar(im, cax=cbar_ax)
    cb.set_label('Reference rupture velocity (km/s)')


def source_time_function(rupt,epicenter):
    '''
    Plot source time function of complete ru
    '''
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from numpy import genfromtxt,unique,log10,where,floor
    from mudpy.forward import get_source_time_function,add2stf
    from mudpy.inverse import d2epi
    
    f=genfromtxt(rupt)
    num=f[:,0]
    #Get slips
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    nfault=len(unum)
    #Count number of windows
    nwin=len(where(num==unum[0])[0])
    #Get rigidities
    mu=f[0:len(unum),13]
    #Get rise times
    rise_time=f[0:len(unum),7]
    #Get areas
    area=f[0:len(unum),10]*f[0:len(unum),11]
    #Get coordinates and compute distances
    source=f[0:len(unum),1:4]
    d=d2epi(epicenter,source)
    #Loop over subfaults
    for kfault in range(nfault):
        if kfault%10==0:
            print '... working on subfault '+str(kfault)+' of '+str(nfault)
        #Get rupture times for subfault windows
        i=where(num==unum[kfault])[0]
        trup=f[i,12]
        #Get slips on windows
        ss=all_ss[i]
        ds=all_ds[i]
        #Add it up
        slip=(ss**2+ds**2)**0.5
        if kfault==0:#Get first source time function
            t1,M1=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip[0])
        #Loop over windows
        for kwin in range(nwin-1):
            #Get next source time function
            t2,M2=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip[kwin+1])
            #Add the soruce time functions
            t1,M1=add2stf(t1,M1,t2,M2)
    #Get power
    exp=floor(log10(M1.max()))
    M1=M1/(10**exp)
    plt.figure()
    plt.fill(t1,M1,'b',alpha=0.5)
    plt.plot(t1,M1,color='k')
    plt.grid()
    plt.xlabel('Time(s)')
    plt.ylabel('Moment Rate ('+r'$\times 10^{'+str(int(exp))+r'}$Nm/s)')
    


def Rm(G,lambda_spatial,lambda_temporal,Ls,Lt,bounds,nstrike,ndip,maxR=0.2):
    '''
    Plot model resolution matrix
    '''
    
    from numpy import diag,zeros,arange
    import matplotlib.pyplot as plt
    from numpy.linalg import inv
    
    #Compute model resolution matrix
    Gs=G.transpose().dot(G)+(lambda_spatial**2)*Ls.transpose().dot(Ls)+(lambda_temporal**2)*Lt.transpose().dot(Lt)
    R=inv(Gs).dot(G.transpose()).dot(G)
    #Keep diagonals only
    r=diag(R)
    ids=arange(1,len(r),2)
    r=r[ids]
    #Get indices for plot
    istrike=zeros(nstrike*ndip)
    idip=zeros(nstrike*ndip)
    k=0
    for i in range(ndip):
         for j in range(nstrike):
             istrike[k]=nstrike-j
             idip[k]=ndip-i
             k+=1           
    #Plot
    plt.figure()
    plt.scatter(istrike,idip,marker='o',c=r,s=250,cmap=plt.cm.PuBuGn,vmax=maxR)
    plt.ylabel('Along-dip index')
    plt.xlabel('Along-strike index')
    cb=plt.colorbar()
    cb.set_label('Diagonal Value of R')
    plt.axis('equal')
    plt.xlim(istrike.min()-1,istrike.max()+1)
    plt.ylim(idip.min()-1,idip.max()+1)
    plt.grid()
    plt.title('Model Resolution')
    plt.show()
    return R

def tslice(rupt,out,dt,cumul):
    '''
    Quick and dirty plot of a .rupt file
    '''
    
    from numpy import genfromtxt,unique,where,zeros,arange,intersect1d,trapz
    import matplotlib.pyplot as plt
    from string import rjust
    
    delta_t=0.01
    f=genfromtxt(rupt)
    trupt=f[:,12]
    trise=f[:,7]
    all_ss=f[:,8]
    all_ds=f[:,9]
    num=f[:,0]
    #Get other parameters
    #lon=f[0:len(unum),1]
    #lat=f[0:len(unum),2]
    #strike=f[0:len(unum),4]
    #Decide on time vector
    tslice=arange(0,trupt.max()+dt,dt)
    #Determine triangle height at all subfaults
    hss=2*all_ss/trise
    hds=2*all_ds/trise
    #Cumulative
    ss_cumul=zeros(len(f))
    ds_cumul=zeros(len(f))
    #Determine time series fo each triangle
    t=arange(0,trupt.max()+trise[0],delta_t)
    for kslice in range(len(tslice-2)):
        print str(kslice)+'/'+str(len(tslice)-1)
        #And initalize slice vectors
        ss_slice=zeros(len(f))
        ds_slice=zeros(len(f))
        for kfault in range(len(f)):
            yss=zeros(t.shape)
            yds=zeros(t.shape)
            #Up going
            i1=where(t>=trupt[kfault])[0]
            i2=where(t<=(trupt[kfault]+trise[0]/2))[0] #Ascending triangle
            i=intersect1d(i1,i2)
            yss[i]=(2*hss[kfault]/trise[0])*t[i]-(2*hss[kfault]*trupt[kfault]/trise[0])
            yds[i]=(2*hds[kfault]/trise[0])*t[i]-(2*hds[kfault]*trupt[kfault]/trise[0])
            #Down going
            i1=where(t>(trupt[kfault]+trise[0]/2))[0]
            i2=where(t<=(trupt[kfault]+trise[0]))[0] #Ascending triangle
            i=intersect1d(i1,i2)
            yss[i]=(-2*hss[kfault]/trise[0])*t[i]+(2*hss[kfault]/trise[0])*(trupt[kfault]+trise[0])
            yds[i]=(-2*hds[kfault]/trise[0])*t[i]+(2*hds[kfault]/trise[0])*(trupt[kfault]+trise[0])
            #Now integrate slip at pertinent time interval
            i1=where(t>=tslice[kslice])[0]
            i2=where(t<=tslice[kslice+1])[0]
            i=intersect1d(i1,i2)
            ss_slice[kfault]=trapz(yss[i],t[i])
            ds_slice[kfault]=trapz(yds[i],t[i])
        #Combine into single model for that time slice
        ss_cumul=ss_cumul+ss_slice
        ds_cumul=ds_cumul+ds_slice
        unum=unique(num)
        lon=f[0:len(unum),1]
        lat=f[0:len(unum),2]
        strike=f[0:len(unum),4]
        ss=zeros(len(unum))
        ds=zeros(len(unum))
        for k in range(len(unum)):
            if cumul==0:
                i=where(unum[k]==num)
                ss[k]=ss_slice[i].sum()
                ds[k]=ds_slice[i].sum()    
            else:
                i=where(unum[k]==num)
                ss[k]=ss_cumul[i].sum()
                ds[k]=ds_cumul[i].sum()        
        slip=(ss**2+ds**2)**0.5
        #Plot
        #Get projection of rake vector
        x,y=slip2geo(ss,ds,strike)
        #Plot
        plt.figure()
        plt.scatter(lon,lat,marker='o',c=slip,s=250,cmap=plt.cm.gnuplot2_r,vmin=0,vmax=35)
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        cb=plt.colorbar()
        plt.quiver(lon,lat,x,y,color='green',width=0.0013)
        plt.grid()
        if cumul==0:
            cb.set_label('Slip (m)')
            plt.title('t = '+str(tslice[kslice])+'s to '+str(tslice[kslice+1])+'s') 
            plt.savefig(out+rjust(str(kslice),4,'0')+'.kin_slice.png')
        else:
            cb.set_label('Cumulative Slip (m)')
            plt.title('t = '+str(tslice[kslice+1])+'s')
            plt.savefig(out+rjust(str(kslice),4,'0')+'.kin_cumulative.png')
        plt.close("all")
    
    
def synthetics(home,project_name,run_name,run_number,gflist,vord,decimate,lowpass,t_lim,sort):
    '''
    Plot synthetics vs real data
    
    gflist: The GF control fiel that decides what to plot/not plot
    datapath
    '''
    from obspy import read
    from numpy import genfromtxt,where,argsort
    import matplotlib.pyplot as plt
    import matplotlib
    from mudpy.green import stdecimate 
    from mudpy.forward import lowpass as lfilter
    
    matplotlib.rcParams.update({'font.size': 14})
    #Decide what to plot
    sta=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=0,dtype='S')
    lon=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[1],dtype='f')
    lat=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[2],dtype='f')
    gf=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[4,5],dtype='f')
    datapath=home+project_name+'/data/waveforms/'
    synthpath=home+project_name+'/output/inverse_models/waveforms/'
    if vord.lower()=='d':
        kgf=0 #disp
        datasuffix='kdisp'
        synthsuffix='disp'
    elif vord.lower()=='v':
        kgf=1 #disp
        datasuffix='kvel'
        synthsuffix='vel'
    #Decide on sorting
    i=where(gf[:,kgf]==1)[0]  
    if sort.lower()=='lon':
        j=argsort(lon[i])[::-1]
        i=i[j]
    elif sort.lower()=='lat':
        j=argsort(lat[i])[::-1] 
        i=i[j]
    nsta=len(i)
    fig, axarr = plt.subplots(nsta, 3)  
    for k in range(len(i)):
        n=read(datapath+sta[i[k]]+'.'+datasuffix+'.n')
        e=read(datapath+sta[i[k]]+'.'+datasuffix+'.e')
        u=read(datapath+sta[i[k]]+'.'+datasuffix+'.u')
        if lowpass!=None:
            fsample=1./e[0].stats.delta
            e[0].data=lfilter(e[0].data,lowpass,fsample,10)
            n[0].data=lfilter(n[0].data,lowpass,fsample,10)
            u[0].data=lfilter(u[0].data,lowpass,fsample,10)
        if decimate!=None:
            n=stdecimate(n,decimate)
            e=stdecimate(e,decimate)
            u=stdecimate(u,decimate)
        ns=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.n.sac')
        es=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.e.sac')
        us=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.u.sac')
        #Make plot
        axn=axarr[k,0]
        axe=axarr[k,1]
        axu=axarr[k,2]
        axn.plot(n[0].times(),n[0].data,'k',ns[0].times(),ns[0].data,'r')
        axn.grid(which='both')
        axe.plot(e[0].times(),e[0].data,'k',es[0].times(),es[0].data,'r')
        axe.grid(which='both')
        axu.plot(u[0].times(),u[0].data,'k',us[0].times(),us[0].data,'r')
        axu.grid(which='both')
        axe.yaxis.set_ticklabels([])
        axu.yaxis.set_ticklabels([])
        axe.set_xlim(t_lim)
        axn.set_xlim(t_lim)
        axu.set_xlim(t_lim)
        axn.yaxis.set_ticklabels([])
        axe.yaxis.set_ticklabels([])
        axu.yaxis.set_ticklabels([])
        axn.yaxis.grid(False)
        axe.yaxis.grid(False)
        axu.yaxis.grid(False)
        axn.yaxis.set_ticks([])
        axe.yaxis.set_ticks([])
        axu.yaxis.set_ticks([])
        
        #Annotations
        trange=t_lim[1]-t_lim[0]
        sign=1.
        if abs(min(n[0].data))>max(n[0].data):
            sign=-1. 
        nmax='%.3f' % (sign*max(abs(n[0].data)))
        sign=1.
        if abs(min(ns[0].data))>max(ns[0].data):
            sign=-1. 
        nsmax='%.3f' % (sign*max(abs(ns[0].data)))
        sign=1.
        nlims=axn.get_ylim()
        nrange=nlims[1]-nlims[0]
        
        if abs(min(e[0].data))>max(e[0].data):
            sign=-1.         
        emax='%.3f' % (sign*max(abs(e[0].data)))
        sign=1.
        if abs(min(es[0].data))>max(es[0].data):
            sign=-1. 
        esmax='%.3f' % (sign*max(abs(es[0].data)))
        sign=1.
        elims=axe.get_ylim()
        erange=elims[1]-elims[0]
        
        if abs(min(u[0].data))>max(u[0].data):
            sign=-1. 
        umax='%.3f' % (sign*max(abs(u[0].data)))
        sign=1.
        if abs(min(us[0].data))>max(us[0].data):
            sign=-1 
        usmax='%.3f' % (sign*max(abs(us[0].data)))
        sign=1.
        ulims=axu.get_ylim()
        urange=ulims[1]-ulims[0]
        
        axn.annotate(nmax,xy=(t_lim[0]+0.02*trange,nlims[0]+0.02*nrange),fontsize=12)
        axe.annotate(emax,xy=(t_lim[0]+0.02*trange,elims[0]+0.02*erange),fontsize=12)
        axu.annotate(umax,xy=(t_lim[0]+0.02*trange,ulims[0]+0.02*urange),fontsize=12)
        axn.annotate(nsmax,xy=(t_lim[0]+0.02*trange,nlims[0]+0.7*nrange),fontsize=12,color='red')
        axe.annotate(esmax,xy=(t_lim[0]+0.02*trange,elims[0]+0.7*erange),fontsize=12,color='red')
        axu.annotate(usmax,xy=(t_lim[0]+0.02*trange,ulims[0]+0.7*urange),fontsize=12,color='red')
        #Station name
        axn.set_ylabel(sta[i[k]],rotation=0)
        if k==0:
            axn.set_title('North (m)')
            axe.set_title('East (m)')
            axu.set_title('Up (m)')
        if k!=len(i)-1:
            axn.xaxis.set_ticklabels([])
            axe.xaxis.set_ticklabels([])
            axu.xaxis.set_ticklabels([])
            xtick=axn.xaxis.get_majorticklocs()
            xtick=xtick[1:]
        if k==len(i)-1: #Last plot
            axe.set_xlabel('Time (s)')
            axn.xaxis.set_ticks(xtick)
            axe.xaxis.set_ticks(xtick)
            axu.xaxis.set_ticks(xtick)
    plt.subplots_adjust(left=0.2, bottom=0.05, right=0.8, top=0.95, wspace=0, hspace=0)
                

def ABIC(home,project_name,run_name):
    '''
    plot values of ABIC vs smoothing parameter for model selection
    '''
    from glob import glob
    from numpy import zeros,argmin
    import matplotlib.pyplot as pl
    
    #Text rendering
    pl.rc('font',family='serif')
    #Get list of log files
    outdir=home+project_name+'/output/inverse_models/models/'
    plotdir=home+project_name+'/plots/'
    logs=glob(outdir+'*'+run_name+'.????.log')
    ABIC=zeros(len(logs))
    ls=zeros(len(logs))
    print 'Gathering statistics for '+str(len(logs))+' inversions...'
    for k in range(len(logs)):
        with open(logs[k]) as f:
            for line in f:
                if 'ABIC' in line:
                    ABIC[k]=float(line.split('=')[1])
                if 'lambda_spatial' in line:
                    ls[k]=float(line.split('=')[1])
    #Get the minimum
    imin=argmin(ABIC)
    #Plot the thing
    pl.figure()
    pl.semilogx(ls,ABIC,'k',linewidth=2)
    pl.grid(which='both')
    pl.semilogx(ls[imin],ABIC[imin],marker='*',color='r',markersize=14)
    pl.xlabel(r'$\lambda$',fontsize=14)
    pl.ylabel('ABIC',fontsize=14)
    pl.annotate(r'$\lambda$'+'='+str(ls[imin]),xy=(ls[imin],ABIC[imin]),xytext=(ls[imin],ABIC[imin]-0.05*(max(ABIC)-ABIC[imin])))
    pl.title('Run Name: '+run_name)
    pl.savefig(plotdir+run_name+'.ABIC.png')
    print 'ABIC is minimized at inversion '+logs[imin]
    print '... lambda = '+repr(ls[imin])
    pl.show()
    
def ABIC2D(home,project_name,run_name,(ABICmin,ABICmax)):
    '''
    plot 2D values of ABIC vs smoothing parameter for model selection
    '''
    from glob import glob
    from numpy import zeros,argmin,log10,linspace
    import matplotlib.pyplot as plt
    from matplotlib import mlab as ml
    
    #Text rendering
    plt.rc('font',family='serif')
    #Get list of log files
    outdir=home+project_name+'/output/inverse_models/models/'
    for k in range(len(run_name)):
        if k==0:
            logs=glob(outdir+'*'+run_name[k]+'.????.log')
        else:
            logs+=glob(outdir+'*'+run_name[k]+'.????.log') 
    ABIC=zeros(len(logs))
    ls=zeros(len(logs))
    lt=zeros(len(logs))
    print 'Gathering statistics for '+str(len(logs))+' inversions...'
    for k in range(len(logs)):
        with open(logs[k]) as f:
            for line in f:
                if 'ABIC' in line:
                    ABIC[k]=float(line.split('=')[1])
                if 'lambda_spatial' in line:
                    ls[k]=log10(float(line.split('=')[1]))
                if 'lambda_temporal' in line:
                    lt[k]=log10(float(line.split('=')[1]))
    #Get the minimum
    imin=argmin(ABIC)
    #Grid
    ABIC=ABIC/1000
    lsi=linspace(ls.min(),ls.max(),100)
    lti=linspace(ls.min(),lt.max(),100)
    ABICi=ml.griddata(ls,lt,ABIC,lsi,lti)
    #Plot the thing
    plt.figure()
    plt.pcolormesh(lsi,lti,ABICi,cmap=plt.cm.spectral_r,vmin=ABICmin,vmax=ABICmax)
    cb=plt.colorbar()
    plt.scatter(ls,lt,c='w',marker='o')
    plt.ylabel(r'$\log(\lambda_t$)',fontsize=18)
    plt.xlabel(r'$\log(\lambda_s$)',fontsize=18)
    cb.set_label(r'ABIC $(\times10^3)$')
    plt.scatter(ls[imin],lt[imin],marker='*',color='r',s=100)
    plt.title(r'$\lambda_s^{min} = $%.4e , $\lambda_t^{min} = $%.4e' % (10**ls[imin],10**lt[imin]))
    plt.show()
    print 'ABIC is minimized at inversion '+logs[imin]
    print '... ls = '+repr(10**ls[imin])+' , lt = '+repr(10**lt[imin])
        
    
def coherence(home,project_name,run_name,run_number,GF_list,vord,f_lims):
    '''
    Plot coherences
    '''
    import matplotlib.pyplot as plt
    from numpy import load,where,genfromtxt,array,log10
    
    force_lims=array([1./200,0.5])
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    gf=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[4,5],dtype='f')
    datapath=home+project_name+'/analysis/frequency/'
    if vord.lower()=='d':
        kgf=0 #disp
        suffix='disp'
    elif vord.lower()=='v':
        kgf=1 #disp
        suffix='vel'
    i=where(gf[:,kgf]==1)[0]  
    #Initalize canvas
    fig, axarr = plt.subplots(len(i), 3)  
    for k in range(len(i)):
        #Read coherences
        coh=load(datapath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+suffix+'.coh.npz')
        fn=coh['fn']
        fe=coh['fe']
        fu=coh['fu']
        cn=coh['cn']
        ce=coh['ce']
        cu=coh['cu']
        #Let's plot them
        #get current axis
        axn=axarr[k,0]
        axe=axarr[k,1]
        axu=axarr[k,2]
        #Plot
        axn.semilogx(fn,cn)
        axe.semilogx(fe,ce,'g')
        axu.semilogx(fu,cu,'r')
        axn.grid(which='both')
        axe.grid(which='both')
        axu.grid(which='both')
        #Arrange axes
        axn.set_xlim(f_lims)
        axe.set_xlim(f_lims)
        axu.set_xlim(f_lims)
        axn.set_ylim([0,1])
        axe.set_ylim([0,1])
        axu.set_ylim([0,1])
        #Text labels
        axn.yaxis.set_ticks(array([0.25,0.5,0.75,1.0]))
        axe.yaxis.set_ticks(array([0.25,0.5,0.75,1.0]))
        axu.yaxis.set_ticks(array([0.25,0.5,0.75,1.0]))
        axe.yaxis.set_ticklabels([])
        axu.yaxis.set_ticklabels([])
        axn.yaxis.set_ticklabels(['','0.5','','1.0'])
        if k!=len(i)-1:
            axn.xaxis.set_ticklabels([])
            axe.xaxis.set_ticklabels([])
            axu.xaxis.set_ticklabels([])
        if k==0: #First plot add some labels
            axn.set_title('North')
            axe.set_title('East')
            axu.set_title('Up')
        if k==len(i)-1: #Last plot
            axe.set_xlabel('Frequency (Hz)')
        #Annotate with station name
        xyannot=(axn.get_xlim()[0]+0.01*log10((log10(axn.get_xlim()[1])-log10(axn.get_xlim()[0]))),axn.get_ylim()[0]+0.05)
        axn.annotate(sta[i[k]], xy=xyannot)
    plt.subplots_adjust(left=0.25, bottom=0.05, right=0.75, top=0.95, wspace=0, hspace=0)
            
def average_coherence(home,project_name,run_name,run_number,GF_list,vord,num_components):
    '''
    Plot coherences
    '''
    import matplotlib.pyplot as plt
    from numpy import load,where,genfromtxt,array,zeros,interp
    
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    gf=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[4,5],dtype='f')
    datapath=home+project_name+'/analysis/frequency/'
    if vord.lower()=='d':
        kgf=0 #disp
        suffix='disp'
    elif vord.lower()=='v':
        kgf=1 #disp
        suffix='vel'
    i=where(gf[:,kgf]==1)[0]  
    for k in range(len(i)):
        #Read coherences
        coh=load(datapath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+suffix+'.coh.npz')
        fn=coh['fn']
        fe=coh['fe']
        fu=coh['fu']
        cn=coh['cn']
        ce=coh['ce']
        cu=coh['cu'] 
        if k==0:
            f=fn
            c1=zeros(cn.shape)
            c1e=zeros(cn.shape)
            c1n=zeros(cn.shape) 
            c1u=zeros(cn.shape)    
        if num_components==1: #Average all
            try:
                c1+=cn
                c1+=ce
                c1+=cu
            except: #Coherence is the wrong size
                cn=interp(f,fn,cn)
                ce=interp(f,fn,ce)
                cu=interp(f,fn,cu)
                c1+=cn
                c1+=ce
                c1+=cu
        else:
            try:
                c1e+=ce
                c1n+=cn
                c1u+=cu
            except: #Coherence is the wrong size
                cn=interp(f,fn,cn)
                ce=interp(f,fn,ce)
                cu=interp(f,fn,cu)
                c1n+=cn
                c1e+=ce
                c1u+=cu
    #Normalize
    c1=c1/(3*len(i))
    c1e=c1e/len(i)
    c1n=c1n/len(i)
    c1u=c1u/len(i)
    #Plot
    plt.figure()
    if num_components==1:
        plt.semilogx(fn,c1)
        plt.fill_between(fn,y1=0,y2=c1,color='r',alpha=0.5)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Mean Coherence')
        plt.grid(which='both')
        plt.xlim(f.min(),f.max())
    else:
        plt.semilogx(fn,c1n,fe,c1e,fu,c1u)
        plt.legend(['North','East','Up'])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Mean Coherence')
        plt.grid(which='both')
        plt.xlim(f.min(),f.max())

#########                  Supporting tools                       ##############

def slip2geo(ss,ds,strike):
    '''
    Determine geogrpahical orientation of rake vector
    '''
    from numpy import deg2rad,sin,cos
    
    #Normalize slips
    ds=ds/((ds**2+ss**2)**0.5)
    ss=ss/((ds**2+ss**2)**0.5)
    #determine contribution of ds and ss slips
    xds=ds*sin(deg2rad(strike-90))
    yds=ds*cos(deg2rad(strike-90))
    xss=ss*sin(deg2rad(strike))
    yss=ss*cos(deg2rad(strike))
    #Add em up
    x=xss+xds
    y=yss+yds
    return x,y