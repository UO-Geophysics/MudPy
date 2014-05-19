# -*- coding: utf-8 -*-
'''
D.Melgar
04/2014

Some routines to make quick and not so quick plots of the forward modeling and 
inversion results
'''

def quick_model_plot(rupt):
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
    plt.scatter(lon,lat,marker='o',c=slip,s=250,cmap=plt.cm.gnuplot_r)
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    cb=plt.colorbar()
    cb.set_label('Slip (m)')
    plt.quiver(lon,lat,x,y,color='green',width=0.0013)
    plt.grid()
    plt.title(rupt)
    plt.show()
    
def quick_static_plot(gflist,datapath,run_name,run_num,c):
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
    
def tile_plot(rupt,nstrike,ndip):
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
    plt.scatter(istrike,idip,marker='o',c=slip,s=250,cmap=plt.cm.gnuplot_r)
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

def model_tslice(rupt,out,dt,cumul):
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
    
    
def plot_synthetics(gflist,datapath,datasuffix,synthpath,synthsuffix,vord):
    '''
    Plot synthetics vs real data
    '''
    from obspy import read
    from numpy import genfromtxt,where
    import matplotlib.pyplot as plt
    import matplotlib
    
    matplotlib.rcParams.update({'font.size': 16})
    #Decide what to plot
    sta=genfromtxt(gflist,usecols=0,dtype='S')
    gf=genfromtxt(gflist,usecols=[4,5],dtype='f')
    kgf=0 #disp
    i=where(gf[:,kgf]==1)[0]
    if gf[i,kgf].sum()>0:
        #Initalize the plot canvas
        plt.figure()
        nsta=len(i)
        left=0.05
        width=0.28
        bottom=0.05
        height=0.75/nsta
        Dy=height+0.03
        for k in range(len(i)):
            n=read(datapath+sta[i[k]]+'.'+datasuffix+'.n')
            e=read(datapath+sta[i[k]]+'.'+datasuffix+'.e')
            u=read(datapath+sta[i[k]]+'.'+datasuffix+'.u')
            ns=read(synthpath+synthsuffix+'.'+sta[i[k]]+'.disp.n.sac')
            es=read(synthpath+synthsuffix+'.'+sta[i[k]]+'.disp.e.sac')
            us=read(synthpath+synthsuffix+'.'+sta[i[k]]+'.disp.u.sac')
            #Make plot
            dy=Dy*k
            rect=[left,bottom+dy,width,height]
            axn=plt.axes(rect)
            axn.plot(n[0].times(),n[0].data,'k',ns[0].times(),ns[0].data,'r')
            axn.grid(which='both')
            axn.set_ylabel(sta[i[k]])
            rect=[left+width+0.03,bottom+dy,width,height]
            axe=plt.axes(rect)
            axe.plot(e[0].times(),e[0].data,'k',es[0].times(),es[0].data,'r')
            axe.grid(which='both')
            rect=[left+2*width+0.06,bottom+dy,width,height]
            axz=plt.axes(rect)
            axz.plot(u[0].times(),u[0].data,'k',us[0].times(),us[0].data,'r')
            axz.grid(which='both')
            if k==0:
                axn.set_xlabel('Time (s)')
                axe.set_xlabel('Time (s)')
                axz.set_xlabel('Time (s)')
            if k!=0:
                axn.get_xaxis().set_ticklabels([])
                axe.get_xaxis().set_ticklabels([])
                axz.get_xaxis().set_ticklabels([])
            if k==nsta-1:
                axn.set_title('North (m)')
                axe.set_title('East (m)')
                axz.set_title('Up (m)')
                axn.legend(['Observed','Inversion'])
                

def plotABIC(home,run_name):
    '''
    plot values of ABIC vs smoothing parameter for model selection
    '''
    from glob import glob
    from numpy import zeros,argmin
    import matplotlib.pyplot as pl
    
    #Text rendering
    pl.rc('font',family='serif')
    #Get list of log files
    outdir=home+'output/inverse_models/models/'
    plotdir=home+'/plots/'
    logs=glob(outdir+'*'+run_name+'*.log')
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