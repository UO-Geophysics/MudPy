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
    
    from numpy import genfromtxt
    import matplotlib.pyplot as plt
    
    f=genfromtxt(rupt)
    lon=f[:,1]
    lat=f[:,2]
    strike=f[:,4]
    ss=f[:,8]
    ds=f[:,9]
    slip=(ss**2+ds**2)**0.5
    #Get projection of rake vector
    x,y=slip2geo(ss,ds,strike)
    #Plot
    plt.scatter(lon,lat,marker='o',c=slip,s=250,cmap=plt.cm.gnuplot_r)
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    cb=plt.colorbar()
    cb.set_label('Slip (m)')
    plt.quiver(lon,lat,x,y,color='green',width=0.002)
    plt.grid()
    plt.title(rupt)
    plt.show()
    
def quick_static_plot(gflist,datapath,h_or_u,c):
    '''
    Make quick quiver plot of static fields
    
    IN:
        gflist: Tod ecide which stations to plot
        datapath: Where are the data files
    '''
    import matplotlib.pyplot as plt
    from numpy import genfromtxt,where,zeros
    
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
    for k in range(len(i)):
        neu=genfromtxt(datapath+sta[i[k]]+'.static.neu')
        n[k]=neu[0]
        e[k]=neu[1]
        u[k]=neu[2]
    #Plot
    if h_or_u.lower()=='u': #plot verticals
        e=zeros(e.shape)
        n=u
    Q=plt.quiver(lon,lat,e,n,width=0.0035,color=c)
    qscale_en=0.1*max((n**2+e**2)**0.5)
    plt.quiverkey(Q,X=0.1,Y=0.9,U=qscale_en,label=str(qscale_en)+'m')
    
#########                  Supporting tools                       ##############

def slip2geo(ss,ds,strike):
    '''
    Determine geogrpahical orientation of rake vector
    '''
    from numpy import deg2rad,sin,cos
    
    #determine contribution of ds and ss slips
    xds=ds*sin(deg2rad(strike-90))
    yds=ds*cos(deg2rad(strike-90))
    xss=ss*sin(deg2rad(strike))
    yss=ss*cos(deg2rad(strike))
    #Add em up
    x=xss+xds
    y=yss+yds
    return x,y