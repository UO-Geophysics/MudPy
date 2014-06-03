def make_dtopo():
    '''
    Make geoclaw dtopo file
    '''
    from numpy import genfromtxt,zeros
    
    #Run params
    f='/Users/dmelgarm/Research/Slip_Inv/tohoku_tsunami/'
    stafile='tohoku.sta'
    dlon=0.033333
    dlat=0.033333
    dt=5
    stat_or_dyn='s'
    
    #Get station list
    
    sta=genfromtxt(f+'data/station_info/'+stafile,usecols=0,dtype='S4')
    s=genfromtxt(f+'data/station_info/'+stafile,usecols=[1,2])
    lon=s[:,0]
    lat=s[:,1]
    if stat_or_dyn.lower()=='s':
        n=zeros(len(sta))
        e=n.copy()
        u=n.copy()
        for ksta in range(len(sta)):
            print ksta
            neu=genfromtxt(f+'output/forward_models/'+str(sta[ksta])+'.static.neu')
            n[ksta]=neu[0]
            e[ksta]=neu[1]
            u[ksta]=neu[2]
            print neu[2]

def make_grid(lon1,lon2,lat1,lat2,dlon,dlat,outfile):
    '''
    Make a grid of points where you will request coseismic deformation
    '''
    
    from numpy import arange
    from string import rjust
    
    lon=arange(lon1,lon2+dlon,dlon)
    lat=arange(lat1,lat2+dlat,dlat)
    lat=lat[::-1]
    #Turn into pseudo station file
    sta=arange(0,len(lon)*len(lat),1)
    ksta=0
    fid=open(outfile,'w')
    for klat in range(len(lat)):
        for klon in range(len(lon)):
            outsta=rjust(str(sta[ksta]),4,'0')
            outlon='%.4f' % lon[klon]
            outlat='%.4f' % lat[klat]
            ksta+=1
            out=outsta+'\t'+outlon+'\t'+outlat+'\n'
            fid.write(out)
    fid.close()    