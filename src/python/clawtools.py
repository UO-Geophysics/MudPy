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
        