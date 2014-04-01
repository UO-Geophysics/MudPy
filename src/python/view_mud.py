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
    slip=(f[:,8]**2+f[:,9]**2)**0.5
    #Plot
    plt.figure()
    plt.scatter(lon,lat,marker='o',c=slip,s=250,cmap=plt.cm.gnuplot_r)
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    cb=plt.colorbar()
    cb.set_label('Slip (m)')
    plt.grid()
    plt.title(rupt)
    plt.show()