'''
InSAR tools gop in here
Diego Melgar
UC Berkeley
05/2015
'''

def quadtree2mudpy(home,project_name,quadtree_file,out_file,prefix):
    '''
    Convert quadtree Matlab generated file into a .sta file with station codes 
    and locations. Also generate individaul .neu files for each los displacement
    '''
    
    from numpy import genfromtxt,c_,savetxt
    from string import rjust
    from matplotlib import pyplot as plt
    from matplotlib import cm
    
    insar=genfromtxt(quadtree_file)
    plt.figure()
    #Write to file
    f=open(out_file,'w')
    for k in range(len(insar)):
        print k
        sta=prefix+rjust(str(k),4,'0')
        out='%s\t%.6f\t%.6f\n' %(sta,insar[k,0],insar[k,1])
        f.write(out)
        #Now make .los file
        los_c=insar[k,2]
        N=insar[k,4]
        E=insar[k,3]
        U=insar[k,5]
        los=c_[los_c,N,E,U]
        #Make plot as you go to verify
        plt.scatter(insar[k,0],insar[k,1],c=los_c,cmap=cm.seismic,s=90,vmin=-1.2,vmax=1.2)
        plt.legend(['LOS'])
        savetxt(home+project_name+'/data/statics/'+sta+'.los',los,header='LOS(m),los unit vector (positive towards satellite) n,e,u')
    plt.colorbar()
    f.close()    
    plt.show()  
    
    
    
def un_nanify(infile,outfile):
    from numpy import genfromtxt,savetxt,where,nan
    insar=genfromtxt(infile)
    i=where(infile[:,2]!=nan)[0] 
    insar=insar[i,:]
    savetxt(outfile,insar,fmt='%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f')