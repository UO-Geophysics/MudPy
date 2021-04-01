'''
InSAR tools go in here
Diego Melgar
University of Oregon
'''

def quadtree2mudpy(home,project_name,quadtree_file,gflist_file,prefix):
    '''
    Convert quadtree Matlab generated file into a .sta file with station codes 
    and locations. Also generate individaul .neu files for each los displacement
    '''
    
    from numpy import genfromtxt,c_,savetxt
    from matplotlib import pyplot as plt
    from matplotlib import cm
    
    insar=genfromtxt(quadtree_file)
    plt.figure()
    vmin=insar[:,2].min()
    vmax=insar[:,2].max()
    #Write to file
    gflist=open(gflist_file,'w')
    #Write gflist header
    gflist.write('#station	lat	lon	static	disp	vel	tsun	strain	Static file	 displacement file	 velocity file	tsunami file	strain file	static sigmas(n	e	u)	displacement sigmas(n	e	u)	velocity sigmas(n	e	u)	tsunami sigma	strain sigmas(5 components?)\n')																						
    for k in range(len(insar)):
        sta=prefix+str(k).rjust(4,'0')
        #Now make .los file
        los_c=insar[k,2]
        N=insar[k,4]
        E=insar[k,3]
        U=insar[k,5]
        los=c_[los_c,N,E,U]
        savetxt(home+project_name+'/data/statics/'+sta+'.los',los,header='LOS(m),los unit vector (positive towards satellite) n,e,u')
        #Generate gflist file as well
        gflist.write('%s\t%10.6f\t%10.6f\t0\t0\t0\t0\t1\t/foo/bar\t/foo/bar\t/foo/bar\t/foo/bar\t%s\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\n'%(sta,insar[k,0],insar[k,1],home+project_name+'/data/statics/'+sta+'.los'))
    
    #Make plot  to verify
    plt.scatter(insar[:,0],insar[:,1],c=insar[:,2],cmap=cm.magma,s=80,lw=0,vmin=vmin,vmax=vmax)
    plt.legend(['LOS'])
    plt.colorbar()
    gflist.close() 
    plt.show()  
    
    
    
def un_nanify(infile,outfile):
    from numpy import genfromtxt,savetxt,where,nan
    insar=genfromtxt(infile)
    i=where(infile[:,2]!=nan)[0] 
    insar=insar[i,:]
    savetxt(outfile,insar,fmt='%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f')