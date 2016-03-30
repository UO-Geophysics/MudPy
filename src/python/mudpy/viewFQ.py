'''
Module for plotting FakeQuakes related studd
'''

def plot_LW_scaling(home,project_name,run_name):
    '''
    Plot a comparisson between used length/width and Blasser scaling laws
    '''
    
    from glob import glob
    from numpy import zeros,arange,log10
    from matplotlib import pyplot as plt
    
    logs=glob(home+project_name+'/output/ruptures/'+run_name+'*.log')
    L=zeros(len(logs))
    W=zeros(len(logs))
    Mw=zeros(len(logs))
    for k in range(len(logs)):
        f=open(logs[k],'r')
        while True:
            line=f.readline()
            if 'Lmax' in line:
                L[k]=float(line.split(':')[-1].split()[0])
            if 'Wmax' in line:
                W[k]=float(line.split(':')[-1].split()[0])
            if 'Actual magnitude' in line:
                Mw[k]=float(line.split()[-1])
                break
    
    #Make plot
    plt.figure(figsize=(9,4))
    
    plt.subplot(121)
    Mw_synth=arange(7.8,9.3,0.1)
    plt.plot(Mw_synth,-2.37+0.57*Mw_synth,c='k',lw=2)
    plt.scatter(Mw,log10(L),marker='+')
    plt.xlim([7.75,9.35])
    plt.xlabel('Actual Mw')
    plt.ylabel('log (L) [km]')
    plt.annotate(r'$\log (L)=-2.37+0.57M_w$',xy=(8.0,3.17))
    
    plt.subplot(122)
    plt.plot(Mw_synth,-1.86+0.46*Mw_synth,c='k',lw=2)
    plt.scatter(Mw,log10(W),marker='+')
    plt.xlim([7.75,9.35])
    plt.xlabel('Actual Mw')
    plt.ylabel('log (W) [km]')
    plt.annotate(r'$\log (W)=-1.86+0.46M_w$',xy=(8.0,2.8))
    plt.show()
    
    plt.subplots_adjust(bottom=0.15)
    