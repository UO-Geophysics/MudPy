from matplotlib import pyplot as plt
from numpy import genfromtxt,arange

f=genfromtxt('/Users/dmelgar/Slip_inv/test_okada/output/forward_models/_test_okada.grid')
ux=f[:,4]
uy=f[:,3]
uz=f[:,5]
xout=f[:,1]
yout=f[:,2]

plt.figure(figsize=(16,16))
h=(ux**2+uy**2)**0.5
plt.scatter(xout,yout,c=h,lw=0,s=300,vmin=0.0,vmax=0.3)
plt.colorbar()
i=arange(0,len(h),1)
plt.quiver(xout[i],yout[i],ux[i]/h[i],uy[i]/h[i],pivot='mid',linewidths=0.01, edgecolors=('k'),scale=50)
plt.grid()
plt.title('MudPy solution')

plt.figure(figsize=(20,5))

plt.subplot(131)
plt.scatter(xout,yout,c=ux,lw=0,s=200,vmin=-0.3,vmax=0.3)
plt.colorbar()
plt.title('ux')

plt.subplot(132)
plt.scatter(xout,yout,c=uy,lw=0,s=200,vmin=-0.3,vmax=0.3)
plt.colorbar()
plt.title('uy')

plt.subplot(133)
plt.scatter(xout,yout,c=uz,lw=0,s=200,vmin=-0.3,vmax=0.3)
plt.colorbar()
plt.title('uz')

plt.show()