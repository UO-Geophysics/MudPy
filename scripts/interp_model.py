from numpy import genfromtxt,sin,cos,deg2rad,rad2deg,zeros,arange,arctan,arccos,ones,savetxt

f=genfromtxt('/Volumes/Kanagawa/Slip_Inv/tohoku_tsunami/data/model_info/tohoku_dense2.fault')
lon=f[:,1]
lat=f[:,2]
depth=f[:,3]
st=f[:,4]
dip=f[:,5]

nstrike=41
ndip=17
ninterp=2
Ndip=2*ndip-1
Nstrike=2*nstrike-1

#Convert everything to xyz
R=6371
x=(R-depth)*sin(deg2rad(90-lat))*cos(deg2rad(lon))
y=(R-depth)*sin(deg2rad(90-lat))*sin(deg2rad(lon))
z=(R-depth)*cos(deg2rad(90-lat))

xout1=zeros(Nstrike*ndip)
yout1=xout1.copy()
zout1=xout1.copy()
stout1=xout1.copy()
dipout1=xout1.copy()
xout2=zeros(Nstrike*Ndip)
yout2=xout2.copy()
zout2=xout2.copy()
stout2=xout2.copy()
dipout2=xout2.copy()
#traverse along strike
inode=arange(nstrike)*2
imid=arange(nstrike-1)*2+1
iread=arange(nstrike)
for k in range(ndip):
    xout1[inode]=x[iread]
    yout1[inode]=y[iread]
    zout1[inode]=z[iread]
    stout1[inode]=st[iread]
    dipout1[inode]=dip[iread]
    xout1[imid]=(x[iread[0:-1]]+x[iread[1:]])/2
    yout1[imid]=(y[iread[0:-1]]+y[iread[1:]])/2
    zout1[imid]=(z[iread[0:-1]]+z[iread[1:]])/2
    stout1[imid]=(st[iread[0:-1]]+st[iread[1:]])/2
    dipout1[imid]=(dip[iread[0:-1]]+dip[iread[1:]])/2
    inode+=Nstrike
    iread+=nstrike
    imid+=Nstrike
#Traverse along dip
iread=arange(Nstrike)
inodes=arange(Nstrike)
imid=arange(Nstrike)+Nstrike
xout2[inodes]=xout1[iread]
yout2[inodes]=yout1[iread]
zout2[inodes]=zout1[iread]
stout2[inodes]=stout1[iread]
dipout2[inodes]=dipout1[iread]
for k in range(ndip-1):
    inodes+=Nstrike*2
    xout2[imid]=(xout1[iread]+xout1[iread+Nstrike])/2
    yout2[imid]=(yout1[iread]+yout1[iread+Nstrike])/2
    zout2[imid]=(zout1[iread]+zout1[iread+Nstrike])/2
    stout2[imid]=(stout1[iread]+stout1[iread+Nstrike])/2
    dipout2[imid]=(dipout1[iread]+dipout1[iread+Nstrike])/2
    imid+=Nstrike*2
    
    xout2[inodes]=xout1[iread+Nstrike]
    yout2[inodes]=yout1[iread+Nstrike]
    zout2[inodes]=zout1[iread+Nstrike]
    stout2[inodes]=stout1[iread+Nstrike]
    dipout2[inodes]=dipout1[iread+Nstrike]
    iread+=Nstrike

#Convert back to geographical
lonout=180+rad2deg(arctan(yout2/xout2))
latout=90-rad2deg(arccos(zout2/((xout2**2+yout2**2+zout2**2)**0.5)))
zout=R-(xout2**2+yout2**2+zout2**2)**0.5
#Save
fout=zeros((len(zout),10))
fout[:,0]=arange(len(zout))+1
fout[:,1]=lonout
fout[:,2]=latout
fout[:,3]=zout
fout[:,4]=stout2
fout[:,5]=dipout2
fout[:,6]=ones(len(zout))*0.5
fout[:,7]=ones(len(zout))*10
fout[:,8]=ones(len(zout))*12500
fout[:,9]=ones(len(zout))*12500

savetxt('/Volumes/Kanagawa/Slip_Inv/tohoku_tsunami/data/model_info/tohoku_dense4.fault',fout,fmt='%i\t%8.4f\t%8.4f\t%10.6f\t%.2f\t%.2f\t%.4f\t%.3f\t%.4f\t%.4f')