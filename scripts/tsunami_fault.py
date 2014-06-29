from numpy import genfromtxt,r_,zeros,arange,ones,savetxt,array,tile

nstrike=81
ndip=33
lengthin=25000.
widthin=25000.

f=genfromtxt('/Volumes/Kanagawa/Slip_Inv/tohoku_tsunami/data/model_info/tohoku_dense4.fault')
fref=genfromtxt('/Volumes/Kanagawa/Slip_Inv/tohoku_tsunami/data/model_info/tohoku.fault')

lon=f[:,1]
lat=f[:,2]
depth=f[:,3]
st=f[:,4]
dip=f[:,5]
length=ones(lon.shape)*lengthin
width=ones(lon.shape)*widthin
#Make first row
#First corner
idense=array([1])
iref=array([1])
lonout=lon[0]
latout=lat[0]
depthout=depth[0]
stout=st[0]
dipout=dip[0]
lengthout=length[0]
widthout=width[0]
#Rest of first row
i=arange(3,78,2)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
l1=arange(iref[-1]+1,(iref[-1]+len(i)+1)/2+1)
iref=r_[iref,array([val for pair in zip(l1, l1) for val in pair])]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]/2.]
widthout=r_[widthout,width[i]]
#Last corner
i=80
idense=r_[idense,idense[-1]+1]
iref=r_[iref,iref[-1]+1]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]

#Second row
#First corner
i=324
idense=r_[idense,idense[-1]+1]
iref=r_[iref,iref[-1]+1]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]
# 1st split
i=arange(246,246+nstrike-6,2)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
l1=arange(iref[-1]+1,iref[-1]+(len(i)+1)/2+1)
iref=r_[iref,array([val for pair in zip(l1, l1) for val in pair])]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]/2.]
widthout=r_[widthout,width[i]/2.]
# 2nd split
i=arange(408,408+nstrike-6,2)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
iref=r_[iref,array([val for pair in zip(l1, l1) for val in pair])]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]/2.]
widthout=r_[widthout,width[i]/2.]
#Last corner
i=324+nstrike-1
idense=r_[idense,idense[-1]+1]
iref=r_[iref,iref[-1]+1]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]

#Third row
#First corner
i=648
idense=r_[idense,idense[-1]+1]
iref=r_[iref,iref[-1]+1]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]
# 1st split
i=arange(570,570+nstrike-6,2)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
l1=arange(iref[-1]+1,iref[-1]+(len(i)+1)/2+1)
iref=r_[iref,array([val for pair in zip(l1, l1) for val in pair])]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]/2.]
widthout=r_[widthout,width[i]/2.]
# 2nd split
i=arange(732,732+nstrike-6,2)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
iref=r_[iref,array([val for pair in zip(l1, l1) for val in pair])]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]/2.]
widthout=r_[widthout,width[i]/2.]
#Last corner
i=648+nstrike-1
idense=r_[idense,idense[-1]+1]
iref=r_[iref,iref[-1]+1]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]

#Fourth row
#First corner
i=972
idense=r_[idense,idense[-1]+1]
iref=r_[iref,iref[-1]+1]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]
# 1st split
i=arange(894,894+nstrike-6,2)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
l1=arange(iref[-1]+1,iref[-1]+(len(i)+1)/2+1)
iref=r_[iref,array([val for pair in zip(l1, l1) for val in pair])]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]/2.]
widthout=r_[widthout,width[i]/2.]
# 2nd split
i=arange(1056,1056+nstrike-6,2)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
iref=r_[iref,array([val for pair in zip(l1, l1) for val in pair])]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]/2.]
widthout=r_[widthout,width[i]/2.]
#Last corner
i=972+nstrike-1
idense=r_[idense,idense[-1]+1]
iref=r_[iref,iref[-1]+1]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]

#Fifth row
i=arange(1296,1296+nstrike,4)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
iref=r_[iref,arange(iref[-1]+1,iref[-1]+len(i)+1)]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]

#6th row
i=arange(1620,1620+nstrike,4)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
iref=r_[iref,arange(iref[-1]+1,iref[-1]+len(i)+1)]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]

#7th row
i=arange(1944,1944+nstrike,4)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
iref=r_[iref,arange(iref[-1]+1,iref[-1]+len(i)+1)]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]

#8th row
i=arange(2268,2268+nstrike,4)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
iref=r_[iref,arange(iref[-1]+1,iref[-1]+len(i)+1)]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]

#9th row
i=arange(2592,2592+nstrike,4)
idense=r_[idense,arange(idense[-1]+1,idense[-1]+len(i)+1)]
iref=r_[iref,arange(iref[-1]+1,iref[-1]+len(i)+1)]
lonout=r_[lonout,lon[i]]
latout=r_[latout,lat[i]]
depthout=r_[depthout,depth[i]]
stout=r_[stout,st[i]]
dipout=r_[dipout,dip[i]]
lengthout=r_[lengthout,length[i]]
widthout=r_[widthout,width[i]]

##Save
fout=zeros((len(latout),10))
fout[:,0]=arange(len(latout))+1
fout[:,1]=lonout
fout[:,2]=latout
fout[:,3]=depthout
fout[:,4]=stout
fout[:,5]=dipout
fout[:,6]=ones(len(latout))*0.5
fout[:,7]=ones(len(latout))*10
fout[:,8]=lengthout
fout[:,9]=widthout
#savetxt('/Volumes/Kanagawa/Slip_Inv/tohoku_tsunami/data/model_info/tohoku_tsunami.fault',fout,fmt='%i\t%8.4f\t%8.4f\t%10.6f\t%.2f\t%.2f\t%.4f\t%.3f\t%.4f\t%.4f')

ruptout=zeros((len(latout),14))
ruptout[:,0]=arange(len(latout))+1
ruptout[:,1]=lonout
ruptout[:,2]=latout
ruptout[:,3]=depthout
ruptout[:,4]=stout
ruptout[:,5]=dipout
ruptout[:,6]=ones(len(latout))*0.5
ruptout[:,7]=ones(len(latout))*10
ruptout[:,8]=zeros(len(latout))
ruptout[:,9]=zeros(len(latout))
ruptout[:,10]=lengthout
ruptout[:,11]=widthout
ruptout[:,12]=zeros(len(latout))
ruptout[:,13]=zeros(len(latout))
ruptout=tile(ruptout,(20,1))
#Read rupture modelr esults
rupt=genfromtxt( u'/Volumes/Kanagawa/Slip_Inv/tohoku_10s/output/inverse_models/models/20win_42_fine2.0000.inv')
#Now replace stuff
iref_out=iref.copy()
for k in range(19):
    iref_out=r_[iref_out,iref+189*(k+1)]
iref_out=iref_out-1
ruptout[:,8]=rupt[iref_out,8]
ruptout[:,9]=rupt[iref_out,9]
ruptout[:,12]=rupt[iref_out,12]
ruptout[:,13]=rupt[iref_out,13]
savetxt('/Volumes/Kanagawa/Slip_Inv/tohoku_tsunami/data/model_info/tohoku_dvt.rupt',ruptout,fmt='%i\t%8.4f\t%8.4f\t%10.6f\t%.2f\t%.2f\t%.4f\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%9.4f\t%.4e')


