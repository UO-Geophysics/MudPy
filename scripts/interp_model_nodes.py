from numpy import sin,cos,deg2rad,rad2deg,zeros,arange,arctan,arccos,ones,savetxt,array,vstack
from numpy.linalg import norm
from scipy.io import loadmat

split=2
f=loadmat(u'/Users/dmelgarm/Research/Data/Tohoku/RTOkada/Romano_etal_fault.mat')['F']
R=6371.
L=25.
Lout=L/split
depth=f[:,8]
strike=f[:,9]
dip=f[:,10]
fout=zeros((len(f)*split*split,10))
kout=0
iout=arange(kout,kout+split**2)
#Work one subfault at a time
for kfault in range(len(f)):
    print kfault
    print iout
    #Get strike and dip
    stout=strike[kfault]
    dipout=dip[kfault]
    #Convert 4 coordinates to xyz
    #Top right
    lon=f[kfault,0]
    lat=f[kfault,1]
    z=depth[kfault]
    x1=(R-z)*sin(deg2rad(90-lat))*cos(deg2rad(lon))
    y1=(R-z)*sin(deg2rad(90-lat))*sin(deg2rad(lon))
    z1=(R-z)*cos(deg2rad(90-lat))
    #Top left
    lon=f[kfault,2]
    lat=f[kfault,3]
    z=depth[kfault]
    x2=(R-z)*sin(deg2rad(90-lat))*cos(deg2rad(lon))
    y2=(R-z)*sin(deg2rad(90-lat))*sin(deg2rad(lon))
    z2=(R-z)*cos(deg2rad(90-lat))
    #bottom left
    lon=f[kfault,4]
    lat=f[kfault,5]
    z=depth[kfault]+L*sin(deg2rad(dip[kfault]))
    x3=(R-z)*sin(deg2rad(90-lat))*cos(deg2rad(lon))
    y3=(R-z)*sin(deg2rad(90-lat))*sin(deg2rad(lon))
    z3=(R-z)*cos(deg2rad(90-lat))
    #bottom right
    lon=f[kfault,6]
    lat=f[kfault,7]
    z=depth[kfault]+L*sin(deg2rad(dip[kfault]))
    x4=(R-z)*sin(deg2rad(90-lat))*cos(deg2rad(lon))
    y4=(R-z)*sin(deg2rad(90-lat))*sin(deg2rad(lon))
    z4=(R-z)*cos(deg2rad(90-lat))
    #Define lines and begin splitting
    P1=array([x1,y1,z1])
    P2=array([x2,y2,z2])
    P3=array([x3,y3,z3])
    P4=array([x4,y4,z4])
    L12=array([x1,y1,z1])
    L23=array([x2,y2,z2])
    L34=array([x3,y3,z3])
    L41=array([x4,y4,z4])
    for k in range(split-1):
        #L12
        L=norm(P2-P1)
        deltaL=L/split
        r=(P2-P1)/L
        temp=P1+(deltaL*(k+1))*r
        L12=vstack((L12,temp))
        #L23
        L=norm(P3-P2)
        deltaL=L/split
        r=(P3-P2)/L
        temp=P2+(deltaL*(k+1))*r
        L23=vstack((L23,temp))
        #L34
        L=norm(P4-P3)
        deltaL=L/split
        r=(P4-P3)/L
        temp=P3+(deltaL*(k+1))*r
        L34=vstack((L34,temp))
        #L41
        L=norm(P1-P4)
        deltaL=L/split
        r=(P1-P4)/L
        temp=P4+(deltaL*(k+1))*r
        L41=vstack((L41,temp))
    #Close the quadrangle
    L12=vstack((L12,P2))
    L23=vstack((L23,P3))
    L34=vstack((L34,P4))
    L41=vstack((L41,P1))
    #Get midpoints for all these bad boys
    for i in range(len(L12)-1): #Moves along strike
        for j in range(len(L12)-1): #Moves down dip
            #These are the midpoints along the strike-slip edge
            tempx=(L12[i,0]+L12[i+1,0])/2
            tempy=(L12[i,1]+L12[i+1,1])/2
            tempz=(L12[i,2]+L12[i+1,2])/2
            #This is how we need to move down dip to get to the midpoint
            projectx=(L23[j,0]+L23[j+1,0])/2-L23[0,0]
            projecty=(L23[j,1]+L23[j+1,1])/2-L23[0,1]
            projectz=(L23[j,2]+L23[j+1,2])/2-L23[0,2]
            #Put em together
            mid=array([tempx+projectx,tempy+projecty,tempz+projectz])
            if i==0 and j==0:
                out=mid
            else:
                out=vstack((out,mid))
    #COnvert back to geogrpahical
    lonout=180+rad2deg(arctan(out[:,1]/out[:,0]))
    latout=90-rad2deg(arccos(out[:,2]/((out[:,0]**2+out[:,1]**2+out[:,2]**2)**0.5)))
    zout=R-(out[:,0]**2+out[:,1]**2+out[:,2]**2)**0.5
    #Append to outvarairfare is about 1200iable
    fout[iout,0]=iout+1
    fout[iout,1]=lonout
    fout[iout,2]=latout
    fout[iout,3]=zout
    fout[iout,4]=ones(len(iout))*stout
    fout[iout,5]=ones(len(iout))*dipout
    fout[iout,6]=ones(len(iout))*0.5
    fout[iout,7]=ones(len(iout))*10
    fout[iout,8]=ones(len(iout))*Lout
    fout[iout,9]=ones(len(iout))*Lout
    kout+=1
    iout+=split**2
savetxt('/Users/dmelgarm/Research/Slip_Inv/tohoku_tsunami/data/model_info/tohoku_dense_2split.fault',fout,fmt='%i\t%8.4f\t%8.4f\t%10.6f\t%.2f\t%.2f\t%.4f\t%.3f\t%.4f\t%.4f')