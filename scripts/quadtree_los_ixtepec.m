format compact

fout='/Users/dmelgar/Ixtepec2017/InSAR/T107_Ascending/asc_quadtree.los'
insar=textread('/Users/dmelgarm/Turkey2020/InSAR/asc_los_detrended.txt');
mindim=32
maxdim=256
thresh=0.04 

plot_lims=[-0.3,0.3]

lon=insar(:,1);
lat=insar(:,2);
los=insar(:,7)/1000;
lookE=insar(:,4);
lookN=insar(:,5);
lookU=insar(:,6);

%Prepare grid
dlon=8.34e-4;
dlat=8.34e-4;

lon_i=linspace(min(lon),max(lon),2048);
lat_i=linspace(min(lat),max(lat),2048);


[X,Y]=meshgrid(lon_i,lat_i);
los_interp_sharp = griddata(lon,lat,los,X,Y);

los_interp=los_interp_sharp;

display('QT decomp')
s=qtdecomp(los_interp,thresh,[mindim,maxdim]);
%Calcualte maximum dimension
maxdim=max(max(full(s)));
%get power of 2 of possible block values
idim=(log(mindim)/log(2)):1:(log(maxdim)/log(2));
%Now loop through
los_out=[];
lon_out=[];
lat_out=[];
for k=1:length(idim)
    [vals, r, c] = qtgetblk(los_interp_sharp, s,2^idim(k));
    icurrent=find(s==2^idim(k));
    %Now get mean and cellc enter of each grid
    for kgrid=1:length(icurrent)
       %Get values
       new_los=mean(mean(vals(:,:,kgrid)));
       if ~isnan(new_los)  %Add to the list
           los_out=[los_out,mean(mean(vals(:,:,kgrid)))];
           %get indices of center as upepr left plus half the size of the
           %block
           r1=r(kgrid)+(2^idim(k))/2;
           r2=r(kgrid)+(2^idim(k))/2+1;
           c1=c(kgrid)+(2^idim(k))/2;
           c2=c(kgrid)+(2^idim(k))/2+1;
           %Now figure out coordiantes of the center
           lon_out=[lon_out,0.5*(X(r1,c1)+X(r2,c2))];
           lat_out=[lat_out,0.5*(Y(r1,c1)+Y(r2,c2))];
       end
    end
end



%Find closest look direction vector
lookE_out=[];
lookN_out=[];
lookU_out=[];
for k=1:length(lon_out)
   d=sqrt((lon_out(k)-lon).^2+(lat_out(k)-lat).^2);
   [a,i]=min(d);
   lookE_out=[lookE_out,lookE(i)];
   lookN_out=[lookN_out,lookN(i)];
   lookU_out=[lookU_out,lookU(i)];
end

%Write to file
out=[lon_out' lat_out' los_out' lookE_out' lookN_out' lookU_out'];
dlmwrite(fout,out,'delimiter','\t','precision',6)

%Plot
figure
subplot(1,2,1)
mesh(X,Y,los_interp_sharp)
colorbar
view([0,90])
hold on
scatter3(lon_out,lat_out,100000*ones(size(lat_out)),'kx')
xlabel('Longitude','FontSize',16)
ylabel('Latitude','FontSize',16)

title('LOS (m)','FontSize',16)
caxis(plot_lims)
colormap('jet')
set(gca,'FontSize',14)
subplot(1,2,2)
scatter(lon_out,lat_out,40,los_out,'filled')
colorbar
view([0,90])
caxis(plot_lims)
colormap('jet')
title('QuadTree Resampled LOS (m)','FontSize',16)
xlabel('Longitude','FontSize',16)
ylabel('Latitude','FontSize',16)
grid on
set(gca,'FontSize',14)
display([int2str(length(lat_out)) ' grids created'])


