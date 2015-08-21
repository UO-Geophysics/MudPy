% Script to downsample Pisco InSAR data using the quadtree algorithm

format compact

%read data
insar=csvread('/Users/dmelgar/Pisco2007/InSAR/DataFromSAR_latlon.csv');
idecimate=1:100000:length(insar);

lon =insar(idecimate,2);
lat=insar(idecimate,1);
los=insar(idecimate,3);

%QuadTree parameters, these are set by trial and error
mindim=16
maxdim=256
thresh=0.02

lon_i=linspace(min(lon),max(lon),2048);
lat_i=linspace(min(lat),max(lat),2048);

[X,Y]=meshgrid(lon_i,lat_i);
los_interp_sharp = griddata(lon,lat,los,X,Y);

%Gaussian filter
G = fspecial('gaussian',[50 50],50);
los_interp = imfilter(los_interp_sharp,G,'same');

%Now turn things outside the domain to NaN

display('QT decomp')
s=qtdecomp(los_interp,thresh,[mindim,maxdim]);
%Calcualte maximum dimension
maxdim=max(max(full(s)));
%get power of 2 of possible block values
idim=(log(mindim)/log(2)):1:(log(maxdim)/log(2));
%Now loop through
los_out=[]
lon_out=[]
lat_out=[]
for k=1:length(idim)
    [vals, r, c] = qtgetblk(los_interp, s,2^idim(k));
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


%Make look direction vector
lookE_out=ones(length(lon_out))*0.9613;
lookN_out=ones(length(lon_out))*0.2756;
lookU_out=ones(length(lon_out))*0.4528;

%Read path files and keep points only inside the polygons
poly1=dlmread('/Users/dmelgar/Slip_inv/Pisco_insar/data/region_path.txt');
poly2=dlmread('/Users/dmelgar/Slip_inv/Pisco_insar/data/region_path2.txt');

in1=inpolygon(lon_out,lat_out,poly1(:,1),poly1(:,2));
in2=inpolygon(lon_out,lat_out,poly2(:,1),poly2(:,2));
%Make union of two sets of points
in_all=logical(in1+in2);

lon_out=lon_out(in_all);
lat_out=lat_out(in_all);
los_out=los_out(in_all);
lookN_out=lookN_out(in_all);
lookE_out=lookE_out(in_all);
lookU_out=lookU_out(in_all);

%Write to file
out=[lon_out' lat_out' los_out' lookE_out' lookN_out' lookU_out'];
dlmwrite('/Users/dmelgar/Slip_inv/Pisco_insar/data/los_qtree.txt',out,'delimiter','\t','precision',6)

%Plot
figure
subplot(1,2,1)
mesh(X,Y,los_interp)
colorbar
view([0,90])
hold on
scatter3(lon_out,lat_out,100000*ones(size(lat_out)),'kx')
plot3(poly1(:,1),poly1(:,2),ones(length(poly1))*100)
plot3(poly2(:,1),poly2(:,2),ones(length(poly2))*100)
xlabel('Longitude','FontSize',16)
ylabel('Latitude','FontSize',16)
title('LOS (m)','FontSize',16)
caxis([-0.6,0.6])
colormap('jet')
set(gca,'FontSize',14)
subplot(1,2,2)
scatter(lon_out,lat_out,40,los_out,'filled')
colorbar
view([0,90])
caxis([-0.7,0.7])
colormap('jet')
title('QuadTree Resampled LOS (m)','FontSize',16)
xlabel('Longitude','FontSize',16)
ylabel('Latitude','FontSize',16)
%xlim([85.7,86.5])
%ylim([26.9,28.1])
grid on
set(gca,'FontSize',14)
display([int2str(length(lat_out)) ' grids created'])


