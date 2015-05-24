format compact

%insar=textread('/Users/dmelgar/Nepal2015/GPS/alos/data.lldtn');
%insar=textread('/Users/dmelgar/Nepal2015/Insar/alos/t157_lltenuds.dat');
insar=textread('/Users/dmelgar/Nepal2015/Insar/raw/t157_resamp_blockmedian_lltenuds.dat');
mindim=32
maxdim=256
thresh=0.19 %0.12 for t047 0.225 for t157
% t047
% lon=insar(:,1);
% lat=insar(:,2);
% los=insar(:,3);
% lookE=insar(:,5);
% lookN=insar(:,6);
% lookU=insar(:,7);

% t147
lon=insar(:,1);
lat=insar(:,2);
los=-insar(:,7)/1000;
lookE=insar(:,4);
lookN=insar(:,5);
lookU=insar(:,6);
%Prepare grid
dlon=8.34e-4;
dlat=8.34e-4;
%Pre filter (T157)
i=find(lat<28.4 & lat>27.3);
lat=lat(i);
lon=lon(i);
los=los(i);
lookN=lookN(i);
lookE=lookE(i);
lookU=lookU(i);
%Pre filter (T048)
% i=find(lat<28.4 & lat>27.3 & lon<86.3 & lon>84.6);
% lat=lat(i);
% lon=lon(i);
% los=los(i);
% lookN=lookN(i);
% lookE=lookE(i);
% lookU=lookU(i);


% lon_i=min(lon):dlon:max(lon);
% lat_i=min(lat):dlat:max(lat);
lon_i=linspace(min(lon),max(lon),2048);
lat_i=linspace(min(lat),max(lat),2048);

display('Gaussian filter')
[X,Y]=meshgrid(lon_i,lat_i);
los_interp_sharp = griddata(lon,lat,los,X,Y);
%Gaussian filter
G = fspecial('gaussian',[30 30],5);
% Filter it
%los_interp = imfilter(los_interp_sharp,G,'same');

%No nan's for detrend
% inan=find(isnan(los_interp_sharp)==1);
% los_nonan=los_interp_sharp;
% los_nonan(inan)=0;
% %Detrend
% los_detrend=detrend_2d(los_nonan);
%Save detrended file
% xsave=reshape(X,2048*2048,1);
% ysave=reshape(Y,2048*2048,1);
% lossave=reshape(los_detrend,2048*2048,1);
% save=[xsave,ysave,lossave];
% dlmwrite('/Users/dmelgar/Nepal2015/Insar/t048_C02_detrended.txt',save,'delimiter','\t','precision',6)
%
los_interp=los_interp_sharp;

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

%Apply line filter
% x1=86.15;
% y1=27.02;
% x2=86.29;
% y2=28.04;
% m=(y1-y2)/(x1-x2);
% b=y1-(m*x1);
% If to the left then keep
% ytest=(m*lon_out)+b;
% i=find(ytest<lat_out);
% lat_out=lat_out(i);
% lon_out=lon_out(i);
% los_out=los_out(i);

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
dlmwrite('/Users/dmelgar/Nepal2015/Insar/t157_qtree.txt',out,'delimiter','\t','precision',6)

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
%xlim([85.7,86.5])
%ylim([26.9,28.1])
title('LOS (m)','FontSize',16)
caxis([-0.7,0.7])
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


