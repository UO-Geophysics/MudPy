%Make regions for filtering out insar points

%Specify path and filename for saving
fout='/Users/dmelgar/Slip_inv/Pisco_insar/data/region_path2.txt'

%read data
%insar=csvread('/Users/dmelgar/Pisco2007/InSAR/DataFromSAR_latlon.csv');
idecimate=1:100:length(insar);
lon =insar(idecimate,2);
lat=insar(idecimate,1);
los=insar(idecimate,3);

%Plot
scatter(lon,lat,10,los,'filled')

%Start interactive tool
h=impoly()

%Get selected path coordiantes
nodes=getPosition(h);

%Save to file
dlmwrite(fout,nodes,'delimiter','\t','precision',6)