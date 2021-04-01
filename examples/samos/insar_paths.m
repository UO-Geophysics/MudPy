%Make regions for filtering out insar points

%Specify path and filename for saving
fout='/Users/dmelgarm/Samos2020/InSAR/T131/samos_path.txt'

%read data
insar=textread('/Users/dmelgarm/Samos2020/InSAR/T131/T131_asc.lltnde');
idecimate=1:4:length(insar);
lon =insar(idecimate,1);
lat=insar(idecimate,2);
los=insar(idecimate,7);

%Plot
scatter(lon,lat,10,los,'filled')
xlim([26.4,27.4])
ylim([37.5,37.9])
colorbar

%Start interactive tool
h=impoly()

%Get selected path coordiantes
nodes=getPosition(h);

%Save to file
dlmwrite(fout,nodes,'delimiter','\t','precision',6)