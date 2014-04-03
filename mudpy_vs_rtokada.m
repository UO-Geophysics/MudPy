function mudpy_vs_rtokada(tohoku,mudpy,gps)

%Read Mudpy station and find it's corresponding gps name, then see if it
%exists int he rtokada run

for k=1:length(mudpy)
    accel=mudpy{k};
    i=find(strcmp(tohoku.aname,accel));
    mudpy_gps=tohoku.gname(i);
    exists=find(mudpy_gps==gps);
    if isempty(exists) %Found a station that doesn't exist
        display(['Accel station ' accel ' and GPS station ' num2str(mudpy_gps) ' were NOT processed in RTOkada'])
    else
end
    