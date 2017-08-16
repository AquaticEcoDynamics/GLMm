function runLA_METAB_halfHour

% runs on top of Run_LA to setup file runs
root = [getenv('USERPROFILE') ...
    '\Desktop\Current Projects\Physics Odom\ChkLakes\highWind\'];
%root = 'C:\Users\Jordan Read\Desktop\Lake Analyzer\G12_ClimPhys\';
lakesAvail = dir(fullfile([root '*bth']));

%lakesAvail = struct('name',{'TroutBog.wtr','CrystalBog.wtr','NorthSparklingBog.wtr'});
numLakes = length(lakesAvail);

startLk = 1;
for j = startLk:numLakes
    charT = lakesAvail(j).name;
    LakeN = charT(1:length(charT)-4);
    Run_LA( LakeN,root(1:end-1),true)
end


end

