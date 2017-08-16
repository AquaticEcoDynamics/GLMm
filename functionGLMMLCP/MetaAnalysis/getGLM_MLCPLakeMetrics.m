function Lake = getGLM_MLCPLakeMetrics(basedir,LakeName)
%function Lake = getGLM_MLCPLakeMetrics(basedir,LakeName)
%
% Inputs:
%      basedir :  base directory
%      LakeName:  Lake name
%
% Outputs:
%     Lake: MATLAB structure containing lake metrics for LakeName such as
%     lake volume, depth, residence time, Kw etc
%
% Uses:
%      
%
% Written by L. Bruce 18 December 2013
%
% Reads GLM input files to get a series of lake characteristic metrics

%-------------------------------------------------------------------------%
%Read in GLM name list information----------------------------------------%
%-------------------------------------------------------------------------%
Lake.nml = getGLMnml([basedir,LakeName,'/nml/glm2_init.nml']);

%-------------------------------------------------------------------------%
%Read in GLM meteorological information-----------------------------------%
%-------------------------------------------------------------------------%
Lake.met = importdata([basedir,LakeName,'/sim/',Lake.nml.metfile],',',1);
Lake.met.time = datenum(Lake.met.textdata(2:end,1));
%Get simulation time index
met_start_i = find(Lake.met.time >= Lake.nml.start_date,1,'first');
met_end_i = find(Lake.met.time <= Lake.nml.end_date,1,'last');

Lake.met.varNames = Lake.met.textdata(1,2:end);
for var_i = 1:length(Lake.met.varNames)
    Lake.met.varNames{var_i} = regexprep(Lake.met.varNames{var_i},' ','');
    Lake.met.varNames{var_i} = regexprep(Lake.met.varNames{var_i},'\t','');
    Lake.met.(Lake.met.varNames{var_i}) = Lake.met.data(met_start_i:met_end_i,var_i);
end

%-------------------------------------------------------------------------%
%Read in GLM flow and volume information----------------------------------%
%-------------------------------------------------------------------------%
Lake.daily = importdata([basedir,LakeName,'/sim/',Lake.nml.lakefile,'.csv'],',',1);
Lake.daily.time = datenum(Lake.daily.textdata(2:end,1));

%create lake structure from data in lake.csv file
varnames = Lake.daily.textdata(1,2:end);
varnames = regexprep(varnames,' ','_');
varnames = regexprep(varnames,'/','_');
for ii = 1:length(varnames)
    Lake.daily.(varnames{ii}) = Lake.daily.data(:,ii);
end

%-------------------------------------------------------------------------%
% 1. Size metrics
%-------------------------------------------------------------------------%

%1.1 Volume
Lake.Volume = mean(Lake.daily.Volume);

%1.2 Surface Area
Lake.Area = mean(Lake.daily.Surface_Area);

%1.3 Depth
Lake.Depth = mean(Lake.daily.Lake_Level);


%1.4 Area/Depth
Lake.AonD = Lake.nml.A(end) / Lake.Depth;

%1.5 Length/Width
Lake.LonW = Lake.nml.bsn_len/Lake.nml.bsn_wid;

%-------------------------------------------------------------------------%
% 2. Flow metrics
%-------------------------------------------------------------------------%

%2.1 Mean inflow
Lake.Inflow = mean(Lake.daily.Tot_Inflow_Vol);

%2.2 Mean residence time (days)
%Lake.ResTime = mean(Lake.daily.Volume(Lake.daily.Tot_Inflow_Vol>0) ./Lake.daily.Tot_Inflow_Vol(Lake.daily.Tot_Inflow_Vol>0));
% L. Bruce 16 February 2016 - changed to average ResTime (mean vol/mean
% inflow)
%If zero flow make NaN
if Lake.Inflow == 0
    Lake.ResTime = NaN;
else
    Lake.ResTime = Lake.Volume/Lake.Inflow;
end


%-------------------------------------------------------------------------%
% 3. Climate metrics
%-------------------------------------------------------------------------%

%3.1 Mean short wave radiation (W/m2)
Lake.ShortWave = mean(Lake.met.ShortWave);

%3.2 Mean air temperature (oC)
Lake.AirTemp = mean(Lake.met.AirTemp);

%3.2 Mean wind speed (m/2)
Lake.WindSpeed = mean(Lake.met.WindSpeed) * Lake.nml.wind_factor;

%-------------------------------------------------------------------------%
% 4. Miscelaneous
%-------------------------------------------------------------------------%

%4.1 Extinction coefficient
Lake.Kw = Lake.nml.Kw;

%4.2 Latitude
Lake.Lat = Lake.nml.Lat;

%4.3 Lake Number
%Remove negative or extra large Lake Numbers
LakeNumber = Lake.daily.LakeNumber;
LakeNumber = LakeNumber(LakeNumber>0);
LakeNumber = LakeNumber(LakeNumber<1e5);
Lake.LN = mean(LakeNumber);



