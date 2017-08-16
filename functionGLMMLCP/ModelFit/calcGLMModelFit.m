function fitparams = calcGLMModelFit(fld_data,simfile,varname,spin_up)
% function fitparams = calcGLMModelFit(fldfile,simfile,varname,spin_up)
%
% Inputs:
%       fld_data   : MATLAB data structure containing observed data
%		simfile    : filename of GLM output
%       varname    : variable names to make calculations of model
%       fit parameters fld = name given to field variable, sim =
%       corresponding name given to GLM simulated data
%       spin_up    : GLM model spin up time to remove from analysis
%
% Outputs
%
% Uses:
%      readGLMnetcdf.m
%
% Written by L. Bruce 19 March 2013
% Takes GLM simulated output and compares against field data for
% variable "varname".
% If plotdata = 'true' plots model against observed for each lake metric
% Calculates the following measures of best fit for varname:
%    Sample size for observed data (N)
%    Mean of observed data (mean)
%    Variance for observed data (var)
%    Normalised mean absolute error (NMAE)
%    Root square mean error (RSME)
%    Normalised mean error (NME)
%    Percent relative error (PRE)
%    Normalised root square mean error (NRSME)
%    Nash-Sutcliffe Model Efficiency (NSE)
%    Correlation coefficient (r), offsett and slope
%  Saves results to .csv file


%Define data subsets to calculate measures of model fit
Data_Subsets = {'all','epi','hyp','thermoD','St'};

%First read field data
%fld_data = readLAwtr(fldfile);
bthA = fld_data.bthA;
bthD = fld_data.bthD;

%Some constant parameters for calculating stratification metrics
drhDz = 0.1;   %min slope for metalimnion (drho/dz per m)
Tdiff = 0.5;    %mixed temp differential (oC)
Tdiff_meta = 0.5; %minimum temperature differential for metalimnion (oC)
Smin = 0.1;      %minimum Salinity
mix_depth_pc = 0.85;  %Percent maximum depth for mixis (if thermoD > mix_depth_pc*max_depth)

%During mixis set all thermoD to maximum depth and make mixed
%fld_data.fld_metrics.mixed(fld_data.fld_metrics.thermoD < 0.01*max(fld_data.depth)) = 1;
fld_data.fld_metrics.mixed(fld_data.fld_metrics.thermoD > 0.9*max(fld_data.depth)) = 1;
%fld_data.fld_metrics.thermoD(fld_data.fld_metrics.thermoD < 0.01*max(fld_data.depth)) = max(fld_data.depth);
fld_data.fld_metrics.thermoD(fld_data.fld_metrics.thermoD > 0.9*max(fld_data.depth)) = max(fld_data.depth);

%Create structure of GLM variable
sim_data = readGLMnetcdf(simfile,{'temp','rho'});

%Save to fit parameters structure
fitparams.(varname).sim_data = sim_data;
fitparams.(varname).fld_data = fld_data;

%Return NaNs if simulation did not run
if length(sim_data.time) <= 1
    for ii = 1:length(Data_Subsets)
        fitparams.(varname).(Data_Subsets{ii}).RMSE = NaN;
        fitparams.(varname).(Data_Subsets{ii}).NSE = NaN;
        fitparams.(varname).(Data_Subsets{ii}).r2 = NaN;
        fitparams.(varname).(Data_Subsets{ii}).PRE = NaN;
    end
else

%Convert simulated z (height from bottom) to depths (mid point for each
%layer)
sim_data.depth = 0.0*sim_data.z - 999;
for time_i = 1:length(sim_data.time)
    max_depth = sim_data.z(time_i,sim_data.NS(time_i));
    sim_data.depth(time_i,1) = max_depth - (sim_data.z(time_i,1))/2;
    for depth_i = 2:sim_data.NS(time_i)
        sim_data.depth(time_i,depth_i) = max_depth - ...
                        (sim_data.z(time_i,depth_i) + sim_data.z(time_i,depth_i-1))/2;
    end
end

%Remove spin up time from analysis
sim_data.startTime = sim_data.startTime + spin_up;
spin_up_idx = find(sim_data.time > sim_data.time(1) + spin_up,1,'first');
sim_data.time = sim_data.time(spin_up_idx:end);
sim_data.NS = sim_data.NS(spin_up_idx:end);

varnames = {'temp','rho','z','depth'};
for ii = 1:length(varnames)
    sim_data.(varnames{ii}) = sim_data.(varnames{ii})(spin_up_idx:end,:);
end
sim_data

%Save to fit parameters structure
fitparams.(varname).sim_data = sim_data;
fitparams.(varname).fld_data = fld_data;

%First determine time array to match field and simulated data
%We have subsampled field data to match end points to simulated data so 
fitparams.(varname).start_date = fld_data.fld_metrics.time(1);
fitparams.(varname).end_date = fld_data.fld_metrics.time(end);    

%Now loop through sim time arrays to get closest match to field data
s_ii=0;
for ii=1:length(fld_data.fld_metrics.time)
    s_ii = s_ii +1;
    [~, sim_i(s_ii)] = min(abs(sim_data.time - fld_data.fld_metrics.time(ii)));
end

%For when field data is measured at a higher frequency than the simulated
%data we need to subsample the field data to match the simulated data
%frequency
sim_i_unique = unique(sim_i);
for ii = 1:length(sim_i_unique)
    [~, fld_i_unique(ii)] = min(abs(sim_data.time(sim_i_unique(ii)) - fld_data.fld_metrics.time));
end

%Number of years
num_years = floor((max(sim_data.time) - min(sim_data.time))/364);
%Number of time and depth samples match simulated data
numDates = length(fld_i_unique);
numDepths = length(fld_data.depth);
%Max depth for simulated and observed data
max_depth = max(fld_data.depth);

%Determine unique field data set
fitparams.(varname).time = fld_data.fld_metrics.time(fld_i_unique);
for kk = 1:length(fld_data.depth)
    fitparams.(varname).data(:,kk) = fld_data.(varname)(fld_i_unique,kk);
end
fitparams.(varname).depth = fld_data.depth;

%Sub set field data to match simulated data
obs_data.surf = fld_data.fld_metrics.surf(fld_i_unique);
obs_data.bot = fld_data.fld_metrics.bot(fld_i_unique);
obs_data.mixed = fld_data.fld_metrics.mixed(fld_i_unique);
obs_data.thermoD = fld_data.fld_metrics.thermoD(fld_i_unique);
obs_data.meta_top = fld_data.fld_metrics.meta_top(fld_i_unique);
obs_data.meta_bot = fld_data.fld_metrics.meta_bot(fld_i_unique);
obs_data.thermoInd = fld_data.fld_metrics.thermoInd(fld_i_unique);
obs_data.rho = fld_data.fld_metrics.rho(fld_i_unique);
obs_data.epi = fld_data.fld_metrics.epi(fld_i_unique);
obs_data.hyp = fld_data.fld_metrics.hyp(fld_i_unique);
obs_data.St = fld_data.fld_metrics.St(fld_i_unique);

%Then get simulated data to match field data

%Initialise varaibles
obs_data.all = [];
sim_data.all = [];
sim_data.surf = [];
sim_data.bot = [];
sim_data.mixed = fld_data.fld_metrics.mixed;
sim_data.thermoD = ones(numDates,1)*max_depth;
sim_data.meta_top = sim_data.thermoD;
sim_data.meta_bot = sim_data.meta_top;
sim_data.thermoInd = ones(numDates,1);
sim_data.rho = zeros(numDates,numDepths);
salFld = zeros(1,numDepths);
sim_data.St = [];

data_i = 1;
for time_i = 1:numDates
    %time_i
    %datestr(fitparams.(varname).time(time_i))
    sim_NS = sim_data.NS(sim_i_unique(time_i));
    %Surface data--------------------------------------------------------%
    %fprintf('Surface layer');
    [~, sim_i_depth] = min(abs(sim_data.depth(sim_i_unique(time_i),:) - fld_data.depth(1)));
    sim_data.surf(time_i) = sim_data.(varname)(sim_i_unique(time_i),sim_i_depth);
    
    %Bottom data--------------------------------------------------------%
    %fprintf('Bottom layer');
    [~, sim_i_depth] = min(abs(sim_data.depth(sim_i_unique(time_i),:) - fld_data.depth(end)));
    sim_data.bot(time_i) = sim_data.(varname)(sim_i_unique(time_i),sim_i_depth);
    
    %All data------------------------------------------------------------%
    %fprintf('All data');
    for kk = 1:length(fld_data.depth)
        obs_data.all(data_i) = fld_data.fld_metrics.all(fld_i_unique(time_i),kk);
        [~, sim_i_depth] = min(abs(sim_data.depth(sim_i_unique(time_i),:) - fld_data.depth(kk)));
        sim_data.all(data_i) = sim_data.(varname)(sim_i_unique(time_i),sim_i_depth);
        %If the difference between observed and simulated depths is <2m
        %then keep data else discard
        if abs(sim_data.depth(sim_i_unique(time_i),sim_i_depth) - fld_data.depth(kk)) < 2
            data_i = data_i + 1;
        end
    end
    
    %---------------------------------------------------------------------%
    %Get measures of temperature, salinity, density to calculate
    %stratification metrics
        
    %Simulated data
    %At this stage assume fresh water as this is how field data metrics are calculated!!!!
    salSim = zeros(1,sim_NS);
    %Water temperature profile
    sim_wtrT = sim_data.temp(sim_i_unique(time_i),1:sim_NS);
    %Water density
    sim_data.rho(time_i,1:length(sim_wtrT)) = waterDensity(sim_wtrT,salSim);
    % test shallowest depth with deepest depth (exclude NaNs)
    sim_wtrT = sim_wtrT(~isnan(sim_wtrT));
    % remove NaNs, need at least 3 values
    sim_rhoT = sim_data.rho(time_i,1:sim_NS); 
    sim_depT = sim_data.depth(sim_i_unique(time_i),1:sim_NS);
    sim_depT(isnan(sim_rhoT)) = [];
    sim_rhoT(isnan(sim_rhoT)) = [];
    sim_wtrT(isnan(sim_rhoT)) = [];
    %Flip simulated arrays so that they go from surface to bottom
    sim_depT = fliplr(sim_depT);
    sim_wtrT = fliplr(sim_wtrT);
    sim_rhoT = fliplr(sim_rhoT);
    %Get layer heights
    sim_hgtT = sim_depT(2:end) - sim_depT(1:end-1);
    sim_hgtT = [sim_depT(1) sim_hgtT];
    
    %Thermocline depth----------------------------------------------------%
    %fprintf('Finding thermal layers');
    
    %Simulated data
    if abs(sim_wtrT(1)-sim_wtrT(end)) > Tdiff % not mixed... % GIVES mixed if NaN!!!!
         if length(sim_depT)>2
             if sim_wtrT(1)>sim_wtrT(end)
                %Tdepth
                [sim_data.thermoD(time_i),~,drho_dz] = FindThermoDepth(sim_rhoT,sim_depT,Smin);
                sim_data.meta_top(time_i) = FindMetaTop(drho_dz,sim_data.thermoD(time_i),sim_depT,drhDz);
                sim_data.meta_bot(time_i) = FindMetaBot(drho_dz,sim_data.thermoD(time_i),sim_depT,drhDz);
             end
         end % or else, keep as NaN
    else
         %thermoD, meta_top and meta_bot stay as max_depth
         sim_data.mixed(time_i) = 1;
    end
    
    %Determine if a real thermocline depth, i.e. clear stratification with
    %delta T > Tdiff
    %First get temperature at thermocline depth
    thermoT = layerTemp(sim_data.thermoD(time_i),sim_data.thermoD(time_i),sim_wtrT,sim_depT,bthA,bthD);
    %Determine if stratification or reverse stratification and check if the
    %difference in temperature across the metalimnion is >
    %minimum temprature differential i.e. Tdiff
    if sim_wtrT(1) > thermoT %ordinary stratification
        if abs(thermoT - sim_wtrT(1)) < Tdiff_meta
            sim_data.mixed(time_i) = 1;
            sim_data.thermoD(time_i) = max_depth;
        end
    else %reverse stratification
        if abs(thermoT - sim_wtrT(end)) < Tdiff_meta
            sim_data.mixed(time_i) = 1;
            sim_data.thermoD(time_i) = max_depth;
        end
    end        
   
    %During mixis set all thermoD to maximum depth and make mixed
    %if (sim_data.thermoD(time_i) < 0.01*max_depth) || (sim_data.thermoD(time_i) > 0.9*max_depth)
    %    sim_data.mixed(time_i) = 1;
    %    sim_data.thermoD(time_i) = max_depth;
    %end
        
    %sim_data.mixed(sim_data.thermoD < 0.01*max_depth) = 1;
    %sim_data.mixed(sim_data.thermoD > 0.9*max_depth) = 1;
    %sim_data.thermoD(sim_data.thermoD < 0.01*max_depth) = max_depth;
    %sim_data.thermoD(sim_data.thermoD > 0.9*max_depth) = max_depth;
    
    %Epilimnion temperature (depth average)--------------------------%
    %Include all layers from surface to meta_top layer
    
    %Simulated Data
    if sim_data.mixed(time_i) == 1 %lake mixed
        sim_data.epi(time_i) = layerTemp(0,max(sim_depT),sim_wtrT,sim_depT,bthA,bthD);
    else
        sim_data.epi(time_i) = layerTemp(0,sim_data.meta_top(time_i),sim_wtrT,sim_depT,bthA,bthD);
    end


    %Hypolimnion temperature (depth average)--------------------------%
    %Include all layers from meta_bot layer to bottom layer
    
    %Simulated Data
    if sim_data.mixed(time_i) == 1 %lake mixed
        sim_data.hyp(time_i) = layerTemp(0,max(sim_depT),sim_wtrT,sim_depT,bthA,bthD);
    elseif (sim_wtrT(1) < sim_wtrT(end)) %Inverse stratification, consider mixed
        sim_data.hyp(time_i) = layerTemp(0,max(sim_depT),sim_wtrT,sim_depT,bthA,bthD);
    else
        sim_data.hyp(time_i) = layerTemp(sim_data.meta_bot(time_i),max(sim_depT),sim_wtrT,sim_depT,bthA,bthD);
    end                 
                         
    %Schmidt Stability ---------------------------------------------------------%
    %    fprintf('Calculating Schmidt Stability');
    
    %Simulated Data
    sim_data.St(time_i) = NaN;
    if sim_data.mixed == 1
        sim_data.St(time_i) = NaN;
    elseif length(sim_wtrT) > 2
        sim_data.St(time_i) = schmidtStability(sim_wtrT,sim_depT,fld_data.bthA,fld_data.bthD,salSim);
    end % else keep as NaN
    
end

%Avoid errors when lake mixed and Schmidt number is close to zero
sim_data.St(sim_data.St<=0.01) = 0.01;

sim_data.thermoD = sim_data.thermoD';

%Calculate mixing
sim_data.mixed(sim_data.thermoD > mix_depth_pc * max_depth) = 1;

%Stratification period ---------------------------------------------------------%
%    fprintf('Calculating Stratification Period');
    
% Note that this section assumes that the simulation period begins
% during a period of mixis.
% Stratification is defined as a period of more than 30 consecutive days
% of non-mixis
% If lake does not mix before the end of the simulation period then the
% final mixing date is given as the last date of the simulation.
    
%Observed Data
if 0
strat_i = find(obs_data.mixed == 0);
strat_start(1) = fitparams.(varname).time(strat_i(1));
st_i = 1;
for ii = 2:length(strat_i)
    if strat_i(ii) > strat_i(ii-1) + 1 %New stratification period
        strat_fini(st_i) = fitparams.(varname).time(strat_i(ii-1));
        st_i = st_i +1;
        strat_start(st_i) = fitparams.(varname).time(strat_i(ii));
    end
end

%If stratification period continues to the end of the observed data then
%finish on date of last observed data
if length(strat_fini) < length(strat_start)
    strat_fini(length(strat_start)) = fitparams.(varname).time(end);
end


%Now remove stratified periods if duration of stratification < 1 month
full_strat = find(strat_fini - strat_start(1:length(strat_fini)) > 30);
strat_start = strat_start(full_strat);
strat_fini  = strat_fini(full_strat);

%For each year find initial and final date of stratification
for year_i = 1:num_years
    start_d = fitparams.(varname).time(find(start_date_sim + 365*(year_i-1) <= fitparams.(varname).time,1,'first'));
    fini_d  = fitparams.(varname).time(find(fitparams.(varname).time < start_date_sim + 365*(year_i),1,'last'));
    obs_data.strat_start(year_i) = strat_start(find(strat_start>start_d,1,'first'));
    obs_data.strat_fini(year_i) = strat_fini(find(strat_fini<=fini_d,1,'last'));
end

obs_data.strat_period = obs_data.strat_fini - obs_data.strat_start;
obs_data.strat = [obs_data.strat_start - fitparams.(varname).time(1) obs_data.strat_fini - fitparams.(varname).time(1) obs_data.strat_period];

%Simulated Data
clear strat_fini strat_start
strat_i = find(sim_data.mixed == 0);
strat_start(1) = fitparams.(varname).time(strat_i(1));
st_i = 1;
for ii = 2:length(strat_i)
    if strat_i(ii) > strat_i(ii-1) + 1 %New stratification period
        strat_fini(st_i) = fitparams.(varname).time(strat_i(ii-1));
        st_i = st_i +1;
        strat_start(st_i) = fitparams.(varname).time(strat_i(ii));
    end
end

%Temporary code for GLM_MLCP where lakes simulated 2 years
%If still stratified at the end of the simulation then use the last date as
%the finish of stratification
if length(strat_fini) < length(strat_start); strat_fini(st_i) = fitparams.(varname).time(end); end
    

%Now remove stratified periods if duration of stratification < 1 month
full_strat = find(strat_fini - strat_start(1:length(strat_fini)) > 30);
strat_start = strat_start(full_strat);
strat_fini  = strat_fini(full_strat);

%For each year find initial and final date of stratification
for year_i = 1:num_years
    start_d = fitparams.(varname).time(find(start_date_sim + 365*(year_i-1) <= fitparams.(varname).time,1,'first'));
    fini_d  = fitparams.(varname).time(find(fitparams.(varname).time < start_date_sim + 365*(year_i),1,'last'));
    sim_data.strat_start(year_i) = strat_start(find(strat_start>start_d,1,'first'));
    sim_data.strat_fini(year_i) = strat_fini(find(strat_fini<=fini_d,1,'last'));
end

sim_data.strat_period = sim_data.strat_fini - sim_data.strat_start;
sim_data.strat = [sim_data.strat_start - fitparams.(varname).time(1) sim_data.strat_fini - fitparams.(varname).time(1) sim_data.strat_period];        
end %If calculating simulation period

%Remove any NaN's
noNaN_all = find(~isnan(obs_data.all));
obs_data.all = obs_data.all(noNaN_all);
sim_data.all = sim_data.all(noNaN_all);
noNaN_surf = find(~isnan(obs_data.surf));
obs_data.surf = obs_data.surf(noNaN_surf);
sim_data.surf = sim_data.surf(noNaN_surf);
noNaN_bot = find(~isnan(obs_data.bot));
obs_data.bot = obs_data.bot(noNaN_bot);
sim_data.bot = sim_data.bot(noNaN_bot);

%Remove and zero observed data (assume measurement not working)
%obs_data.all = obs_data.all(obs_data.all > 0.0001);
%sim_data.all = sim_data.all(obs_data.all > 0.0001);
%obs_data.surf = obs_data.surf(obs_data.surf > 0.0001);
%sim_data.surf = sim_data.surf(obs_data.surf > 0.0001);
%obs_data.bot = obs_data.bot(obs_data.bot > 0.0001);
%sim_data.bot = sim_data.bot(obs_data.bot > 0.0001);


fitparams.(varname).obs_data = obs_data;
fitparams.(varname).sim_data = sim_data;

%Calculate measure of model fit------------------------------------------%

for ii = 1:length(Data_Subsets)
   % ii
    
   %First check maximum size is first element in size top avoid large arrays
   %of measures of model fit
   size_data = size(obs_data.(Data_Subsets{ii}));
   if size_data(2) > size_data(1)
       obs_data.(Data_Subsets{ii}) = obs_data.(Data_Subsets{ii})';
   end
   size_data = size(sim_data.(Data_Subsets{ii}));
   if size_data(2) > size_data(1)
       sim_data.(Data_Subsets{ii}) = sim_data.(Data_Subsets{ii})';
   end
    
   %Number of samples
   fitparams.(varname).(Data_Subsets{ii}).N = length(obs_data.(Data_Subsets{ii}));
    
   %Calculate mean
   fitparams.(varname).(Data_Subsets{ii}).mean = mean(obs_data.(Data_Subsets{ii}));

   %Calculate variance
   fitparams.(varname).(Data_Subsets{ii}).var = var(obs_data.(Data_Subsets{ii}));

   %Calculate NMAE
   fitparams.(varname).(Data_Subsets{ii}).NMAE = sum(abs(sim_data.(Data_Subsets{ii}) - obs_data.(Data_Subsets{ii}))) ...
                                        /length(obs_data.(Data_Subsets{ii})) /mean(obs_data.(Data_Subsets{ii}));
   
   %Calculate RMSE
   fitparams.(varname).(Data_Subsets{ii}).RMSE = rms(sim_data.(Data_Subsets{ii}) - obs_data.(Data_Subsets{ii}), 1);

   %Calculate NME
   fitparams.(varname).(Data_Subsets{ii}).NME = mean(sim_data.(Data_Subsets{ii}) - obs_data.(Data_Subsets{ii})) /mean(obs_data.(Data_Subsets{ii}));

   %Calculate percent relative error PRE
   fitparams.(varname).(Data_Subsets{ii}).PRE = mean((sim_data.(Data_Subsets{ii}) - obs_data.(Data_Subsets{ii}))) /mean(obs_data.(Data_Subsets{ii}))*100;

   %Calculate NRSME
   fitparams.(varname).(Data_Subsets{ii}).NRSME = 1 - sum((sim_data.(Data_Subsets{ii}) - obs_data.(Data_Subsets{ii})).^2) /sum((obs_data.(Data_Subsets{ii}) - mean(obs_data.(Data_Subsets{ii}))).^2);

   %Calculate Nash-Sutcliffe Model Efficiency NSE
   fitparams.(varname).(Data_Subsets{ii}).NSE = 1 - mean((sim_data.(Data_Subsets{ii}) - obs_data.(Data_Subsets{ii})).^2) / ...
                                                    mean((obs_data.(Data_Subsets{ii}) - mean(obs_data.(Data_Subsets{ii}))).^2);
                
   %Calculate offset and slope
   [fitparams.(varname).(Data_Subsets{ii}).r2 fitparams.(varname).(Data_Subsets{ii}).slope ...
       fitparams.(varname).(Data_Subsets{ii}).offset] = regression( ...
       obs_data.(Data_Subsets{ii})',sim_data.(Data_Subsets{ii})');


end %Loop through data subsets

end %If simulation crashed

%Print table of lake metrics for each measure of model fit


% Open file
% ---------------------------- DEFINE THE FILE --------------------------- %

% Create .csv file.
fid = fopen('GLM_fitparams.csv','w');

% Write header information

header_line = 'Lake_Metric, RMSE, NSE, r2, PRE, NMAE';

fprintf(fid,'%s \n', header_line);


% ---------------------------- STORE THE DATA ---------------------------- %

format_line = '%s';
for ii = 1:5
    format_line = [format_line,', %f '];
end
format_line = [format_line,'\n'];

for ii = 1:length(Data_Subsets)
    %ii
    name_str = Data_Subsets{ii};
    data_line = [fitparams.(varname).(Data_Subsets{ii}).RMSE, ...
                 fitparams.(varname).(Data_Subsets{ii}).NSE, ...
                 fitparams.(varname).(Data_Subsets{ii}).r2, ...
                 fitparams.(varname).(Data_Subsets{ii}).PRE, ...
                 fitparams.(varname).(Data_Subsets{ii}).NMAE];
    fprintf(fid,format_line,name_str,data_line);
end 

fclose(fid);                                      % Close the file.


% ------------------------------------------------------------------------ %
