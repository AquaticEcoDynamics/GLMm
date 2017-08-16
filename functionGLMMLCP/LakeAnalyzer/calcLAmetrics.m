function LA_metrics = calcLAmetrics(fld_data,varname,GLM_nml,paths)
% function LA_metrics = calcLAmetrics(fld_data,varname,GLM_nml,paths)
%
% Inputs:
%       fld_data   : MATLAB data structure containing observed data
%       varname    : variable names to make calculations of observed Lake
%       Analyser etrics
%       GLM_nml    : MATLAB structure containing name list information for
%       the GLM run including start and finish dates
%
% Outputs:
%       LA_metrics  : MATLAB structure containing series of LA metrics for
%       observed data
%
% Uses:
%      Various Lake Analyzer scripts
%
% Written by L. Bruce 27 February 2014
% Takes field data for variable "varname" and calculated various Lake
% Analyzer metrics of stratification and mixing.
%
%  Saves results to MATLAB file

%Define set of Lake Analyzer metrics to calculate
LA_metrics.LA_metrics = {'all','epi','hyp','thermoD','St','LA'};


%Simulation data
start_date = GLM_nml.start_date;
end_date = GLM_nml.end_date;
start_i = find(fld_data.time >= start_date,1,'first');
end_i = find(fld_data.time <= end_date,1,'last');
LA_metrics.time = fld_data.time(start_i:end_i);


%Hypsographic curve
bthA = fld_data.bthA;
bthD = fld_data.bthD;
LA_metrics.bthA = bthA;
LA_metrics.bthD = bthD;
LA_metrics.depth = fld_data.depth;
windH = 1.0;

%Wind data
%-------------------------------------------------------------------------%
%Read in GLM meteorological information-----------------------------------%
%-------------------------------------------------------------------------%
metData = importdata([paths.working_dir,paths.sim_dir,GLM_nml.metfile],',',1);
metData.time = datenum(metData.textdata(2:end,1));

met_varNames = metData.textdata(1,2:end);
met_varNames = regexprep(met_varNames,' ','');
met_varNames = regexprep(met_varNames,'\t','');
wind_indx = find(strcmp(met_varNames,'WindSpeed')==1);
WindSpeed = metData.data(:,wind_indx);


%Some constant parameters for calculating stratification metrics
drhDz = 0.1;   %min slope for metalimnion (drho/dz per m)
Tdiff = 0.5;    %mixed temp differential (oC)
Tdiff_meta = 0.5; %minimum temperature differential for metalimnion (oC)
Smin = 0.1;      %minimum Salinity
mix_depth_pc = 0.85;  %Percent maximum depth for mixis (if thermoD > mix_depth_pc*max_depth)


%Cut off data to match period of simulation
numDates = end_i - start_i + 1;
numDepths = length(fld_data.depth);
LA_metrics.depth = fld_data.depth;
%Max depth for simulated and observed data
max_depth = max(fld_data.depth);

%Then get simulated data to match field data

%Initialise varaibles
LA_metrics.all = [];
LA_metrics.surf = [];
LA_metrics.bot = [];
LA_metrics.mixed = zeros(numDates,1);
LA_metrics.thermoD = ones(numDates,1)*max_depth;
LA_metrics.meta_top = LA_metrics.thermoD;
LA_metrics.meta_bot = LA_metrics.meta_top;
LA_metrics.thermoInd = ones(numDates,1);
LA_metrics.rho = zeros(numDates,numDepths);
salFld = zeros(1,numDepths);
LA_metrics.St = [];
LA_metrics.Ustar = [];
LA_metrics.LN = [];

for time_i = 1:numDates
    fld_i = time_i + start_i - 1;
    %time_i
    %datestr(LA_metrics.time(time_i))
    %Surface data--------------------------------------------------------%
    %disp('Surface layer');
    LA_metrics.surf(time_i) = fld_data.(varname)(fld_i,1);
    
    %Bottom data--------------------------------------------------------%
    %disp('Bottom layer');
    LA_metrics.bot(time_i) = fld_data.(varname)(fld_i,end);
    
    %All data------------------------------------------------------------%
    %disp('All data');
    for kk = 1:length(fld_data.depth)
        LA_metrics.all(time_i,kk) = fld_data.(varname)(fld_i,kk);
    end
    
    %---------------------------------------------------------------------%
    %Get measures of temperature, salinity, density to calculate
    %stratification metrics
    
    %Observed data
    obs_wtrT = fld_data.temp(fld_i,:);
    LA_metrics.rho(time_i,:) = waterDensity(obs_wtrT,salFld);
    % test shallowest depth with deepest depth (exclude NaNs)
    obs_wtrS = salFld(~isnan(obs_wtrT));
    obs_wtrT = obs_wtrT(~isnan(obs_wtrT));
    % remove NaNs, need at least 3 values
    obs_rhoT = LA_metrics.rho(time_i,:); 
    obs_depT = fld_data.depth';
    obs_depT(isnan(obs_rhoT)) = [];
    obs_rhoT(isnan(obs_rhoT)) = [];
    %Get layer heights
    obs_hgtT = obs_depT(2:end) - obs_depT(1:end-1);
    obs_hgtT(end+1) = obs_hgtT(end);
        
    %Thermocline depth----------------------------------------------------%
    %disp('Finding thermal layers');
    %Observed data
    if abs(obs_wtrT(1)-obs_wtrT(end)) > Tdiff % not mixed... % GIVES mixed if NaN!!!!
         if length(obs_depT)>2
             [LA_metrics.thermoD(time_i),~,drho_dz] = FindThermoDepth(obs_rhoT,obs_depT,Smin);
             LA_metrics.meta_top(time_i) = FindMetaTop(drho_dz,LA_metrics.thermoD(time_i),obs_depT,drhDz);
             LA_metrics.meta_bot(time_i) = FindMetaBot(drho_dz,LA_metrics.thermoD(time_i),obs_depT,drhDz);
         end % or else, keep as NaN
     else
         LA_metrics.mixed(time_i) = 1;
    end
    
    %Determine if a real thermocline depth, i.e. clear stratification with
    %delta T > Tdiff
    %First get temperature at thermocline depth
    thermoT = layerTemp(LA_metrics.thermoD(time_i),LA_metrics.thermoD(time_i),obs_wtrT,obs_depT,bthA,bthD);
    %Determine if stratification or reverse stratification and check if the
    %difference in temperature across the metalimnion is >
    %minimum temprature differential i.e. Tdiff
    if obs_wtrT(1) > thermoT %ordinary stratification
        if abs(thermoT - obs_wtrT(1)) < Tdiff_meta
            LA_metrics.mixed(time_i) = 1;
            LA_metrics.thermoD(time_i) = max_depth;
        end
    else %reverse stratification
        if abs(thermoT - obs_wtrT(end)) < Tdiff_meta
            LA_metrics.mixed(time_i) = 1;
            LA_metrics.thermoD(time_i) = max_depth;
        end
    end        
    
    %Epilimnion temperature (depth average)--------------------------%
    %Include all layers from surface to meta_top layer
    
    %Observed Data
    if LA_metrics.mixed(time_i) == 1 %lake mixed
        LA_metrics.epi(time_i) = layerTemp(0,max(obs_depT),obs_wtrT,obs_depT,bthA,bthD);
    else
        LA_metrics.epi(time_i) = layerTemp(0,LA_metrics.meta_top(time_i),obs_wtrT,obs_depT,bthA,bthD);
    end

    %Hypolimnion temperature (depth average)--------------------------%
    %Include all layers from meta_bot layer to bottom layer
    
    %Observed Data
    if LA_metrics.mixed(time_i) == 1 %lake mixed
        LA_metrics.hyp(time_i) = layerTemp(0,max(obs_depT),obs_wtrT,obs_depT,bthA,bthD);
    else
        LA_metrics.hyp(time_i) = layerTemp(LA_metrics.meta_bot(time_i),max(obs_depT),obs_wtrT,obs_depT,bthA,bthD);
    end
    
    %Schmidt Stability ---------------------------------------------------------%
    %    disp('Calculating Schmidt Stability');
    
    %Observed Data
    LA_metrics.St(time_i) = NaN;
    if length(obs_wtrT) > 2
        LA_metrics.St(time_i) = schmidtStability(obs_wtrT,obs_depT,bthA,bthD,salFld);
    end % else keep as NaN
    
    %U* U_star ---------------------------------------------------------%
    %    disp('Calculating Ustar');
    
    %Observed Data
    [~,wind_i] = min(abs(metData.time - fld_data.time(fld_i)));
    wnd = WindSpeed(wind_i);
    LA_metrics.Ustar(time_i) = NaN;
    if length(obs_wtrT) > 2
        AvEp_rho = layerDensity(obs_depT(1),LA_metrics.meta_top(time_i),obs_wtrT,obs_depT,bthA,bthD); 
        LA_metrics.Ustar(time_i) = uStar(wnd,windH,AvEp_rho);
    end % else keep as NaN       
    
    
    %Lake Number ---------------------------------------------------------%
    %    disp('Calculating Lake Number');
    
    %Observed Data
    LA_metrics.LN(time_i) = NaN;
    if length(obs_wtrT) > 2
        AvHyp_rho = layerDensity(LA_metrics.meta_bot(time_i),max(bthD),obs_wtrT,obs_depT,...
                    bthA,bthD,obs_wtrS);
        LA_metrics.LN(time_i) = lakeNumber(bthA,bthD,LA_metrics.Ustar(time_i),LA_metrics.St(time_i), ...
                                  LA_metrics.meta_top(time_i),LA_metrics.meta_bot(time_i),AvHyp_rho);
    end % else keep as NaN


end

%During mixis set all thermoD to maximum depth
%LA_metrics.thermoD(LA_metrics.thermoD < 0.01*max_depth) = max_depth;
LA_metrics.thermoD(LA_metrics.thermoD > 0.9*max_depth) = max_depth;

%Avoid errors when lake mixed and Schmidt number is close to zero
LA_metrics.St(LA_metrics.St<=0.01) = 0.01;
LA_metrics.thermoD = LA_metrics.thermoD';

%Calculate mixing
LA_metrics.mixed(LA_metrics.thermoD > mix_depth_pc * max_depth) = 1;

%Stratification period ---------------------------------------------------------%
%    disp('Calculating Stratification Period');
    
% Note that this section assumes that the simulation period begins
% during a period of mixis.
% Stratification is defined as a period of more than 30 consecutive days
% of non-mixis
% If lake does not mix before the end of the simulation period then the
% final mixing date is given as the last date of the simulation.
    
%Observed Data
if 0
strat_i = find(obs_data.mixed == 0);
strat_start(1) = LA_metrics.(varname).time(strat_i(1));
st_i = 1;
for ii = 2:length(strat_i)
    if strat_i(ii) > strat_i(ii-1) + 1 %New stratification period
        strat_fini(st_i) = LA_metrics.(varname).time(strat_i(ii-1));
        st_i = st_i +1;
        strat_start(st_i) = LA_metrics.(varname).time(strat_i(ii));
    end
end

%If stratification period continues to the end of the observed data then
%finish on date of last observed data
if length(strat_fini) < length(strat_start)
    strat_fini(length(strat_start)) = LA_metrics.(varname).time(end);
end


%Now remove stratified periods if duration of stratification < 1 month
full_strat = find(strat_fini - strat_start(1:length(strat_fini)) > 30);
strat_start = strat_start(full_strat);
strat_fini  = strat_fini(full_strat);

%For each year find initial and final date of stratification
for year_i = 1:num_years
    start_d = LA_metrics.(varname).time(find(start_date_sim + 365*(year_i-1) <= LA_metrics.(varname).time,1,'first'));
    fini_d  = LA_metrics.(varname).time(find(LA_metrics.(varname).time < start_date_sim + 365*(year_i),1,'last'));
    obs_data.strat_start(year_i) = strat_start(find(strat_start>start_d,1,'first'));
    obs_data.strat_fini(year_i) = strat_fini(find(strat_fini<=fini_d,1,'last'));
end

obs_data.strat_period = obs_data.strat_fini - obs_data.strat_start;
obs_data.strat = [obs_data.strat_start - LA_metrics.(varname).time(1) obs_data.strat_fini - LA_metrics.(varname).time(1) obs_data.strat_period];
end

%Remove any NaN's
%LA_metrics.all = LA_metrics.all(~isnan(LA_metrics.all));
%LA_metrics.surf = LA_metrics.surf(~isnan(LA_metrics.surf));
%LA_metrics.bot = LA_metrics.bot(~isnan(LA_metrics.bot));

%LA_metrics.(varname).obs_data = obs_data;


%Save results
save([paths.working_dir,paths.fld_dir,'LA_metrics.mat'],'LA_metrics')