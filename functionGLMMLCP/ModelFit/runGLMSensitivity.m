
function sens_anal = runGLMSensitivity(lakename,conf_file)
% Inputs:
%       fld_data   : MATLAB data structure containing observed data
%		params     : list of parameters with initial values
%
% Outputs:
%       sens_anal = model sensitivity to measure of model fit for
%                    observed vs simulated data of data subset (eg all, 
%                    surface,  hypolimnion temp, Schmidt stability etc
%
% Uses:
%     runGLMModelFit.m
%
% Written by L. Bruce 22 April 2013
% Runs through a set of physical parameters, adjusting each +/- 20% then
% running GLM and comparing how sensitive measures of model fit are to 
% changes in GLM parameters.
%

%Start counting run time
tic;

%CONFIGURATION INFORMATION------------------------------------------------%

conf = readGLMconfig(conf_file);

%Working directory to include Lake Name
conf.paths.working_dir = [conf.paths.base_dir,lakename,'/',conf.paths.working_dir];
base_dir = conf.paths.base_dir;
working_dir = conf.paths.working_dir;

%Config for this lake
conf.paths.lake_dir = working_dir;

%Create Output and Results folders
mkdir([working_dir,'/Output']);
mkdir([working_dir,'/Results']);
mkdir([working_dir,'/Results/Plots']);

%Move to working directory
paths = conf.paths;
cd (working_dir)

%If running remotely then add ganymed-ssh2-build250 directory'
if conf.config.remote_flag == 1
    
   ganymed_dir = '/Users/hydro/Work/MATLAB/functions/functionMisc/ganymed-ssh2-build250/';
    copyfile(ganymed_dir,[working_dir,'/ganymed-ssh2-build250/']);

   %Some information for ssh connection
   conf.run_info.host_name = 'aqua';
   conf.run_info.usr_name = 'hydro';
   conf.run_info.password = 'modeller';
   conf.run_info.remote_dir = ['Louise/GLM_MLCP/v2.0.0/',lakename,'/sim'];
   conf.run_info.run_glm = ['cd ',conf.run_info.remote_dir,'; ./glm'];
   conf.run_info.output_file = [conf.run_info.remote_dir,'/output.nc'];
end

%Other MCMC configuration information
varname = conf.config.varname;

%GLM configuration for base simulation
glm_nml = getGLMnml([working_dir,'nml/glm2_init.nml']);
sim_start = glm_nml.start_date;
sim_stop = glm_nml.end_date;
sim_time = sim_start:glm_nml.out_freq:sim_stop;
max_depth_plots = glm_nml.max_depth_plots;


%List of data subsets to calculate sensitivity to
Data_Subsets = conf.dataset.Data_Subsets;
%List of measures of model fit
Model_Fit = conf.dataset.Model_Fit;


%FIELD DATA--------------------------------------------------------------%

%Field data file
fld_temp_file = [working_dir,conf.paths.fld_dir,lakename,'_Fld_temp.wtr'];

%First read in field data
fld_data = readLAwtr(fld_temp_file);

%Get bathymetry data from glm namelist HA curve
fld_data.bthD = flipud(glm_nml.H(end) - glm_nml.H)';
fld_data.bthA = flipud(glm_nml.A)';

%Now plot field data as colour contour
figure
varInformation
x_tick = sim_start:round((sim_start-sim_stop)/4):sim_stop;
x_tick_date = datestr(x_tick,'dd-mmm-yy');
pcolor(fld_data.time,fld_data.depth,fld_data.temp');
 axis ij
shading interp
 xlabel('Date')
 ylabel('Depth (m)')
set(gca,'XTick',x_tick,'XTickLabel',x_tick_date)
var_lim = varIndex.temp.caxis;
var_title = varIndex.temp.title;
 title(var_title)
 caxis(var_lim)
 colorbar

%Save plot of field temperature into Field file
xlim([sim_start sim_stop])
ylim([0 max_depth_plots])
fig_name = [working_dir,'Results/Plots/Field_Temp','.png'];
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'paperposition',[0.1  0.1 15 8]);
print(gcf,'-dpng',fig_name,'-opengl');

close(gcf);

%Calculate thermodynamics metrics using Lake Analyzer scripts for observed
%data.  Note includes metrics such as Lake Number that do not appear in the
%data subset listed in the MCMC configuration file.  Also calculated for
%all observed data bounded by the simulation specified start and stop
%dates.
fld_data.fld_metrics = calcLAmetrics(fld_data,varname,glm_nml,paths);

%INITIAL RUN-------------------------------------------------------------%

%List of parameters and their initial values
pars = fieldnames(conf.params);
for i = 1:length(pars)
    params{i} = {pars{i},conf.params.(pars{i})};
end

theta_initial = zeros(length(params),1);
for param_i = 1:length(params)
    theta_initial(param_i) = params{param_i}{2};
end

%Calculate GLM model fit parameters using initial run
%And save results
fit_params_init = runGLMModelFit(theta_initial,fld_data,params,'init',conf);

%Calculate Thermodynamic Metrics for simulation file
sens_anal.init.sim_metrics = calcLaGLMmetrics('Output/output_init.nc',conf);


%Plot initial run and comparison against field - save results
plotLake(paths,'output_init.nc')
xlim([sim_start sim_stop]);
%x_tick = [sim_start:round((sim_start-sim_stop)/4):sim_stop];
x_tick_i = get(gca,'XTick');
x_tick = sim_start:(sim_stop - sim_start)/4:sim_stop;
set(gca,'XTick',x_tick,'XTickLabel',datestr(x_tick,'mm-yy'))

ylim([0 max_depth_plots])
fig_name = [working_dir,'Results/Plots/Sim_Temp','.png'];
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'paperposition',[0  0 15 8]);
print(gcf,'-dpng',fig_name,'-opengl');


%INITIAL RUN COMPARISON PLOTS --------------------------------------------%

%Now plot scatter and time series comparisons

plotGLMModelFit(fit_params_init.(varname),[paths.working_dir,'Results/Plots/'])

close all;


%SENSITIVITY ANALYSIS-----------------------------------------------------%

%Now run through list of parameters and calculate indicators of model fit

%Set parameter bounds as +/- 20% of initial value
for param_i = 1:length(params)
    params{param_i}{3} = 0.8*params{param_i}{2};
    params{param_i}{4} = 1.2*params{param_i}{2};
end

%Loop through parameters and determine sensitivity to model performance for each
%parameter

%Note that this loop used to save all output netcdf files, changed to init
%so overide and not fill disk

for param_i = 1:length(params)
    %Initialies parameters to initial
    theta_lower = theta_initial;
    theta_upper = theta_initial;
    param_name = params{param_i}{1};
    disp(['Parameter Name:  ',param_name]);
    theta_lower(param_i) = params{param_i}{3};
    theta_upper(param_i) = params{param_i}{4};
    %sens_anal.(param_name).fitparams_lower = runGLMModelFit(theta_lower,fld_data,params,[param_name,'_lower'],conf);
    %sens_anal.(param_name).fitparams_upper = runGLMModelFit(theta_upper,fld_data,params,[param_name,'_upper'],conf);
    sens_anal.(param_name).fitparams_lower = runGLMModelFit(theta_lower,fld_data,params,'new',conf);
    %Calculate Thermodynamic Metrics for simulation file
    sens_anal.(param_name).sim_metrics_lower = calcLaGLMmetrics('Output/output_new.nc',conf);
    sens_anal.(param_name).fitparams_upper = runGLMModelFit(theta_upper,fld_data,params,'new',conf);
    %Calculate Thermodynamic Metrics for simulation file
    sens_anal.(param_name).sim_metrics_upper = calcLaGLMmetrics('Output/output_new.nc',conf);
end

%Tempory code for mac users when runs fail with NaN results
num_nan=1;
while num_nan > 0
    num_nan = 0;
    for param_i = 1:length(params)
        %Initialies parameters to initial
        theta_lower = theta_initial;
        theta_upper = theta_initial;
        param_name = params{param_i}{1};
        theta_lower(param_i) = params{param_i}{3};
        theta_upper(param_i) = params{param_i}{4};
        if isnan(sens_anal.(param_name).fitparams_lower.(varname).all.RMSE)
            sens_anal.(param_name).fitparams_lower = runGLMModelFit(theta_lower,fld_data,params,[param_name,'_lower'],conf);
        end
        if isnan(sens_anal.(param_name).fitparams_upper.(varname).all.RMSE)
            sens_anal.(param_name).fitparams_upper = runGLMModelFit(theta_upper,fld_data,params,[param_name,'_upper'],conf);
        end
    end
    for param_i = 1:length(params)
        if isnan(sens_anal.(param_name).fitparams_lower.(varname).all.RMSE)
            num_nan = num_nan +1;
        end
        if isnan(sens_anal.(param_name).fitparams_upper.(varname).all.RMSE)
            num_nan = num_nan +1;
        end
    end
end

%Finish by running through with initial values so that lake.csv is in sim
%file
runGLMModelFit(theta_initial,fld_data,params,'init',conf);


%Determine sensitivity to each parameter for each data set using PRE
for param_i = 1:length(params)
    param_name = params{param_i}{1};
    for ii = 1:length(Data_Subsets)
        sens_anal.(param_name).(Data_Subsets{ii}).sens_lower = ...
            (sens_anal.(param_name).fitparams_lower.(varname).(Data_Subsets{ii}).PRE - ...
                fit_params_init.(varname).(Data_Subsets{ii}).PRE)/0.2;
        sens_anal.(param_name).(Data_Subsets{ii}).sens_upper = ...
            (sens_anal.(param_name).fitparams_upper.(varname).(Data_Subsets{ii}).PRE - ...
                fit_params_init.(varname).(Data_Subsets{ii}).PRE)/0.2;
    end
end


% ------------------------- PRINT RESULTS --------------------------------%

%Print table of lake metrics for each measure of model fit


% Open file
% ---------------------------- DEFINE THE FILE --------------------------- %

% Create .csv file.
fid = fopen('GLM_sensitivity.csv','w');

% Write header information

%Header of parmaeter names
header_line = '           ';
header_line = [header_line, ', Initial'];
for ii = 1:length(params)
    header_line = [header_line,', ',params{ii}{1}, ', '];
end
fprintf(fid,'%s\n', header_line);

%Header of uppers and lowers
header_line = 'Lake_Metric';
header_line = [header_line, ', Initial'];
for ii = 1:length(params)
    header_line = [header_line,', Lower', ', Upper'];
end
fprintf(fid,'%s\n', header_line);

% ---------------------------- STORE THE DATA ---------------------------- %

format_line = '%s, %f';
for ii = 1:length(params)
    format_line = [format_line,', %f , %f'];
end
format_line = [format_line,'\n'];

for fit_i = 1:length(Model_Fit)
    for data_i = 1:length(Data_Subsets)
        name_str = [Model_Fit{fit_i},'_',Data_Subsets{data_i}];
        data_line = fit_params_init.(varname).(Data_Subsets{data_i}).(Model_Fit{fit_i});
        for param_i = 1:length(params)
            data_line = [data_line sens_anal.(params{param_i}{1}).fitparams_lower.(varname).(Data_Subsets{data_i}).(Model_Fit{fit_i}), ...
                                   sens_anal.(params{param_i}{1}).fitparams_upper.(varname).(Data_Subsets{data_i}).(Model_Fit{fit_i})];
        end
        fprintf(fid,format_line,name_str,data_line);
    end
end
%Finish table with sensitivity analysis using slope of percentage relative
%error (PRE)
for data_i = 1:length(Data_Subsets)
    name_str = ['SensAnal_PRE_',Data_Subsets{data_i}];
    data_line = 0.0;
    for param_i = 1:length(params)
        data_line = [data_line sens_anal.(params{param_i}{1}).(Data_Subsets{data_i}).sens_lower, ...
                               sens_anal.(params{param_i}{1}).(Data_Subsets{data_i}).sens_upper];
    end
    fprintf(fid,format_line,name_str,data_line);
end

fclose(fid);                                      % Close the file.

save([lakename,'SensAnal.mat'],'sens_anal', '-v7.3')

% ------------------------------------------------------------------------ %
%If running remotely then remove ganymed-ssh2-build250 directory
if conf.config.remote_flag == 1
    rmdir('ganymed-ssh2-build250/','s');
end
cd ([base_dir,'MetaAnalysis/SensitivityAnalysis'])

% ------------------------------------------------------------------------ %
%Display time message
runTime     = toc;                %stop counting run time
hours       = runTime / 60;     %convert runTime from seconds to minute
timeMessage = ['Total Run Time = ', num2str(hours), ' minutes'];
disp(timeMessage)

    