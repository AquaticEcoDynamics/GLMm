function fit_params = rerunGLMinit(lakename,conf_file)
%function fit_params = rerunGLMinit(lakename,conf_file)
%
% Inputs:
%      lakenames:  name of lake
%      conf_file:  name of configuration file specifying paths and
%      parameters
%
% Outputs:
%
% Uses:
%      runGLMModelFit.m
%
% Written by L. Bruce 19 May 2014
%
% Re runs the intitial GLM for the lake, lakename and calculated model fit
% parameters

%CONFIGURATION INFORMATION------------------------------------------------%

conf = readGLMconfig(conf_file);

%Working directory to include Lake Name
conf.paths.working_dir = [conf.paths.base_dir,lakename,'/',conf.paths.working_dir];
base_dir = conf.paths.base_dir;
working_dir = conf.paths.working_dir;

%Create Output and Results folders
mkdir([working_dir,'/Output']);
mkdir([working_dir,'/Results']);
mkdir([working_dir,'/Results/Plots']);

%Move to working directory
paths = conf.paths;
cd (working_dir)


%If running remotely then add ganymed-ssh2-build250 directory'
if conf.config.remote_flag == 1
    
   ganymed_dir = 'C:\Louise\MATLAB\function\functionMisc\ganymed-ssh2-build250';
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
fit_params = runGLMModelFit(theta_initial,fld_data,params,'init',conf);

%Plot initial run and comparison against field - save results
plotLake(paths,'output_init.nc')
xlim([sim_start sim_stop]);
%x_tick = [sim_start:round((sim_start-sim_stop)/4):sim_stop];
x_tick_i = get(gca,'XTick');
x_tick = sim_start:(sim_stop - sim_start)/4:sim_stop;
set(gca,'XTick',x_tick,'XTickLabel',datestr(x_tick,'mm-yy'))

ylim([0 max_depth_plots])
fig_name = [paths.working_dir,'Results/Plots/Sim_Temp','.png'];
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'paperposition',[0  0 15 8]);
print(gcf,'-dpng',fig_name,'-opengl');


%INITIAL RUN COMPARISON PLOTS --------------------------------------------%

%Now plot scatter and time series comparisons

plotGLMModelFit(fit_params.(varname),[paths.working_dir,'Results/Plots/'])

%Plot temperature comparison figure for MLCP
temp_range = [0 35];
figure
set(gcf,'defaultAxesFontSize', 6) 
plot(fit_params.(varname).obs_data.all,fit_params.(varname).sim_data.all,'*');
xlim(temp_range);
ylim(temp_range);
hold on
plot(temp_range,temp_range,'k')
hold off
xlabel('Measured Temperature (^oC)','FontSize', 6)
ylabel('Simulated Temperature (^oC)','FontSize', 6)
title(lakename,'FontSize', 6)
fig_name = [base_dir,'MetaAnalysis/ModelFit/Plots/',lakename,'_alltemp.png'];
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'paperposition',[0.1  0.1 4 4]);
print(gcf,'-dpng',fig_name,'-opengl');

close all;

% ------------------------------------------------------------------------ %
%If running remotely then remove ganymed-ssh2-build250 directory
if conf.config.remote_flag == 1
    rmdir('ganymed-ssh2-build250/','s');
end
cd (base_dir)