function plotGLMModelFit(fitparams,plot_dir)
% Inputs:
%       varname    : variable name
%       obs_data   : MATLAB data structure containing observed data
%       sim_data   : MATLAB data structure containing simulated data
%       plot_dir   : directory to save plots to
%
% Outputs
%
% Uses:
%      plot.m
%
% Written by L. Bruce 15 August 2013
% Takes GLM simulated output and compares against field data for
% variable "varname".
% Plots the following lake metrics for observed against simulated data:
%    All data from observed profiles (all)
%    Volumetric average epilimnion temperature (epi)
%    Volumetric average hypolimnion temperature (hyp)
%    Schmidt Number (St)

%Get observed and simulated data
obs_data = fitparams.obs_data;
sim_data = fitparams.sim_data;


%First determine x date lables for time series plots
%Determine x tick marks
xmin = min(fitparams.time);
xmax = max(fitparams.time);
x_lim = [xmin xmax];
x_tick = [xmin:round((xmax-xmin)/4):xmax];
x_tick_date = datestr(x_tick,'dd-mmm-yy');
max_depth = max(fitparams.depth);

%Temperature range (should change)
temp_range = [0 30];


%All data
figure
plot(obs_data.all,sim_data.all,'*');
xlim(temp_range);
ylim(temp_range);
hold on
plot(temp_range,temp_range,'k')
hold off
xlabel('Measured Temperature (^oC)')
ylabel('Simulated Temperature (^oC)')
title('All Measured Data')
fig_name = [plot_dir,'All_Scatter','.png'];
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'paperposition',[0.1  0.1 15 8]);
print(gcf,'-dpng',fig_name,'-opengl');


%Epilimnion and Hypolimnion Data
figure
plot(obs_data.epi,sim_data.epi,'r*');
hold on
plot(obs_data.hyp,sim_data.hyp,'b*');
xlim(temp_range);
ylim(temp_range);
hold on
plot(temp_range,temp_range,'k')
hold off
xlabel('Measured Temperature (^oC)')
ylabel('Simulated Temperature (^oC)')
title('Epilimnion: red,  Hypolimnion: blue')
fig_name = [plot_dir,'EpiHyp_Scatter','.png'];
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'paperposition',[0.1  0.1 15 8]);
print(gcf,'-dpng',fig_name,'-opengl');
%Time series
figure
plot(fitparams.time,obs_data.hyp,'b*', ...
     fitparams.time,sim_data.hyp,'b', ...
     fitparams.time,obs_data.epi,'r*', ...
     fitparams.time,sim_data.epi,'r');
ylim(temp_range);
ylabel('Epilimnion and Hypolimnion Temperature (^oC)')
xlabel('Date')
xlim(x_lim)
set(gca,'XTick',x_tick,'XTickLabel',x_tick_date)
legend('Observed Hyp','Simulated Hyp','Observed Epi','Simulated Epi')
fig_name = [plot_dir,'EpiHyp_TS','.png'];
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'paperposition',[0.1  0.1 15 8]);
print(gcf,'-dpng',fig_name,'-opengl');

%Thermocline depth
figure
plot(obs_data.thermoD,sim_data.thermoD,'*');
xlim([0 max_depth]);
ylim([0 max_depth]);
hold on
plot([0 max_depth],[0,max_depth],'k')
hold off
xlabel('Measured Depth (m)')
ylabel('Simulated Depth (m)')
title('Thermocline depth')
fig_name = [plot_dir,'ThermoD_Scatter','.png'];
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'paperposition',[0.1  0.1 15 8]);
print(gcf,'-dpng',fig_name,'-opengl');
%Time series
figure
plot(fitparams.time,obs_data.thermoD,'*', ...
     fitparams.time,sim_data.thermoD,'b');
axis ij
ylim([0 max_depth]);
ylabel('Thermocline depth (m)')
xlabel('Date')
xlim(x_lim)
set(gca,'XTick',x_tick,'XTickLabel',x_tick_date)
legend('Observed','Simulated')
fig_name = [plot_dir,'ThermoD_TS','.png'];
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'paperposition',[0.1  0.1 15 8]);
print(gcf,'-dpng',fig_name,'-opengl');

%Schmidt Number
if min(isnan(sim_data.St)) == 0
figure
max_St = ceil(max(max(obs_data.St),max(sim_data.St))/1000)*1000;
plot(obs_data.St,sim_data.St,'*');
xlim([0 max_St]);
ylim([0 max_St]);
hold on
plot([0 max_St],[0,max_St],'k')
hold off
xlabel('Measured St')
ylabel('Simulated St')
title('Schmidt Stability')
fig_name = [plot_dir,'Schmidt_Scatter','.png'];
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'paperposition',[0.1  0.1 15 8]);
print(gcf,'-dpng',fig_name,'-opengl');
%Time series
figure
plot(fitparams.time,obs_data.St,'*', ...
     fitparams.time,sim_data.St,'b');
ylim([0 max_St]);
ylabel('Schmidt Stability')
xlabel('Date')
xlim(x_lim)
set(gca,'XTick',x_tick,'XTickLabel',x_tick_date)
legend('Observed','Simulated')
fig_name = [plot_dir,'Schmidt_TS','.png'];
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'paperposition',[0.1  0.1 15 8]);
print(gcf,'-dpng',fig_name,'-opengl');
end %if St not all NaNs


%%%%--------------------MEGA PLOT---------------------------------%%%%%%%

figure
%nml=getGLMnml;
x_tick = [xmin:round((xmax-xmin)/2):xmax];
x_tick_date = datestr(x_tick,'dd-mmm-yy');

%Epi and Hyp Temperature
subplot(3,2,1)
plot(obs_data.epi,sim_data.epi,'r*');
hold on
plot(obs_data.hyp,sim_data.hyp,'b*');
xlim(temp_range);
ylim(temp_range);
hold on
plot(temp_range,temp_range,'k')
hold off
xlabel('Measured Temp (^oC)')
ylabel('Simulated Temp (^oC)')
title('Epilimnion: red,  Hypolimnion: blue')
%Time series
subplot(3,2,2)
plot(fitparams.time,obs_data.hyp,'b*', ...
     fitparams.time,sim_data.hyp,'b', ...
     fitparams.time,obs_data.epi,'r*', ...
     fitparams.time,sim_data.epi,'r');
ylim(temp_range);
ylabel('Temperature (^oC)')
xlabel('Date')
xlim(x_lim)
set(gca,'XTick',x_tick,'XTickLabel',x_tick_date)

%Thermocline depth
subplot(3,2,3)
plot(obs_data.thermoD,sim_data.thermoD,'*');
xlim([0 max_depth]);
ylim([0 max_depth]);
hold on
plot([0 max_depth],[0,max_depth],'k')
hold off
xlabel('Measured Depth (m)')
ylabel('Simulated Depth (m)')
title('Thermocline depth')
%Time series
subplot(3,2,4)
plot(fitparams.time,obs_data.thermoD,'*', ...
     fitparams.time,sim_data.thermoD,'b');
axis ij
ylim([0 max_depth]);
ylabel('Thermocline depth (m)')
xlabel('Date')
xlim(x_lim)
set(gca,'XTick',x_tick,'XTickLabel',x_tick_date)


%Schmidt Number
if min(isnan(sim_data.St)) == 0
subplot(3,2,5)
max_St = ceil(max(max(obs_data.St),max(sim_data.St))/1000)*1000;
plot(obs_data.St,sim_data.St,'*');
xlim([0 max_St]);
ylim([0 max_St]);
hold on
plot([0 max_St],[0,max_St],'k')
hold off
xlabel('Measured St')
ylabel('Simulated St')
title('Schmidt Stability')
%Time series
subplot(3,2,6)
plot(fitparams.time,obs_data.St,'*', ...
     fitparams.time,sim_data.St,'b');
ylim([0 max_St]);
ylabel('Schmidt Stability')
xlabel('Date')
xlim(x_lim)
set(gca,'XTick',x_tick,'XTickLabel',x_tick_date)
end %if St not all NaNs

%Print Mega Plot
fig_name = [plot_dir,'ModelFit','.png'];
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperUnits', 'centimeters');
print(gcf,'-dpng',fig_name,'-opengl');

