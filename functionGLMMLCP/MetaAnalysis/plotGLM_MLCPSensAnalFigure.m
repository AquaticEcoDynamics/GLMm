%plotGLM_MLCPSensAnalFigure

%Script to plot model fit correlations for the first MLCP paper.

%Created by L. Bruce 15 January 2016

%Clean up
clear all
close all

%Must have loaded MLCP_ModelFit and MLCP Lake
load('MLCP_SensAnal_CV.mat');
load('MLCP_SensAnal_v2.2.0.mat');

%Base directory
base_dir = 'C:/Louise/GLM/GLM_v2.2.0/';

%Folder to save plots to
dirName = [base_dir,'MetaAnalysis\SensitivityAnalysis\Plots\'];


%-------------------------------------------------------------------------%
%Set up metric and lake names
MetricNames = {'all','epi','hyp','thermoD','St'};
MetricNamesFig = {'All Temperature (^oC)','Epilimnion Temperature (^oC)','Hypolimnion Temperature (^oC)','Thermocline Depth (m)','Schmidt Number'};
FigNumber = {'9','10','11','12','13'};

LakeNames = [{'Alexandrina'},'Ammersee','Blelham','Bourget','Cannonsville',...
    'Como','Constance','ElGergal','Emaiksoun','Esthwaite','Feeagh', ...
    'Geneva01','Geneva03','GrosseDhunn','Harp','Iseo','Kinneret03', ...
    'Kinneret97','Mendota','MtBold','Muggelsee','NamCo','Oneida', ...
    'Pusiano','Rappbode', 'Rassnitzersee','Ravn','Rotorua', ...
    'Stechlin','Tarawera','Toolik','Windermere','Woods','Zurich'];
LakeInitials = [{'AL'},'AM','BL','BO','CA', ...
                 'CO','CN','EG','EM','ES','FE', ...
                 'G1','G3','GD','HA','IS','K3', ...
              'K7','ME','MB','MG','NM','ON', ...
              'PU','RP', 'RS','RV','RO', ...
              'ST','TA','TO','WI','WO','ZU'];
          
MixedLakes = [{'Alexandrina'},'Emaiksoun','Muggelsee','Woods'];
          
LakeCharNames = [{'Volume'},'Depth','Area','AonD','LonW','Inflow','ResTime','ShortWave','AirTemp','WindSpeed','Kw','Lat','LN'];
LakeCharNamesFig = [{'Volume (m^3)'},'Depth (m)','Area (m^2)','Area on Depth (m)','Length on Width','Inflow (m^3/s)','Residence Time (days)','Short Wave Radiation (W/m^2)','Air Temperature (^oC)','Wind Speed (m/s)','Extinction Coefficient','Latitude','Lake Number'];

ParamNames = [{'coef_mix_conv'},'coef_wind_stir','coef_mix_shear','coef_mix_turb', ...
               'coef_mix_KH','coef_mix_hyp','ce','ch','cd']; 
ParamNamesFig = [{'C_c'},'C_w','C_s','C_t','C_K_H','C_h_y_p','C_e','C_h','C_d'];
ParamNamesFull = {'Convective','Wind stir','Shear','Turbulence', ...
               'Kelvin Helmholtz','Hypolimnion','Evaporative','Sensible heat','Wind drag'}; 

numLakes = length(LakeNames);
numMetrics = length(MetricNames);
numLakeChars = length(LakeCharNames);
numParams = length(ParamNames);

%-------------------------------------------------------------------------%
%Set up figure settings

%Axes position
%3 across  down
axes_pos = zeros(6,4);
for ii = 1:3
    axes_pos(ii,:) = [0.08 + (ii-1)*0.33 0.2 0.25 0.75];
end
% Font size setting for figure
figure_fontsize = 8;
xlabel_fontsize = 8;
ylabel_fontsize = 8;

%X and Y limits
% Lake characteristics
xlims.Depth = [1e0 1e3];
xlims.Area = [1e5 1e10];
xlims.WindSpeed = [0 8];
xlims.ResTime = [1e0 1e6];
xlims.ShortWave = [0 300];
xlims.LonW = [0 40];
xlims.AirTemp = [-10 25];
xlims.Volume = [1e5 1e15];
xlims.Kw = [0 1.5];
xlims.Inflow = [1e2 1e10];
xlims.AonD = [1e2 1e10];

%Model fit
ylims.all = [0 0.4];
ylims.epi = [0 0.2];
ylims.hyp = [0 0.4];
ylims.thermoD = [0 1.5];
ylims.St = [0 1];
%-------------------------PLOT FIGURES------------------------------------%

%Set which parameters to plot---------------------------------------------%
fig_an = {'(a)','(b)','(c)'};
%All temperatures
lchars(1,:) = [3 2 9]; %Area, Depth, Temp 
params(1,:) = [5 7 8]; %kh, ce, ch
%Epi temperatures
lchars(2,:) = [9 2 7]; %Temp, Depth, ResTime 
params(2,:) = [3 5 9]; %shear,kh,cd
%Hypolimnion temperatures
lchars(3,:) = [9 8 6]; %Temp, shortwave, Inflow 
params(3,:) = [5 7 8]; %KH, ce, ch
%Thermocline depth
lchars(4,:) = [2 2 2]; %Depth, Depth, Depth
params(4,:) = [1 2 5]; %conv,wind,kh
%Schmidt number
lchars(5,:) = [10 10 10]; %wind LN
params(5,:) = [7 7 7]; %ce

for metric_i = 1:numMetrics
figure
for ii = 1:3
    lchar_i = lchars(metric_i,ii);
    param_i = params(metric_i,ii);
    if sum(find([5 8 9 10 11 12] == lchars(metric_i,ii))) %Linear plot
        lin_char = 1;
    else
        lin_char = 0;
    end
    %Plot figure including p value and r2
    x = MLCP_SA.(ParamNames{param_i}).(MetricNames{metric_i}).(LakeCharNames{lchar_i}).x;
    y = MLCP_SA.(ParamNames{param_i}).(MetricNames{metric_i}).(LakeCharNames{lchar_i}).y;
    %--------Plot figure-----------------------%
    %Set position of subplot 
    axes('position',axes_pos(ii,:))
    set(gcf,'defaultAxesFontSize', figure_fontsize)
    xlim(xlims.(LakeCharNames{lchar_i}))
    ylim(ylims.(MetricNames{metric_i}))
    x_lim=xlim; y_lim = ylim;
    if lin_char %Linear plot
        plot(x,y,'w*')
        plotr2andpval %Plot significance of relationship
    else %Log plot
        semilogx(x,y,'w*')
        plotr2andpval_semilog %Plot significance of relationship
    end

    %Plot lakes as initials
    for lake_i = 1:numLakes
        text(x(lake_i),y(lake_i),LakeInitials{lake_i}, ...
            'HorizontalAlignment','center','Fontsize',6,'Color','blue')
    end
    
    %Plot line of best fit
    hold on
    if lin_char %Linear plot
        plot(x_lim,ab(2) + ab(1)*(x_lim),'k')
    else
        plot(x_lim,ab(2) + ab(1)*log(x_lim),'k')
    end
    plot(x_lim,[0 0],'k--')
    hold off

    %Labels
    setXYLabelFontsize
    xlabel(LakeCharNamesFig{lchar_i})
    ylabel(['SI for ',ParamNamesFig{param_i}])
    %Figure annotation
   if lin_char
       text_x = x_lim(1)+0.9*(x_lim(2) - x_lim(1));
       text_y = y_lim(1)+0.05*(y_lim(2) - y_lim(1));
   else
       text_x = exp(log(x_lim(1))+0.9*(log(x_lim(2)) - log(x_lim(1))));
       text_y = y_lim(1)+0.05*(y_lim(2) - y_lim(1));
   end
    text(text_x,text_y,fig_an{ii},'FontSize',figure_fontsize)%,'Fontweight','bold')
    xlim(x_lim); ylim(y_lim);
end
if metric_i == 3
    set(gca,'XTick',[1e2 1e5 1e8])
end
%--------Save figure-----------------------%
dirandnametosave = [dirName, 'Figures\Figure',FigNumber{metric_i},'_SA_sig',MetricNames{metric_i}];
% Set up figure properties for printed version----------------------------%
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
xSize = 16;
ySize = 16/3; %http://cdn.elsevier.com/assets/pdf_file/0010/109963/Artwork.pdf
xLeft = (21-xSize)/2;
yTop = (30-ySize)/2;
set(gcf,'paperposition',[0 0 xSize ySize])
%-------------------------------------------------------------------------%
% Save figure as a png file
saveas(gcf, dirandnametosave,'png'); 
%Save figure as eps
print(gcf,'-depsc2',[dirandnametosave,'.eps'],'-painters'); 


end

close all

%Figure for colour SA colour map - 5 strat metrics
figure('Units','centimeters','Position',[15 5 14 14*3/2])
%Axes position
%2 across 3 down
axes_pos = zeros(6,4);
for ii = 1:2
    for jj = 1:3
       axes_pos((jj-1)*2+ii,:) = [0.05 + (ii-1)*0.5 0.07 + (3-jj)*0.3 0.35 0.25];
    end
end
%Loop through metrics
for strat_i = 1:numMetrics
    clear plot_matrix
    for lake_i = 1:numLakes
        for param_i = 1:numParams
            plot_matrix(lake_i,param_i) = Lake.(LakeNames{lake_i}).SA.Sens_Coef.(ParamNames{param_i}).(MetricNames{strat_i});
        end
    end
    %Reverse lake order so first is at top of plot then read down
    plot_matrix = flipud(plot_matrix);
    %Remove NaNs
    plot_matrix(isnan(plot_matrix)) = 0.0;
    %For plotting with pcolor need to extend matix by 1 value
    plot_matrix(:,numParams+1) = NaN;
    plot_matrix(numLakes+1,:) = NaN;
    axes('position',axes_pos(strat_i,:))
    set(gcf,'defaultAxesFontSize', 6)
    pcolor(plot_matrix)
    caxis([0 1])
    set(gca,'FontSize',6)
    set(gca,'XTick',1.5:1:numParams+0.5,'XTickLabel',[]);%ParamNamesFig);
    for param_i = 1:numParams
        text(0.2+param_i,-0.7,ParamNamesFig{param_i},'FontSize',8)
    end
    set(gca,'YTick',1.5:1:numLakes+0.5 ,'YTickLabel',fliplr(LakeInitials),'FontSize',6);%, 'TickLength',[0 0.025] );

    title(MetricNamesFig{strat_i},'FontSize',8)
end
axes('position',axes_pos(strat_i+1,:))
colorbar
caxis([0 1])
%Figure legend
for param_i = 1:numParams
    text(0,1-0.1*param_i,[ParamNamesFig{param_i},' = ',ParamNamesFull{param_i}],'FontSize',8)
end
set(gca,'Visible','off')
%--------Save figure-----------------------%
dirandnametosave = [dirName, 'SA_Tables\SA_AllFigure2'];
% Set up figure properties for printed version----------------------------%
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
xSize = 16;
ySize = 16/2*3; %http://cdn.elsevier.com/assets/pdf_file/0010/109963/Artwork.pdf
xLeft = (21-xSize)/2;
yTop = (30-ySize)/2;
set(gcf,'paperposition',[0 0 xSize ySize])
%-------------------------------------------------------------------------%
saveas(gcf, dirandnametosave,'png');
%Save figure as eps
print(gcf,'-depsc2',[dirandnametosave,'.eps'],'-painters'); 


%----------Plot Depth vs Area for MLCP----------------------------------%

figure
for lake_i = 1:numLakes
    x(lake_i) = Lake.(LakeNames{lake_i}).Area;
    y(lake_i) = Lake.(LakeNames{lake_i}).Depth;
end
    %--------Plot figure-----------------------%
    %Set position of subplot 
    set(gcf,'defaultAxesFontSize', figure_fontsize)
    xlim(xlims.Area)
    ylim(xlims.Depth)
    x_lim=xlim; y_lim = ylim;
    loglog(x,y,'b.')

    %Plot lakes as initials
    for lake_i = 1:numLakes
        text(x(lake_i),y(lake_i),LakeInitials{lake_i}, ...
            'HorizontalAlignment','left','Fontsize',8,'Color','blue')
    end
    
    %Labels
    setXYLabelFontsize
    xlabel('Area (m^2)')
    ylabel('Depth (m)')
    %Figure annotation
%--------Save figure-----------------------%
dirandnametosave = [dirName, 'Figures\Bathy_loglog'];
% Set up figure properties for printed version----------------------------%
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
xSize = 16;
ySize = 16; %http://cdn.elsevier.com/assets/pdf_file/0010/109963/Artwork.pdf
xLeft = (21-xSize)/2;
yTop = (30-ySize)/2;
set(gcf,'paperposition',[0 0 xSize ySize])
%-------------------------------------------------------------------------%
% Save figure as a png file
saveas(gcf, dirandnametosave,'png');    
