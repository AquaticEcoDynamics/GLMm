%plotGLM_MLCPModelFitFigure

%Script to plot model fit correlations for the first MLCP paper.

%Created by L. Bruce 14 January 2016

%Clean up
close all
clear all

%Must have loaded MLCP_ModelFit and MLCP Lake
load('MLCP_modelfit_pval_v2.2.0_MLCP.mat')
load('MLCP_modelfit_v2.2.0_MLCP.mat');

%Base directory
base_dir = 'C:/Louise/GLM/GLM_v2.2.0_MLCP/';

%Folder to save plots to
dirName = [base_dir,'MetaAnalysis\ModelFit\Plots\'];

%Create plotting directory
if ~exist(dirName,'dir'); mkdir(dirName); end
if ~exist([dirName '/Figures'],'dir'); mkdir([dirName '/Figures']); end


%-------------------------------------------------------------------------%
%Set up metric and lake names
MetricNames = {'all','epi','hyp','thermoD','St'};
MetricNamesFig = {'Full Profile Temperature','Epilimnion Temperature','Hypolimnion Temperature','Thermocline Depth','Schmidt Number'};

ModelFitNames = {'RMSE','NSE','r2','PRE','NMAE'};
ModelFitNamesFig = {'RMSE (o^C)','MEFF','r','PRE (%)','NMAE'};


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

LakeCharNames = {'Volume','Area','Depth','AonD','LonW','Inflow','ResTime', ...
                 'ShortWave','AirTemp','WindSpeed','Kw','Lat','LN','LN_strat','pcLN_lt1'};
LakeCharNamesFig = {'Volume (m^3)','Area (m^2)','Depth (m)','Area on Depth (m)', ...
                    'Length on Width','Inflow (m^3/day)','Residence Time (days)', ...
                    'Short Wave Radiation (W/m^2)','Air Temperature (^oC)', ...
                    'Wind Speed (m/s)','Extinction Coefficient','Latitude','Lake Number','Lake Number','%LN<1'};

numLakes = length(LakeNames);
numMetrics = length(MetricNames);
numLakeChars = length(LakeCharNames);
numModelFitNames = length(ModelFitNames);

%Create list of stratified (i.e. not mixed) lakes
StratLakes = [];
for lake_i = 1:numLakes
    if max(strcmp(LakeNames{lake_i},MixedLakes)) == 0
        StratLakes = [StratLakes lake_i];
    end
end

%Lake characteristics linear plot, all others log
linear_chars =  [5 8 9 10 11 12 15];

%-------------------------------------------------------------------------%
%Set up figure settings

%Axes position
%3 across 2 down
axes_pos = zeros(6,4);
for ii = 1:3
    axes_pos(ii+3,:) = [0.07 + (ii-1)*0.33 0.13 0.25 0.35];
    axes_pos(ii,:) = [0.07 + (ii-1)*0.33 0.6 0.25 0.35];
end

% Font size setting for figure
figure_fontsize = 8;
xlabel_fontsize = 8;
ylabel_fontsize = 8;

%X and Y limits
% Lake characteristics
xlims.Depth = [1e0 1e3];
xlims.WindSpeed = [0 6];
xlims.ResTime = [1e0 1e5];
xlims.ShortWave = [0 300];
xlims.LonW = [0 40];
xlims.AirTemp = [-10 30];
xlims.Volume = [1e5 1e15];
xlims.Kw = [0 1.5];
xlims.LN = [1e-1 1e5];
xlims.Inflow = [1e0 1e10];
xlims.LN_strat = [1e-2 1e10];
xlims.pcLN_lt1 = [0 100];

%Model fit
ylims.all.RMSE = [0 4];
ylims.all.PRE = [-20 20];
ylims.epi.RMSE = [0 5];
ylims.epi.PRE = [-20 20];
ylims.hyp.RMSE = [0 4];
ylims.hyp.NSE = [-4 1];
ylims.hyp.PRE = [-40 40];
ylims.thermoD.NSE = [0 1];
ylims.thermoD.PRE = [-100 100];
ylims.thermoD.NMAE = [0 1.5];
ylims.thermoD.r2 = [0.2 1];
ylims.St.NSE = [0 1];
ylims.St.PRE = [-100 100];
ylims.St.NMAE = [0 1];
ylims.St.r2 = [0.6 1];
%-------------------------PLOT FIGURES------------------------------------%

%Temperature figure-------------------------------------------------------%
%2 down by 3 across (2 x epi, 2 x hyp, 1 x epi) (4 x RMSE & PRE)
metrics = [1 3 3 1 2 3]; %2 x all, 2 x epi,  2 x hyp
lchars = [11 11 3 6 8 7]; %Kw, Depth, inflow/ inflow, SW, Res time 
mfs = [1 1 1 4 4 4]; %RMSE x 3, PRE x 3
fig_an = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};

figure
for ii = 1:6
    lchar_i = lchars(ii);
    if sum(find(linear_chars == lchars(ii))) %Linear plot
        lin_char = 1;
    else
        lin_char = 0;
    end
    %Plot figure including p value and r2
    x = MLCP_ModelFit.(MetricNames{metrics(ii)}).(LakeCharNames{lchars(ii)}).(ModelFitNames{mfs(ii)}).x;
    y = MLCP_ModelFit.(MetricNames{metrics(ii)}).(LakeCharNames{lchars(ii)}).(ModelFitNames{mfs(ii)}).y;
    %--------Plot figure-----------------------%
    %Set position of subplot (1:4)
    axes('position',axes_pos(ii,:))
    set(gcf,'defaultAxesFontSize', figure_fontsize)
    xlim(xlims.(LakeCharNames{lchars(ii)}))
    ylim(ylims.(MetricNames{metrics(ii)}).(ModelFitNames{mfs(ii)}))
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
    xlabel(LakeCharNamesFig{lchars(ii)})
    ylabel([MetricNamesFig{metrics(ii)},' ',ModelFitNamesFig{mfs(ii)}])
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
%--------Save figure-----------------------%
dirandnametosave = [dirName, 'Figures\Figure5_mf_sigTemp'];
% Set up figure properties for printed version----------------------------%
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
xSize = 16;
ySize = 16/3*2; %http://cdn.elsevier.com/assets/pdf_file/0010/109963/Artwork.pdf
xLeft = (21-xSize)/2;
yTop = (30-ySize)/2;
set(gcf,'paperposition',[0 0 xSize ySize])
%-------------------------------------------------------------------------%
% Save figure as a png file
saveas(gcf, dirandnametosave,'png'); 
%Save figure as eps
print(gcf,'-depsc2',[dirandnametosave,'.eps'],'-painters'); 

%-------------------------------------------------------------------------%
%Thermocline depth and St figure------------------------------------------%
%-------------------------------------------------------------------------%

%3 across ((2 x MEFF & PRE), 2 down (Thermo D & St)
metrics = [5 5 5 4 5 4]; %4 x St, 2 x Td)
lchars = [3 7 11 3 8 9]; %Kw, Depth, ResTime, Kw, SW, LN 
mfs = [5 5 5 5 4 4]; %NMAE x 5, PRE x 1

figure
for ii = 1:6
    lchar_i = lchars(ii);
    if sum(find(linear_chars == lchars(ii))) %Linear plot
        lin_char = 1;
    else
        lin_char = 0;
    end
    %Plot figure including p value and r2
    x = MLCP_ModelFit.(MetricNames{metrics(ii)}).(LakeCharNames{lchars(ii)}).(ModelFitNames{mfs(ii)}).x;
    y = MLCP_ModelFit.(MetricNames{metrics(ii)}).(LakeCharNames{lchars(ii)}).(ModelFitNames{mfs(ii)}).y;
    %--------Plot figure-----------------------%
    %Set position of subplot (1:4)
    axes('position',axes_pos(ii,:))
    set(gcf,'defaultAxesFontSize', figure_fontsize)
    xlim(xlims.(LakeCharNames{lchars(ii)}))
    ylim(ylims.(MetricNames{metrics(ii)}).(ModelFitNames{mfs(ii)}))
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
    xlabel(LakeCharNamesFig{lchars(ii)})
    ylabel([MetricNamesFig{metrics(ii)},' ',ModelFitNamesFig{mfs(ii)}])
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
%--------Save figure-----------------------%
dirandnametosave = [dirName, 'Figures\Figure6_mf_sigTdSt'];
% Set up figure properties for printed version----------------------------%
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
xSize = 16;
ySize = 16/3*2; %http://cdn.elsevier.com/assets/pdf_file/0010/109963/Artwork.pdf
xLeft = (21-xSize)/2;
yTop = (30-ySize)/2;
set(gcf,'paperposition',[0 0 xSize ySize])
%-------------------------------------------------------------------------%
% Save figure as a png file
saveas(gcf, dirandnametosave,'png'); 
%Save figure as eps
print(gcf,'-depsc2',[dirandnametosave,'.eps'],'-painters'); 

%-------------------------------------------------------------------------%
%Lake Number figure------------------------------------------%
%-------------------------------------------------------------------------%

%2 across 4 down
axes_pos = zeros(8,4);
for ii = 1:2:8
    axes_pos(ii,:) = [0.07 0.1+(7-ii)*0.115 0.4 0.185];
    axes_pos(ii+1,:) = [0.07+0.5 0.1+(7-ii)*0.115 0.4 0.185];
end

%2 across LN & %LN<1, 4 down epi, hyp, Td, St 
metrics = [2 2 3 3 4 4 5 5 6 6]; 
lchars = [14 15 14 15 14 15 14 15]; %LN x 3, %LN<1 x 3
mfs = [1 1 2 2 3 3 3 3]; 

figure
for ii = 1:8
    lchar_i = lchars(ii);
    if sum(find(linear_chars == lchars(ii))) %Linear plot
        lin_char = 1;
    else
        lin_char = 0;
    end
    %Plot figure including p value and r2
    x = MLCP_ModelFit.(MetricNames{metrics(ii)}).(LakeCharNames{lchars(ii)}).(ModelFitNames{mfs(ii)}).x;
    y = MLCP_ModelFit.(MetricNames{metrics(ii)}).(LakeCharNames{lchars(ii)}).(ModelFitNames{mfs(ii)}).y;
    %--------Plot figure-----------------------%
    %Set position of subplot (1:4)
    axes('position',axes_pos(ii,:))
    set(gcf,'defaultAxesFontSize', figure_fontsize)
    xlim(xlims.(LakeCharNames{lchars(ii)}))
    ylim(ylims.(MetricNames{metrics(ii)}).(ModelFitNames{mfs(ii)}))
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
    xlabel(LakeCharNamesFig{lchars(ii)})
    ylabel([MetricNamesFig{metrics(ii)},' ',ModelFitNamesFig{mfs(ii)}])
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
%--------Save figure-----------------------%
dirandnametosave = [dirName, 'Figures\Figure7_mf_sigLN'];
% Set up figure properties for printed version----------------------------%
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
xSize = 16;
ySize = 16/2*4; %http://cdn.elsevier.com/assets/pdf_file/0010/109963/Artwork.pdf
xLeft = (21-xSize)/2;
yTop = (30-ySize)/2;
set(gcf,'paperposition',[0 0 xSize ySize])
%-------------------------------------------------------------------------%
% Save figure as a png file
saveas(gcf, dirandnametosave,'png'); 
%Save figure as eps
print(gcf,'-depsc2',[dirandnametosave,'.eps'],'-painters'); 




