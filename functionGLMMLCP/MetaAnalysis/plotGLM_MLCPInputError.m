%function MLCPmodel_fit = plotGLM_MLCPInputError%(LakeNames,MetricNames)%,ParamNames)
%function MLCPmodel_fit = plotGLM_MLCPInputError%(LakeNames,MetricNames)%,ParamNames)
%
% Inputs:
%      LakeNames:  Array of strings corresponding to each of the lakes used
%      in the MLCP sensitivity analysis
%      MetricNames:  Names of Lake Metrics to plot sensitivity to
%      ParamNames: List of parameters used in the sensitivity analysis
%
% Outputs:
%      MLCPsens_anal: MATLAB structure containing percent average error
%      (PRE) for each of the lakes listed in LakeNames for each of the
%      parameters listed in ParamNames.
%
% Uses:
%      plot.m
%
% Written by L. Bruce 29 January 2016
%
% Plots results of calcGLMModelFit
% Takes the Model Fit calculations and correlates with input error

%Clean up
close all
clear all

%Base directory
base_dir = 'C:/Louise/GLM/GLM_v2.2.0_MLCP/';

%Folder to save plots to
dirName = [base_dir,'MetaAnalysis\ModelFit\Plots\InputError\'];

%Create directory for plots
if ~exist(dirName,'dir') mkdir(dirName); end

%Input data error file
input_file = 'C:\Louise\GLM\Puplications\Multi-Lake_part1\Paper\Tables\MLCP_InputData.xlsx';

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

InputNames = {'morph','distmet','freqmet','flow','Kw','obs','mean'};
InputNamesXls = {'Morphometry','Dist. Met','Freq. Met','Flow','Kw','Obs Data','Mean'};
InputNamesFig = {'Morphometry','Meteorology Distance','Meteorology Frequency','Flows','Light Extinction','Observed Data','All Input Error'};

numLakes = length(LakeNames);
numMetrics = length(MetricNames);
numInputNames = length(InputNames);
numModelFitNames = length(ModelFitNames);

%Create list of stratified (i.e. not mixed) lakes
StratLakes = [];
for lake_i = 1:numLakes
    if max(strcmp(LakeNames{lake_i},MixedLakes)) == 0
        StratLakes = [StratLakes lake_i];
    end
end

%Load previously saved model fit and lake characteristic data
load([base_dir,'MetaAnalysis\ModelFit\MLCP_modelfit_v2.2.0_MLCP.mat']);

%Read in input error metrics
[Input_Error,InputDataNames,~] = xlsread(input_file,'Input Error','AG2:AM36');
for ii = 1:numInputNames
    input_idx = find(strcmp(InputDataNames,InputNamesXls{ii}));
    for lake_i = 1:numLakes
        Lake.(LakeNames{lake_i}).(InputNames{ii}) = Input_Error(lake_i,input_idx);
    end
end

%-------------------------------------------------------------------------%
%Produce scatter plots for model fit against input error

%-------------------------------------------------------------------------%
% Font size setting for figure
xlabel_fontsize = 10;
ylabel_fontsize = 10;
figure_fontsize = 10;

for mf_i = 1:numModelFitNames%[1 2 4]

%Loop through input error metrics and plot for each of 5
%stratification metrics combining on single plot

%Plot for each model input error metric
for inerr_i = 1:numInputNames
    figure
    clear x
    for lake_i = 1:numLakes
        x(lake_i) = Lake.(LakeNames{lake_i}).(InputNames{inerr_i});
    end
    %Loop through 5 lake thermodynamic metrics
    for metric_i = 1:numMetrics
        clear y
        if metric_i <= 2 %Include all lakes
            for lake_i = 1:numLakes
                y(lake_i) = Lake.(LakeNames{lake_i}).(MetricNames{metric_i}).(ModelFitNames{mf_i});
            end
        else %Only include stratified lakes
            y = NaN*ones(numLakes,1)';
            for lake_i = StratLakes
                 y(lake_i) = Lake.(LakeNames{lake_i}).(MetricNames{metric_i}).(ModelFitNames{mf_i});
            end
        end
            
        axes('position',[0.06 + (metric_i-1)*0.19 0.20 0.2 0.65])
        plot(x,y,'w*')
        for lake_i = 1:numLakes
            text(x(lake_i),y(lake_i),LakeInitials{lake_i},'HorizontalAlignment','center','Fontsize',6)
        end
        %Plot significance of relationship
        plotr2andpval
        if metric_i == 1
            ylabel(ModelFitNamesFig{mf_i},'Fontsize',8)
        end
        xlabel(InputNames{inerr_i},'Fontsize',8)
        setXYLabelFontsize
        title(MetricNames{metric_i},'Fontsize',8)
        
        %Save significance and correlation
        MLCP_ModelErr.(MetricNames{metric_i}).(InputNames{inerr_i}).(ModelFitNames{mf_i}).r2 = r2;
        MLCP_ModelErr.(MetricNames{metric_i}).(InputNames{inerr_i}).(ModelFitNames{mf_i}).pval = pval;
        MLCP_ModelErr.(MetricNames{metric_i}).(InputNames{inerr_i}).(ModelFitNames{mf_i}).x = x;
        MLCP_ModelErr.(MetricNames{metric_i}).(InputNames{inerr_i}).(ModelFitNames{mf_i}).y = y;

    end %Loop through lake thermodynamic properties
    
   
    %-------------------------------------------------------------------------%
    % Set up figure properties for printed version
    %--% Paper Size
                   set(gcf, 'PaperPositionMode', 'manual');
                   set(gcf, 'PaperUnits', 'centimeters');
                   xSize = 16;
                   ySize = 6; %http://cdn.elsevier.com/assets/pdf_file/0010/109963/Artwork.pdf
                   xLeft = (21-xSize)/2;
                   yTop = (30-ySize)/2;
                   set(gcf,'paperposition',[0 0 xSize ySize])
    %-------------------------------------------------------------------------%
    dirandnametosave = [dirName, 'MLCP_ModelErr_',InputNames{inerr_i},'_',ModelFitNames{mf_i}];
    saveas(gcf, dirandnametosave,'png'); 
    close
end
end %model fit metrics

%Save significance and r2
 save('MLCP_modelerr_pval_v2.2.0_MLCP.mat','MLCP_ModelErr')

%----------------------SAVE SIG PVALS TO FILE-----------------------------%
%Save input error, metrics (model fit and lake), pval and r2 to text file
%-------------------------------------------------------------------------%

    fid = fopen('MLCP_ModelErr_pval_v2.2.0_MLCP.csv','w');

    %Headers
    header = ['Number, Name'];
    for ii = 1:length(MetricNames)
        for jj = 1:numModelFitNames%length([1 2 4])
            header = [header,',',MetricNamesFig{ii}];
        end
    end

    fprintf(fid,'%s \n',header);

    header = ['Number, Name'];
    for ii = 1:length(MetricNames)
        header = [header,',RMSE,   NSE,   r2,    PRE,   NMAE'];
    end

    fprintf(fid,'%s \n',header);

    %Print model fit r2 for each lake characteristic
    for inerr_i = 1:numInputNames
        mf_string = [num2str(inerr_i),',',InputNames{inerr_i}];
        for metric_i = 1:length(MetricNames)
            for mf_i = 1:numModelFitNames%[1 2 4]
                mf_r2 = MLCP_ModelErr.(MetricNames{metric_i}).(InputNames{inerr_i}).(ModelFitNames{mf_i}).r2;
                mf_string = [mf_string,  ',',num2str(mf_r2)];
            end
        end
        fprintf(fid,'%s \n',mf_string);
    end

    %Print model fit pval for each lake characteristic
    for inerr_i = 1:numInputNames
        mf_string = [num2str(inerr_i),',',InputNames{inerr_i}];
        for metric_i = 1:length(MetricNames)
            for mf_i = 1:numModelFitNames%[1 2 4]
                mf_pval = MLCP_ModelErr.(MetricNames{metric_i}).(InputNames{inerr_i}).(ModelFitNames{mf_i}).pval;
                mf_string = [mf_string,  ',',num2str(mf_pval)];
            end
        end
        fprintf(fid,'%s \n',mf_string);
    end
    fclose(fid);

%-------------------------------------------------------------------------%
%-----------------CREATE SIGNIFICANT PLOTS FOR FIGURE---------------------%
%-------------------------------------------------------------------------%

%Set up figure settings

%Axes position
%3 across 2 down
axes_pos = zeros(6,4);
for ii = 1:3
    axes_pos(ii,:) = [0.07 + (ii-1)*0.33 0.17 0.25 0.7];
end

% Font size setting for figure
figure_fontsize = 8;
xlabel_fontsize = 8;
ylabel_fontsize = 8;

%X and Y limits
% Lake characteristics
xlims.InputErr = [0 3];

%Model fit
ylims.all.RMSE = [0 4];
ylims.all.PRE = [-20 20];
ylims.all.r2 = [0.75 1];
ylims.all.NMAE = [0 0.4];
ylims.epi.RMSE = [0 4];
ylims.epi.PRE = [-25 25];
ylims.epi.r2 = [0.75 1];
ylims.hyp.RMSE = [0 4];
ylims.hyp.PRE = [-40 30];
ylims.hyp.NMAE = [0 0.4];
ylims.thermoD.NSE = [-0.5 1];
ylims.thermoD.PRE = [-100 100];
ylims.thermoD.NMAE = [0 1.2];
ylims.St.NSE = [-0.5 1];
ylims.St.PRE = [-100 50];
ylims.St.NMAE = [0 1.5];
%-------------------------PLOT FIGURES------------------------------------%

%Create plotting directory
if ~exist([dirName, 'Figures'],'dir') mkdir([dirName, 'Figures']); end

%PRE figure-------------------------------------------------------%
%3 across, 1 down(3 x all, 3 St)
metrics = [1 1 3]; %2 x all, 1Td, 1Hyp 2 x St
inerr = [2 7 3]; %mean, mean, dist met 
mfs = [4 4 4]; %PRE
fig_an = {'(a)','(b)','(c)','(d)','(e)','(f)'};

figure
for ii = 1:3
    %Plot figure including p value and r2
    x = MLCP_ModelErr.(MetricNames{metrics(ii)}).(InputNames{inerr(ii)}).(ModelFitNames{mfs(ii)}).x;
    y = MLCP_ModelErr.(MetricNames{metrics(ii)}).(InputNames{inerr(ii)}).(ModelFitNames{mfs(ii)}).y;
    %--------Plot figure-----------------------%
    %Set position of subplot (1:4)
    axes('position',axes_pos(ii,:))
    set(gcf,'defaultAxesFontSize', figure_fontsize)
    xlim(xlims.InputErr)
    ylim(ylims.(MetricNames{metrics(ii)}).(ModelFitNames{mfs(ii)}))
    x_lim=xlim; y_lim = ylim;
    plot(x,y,'w*')
    plotr2andpval %Plot significance of relationship

    %Plot lakes as initials
    for lake_i = 1:numLakes
        text(x(lake_i),y(lake_i),LakeInitials{lake_i}, ...
            'HorizontalAlignment','center','Fontsize',6,'Color','blue')
    end
    
    %Plot line of best fit
    hold on
    plot(x_lim,ab(2) + ab(1)*x_lim,'k')
    plot(x_lim,[0 0],'k--')
    hold off

    %Labels
    setXYLabelFontsize
    xlabel(InputNamesFig{inerr(ii)})
    ylabel([MetricNamesFig{metrics(ii)},' ',ModelFitNamesFig{mfs(ii)}])
    %Xtick lables
    set(gca,'XTick',[0 1 2 3])
    set(gca,'XTickLabel',{'ideal','low','medium','high'})
    %Figure annotation
    text_x = x_lim(1)+0.9*(x_lim(2) - x_lim(1));
    text_y = y_lim(1)+0.05*(y_lim(2) - y_lim(1));
    text(text_x,text_y,fig_an{ii},'FontSize',figure_fontsize)%,'Fontweight','bold')
    xlim(x_lim); ylim(y_lim);
end
%--------Save figure-----------------------%
dirandnametosave = [dirName, 'Figures\Fig4abcinput_error_PRE'];
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
saveas(gcf, dirandnametosave,'png')
%Save figure as eps
print(gcf,'-depsc2',[dirandnametosave,'.eps'],'-painters'); 

%-------------------------------------------------------------------------%
%-----------------CREATE 2 MAIN SIGNIFICANT PLOTS FOR FIGURE--------------%
%-------------------------------------------------------------------------%

%Set up figure settings for NMAE

%Axes position
%3 across 1 down
axes_pos = zeros(6,4);
for ii = 1:3
    axes_pos(ii,:) = [0.07 + (ii-1)*0.33 0.17 0.25 0.7];
end

%-------------------------PLOT FIGURES------------------------------------%

%Temperature figure-------------------------------------------------------%
%2 across, 1 down(1 x all, 1 Td)
metrics = [1 2 3]; %1 x all, 2 x Td
inerr = [5 7 6]; %Freq met, obs
mfs = [3 3 5]; %NMAE
fig_an = {'(d)','(e)','(f)'};

figure
for ii = 1:3
    %Plot figure including p value and r2
    x = MLCP_ModelErr.(MetricNames{metrics(ii)}).(InputNames{inerr(ii)}).(ModelFitNames{mfs(ii)}).x;
    y = MLCP_ModelErr.(MetricNames{metrics(ii)}).(InputNames{inerr(ii)}).(ModelFitNames{mfs(ii)}).y;
    %--------Plot figure-----------------------%
    %Set position of subplot (1:4)
    axes('position',axes_pos(ii,:))
    set(gcf,'defaultAxesFontSize', figure_fontsize)
    xlim(xlims.InputErr)
    ylim(ylims.(MetricNames{metrics(ii)}).(ModelFitNames{mfs(ii)}))
    x_lim=xlim; y_lim = ylim;
    plot(x,y,'w*')
    plotr2andpval %Plot significance of relationship

    %Plot lakes as initials
    for lake_i = 1:numLakes
        text(x(lake_i),y(lake_i),LakeInitials{lake_i}, ...
            'HorizontalAlignment','center','Fontsize',6,'Color','blue')
    end
    
    %Plot line of best fit
    hold on
    plot(xlims.InputErr,ab(2) + ab(1)*xlims.InputErr,'k')
    hold off

    %Labels
    setXYLabelFontsize
    xlabel(InputNamesFig{inerr(ii)})
    ylabel([MetricNamesFig{metrics(ii)},' ',ModelFitNamesFig{mfs(ii)}])
    %Xtick lables
    set(gca,'XTick',[0 1 2 3])
    set(gca,'XTickLabel',{'ideal','low','medium','high'})
    %Figure annotation
    text_x = x_lim(1)+0.9*(x_lim(2) - x_lim(1));
    text_y = y_lim(1)+0.05*(y_lim(2) - y_lim(1));
    text(text_x,text_y,fig_an{ii},'FontSize',figure_fontsize)%,'Fontweight','bold')
    xlim(x_lim); ylim(y_lim);
end
%--------Save figure-----------------------%
dirandnametosave = [dirName, 'Figures\Fig4definput_error_NMAE'];
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
