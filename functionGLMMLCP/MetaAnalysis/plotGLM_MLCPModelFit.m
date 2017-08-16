%function MLCPmodel_fit = plotGLM_MLCPModelFit%(LakeNames,MetricNames)%,ParamNames)
%function MLCPmodel_fit = plotGLM_MLCPModelFit%(LakeNames,MetricNames)%,ParamNames)
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
% Written by L. Bruce 18 December 2013
%
% Plots results of calcGLMModelFit
% Takes the PRE's of each parameter to calculate a probability
% distribution

%Clean up
close all
clear all

%Base directory
base_dir = 'C:/Louise/GLM/GLM_v2.2.0_MLCP/';

%Folder to save plots to
dirName = [base_dir,'MetaAnalysis\ModelFit\Plots\'];

%CONFIGURATION INFORMATION------------------------------------------------%
conf = readGLMconfig([base_dir,'MetaAnalysis\SensitivityAnalysis\SA_config.nml']);
working_dir = conf.paths.working_dir;

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
LakeFont = {'normal','bold','normal','normal','bold', ...
            'normal','normal','normal','normal','normal', ...
            'normal','bold','bold','normal','bold', ...
            'normal','normal','normal','normal','bold', ...
            'normal','normal','normal','normal','bold', ...
            'bold','normal','normal','bold','normal', ...
            'normal','normal','normal','normal'};
LakeColour = {'black','red','black','black','red', ...
            'green','black','black','black','black', ...
            'green','red','red','black','red', ...
            'black','green','green','green','red', ...
            'black','black','black','black','red', ...
            'red','black','black','red','black', ...
            'black','black','black','black'};
          
MixedLakes = [{'Alexandrina'},'Emaiksoun','Muggelsee','Woods'];

LakeCharNames = {'Volume','Area','Depth','AonD','LonW','Inflow','ResTime', ...
                 'ShortWave','AirTemp','WindSpeed','Kw','Lat','LN','LN_strat','pcLN_lt1'};
LakeCharNamesExtra = {'MaxDepth','Max_Area','Crest Elevation','Start Date'};
LakeCharNamesFig = {'Volume (m^3)','Area (m^2)','Depth (m)','Area on Depth (m)', ...
                    'Length on Width','Inflow (m^3/day)','Residence Time (days)', ...
                    'Short Wave Radiation (W/m^2)','Air Temperature (^oC)', ...
                    'Wind Speed (m/s)','Extinction Coefficient','Latitude','Lake Number','LN_strat','%LN<1'};

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

if 1 %Invoke when need to collect new data
%-------------------------------------------------------------------------%
%Loop through lakes and:
%     1)Get Lake characterisation metrics for each lake
%     2) read PRE values for each parameter
 for lake_i = 1:numLakes
    disp(['Loading:  ',LakeNames{lake_i}]);
    %Working directory to include Lake Name
    conf.paths.working_dir = [conf.paths.base_dir,LakeNames{lake_i},'/',working_dir];
    %Load in resuts from GLM simulation such as depth, inflow etc
    Lake.(LakeNames{lake_i}) = getGLM_MLCPLakeMetrics(base_dir,LakeNames{lake_i});
    %Load in Lake Analyzer simulation metrics
    LA_Metrics = calcLaGLMmetrics([conf.paths.working_dir,'Output/output_init.nc'],conf);
     %Make various LN metrics
     if max(lake_i == StratLakes) == 0
         Lake.(LakeNames{lake_i}).LN = NaN; %Mixed lakes no LN
         Lake.(LakeNames{lake_i}).LN_strat = NaN;
         Lake.(LakeNames{lake_i}).pcLN_lt1 = NaN;
     else
         LN = LA_Metrics.LN;
         LN(LN<0) = 0;
         LN(isinf(LN)) = 0;
         LN_strat = LN(~LA_Metrics.mixed);
         LN_lt_1 = LN_strat(LN_strat<1);
         Lake.(LakeNames{lake_i}).LN = mean(LN);
         Lake.(LakeNames{lake_i}).LN_strat = mean(LN_strat);
         Lake.(LakeNames{lake_i}).pcLN_lt1 = length(LN_lt_1)/length(LN_strat)*100;
     end
    %Load in model fit data
    filename = [base_dir,LakeNames{lake_i},'/Results/GLM_fitparams_init.csv'];
    mf_data = importdata(filename,',',1);
    for metric_i = 1:numMetrics
        for mf_i = 1:numModelFitNames
            %For well mixed lakes like Emaliksoun all Model Fit metrics for
            %Schmidt Stability (St) are NaN therefore set to NaN
            if size(mf_data.data,1) >= metric_i
                Lake.(LakeNames{lake_i}).(MetricNames{metric_i}).(ModelFitNames{mf_i}) = ...
                    mf_data.data(metric_i,mf_i);
            else
               Lake.(LakeNames{lake_i}).(MetricNames{metric_i}).(ModelFitNames{mf_i}) = NaN;
            end
        end
    end
    %To compare distance from equator make latitude absolute
    Lake.(LakeNames{lake_i}).Lat = abs(Lake.(LakeNames{lake_i}).Lat);
 end
 
 save('MLCP_modelfit_v2.2.0_MLCP.mat','Lake')
 
 %Display Lake Metrics to check if all OK
 for lake_i = 1:numLakes
     disp(['Current Lake:  ',LakeNames{lake_i}]);
     Lake.(LakeNames{lake_i})
 end

elseif 0 %Just load previously saved data
    load('MLCP_modelfit_v2.2.0_MLCP.mat');
end

%-------------------------------------------------------------------------%

%Save Lake Names, and metrics to text file
%-------------------------------------------------------------------------%

fid = fopen('MLCP_LakeSummary_2.2.0_MLCP.csv','w');

%Headers
header = ['Abbreviation, Lake Name'];
for ii = 1:numLakeChars
    header = [header,',',LakeCharNames{ii}];
end
for ii = 1:length(LakeCharNamesExtra)
    header = [header,',',LakeCharNamesExtra{ii}];
end

fprintf(fid,'%s \n',header);

%Print lake metrics for each lake
for lake_i = 1:numLakes
    lc_string = [LakeInitials{lake_i},',',LakeNames{lake_i}];
    for ii = 1:numLakeChars
        lc_string = [lc_string,  ',',num2str(Lake.(LakeNames{lake_i}).(LakeCharNames{ii}))];
    end
    lc_string = [lc_string,',',num2str(Lake.(LakeNames{lake_i}).nml.H(end) - Lake.(LakeNames{lake_i}).nml.H(1)),',', ...
                               num2str(Lake.(LakeNames{lake_i}).nml.A(end)),',', ...
                               num2str(Lake.(LakeNames{lake_i}).nml.crest_elev),',', ...
                               datestr(Lake.(LakeNames{lake_i}).nml.start_date)];
    fprintf(fid,'%s \n',lc_string);
end
LakeCharNamesExtra = {'MaxDepth','Max_Area','Crest Elevation','Start Date'};

fclose(fid);


%-------------------------------------------------------------------------%
%Produce scatter plots for model fit against Lake Metrics

%-------------------------------------------------------------------------%
% Font size setting for figure
xlabel_fontsize = 10;
ylabel_fontsize = 10;
figure_fontsize = 10;

for mf_i = 1:numModelFitNames

%Loop through lake characterisation metrics and plot for each of 5
%stratification metrics combining on single plot

%For the following lake characteristics use a linear plot and for the rest
%use a log plot
linear_chars =  [5 8 9 10 11 12 15];
for lchar_i = 1:numLakeChars
    figure
    clear x
    clear y
    for lake_i = 1:numLakes
        x(lake_i) = Lake.(LakeNames{lake_i}).(LakeCharNames{lchar_i});
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
        if ~isempty(find(linear_chars == lchar_i))
            plot(x,y,'w*')
        else
            semilogx(x,y,'w*')
        end
        for lake_i = 1:numLakes
            text(x(lake_i),y(lake_i),LakeInitials{lake_i},'HorizontalAlignment','center','Fontsize',6)
        end
        %Plot significance of relationship
        if  ~isempty(find(linear_chars == lchar_i))
            plotr2andpval
        else
            plotr2andpval_semilog
        end
        if metric_i == 1
            ylabel(ModelFitNamesFig{mf_i},'Fontsize',8)
        end
        xlabel(LakeCharNames{lchar_i},'Fontsize',8)
        setXYLabelFontsize
        title(MetricNames{metric_i},'Fontsize',8)
        
        %Save significance and correlation
        MLCP_ModelFit.(MetricNames{metric_i}).(LakeCharNames{lchar_i}).(ModelFitNames{mf_i}).r2 = r2;
        MLCP_ModelFit.(MetricNames{metric_i}).(LakeCharNames{lchar_i}).(ModelFitNames{mf_i}).pval = pval;
        MLCP_ModelFit.(MetricNames{metric_i}).(LakeCharNames{lchar_i}).(ModelFitNames{mf_i}).x = x;
        MLCP_ModelFit.(MetricNames{metric_i}).(LakeCharNames{lchar_i}).(ModelFitNames{mf_i}).y = y;

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
    dirandnametosave = [dirName, 'MLCP_ModelFit_',LakeCharNames{lchar_i},'_',ModelFitNames{mf_i}];
    saveas(gcf, dirandnametosave,'png'); 
    close
end
end %model fit metrics

%Save significance and r2
 save('MLCP_modelfit_pval_v2.2.0_MLCP.mat','MLCP_ModelFit')

%-------------------------------------------------------------------------%
%Save Lake Names, numbers and model fit values to text file
%-------------------------------------------------------------------------%

fid = fopen('MLCP_ModelFit_v2.2.0_MLCP.csv','w');

%Headers
header = ['Number, Name'];
for ii = 1:length(MetricNames)
    for jj = 1:length(ModelFitNames)
        header = [header,',',MetricNamesFig{ii}];
    end
end

fprintf(fid,'%s \n',header);

header = ['Number, Name'];
for ii = 1:length(MetricNames)
    header = [header,',  RMSE, MEFF,     r,    PRE,   NMAE'];
end

fprintf(fid,'%s \n',header);

%Print model fit values for each lake
for lake_i = 1:numLakes
    mf_string = [num2str(lake_i),',',LakeNames{lake_i}];
    lake_mf = Lake.(LakeNames{lake_i});
    for ii = 1:length(MetricNames)
        for jj = 1:length(ModelFitNames)
            mf_string = [mf_string,  ',',num2str(lake_mf.(MetricNames{ii}).(ModelFitNames{jj}))];
        end
    end
    fprintf(fid,'%s \n',mf_string);
end

%Calculate mean and print for each model fit metric
mf_string = 'All  , Mean';
for ii = 1:length(MetricNames)
    for jj = 1:length(ModelFitNames)
        mf_array = [];
        for lake_i = 1:numLakes
            lake_mf = Lake.(LakeNames{lake_i});
            mf_array = [mf_array lake_mf.(MetricNames{ii}).(ModelFitNames{jj})];
        end
        mf_string = [mf_string,  ',',num2str(mean(mf_array(~isnan(mf_array))))];
    end
end
    fprintf(fid,'%s \n',mf_string);


fclose(fid);

%-------------------------------------------------------------------------%

%Save Lake Characteristics, metrics (model fit and lake), pval and r2 to text file
%-------------------------------------------------------------------------%

fid = fopen('MLCP_ModelFit_pval_v2.2.0_MLCP.csv','w');

%Headers
header = ['Number, Name'];
for ii = 1:length(MetricNames)
    for jj = 1:length(ModelFitNames)
        header = [header,',',MetricNamesFig{ii}];
    end
end

fprintf(fid,'%s \n',header);

header = ['Number, Name'];
for ii = 1:length(MetricNames)
    header = [header,',  RMSE, MEFF,   r2,   PRE, NMAE'];
end

fprintf(fid,'%s \n',header);

%Print model fit r2 for each lake characteristic
for lchar_i = 1:numLakeChars
    mf_string = [num2str(lchar_i),',',LakeCharNames{lchar_i}];
    for metric_i = 1:numMetrics
        for mf_i = 1:numModelFitNames
            mf_r2 = MLCP_ModelFit.(MetricNames{metric_i}).(LakeCharNames{lchar_i}).(ModelFitNames{mf_i}).r2;
            mf_string = [mf_string,  ',',num2str(mf_r2)];
        end
    end
    fprintf(fid,'%s \n',mf_string);
end

%Print model fit pval for each lake characteristic
for lchar_i = 1:numLakeChars
    mf_string = [num2str(lchar_i),',',LakeCharNames{lchar_i}];
    for metric_i = 1:numMetrics
        for mf_i = 1:numModelFitNames
            mf_pval = MLCP_ModelFit.(MetricNames{metric_i}).(LakeCharNames{lchar_i}).(ModelFitNames{mf_i}).pval;
            mf_string = [mf_string,  ',',num2str(mf_pval)];
        end
    end
    fprintf(fid,'%s \n',mf_string);
end
fclose(fid);