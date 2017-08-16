function MLCPsens = runGLM_MLCPSensitivityAnalysis%(LakeNames)%,MetricNames)%,ParamNames)
%function plotGLMSensitivityAnalysis(LakeNames)
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
% Written by L. Bruce 19 August 2013
%
% Plots results of runGLMSensitivity.
% Takes the PRE's of each parameter to calculate a probability distribution

%Start by cleaning up
clear all
close all

%Add path
addpath(genpath('C:\Louise\MATLAB\function\functionGLMMLCP'))

%GLM_MLCP Sensitivity analysis configuration file
config_file = 'C:\Louise\GLM\GLM_v2.2.0_MLCP\MetaAnalysis\SensitivityAnalysis\SA_config.nml';

ParamNames = [{'coef_mix_conv'},{'coef_wind_stir'},{'coef_mix_shear'}, ...      
    {'coef_mix_turb'},{'coef_mix_KH'},{'coef_mix_hyp'}, ...
    {'ce'},{'ch'},{'cd'}];

MetricNames = {'all','epi','hyp','thermoD','St'};

LakeNames = [{'Alexandrina'},'Ammersee','Blelham','Bourget','Cannonsville',...
    'Como','Constance','ElGergal','Emaiksoun','Esthwaite','Feeagh', ...
    'Geneva01','Geneva03','GrosseDhunn','Harp','Iseo','Kinneret03', ...
    'Kinneret97','Mendota','MtBold','Muggelsee','NamCo','Oneida', ...
    'Pusiano','Rappbode', 'Rassnitzersee','Ravn','Rotorua', ...
    'Stechlin','Tarawera','Toolik','Windermere','Woods','Zurich'];

numLakes = length(LakeNames);
numMetrics = length(MetricNames);
numParams = length(ParamNames);

%Loop through lakes and run GLM Sensitivity Analysis for each lake
for lake_i = 18:numLakes
    lakename = LakeNames{lake_i};
    disp(['Current Lake:  ',lakename]);
    MLCPsens.(lakename) = runGLMSensitivity(lakename,config_file);
    close all
end

%save('MLCPsens_2.2.0.mat','MLCPsens','-v7.3')

% Remove unwanted output files
if 0
    conf = readGLMconfig(conf_file);
    for lake_i = 1:numLakes
        working_dir = [conf.paths.base_dir,lakename,'/',conf.paths.working_dir];
        cd ([working_dir,'Output'])
        delete('*.nc')
    end
end