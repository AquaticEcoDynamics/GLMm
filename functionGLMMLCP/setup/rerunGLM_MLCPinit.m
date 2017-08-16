%function rerunGLM_MLCPinit%(LakeNames)%,MetricNames)%,ParamNames)
%function rerunGLM_MLCPinit(LakeNames)
%
% Inputs:
%      LakeNames:  Array of strings corresponding to each of the lakes used
%      in the MLCP sensitivity analysis
%
% Outputs:
%
% Uses:
%      rerunGLMinit.m
%
% Written by L. Bruce 19 May 2014
%
% Runs a series of GLM initial runs (standard parameter set) for the lakes in
% the Multi Lake Comparison Project (MLCP)

%Add path to GLM MLCP scripts
addpath('C:\Louise\MATLAB\function\functionGLMMLCP\')

%Start by cleaning up
clear all
close all

base_dir = pwd;

%GLM_MLCP configuration file
config_file = 'C:\Louise\GLM\GLM_v2.2.0_MLCP\MetaAnalysis\ModelFit\rerun_config.nml';

LakeNames = [{'Alexandrina'},'Ammersee','Blelham','Bourget','Cannonsville',...
    'Como','Constance','ElGergal','Emaiksoun','Esthwaite','Feeagh', ...
    'Geneva01','Geneva03','GrosseDhunn','Harp','Iseo','Kinneret03', ...
    'Kinneret97','Mendota','MtBold','Muggelsee','NamCo','Oneida', ...
    'Pusiano','Rappbode', 'Rassnitzersee','Ravn','Rotorua', ...
    'Stechlin','Tarawera','Toolik','Windermere','Woods','Zurich'];

numLakes = length(LakeNames);

%Loop through lakes and re run the initial configuration for each lake
for lake_i = 1:numLakes
    lakename = LakeNames{lake_i};
        disp(['Current Lake:  ',lakename]);
     lake.(LakeNames{lake_i}).fit_params = rerunGLMinit(lakename,config_file);
     disp(['Model Fit All RMSE : ',num2str(lake.(LakeNames{lake_i}).fit_params.temp.all.RMSE)])
     disp(['Model Fit All PRE : ',num2str(lake.(LakeNames{lake_i}).fit_params.temp.all.PRE)])
    close all
    cd(base_dir)
end

save('MLCP_v2.2.0_MLCP.mat','lake')
