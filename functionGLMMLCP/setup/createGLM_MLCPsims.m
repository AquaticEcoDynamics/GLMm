%createGLM_MLCPsims
%
% Inputs:
%      LakeNames:  Array of strings corresponding to each of the lakes used
%      in the MLCP meta analysis
%
% Outputs:
%
% Uses:
%
%
% Written by L. Bruce 24 February 2014
% Modified by L. Bruce 14 April 2014 for MCMC Set up
%
% Creates a series of working sim directories for each of the lakes in the
% GLM Mulit-Lake Comparison Project (GLM_MLCP).  Using the latest
% input files and configuration from the previous runs.

%Clean up
close all
clear all

LakeNames = [{'Alexandrina'},'Ammersee','Blelham','Bourget','Cannonsville',...
    'Como','Constance','ElGergal','Emaiksoun','Esthwaite','Feeagh', ...
    'Geneva01','Geneva03','GrosseDhunn','Harp','Iseo','Kinneret03', ...
    'Kinneret97','Mendota','MtBold','Muggelsee','NamCo','Oneida', ...
    'Pusiano','Rappbode', 'Rassnitzersee','Ravn','Rotorua', ...
    'Stechlin','Tarawera','Toolik','Windermere','Woods','Zurich'];
numLakes = length(LakeNames);

%Location of original project directory
path_org = 'C:\Louise\GLM\GLM_v2.2.0\';

%Location of new project direcory
path_new = 'C:\Louise\GLM\GLM_v2.2.0_MLCP\';

%Run through lakes:
%   1) Create directory named lakename
%   2) Copy glm_init.nml from lakename/nml directory
%   3) Copy all input files *.csv from lakename/sim directory
%   4) Create Output and Results folders
%   4) Test run ./glm.exe

for lake_i = 4:numLakes
    
    lakename = LakeNames{lake_i};
    
    disp(['Current Lake:  ',lakename]);

    %Create LakeName directory
    mkdir([path_new lakename]);
    
    %Create subdirectories
     mkdir([path_new lakename '/nml']);
     mkdir([path_new lakename '/sim']);
     mkdir([path_new lakename '/Field']);
     mkdir([path_new lakename '/Ouput']);
     mkdir([path_new lakename '/Results']);
     
    %Copy across initial glm.nml file
    copyfile([path_org lakename '/nml/glm2_init.nml'],[path_new lakename '/nml/glm2_init.nml']);
    copyfile([path_org lakename '/nml/glm2_init.nml'],[path_new lakename '/sim/glm2.nml']);
    
    %Copy across input files
    copyfile([path_org lakename '\sim\*.csv'],[path_new lakename '\sim\']);
    
    %Copy field file
    copyfile([path_org lakename '\Field\' lakename '_Fld_temp.wtr'],[path_new lakename '\Field\']);
    
    %Test run
    orgpth=pwd;
    cd([path_new lakename '\sim'])
    !C:\Louise\GLM\glm_2.2.0rc4_win64\glm-bin\glm.exe
    cd(orgpth) 
end