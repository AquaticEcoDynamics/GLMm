function plotLake(paths,filename)
%
%
% Simple script to load and plot data from a GLM output.nc file.
% 
% 
% Uses: plotGLM
%       readGLMnetcdf
%
% Created by Brendan Busch
% Modified Louise Bruce 5/3/13 for GLM v1.0.2


foldername = [paths.working_dir,'Output/'];
outname = [paths.working_dir,'Output/Plots'];
mkdir(outname);

%Extract information from netcdf file
finfo = ncinfo([foldername,filename]);

%Extract list of variable names
varNames = {finfo.Variables.Name}';

%if length(varNames)<=19
    %validNames = varNames([15,16,17]);
%else
    %validNames = varNames([15:20,25:27,29:35,38,41,44,47,50:52]);
%end
%testsim1
%validNames = varNames([15:20,25:27,29:35,38,41,44])
%testsim2
%validNames = varNames([15:20,25:27,29:35,38,41])
%validNames = {'aed_phytoplankton_green'};
%validNames = varNames([15:20,25:27,29:34,35,38,41,44])
%validNames = varNames([16,18,19,34,51,52,92]) %Mount Bold
%validNames = varNames([16,18,19,34,35,38,41]) %Lake Alex
validNames = varNames(16);

data = readGLMnetcdf([foldername,filename],['time','NS','z',validNames']);

for ii = 1:length(validNames)
    finalOut = [outname,'/',validNames{ii},'.png'];
    newFig = plotGLM(validNames{ii},data);
    
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'paperposition',[0.1  0.1 15 8]);
    
    %-----------------------------------------------------------------------------
    %          Finally save figure
    % ------------------------------------ 
    %print(gcf,'-dpng',finalOut,'-opengl');
    %closes the figure
    %close all
    
end