function fld_data = readLAwtr(filename)
% Usage: fld_data = readALwtr(filename)
%
% Inputs:
%		filename     : Lake Analyser field data .wtr file name
%                     
% Outputs:
%       fld_data     : field data structure
%
% Uses: 
%       importdata
%
% Simple script to load field data from a Lake Analyser field data .wtr file.
% 
% 
% Created by Louise Bruce 23/4/2013
% Modified for Windows 2/2/2014 importdata doesn't work for Como MATLAB2014

%import data from field and bathymetry files
data = importdata_bb(filename);
%bathdata = importdata([filename(1:end-9) '.bth']);

%Get bathymetry data - use glm2.nml so no longer needed
%fld_data.bthD   = bathdata.data(:,1)';
%fld_data.bthA   = bathdata.data(:,2)';

%First text data column is time
fld_data.time = data.mdates;

%Depth from surface level
fld_data.depth = data.headers;

%Temperature matix
fld_data.temp = data.data;