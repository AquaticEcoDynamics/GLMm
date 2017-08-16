function [data] = readGLMnetcdf(filename,namelist)
% function [data] = nldnc(filename,namelist)
%
% Inputs:
%		filename   :  a netcdf file to take slice of variables from
%       namelist   :  variable names to load eg: 'TEMPERATURE'
% Outputs
%		data    : a matlab structure that contains the data from the
%		variables in namelist.
%
% Uses: netcdf, ncinfo, names_netcdf
%
% Written by L. Bruce 5 March 2013
% Based on nldncGLM B. Busch 2011
%

%Open file and create netcdf file identifier ncid
ncid = netcdf.open(filename,'nc_nowrite');

%Extract information from netcdf file
finfo = ncinfo(filename);

%Extract list of variable names
vnamesfull = {finfo.Variables.Name}';

%Determine varaibles to extract from file
if ~exist('namelist','var') || isempty(namelist); %Default extract all varaibles from file
	vnames = vnamesfull;
else
	if iscell(namelist)
        vnames = namelist;
    else
		vnames = {namelist};
	end
end

for ii = 1:length(vnames)
  found=strcmp(vnames(ii),vnamesfull);
  if ~isempty(found)
	  disp(['loading ',char(vnames(ii))]);
	  varid = netcdf.inqVarID(ncid,char(vnames(ii)));
      %Need to reverse order of dimensions
      data.(vnames{ii}) = permute(netcdf.getVar(ncid,varid),length(size(netcdf.getVar(ncid,varid))):-1:1);
  else %File does not contain variable ii
	found=strcmp(vnames(ii),namelist);
	if ~isempty(found)
	  disp(['Error loading ',char(vnames(ii)),' Variable not in file']);
	end
  end
end

%Determine start time date hours measured from
if ~strcmp('time',vnames) %then need to get time
	disp('loading time');
	varid = netcdf.inqVarID(ncid,'time');
    %Need to reverse order of dimensions
    data.time = permute(netcdf.getVar(ncid,varid),length(size(netcdf.getVar(ncid,varid))):-1:1);
end

time_i = find(strcmp('time',vnamesfull)==1);
data.startTime = datenum(finfo.Variables(time_i).Attributes.Value(end-19:end),'yyyy-mm-dd HH:MM');

%Convert data.time from hours from startTime to matlab date number
data.time = data.startTime + data.time/24;

%Get NS if not included in vnames
if ~strcmp('NS',vnames) %then need to get NS
	disp('loading NS');
	varid = netcdf.inqVarID(ncid,'NS');
    %Need to reverse order of dimensions
    data.NS = permute(netcdf.getVar(ncid,varid),length(size(netcdf.getVar(ncid,varid))):-1:1);
end

%Get depths if not included in vnames
if ~strcmp('z',vnames) %then need to get time
	disp('loading z');
	varid = netcdf.inqVarID(ncid,'z');
    %Need to reverse order of dimensions
    data.z = permute(netcdf.getVar(ncid,varid),length(size(netcdf.getVar(ncid,varid))):-1:1);
end

netcdf.close(ncid)