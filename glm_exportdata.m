function [data,XX,YY,ZZ,mTime,surf,bot,datum] = glm_exportdata(ncfile,lakefile,surface,bottom,useDatum)

% If use datum is 1, then Z will be coverted to the bottom depth as
% specified in the gml2.nml file. File must be in the same directory as the
% ncfile.

% ncfile = 'output.nc';
% lakefile 'lake.csv';
% varname = 'temp';
% 
% bottom = [];
% surface = [];
% 
% useDatum = 1;


vars = glm_infonetcdf(ncfile);


first_var = 1;

lake_file = regexprep(ncfile,'output.nc','lake.csv');

import_lake_file(lake_file);

datum = 0;


if useDatum
    
    glm_file = 'glm2.nml';
    
    fid = fopen(glm_file,'rt');
    
    while ~feof(fid)
        fline = fgetl(fid);
        str = regexprep(fline,' ','');
        
        spt = strsplit(str,'=');
        if strcmpi(spt{1},'base_elev')
            datum = str2num(regexprep(spt{2},',',''));
        end
    end
end





% Just an array - ignore for now
alldepths = [-2500:10:2500];


% Loop through available vars;

for bb = 15:length(vars)
    
    if first_var
        
        % Import in the data
        data = readGLMnetcdf(ncfile,vars{bb});
        
        %Number of dates simulated output
        numDates = length(data.time);
        
        %Get time and depths
        mTime = data.time;
        mDepth = data.z;
        mNS = data.NS;
        
        mDepth(mDepth > 9.96920996838687e+30) = NaN;
        
        % Routine to find the max depth and build the array for the pcolor plots
        for i = 1:length(mTime)
            mdepth(i) = mDepth(i,mNS(i));
        end
        
        max_depth = max(mdepth);
        
        [~,ind] = min(abs(alldepths-max_depth));
        
        if alldepths(ind) < max_depth
            ind = ind + 1;
        end
        
        nd = [0:0.01:alldepths(ind)];
        
        [XX,YY] = meshgrid(mTime,nd);
        
        % Now get the data for the pcolor plot
        
        mLine.(vars{bb}).surf = [];
        mLine.(vars{bb}).bot = [];
        
        if isempty(surface)
            surface = [max(max(mDepth)) max(max(mDepth))-10];
            disp(['No surface values specifed, using ',num2str(surface(1)),' to ',num2str(surface(2))]);
        end
        
        if isempty(bottom)
            bottom = [0 10];
            disp(['No bottom values specifed, using ',num2str(bottom(1)),' to ',num2str(bottom(2))]);
        end
        
        
        
        
        for i = 1:length(mTime)
            z(1) = 0;
            z(2:(length(1:mNS(i))+1),1) = mDepth(i,1:mNS(i))';
            x(1) = data.(vars{bb})(i,1);
            x(2:(length(1:mNS(i))+1),1) = data.(vars{bb})(i,1:mNS(i))';
            if length(z) > 1
                ZZ(:,i) = interp1(z,x,nd);
            else
                ZZ(1:length(nd),i) = NaN;
            end
            
            % Data for line plot.
            
                
                sInd = find(z >= surface(1) & z <= surface(2));
                mLine.(vars{bb}).surf(i) = mean(x(sInd));
            
                bInd = find(z >= bottom(1) & z <= bottom(2));
                mLine.(vars{bb}).bot(i) = mean(x(bInd));
            clear z x;
        end
        
        
        if useDatum
            YY = YY + datum;
        end
        
        mData.(vars{bb}) = ZZ;
        
        first_var = 0;
        
        clear ZZ;
        
    else
        
        % Import in the data
        data = readGLMnetcdf(ncfile,vars{bb});
        mLine.(vars{bb}).surf = [];
        mLine.(vars{bb}).bot = [];
        for i = 1:length(mTime)
            z(1) = 0;
            z(2:(length(1:mNS(i))+1),1) = mDepth(i,1:mNS(i))';
            x(1) = data.(vars{bb})(i,1);
            x(2:(length(1:mNS(i))+1),1) = data.(vars{bb})(i,1:mNS(i))';
            if length(z) > 1
                ZZ(:,i) = interp1(z,x,nd);
            else
                ZZ(1:length(nd),i) = NaN;
            end
            
            % Data for line plot.
            
            if ~isempty(surface)
                
                sInd = find(z >= surface(1) & z <= surface(2));
                mLine.(vars{bb}).surf(i) = mean(x(sInd));
            end
            
            if ~isempty(bottom)
                bInd = find(z >= bottom(1) & z <= bottom(2));
                mLine.(vars{bb}).bot(i) = mean(x(bInd));
            end
            clear z x;
        end
        
        mData.(vars{bb}) = ZZ;
        
        
        clear ZZ;
        
    end
    
end

depth_bin.surface = surface;
depth_bin.bottom = bottom;

save mData.mat mData XX YY mTime mLine datum depth_bin -mat;
end

function [varnames,dimnames,data] = glm_infonetcdf(filename);
% Simple function to get all variable names from a Tuflow netcdf file
data = [];
ncid = netcdf.open(filename,'NC_NOWRITE');
[ndims,nvars,~,unlimdimid] = netcdf.inq(ncid);
dimids = (0:ndims-1)';
dimnames = cell(ndims,1);
dimlen = zeros(1,ndims);
for i = 1 : ndims
    [dimnames{i},dimlen(i)] = netcdf.inqDim(ncid,dimids(i));
end
varid = (0:nvars-1)';
varnames = cell(nvars,1);
xtype = zeros(nvars,1);
vardimids = cell(nvars,1);
varunlimdim = cell(nvars,1);
for i = 1 : nvars
    [varnames{i},xtype(i),vardimids{i}] = netcdf.inqVar(ncid,varid(i));
    varunlimdim{i} = find(vardimids{i}==unlimdimid,1,'first');
end

end
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

netcdf.close(ncid);
end
function lake = import_lake_file(filename)

disp(['Importing ',filename]);

fid = fopen(filename,'rt');

fline = fgetl(fid);

headers = strsplit(fline,',');

headers = regexprep(headers,' ','_');
headers = regexprep(headers,'/','_');


frewind(fid)
x  = length(headers);
textformat = [repmat('%s ',1,x)];

datacell = textscan(fid,textformat,'Headerlines',3,'Delimiter',',');
fclose(fid);

lake.Date(:,1) = datenum(datacell{1},'yyyy-mm-dd HH:MM:SS');

for i = 2:length(headers)
    lake.(headers{i}) = str2double(datacell{i});
end

save lake.mat lake -mat


end
