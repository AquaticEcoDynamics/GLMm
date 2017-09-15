function glm_6_panel
% A very simple plot to veiw 3 variables, across 6 axis. Doesn't require
% any other functions.

% User input required below.....

datearray = datenum(2004,05:06:29,01);

image_name = 'TSO.png';

ncfile = 'output.nc';

surface_range = [65 75];
bottom_range = [0 10];

% Top Plot
var1.name = 'temp'; % GLM Variable Name
var1.caxis = [10 25]; % Axis limits for both plots
var1.Label = 'Temperature (C)';
var1.conversion = 1;
var1.legend_location = 'northwest';

% Middle Plot
var2.name = 'salt';% GLM Variable Name
var2.caxis = [0 2];
var2.Label = 'Salinity (psu)';
var2.conversion = 1;
var2.legend_location = 'southwest';

% Bottom Plot
var3.name = 'OXY_oxy';% GLM Variable Name
var3.caxis = [180 350];
var3.Label = 'Oxygen (mmol/m^3)';
var3.conversion = 1;
var3.legend_location = 'northeast';



% End of User Input_______________________________________________________


% Loading the external data;



% Shouldn't need to change anything under here...

fig = figure('DefaultAxesFontSize',7);

% Load Variable 1 (Top plot).

[data,XX,YY,ZZ,mTime,surf,bot] = glm_exportdata(ncfile,var1.name,surface_range,bottom_range);

axes('position',[0.05 0.65 0.4 0.25]); % Top Left

pcolor(XX,YY,ZZ);shading flat;
caxis(var1.caxis);
title(var1.Label,'fontsize',7,'fontweight','bold');
cb = colorbar;
set(cb,'position',[0.46 0.65 0.01 0.25]);
set(gca,'xtick',datearray,'xticklabel',[]);

axes('position',[0.55 0.65 0.4 0.25]); % Top Right

yt = get(gca,'ytick');
set(gca,'ytick',yt,'fontsize',7);

plot(mTime,surf,'r','displayname','GLM Surface');hold on
plot(mTime,bot,'--k','displayname','GLM Bottom');hold on




set(gca,'xtick',datearray,'xticklabel',[]);
%title(var1.Label,'fontsize',8,'fontweight','bold');
leg = legend('location',var1.legend_location);
set(leg,'fontsize',6);
ylim(var1.caxis);
xlim([datearray(1) datearray(end)]);
grid on;

%________________________________________________________________________
% Load Variable 2 (Middle plot).

[data,XX,YY,ZZ,mTime,surf,bot] = glm_exportdata(ncfile,var2.name,surface_range,bottom_range);

axes('position',[0.05 0.35 0.4 0.25]); % Mid Left

pcolor(XX,YY,ZZ);shading flat;
caxis(var2.caxis);
title(var2.Label,'fontsize',7,'fontweight','bold');
cb = colorbar;
set(cb,'position',[0.46 0.35 0.01 0.25]);
set(gca,'xtick',datearray,'xticklabel',[]);

axes('position',[0.55 0.35 0.4 0.25]); % Mid Right

plot(mTime,surf,'r','displayname','GLM Surface');hold on
plot(mTime,bot,'--k','displayname','GLM Bottom');hold on

%________________________________________________________________________
set(gca,'xtick',datearray,'xticklabel',[]);
%title(var2.Label,'fontsize',8,'fontweight','bold');
leg = legend('location',var2.legend_location);
set(leg,'fontsize',6);
ylim(var2.caxis);
xlim([datearray(1) datearray(end)]);
grid on;
% Load Variable 3 (Bottom plot).

[data,XX,YY,ZZ,mTime,surf,bot] = glm_exportdata(ncfile,var3.name,surface_range,bottom_range);

axes('position',[0.05 0.05 0.4 0.25]); % Bot Left

pcolor(XX,YY,ZZ);shading flat;
caxis(var3.caxis);
title(var3.Label,'fontsize',7,'fontweight','bold');
cb = colorbar;
set(cb,'position',[0.46 0.05 0.01 0.25]);
set(gca,'xtick',datearray,'xticklabel',datestr(datearray,'dd/mm/yy'),'fontsize',7);

axes('position',[0.55 0.05 0.4 0.25]); % Bot Right

plot(mTime,surf,'r','displayname','GLM Surface');hold on
plot(mTime,bot,'--k','displayname','GLM Surface');hold on

set(gca,'xtick',datearray,'xticklabel',datestr(datearray,'dd/mm/yy'),'fontsize',7);
%title(var3.Label,'fontsize',8,'fontweight','bold');
leg = legend('location',var3.legend_location);
set(leg,'fontsize',6);
ylim(var3.caxis);
xlim([datearray(1) datearray(end)]);
grid on;

%--% Paper Size
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
xSize = 21;
ySize = 24;
xLeft = (21-xSize)/2;
yTop = (30-ySize)/2;
set(gcf,'paperposition',[0 0 xSize ySize])

print(gcf,'-dpng',image_name,'-opengl');

close all;
end

function [data,XX,YY,ZZ,mTime,surf,bot] = glm_exportdata(ncfile,varname,surface,bottom)

% Just an array - ignore for now
alldepths = [0:1:2500];


% Import in the data
data = readGLMnetcdf(ncfile,varname);

%Number of dates simulated output
numDates = length(data.time);

%Get time and depths
mTime = data.time;
mDepth = data.z;
mNS = data.NS;

% Routine to find the max depth and build the array for the pcolor plots
for i = 1:length(mTime) 
    mdepth(i) = mDepth(i,mNS(i));
end

max_depth = max(mdepth);

[~,ind] = min(abs(alldepths-max_depth));

if alldepths(ind) < max_depth
    ind = ind + 1;
end

nd = [0:0.1:alldepths(ind)];

[XX,YY] = meshgrid(mTime,nd);

% Now get the data for the pcolor plot

for i = 1:length(mTime)    
    z = mDepth(i,1:mNS(i))';
    x = data.(varname)(i,1:mNS(i))';
    ZZ(:,i) = interp1(z,x,nd);
    
    % Data for line plot.
    sInd = find(z >= surface(1) & z <= surface(2));
    surf(i) = mean(x(sInd));
    bInd = find(z >= bottom(1) & z <= bottom(2));
    bot(i) = mean(x(bInd));
    
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

netcdf.close(ncid)
end