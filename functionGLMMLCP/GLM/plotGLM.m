function H = plotGLM(varname,data)
% Usage: plotGLM(variablename,data)
%
% Inputs:
%		varname     : GLM variable name
%       data        : a matlab structure that contains all the data from
%                     the GLM output.nc file
% Outputs:
%		H           : handle to figure
%
% Uses: 
%       pcolor
%
% Simple script to load and plot data from a GLM output.nc file.
% 
% 
% Created by Brendan Busch
% Modified Louise Bruce 5/9/11
% Modified Louise Bruce 5/3/13 for GLM v1.0.2

arraySize = size(data.z);
%Resize variable matrix so that bottom values = lowest layer
% Needs to remove dud values - May not be required in future versions of
% the code
for ii = 1:length(data.time)
    plotVar(ii,1) = data.(varname)(ii,1);
    plotVar(ii,2:arraySize(2)+1) = data.(varname)(ii,:);
    plotVar(ii,data.NS(ii)+2) = data.(varname)(ii,data.NS(ii));
    plotVar(ii,(data.NS(ii)+3:end)) = NaN;
    height(ii,1) = 0;
    height(ii,2:arraySize(2)+1) = data.z(ii,:);
    height(ii,data.NS(ii)+2) = data.z(ii,data.NS(ii));
    height(ii,(data.NS(ii)+3:end)) = NaN;
    [maxHeight(ii) surf_i(ii)] = max(height(ii,:));
    surf_var(ii) = plotVar(ii,surf_i(ii));
end
disp('Finished doing the heights fix')

for ii = 1:length(data.time)
    NaNGuys = ~isnan(plotVar(ii,:));
    ss = find(NaNGuys >0);
    if(length(ss) > 3)
    	newHeight(ii,:) = 0:max(height(ii,:)/100):max(height(ii,:));
    	newVar(ii,:) = interp1(height(ii,ss(1:end-1)),plotVar(ii,ss(1:end-1)),newHeight(ii,:),'linear');
	else
	    newHeight(ii,:) = 0:max(height(ii,:)/100):max(height(ii,:));
		newVar(ii,1:length(newHeight(ii,:))) = NaN;
    end
    
    plotTime(ii,1) = data.time(ii);
end
disp('finished doing the interp bit..')

%Determine x tick marks
xmin = min(data.time);
xmax = max(data.time);
x_lim = [xmin xmax];
x_tick = [xmin:round((xmax-xmin)/2):xmax];
x_tick_date = datestr(x_tick,'dd-mmm-yy');

figure
varInformation
H = pcolor(plotTime,newHeight',newVar'); hold on;shading flat
  plot(plotTime,maxHeight,'k');hold on
  set(gca,'box','on')
 xlabel('Date')
 ylabel('Height (m)')
 xlim(x_lim)
set(gca,'XTick',x_tick,'XTickLabel',x_tick_date)
var_lim = eval(['varIndex.',varname,'.caxis']);
var_title = eval(['varIndex.',varname,'.title']);
 title(var_title)
 caxis(var_lim)
 colorbar
 

 


