%plotGLM_MLCP1x4figure

%Script to plot a scatter plot in a figure with with 1 row of 1 square subplots
%Used in the GLM MLCP Sensitivity Analysis for correlation between lake
%characteristics and sensitivity coefficients for physical parameters.

%Created by L. Bruce 15th December 2014
%Settings pulled from B. Busch circa 2013
%Added box plot to right for MLCP Paper 29th September 2015

%Add box plot to right of scatter plot
axes('position',[0.85 0.1 0.12 0.85])
boxplot(y)

%Main plot
axes('position',[0.1 0.1 0.7 0.85])
if ~isempty(find(linear_chars == lchar_i))
    plot(x,y,'w*')
    %Plot significance of relationship
    plotr2andpval
else
    semilogx(x,y,'w*')
    %Plot significance of relationship
    plotr2andpval_semilog
end

%Plot lakes as 2 letter initials
for lake_i = 1:numLakes
    text(x(lake_i),y(lake_i),LakeInitials{lake_i}, ...
        'HorizontalAlignment','center','Fontsize',6, ...
        'Color',LakeColour{lake_i},'Fontweight',LakeFont{lake_i})
end

%Figure annotation
x_lim=xlim; y_lim = ylim;
if ~isempty(find(linear_chars == lchar_i))
   text_x = x_lim(1)+0.95*(x_lim(2) - x_lim(1));
   text_y = y_lim(1)+0.05*(y_lim(2) - y_lim(1));
else
   text_x = exp(log(x_lim(1))+0.95*(log(x_lim(2)) - log(x_lim(1))));
   text_y = y_lim(1)+0.05*(y_lim(2) - y_lim(1));
end

