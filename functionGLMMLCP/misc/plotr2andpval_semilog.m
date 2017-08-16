%plotr2andpval

%Handy script to plot the r2 and pval in the top right hand corner of the
%graph

%Created by L. Bruce 19 December 2013

%Log can't handle negatives so use
y1 = y(x>0);
x1 = x(x>0);

[R,P] = corrcoef(log(x1(~isnan(y1))),y1(~isnan(y1)));
r2=R(1,2);
pval=P(1,2);

ab = polyfit(log(x1(~isnan(y1))),y1(~isnan(y1)),1);


if r2 > 0
    text(0.07,0.95,['r=',num2str(r2,'%1.2f')],'Units','Normalized','Fontsize',figure_fontsize);%,'Color',[.7 .9 .7]);
    text(0.07,0.85,['p=',num2str(pval,'%1.2f')],'Units','Normalized','Fontsize',figure_fontsize);%,'Color',[.7 .9 .7]);
else
    text(0.7,0.95,['r=',num2str(r2,'%1.2f')],'Units','Normalized','Fontsize',figure_fontsize);%,'Color',[.7 .9 .7]);
    text(0.7,0.85,['p=',num2str(pval,'%1.2f')],'Units','Normalized','Fontsize',figure_fontsize);%,'Color',[.7 .9 .7]);
end