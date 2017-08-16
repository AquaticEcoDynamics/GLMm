%Handy script to set the font size of the current plot

%Needs preset variable: xlabel_fontsize & y_label_fontsize

%Written by B. Busch 20th March 2014

ylab = get(gca,'YLabel');
set(gca,'YLabel',ylab,'fontsize',xlabel_fontsize);

xlab = get(gca,'XLabel');
set(gca,'XLabel',xlab,'fontsize',ylabel_fontsize);