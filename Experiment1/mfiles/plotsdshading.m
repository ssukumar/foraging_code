function h=plotsdshading(xaxis,meandata,sddata,c)
%PLOTSDSHADING  Creates shaded region around plotted lines
%  Created by Alaa Ahmed on 22-Apr-2007
%  Last modified  22-Apr-2007

upperline=meandata+sddata;
lowerline=meandata-sddata;
xdata=[xaxis fliplr(xaxis) xaxis(1)];
ydata=[upperline' fliplr(lowerline') upperline(1)];

%Plot shaded area

h(1)=patch(xdata,ydata,[0,0,0],'FaceColor',c);

set(h(1),'FaceAlpha',0.2)
hold on

%Plot mean data
h(2)=plot(xaxis,meandata,'linewidth',1.5,'Color',c);

