%% beautify fig
set(gcf,'color',[ 1 1 1])
allhandles=findall(gcf);
allaxeshandles=findall(gcf,'type','axes');
alltitlehandles=get(allaxeshandles,'title');
allxlabelhandles=get(allaxeshandles,'xlabel');
allylabelhandles=get(allaxeshandles,'ylabel');

set(allaxeshandles,'tickdir','out','box','off','color','w','fontname','arial','fontsize',12,'fontweight','bold');

% if iscell(alltitlehandles) 
%     alltitlehandles=cell2mat(alltitlehandles); 
% end
% set(alltitlehandles,'fontname','arial','fontsize',14,'fontweight','bold');
% 
% if iscell(allxlabelhandles) allxlabelhandles=cell2mat(allxlabelhandles); end
% if iscell(allylabelhandles) allylabelhandles=cell2mat(allylabelhandles); end
% 
% set([allxlabelhandles allylabelhandles],'fontname','arial','fontsize',12,'fontweight','bold');
