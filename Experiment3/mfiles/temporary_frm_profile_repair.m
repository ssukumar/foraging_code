function temporary_frm_profile_repair()
% Function to beautify the velocity profiles as much as possible before
% upload to manuscript
close all
addpath(genpath('../plots'));
fig1 = openfig('FRM_profile_env1_exp_1.fig');

ylabel('Max Grip Force Rate (N/s)');
beautifyfig
axObjs = fig1.Children;
dataObjs = axObjs.Children;
% xprb_dat = dataObjs(1).XData;
subplt = 0;
if ~subplt
    h=figure;
    ylow = dataObjs(3).YData(275:725);
    xlow = 1:length(ylow);
    ylow_ser = dataObjs(4).YData(275:725)' - dataObjs(3).YData(275:725);
    
    gcf;
    % plot(xlow,ylow,'r', 'LineWidth',1.5);
    plotsdshading(xlow, ylow', ylow_ser','#F8766D');
    
    hold on
    yprb_low = dataObjs(1).YData(200:650);
    xprb_low = (1:length(yprb_low))+ 500;
    yprb_low_ser = dataObjs(2).YData(200:650)' - yprb_low;
    % plot(xprb_low,yprb_low,'b', 'LineWidth',1.5);
    plotsdshading(xprb_low, yprb_low', yprb_low_ser','#A7ACD7');
    fig2 = openfig('FRM_profile_env2_exp_1.fig');
    axObjs2 = fig2.Children;
    dataObjs2 = axObjs2.Children;
    figure(h)
    hold on
    yhigh = dataObjs2(3).YData(740:1300);
    yhigh_ser = dataObjs2(4).YData(740:1300)' - yhigh;
    yprb_high = dataObjs2(1).YData(440:900);
    yprb_high_ser = dataObjs2(2).YData(440:900)' - yprb_high;
    
    xprb_high = (1:length(yprb_high))+1000;
    xhigh = (1:length(yhigh))+1500;
    % plot(xprb_high, yprb_high,'b','LineWidth',1.5);
    plotsdshading(xprb_high, yprb_high', yprb_high_ser', '#A7ACD7');
    
    hold on
    % plot(xhigh, yhigh,'g','LineWidth',1.5);3
    plotsdshading(xhigh, yhigh', yhigh_ser', '#00BFC4');
    % set(gca,'xtick',[])
    
    ylabel('Max Grip Force Rate (N/s)');
    beautifyfig
    
    keyboard
    
else
    fig_h = figure;
    ylow = dataObjs(3).YData(275:725);
    xlow = 1:length(ylow);
    ylow_ser = dataObjs(4).YData(275:725)' - dataObjs(3).YData(275:725);
    hold on 
    plotsdshading(xlow, ylow', ylow_ser','#F8766D');
    
    yprb_low = dataObjs(1).YData(200:650);
    xprb_low = (1:length(yprb_low)) + 500;
    yprb_low_ser = dataObjs(2).YData(200:650)' - yprb_low;
    plotsdshading(xprb_low, yprb_low', yprb_low_ser','#A7ACD7');
    hold on
    
    fig2 = openfig('FRM_profile_env2_exp_1.fig');
    axObjs2 = fig2.Children;
    dataObjs2 = axObjs2.Children;
    hold on
    yhigh = dataObjs2(3).YData(740:1300);
    yhigh_ser = dataObjs2(4).YData(740:1300)' - yhigh;
    yprb_high = dataObjs2(1).YData(440:900);
    yprb_high_ser = dataObjs2(2).YData(440:900)' - yprb_high;
    xprb_high = (1:length(yprb_high))+1000;
    xhigh = (1:length(yhigh))+1500;
    figure(fig_h)
    hold on
    plotsdshading(xprb_high, yprb_high', yprb_high_ser', '#A7ACD7');
    hold on
    plotsdshading(xhigh, yhigh', yhigh_ser', '#00BFC4');
    xticks([(xlow(1:100:end)-1), (xprb_low(1:100:end)-1),...
        (xprb_high(1:100:end)-1), (xhigh(1:100:end)-1)])
    xticks_labels = string(5*[(xlow(1:100:end)-1), ...
        (xprb_low(1:100:end)-501),(xprb_high(1:100:end)-1001) ,...
        (xhigh(1:100:end)-1501)]);
    xticklabels(cellstr(xticks_labels))
    ylabel('Max Grip Force Rate (N/s)');
    beautifyfig

end