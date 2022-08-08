close all

fig1 = openfig('FRM_profile_env1_exp_1.fig');
ylabel('Max Grip Force Rate (N/s)');
xlabel('Sample #');
title('Low Effort Env','Fontsize',13);
beautifyfig


axObjs = fig1.Children;
dataObjs = axObjs.Children;
% xprb_dat = dataObjs(1).XData;
subplt = 0;
len_ =150;
plot_spacing = 170;

if ~subplt
    h=figure;
    ylow = dataObjs(3).YData(1:len_);
    xlow = 1:length(ylow);
    ylow_ser = dataObjs(4).YData(1:length(ylow))' - ylow;
    gcf;
    plotsdshading(xlow, ylow', ylow_ser','#F8766D');
    
    hold on
    yprb_low = dataObjs(1).YData(1:len_);
    xprb_low = (1:length(yprb_low))+ plot_spacing;
    yprb_low_ser = dataObjs(2).YData(1:length(yprb_low))' - yprb_low;
    plotsdshading(xprb_low, yprb_low', yprb_low_ser','#A7ACD7');
    
    fig2 = openfig('FRM_profile_env2_exp_1.fig');
    axObjs2 = fig2.Children;
    dataObjs2 = axObjs2.Children;
    figure(h)
    hold on
    yhigh = dataObjs2(3).YData(1:len_);
    yhigh_ser = dataObjs2(4).YData(1:len_)' - yhigh;
    yprb_high = dataObjs2(1).YData(1:len_);
    yprb_high_ser = dataObjs2(2).YData(1:len_)' - yprb_high;
    xprb_high = (1:length(yprb_high))+ 2*plot_spacing;
    xhigh = (1:length(yhigh))+ 3*plot_spacing;
    
    plotsdshading(xprb_high, yprb_high', yprb_high_ser', '#A7ACD7');
    hold on
    plotsdshading(xhigh, yhigh', yhigh_ser', '#00BFC4');

    ylabel('Max Grip Force Rate (N/s)');
    beautifyfig
    
    keyboard
    
end

fig2 = openfig('FRM_profile_env2_exp_1.fig');
ylabel('Max Grip Force Rate (N/s)');
xlabel('Sample #');
title('High Effort Env','Fontsize',13);
beautifyfig
