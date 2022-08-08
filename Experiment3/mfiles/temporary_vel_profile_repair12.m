function temporary_vel_profile_repair12()
% Function to beautify the velocity profiles as much as possible before
% upload to manuscript
close all
addpath(genpath('../plots'));
fig1 = openfig('vel_profile_temp_low_YRT.fig');
ylabel('Velocity (m/s)');
beautifyfig
axObjs = fig1.Children;
dataObjs = axObjs.Children;
% xprb_dat = dataObjs(1).XData;
subplt = 1;
sidebyside = 0;
window_size_half = 250;

if ~subplt
    h=figure;
    
    idxs = 921:1449;
    ylow = dataObjs(3).YData(idxs);
    xlow = 1:length(ylow);
    ylow_ser = dataObjs(4).YData(idxs)' - dataObjs(3).YData(idxs);
    
    gcf;
    % plot(xlow,ylow,'r', 'LineWidth',1.5);
    plotsdshading(xlow, ylow', ylow_ser','#F8766D');
    
    hold on
    yprb_low = dataObjs(1).YData(200:650);
    xprb_low = (1:length(yprb_low))+ 500;
    yprb_low_ser = dataObjs(2).YData(200:650)' - yprb_low;
    % plot(xprb_low,yprb_low,'b', 'LineWidth',1.5);
    plotsdshading(xprb_low, yprb_low', yprb_low_ser','#A7ACD7');
    fig2 = openfig('vel_profile_temp_high_YRT.fig');
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
    
    ylabel('Velocity (m/s)');
    beautifyfig
    
    keyboard
    
elseif sidebyside
    fig_h = figure;
    
    max_idx_low = 824;
    idx_low =max_idx_low-window_size_half:max_idx_low+window_size_half;
    ylow = dataObjs(3).YData(idx_low);
    xlow = 1:length(ylow);
    ylow_ser = dataObjs(4).YData(idx_low)' - dataObjs(3).YData(idx_low);
    hold on 
    plotsdshading(xlow, ylow', ylow_ser','#F8766D');
    
    max_idx_prb_low = 893;
    idxs_prb_low = max_idx_prb_low-window_size_half:max_idx_prb_low+window_size_half;
    yprb_low = dataObjs(1).YData(idxs_prb_low);
    xprb_low = (1:length(yprb_low)) + 800;
    yprb_low_ser = dataObjs(2).YData(idxs_prb_low)' - yprb_low;
    plotsdshading(xprb_low, yprb_low', yprb_low_ser','#A7ACD7');
    hold on
    
    fig2 = openfig('vel_profile_temp_high_YRT.fig');
    axObjs2 = fig2.Children;
    dataObjs2 = axObjs2.Children;
    hold on
    
    idx_max_high = 781;
    idx_high = idx_max_high-window_size_half:idx_max_high+window_size_half;
    yhigh = dataObjs2(3).YData(idx_high);
    yhigh_ser = dataObjs2(4).YData(idx_high)' - yhigh;
    
    idx_max_prb_high = 522;
    idx_prb_high = idx_max_prb_high-window_size_half:idx_max_prb_high+window_size_half;
    yprb_high = dataObjs2(1).YData(idx_prb_high);
    yprb_high_ser = dataObjs2(2).YData(idx_prb_high)' - yprb_high;
    xprb_high = (1:length(yprb_high))+1600;
    xhigh = (1:length(yhigh))+2400;
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
    ylabel('Velocity (m/s)');
    beautifyfig

    
else 
    % still need to figure out indexing for dataObj and note it down this
    % time
    
    % ## Low effort
    [~,max_idx_low] = max(dataObjs(3).YData);
    idx_low =max_idx_low-window_size_half:max_idx_low+window_size_half;
    ylow = dataObjs(3).YData(idx_low);
    [~,new_max_ix_ylow] = max(ylow);
    ylow_ser = dataObjs(4).YData(idx_low)' - dataObjs(3).YData(idx_low);
    
    
    [~,max_idx_prb_low] = max(dataObjs(2).YData);
    idxs_prb_low = max_idx_prb_low-window_size_half:max_idx_prb_low+window_size_half;
    yprb_low = dataObjs(1).YData(idxs_prb_low);
    [~,new_max_ix_yprb_low] = max(yprb_low);
    yprb_low_ser = dataObjs(2).YData(idxs_prb_low)' - yprb_low;
    figure;
    
    ylow_diff =[];
    
    xlow = 1:length(ylow);
    xprb_low = (1:length(yprb_low));
    
    
    if new_max_ix_ylow <= new_max_ix_yprb_low
        xlow = xlow + (new_max_ix_yprb_low- new_max_ix_ylow);
        ylowhigh_diff = -1*yprb_low + [zeros(1,(new_max_ix_yprb_low- new_max_ix_ylow)),...
            ylow];
        ylowhigh_ser = sqrt(yprb_low_ser.^2 + [zeros(1,(new_max_ix_yprb_low- ...
            new_max_ix_ylow)), ylow_ser].^2);
    else 
        xprb_low = xprb_low + (new_max_ix_ylow- new_max_ix_yprb_low);
        ylowhigh_diff = -1*[zeros(1,(new_max_ix_yprb_low- new_max_ix_ylow)),yhigh]...
            + ylow;
        ylowhigh_ser = sqrt([zeros(1,(new_max_ix_yprb_low- new_max_ix_ylow))...
            ,yprb_low_ser].^2 + ylow_ser.^2);
    end
    
    
    
    hold on 
    plotsdshading(xlow, ylow', ylow_ser','#F8766D');
    hold on
    plotsdshading(xprb_low, yprb_low', yprb_low_ser','#00BFC4');
    
    hold on 
    plotsdshading(xlow, ylowhigh_diff', ylowhigh_ser','k');
    
    ylabel('Velocity (m/s)');
    title('Low Effort');
    beautifyfig
    
    %## high effort 
    fig2 = openfig('high_reward_vel_prof_subj14_2.fig');
    axObjs2 = fig2.Children;
    dataObjs2 = axObjs2.Children;
    hold on
    [~,idx_max_prb_high] = max(dataObjs2(1).YData);
    idx_prb_high = idx_max_prb_high-window_size_half:idx_max_prb_high+window_size_half;
    yprb_high = dataObjs2(1).YData(idx_prb_high);
    [~,new_max_ix_yprb_high] = max(yprb_high);
    yprb_high_ser = dataObjs2(2).YData(idx_prb_high)' - yprb_high;
    figure;
    
    [~,idx_max_high] = max(dataObjs2(3).YData);
    idx_high = idx_max_high-window_size_half:idx_max_high+window_size_half;
    yhigh = dataObjs2(3).YData(idx_high);
    yhigh_ser = dataObjs2(4).YData(idx_high)' - yhigh;
    figure;
    
    ylow_diff =[];
    [~,new_max_ix_yhigh] = max(yhigh);
    xprb_high = 1:length(yprb_high);
    xhigh = (1:length(yhigh));
    
    
    if new_max_ix_ylow <= new_max_ix_yhigh
        xprb_high = xprb_high + (new_max_ix_yhigh- new_max_ix_yprb_high);
        ylowhigh_diff =  yhigh - [zeros(1,(new_max_ix_yhigh-...
            new_max_ix_yprb_high)),...
            yprb_high];
        ylowhigh_ser = sqrt(yhigh_ser.^2 + [zeros(1,(new_max_ix_yhigh- ...
            new_max_ix_yprb_high)), yprb_high_ser].^2);
    else 
        xhigh = xhigh + (new_max_ix_yprb_high- new_max_ix_yhigh);
        ylowhigh_diff = [zeros(1,(new_max_ix_yhigh- new_max_ix_yprb_high)),yhigh]...
            - yprb_high;
        ylowhigh_ser = sqrt([zeros(1,(new_max_ix_yhigh- new_max_ix_yprb_high))...
            ,yhigh_ser].^2 + yprb_high_ser.^2);
    end
    
    
    
    hold on 
    plotsdshading(xprb_high, yprb_high', yprb_high_ser','#F8766D');
    hold on
    plotsdshading(xhigh, yhigh', yhigh_ser','#00BFC4');
    
    hold on 
    plotsdshading(xhigh, ylowhigh_diff', ylowhigh_ser','k');
    
    ylabel('Velocity (m/s)');
    title('High Effort');
    beautifyfig
    
    % ### Probe 
    
    [~,max_idx_prb_low] = max(dataObjs(2).YData);
    idxs_prb_low = max_idx_prb_low-window_size_half:max_idx_prb_low+window_size_half;
    yprb_low = dataObjs(1).YData(idxs_prb_low);
    [~,new_max_ix_yprb_low] = max(yprb_low);
    yprb_low_ser = dataObjs(2).YData(idxs_prb_low)' - yprb_low;
    
    [~,idx_max_prb_high] = max(dataObjs2(1).YData);
    idx_prb_high = idx_max_prb_high-window_size_half:idx_max_prb_high+window_size_half;
    yprb_high = dataObjs2(1).YData(idx_prb_high);
    figure;
    
    [~,new_max_ix_yprb_high] = max(yprb_high);
    yprb_high_ser = dataObjs2(2).YData(idx_prb_high)' - yprb_high;
    xprb_low = 1:length(yprb_low);
    xprb_high = (1:length(yprb_high));
    if new_max_ix_yprb_low <= new_max_ix_yprb_high
        xprb_low = xprb_low + (new_max_ix_yprb_high- new_max_ix_yprb_low);
        yprb_lowhigh_diff = -1*yprb_high + [zeros(1,(new_max_ix_yprb_high- new_max_ix_yprb_low)),...
            yprb_low];
        yprb_lowhigh_ser = sqrt(yprb_high_ser.^2 + [zeros(1,...
            (new_max_ix_yprb_high- new_max_ix_yprb_low)), yprb_low_ser].^2);
    else 
        xprb_high = xprb_high + (new_max_ix_yprb_low- new_max_ix_yprb_high);
        yprb_lowhigh_diff = -1*yprb_high + [zeros(1,(new_max_ix_yprb_high- new_max_ix_yprb_low)),...
            yprb_low];
        yprb_lowhigh_ser = sqrt(yprb_high_ser.^2 + [zeros(1,...
            (new_max_ix_yprb_high- new_max_ix_yprb_low)), yprb_low_ser].^2);
    end
    
    
    
    hold on 
    plotsdshading(xprb_low, yprb_low', yprb_low_ser','#F8766D');
    
    plotsdshading(xprb_high, yprb_high', yprb_high_ser','#00BFC4');
    hold on
    plotsdshading(xprb_high, yprb_lowhigh_diff', yprb_lowhigh_ser','k');
    
    ylabel('Velocity (m/s)');
    title('Probe');
    beautifyfig

end