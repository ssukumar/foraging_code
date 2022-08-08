function simulate_jbar_change()
close all
% First start by simulating the energy expenditure wrt mass

mh = 4.5; % Base Mass =1 
ml = 1.0;

tm = 0:0.01:5; % assume seconds
tm = tm+0.45;

subj_low = [15,77];
subj_high = [45, 277];

num_subjs = 2;
subjs = cell(1,num_subjs);
subjs{1} = subj_low;
subjs{2} = subj_high;
figure
subj_lines = {'-','--'};
clrs = {'b','g'};

for s = 1:num_subjs
    fprintf('\n############ Subj # %d: \t',s),
    % low mass 
    
    u_ml = energy_u(subjs{s}(1), subjs{s}(2),ml,tm);
    plot(tm, u_ml,strcat(clrs{s},subj_lines{1}), 'LineWidth',1.5);
    hold on
    [~,tm_min] = min(u_ml);
    xline(tm(tm_min),strcat(clrs{s},subj_lines{1}));
    tm_low = tm(tm_min);
%     plot(tm(tm_min), 0, '+','MarkerSize',ml+5);
    fprintf('Low Mass : tm = %.3f;\t',tm(tm_min))
    
    %high mass
    u_mh = energy_u(subjs{s}(1), subjs{s}(2),mh,tm);
    plot(tm, u_mh,strcat(clrs{s},subj_lines{2}), 'LineWidth',1.5);
    hold on
    [~,tm_min] = min(u_mh);
    xline(tm(tm_min),strcat(clrs{s},subj_lines{2}));
%     plot(tm(tm_min), 0, '+','MarkerSize',mh+5);
    hold on
    fprintf('High Mass : tm = %.3f\n',tm(tm_min))
    tm_high = tm(tm_min);
    delta_tm = tm_high - tm_low;
    
    fprintf('############ \t\t \\Delta Travel Duration = %.3f\n',delta_tm);
    
end
xlabel('Movement Duration');
ylabel('Movement Expenditure');
legend('High Effort','Low Effort');

function u= energy_u (a,b, m, t)
% parameters from Shadmehr, Huang, Ahmed Curr Biol. (2016)
% Supplemental materials Table S1 : Movement Utility Row


%a = 15; b = 77; 
d = 0.30; % 30 cm
i = 1.1; j =2; k=.7;

u = a*t + b*(m^k)*(d^i)./(t.^j);
