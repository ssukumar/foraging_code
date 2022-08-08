function berry_noise_determination ()
% Function to determine how far off the algorithm is in dispensing berries

% Simulate the berry dispensing
exp_no = 1;
% read in the data_subj array
if exp_no == 1
    filename = 'subj_order.txt';
    table_name= 'table_R.csv';
else
    filename = 'subj_order_e2.txt';
    table_name='table_R_e2.csv';
end

%Data directory
folder = '..';
f = pwd;
cd(folder);
dir_file = pwd;
cd(f);

close all



filename = fullfile(folder, filename);
csv_filename = fullfile(folder, table_name);

subj_order = readtable(filename, 'Delimiter',' ',...
    'ReadVariableNames',0);
subj_order = table2cell(subj_order);
num_subjs = size(subj_order,1);
num_trials = 200;
num_env=2;
max_berries = 16;
time_to_berry_env = cell(1,num_env);
avg_per_subj_gut = zeros(num_subjs,num_env);
clrs = {'b','r'};
for fn = 1:num_env
%     figure(fn);
    
    env = {'Low Effort', 'High Effort'};
    subj_ct = 1;
    time_to_berry_env{fn} = repmat({nan(1,max_berries)}, num_subjs,2);
    diff_time_env{fn} = repmat({nan(1,max_berries)}, num_subjs,2);
    norm_GUT = nan(num_subjs,num_trials);
    
    for i = [cell2mat(subj_order(:,1))]'
        subj_ct = i;
        subj_name = subj_order{subj_ct,2};
        order = subj_order{subj_ct,3};
        load(fullfile(folder, 'matfiles_data',subj_name,...
            'data_subj.mat'));
        time_to_berry_ = nan(num_trials ,max_berries);
        diff_time_to_berry = nan(num_trials ,max_berries);
        for tr = 1:num_trials
            if length(data_subj{fn}.time_to_berry{tr})>0
                time_to_berry_(tr,:) = data_subj{fn}.time_to_berry{tr} - ...
                    data_subj{fn}.grip_drop_berry{tr};
                diff_time_to_berry(tr,:) = [time_to_berry_(tr,1), ...
                    diff(time_to_berry_(tr,:))];
            end
            
%             if data_subj{fn}.berries(tr)>0
%                 norm_GUT(tr) = data_subj{fn}.giving_up_time(tr)/...
%                     diff_time_to_berry(data_subj{fn}.berries(tr));
%             end
        end
%         keyboard
        diff_time_env{fn}{i,1} = nanmean(diff_time_to_berry ,1)./1000;
        diff_time_env{fn}{i,2} = nanstd(diff_time_to_berry ,0,1)./1000;
        time_to_berry_env{fn}{i,1} = nanmean(time_to_berry_,1)./1000;
        time_to_berry_env{fn}{i,2} = nanstd(time_to_berry_,0,1)./1000;
        
%         if i <= 10
%             subplot(2,1,1);
%             title('\beta = 0.5','Interpreter','Tex');
%         else
%             subplot(2,1,2);
%             title('\beta = 0.2','Interpreter','Tex');
%         end
        
%         errorbar(1:max_berries, time_to_berry_env{fn}{i,1}, ...
%              time_to_berry_env{fn}{i,2});
%         hold on
%         plot(1:max_berries, time_to_berry_env{fn}{i,1},'-o');
%         hold on
        
%         dur_ret = simulate_berry_dur(fn,i,time_to_berry_env{fn}{i,1});
%         plot(1:max_berries, dur_ret);
        xlabel('Berries','Fontsize',11);
        ylabel('Time to the berry','Fontsize',11);
        drawnow
%         keyboard
        time_avg_trial = nan(1,num_trials);
        
        for b = 2:max_berries
            time_avg_trial(data_subj{fn}.berries == b) = ...
                diff_time_env{fn}{i,1}(b)*1000;
        end
        
        norm_GUT(i,:) = data_subj{fn}.giving_up_time./time_avg_trial;
%         norm_GUT(i,isinf(norm_GUT(i,:))) = NaN;
        % create a new field
        data_subj{fn}.gut_norm = norm_GUT(i,:);
        fprintf('Saving %s summary data....\n',subj_name);
        save(fullfile(dir_file,'matfiles_data', subj_name,'data_subj.mat'),'data_subj');
    end
    
    
%     keyboard
    avg_per_subj_gut(:,fn) = nanmean(norm_GUT,2);
    
%     simulate_berry_dur(fn,18,'Fontsize',15);
%     simulate_berry_dur(fn,11);

%     if fn==1
%         ttl = 'Low Effort Environment';
%     else 
%         ttl = 'High Effort Environment';
%     end
%     

%     subplot(2,1,1);
%     title('\beta = 0.5','Interpreter','Tex');
%     env_time_to_berry_means = cell2mat({time_to_berry_env{fn}{1:10,1}}');
%     env_mean = nanmean(env_time_to_berry_means,1);
%     marker_ = strcat('-o',clrs{fn});
%     plot(1:max_berries, env_mean, marker_, 'LineWidth',2);
%     hold on
% 
%     subplot(2,1,2);
%    
%     title('\beta = 0.2','Interpreter','Tex');
% 
%     env_time_to_berry_means = cell2mat({time_to_berry_env{fn}{11:18,1}}');
%     env_mean = nanmean(env_time_to_berry_means,1);
%     marker_ = strcat('-o',clrs{fn});
%     plot(1:max_berries, env_mean, marker_, 'LineWidth',2);
%     
% %     sgtitle(ttl,'Fontsize',15);
%     hold on
end

bar(avg_per_subj_gut);
figure
subj_avg_gut = nanmean(avg_per_subj_gut,1);
subj_sd_gut = nanstd(avg_per_subj_gut,0,1);
bar(subj_avg_gut,'grouped');
hold on
errorbar(subj_avg_gut,subj_sd_gut);
keyboard
tt_subj_avg = diff(avg_per_subj_gut,1,2);
keyboard


function ret_dur= simulate_berry_dur(fn, i, t, alpha, beta)

%max berries

if nargin < 4
    % error('Need to provide time');
    % elseif nargin<2
    
    if i <=10
        g =1;
        beta = 0.9;
    else
        g=2;
        beta = .6;
    end
    alpha =  50;
end

dur = @(t) alpha*beta./((1+ beta*t).^2 );
% 
ret_dur = 1./dur(t);

fg = zeros(1,1000)+30; 
fg(end+1) = 9;
t_cur = 0;
t_next = t_cur;
berries =[];
fg_cur = fg(1);
ct=1;
t_h = [];
while fg_cur > 27
    if t_cur >= t_next 
        
        berries(end+1)= 1;
        t_h(end+1) = t_cur;
        t_next = 1/dur(t_cur);
    end
    t_cur = t_cur +10;
    ct = ct+1;
    fg_cur = fg(ct);
end 

% keyboard
% figure(fn);
% subplot(2,1,g);
% plot(cumsum(berries),t_h );


