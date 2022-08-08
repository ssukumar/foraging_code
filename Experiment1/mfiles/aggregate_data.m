function aggregate_data(exp_no)
% Function to aggregate main metrics/summary metrics across subjects and
% compile it in a table to be used by R


if nargin<1
    exp_no=1;
end

if exp_no == 1
    % pilot
    filename = 'subj_order.txt';
    table_name= 'foraging_reward_main.csv';
else
    filename = 'subj_order_e2.txt';
    table_name='foraging_reward_data.csv';
end
% folder = fullfile('..','pilot');
folder = '..';

filename = fullfile(folder, filename);
csv_filename = fullfile(folder, table_name);

subj_order = readtable(filename, 'Delimiter',' ',...
    'ReadVariableNames',0);
subj_order = table2cell(subj_order);
num_subj = size(subj_order,1);
num_trials = 200;
num_env=2;
%%
reqd_vars = {'ntrials','pk_vel', 'movement_duration','travel_duration',...
    'harvest_duration','pk_acc_rep','max_distance','avg_pk_fr','avg_hrt',...
    'grip_count','berries','alpha_vals'};


table_cols = {'Subj_id','Trial','Env','Order','Probe','pk_vel',...
    'movement_duration', 'travel_duration','harvest_duration','pk_acc_rep','max_distance',...
    'avg_pk_fr','avg_hrt','grip_count','berries','alpha_vals'};


% csv_arr = zeros(num_subj*num_env*num_trials,numel(table_cols));
csv_arr= [];
% create csv arr per subject per env and append it
m = 1;
for env =1:2
    data_all(env) = cell2struct(cell(length(...
        reqd_vars),1),reqd_vars);
end

if 1
    for s = 1:num_subj
        
        subj_name = subj_order{s,2};
        order = subj_order{s,3};
        load(fullfile(folder, 'matfiles_data',subj_name,...
            'data_subj.mat'));
        total_trials = data_subj{1}.ntrials + data_subj{2}.ntrials;
        subj_arr = zeros(total_trials, numel(table_cols));
        
        for env = 1:2
            prv_num_trials = num_trials;
            num_trials = data_subj{env}.ntrials;
            
            subj_id_ix =  strfind(table_cols,'Subj_id');
            % column index for subj_id
            subj_id_ix = not(cellfun('isempty',subj_id_ix));
            rows_ix_subj = (env-1)*prv_num_trials + (1:num_trials);
            subj_arr(rows_ix_subj, subj_id_ix) = zeros(length(rows_ix_subj),1)...
                + s ;
            trial_ix = strfind(table_cols, 'Trial');
            trial_ix = not(cellfun('isempty',trial_ix));
            subj_arr(rows_ix_subj, trial_ix) = [1:num_trials]';
            
            env_ix = strfind(table_cols,'Env');
            env_ix = not(cellfun('isempty', env_ix));
            subj_arr(rows_ix_subj, env_ix) = zeros(length(rows_ix_subj),1)...
                + (env - 1) ;
            
            order_ix = strfind(table_cols,'Order');
            order_ix = not(cellfun('isempty', order_ix));
            subj_arr(rows_ix_subj, order_ix) = zeros(length(rows_ix_subj),1)...
                + order ;
            
            if s ==3 && env==2
                keyboard
            end
            probe_ix = strfind(table_cols, 'Probe');
            probe_ix = not(cellfun('isempty', probe_ix));
            subj_arr(rows_ix_subj, probe_ix) = zeros(length(rows_ix_subj),1);
            subj_arr((env-1)*prv_num_trials + ...
                data_subj{env}.probe_indices, probe_ix) = 1;
            
            for k = 1:numel(reqd_vars)
                tab_ix = strfind(table_cols,reqd_vars{k});
                if strcmp(reqd_vars{k},'alpha_vals')
                    fprintf('%s \n',subj_name)
%                     keyboard;
                end
                tab_ix = not(cellfun('isempty',tab_ix));
                num_trials = length(data_subj{env}.(reqd_vars{k}));
                %             row_ix_add = (env-1)*prv_num_trials + (1:num_trials);
                subj_arr(rows_ix_subj ,tab_ix) = data_subj{env}.(reqd_vars{k});
                data_all(env).(reqd_vars{k}) = [data_all(env).(reqd_vars{k}),...
                    nanmean(data_subj{env}.(reqd_vars{k}))];
            end
            
        end
        csv_arr = [csv_arr; subj_arr];
    end
    
    % Creating table in long format in case analysis needs to be done in R
    table_R = array2table(csv_arr, 'VariableNames', table_cols);
    writetable(table_R, csv_filename, 'WriteVariableNames',1);
end

keyboard
%% Quantifying harvest vigor per subject
% 1. Looking at the trend in pk force rate over the course of the trial
% 2. Looking at the hazard rate (probability of leaving) after each berry

for env =1:2
    data_pfr(env) = cell2struct(cell(num_subj,1), [subj_order(:,2)']);
    data_prb_pfr(env) = cell2struct(cell(num_subj,1), [subj_order(:,2)']);
end
keyboard

clrs={'b-','r-'};
ct = 1; 


for s = 1:num_subj
    
    subj_name = subj_order{s,2};
    order = subj_order{s,3};
    load(fullfile(folder, 'matfiles_data',subj_name,...
        'data_subj.mat'));
    
    
    for env = 1:2
        % 1 is low reward and 2 is high reward]
        pk_rate_vals_arr =[];
        probe_ix = strfind(table_cols, 'Probe');
        probe_ix = not(cellfun('isempty', probe_ix));
        for tr = 1:data_subj{env}.ntrials-1
            fprintf('Subj : %d; Env:%d, Trial: %d\n', s, env, tr); 
            
            if isempty(pk_rate_vals_arr)
                pk_rate_vals_arr = data_subj{env}.peak_force_rate{tr};
            else
                
                if length(data_subj{env}.peak_force_rate{tr}) >=...
                        size(pk_rate_vals_arr,2)
                    pk_rate_vals_arr = [pk_rate_vals_arr,...
                        nan(size(pk_rate_vals_arr,1),...
                        length(data_subj{env}.peak_force_rate{tr}) -...
                        size(pk_rate_vals_arr,2))];
                    pk_rate_vals_arr = [pk_rate_vals_arr;...
                        data_subj{env}.peak_force_rate{tr}];
                else
                    pk_rate_vals_arr = [pk_rate_vals_arr;...
                        [data_subj{env}.peak_force_rate{tr},...
                        nan(1,size(pk_rate_vals_arr,2) -...
                        length(data_subj{env}.peak_force_rate{tr}))]];
                end
            end
        end
        data_pfr(env).(subj_name)= nanmean(pk_rate_vals_arr, 1);
        subplot(2,2,ct);
        plot(data_pfr(env).(subj_name),clrs{env});
        hold on
        errorbar(data_pfr(env).(subj_name), nanstd(pk_rate_vals_arr,0, 1),...
            clrs{env});
        title(sprintf('Subj # %d',s), 'Fontsize',14);
        xlabel('Grip Count #');
        ylabel('Peak Force Rate (N/s)');
        hold on
    end
    ct = ct+1; 
    ct = mod(ct,5);
    if ct == 0 
        ct =1;
        figure;
    end
end

