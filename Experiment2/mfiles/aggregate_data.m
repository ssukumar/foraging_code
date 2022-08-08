function aggregate_data(exp_no)
% Function to aggregate main metrics/summary metrics across subjects and
% compile it in a table to be used by R

if nargin<1
    exp_no=1;
end

if exp_no == 1
    % pilot
    filename = 'subj_order_new.txt';
    table_name= 'foraging_reward_new.csv';
% else
%     filename = 'subj_order_e2.txt';
%     table_name='foraging_reward_data.csv';
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
    'alpha_vals','grip_count','berries','mid_15_trial',...
    'time_bw_patch_harvests', 'time_last_berry'};


table_cols = {'Subj_id','Trial','Env','Order','Probe','pk_vel',...
    'movement_duration', 'travel_duration','harvest_duration','pk_acc_rep','max_distance',...
    'avg_pk_fr','avg_hrt','alpha_vals','grip_count','berries','mid_15_trial',...
    'time_bw_patch_harvests', 'time_last_berry'};


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
        if isnumeric(subj_order{s,2})
            subj_num = subj_order{s,2}; 
            subj_name = num2str(subj_order{s,2});
        else
            subj_name = subj_order{s,2};
        end
        order = subj_order{s,3};
        load(fullfile(folder, 'matfiles_data',subj_name,...
            'data_subj.mat'));
        total_trials = data_subj{1}.ntrials + data_subj{2}.ntrials; 
        % variable number of trials due to fixed time
        
        subj_arr = zeros(total_trials, numel(table_cols));
        
        for env = 1:2
            prv_num_trials = num_trials; % using deafault val of 200 but multiplied by zero for first env
            num_trials = data_subj{env}.ntrials;
            
            subj_id_ix =  strfind(table_cols,'Subj_id');
            % column index for subj_id
            subj_id_ix = not(cellfun('isempty',subj_id_ix));
            rows_ix_subj = (env-1)*prv_num_trials + (1:num_trials);
            subj_arr(rows_ix_subj, subj_id_ix) = zeros(length(rows_ix_subj),1)...
                + subj_num ;
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
            

            probe_ix = strfind(table_cols, 'Probe');
            probe_ix = not(cellfun('isempty', probe_ix));
            subj_arr(rows_ix_subj, probe_ix) = zeros(length(rows_ix_subj),1);
            subj_arr((env-1)*prv_num_trials + ...
                data_subj{env}.probe_indices, probe_ix) = 1;
            
            for k = 1:numel(reqd_vars)
                tab_ix = strfind(table_cols,reqd_vars{k});
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



