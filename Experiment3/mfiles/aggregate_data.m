function aggregate_data(exp_no)
% Function to aggregate main metrics/summary metrics across subjects and
% compile it in a table to be used by R


if nargin<1
    exp_no=1;
end

if exp_no == 1
    filename = 'subj_order.txt';
    table_name= 'table_R.csv';
    table_new_name= 'table_R_grip.csv';
    struct_mat_file = 'forRazan_grip_shruthi.mat';
else 
    filename = 'subj_order_e2.txt';
    table_name='table_R_e2.csv';
end

folder = '..';

filename = fullfile(folder, filename);
csv_filename = fullfile(folder, table_name);
csv_new_filename = fullfile(folder, table_new_name);
subj_order = readtable(filename, 'Delimiter',' ',...
    'ReadVariableNames',0);
subj_order = table2cell(subj_order);
num_subj = size(subj_order,1);
num_trials = 200;    
num_env=2;
reqd_vars = {'pk_vel', 'movement_duration','travel_duration',...
        'max_distance','harvest_duration','patch_residence_duration',...
        'force_rate_max', 'force_per_berry','berries','grip_ramp_up_dur'...
        'grip_ramp_up2','harvest_reaction_time', 'grip_drop_time','grip_release_time',...
        'avg_rest_dur','avg_grip_bias','avg_grip_variance',...
        'grip_overshoot','time_to_last_berry','grip_drop_last_berry',...
        'giving_up_time','gut_norm','norm_dist','norm_pv'};
    
    
table_cols = {'Subj_id','Trial','Env','Order','Probe','pk_vel',...
    'movement_duration', 'travel_duration','max_distance',...
    'harvest_duration','patch_residence_duration', 'force_rate_max',...
    'force_per_berry','berries','grip_ramp_up_dur','grip_ramp_up2',...
    'harvest_reaction_time', 'grip_drop_time', 'grip_release_time',...
    'avg_rest_dur','avg_grip_bias', 'avg_grip_variance','grip_overshoot',...
    'time_to_last_berry','grip_drop_last_berry', 'giving_up_time',...
    'gut_norm','norm_dist','norm_pv'};

table_new_cols = {'Exp','Subj_id','Env','Probe','Trial', 'target_force'};
forRazan = struct([]);
csv_arr = zeros(num_subj*num_env*num_trials,numel(table_cols));
duration_mat_avg = cell(1,2);
m = 1;
csv_new_arr = zeros(num_subj*num_env*num_trials, numel(table_new_cols));
force_col = cell(num_subj*num_env*num_trials,1);
vel_col = cell(num_subj*num_env*num_trials,1);
% temporary code for probability matrix
frac_begin = round(0:0.1:0.9,1,'decimals');
frac_duration = round(0:0.01:0.39,2,'decimals');

csv_new_arr(:,1) = ones(num_subj*num_env*num_trials,1);
csv_new_arr(:,6) = zeros(num_subj*num_env*num_trials,1)+30;

for fn = 1:num_env
%     csv_arr = zeros(num_subj*num_env*num_trials,numel(table_cols));
    data_all(fn) = cell2struct(repmat({nan(num_subj,num_trials)},length(...
        reqd_vars),1),reqd_vars);
    
    duration_mat_avg{fn} = zeros(40,10);
%     csv_new_arr{};
    env = {'Low Effort', 'High Effort'};
    env_ix = strfind(table_cols, 'Env');
    env_ix = not(cellfun('isempty',env_ix));
    rows_ix_env = (fn-1) * num_trials*num_subj + (1:(num_subj*num_trials));
    csv_arr(rows_ix_env, env_ix) = zeros(length(rows_ix_env),1) + fn-1 ;
    env_new_ix = strfind(table_new_cols, 'Env');
    env_new_ix = not(cellfun('isempty',env_new_ix));
    csv_new_arr(rows_ix_env,env_new_ix) = zeros(length(rows_ix_env),1)...
            + fn-1 ;
    subj_ct = 1;
    for i = [cell2mat(subj_order(:,1))]'
        
        subj_name = subj_order{subj_ct,2};
        order = subj_order{subj_ct,3};
        load(fullfile(folder, 'matfiles_data',subj_name,...
            'data_subj.mat'));
        % Aggregate peak velocity
        
        duration_mat_avg{fn} = duration_mat_avg{fn} + data_subj{fn}.duration_matrix;
       
        subj_id_ix =  strfind(table_cols,'Subj_id');
        subj_id_ix = not(cellfun('isempty',subj_id_ix));
        rows_ix_subj = (fn-1)*num_trials*num_subj + (subj_ct-1)*num_trials + ...
            (1:num_trials);
        csv_arr(rows_ix_subj, subj_id_ix) = zeros(length(rows_ix_subj),1)...
            + i ;
        
        subj_id_new_ix =  strfind(table_new_cols,'Subj_id');
        subj_id_new_ix = not(cellfun('isempty',subj_id_new_ix));
        csv_new_arr(rows_ix_subj, subj_id_new_ix) = zeros(length(rows_ix_subj),1)...
            + i ;
        
        trial_ix = strfind(table_cols, 'Trial');
        trial_ix = not(cellfun('isempty',trial_ix));
        csv_arr(rows_ix_subj, trial_ix) = [1:num_trials]';
        trial_new_ix = strfind(table_new_cols, 'Trial');
        trial_new_ix = not(cellfun('isempty',trial_new_ix));
        csv_new_arr(rows_ix_subj, trial_new_ix) = [1:num_trials]';
        
        order_ix = strfind(table_cols,'Order');
        order_ix = not(cellfun('isempty', order_ix));
        csv_arr(rows_ix_subj, order_ix) = zeros(length(rows_ix_subj),1)...
            + order ;
%         csv_new_arr(rows_ix_subj, order_ix) = zeros(length(rows_ix_subj),1)...
%             + order ;
        
        probe_ix = strfind(table_cols, 'Probe');
        probe_ix = not(cellfun('isempty', probe_ix));
        probe_indices = data_subj{fn}.probe_trials;
        rows_ix_probe = (fn-1) * num_trials*num_subj + (subj_ct-1)*num_trials+...
            probe_indices; 
        csv_arr(rows_ix_probe, probe_ix) = ones(length(rows_ix_probe),1);
        probe_new_ix = strfind(table_new_cols, 'Probe');
        probe_new_ix = not(cellfun('isempty', probe_new_ix));
        csv_new_arr(rows_ix_probe, probe_new_ix) = ones(length(rows_ix_probe),1);
        
        for k = 1:numel(reqd_vars)
            tab_ix = strfind(table_cols,reqd_vars{k});
            tab_ix = not(cellfun('isempty',tab_ix));
            num_trials = length(data_subj{fn}.(reqd_vars{k}));
            row_ix_add = (fn-1)*num_trials*num_subj + (subj_ct-1)*num_trials; 
            csv_arr(row_ix_add + (1:num_trials) ,tab_ix) = ...
                data_subj{fn}.(reqd_vars{k});
            data_all(fn).(reqd_vars{k})(subj_ct,:) = data_subj{fn}.(reqd_vars{k});
        end
        
        force_col(row_ix_add + (1:num_trials)) = ...
            data_subj{fn}.grip_force';
        vel_col(row_ix_add + (1:num_trials)) = ...
            data_subj{fn}.vel_profile_rep';
        
        if ~isempty(data_subj{fn}.scatter_mat)
% %             figure(fn);
% %             hold on
% %             scatter(data_subj{fn}.scatter_mat(:,1),data_subj{fn}.scatter_mat(:,2));
% %             xlabel('Fraction of grip duration at which dip/rest begins');
% %             ylabel('Duration of rest as a fraction of grip duration');
% %             title(sprintf('Environment : %s', env{fn}));
        end
        subj_ct = subj_ct+1;
        
    end
    
%     duration_mat_avg{fn} = duration_mat_avg{fn}./length([cell2mat(subj_order(:,1))]);
%     figure(m);
%     imagesc(frac_begin, frac_duration, duration_mat_avg{fn});
%     colorbar;
%     ylabel('Duration of rest period as a fraction of total');
%     xlabel('Fraction of grip duration at which rest began');
%     title(' Probability of rest occuring ');
%     keyboard

    
    m=fn+1;
end

keyboard
% compute_hazard_probability(num_env, num_trials, num_subj, data_all)

% Creating table in long format in case analysis needs to be done in R
table_R = array2table(csv_arr, 'VariableNames', table_cols);
writetable(table_R, csv_filename, 'WriteVariableNames',1);

% table_new_R = array2table(csv_new_arr, 'VariableNames', table_new_cols);
% table_new_R.grip_profile = force_col;
% table_new_R.vel_profile = vel_col;
% forRazan = table2struct(table_new_R); 
% writetable(table_new_R, csv_new_filename, 'WriteVariableNames',1);
% save(struct_mat_file,'forRazan');
% Correlation between travel duration and peak velocity
clrs = {'b','r'};
% m=2;
if (1)
    m_ix = m;
    for env = 1:2
        for m = 1:3
            figure(m+m_ix);
            xx='Peak Velocity';
            x = data_all(env).pk_vel(:,2:200);
            x=x(:);
            if m== 1
                y = data_all(env).travel_duration(:,2:200);
                y=y(:);
                ind = ~isnan(x) & ~isnan(y);
                y=y(ind);
                yy= 'travel duration';
                
            elseif m==2
                y=  data_all(env).movement_duration(:,2:200);
                y=y(:);
                ind = ~isnan(x) & ~isnan(y);
                y=y(ind);
                yy='movement duration';
            else
                x= data_all(env).force_rate_max(:,2:200);
                x=x(:);
                xx='PeakF Force Rate';
                y= data_all(env).grip_ramp_up_dur(:,2:200);
                y=y(:);
                ind = ~isnan(x) & ~isnan(y);
                y=y(ind);
                yy='Grip Ramp Up Duraition';
            end
            
            
            x = x(ind);
            scatter(x,y, clrs{env});
            [coeffs,S] = polyfit(x', y', 1);
            fittedX = linspace(min(x), max(x), length(x));
            fittedY = polyval(coeffs, fittedX);
            r_squared = 1 - (S.normr/norm(y - mean(y)))^2;
            fprintf('\n R squared for %s correlated with %s : %d\n', yy, xx,...
                r_squared);
            % Plot the fitted line
            hold on;
            plot(fittedX, fittedY, 'r-', 'LineWidth', 3);
            ylabel(yy);
            xlabel(xx);
            title ('Correlation');
        end
    end
    m = m_ix+m;
end

% Correlation between "harvest vigor" 
if (1)
    figure(m);
    for env = 1:num_env
        xx = 'Berries';
        x = data_all(env).berries;
        x = x(:);
        
        yy = 'Harvest Vigor/ Max force rate';
        y = data_all(env).force_rate_max;
        y = y(:);
        ind = ~isnan(x) & ~isnan(y);
        y = y(ind);
        x = x(ind);
        scatter(x,y, clrs{env});
        [coeffs,S] = polyfit(x', y', 1);
        fittedX = linspace(min(x), max(x), length(x));
        fittedY = polyval(coeffs, fittedX);
        r_squared = 1 - (S.normr/norm(y - mean(y)))^2;
        fprintf('\n R squared for %s correlated with %s : %d\n', yy, xx,...
            r_squared);
        % Plot the fitted line
        hold on;
        plot(fittedX, fittedY, 'r-', 'LineWidth', 3);
        ylabel(yy);
        xlabel(xx);
        title ('Correlation');
    end
    m=m+1;
end

%plot the population means for different metrics
for env = 1:num_env
    k_ix = m+1;
    m_ix = 1; % subplot tile number

    for k = 1:numel(reqd_vars)
        figure(k_ix)
        subplot(2,2,m_ix)
        plot_dat = data_all(env).(reqd_vars{k});
        plot_mean = nanmean(plot_dat,1);
        
        plot_ster = nanstd(plot_dat,0,1)./sqrt(num_subj);
        if strcmp(reqd_vars{k},'travel_duration') || strcmp(reqd_vars{k},'max_distance')
           xdata=[2:num_trials];
           plot_mean=plot_mean(2:num_trials);
           plot_ster = plot_ster(2:num_trials);
        else
            xdata = [1:num_trials];
        end
        plotsdshading(xdata, plot_mean', plot_ster', clrs{env});
        title(reqd_vars{k});
        hold on;
        if env==2
            patch([probe_indices(1) probe_indices(1) probe_indices(end/4)...
                probe_indices(end/4)] ,[0.75*min(plot_mean) 1.25*max(plot_mean)...
                1.25*max(plot_mean) 0.75*min(plot_mean)] , 'r');
            alpha(0.2)
            patch([probe_indices(end/4+1) probe_indices(end/4+1) probe_indices(end/2)...
                probe_indices(end/2)] ,[0.75*min(plot_mean) 1.25*max(plot_mean)...
                1.25*max(plot_mean) 0.75*min(plot_mean)] , 'r');
            alpha(0.2)
            patch([probe_indices(end/2+1) probe_indices(end/2+1) probe_indices(end/4)+100 ...
                probe_indices(end/4)+100] ,[0.75*min(plot_mean) 1.25*max(plot_mean)...
                1.25*max(plot_mean) 0.75*min(plot_mean)] , 'r');
            alpha(0.2)
            patch([probe_indices(end/4+1)+100 probe_indices(end/4+1)+100 probe_indices(end) ...
                probe_indices(end)] ,[0.75*min(plot_mean) 1.25*max(plot_mean)...
                1.25*max(plot_mean) 0.75*min(plot_mean)] , 'r');
            alpha(0.2)
        end
        m_ix = m_ix +1;
        
        if mod(m_ix, 4) == 1
            k_ix = k_ix + 1;
           m_ix = 1;
        end
        
    end
end

hold off
m = k_ix;
% Plot only probe trials metrics
for env = 1:num_env
    k_ix = m + 1;
    m_ix = 1;
    
    for k = 1:numel(reqd_vars)
        plot_dat = data_all(env).(reqd_vars{k});
        plot_dat = plot_dat(:,probe_indices);
        plot_mean = nanmean(plot_dat,1);
        figure(k_ix);
        subplot(2, 2, m_ix);
        plot_ster = nanstd(plot_dat, 0, 1)/sqrt(num_subj);
        xdata = [1:length(probe_indices)];
        plotsdshading(xdata, plot_mean', plot_ster', clrs{env});
        title(reqd_vars{k})
        hold on;

        m_ix = m_ix + 1;
        if mod(m_ix, 4) == 1 
            k_ix = k_ix + 1;
            m_ix = 1;
        end
        
    end
end

m = k_ix;

% Compute and plot hazard rate for berries and for time

function compute_hazard_probability(num_env, num_trials, num_subj, data_all)
% This function computes the hazard probability associated with each berry
% and each binned time step

% [max_hazard_time, max_berries] =  compute_time_berries(alpha, beta, 0.8);
alpha = 20;
frac = 0.7;
max_berries = frac*alpha;
hazard_matrix_tr = cell(1,num_env);
hazard_matrix_tr{1} = zeros(num_trials, max_berries);
hazard_matrix_tr{2} = zeros(num_trials, max_berries);

hazard_prob = zeros(1,max_berries);
env_name= {'Low eff', 'high eff'};

for env = 1:num_env
    for tr = 1:num_trials
        max_stay = max(data_all(env).berries(:,tr));
        for b = 1:max_berries
            if b <= max_stay
                hazard_matrix_tr{env}(tr, b) = sum(data_all(env).berries(:,tr)...
                    == b)/sum(data_all(env).berries(:,tr)>=b);
            else
                hazard_matrix_tr{env}(tr,b) = Inf;
            end
        end
    end
    subplot(1,2,env);
    imagesc(hazard_matrix_tr{env});
    xlabel('Berries');
    ylabel('Trial');
    title(env_name{env})
end

figure
hazard_matrix_subj= cell(1,num_env);
hazard_matrix_subj{1}= zeros(num_subj, max_berries);
hazard_matrix_subj{2}= zeros(num_subj, max_berries);
for env = 1:num_env
    for s = 1:num_subj
        max_stay = max(data_all(env).berries(s,:));
        for b = 1:max_berries
%             max_stay = max(data_all(env).berries(s,:));
            if b <= max_stay
                hazard_matrix_subj{env}(s, b) = sum(data_all(env).berries(s,:)...
                    == b)/sum(data_all(env).berries(s,:) >= b);
            else
                hazard_matrix_subj{env}(s, b) = Inf;
            end
        end
    end
    subplot(1,2,env)
    imagesc(hazard_matrix_subj{env});
    xlabel('Berries');
    ylabel('Subj');
    title(env_name{env})
end
