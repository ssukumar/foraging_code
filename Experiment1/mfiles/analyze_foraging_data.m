function analyze_foraging_data(pilot, subj, exp, test)
% Code that is created to analyze the data for all subjects that have done
% the two foraging experiments designed by Shruthi Sukumar
% Code accepts 3-letter subject code.

% Make sure all residual graphs are closed
close all;
clear all; 
if nargin<1 
    pilot = 0;
    filename_low = 'foraging_DMM1_low_rewd1';
    filename_high = 'foraging_DMM_high_rewd1';
    filenames = {filename_low, filename_high};
    subj = 'DMM';
    exp =1;
    test = 1;
elseif pilot == 1
    filename_low = input('\nEnter the name of the filename for low effort environment: ');
    filename_high = input('\nEnter the name of the filename for high effort environment: ');
    filenames = {filename_low, filename_high};
    if nargin<2
        subj = input('\nEnter the subject identifier for this: ');
        exp = 1;
        test = 1;
    end
elseif pilot == 0
    
    if nargin<2
        subj = 'ZH1';
        exp = 1;
        test = 1;
    end
    
    conditions = {'low','high'};
    
    if exp == 1
        exp_name = 'foraging_';
    else
        exp_name = 'harvest_';
    end
    
    fprintf('Running throught %s`s data \n', subj);
    
    % Code for robot filename
    
    ret_subj = @(subj,exp_name,cond) ( strcat(exp_name,subj,'20_',cond,...
        num2str(test+1)));
    
    filenames = {ret_subj(subj,exp_name, conditions{1}), ret_subj(subj, ...
        exp_name,conditions{2}) };
end
%Data directory
f = pwd;
cd('..');
if pilot
    dir_file = fullfile(pwd,'pilot');
else
    dir_file = pwd;
end

if exist(fullfile(dir_file,'plots',subj)) ==0
    mkdir(fullfile(dir_file,'plots',subj));
end
cd(f);

% When fn == 1 the contents of the matfile are associated with the low
% effort environment in both experiments. This is because of how filenames
% are declared above.

global data_subj avg_profiles
clrs = {'b','r'};
fig_ix =1;

avg_profiles = cell(2,2);

for fn = 1:length(filenames)
%     clear total_score;
    % Put in checks to ensure that the environment order is right and that
    % the trials are probe/env according to the predetermined order
    % Read data from the existing robot data folder
    
    filename = filenames{fn};
    filedir = fullfile(dir_file,'robot_data', filename);
    %     ix_1 = length(exp_name)+1;
    subj_mat_filename = strcat(filename(1:end),'.mat');
    discard_mat_name = strcat(filename(1:end),'_discard', '.mat');
    subj_dir = fullfile(dir_file, 'matfiles_data',subj);
    if exp ==2
        subj_dir = fullfile(subj_dir, 'harvest');
    end
    if ~isdir(subj_dir)
        mkdir(subj_dir);
    end
    mat_filename = fullfile(subj_dir, subj_mat_filename);
    discard_mat_name = fullfile(subj_dir, discard_mat_name);
    
    logfnid = 15; % Log function ID used for robotconv.m
    f= pwd;
    
    if (exist(mat_filename, 'file') == 2)
        load(mat_filename);
        load(discard_mat_name);
    else
        T = robotdataread(filename,filedir, logfnid);
        discard_trials = [];
        save(discard_mat_name, 'discard_trials');
        save(mat_filename, 'T');
    end
    cd(f);
    
    ntrials = T.trial_arr'; % Trials that actually ocurred
    num_trials = T.trials;
    
    % Initialize all required variables: metrics of interest
    % scalar metrics where there is one value per metric per trial
    params_env = {'ntrials'};
    scalar_metrics = {'harvest_duration', 'travel_duration',...
        'movement_duration','pk_vel', 'pk_vel_ix', 'pk_vel_ix_rep',...
        'pk_acc','pk_acc_rep','pk_acc_ix_rep','berries','grip_count',...
        'max_distance','harvest_reaction_time','grip_overshoot',...
        'hazard_time','giving_up_time','avg_pk_fr','avg_hrt',...
        'avg_tbh','norm_pv', 'norm_td', 'norm_hd','alpha_vals'};

    % array metrics where each metric is an array of values per trial
    arr_metrics = {'vel_profile','vel_profile_rep','x_profile',...
        'x_profile_rep','time_ms','time_ms_rep', 'Probe','pk_force_rate',...
        'berry_diff','pk_force_rate_ix', 'grip_force','grip_begin_ix',...
        'harvest_reaction_times','time_bw_harvests'}; 
    
    % all metrics
    reqd_vars =  cat( 2,params_env,scalar_metrics, arr_metrics);
    cell_rep =cat(1, num_trials,repmat({nan(1,num_trials)},length(scalar_metrics),1),...
        repmat({cell(1,num_trials)},length(arr_metrics),1));
    
    % creating the structure per environment type (low-vs-high)
    data_subj{fn} = cell2struct(cell_rep ,reqd_vars);
    data_subj{fn}.probe_indices = find(T.flag_probe);
    probe_indices = data_subj{fn}.probe_indices(data_subj{fn}.probe_indices<=num_trials);
    data_subj{fn}.probe_indices = probe_indices;
    outtarget_time = 0;
    outtarget_ix_prev = 1; % for the first trial, it has to be the first index
    time_ms_prev = [];
    alpha_arr = T.alpha(ntrials)';
    data_subj{fn}.alpha_vals = alpha_arr;
%     keyboard;
    for tr = 1:num_trials
        
        if ~isempty(find(tr==ntrials, 1))
            
            data = T.framedata(tr);
            
            x= data.x;
            y= data.y;
            vx = data.vx;
            vy = data.vy;
            
            % Robot has a problem where the kinematic data spikes
            
            data1 = find_discont(x,y,vx,data,0,1);
            vvx = diff23f5(data1.x,1/200,10);
            vvy = diff23f5(data1.y, 1/200, 10);
            
            % We use the differentiation of the displacement vector
            data_subj{fn}.vel_profile{tr} = vvx(:,2);
            data_subj{fn}.x_profile{tr} = x;
            data_subj{fn}.y_profile{tr} = y;
            
            % harvest behavior variables
            alpha = T.alpha(tr);
            time_bw_harvests = T.time_bw_harvests(tr);
            score = data.score;
            time_ms = data.time_ms;
            data_subj{fn}.time_ms{tr} = time_ms;
            grip_force = data.grip_force;
            data_subj{fn}.grip_force{tr} = grip_force;
            % Getting the Indices for time stamps of when the entered and
            % left patches
            attarget_ix = find(data.statenumber==6, 1, 'first');
            outtarget_ix = find(data.statenumber==6, 1,'last');
            
            if isempty(find(tr == discard_trials))
                
%                 if discard_prev_flag
%                     
%                     x_pos = -1*sign(x(1))*x;
%                     attarget_ix = find(x_pos > 0.165 , 1 ,'first');
%                     discard_prev_flag = 0;
%                 end
                % #########Find harvest duration#######
                fl =0;
                for j = 1:length(data.frame)
                    if data.frame(j) ~= j
                        discard_trials(end+1) = tr;
                        fl = 1;
                        break;
                    end
                end
                
                if fl ==1 
                    time_ms_prev = time_ms(1:j-1);
                    continue;
                end
                
                
                fl_nb = harvest_behavior(fn, tr, data.flag_berry, score,...
                    grip_force, time_ms, time_bw_harvests,...
                    30, attarget_ix, outtarget_ix);
                if fl_nb == 1
                    discard_trials(end+1) = tr;
                end
%                 % determine if trials need to be discarded 
%                 if ~isempty(data_subj{fn}.grip_valid_ix{tr})
                    
                return_pk_vel(fn, tr,0, vvx, x, y, time_ms, time_ms_prev);

                get_mvt_duration(fn, tr, time_ms, attarget_ix,...
                    outtarget_time, outtarget_ix_prev);
%                 plot_diagnostics(1, fn, tr, vvx);
%                 else
%                     discard_trials = [discard_trials, tr];
%                     save(discard_mat_name, 'discard_trials');
%                     continue;
%                 end
                save(discard_mat_name, 'discard_trials');
                time_ms_prev = time_ms;
            else 
                discard_prev_flag = 1;
                for j = 1:length(data.frame)
                    if data.frame(j) ~= j
                        break;
                    end
                end
                time_ms_prev = time_ms(1:j-1);
%                 % determining what happens on those discarded trials
            end
        end
        
        % current val to be used in next trial
        outtarget_time = time_ms(outtarget_ix);
        outtarget_ix_prev = outtarget_ix;
        
        
    end
    
%     fig_ix = avg_vel_profiles(fn, data_subj{fn}.probe_indices, clrs{fn}, fig_ix);
%     keyboard

end
% fig_ix = berry_raster_plot(clrs, fig_ix);
% keyboard
fig_ix = plot_comparisons_env(fig_ix);

% plot_mvt_traces(fig_ix);
fprintf('Saving %s summary data....\n',subj);
save(fullfile(dir_file,'matfiles_data', subj,'data_subj.mat'),'data_subj');

save(fullfile(dir_file,'matfiles_data', subj,'vel_profile.mat'),'avg_profiles');


% Function to determine time between harvests - testing tcl timing
function fl_no_berry = harvest_behavior(fn, tr, flag_berry, score, grip_force, time_ms, ...
    grip_threshold, time_bw_harvests, attarget_ix, outtarget_ix)
% compute the time between collecting berry and next time the flag_berry
% gets set
global data_subj
persistent total_score; 

fl_no_berry = 0;
if score(end) == 0
    fl_no_berry = 1;
    return
end 

begin_ix = attarget_ix; 

harvest_ok_ix = begin_ix -1 + find(flag_berry(begin_ix:end),...
    1,'first');
first_har_ok_ix = harvest_ok_ix; 
harvest_no_ix = harvest_ok_ix - 1 + find(~flag_berry(harvest_ok_ix:end),1,...
    'first');
last_berry =0;
% variables to determine
grip_count = 0; 
berries = 0; 
% number of grip pulses resulting in a reward (for now this ignores unrewarded grips)
pk_rate_ix_arr = [];
pk_rate_vals = [];
harvest_reaction_times = [];
time_bw_=[]; 

force_rate = diff23f5(grip_force, 0.005, 10);
force_rate = force_rate(:,2);
force_rate_thresh = 10;

if tr == 1
    total_score = 0;
end

if isempty(total_score) 
    total_score = 0;
end
total_score= total_score + score(end);
harvest_ok_ixs = first_har_ok_ix;
curr_berry_ct = 0;
while ~last_berry
    time_bw_(end+1) = time_ms(harvest_ok_ix) - time_ms(begin_ix);
    [~,grip_pk_ix] = max(grip_force(harvest_ok_ix:harvest_no_ix));
    grip_pk_ix = grip_pk_ix + harvest_ok_ix - 1;
%     data_subj{fn}.grip_begin_ix{tr}(end+1) = grip_begin_ix;
    if score(harvest_no_ix) > curr_berry_ct
        
        curr_berry_ct = score(harvest_no_ix );
        grip_count = grip_count + 1;
    end
%     grip_count = grip_count + isempty(find(grip_force(harvest_ok_ix:end)...
%         >= grip_threshold,1,'first'));
%     force_rate_flip_ix = grip_pk_ix - 1 +...
%         find(force_rate(grip_pk_ix:end) <0 ,1,'first');
    
    if isempty(grip_pk_ix) 
        [pk_force_rate, pk_rate_ix] = max(force_rate(harvest_ok_ix:end));
    else
        [pk_force_rate, pk_rate_ix] = max(force_rate(harvest_ok_ix:grip_pk_ix));
    end
    if isempty(pk_force_rate)
        keyboard;
    end
    
    pk_rate_vals(end+1) = pk_force_rate;
%     if isempty(pk_ix)
%         keyboard;
%     end
    pk_rate_ix = harvest_ok_ix - 1 + pk_rate_ix;
    pk_rate_ix_arr(end+1) = pk_rate_ix;
    grip_begin_ix_flip = find(force_rate(pk_rate_ix:-1:harvest_ok_ix) <= ...
        0.05* pk_force_rate, 1, 'first');
    grip_begin_ix = length(force_rate(pk_rate_ix:-1:harvest_ok_ix)) - grip_begin_ix_flip;
    grip_begin_ix = grip_begin_ix + harvest_ok_ix - 1;
%     eta = 0.08;
%     while isempty(grip_begin_ix)
% %         keyboard
% 
%         if force_rate(harvest_ok_ix 
% 
%         grip_begin_ix_flip = find(force_rate(pk_rate_ix:-1:harvest_ok_ix) <= ...
%            	eta* pk_force_rate, 1, 'first');
%         grip_begin_ix = length(force_rate(pk_rate_ix:-1:harvest_ok_ix)) - grip_begin_ix_flip;
%         grip_begin_ix = grip_begin_ix + harvest_ok_ix - 1;
%         
%         eta = eta*1.25;
%     end
    if isempty(grip_begin_ix)
        grip_begin_ix = harvest_ok_ix;
    end

    data_subj{fn}.grip_begin_ix{tr}(end+1) = grip_begin_ix;
    harvest_reaction_times(end+1) = time_ms(grip_pk_ix) - ...
        time_ms(grip_begin_ix);
    begin_ix = harvest_no_ix; 
    harvest_ok_ix = begin_ix -1 + find(flag_berry(begin_ix:end),...
        1,'first');
    harvest_no_ix = find(~flag_berry(harvest_ok_ix:end),1,'first') +...
        harvest_ok_ix;

    if isempty(harvest_no_ix)
        harvest_no_ix = length(flag_berry);
    end
    if isempty(find(grip_force(harvest_ok_ix:end)>30)) || (score(harvest_ok_ix) == score(end))
        last_berry = 1;
        break
    end
%     ct = ct + 1;
    harvest_ok_ixs(end+1) = harvest_ok_ix;
end

data_subj{fn}.berries(tr) = curr_berry_ct; 
data_subj{fn}.grip_count(tr) = grip_count;
data_subj{fn}.peak_force_rate{tr} = pk_rate_vals;
data_subj{fn}.berry_diff{tr} = diff(score)';
data_subj{fn}.pk_force_rate_ix{tr} = pk_rate_ix_arr;
data_subj{fn}.harvest_reaction_times{tr} = harvest_reaction_times;
data_subj{fn}.avg_pk_fr(tr)= mean(pk_rate_vals);
data_subj{fn}.avg_hrt(tr) = mean(harvest_reaction_times);

data_subj{fn}.harvest_duration(tr) = time_ms(harvest_no_ix) - time_ms(first_har_ok_ix);

data_subj{fn}.leave_patch_ix(tr) = outtarget_ix; 
data_subj{fn}.enter_patch_ix(tr) = attarget_ix;
% if data_subj{fn}.avg_pk_fr(tr) <= 0 || data_subj{fn}.avg_hrt(tr)
%     keyboard;
% end
% clf
% plot(force_rate);
% hold on
% plot(pk_rate_ix_arr, force_rate(pk_rate_ix_arr),'b*');
% hold on
% yyaxis right;

time_samp_x = (1:length(grip_force))*0.005; % convert to seconds
% plot(time_samp_x, grip_force);
% yline(30,'--');
% xlabel('Time');
% 
% hold on 
% % hold on
% % plot(score)
% % keyboard
% 
% plot(time_samp_x(harvest_ok_ixs), grip_force(harvest_ok_ixs),'k+')
% 
% % plot(time_bw_, 'b*');
% % hold on
% % yline(time_bw_harvests,':');
% keyboard;
% hold off
data_subj{fn}.time_bw_harvests{tr} = time_bw_;
data_subj{fn}.avg_tbh(tr) = mean(time_bw_);
% determine the number of berries collected and the number of grip pulses


function fig_ix = avg_vel_profiles(fn, probe_indices, c, fig_ix)
% function to generate average grip force profiles

global data_subj avg_profiles

figure(fig_ix)

max_ix_mat = [data_subj{fn}.pk_vel_ix_rep]';
avg_vel_ = align_average(data_subj{fn}.vel_profile_rep, max_ix_mat,c,...
    1, probe_indices);
hold on

avg_profiles{fn,1} = avg_vel_{1};
avg_profiles{fn,2} = avg_vel_{2};
keyboard
fig_ix = fig_ix + 1;

function fig_ix = berry_raster_plot( clrs, fig_ix)
% function to generate average grip force profiles

global data_subj 

figure(fig_ix)
ber_raster = plot_raster_berries(data_subj,clrs);
hold on

keyboard
fig_ix = fig_ix + 1;

% function fig_ix = avg_grip_profiles(fn, probe_indices, c, fig_ix)
% % function to generate average grip force profiles
% 
% global data_subj
% 
% % begin_ix_mat = cell2mat(data_subj{fn}.grip_begin_ix);\
% 
% figure(fig_ix)
% 
% begin_ix_mat = [data_subj{fn}.grip_begin_ix]';
% align_average(data_subj{fn}.force_profile, begin_ix_mat,c,...
%     1, probe_indices);
% hold on
% 
% fig_ix = fig_ix + 1;
% 
% function fig_ix = avg_force_rate_profiles(fn, probe_indices, c, fig_ix)
% % function to generate average grip force profiles
% 
% global data_subj
% 
% % begin_ix_mat = cell2mat(data_subj{fn}.grip_begin_ix);\
% 
% figure(fig_ix)
% % for p =1:2
% % 
% %     if p ==1
% %         begin_ix_mat = [data_subj{fn}.grip_begin_ix(non_probe_indices)]';
% %         align_average(data_subj{fn}.grip_force(non_probe_indices), ...
% %             begin_ix_mat,c);
% %         hold on;
% %     else
% %         begin_ix_mat = [data_subj{fn}.grip_begin_ix(probe_indices)]';
% %         align_average(data_subj{fn}.grip_force(probe_indices), ...
% %             begin_ix_mat,'g');
% %         hold off;
% %     end
% %     
% % end
% 
% begin_ix_mat = [data_subj{fn}.force_rate_max_ix]';
% align_average(data_subj{fn}.force_rate, begin_ix_mat,c,...
%     1, probe_indices);
% hold on
% 
% fig_ix = fig_ix + 1;
% 
% function fig_ix = avg_vel_profiles(fn, probe_indices, c, fig_ix)
% % function to generate average grip force profiles
% 
% global data_subj
% 
% figure(fig_ix)
% 
% max_ix_mat = [data_subj{fn}.pk_vel_ix_rep]';
% align_average(data_subj{fn}.vel_profile_rep, max_ix_mat,c,...
%     1, probe_indices);
% hold on
% 
% fig_ix = fig_ix + 1;

function return_pk_vel(fn,tr,is_discard,vvx, x, y, time_ms, time_ms_prev)
% Returns index of peak velocity and the value
% Also repairs the velocity profiles : trial i vel profile is the movement
% to patch # i

global data_subj;

if data_subj{fn}.grip_count(tr) > 0
    start_grip_ix = data_subj{fn}.grip_begin_ix{tr}(1);
    [data_subj{fn}.pk_vel(tr),data_subj{fn}.pk_vel_ix(tr)] = ...
        max(abs(vvx(1:start_grip_ix,2)));
else
    [data_subj{fn}.pk_vel(tr),data_subj{fn}.pk_vel_ix(tr)] = ...
    max(abs(vvx(1:end,2)));
end


% Repair velocity profiles

if tr==1 
    data_subj{fn}.vel_profile_rep{tr} = vvx(:,2)';
    data_subj{fn}.x_profile_rep{tr} = x';
    data_subj{fn}.y_profile_rep{tr} = y';
    data_subj{fn}.pk_vel_ix_rep(tr) = data_subj{fn}.pk_vel_ix(tr);
    data_subj{fn}.time_ms_rep{tr} = time_ms';
    data_subj{fn}.profile_shift_rep(tr)=  0;
else
    if is_discard == 1
        profile_end_ix = length(data_subj{fn}.vel_profile{tr});
    else
        profile_end_ix = round(3*length(data_subj{fn}.vel_profile{tr})/4);
    end
    vel_prof_rep = [data_subj{fn}.vel_profile{tr-1}(end-round(length(...
        data_subj{fn}.vel_profile{tr-1})/3):end)' ...
        data_subj{fn}.vel_profile{tr}(1:profile_end_ix)' ];
    x_profile_rep = [data_subj{fn}.x_profile{tr-1}(end-round(length(...
        data_subj{fn}.x_profile{tr-1})/3):end)' ...
        data_subj{fn}.x_profile{tr}(1:profile_end_ix)' ];
    y_profile_rep = [data_subj{fn}.y_profile{tr-1}(end-round(length(...
        data_subj{fn}.y_profile{tr-1})/3):end)' ...
        data_subj{fn}.y_profile{tr}(1:profile_end_ix)' ];
    time_stitch = [time_ms_prev(end-round(length(data_subj{fn...
        }.vel_profile{tr-1})/3):end)',time_ms(1:profile_end_ix)'];
    t_miss = time_ms(1) - time_ms_prev(end);
    t_vec = linspace(time_ms_prev(end)+5, time_ms(1), t_miss/5 +1 );
    vel_bet = interp1(time_stitch, vel_prof_rep, t_vec);
    x_bet = interp1(time_stitch, x_profile_rep, t_vec);
    y_bet = interp1(time_stitch, y_profile_rep, t_vec);
    vel_prof_rep = [data_subj{fn}.vel_profile{tr-1}(end-round(length(...
        data_subj{fn}.vel_profile{tr-1})/2):end)', vel_bet, ...
        data_subj{fn}.vel_profile{tr}(1:profile_end_ix)' ];
    x_profile_rep =  [data_subj{fn}.x_profile{tr-1}(end-round(length(...
        data_subj{fn}.x_profile{tr-1})/2):end)', x_bet, ...
        data_subj{fn}.x_profile{tr}(1:profile_end_ix)' ];
    y_profile_rep =  [data_subj{fn}.y_profile{tr-1}(end-round(length(...
        data_subj{fn}.y_profile{tr-1})/2):end)', y_bet, ...
        data_subj{fn}.y_profile{tr}(1:profile_end_ix)' ];
    data_subj{fn}.vel_profile_rep{tr} = vel_prof_rep;
    data_subj{fn}.x_profile_rep{tr} = x_profile_rep;
    data_subj{fn}.y_profile_rep{tr} = y_profile_rep;
    data_subj{fn}.pk_vel_ix_rep(tr) = length([data_subj{fn}.vel_profile{tr-1}(...
        end-round(length(data_subj{fn}.vel_profile{tr-1})/2):end)',...
            vel_bet]) + data_subj{fn}.pk_vel_ix(tr) - 1;
    time_ms_rep =  [time_ms_prev(end-round(length(data_subj{fn...
        }.vel_profile{tr-1})/2):end)', t_vec, time_ms(1:profile_end_ix)'];
    data_subj{fn}.time_ms_rep{tr} = time_ms_rep;
    data_subj{fn}.profile_shift_rep(tr)=  length([data_subj{fn}.vel_profile{tr-1}(...
        end-round(length(data_subj{fn}.vel_profile{tr-1})/2):end)',...
            vel_bet]);
    
end



function get_mvt_duration(fn, tr, time_ms, attarget_ix, outtarget_time, ...
    outtarget_ix_prev)
% Computes the movement duration of the reach between two patches 
% This is based on onset time estimated from the previous trial and offset
% time estimated from current trial. Therefore, movement duration on trial
% i corresponds to movement to patch i.

global data_subj;

% First compute the travel durations, aka , the entire duration spent
% between when the subject left the last patch to when they entered the
% current patch

if tr ~=1
    if isempty( time_ms(attarget_ix) - outtarget_time) ||...
            isnan( time_ms(attarget_ix) - outtarget_time)
        keyboard;
    end
    data_subj{fn}.travel_duration(tr) = time_ms(attarget_ix) - outtarget_time;
else 
    data_subj{fn}.travel_duration(tr) = nan;
end

% Check if the travel duration computation is right

% Compute movement duration based on repaired velocity profiles and
% repaired time_ms matrix

% free parameters that can be changed 
sd_win = 10; % Window for computing standard deviation
vel_thresh = 0.005;
sd_thresh = 0.001;
vx = data_subj{fn}.vel_profile_rep{tr};
sd_vec = [];
for i =1:length(vx)-sd_win
    sd_vec(end+1) = std(vx(i:i+sd_win));
end

pk_vel_ix_rep = data_subj{fn}.pk_vel_ix_rep(tr);


if data_subj{fn}.vel_profile_rep{tr}(pk_vel_ix_rep) < 0 
    vel_profile_rep = -1 * data_subj{fn}.vel_profile_rep{tr};
else
    vel_profile_rep = data_subj{fn}.vel_profile_rep{tr};
end
acc_profile_rep = diff23f5(vel_profile_rep', 1/200, 10);
acc_profile_rep = acc_profile_rep(:,2);
[data_subj{fn}.pk_acc_rep(tr), data_subj{fn}.pk_acc_ix_rep(tr)] =...
        max(acc_profile_rep(1:pk_vel_ix_rep));
time_ms_rep = data_subj{fn}.time_ms_rep{tr};

vel_flipped = vel_profile_rep(end:-1:1);
acc_flipped = acc_profile_rep(end:-1:1);
pk_vel_ix_flipped = length(vel_flipped) - pk_vel_ix_rep + 1;

% determining indices on the repaired array (REMEMBER)
onset_ix_flipped = pk_vel_ix_flipped + find(abs(vel_flipped(...
    pk_vel_ix_flipped:end)) < vel_thresh,1,'first');


sd_flipped = sd_vec(end:-1:1); 
sd_drop_flipped = find(sd_flipped(pk_vel_ix_flipped:end)<sd_thresh, 1,...
    'first') + pk_vel_ix_flipped;
sd_ix_rep = length(sd_vec) - sd_drop_flipped + 1;

if onset_ix_flipped > sd_drop_flipped 
    onset_ix_rep = length(vel_flipped)- onset_ix_flipped + 1;
else
    onset_ix_rep = sd_ix_rep;
end

if isempty(onset_ix_rep)
   onset_ix_rep = find(vel_profile_rep(1:pk_vel_ix_rep) < 0 ,1,'last');
    
end

onset_time = time_ms_rep(onset_ix_rep);

vel0_ix_rep = pk_vel_ix_rep + find(vel_profile_rep(pk_vel_ix_rep:end) <0 ...
    , 1, 'first')-1;

sd0_ix_rep = pk_vel_ix_rep + find(sd_vec(pk_vel_ix_rep:end) < ...
    sd_thresh,1,'first')-1;

if ~isempty(vel0_ix_rep) && ~isempty(sd0_ix_rep)
    if sd0_ix_rep > vel0_ix_rep
        offset_ix_rep = sd0_ix_rep;
    else
        offset_ix_rep = vel0_ix_rep;
    end
else
    
    offset_ix_rep = pk_vel_ix_rep + find(vel_profile_rep(pk_vel_ix_rep:end)...
        == min(vel_profile_rep(pk_vel_ix_rep:end)),1,'first')-1;
end

% 
% grip_valid_ix_rep = data_subj{fn}.grip_valid_ix{tr} + ...
%     data_subj{fn}.profile_shift_rep(tr) -1;
% 
% if ~isempty(grip_valid_ix_rep)
%     if grip_valid_ix_rep(1) < offset_ix_rep
%         offset_ix_rep = grip_valid_ix_rep(1);
%     end
% end

offset_time = time_ms_rep(offset_ix_rep);
x_profile_rep = data_subj{fn}.x_profile_rep{tr};

data_subj{fn}.max_distance(tr) = max(x_profile_rep) - min(x_profile_rep);

if isempty(offset_time) || isempty(onset_time)
    data_subj{fn}.movement_duration(tr) = nan;
    keyboard;
else
    data_subj{fn}.movement_duration(tr) = offset_time - ...
        onset_time;
end
% Just for checking; remove once check complete and port to diagnostics fn
% h=figure(1)
% % WinOnTop(h,true)
% plot(data_subj{fn}.x_profile_rep{tr})
% hold on
% drawnow
% plot(vel_profile_rep);
% hold on
% plot(sd_vec,'k')
% hold on
% plot(onset_ix_rep, vel_profile_rep(onset_ix_rep),'r^');
% hold on
% plot(offset_ix_rep ,vel_profile_rep(offset_ix_rep),'g^');
% hold on
% plot(pk_vel_ix_rep ,vel_profile_rep(pk_vel_ix_rep),'k*');
% 
% yyaxis right
% 
% plot(acc_profile_rep)
% hold on 
% plot(data_subj{fn}.pk_acc_ix_rep(tr), acc_profile_rep(...
%     data_subj{fn}.pk_acc_ix_rep(tr)),'o');
%  
% hold off 
% keyboard;
% clf


function plot_diagnostics(diag_no,fn, tr, vvx)
% Function to check whether the computed values for peak velocity, harvest
% duration and movement duration make sense
% Diag no: 1- Pk Vel, 2- Velocity profile, 3- Harvest Duration
global data_subj
if nargin<1 
    diag_no=1;
    tr = input('Please Enter A Trial Number: ');
end
gcf;

if diag_no==1
    
    plot(vvx(:,2))
    hold on
    plot(data_subj{fn}.pk_vel_ix(tr),data_subj{fn}.pk_vel(tr),'k^')
    hold off
elseif diag_no==2
    
    plot(data_subj{fn}.vel_profile_rep{tr});
    hold on 
    plot(data_subj{fn}.pk_vel_ix_rep(tr), data_subj{fn}.vel_profile_rep{...
        tr}(data_subj{fn}.pk_vel_ix_rep(tr)),'k^');
    hold off
elseif diag_no==3
    
    plot(data_subj{fn}.vel_profile{tr});
    hold on
    if ~isempty(data_subj{fn}.grip_valid_ix{tr})
    
        plot(data_subj{fn}.grip_valid_ix{tr}, ...
            data_subj{fn}.vel_profile{tr}(data_subj{fn}.grip_valid_ix{tr})...
            ,'k^');
    end
    hold off
    
else
    fprintf('the other diagnostics are under construction ... ');
end
keyboard
% 
% function test_timing_diagnostics(metric)
% %function to determine the amount of noise present in prescribed times
% 
% global data_subj
% 



function plot_diag_env(fig_ix)
% Diagnostics at the environment level not the trial level

global data_subj 

if nargin < 1
    fig_ix = 1;
end
clrs = {'b','r'};

for i = 1:3
    
    for env = 1:2
        figure(fig_ix);
        % Correlations between movewment duration and peak velocity
        
        if i ==1 
            y =  data_subj{env}.movement_duration;
            x = data_subj{env}.pk_vel;
            yy = 'Movement Duration';
            xx= 'Peak Velocity';
        elseif i ==2
            y =  data_subj{env}.travel_duration;
            x = data_subj{env}.pk_vel;
            yy = 'Travel Duration';
            xx= 'Peak Velocity';
        else
            y =  data_subj{env}.movement_duration;
            x =  data_subj{env}.travel_duration;
            yy = 'Movement Duration';
            xx = 'Travel Duration';
        end
        
        scatter(x(2:end),y(2:end), clrs{env});
        [coeffs,S] = polyfit(x(2:end), y(2:end), 1);
        fittedX = linspace(min(x(2:end)), max(x(2:end)), 199);
        fittedY = polyval(coeffs, fittedX);
        r_squared = 1 - (S.normr/norm(y(2:end) - mean(y(2:end))))^2;
        fprintf('R squared for correlation of %s to %s - %d', yy, xx,...
            r_squared);
        % Plot the fitted line
        hold on;
        plot(fittedX, fittedY, 'r-', 'LineWidth', 3);
        ylabel(yy);
        xlabel(xx);
        title ('Correlation');
        hold on;
    end
    beautifyfig;
    hold off;
    fig_ix = fig_ix+1;
end


function fig_ix = plot_comparisons_env (fig_ix)
% Plotting within subject differences in the different metrics between
% environments

global data_subj

if nargin<1
    fig_ix = 1;
    
end

% Plot peak velocity differences 

clrs = {'b','r'}; 
metrics = {'pk_vel', 'movement_duration','travel_duration', 'pk_force_rate',...
    'harvest_reaction_time','time_bw_harvests','grip_count'};

m_ix = 1; % metric index
for metric = metrics
    
    figure(fig_ix);
    subplot(2,2,m_ix);
    max_val = -10000;
    min_val = 10000;
    
    for env = 1:2 % env=1 low effort; env = 2 is high effort
        if strcmp(metric,'pk_vel')
            plot_dat = data_subj{env}.pk_vel;
            tt = 'Pk Velocity';
        elseif strcmp(metric,'movement_duration')
            plot_dat = data_subj{env}.movement_duration;
            tt = 'Movement duration';
        elseif strcmp(metric,'travel_duration')
            plot_dat = data_subj{env}.travel_duration;
            tt = 'Travel duration';
        elseif strcmp(metric, 'time_bw_harvests')
            plot_dat = data_subj{env}.avg_tbh;
            tt = 'Time bw harvests';
        elseif strcmp(metric, 'pk_force_rate')
            plot_dat = data_subj{env}.avg_pk_fr;
            tt = 'Peak Force Rate';
        elseif strcmp(metric, 'harvest_reaction_time')
            plot_dat = data_subj{env}.avg_hrt;
            tt= 'Harvest Reaction Time';
        elseif strcmp(metric,'harvest_duration')
            plot_dat = data_subj{env}.harvest_duration;
            tt = 'Harvest duration';
        elseif strcmp(metric, 'patch_residence')
            plot_dat = data_subj{env}.patch_residence_duration;
            tt = 'Patch residence duration';
        elseif strcmp(metric, 'grip_count')
            plot_dat = data_subj{env}.grip_count;
            tt = 'Grip Count';
        elseif strcmp(metric,'avg_grip_bias')
            plot_dat = data_subj{env}.avg_grip_bias;
            tt='Average Grip Bias';
        elseif strcmp(metric,'avg_grip_variance')
            plot_dat = data_subj{env}.avg_grip_variance;
            tt='Average Grip Variance';
        end
        probe_indices = data_subj{env}.probe_indices;
        plot(1:length(plot_dat), plot_dat, clrs{env});
        hold on
        title(tt,'Fontsize',15);
%         legend('low','high');
        ylabel(tt);
        xlabel('Trial #');
        
        if max(plot_dat) > max_val
            max_val = max(plot_dat);
        end
        
        if min(plot_dat) < min_val
            min_val = min(plot_dat);
        end

    end
    gcf;
    
%     keyboard;
    
    d = diff(probe_indices); 
    idxs = [find(d~=1);length(probe_indices)]';
    num_probe_blocks = sum(double(d~=1)) + 1;
    beg_ix = 1;
    for i = 1:num_probe_blocks
        patch([probe_indices(beg_ix) probe_indices(beg_ix) probe_indices(idxs(i))...
            probe_indices(idxs(i))] ,[min_val max_val max_val min_val] , 'r');
        alpha(0.2)
        beg_ix = idxs(i)+1;
    end
    
    hold off;
    beautifyfig;
    m_ix = m_ix + 1;
    if mod(m_ix,4) == 1
        fig_ix = fig_ix +1;
        m_ix = 1;
    end
    
end

function fig_ix = plot_mvt_traces(fig_ix)
% Plotting the movement traces
global data_subj;
num_env = 2;
num_trials = 200;
figure(fig_ix);
clrs = {'b--','b-.', 'r--','r-.'};
clr_p = {'k-.'};
hAx=gobjects(2*2,1);        % preallocate for the axes handles

env_tt = {'Low Eff', 'High Eff'};
dir_tt = {'Leftward', 'Rightward'};

for e = 1:num_env
    probe_indices = data_subj{e}.probe_trials;

    e0 = e-1;
    for tr = 1:num_trials
        tr0 = 1 - mod(tr,2);
        m_ix = e0 * 2^1 + tr0 * 2^0 + 1;
        hAx(m_ix) = subplot(2,2,m_ix);
        clr = clrs{m_ix};
        lw=1;
        if ~isempty(find(tr==probe_indices))
            clr = clr_p{1};
            lw = 1.5;
        end
        plot(data_subj{e}.x_profile_rep{tr}, ...
            data_subj{e}.y_profile_rep{tr}, clr,'linewidth',lw)
        hold on
        title(strcat(env_tt{e},' /',dir_tt{tr0+1}));
    end
    
end

gca;

xlim(hAx,[-.25,.25]);
ylim(hAx, [-.1,.1]);

for m_ix = 1:4
    gca;
    subplot(2,2,m_ix);
    circle(-.2,0,.03);
    circle(.2,0,.03);
    xlabel( 'X coordinate');
    ylabel('Y coordinate');
end

beautifyfig;

function h = circle(x,y,r)

hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,'--');
hold off


function normalized_probe_vals(fn)
% Normalize the vector to the average of the mtric in non-probe trials
% Length of metric_norm is the same as probe_indices

global data_subj 

num_trials = 200;
probe_indices = data_subj{fn}.probe_trials;

non_probe_indices = setxor(1:num_trials, probe_indices);
avg_pv_np = nanmean(data_subj{fn}.pk_vel(non_probe_indices));
avg_dist_np = nanmean(data_subj{fn}.max_distance(non_probe_indices));
avg_td_np = nanmean(data_subj{fn}.travel_duration(non_probe_indices));
avg_hd_np = nanmean(data_subj{fn}.harvest_duration(non_probe_indices));

data_subj{fn}.norm_dist(probe_indices) = data_subj{fn}.max_distance(probe_indices)/avg_dist_np;
data_subj{fn}.norm_pv(probe_indices) = data_subj{fn}.pk_vel(probe_indices)/avg_pv_np;
data_subj{fn}.norm_td(probe_indices) = data_subj{fn}.travel_duration(probe_indices)/avg_td_np;
data_subj{fn}.norm_hd(probe_indices) = data_subj{fn}.harvest_duration(probe_indices)/avg_hd_np;


function [fig_ix, hazard_time] = compute_hazard_rate( fig_ix, exp_no)
% Computes the hazard rate averaged across patches in an environment
% grip_force_cell : normalized grip force value divided into two cells
% Cell index 1 : Low effort environment
% Cell index 2 : High effort environment

global data_subj
% 
% if size(grip_norm_cell)~=2
%     error('Size of variable grip_force_cell is not right. One cell is required per environment');
% end


hazard_time = cell(1,length(2));
max_time_hazard = zeros(1,2);
probe_indices = [21:30, 71:80, 121:130, 171:180];
num_drop = cell(1,2);

for env = 1:2
    
    harvest_time{env} = nan+zeros(1,200);
    
%     grip_norm = grip_norm_cell{env};
%     time_ms_ = time_ms_cell{env};

    grip_norm = data_subj{env}.grip_norm;
    time_ms_ = data_subj{env}.time_ms;
    grip0_ix = data_subj{env}.start_grip_ix;
    end_grip_ix = data_subj{env}.end_grip_ix;
    for n = 1:length(grip_norm)
        drop_last =0;
        grip_arr = grip_norm{n};
        time_arr = time_ms_{n};
        
        % determine normalized threshold based on required force
        if exp_no == 1
            thresh = 0.93;
        elseif exp_no ==2 
            if isempty( find(n==probe_indices)) 
                if env==1 
                    thresh = 0.8;
                else 
                    thresh = 0.96;
                end
            else
                thresh = 0.93;
            end
        end
        
        % count 
        num_drop{env}(end+1) = 0;
        
        % Check that counter size updates appropriately
        if length(num_drop{env}) ~= n
            fprintf(' Count size not updated correctly');
            keyboard;
        end
        
        % First point when grip force increases above threshold
%         ix0_har = find(grip_arr >= thresh,1,'first');
        ix0_har = grip0_ix(n);
        
        
        % Find subsequent point when grip force drops below required
        if ~isempty (ix0_har) && ~isnan(ix0_har)
            
            t0_har = time_arr(ix0_har);
            ix_drop = ix0_har - 1 + find(grip_arr(ix0_har:end) < thresh...
                ,1,'first');
            
            if ~isempty(ix_drop)
                ix_drop_curr = ix_drop;
                while (drop_last == 0)
                    
                    ix_inc_next = ix_drop_curr  + find(grip_arr(ix_drop_curr:end) > thresh...
                        ,1,'first');
                    if ~isempty(ix_inc_next)
                        num_drop{env}(n) = num_drop{env}(n)+1;
                        ix_drop_curr=  ix_inc_next + find(grip_arr(ix_inc_next:end) < thresh...
                            ,1,'first');
                    else
                        drop_last = 1;
                    end
                end
                
%                 num_drop{env}(n)
%                 keyboard
                if (0) %num_drop{env}(n) ~=0
                    figure
                    plot(grip_arr);
                    keyboard;
                end
                time_drop = time_arr(ix_drop);
                time_hazard = time_drop- t0_har;
                hazard_time{env}(n) = time_hazard;
                if time_hazard > max_time_hazard(env)
                    max_time_hazard(env)=time_hazard;
                end
            else
                hazard_time{env}(n) = nan;
            end
        else
            hazard_time{env}(n)=nan;
        end
        
    end
    
end

max([num_drop{1},num_drop{2}])
% keyboard

prob_hazards = cell(1,2);
time_bin_cell = cell(1,2);
n_bins=40;
clrs ={'b','r'};
for env = 1:2
    figure(fig_ix)

    [n,time_bins] = hist(hazard_time{env}, n_bins);    
    denom = cumsum(n,2, 'reverse');
    prob_hazards{env} = n./denom ;
    time_bin_cell{env} = time_bins;
    plot(time_bins,prob_hazards{env},strcat(clrs{env},'-o') );
    
    xlabel('Time grip dropped');
    ylabel('Hazard probability');
    
%     keyboard;
    hold on
end
fig_ix = fig_ix + 1;


function fig_ix = grip_drop_dyn (fn, discard_trials, exp, fig_ix)
% Function to plot when grip force drops below the required force on
% average and for how long

global data_subj

num_trials = length(data_subj{fn}.grip_force);
grip_cells = data_subj{fn}.grip_force_cell; 
probe_indices = data_subj{fn}.probe_trials;

grip_reqd = zeros(1, num_trials);
if exp == 1 
    grip_reqd(:) = 30;
    
else
    if fn ==1
        grip_reqd(:) = 10;
    else
        grip_reqd(:) = 50;
    end
    
    grip_reqd(probe_indices) = 30;
end

grip_reqd= grip_reqd - 2; % tolerance

frac_begin = round(0:0.1:0.9,1,'decimals');
frac_duration = round(0:0.01:0.39,2,'decimals');

duration_matrix = zeros(length(frac_duration), length(frac_begin));
duration_matrix_ct = zeros(length(frac_duration), length(frac_begin));

scatter_mat = [];

duration_cell = cell(1,num_trials);
grip_fracs = zeros(1, num_trials);
idx_fracs = cell(1,num_trials);
ixs_req = cell(length(frac_begin),length(idx_fracs));

for tr = 1:num_trials

    if isempty(find(tr == discard_trials))
        time_ms = data_subj{fn}.time_ms{tr};
        start_grip_ix = data_subj{fn}.start_grip_ix(tr);
        end_grip_ix = data_subj{fn}.end_grip_ix(tr);
        if isnan(end_grip_ix)
            discard_trials = [discard_trials, tr];
        elseif start_grip_ix == end_grip_ix
            discard_trials = [discard_trials, tr];
        else
            grip_fracs(tr) = time_ms(end_grip_ix) - time_ms(start_grip_ix);
            idx_fracs{tr} = (time_ms(start_grip_ix:end_grip_ix) - time_ms(start_grip_ix))./grip_fracs(tr);
            idx_fracs{tr} = round(idx_fracs{tr}, 1, 'decimals')';
        end
    end
end

b_ct = 0;

start_grip_cell = num2cell(data_subj{fn}.start_grip_ix);

for beg = frac_begin
    
    b_ct = b_ct+1;
    
    ixs = cellfun( @(ix_mat,start_grip_ix)( find(ix_mat == beg)+ start_grip_ix-1),...
        idx_fracs, start_grip_cell, 'UniformOutput', false);
    ixs_req(b_ct, :)= ixs;
%     keyboard

end

% Now iterate through all the trials: for each isolate the number of trials
% for which "rest" began in any given fraction. (Note: There may be many
% such beginnings for a single trial".
time_ms = data_subj{fn}.time_ms;
for tr=  1:num_trials

    b_ct_yet = 0;
    ix_drop_yet=[];
    duration_cell{tr} = [];
    if isempty(find(tr == discard_trials))
        for b = frac_begin
            b_ct = find(b==frac_begin);
            rest_flag = 1;
            
            if ~isempty(ixs_req(b_ct,tr)) && length(ixs_req{b_ct,tr})>1
                
                start_ix_beg=  ixs_req{b_ct,tr}(1);
                end_ix_beg = ixs_req{b_ct,tr}(end);
                
                while rest_flag==1
                    
                    if grip_cells{tr}(start_ix_beg) > grip_reqd(tr)
                        
                        ix_drop1 = find(grip_cells{tr}(start_ix_beg+1:end_ix_beg)...
                            <grip_reqd(tr), 1, 'first') + start_ix_beg-1;
                    else
                        ix_back = find(grip_cells{tr}(start_ix_beg+1:end_ix_beg-1)...
                            >= grip_reqd(tr), 1, 'first') + start_ix_beg-1;
                        % if there is a rest period from the prev fraction:
                        if ~isempty(ix_drop_yet) && ~isempty(ix_back)
                            residual_drop =time_ms{tr}(ix_back) -  ...
                                time_ms{tr}(ix_drop_yet);
                            residual_frac = residual_drop/grip_fracs(tr);
                            time_drop_beg = time_ms{tr}(ix_drop_yet) - ...
                                time_ms{tr}(start_grip_cell{tr}(1));
                            duration_cell{tr} = [duration_cell{tr}, residual_drop];
                            time_beg_frac = time_drop_beg/grip_fracs(tr);
                            scatter_mat = [scatter_mat; [time_beg_frac, residual_frac]];
                            residual_frac = round(residual_frac, 1, 'decimals');
                            row_yet_ix = find(residual_frac == frac_duration);
                            duration_matrix(row_yet_ix,b_ct_yet)...
                                = duration_matrix(row_yet_ix, b_ct_yet) + 1;
                        end
                        ix_drop1 = find( grip_cells{tr}(ix_back: end_ix_beg) < ...
                            grip_reqd(tr), 1, 'first')+ ix_back-1;
                    end
                    
                    ct_req = numel(ix_drop1);
                    
                    if ~isempty(ix_drop1)
                        ix_end = find(grip_cells{tr}(ix_drop1+1:end_ix_beg)>=grip_reqd(tr),...
                            1, 'first')+ ix_drop1 -1;
                        
                        if ~isempty(ix_end)
                            
                            time_drop = time_ms{tr}(ix_end) - time_ms{tr}(ix_drop1);
                            time_frac = time_drop/ grip_fracs(tr);
                            duration_cell{tr} = [duration_cell{tr}, time_drop];
                            time_drop_beg = time_ms{tr}(ix_drop1) - ...
                                time_ms{tr}(start_grip_cell{tr}(1));
                            time_beg_frac = time_drop_beg/grip_fracs(tr);
                            
                            scatter_mat=[scatter_mat; [time_beg_frac,time_frac]];
                            
                            time_frac = round(time_frac, 2 ,'decimals');
                            row_ix = find(time_frac==frac_duration, 1 ,'first');
                            duration_matrix(row_ix,b_ct) =...
                                duration_matrix(row_ix,b_ct)+ ct_req;
                        else
                            b_ct_yet= b_ct;
                            ix_drop_yet = ix_drop1;
                        end
                        
                        
                        more_ix = find(grip_cells{tr}(ix_end+1:end_ix_beg) < grip_reqd,...
                            1, 'first')+ ix_end;
                        if isempty(more_ix)
                            rest_flag=0;
                        else
                            start_ix_beg = ix_end+1;
                        end
                        
                    else
                        rest_flag= 0;
                    end
                end
                
            end
        end
        
    end
    if ~isempty(duration_cell{tr})
        data_subj{fn}.avg_rest_dur(tr) = mean(duration_cell{tr});
    else
        data_subj{fn}.avg_rest_dur(tr) = 0;
    end
end

duration_matrix_ct = duration_matrix;
duration_matrix = duration_matrix/ num_trials;
data_subj{fn}.duration_matrix = duration_matrix;

figure(fig_ix);
clf
imagesc(frac_begin, frac_duration, duration_matrix);
% mesh(frac_begin, frac_duration, duration_matrix);
colorbar;
ylabel('Duration of rest period as a fraction of total');
xlabel('Fraction of grip duration at which rest began');
title(' Probability of rest occuring ');

fig_ix = fig_ix+1;
% keyboard

if ~isempty(scatter_mat)
    figure(fig_ix);
    
    [~, sort_ix] = sort(scatter_mat(:,1));
    
    scatter_mat = scatter_mat(sort_ix,:);
    
    scatter(scatter_mat(:,1), scatter_mat(:,2));
    
    fig_ix = fig_ix+1;
    
    xlabel('Fraction of grip duration at which dip/rest begins');
    ylabel('Duration of rest (ms)');
    
end

data_subj{fn}.scatter_mat = scatter_mat;
