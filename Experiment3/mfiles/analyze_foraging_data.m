function analyze_foraging_data(subj, exp, test)
% Code that is created to analyze the data for all subjects that have done
% the two foraging experiments designed by Shruthi Sukumar
% Code accepts 3-letter subject code.

% Make sure all residual graphs are closed
% close all;

if nargin<1
    subj = 'YRT';
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

%Data directory
f = pwd;
cd('..');
dir_file = pwd;
cd(f);

if exist(fullfile('..','plots',subj)) ==0
    mkdir(fullfile('..','plots',subj));
end

% When fn == 1 the contents of the matfile are associated with the low
% effort environment in both experiments. This is because of how filenames
% are declared above.

global data_subj;
clrs = {'b','r'};
fig_ix =1;
for fn = 1:length(filenames)
    
    % Put in checks to ensure that the environment order is right and that
    % the trials are probe/env according to the predetermined order
    % Read data from the existing robot data folder
    
    filename = filenames{fn};
    filedir = fullfile(dir_file,'robot_data', filename);
    ix_1 = length(exp_name)+1;
    subj_mat_filename = strcat(filename(ix_1:end),'.mat');
    discard_mat_name = strcat(filename(ix_1:end),'_discard', '.mat');
        subj_dir = fullfile(dir_file, 'matfiles_data',subj);
    if exp ==2 
        subj_dir = fullfile(subj_dir, 'harvest');
    end
    if ~isdir(subj_dir)
        mkdir(subj_dir);
    end
    mat_filename = fullfile(subj_dir,...
            subj_mat_filename);
    discard_mat_name = fullfile(subj_dir, ...
        discard_mat_name);
    
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
    num_trials = max(T.trialnumber);
    
    %  Deterministic reward parameters: 
    alpha = T.alpha(1); % index is because of the way tcl logged the parameters
    beta = T.beta(1);
    
    % Initialize all required variables
    
    reqd_vars = {'harvest_duration','patch_residence_duration',...
        'grip_drop_duration', 'travel_duration','movement_duration',...
        'pk_vel', 'pk_vel_ix', 'pk_vel_ix_rep','pk_acc_rep','time_ms',...
        'pk_acc_ix_rep','berries', 'grip_norm_avg','grip_valid_ix','grip_force'...
        'grip_begin_ix','vel_profile','vel_profile_rep','x_profile','x_profile_rep',...
        'time_ms_rep', 'probe_trials', 'profile_shift_rep','start_grip_ix', ...
        'end_grip_ix','force_integral', 'force_per_berry','force_rate_max',...
        'force_profile','grip_norm','force_rate','force_rate_max_ix',...
        'max_distance','harvest_reaction_time','grip_ramp_up_dur',...
        'grip_ramp_up2','grip_drop_time', 'berry_ixs','grip_release_time',...
        'avg_rest_dur','avg_grip_bias','avg_grip_variance',...
        'grip_overshoot','hazard_time','time_to_berry',...
        'grip_drop_berry','grip_drop_last_berry','time_to_last_berry',...
        'giving_up_time','norm_dist','norm_pv', 'norm_td', 'norm_hd'};
    
    data_subj{fn} = cell2struct(repmat({nan(1,num_trials)},...
        length(reqd_vars),1),reqd_vars);
    data_subj{fn}.grip_valid_ix = cell(1,num_trials);
    data_subj{fn}.grip_force = cell(1,num_trials);
    data_subj{fn}.force_rate =cell(1,num_trials);
    data_subj{fn}.grip_norm = cell(1,num_trials);
    data_subj{fn}.berry_ixs = cell(1,num_trials);
    data_subj{fn}.vel_profile = cell(1,num_trials);
    data_subj{fn}.vel_profile_rep = cell(1,num_trials);
    data_subj{fn}.x_profile = cell(1,num_trials);
    data_subj{fn}.x_profile_rep = cell(1,num_trials);
    data_subj{fn}.y_profile = cell(1,num_trials);
    data_subj{fn}.y_profile_rep = cell(1,num_trials);
    data_subj{fn}.force_profile = cell(1,num_trials);
    data_subj{fn}.time_ms = cell(1,num_trials);
    data_subj{fn}.time_ms_rep = cell(1,num_trials);
    data_subj{fn}.time_to_berry = cell(1,num_trials);
    data_subj{fn}.grip_drop_berry = cell(1,num_trials);
    
    probe_indices = [21:30, 71:80, 121:130, 171:180];
    
    data_subj{fn}.probe_trials = probe_indices;
    
    grip_tol = 2; % Tolerance of 2 newtons
    time_ms_prev = [];
    outtarget_time = 0;
    outtarget_ix_prev = 1; % for the first trial, it has to be the first index
    discard_prev_flag = 0;
    max_berries =16;
    
    for tr = 1:num_trials
        
        if ~isempty(find(tr==ntrials, 1))
            
            data = T.framedata(tr);
            grip_force = data.grip_force;
            data_subj{fn}.grip_force_cell{tr} = grip_force;
            x= data.x;
            y= data.y;
            vx = data.vx;
            vy = data.vy;
            
            data1 = find_discont(x,y,vx,data,0,1);
            
            vvx = diff23f5(data1.x,1/200,10);
            vvy = diff23f5(data1.y, 1/200, 10);
            
            data_subj{fn}.vel_profile{tr} = vvx(:,2);
            data_subj{fn}.x_profile{tr} = x;
            data_subj{fn}.y_profile{tr} = y;
            score = data.score;
            time_ms = data.time_ms;
            data_subj{fn}.time_ms{tr} = time_ms;
            
            % Getting the Indices for time stamps of when the entered and
            % left patches
            attarget_ix = find(data.statenumber==6, 1, 'first');
            outtarget_ix = find(data.statenumber==6, 1,'last');
            
%             if exp == 1 || ~isempty(find(tr==probe_indices)) 
            grip_req = 30;
%             elseif exp == 2 && fn ==1
%                 grip_req = 10;
%             else 
%                 grip_req = 50;
%             end
            
            data_subj{fn}.grip_norm{tr} = grip_force/grip_req;
            
            if isempty(find(tr == discard_trials))
                
                if discard_prev_flag
                    x_pos = -1*sign(x(1))*x;
                    attarget_ix = find(x_pos > 0.165 , 1 ,'first');
                    discard_prev_flag = 0;
                end
                % #########Find harvest duration#######
                
                compute_harvest_duration(fn ,tr, score, grip_force,...
                    grip_req, grip_tol, time_ms, attarget_ix,...
                    outtarget_ix, alpha, beta, max_berries, 0);
                
                % determine if trials need to be discarded 
                if ~isempty(data_subj{fn}.grip_valid_ix{tr})
                    return_pk_vel(fn, tr,0, vvx, x, y, time_ms,...
                        time_ms_prev);
                    get_mvt_duration(fn, tr, time_ms, attarget_ix,...
                        outtarget_time, outtarget_ix_prev);
                else
                    discard_trials = [discard_trials, tr];
                    save(discard_mat_name, 'discard_trials');
                    continue;
                end
            else 
%                 if fn == 1 && tr == 199
%                     keyboard;
%                 end
                discard_prev_flag = 1;
%                 % determining what happens on those discarded trials
                return_pk_vel(fn, tr,1, vvx,x,y, time_ms, time_ms_prev);
                get_mvt_duration(fn, tr, time_ms, attarget_ix,outtarget_time);
            end
        end
        
        % current val to be used in next trial
        outtarget_time = time_ms(outtarget_ix);
        outtarget_ix_prev = outtarget_ix;
        time_ms_prev = time_ms;
        
    end
    normalized_probe_vals(fn);
    fig_ix = avg_vel_profiles(fn, probe_indices, clrs{fn}, fig_ix);
%     fig_ix = avg_force_rate_profiles(fn, probe_indices, clrs{fn}, fig_ix);
%     saveas(gcf, fullfile('..','plots',subj,strcat('FRM_profile_env',...
%         num2str(fn),'_exp_',num2str(exp))),'png');
%     saveas(gcf, fullfile('..','plots',subj,strcat('FRM_profile_env',...
%         num2str(fn),'_exp_',num2str(exp))),'fig');
%     fig_ix = avg_grip_profiles(fn, probe_indices, clrs{fn}, fig_ix);
%     saveas(gcf, fullfile('..','plots',subj,strcat('grip_profile_env',...
%         num2str(fn),'_exp_',num2str(exp))),'png');
    fig_ix = grip_drop_dyn(fn, discard_trials, exp, fig_ix);
    saveas(gcf, fullfile('..','plots',subj,strcat('grip_drop',...
        num2str(fn),'_exp_',num2str(exp))),'png');
end

fig_ix = compute_hazard_rate(fig_ix, exp);
saveas(gcf, fullfile('..','plots',subj,strcat('hazard_rate',...
        '_exp_',num2str(exp))),'png');
fig_ix = plot_comparisons_env(fig_ix);

fprintf('\nSubj : %s, alpha = %d, beta = %0.2f', ...
            subj, alpha, beta);
% keyboard
% plot_mvt_traces(fig_ix);
fprintf('Saving %s summary data....\n',subj);
save(fullfile(dir_file,'matfiles_data', subj,'data_subj.mat'),'data_subj');

% Function to plot single grip profile and corresponding berry points
function plot_force_(fn, tr, score) 
global data_subj
% to be called inside trial loop afer the foce profiel funciton is defined
f_req = data_subj{fn}.force_profile{tr};
xticks_ = .005*(1:length(f_req));

diff_scores = diff(score);
berry_ixs = data_subj{fn}.berry_ixs{tr};
figure
plot(xticks_, f_req,'b-', 'Linewidth',2);
hold on
plot(xticks_(berry_ixs), f_req(berry_ixs),'ro');

xlabel('Time');
ylabel('Grip Force(N)');
yline(30,'k--');

keyboard

function compare_force_profiles (tr)
global data_subj

f_req1 = data_subj{1}.force_profile{tr};
% cross_ix1 = find(f_req1>=30, 1, 'first');
begin_ix1 = data_subj{1}.grip_begin_ix(tr);
begin_ix1 = find(f_req1(begin_ix1:end) >=10, 1, 'first') +begin_ix1 -1;
f_req2 = data_subj{2}.force_profile{tr};
begin_ix2 = data_subj{2}.grip_begin_ix(tr);
begin_ix2 = find(f_req2(begin_ix2:end) >=10, 1, 'first') + begin_ix2-1;
cross_ix2 = find(f_req2>=30, 1, 'first');
f_req2 = f_req2(begin_ix2:begin_ix2+60);
f_req1 = f_req1(begin_ix1:begin_ix1+60); % ensuring they're the same length

xticks = 0.005*(1:length(f_req2));
plot(xticks, f_req1 , 'b-','Linewidth',1.5);
hold on
% yline(30- f_req1(1), 'b--');
hold on

plot(xticks, f_req2, 'r-','Linewidth',1.5); 
yline(30, 'k--');
xlabel('Time since grip beginning(s)');
ylabel('Force increase from grip beginning(N)');
beautifyfig

f_dif1 = f_req1 - f_req2;
figure,
plot(xticks, f_dif1, 'Linewidth',1.5);
xlabel('Time since grip beginning(s)');
ylabel('Plain Force Difference (N)');
beautifyfig

figure,
fdif2 = (f_req1 - f_req1(1)) - (f_req2 - f_req2(1));
plot(xticks, fdif2, 'Linewidth',1.5);
xlabel('Time since grip beginning(s)');
ylabel('Delta Force Difference (N)');
beautifyfig
keyboard
% Function to compute harvest duration
function compute_harvest_duration(fn ,tr, score, grip_force, grip_req,...
    grip_tol, time_ms, attarget_ix, outtarget_ix, alpha, beta,...
    max_berries, diag01)
% Computes harvest duration

global data_subj;

% Index at which score is reset to zero 
% Sometimes the score isn't reset from the previous trial until a few
% indices after the beginning of a new trial 

score_first0_ix = find(score==0, 1, 'first');

% First index when the subject is in a new patch and starts harvesting
start_grip_ix = score_first0_ix - 1 + ...
    find((grip_force(score_first0_ix:end) >=grip_req-grip_tol) & ...
    score(score_first0_ix:end) >0,1,'first');

%Last point of harvest within patch in this trial
end_grip_ix =  score_first0_ix - 1 + ...
    find(grip_force(score_first0_ix:end) >=grip_req-grip_tol & ...
    score(score_first0_ix:end) >0,1,'last');

grip_drop_ix = start_grip_ix - 1+ find(grip_force(start_grip_ix: ...
    end_grip_ix) < grip_req-grip_tol);

if ~isempty(grip_drop_ix)
    data_subj{fn}.grip_drop_time(tr) = length(grip_drop_ix) * 5; 
    % 5ms is sampling period of robot
else
    data_subj{fn}.grip_drop_time(tr) = 0;
end

if isempty(start_grip_ix)
    data_subj{fn}.harvest_duration(tr) = nan;
    data_subj{fn}.patch_residence_duration(tr) = nan;
    data_subj{fn}.grip_release_time(tr) = nan;
    data_subj{fn}.grip_valid_ix{tr} = [];
    data_subj{fn}.start_grip_ix(tr) = nan;
    data_subj{fn}.end_grip_ix(tr) = nan;
    data_subj{fn}.grip_drop_time(tr) = 0;
else
    data_subj{fn}.harvest_duration(tr) = time_ms(end_grip_ix) - ...
        time_ms(start_grip_ix) - data_subj{fn}.grip_drop_time(tr);
    data_subj{fn}.patch_residence_duration(tr) = time_ms(outtarget_ix)-...
        time_ms(attarget_ix);
%     data_subj{fn}.grip_release_time(tr) = time_ms(outtarget_ix) - ...
%         time_ms(end_grip_ix);
    data_subj{fn}.grip_valid_ix{tr} = setxor([start_grip_ix:end_grip_ix...
        ], grip_drop_ix);
    data_subj{fn}.start_grip_ix(tr) = start_grip_ix;
    data_subj{fn}.end_grip_ix(tr) = end_grip_ix;
end

% Also determine the rate of force generation in a patch
max_force_ix = start_grip_ix - 1 + find(grip_force( start_grip_ix:...
    end_grip_ix) == max(grip_force(start_grip_ix:end_grip_ix)),1,'first');
grip_flipped = grip_force(end:-1:1);
low_thresh = 10; % enforced in tcl code. lower force threshold
start_grip_flipped_ix = length(grip_force) - start_grip_ix + 1;

% find the point when subject begins to grip to get reward
% This compute the rate and max rate immediately after you move into the
% patch until you start harvesting
force_rate = diff23f5(grip_force,1/200,10);
force_rate= force_rate(:,2);

force_rate_flipped = force_rate(end:-1:1);
[force_rate_max, fr_max_ix] = max(force_rate(attarget_ix: start_grip_ix));
fr_max_ix = attarget_ix - 1 + fr_max_ix;
start_ix_rate_flipped = length(force_rate) - start_grip_ix + 1;
fr_max_ix_flipped = length(force_rate) - fr_max_ix + 1;
begin_grip_flipped_ix = find(force_rate_flipped(fr_max_ix_flipped:end) ....
    <= 0.05*force_rate_max, 1, 'first') + fr_max_ix_flipped-1;

eta = 0.08;
while isempty(begin_grip_flipped_ix)
    
    begin_grip_flipped_ix = find(force_rate_flipped(fr_max_ix_flipped:end) ....
        <= eta*force_rate_max, 1, 'first') + fr_max_ix_flipped -1;
    
    eta= eta*1.25;
    fprintf('Reducing eta : %f\n', eta);
end
begin_grip_ix = length(force_rate_flipped) - begin_grip_flipped_ix + 1;

if isempty(begin_grip_ix) || begin_grip_ix<0
    % this is an error
    keyboard
end

data_subj{fn}.force_rate{tr} = force_rate(begin_grip_ix: start_grip_ix)';
% crossing low force threshold at the beginning of gripping
cross_low_thresh_ix = find(grip_force(begin_grip_ix:start_grip_ix)>=low_thresh,...
    1,'first') + begin_grip_ix -1;
% crossing low force threshold at the end of harvest (when grip is
% released)

cross_low_thresh_ix2 = find(grip_force(end_grip_ix:outtarget_ix)>=low_thresh,...
    1,'last') + end_grip_ix -1;

% Computing average bias and variability in force generate for each trial
diff_scores = diff(score);
if begin_grip_ix> attarget_ix
    data_subj{fn}.force_profile{tr} = grip_force(attarget_ix:end)';
    data_subj{fn}.berry_ixs{tr} = find(diff_scores>0)+1 - attarget_ix +1; 
    data_subj{fn}.grip_begin_ix(tr) = begin_grip_ix - attarget_ix + 1;
else
    data_subj{fn}.force_profile{tr} = grip_force(begin_grip_ix:end)';
    data_subj{fn}.berry_ixs{tr} = find(diff_scores>0)+1 - begin_grip_ix +1; 
    data_subj{fn}.grip_begin_ix(tr) = 1 ;
end
data_subj{fn}.grip_force{tr} = grip_force';

avg_grip_force = nanmean(grip_force(start_grip_ix:end_grip_ix));
var_grip_force = nanvar(grip_force(start_grip_ix:end_grip_ix));
data_subj{fn}.avg_grip_bias(tr) = nanmean( grip_force(start_grip_ix:...
    end_grip_ix)- (grip_req-grip_tol) );
data_subj{fn}.avg_grip_variance(tr)= var_grip_force;
% keyboard

grip_reach_thresh_ix = find(grip_force(cross_low_thresh_ix:end) >= ...
    grip_req - grip_tol,...
    1, 'first') + cross_low_thresh_ix - 1;


if ~isempty(grip_reach_thresh_ix)
    [data_subj{fn}.force_rate_max(tr), data_subj{fn}.force_rate_max_ix(tr)] =...
        max(force_rate(cross_low_thresh_ix:grip_reach_thresh_ix));
%     data_subj{fn}.grip_ramp_up_dur(tr) = time_ms(start_grip_ix) - ...
%         time_ms(cross_low_thresh_ix);
    data_subj{fn}.grip_ramp_up_dur(tr) = time_ms(grip_reach_thresh_ix) - ...
        time_ms(cross_low_thresh_ix);
    data_subj{fn}.grip_release_time(tr) = time_ms(cross_low_thresh_ix2) - ...
        time_ms(end_grip_ix);
    data_subj{fn}.grip_ramp_up2(tr) = time_ms(grip_reach_thresh_ix) - ...
        time_ms(begin_grip_ix);
    data_subj{fn}.harvest_reaction_time(tr) = time_ms(cross_low_thresh_ix) - ...
        time_ms(attarget_ix);
else
    data_subj{fn}.force_rate_max(tr) = 0;
end
% Compute the sum of forces generated through the time in patch
data_subj{fn}.force_integral(tr) = sum(grip_force(attarget_ix: ... 
    outtarget_ix))/length(grip_force(attarget_ix:outtarget_ix));

% Compute the integral of force per number of berries
if ~isempty(start_grip_ix)
    data_subj{fn}.berries(tr) = score(end_grip_ix);
    %sum(score(start_grip_ix : end_grip_ix));
    data_subj{fn}.force_per_berry(tr) = sum(grip_force(start_grip_ix:end_grip_ix))...
        ./data_subj{fn}.berries(tr);
else
    data_subj{fn}.berries = 0; 
    data_subj{fn}.force_per_berry(tr) = nan;
end

if diag01
    plot_diagnostics_grip_metrics(grip_force, force_rate, attarget_ix,...
        cross_low_thresh_ix, begin_grip_ix, grip_reach_thresh_ix)
end


% Compute the max overshoot of the subjects when grip begins

% start_grip_ix : index at which they start earning reward
grip_diff = diff(grip_force);

flip_sign_diff_ix = find(force_rate(start_grip_ix:end)<0,1,'first') + ...
    start_grip_ix - 1;

grip_overshoot = max(grip_force(start_grip_ix: flip_sign_diff_ix)) - grip_req;

if ~isempty(grip_overshoot)
    data_subj{fn}.grip_overshoot(tr) = grip_overshoot;
end

% Compute hazard time at each trial
drop_last = 0;
thresh = grip_req - grip_tol;

max_time_hazard = compute_time_berries(alpha, beta, 0.8);

if ~isempty (start_grip_ix) && ~isnan(start_grip_ix)
    
    t0_har = time_ms(start_grip_ix);
    ix_drop = start_grip_ix - 1 + find(grip_force(start_grip_ix:end) < thresh...
        ,1,'first');
    
    if ~isempty(ix_drop)
        ix_drop_curr = ix_drop;
        while (drop_last == 0)
            
            ix_inc_next = ix_drop_curr  + find(grip_force(ix_drop_curr:end) > thresh...
                ,1,'first');
            if ~isempty(ix_inc_next)
                ix_drop_curr=  ix_inc_next + find(grip_force(ix_inc_next:end) < thresh...
                    ,1,'first');
            else
                drop_last = 1;
            end
        end
        time_drop = time_ms(ix_drop);
        time_hazard = time_drop- t0_har;
        data_subj{fn}.hazard_time(tr) = time_hazard;
    else
        data_subj{fn}.hazard_time(tr) = nan;
    end
else
    data_subj{fn}.hazard_time(tr)=nan;
end

% time to each berry
time_to_berry = zeros(1, max_berries) + nan;
grip_drop_berry = zeros(1, max_berries) + nan; 

for b = 1:data_subj{fn}.berries(tr)
    time_stamp_berry_ix = start_grip_ix - 1 + find(score(start_grip_ix:end)...
        ==b,1,'first');
    time_to_berry(b) = time_ms(time_stamp_berry_ix) - time_ms(start_grip_ix);
    grip_drop_berry_ix = start_grip_ix - 1+ find(grip_force(start_grip_ix: ...
    time_stamp_berry_ix) < grip_req-grip_tol);
    grip_drop_berry(b) = length(grip_drop_berry_ix) * 5; 
end

data_subj{fn}.time_to_berry{tr} = time_to_berry;
data_subj{fn}.grip_drop_berry{tr} = grip_drop_berry;
data_subj{fn}.time_to_last_berry(tr) = time_to_berry(data_subj{fn}.berries(tr));
data_subj{fn}.grip_drop_last_berry(tr) = grip_drop_berry(data_subj{fn}.berries(tr));
data_subj{fn}.giving_up_time(tr) = time_ms(end_grip_ix) - ...
    time_ms(time_stamp_berry_ix);


function plot_diagnostics_grip_metrics(grip_force, force_rate, attarget_ix,...
    cross_low_thresh_ix, begin_grip_ix, grip_reach_thresh_ix)
% to be called from inside compute_harvest_duration

figure(1)
clf
h = gcf;
plot(grip_force)
hold on

plot(attarget_ix, grip_force(attarget_ix),'r*')
hold on
% plot(cross_low_thresh_ix, grip_force(cross_low_thresh_ix),'b+')
yline(grip_force(cross_low_thresh_ix),'b:')
hold on
plot(begin_grip_ix, grip_force(begin_grip_ix),'k^')
hold on
xline(begin_grip_ix, ':');
hold on
xline(grip_reach_thresh_ix, ':');
hold on
% plot(start_grip_ix, grip_force(start_grip_ix),'ro')
yline(grip_force(grip_reach_thresh_ix),'k--')
hold on

yyaxis right
plot(force_rate);

hold off
keyboard

% h=figure;
% plot(1:max_berries, time_to_berry,'-o');
% xlabel('Berry #')
% ylabel('Time to berry #');
% keyboard;
% clf;
% close(h)


% Avraging grip profiles

function fig_ix = avg_grip_profiles(fn, probe_indices, c, fig_ix)
% function to generate average grip force profiles

global data_subj

% begin_ix_mat = cell2mat(data_subj{fn}.grip_begin_ix);\

figure(fig_ix)
% for p =1:2
% 
%     if p ==1
%         begin_ix_mat = [data_subj{fn}.grip_begin_ix(non_probe_indices)]';
%         align_average(data_subj{fn}.grip_force(non_probe_indices), ...
%             begin_ix_mat,c);
%         hold on;
%     else
%         begin_ix_mat = [data_subj{fn}.grip_begin_ix(probe_indices)]';
%         align_average(data_subj{fn}.grip_force(probe_indices), ...
%             begin_ix_mat,'g');
%         hold off;
%     end
%     
% end

begin_ix_mat = [data_subj{fn}.grip_begin_ix]';
align_average(data_subj{fn}.force_profile, begin_ix_mat,c,...
    1, probe_indices);
hold on

fig_ix = fig_ix + 1;

function fig_ix = avg_force_rate_profiles(fn, probe_indices, c, fig_ix)
% function to generate average grip force profiles

global data_subj

% begin_ix_mat = cell2mat(data_subj{fn}.grip_begin_ix);\

figure(fig_ix)
% for p =1:2
% 
%     if p ==1
%         begin_ix_mat = [data_subj{fn}.grip_begin_ix(non_probe_indices)]';
%         align_average(data_subj{fn}.grip_force(non_probe_indices), ...
%             begin_ix_mat,c);
%         hold on;
%     else
%         begin_ix_mat = [data_subj{fn}.grip_begin_ix(probe_indices)]';
%         align_average(data_subj{fn}.grip_force(probe_indices), ...
%             begin_ix_mat,'g');
%         hold off;
%     end
%     
% end

begin_ix_mat = [data_subj{fn}.force_rate_max_ix]';
align_average(data_subj{fn}.force_rate, begin_ix_mat,c,...
    1, probe_indices);
hold on

fig_ix = fig_ix + 1;

function fig_ix = avg_vel_profiles(fn, probe_indices, c, fig_ix)
% function to generate average grip force profiles

global data_subj

figure(fig_ix)

max_ix_mat = [data_subj{fn}.pk_vel_ix_rep]';
align_average(data_subj{fn}.vel_profile_rep, max_ix_mat,c,...
    1, probe_indices);
hold on

fig_ix = fig_ix + 1;

function time_since = compute_time_berries ( alpha, beta, frac_no) 
% alpha here represents the maximum number of berries
% beta represents the rate of decay
% frac_no represents the fraction of the total number of berries which will
% be the end time for hazard computation

time_comp = @(frac,alpha,beta)( (1/beta)*(1/ frac - 1));

time_since = time_comp(frac_no,alpha, beta);

function return_pk_vel(fn,tr,is_discard,vvx, x, y, time_ms, time_ms_prev)
% Returns index of peak velocity and the value
% Also repairs the velocity profiles : trial i vel profile is the movement
% to patch # i

global data_subj;

if ~isempty(data_subj{fn}.grip_valid_ix{tr}) 
    start_grip_ix = data_subj{fn}.grip_valid_ix{tr}(1);
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

if tr ~=1 && ~isempty(attarget_ix)
    data_subj{fn}.travel_duration(tr) = time_ms(attarget_ix) - outtarget_time;
elseif isempty(attarget_ix) 
    data_subj{fn}.travel_duration(tr) = time_ms(end) - outtarget_time;
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


grip_valid_ix_rep = data_subj{fn}.grip_valid_ix{tr} + ...
    data_subj{fn}.profile_shift_rep(tr) -1;

if ~isempty(grip_valid_ix_rep)
    if grip_valid_ix_rep(1) < offset_ix_rep
        offset_ix_rep = grip_valid_ix_rep(1);
    end
end

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
% % h=figure(1)
% % % WinOnTop(h,true)
% % % plot(data_subj{fn}.x_profile_rep{tr})
% % % hold on
% % plot(vel_profile_rep);
% % hold on
% % % plot(sd_vec,'k')
% % % hold on
% % plot(onset_ix_rep, vel_profile_rep(onset_ix_rep),'r^');
% % hold on
% % plot(offset_ix_rep ,vel_profile_rep(offset_ix_rep),'g^');
% % hold on
% % plot(pk_vel_ix_rep ,vel_profile_rep(pk_vel_ix_rep),'k*');
% % 
% % 
% % 
% % % yyaxis right
% % % 
% % % plot(acc_profile_rep)
% % % hold on 
% % % plot(data_subj{fn}.pk_acc_ix_rep(tr), acc_profile_rep(...
% % %     data_subj{fn}.pk_acc_ix_rep(tr)),'o');
% % %  
% % hold off 
% % keyboard;
% % clf


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
        
        nan_idx = (isnan(x) | isnan(y));
        scatter(x(~nan_idx),y(~nan_idx), clrs{env});
        [coeffs,S] = polyfit(x(~nan_idx), y(~nan_idx), 1);
        fittedX = linspace(min(x(~nan_idx)), max(x(~nan_idx)), 199);
        fittedY = polyval(coeffs, fittedX);
        r_squared = 1 - (S.normr/norm(y(~nan_idx) - mean(y(~nan_idx))))^2;
        fprintf('R squared for correlation of %s to %s - %d\n\n', yy, xx,...
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
metrics = {'pk_vel', 'harvest_duration','movement_duration',...
    'travel_duration','force_rate_max','grip_ramp_up_dur', ...
    'grip_ramp_up2','harvest_reaction_time', 'avg_grip_bias',...
    'grip_overshoot','avg_grip_variance'};
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
        elseif strcmp(metric, 'harvest_reaction_time')
            plot_dat = data_subj{env}.harvest_reaction_time;
            tt = 'Harvest Rxn Time';
        elseif strcmp(metric, 'force_rate_max')
            plot_dat = data_subj{env}.force_rate_max;
            tt = 'Max Force Rate';
        elseif strcmp(metric, 'grip_ramp_up_dur')
            plot_dat = data_subj{env}.grip_ramp_up_dur;
            tt= 'Grip Ramp Up Duration';
        elseif strcmp(metric, 'grip_ramp_up2')
            plot_dat = data_subj{env}.grip_ramp_up2;
            tt= 'Grip Ramp Up Duration2';
        elseif strcmp(metric,'harvest_duration')
            plot_dat = data_subj{env}.harvest_duration;
            tt = 'Harvest duration';
        elseif strcmp(metric, 'patch_residence')
            plot_dat = data_subj{env}.patch_residence_duration;
            tt = 'Patch residence duration';
        elseif strcmp(metric, 'grip_overshoot')
            plot_dat = data_subj{env}.grip_overshoot;
            tt = 'Overshoot Grip Force';
        elseif strcmp(metric,'avg_grip_bias')
            plot_dat = data_subj{env}.avg_grip_bias;
            tt='Average Grip Bias';
        elseif strcmp(metric,'avg_grip_variance')
            plot_dat = data_subj{env}.avg_grip_variance;
            tt='Average Grip Variance';
        end
        probe_indices = data_subj{env}.probe_trials;
        plot(1:length(plot_dat), plot_dat, clrs{env});
        hold on
        title(tt,'Fontsize',15);
        legend('low','high');
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
    
    patch([probe_indices(1) probe_indices(1) probe_indices(end/4)...
        probe_indices(end/4)] ,[min_val max_val max_val min_val] , 'r');
    alpha(0.2)
    patch([probe_indices(end/4+1) probe_indices(end/4+1) probe_indices(end/2)...
        probe_indices(end/2)] ,[min_val max_val max_val min_val] , 'r');
    alpha(0.2)
    patch([probe_indices(end/2+1) probe_indices(end/2+1) probe_indices(end/4)+100 ...
        probe_indices(end/4)+100] ,[min_val max_val max_val min_val] , 'r');
    alpha(0.2)
    patch([probe_indices(end/4+1)+100 probe_indices(end/4+1)+100 probe_indices(end) ...
        probe_indices(end)] ,[min_val max_val max_val min_val] , 'r');
    alpha(0.2)
    
    hold off;
    beautifyfig;
    m_ix = m_ix + 1;
    if mod(m_ix,4) == 1
        fig_ix = fig_ix +1;
        m_ix = 1;
    end
    
end
% keyboard

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
