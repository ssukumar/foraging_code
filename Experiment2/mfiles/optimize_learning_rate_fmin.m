function optimize_learning_rate_fmin(df)
% New foraging data set optimizing the alpha/learning rates per subject 
% df: dataframe from R 

global df_tbl

% read to table
df_tbl = readtable('df_from_R2.csv','ReadVariableNames',1);
method=1;
% get list of subjects 

subj_list = unique(df_tbl.Subj_id);
n_subjs = length(subj_list);
% keyboard
fprintf('\n');

rate_update= zeros(height(df_tbl),1);
cumul_rate_new= zeros(height(df_tbl),1);
init_alpha = 0.02;
init_beta=0.1;
init_gamma=2;

% init_alph_bet = [init_alpha, init_beta];

if ~ismember('hd_log',df_tbl.Properties.VariableNames)
    hd_log = log10(cellfun(@conv2double,df_tbl.hd_scale));
    df_tbl = addvars(df_tbl, hd_log,'NewVariableNames','hd_log');
end
if method ==1 || method==3
    [alpha_opt, fopt_met] = fminsearchbnd( @(alpha) ensemble_model(alpha,init_beta, init_gamma, subj_list, method),...
        init_alpha,0, 2);
    fprintf('\n\nBest estimate; alpha, %0.3f, AIC, %f\n',...
    alpha_opt, fopt_met);
elseif method ==2 
    [beta_opt, fopt_met] = fminsearchbnd( @(beta) ensemble_model(init_alpha, beta, init_gamma,subj_list,method),...
        init_beta,0, 2);
    fprintf('\n\nBest estimate; beta, %0.3f, AIC, %f\n',...
    beta_opt, fopt_met);
else 
    [gamma_opt, fopt_met] = fminsearchbnd( @(gamma) ensemble_model(init_alpha, init_beta, gamma,...
        subj_list,method), init_gamma ,0, 50);
    fprintf('\n\nBest estimate; gamma %0.3f; AIC, %f\n',...
    gamma_opt, fopt_met);
end
    


keyboard;
figure(1);

for s = 1:n_subjs
    

    subplot(2,5,s)
%     init_alpha = rand;
%     [alpha_opt, fopt_rsq] = fminsearchbnd(@(alpha) optimize_subj(alpha,subj_list(s)),...
%         init_alpha, 0, 10);
%     fprintf('\n\nBest estimate; alpha, %0.3f, R squared, %1.4f',...
%         alpha_opt, -1*fopt_met);
%     [alpha_opt, fopt_aic] = fminsearchbnd(@(alpha) optimize_subj(alpha,subj_list(s)),...
%         init_alpha, 0, 10);
%     fprintf('\n\nBest estimate; alpha, %0.3f, AIC, %f\n',...
%         alpha_opt, fopt_aic);   
    fprintf('\n**********************************************************\n');
    if method ==1
        [~, rate_upd_best, ~,~,~] = optimize_subj(alpha_opt, init_beta,init_gamma, subj_list(s),0,0);

    elseif method==2
        [~,~, rate_upd_best,~,~] = optimize_subj(init_alpha, beta_opt,init_gamma, subj_list(s),0,0);
    elseif method==3
        [~, ~, ~,rate_upd_best,~] = optimize_subj(alpha_opt, init_beta,init_gamma, subj_list(s),0,0);
    else
        [~, ~, ~,~,rate_upd_best] = optimize_subj(init_alpha, init_beta,gamma_opt, subj_list(s),0,0);
    end
     rate_update(strcmp(df_tbl.Subj_id, subj_list(s))) =  rate_upd_best;
    
end
if method==1
    df_tbl = addvars(df_tbl, rate_update,'NewVariableNames','rate_update_new');
elseif method ==2 
    df_tbl = addvars(df_tbl, rate_update,'NewVariableNames','cumul_rate_beta');
elseif method==3
    df_tbl = addvars(df_tbl, rate_update,'NewVariableNames','rate_update_perf');
else
    df_tbl = addvars(df_tbl, rate_update,'NewVariableNames','util_rate');
end
keyboard;

writetable(df_tbl, '../rfiles/df_from_matlab.csv','WriteVariableNames',1);


function [met_avg] = ensemble_model (alpha, beta, gamma, subj_list, method )
% Optimizing across multiple subjects

global df_tbl

if nargin<2
    subset =0;
end

n_subjs = length(subj_list);
% add median split functionality with subset 

met_all = zeros(1,n_subjs); 
% ensembling the models 
% for s = 1:n_subjs 
%     [met_subj, ~] = optimize_subj(alpha, subj_list(s),0,1); 
% %     rate_update = [rate_update; rate_upd_best];
%     met_all(s) = met_subj;
%     
% end

% substituting this loop with LMER model from best performing returned
% model

if method ==1 
    for s = 1: n_subjs
        
        [~, rate_upd_col,~,~,~] = optimize_subj(alpha,beta,gamma, subj_list(s), 1,0);
        df_tbl.rate_update(strcmp(df_tbl.Subj_id, subj_list(s))) =  rate_upd_col;
        
    end
    lmer_all = fitlme(df_tbl, ...
        'hd_log ~ rate_update * reward + trial_bin_ord_c * Order + (1|Subj_id)');
elseif method ==2
    
    for s = 1: n_subjs
        
        [~,~, rate_upd_col,~,~] = optimize_subj(alpha,beta,gamma, subj_list(s), 1,0);
        df_tbl.cumul_rate_new(strcmp(df_tbl.Subj_id, subj_list(s))) = ...
            rate_upd_col;
        
    end
    
    lmer_all = fitlme(df_tbl, ...
        'hd_log ~ cumul_rate_new * reward + trial_bin_ord_c*Order + (1|Subj_id)');
    
elseif method ==3
    for s = 1: n_subjs
        
        [~,~,~,rate_upd_col,~] = optimize_subj(alpha,beta,gamma, subj_list(s), 1,0);
        df_tbl.rate_perf(strcmp(df_tbl.Subj_id, subj_list(s))) =  rate_upd_col;
        
    end
    lmer_all = fitlme(df_tbl, ...
        'hd_log ~ rate_perf * reward + trial_bin_ord_c*Order + (1|Subj_id)');
else
    for s = 1: n_subjs
        
        [~,~,~,~,rate_upd_col] = optimize_subj(alpha,beta,gamma, subj_list(s), 1,0);
        df_tbl.util_rate(strcmp(df_tbl.Subj_id, subj_list(s))) =  rate_upd_col;
        
    end
    lmer_all = fitlme(df_tbl, ...
        'hd_log ~ util_rate * reward + trial_bin_ord_c*Order + (1|Subj_id)');
end


met_avg= -1* lmer_all.LogLikelihood;
lmer_all.Coefficients
% met_avg = sum(met_all);

fprintf('Current estimate; alpha, %0.3f,beta %0.3f ,gamma %0.3f , LogLik %f, AIC, %f\n',...
    alpha, beta, gamma, lmer_all.ModelCriterion.AIC , met_avg);

function [met_subj, rate_update, cumul_rate_new,rate_perf, util_rate] = optimize_subj(alpha,beta, gamma, subj_id,...
    plt_debug, opt_bool)


if nargin<3
    plt_debug=0;
end

global df_tbl

df_tmp = df_tbl(strcmp(df_tbl.Subj_id,subj_id),: );

rate_update = zeros(size(df_tmp,1),1);
cumul_rate_new = zeros(size(df_tmp,1),1);
rate_perf = zeros(size(df_tmp,1),1);
util_rate = zeros(size(df_tmp,1),1);
[trial_bin_ord_c, trial_asc_idx] = sort(df_tmp.trial_bin_ord_c);

% comment the following two lines to check for another metric while fitting
% alpha
pk_vel_out = cellfun(@conv2double,df_tmp.pk_vel_out);
pk_vel_out = pk_vel_out(trial_asc_idx);

hd_scale = cellfun(@conv2double,df_tmp.hd_scale);
hd_scale = hd_scale(trial_asc_idx);
hd_scale(isnan(hd_scale)) = eps;
hd_log = log10(hd_scale);

Order = cellfun(@conv2double,df_tmp.Order);
Order =Order(trial_asc_idx);

reward = df_tmp.reward;
reward = reward(trial_asc_idx);
rate_update(trial_asc_idx(1))=0;
cumul_rate_new(trial_asc_idx(1))=0;
rate_perf(trial_asc_idx(1)) = 0;
util_rate(trial_asc_idx(1)) = 0;

berries =  df_tmp.berries; 
berries = berries(trial_asc_idx); 

travel_scale = cellfun(@conv2double, df_tmp.travel_scale); 
travel_scale = travel_scale(trial_asc_idx); 
travel_scale(isnan(travel_scale))=0;

assert(size(berries,2)==1, 'Wrong dimensions; check math');

t=2;
for tr = trial_bin_ord_c(2:end)'
%     reward = table2array(df_tmp(df_tmp.trial_bin_ord_c==tr,{'reward'}));
    rate_update(trial_asc_idx(t)) = rate_update(trial_asc_idx(t-1)) + ...
        alpha* (reward(t) - rate_update(trial_asc_idx(t-1)) );
    
    reward_vec = berries(1:t-1)'*exp(-1*beta*(t- [1:t-1]'));
    time_vec = (hd_scale(1:t-1) + travel_scale(1:t-1))'*exp(-1*beta*(t- [1:t-1]'));
    
    assert(length(reward_vec) ==1 && length(time_vec)==1,...
        'Something\`s wrong with the math');
    
    cumul_rate_new(trial_asc_idx(t)) = reward_vec/time_vec;
    
    rate_perf(trial_asc_idx(t)) = rate_perf(trial_asc_idx(t-1)) + ...
        alpha* (berries(t)/(hd_scale(t)+travel_scale(t)) - rate_perf(trial_asc_idx(t-1)) );
    
%     travel_scale(travel_scale==0) =NaN;
    
    util_vec = nansum(berries(1:t-1) - gamma./(travel_scale(1:t-1)));
    time_vec_all = nansum(hd_scale(1:t-1) + travel_scale(1:t-1));
    util_rate(trial_asc_idx(t)) = util_vec/time_vec_all;
    
    t=t+1;
end

if plt_debug
    clf
    plot(df_tmp.trial_bin_ord_c, rate_update);
%     keyboard;
    drawnow
end
 
if opt_bool
    tbl = table(rate_update, cumul_rate_new, hd_log, reward, trial_bin_ord_c, Order, ...
        'VariableNames',{'rate_update','cumul_rate_new', 'hd_log', 'reward' , 'trial_bin_ord_c',...
        'Order'});
    % addvars(df_tmp ,rate_update);
    
    mdl_subj = fitlm(tbl,'hd_log ~cumul_rate_new + reward +trial_bin_ord_c' ,...
        'RobustOpts','on');
    
    if plt_debug
        %     plotEffects(mdl_subj)
        plotAdded(mdl_subj, 'rate_update')
        mdl_subj.Coefficients
        drawnow
    end
    
    r_sq= -1*mdl_subj.Rsquared.Ordinary;
    AIC_val = mdl_subj.ModelCriterion.AIC;
    
    met_subj = AIC_val;
else
    met_subj=NaN;
end
% fprintf('Current estimate; alpha, %0.3f, R squared, %1.4f, AIC, %f\n',...
%     alpha,-1*r_sq, AIC_val);
%  
function op = conv2double(x)

if strcmp(x,'NA')
    op = NaN;
else 
    op = str2double(x);
end
