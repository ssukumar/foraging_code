function optimize_learning_rate_fmin2(df)
% New foraging data set optimizing the alpha/learning rates per subject 
% df: dataframe from R 

global df_tbl

% read to table
df_tbl = readtable('df_from_R.csv','ReadVariableNames',1);

% get list of subjects 

subj_list = unique(df_tbl.Subj_id);
n_subjs = length(subj_list);
% keyboard
fprintf('\n');

rate_update= zeros(height(df_tbl),1);
init_alpha = [0.02,0.06];

hd_log = log10(cellfun(@conv2double,df_tbl.hd_scale));
df_tbl = addvars(df_tbl, hd_log,'NewVariableNames','hd_log');

[alpha_opt, fopt_met] = fminsearchbnd( @(alpha) ensemble_model(alpha,subj_list),...
    init_alpha,[0,0], [2,2]);
fprintf('\n\nBest estimate; alpha, %0.3f %0.3f, AIC, %f\n',...
    alpha_opt(1), alpha_opt(2) , fopt_met);

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
    
    [~, rate_upd_best] = optimize_subj(alpha_opt, subj_list(s),0,0);
%     rows_subj = find()
%     keyboard;
    rate_update(strcmp(df_tbl.Subj_id, subj_list(s))) =  rate_upd_best;
    
end

df_tbl = addvars(df_tbl, rate_update,'NewVariableNames','rate_update_new');
keyboard;

writetable(df_tbl, '../rfiles/df_from_matlab2.csv','WriteVariableNames',1);


function [met_avg] = ensemble_model (alpha, subj_list )
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

for s = 1: n_subjs
    
    [~, rate_upd_col] = optimize_subj(alpha, subj_list(s), 1,0);
    df_tbl.rate_update(strcmp(df_tbl.Subj_id, subj_list(s))) =  rate_upd_col;

end
lmer_all = fitlme(df_tbl, ...
    'hd_log ~ rate_update * reward + trial_bin_ord_c*Order + (1|Subj_id)');
met_avg= -1* lmer_all.LogLikelihood;
lmer_all.Coefficients
% met_avg = sum(met_all);

fprintf('Current estimate; alpha, %0.3f %0.3f, LogLik %f, AIC, %f\n',...
    alpha(1), alpha(2) ,lmer_all.ModelCriterion.AIC , met_avg);

function [met_subj, rate_update] = optimize_subj(alpha, subj_id, plt_debug, opt_bool)


if nargin<3
    plt_debug=0;
end

assert(length(alpha)==2, 'Alpha must contain two rate parameters');

global df_tbl

df_tmp = df_tbl(strcmp(df_tbl.Subj_id,subj_id),: );

rate_update = zeros(size(df_tmp,1),1);
[trial_bin_ord_c, trial_asc_idx] = sort(df_tmp.trial_bin_ord_c);

% comment the following two lines to check for another metric while fitting
% alpha
pk_vel_out = cellfun(@conv2double,df_tmp.pk_vel_out);
pk_vel_out = pk_vel_out(trial_asc_idx);

hd_scale = cellfun(@conv2double,df_tmp.hd_scale);
hd_scale = hd_scale(trial_asc_idx);
hd_log = log10(hd_scale);

Order = cellfun(@conv2double,df_tmp.Order);
Order =Order(trial_asc_idx);

reward = df_tmp.reward;
reward = reward(trial_asc_idx);
rate_update(trial_asc_idx(1))=0;

t=2;
for tr = trial_bin_ord_c(2:end)'
    
    k = ((reward(t) - rate_update(trial_asc_idx(t-1)) ) >=0) +1;
    
    rate_update(trial_asc_idx(t)) = rate_update(trial_asc_idx(t-1)) + ...
        alpha(k)* (reward(t) - rate_update(trial_asc_idx(t-1)) );
    t=t+1;
end

if plt_debug
    clf
    plot(df_tmp.trial_bin_ord_c, rate_update);
%     keyboard;
    drawnow
end

if opt_bool
    tbl = table(rate_update, hd_log, reward, trial_bin_ord_c, Order, ...
        'VariableNames',{'rate_update', 'hd_log', 'reward' , 'trial_bin_ord_c',...
        'Order'});
    % addvars(df_tmp ,rate_update);
    
    mdl_subj = fitlm(tbl,'hd_log ~rate_update + reward +trial_bin_ord_c' ,...
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
