library('ggplot2')
require('reshape2')
require(lmerTest)
library('lme4')
library(sjPlot)
library('pwr')
library('dplyr')
library('car')
library(readr)
setwd("C:/Users/shrut/Documents/Research/foraging_reward_new/rfiles")
## read in experiment 3 data

data.exp3 <- read_rds('exp3_dataframe.rds')

data_tmp <- data.exp3 %>% group_by(Subj_id, Env) %>% 
  mutate( trial_bin_ord = if_else( Order!=Env,1+round(Trial/max(Trial,na.rm=TRUE), 1),round(Trial/max(Trial,na.rm=TRUE), 1) ))
data_tmp$trial_bin_ord = as.factor(data_tmp$trial_bin_ord)
data.exp3 <- data_tmp

st.err <- function(x) {
  sd(x)/sqrt(length(x))
}

create_agg <- function(df, form, varlist, fn=mean, fn_spread= sd) {
  
  df.agg = aggregate(formula = form,data=df,FUN= fn)
  cnames = colnames(df.agg)
  df.agg_sd = aggregate(formula = form,data= df, FUN=fn_spread)
  
  varlist = varlist
  
  # varlist
  for (i in 1:length(varlist)) {
    v= varlist[i]
    max_varname = paste(v,'_max', sep="")
    min_varname = paste(v,'_min', sep="")
    # print(min_varname)
    df.agg = cbind.data.frame(df.agg, df.agg[v] - df.agg_sd[v], df.agg[v] + df.agg_sd[v])
    cnames = c(cnames, min_varname, max_varname)
    
  }
  
  colnames(df.agg) <- cnames
  return(df.agg)
}

### subject plots PV

agg_pv_subj = aggregate(data=data.exp3, cbind(pk_vel,pv_norm) ~ Env+ord_labels+Subj_id, mean)
cnames=colnames(agg_pv_subj)
agg_pv_subj_sd =  aggregate(data=data.exp3, cbind(pk_vel,pv_norm) ~ Env+ord_labels+Subj_id, sd)
varlist = c('pk_vel','pv_norm')

varlist
for (i in 1:length(varlist)) {
  v= varlist[i]
  max_varname = paste(v,'_max', sep="")
  min_varname = paste(v,'_min', sep="")
  print(min_varname)
  agg_pv_subj = cbind.data.frame(agg_pv_subj, agg_pv_subj[v] - agg_pv_subj_sd[v], agg_pv_subj[v] + agg_pv_subj_sd[v])
  cnames = c(cnames, min_varname, max_varname)
  
}

colnames(agg_pv_subj) <- cnames

## 

plt_pv_subj <- ggplot(data=agg_pv_subj)+geom_bar(aes(x=Env, y=pv_norm, fill=Env), stat="identity")+
  geom_errorbar(aes(x=Env,ymin=pv_norm_min, ymax=pv_norm_max)) +
  facet_wrap(ord_labels~ Subj_id, scales="free",nrow=2)+
  ylab('Peak Velocity (norm to subj avg.)')+xlab('Environment reward')+
  theme_classic()

plt_pv_subj <- plt_pv_subj + scale_fill_discrete(name="Environemnt",
                                               breaks=c(0,1),
                                               labels=c("Low Reward",'High Reward'))

plt_pv_subj

ggsave('pv_plt_subj.pdf',plt_pv_subj, width=8, height=4)
ggsave('pv_plt_subj.png',plt_pv_subj, width=8, height=4)

## Bar plot for peak velocity as a funciton of environment
agg_pv_all = aggregate(data=agg_pv_subj, cbind(pk_vel,pv_norm) ~ Env, mean)
cnames=colnames(agg_pv_all)
agg_pv_all_sd =  aggregate(data=agg_pv_subj, cbind(pk_vel,pv_norm) ~ Env, st.err)
varlist = c('pk_vel','pv_norm')

varlist
for (i in 1:length(varlist)) {
  v= varlist[i]
  max_varname = paste(v,'_max', sep="")
  min_varname = paste(v,'_min', sep="")
  print(min_varname)
  agg_pv_all = cbind.data.frame(agg_pv_all, agg_pv_all[v] - agg_pv_all_sd[v], 
                                agg_pv_all[v] + agg_pv_all_sd[v])
  cnames = c(cnames, min_varname, max_varname)
  
}

colnames(agg_pv_all) <- cnames

plt_pv_bar <- ggplot()+geom_bar(data=agg_pv_all, aes(x=Env, y=pv_norm, fill=Env), stat="identity")+
  geom_errorbar(data=agg_pv_all, aes(x=Env,ymin=pv_norm_min, ymax=pv_norm_max), width=0) +
  ylab('Peak Velocity (norm to subj avg.)')+xlab('Environment reward')+
  theme_classic()


plt_pv_bar <- plt_pv_bar + geom_line(data=agg_pv_subj, 
                                     aes(x=Env, y= pv_norm, group=Subj_id))
  # geom_errorbar(data=agg_pv_subj, aes(x=Env, ymin= pv_norm_min, ymax = pv_norm_max,
  #                                     color=ord_labels), width=0.0)

plt_pv_bar
ggsave('plt_pv_bar.pdf',plt_pv_bar, width=4, height=3)
ggsave('plt_pv_bar.png',plt_pv_bar, width=4, height=3)

## time dependency of peak velocity 

agg_pv_time = create_agg(df= data.exp3, form=formula(pv_norm ~ trial_bin_ord + Subj_id), 
                         varlis=c('pv_norm'))
agg_pv_time_all = create_agg(df = agg_pv_time, form = formula(pv_norm ~ trial_bin_ord),
                            varlist = c('pv_norm'))

plt_pv_time<- ggplot()+geom_line(data =agg_pv_time_all,aes(x=trial_bin_ord, y=pv_norm, group=1),
                                                      size=2)+
  geom_errorbar(data =agg_pv_time_all, aes(x=trial_bin_ord, ymin =pv_norm_min, ymax=pv_norm_max),width=0.0)+
  theme_classic()
plt_pv_time

ggsave('../plots/plt_pv_time.pdf', plt_pv_time, width=6, height=4)

plt_pv_time <- plt_pv_time + geom_line(data=agg_pv_time ,aes(x=trial_bin_ord, y=pv_norm, 
                                                             group=Subj_id),alpha=0.5)
ggsave('../plots/plt_pv_time_subj.pdf', plt_pv_time, width=6, height=4)

#######Harvest duration


agg_dat_subj = aggregate(data= data.exp3, cbind(hd_scale,hd_norm)~ Env+alpha_vals+Subj_id+ord_labels, mean)
cnames = colnames(agg_dat_subj)

agg_sd_subj = aggregate(data= data.exp3, cbind(hd_scale,hd_norm)~ Env+alpha_vals+Subj_id+ord_labels, sd)

varlist = c('hd_scale','hd_norm')

varlist
for (i in 1:length(varlist)) {
  v= varlist[i]
  max_varname = paste(v,'_max', sep="")
  min_varname = paste(v,'_min', sep="")
  # print(min_varname)
  agg_dat_subj = cbind.data.frame(agg_dat_subj, agg_dat_subj[v] - agg_sd_subj[v], agg_dat_subj[v] + agg_sd_subj[v])
  cnames = c(cnames, min_varname, max_varname)
  
}

colnames(agg_dat_subj) <- cnames

plt_hd_subj <- ggplot(data=agg_dat_subj)+geom_line(aes(x=alpha_vals, y=hd_norm, color=Env, group=Env))+
  geom_errorbar(aes(x=alpha_vals,ymin=hd_norm_min, ymax=hd_norm_max, color=Env), width=0.0) +
  facet_wrap(ord_labels~ Subj_id, scales="free",nrow=2)+
  ylab('Harvest Duration (norm to subj avg.)')+xlab('Initial patch reward')+
  theme_classic()

plt_hd_subj <- plt_hd_subj + scale_color_discrete(name="Environment",
                                                 breaks=c(0,1),
                                                 labels=c("Low Reward",'High Reward'))

plt_hd_subj

ggsave('hd_plt_subj.pdf',plt_hd_subj, width=8, height=4)
ggsave('hd_plt_subj.png',plt_hd_subj, width=8, height=4)

## Get similar trial series plot 
# agg_dat_tr_subj <- create_agg(df=data.exp3, form = formula(cbind(hd_norm, berries, grip_count)~ Trial +Env+alpha_vals +Subj_id), 
#                               varlist=c('hd_norm'), fn=mean, fn_spread= st.err)
# agg_dat_tr_all <- create_agg(df= agg_dat_tr_subj, form = formula(hd_norm~ trial_bin_ord_c +Env+alpha_vals),fn=mean, fn_spread= st.err)

# ggplot(data=agg_dat_tr_subj)+ geom_line(aes(x=Trial, y= hd_norm, color=Env, group= interaction (Env, Subj_id)))+
#   facet_wrap(~alpha_vals, nrow =2, scales="free")+
#   theme_classic()

agg_tmp <- data.exp3 %>% group_by(Subj_id,alpha_vals,Env) %>% mutate(tr_in_alpha=row_number())
data.exp3 = agg_tmp


ggplot(data=agg_dat_tr_subj)+ geom_line(aes(x=tr_in_alpha, y= hd_norm, color=Env, group= interaction (Env, Subj_id)))+
  facet_wrap(~alpha_vals, nrow =2, scales="free")+
  theme_classic()

ggplot(data=agg_dat_tr_subj)+ geom_line(aes(x=tr_in_alpha, y= berries, color=Env, group= interaction (Env, Subj_id)))+
  facet_wrap(~alpha_vals, nrow =2, scales="free")+
  theme_classic()

ggplot(data=agg_dat_tr_subj)+ geom_line(aes(x=tr_in_alpha, y= grip_count, color=Env, group= interaction (Env, Subj_id)))+
  facet_wrap(~alpha_vals, nrow =2, scales="free")+
  theme_classic()


agg_dat_tr_subj = create_agg(df= data.exp3, 
                             form=formula(cbind(hd_norm, berries, grip_count)~ Env+tr_in_alpha+alpha_vals),
                             varlist = c('hd_norm', 'berries', 'grip_count'))

ggplot(data=agg_dat_tr_subj)+ geom_line(aes(x=tr_in_alpha, y= hd_norm, color=Env))+
  geom_errorbar(aes(x=tr_in_alpha, ymin = hd_norm_min, ymax =hd_norm_max, color=Env))+
  facet_wrap(~alpha_vals, nrow =1, scales="free")+
  theme_classic()

agg_dat_tr_prb = agg_dat_tr_subj %>% subset(alpha_vals==30)
plt_hd_prob<- ggplot(data=agg_dat_tr_prb)+ geom_line(aes(x=tr_in_alpha, y= hd_norm, color=Env))+
  geom_errorbar(aes(x=tr_in_alpha, ymin = hd_norm_min, ymax =hd_norm_max, color=Env),
                width=0.0)+
  facet_wrap(~alpha_vals, nrow =1, scales="free")+
  theme_classic()

ggsave('../plots/plt_hd_prob.pdf', plt_hd_prob, width=8, height=4)

agg_dat_tr_subj = create_agg(df= data.exp3, 
                             form=formula(cbind( berries, grip_count)~ Env+tr_in_alpha+alpha_vals),
                             varlist = c('berries', 'grip_count'), fn = median, fn_spread =st.err)

ggplot(data=agg_dat_tr_subj)+ geom_line(aes(x=tr_in_alpha, y= grip_count, color=Env))+
  geom_errorbar(aes(x=tr_in_alpha, ymin = grip_count_min, ymax =grip_count_max, color=Env))+
  facet_wrap(~alpha_vals, nrow =1, scales="free")+
  theme_classic()

data_single_subj = data.exp3 %>% subset(Subj_id == 321)

ggplot(data=data_single_subj)+ geom_line(aes(x=tr_in_alpha, y= grip_count, color=Env))+
  # geom_errorbar(aes(x=tr_in_alpha, ymin = grip_count_min, ymax =grip_count_max, color=Env))+
  facet_wrap(Subj_id~alpha_vals, ncol=3, scales="free")+
  theme_classic()

agg_dat_tr = create_agg(df= data.exp3, 
                             form=formula(cbind( berries, grip_count, hd_norm)~ Env+Trial),
                             varlist = c('berries', 'grip_count', 'hd_norm'), fn = mean, fn_spread =st.err)



ggplot(data=agg_dat_tr)+ geom_line(aes(x=Trial, y= grip_count, color=Env))+
  geom_errorbar(aes(x=Trial, ymin = grip_count_min, ymax =grip_count_max, color=Env))+
  # facet_wrap(Subj_id~alpha_vals, ncol=3, scales="free")+
  theme_classic()


ggplot(data=agg_dat_tr)+ geom_line(aes(x=Trial, y= hd_norm, color=Env))+
  geom_errorbar(aes(x=Trial, ymin = hd_norm_min, ymax =hd_norm_max, color=Env))+
  # facet_wrap(Subj_id~alpha_vals, ncol=3, scales="free")+
  theme_classic()
