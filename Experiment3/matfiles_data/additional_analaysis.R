# code to determine the window of the memory

require(ggplot2)
require(lme4)
require(sjPlot)

subj_data1$trial_scale= subj_data1$Trial/400
subj_data1$Probe <- as.factor(subj_data1$Probe)
subj_data1$Env <- as.factor(subj_data1$Env)
subj_data1$travel_scale = subj_data1$travel_dur/1000
subj_data1$move_scale = subj_data1$Move_dur/1000
subj_data1$Harvest_scale =subj_data1$Harvest_dur/1000
subj_data1$Reward_rate_scale = subj_data1$Reward_rate*1000
subj_data1$grip_norm = subj_data1$grip_force/30


# Split the probe window into halves

subj_data1$probe_blockno = numeric(nrow(subj_data1))

for (s in 1:9) {
    
  probe_blockno = subj_data1$probe_blockno[subj_data1$Subj==s]
  probe_blockno[21:30] = 1  
  probe_blockno[71:80] = 2  
  probe_blockno[121:130] = 3
  probe_blockno[171:180] = 4
  subj_data1$probe_blockno[subj_data1$Subj==s] = probe_blockno

}

eq_pv = Pk_Vel ~ Env*Probe + Env*Order + trial_scale + (1|Subj)
mdf_pv = lmer(eq_pv, data=subj_data1, REML=FALSE)
summary(mdf_pv)
ggplot()+geom_histogram(data=subj_data1, binwidth = 0.01, aes(x=Pk_Vel,fill=factor(Env)))+
  facet_wrap(~Subj, scales="free")

# Measure consistency of peak velocity within subjects

eq_harvest = Harvest_scale~Env*Probe+ Env*Order + Pk_Vel + trial_scale + (1|Subj)
mdf_harvest = lmer(eq_harvest, data=subj_data1)
summary(mdf_harvest)

# Reward rates discrepancy

ggplot()+geom_histogram(data=subj_data1, binwidth=.05, aes(x=Reward_rate_scale,fill=factor(Env)))+
  facet_wrap(~Subj, scales = "free")
eq_r = Reward_rate_scale~Env*Order + Env*Probe + Pk_Vel+ trial_scale+(1|Subj)
mdf_r = lmer( eq_r, data= subj_data1)
summary(mdf_r)
plot_model(mdf_r,type="diag")
plot_model(mdf_r, type="int")
ggplot()+geom_point(data=subj_data1,aes(x=factor(Env), y = Reward_rate_scale, color=factor(Order)))


# Average reward rate

reward_rate = numeric(18)
subj_arr = numeric(18)
rr_type = numeric(18)
order = numeric(18)

for (s in 1:9) {
  
  reward_rate_arr = subj_data1$avg_reward_rate[subj_data1$Subj==s]
  ord = subj_data1$Order[subj_data1$Subj==s]
  reward_rate[(s-1)*2+1] <- (reward_rate_arr[1]*1000)
  reward_rate[(s-1)*2+2] <- (reward_rate_arr[201]*1000)
  subj_arr[(s-1)*2+1] <- s
  subj_arr[(s-1)*2+2] <- s
  rr_type[(s-1)*2+1] <- 0
  rr_type[(s-1)*2+2] <- 1
  order[(s-1)*2+1] <- ord[1]
  order[(s-1)*2+2] <- ord[1]
    
}


s_inc = c(1,2,3,5,6,9)
s_dec= c(4,7,8)



ccc <- cbind(subj_arr, reward_rate, rr_type,order)
colnames(ccc) = c('subj', 'reward_rate','rr_type','order')
ccc <- as.data.frame(ccc)

ggplot()+geom_point(data=ccc, aes(x=factor(rr_type), y=reward_rate,color=factor(subj)))+
  geom_line(data=ccc, aes(x=factor(rr_type), y=reward_rate, group=subj,color=factor(subj), 
                          linetype=factor(order)))

ggplot()+geom_line(data=ccc, aes(x=rr_type, y=reward_rate,linetype=factor(order)))


# grip force normalized

a_gg = aggregate(grip_norm~Trial+Env, data=subj_data1, mean) 
b_gg = aggregate(grip_norm~Trial+Env, data=subj_data1, st.err)
c_gg <-  cbind( a_gg$Trial, a_gg$Env, a_gg$grip_norm, a_gg$grip_norm+ b_gg$grip_norm, a_gg$grip_norm - b_gg$grip_norm)
colnames(c_gg) = c("Trial", "Env", "grip_norm","ymax","ymin")
c_gg <- as.data.frame(c_gg) 

ggplot()+geom_ribbon(data=c_gg, aes(x=Trial, ymin= ymin, ymax=ymax, fill=factor(Env)),alpha=0.3)+
  geom_line(data=c_gg, aes(x=Trial, y= grip_norm,color=factor(Env) ))+
  annotate("rect", xmin=21,xmax=30,ymin=0.5, ymax=1.5,alpha=0.1,fill="red")+
  annotate("rect", xmin=71,xmax=80,ymin=0.5, ymax=1.5,alpha=0.1,fill="red")+
  annotate("rect", xmin=121,xmax=130,ymin=0.5, ymax=1.5,alpha=0.1,fill="red")+
  annotate("rect", xmin=171,xmax=180,ymin=0.5, ymax=1.5,alpha=0.1,fill="red")+
  ylab("Normalized Grip Force ")+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=20, hjust=0.5,face="bold"))


# grip force normalized

a_gg = aggregate(grip_norm~Trial+Env, data=subj_data1, mean) 
b_gg = aggregate(grip_norm~Trial+Env, data=subj_data1, st.err)
c_gg <-  cbind( a_gg$Trial, a_gg$Env, a_gg$grip_norm, a_gg$grip_norm+ b_gg$grip_norm, a_gg$grip_norm - b_gg$grip_norm)
colnames(c_gg) = c("Trial", "Env", "grip_norm","ymax","ymin")
c_gg <- as.data.frame(c_gg) 

ggplot()+geom_ribbon(data=c_gg, aes(x=Trial, ymin= ymin, ymax=ymax, fill=factor(Env)),alpha=0.3)+
  geom_line(data=c_gg, aes(x=Trial, y= grip_norm,color=factor(Env) ))+
  annotate("rect", xmin=21,xmax=30,ymin=0.5, ymax=1.5,alpha=0.1,fill="red")+
  annotate("rect", xmin=71,xmax=80,ymin=0.5, ymax=1.5,alpha=0.1,fill="red")+
  annotate("rect", xmin=121,xmax=130,ymin=0.5, ymax=1.5,alpha=0.1,fill="red")+
  annotate("rect", xmin=171,xmax=180,ymin=0.5, ymax=1.5,alpha=0.1,fill="red")+
  ylab("Normalized Grip Force ")+
  theme(axis.text=element_text(size=18), axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=20, hjust=0.5,face="bold"))+theme_classic()

eq_harvest = Harvest_scale~Env*Probe+ Env*Order +grip_norm + Pk_Vel + trial_scale + (1|Subj)
mdf_harvest = lmer(eq_harvest, data=subj_data1)
summary(mdf_harvest)

  
