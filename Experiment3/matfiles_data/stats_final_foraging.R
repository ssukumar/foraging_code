# Inlcude required libraries
library(lme4)
library(sjPlot)
library(ggplot2)
library(car)

# Stats for experiment 1 or 2 : 

data = subj_data1

# Repeated measures anova 


data$Env <- as.factor(data$Env)
data$Probe <- as.factor(data$Probe)
data$Order <- as.factor(data$Order)
data$trial_scale <- data$Trial/200
data$Trial <- as.factor(data$Trial)
data$Harvest_scale <- data$Harvest_dur/1000
data$probe_blockno <-as.factor(data$probe_blockno)
data$probeonly_blockno <-as.factor(data$probeonly_blockno)


# Peak Velocity / Vigor Results
# Repeated measures anova with factors Probe, Env, Order and trial

aov_agg = aggregate(data=data, Pk_Vel ~ Subj+Env+Probe+probe_blockno, mean )
pv_agg.aov <-with(data=aov_agg , aov(Pk_Vel~Env*Probe*probe_blockno + Error(Subj)))
summary(pv_agg.aov)

pv.aov <- with(data, aov(Pk_Vel ~ Env*Probe*Order +Trial + Error(Subj) ))
summary(pv.aov)


td.aov <- with(data, aov(travel_scale ~ Env*Probe*Order +Trial + Error(Subj) ))
summary(td.aov)

md.aov <- with(data, aov(move_scale ~ Env*Probe*Order +Trial + Error(Subj) ))
summary(md.aov)


pv2.aov <- with(data, aov(Pk_Vel ~ Env*probe_blockno*Probe + Order + Error(Subj) ) )
summary(pv2.aov)
# Lmer 

eq_pv = Pk_Vel ~ Env*Probe + probe_blockno*Probe+ Order  + (1|Subj)
mdf_pv = lmer(eq_pv, data=data)

plot_model(mdf_pv, type="diag")
vif(mdf_pv)
summary(mdf_pv)

eq_ts = travel_scale ~ Env*Probe + probe_blockno*Probe+ Order  + (1|Subj)
mdf_ts = lmer(eq_ts, data=data)

plot_model(mdf_ts,type="diag")
vif(mdf_ts)

summary(mdf_ts)

eq_pv2 = Pk_Vel ~ Env*Probe + probe_blockno*Probe + Order + (1+Env|Subj)
mdf_pv2 = lmer(eq_pv2 , data =data)

ss <-getME(mdf_pv2, c('theta', 'fixef'))
mdf_pv2.conv <- update(mdf_pv2, start=ss)

# diagnostics for the model
vif(mdf_pv2.conv)
plot_model(mdf_pv2.conv, type="diag")

summary(mdf_pv2.conv)

 # Harvest Duration results 
# Repeated measures anova

ggplot()+ geom_histogram(data=data, aes(x=Harvest_scale), binwidth =0.1)+facet_wrap(~Subj)
hs.aov <- with(data, aov(Harvest_scale ~ Env*Probe+ Order + Trial + Error(Subj) ))
summary(hs.aov)

hs1.aov <- with(data, aov(Harvest1_scale ~ Env*Probe*Order + Trial + Error(Subj) ))
summary(hs1.aov)


eq_har = Harvest_scale ~ Env*Order + Env*Probe + trial_scale + (1|Subj)
mdf_har = lmer(eq_har, data=data)

gmd_har = glmer(eq_har, family = Gamma("inverse"), data=data)

vif(mdf_har)
plot_model(mdf_har,type="diag")

summary(mdf_har)

# Grip force data (normalized )
# Repeated measures anova

grp.aov <- with(data, aov(grip_norm ~ Env*Order +Trial + Error(Subj)) )
summary(grp.aov)


rr.aov <- with(data, aov(rr~ Env*Probe*Order + probe_blockno + Error(Subj)))
summary(rr.aov)

eq_rr = rr~ Env*Probe*Order + poly(trial_scale,2)+ (1+Env|Subj)
mdf_rr = lmer(eq_rr, data=data)

vif(mdf_rr)
plot_model(mdf_rr,type="diag")

summary(mdf_rr)

eq_rr1 = rr~ Env*Probe*Order + poly(trial_scale,2)+ (1|Subj)
mdf_rr1 = lmer(eq_rr1, data=data)
vif(mdf_rr1)
plot_model(mdf_rr1, type="diag")

summary(mdf_rr1)

rrt.aov <- with(data, aov(rr_tot~ Env*Probe*Order + probe_blockno + Error(Subj)))
summary(rrt.aov)

eq_rrt = rr_tot~ Env*Probe*Order + poly(trial_scale,2)+ (1+Env|Subj)
mdf_rrt = lmer(eq_rrt, data=data)
ss <-getME(mdf_rrt, c('theta', 'fixef'))
mdf_rrt.conv <- update(mdf_rrt, start=ss)
vif(mdf_rrt.conv)
plot_model(mdf_rrt.conv,type="diag")

summary(mdf_rrt.conv)

eq_rrt1 = rr_tot~ Env*Probe*Order + poly(trial_scale,2)+ (1|Subj)

mdf_rrt1 = lmer(eq_rrt1, data=data)
vif(mdf_rrt1)
plot_model(mdf_rrt1, type="diag")

summary(mdf_rrt1)





