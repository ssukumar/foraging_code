

eq = Pk_Vel~Env*Order + Env*Probe+trial_scale+(1|Subj)
mdf = lmer( eq, data= subj_data1)
a = aggregate(Pk_Vel~Probe+Env+Subj,subj_data1, mean)
b = aggregate(Pk_Vel~Probe+Env+Subj,subj_data1, st.err)
cc = cbind(a$Subj, a$Probe,a$Env,a$Pk_Vel, a$Pk_Vel+b$Pk_Vel,a$Pk_Vel-b$Pk_Vel)
cc<-as.data.frame(cc)
head(cc)
colnames(cc) = c('Subj','Probe','Env','Pk_Vel','ymax','ymin')

cc<-as.data.frame(cc)
#g1=
  ggplot()+
  geom_line(data=cc, aes(x=Env,y=Pk_Vel,linetype=factor(Probe)), size=1.7)+
  scale_x_discrete(limits=c("1","2"),labels=c("Low","High"))+
  scale_linetype_discrete(name="", breaks = c(1,2),labels=c("Environment","Probe") )+
  geom_errorbar(data=cc, aes(x=Env, ymax = ymax, ymin=ymin),width=0.1, size=1.7)+
  facet_wrap(~Subj,scales="free")+
  xlab('Environment Effort Level')+ylab('Peak Velocity (m/s)')+
  theme(text=element_text(face="bold", size=18),
        plot.title = element_text(size=18, hjust=0.5,face="bold"), legend.key.size = unit(0.5,"in"),
        legend.text=element_text( size=18, face="bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
#ggsave(plot = g1,filename="/Volumes/Shruthi/CU/Research/Alaa_lab/foraging_data/facet1.pdf", useDingbats=FALSE)
ggexport(plot = g1,filename="/Volumes/Shruthi/CU/Research/Alaa_lab/foraging_data/facet1.pdf" , width = 15, height = 5, units = "cm")


a_tr1 = aggregate(travel_scale~Env+Probe+Subj, subj_data1,mean)
b_tr1 = aggregate(travel_scale~Env+Probe+Subj, subj_data1, st.err)
c_tr1 = cbind(a_tr1$Env, a_tr1$Probe,  a_tr1$travel_scale,a_tr1$Subj, a_tr1$travel_scale+b_tr1$travel_scale, 
              a_tr1$travel_scale-b_tr1$travel_scale)
head(c_tr1)
colnames(c_tr1)=c('Env','Probe','travel_scale','Subj','ymin','ymax')

c_tr1 <- as.data.frame(c_tr1)

ggplot()+geom_line(data=c_tr1, aes(x=Env, y=travel_scale,linetype=factor(Probe2)),size=1.7)+
  scale_x_discrete(limits=c("1","2"), labels=c("Low","High"))+
  scale_linetype_discrete(name = "", breaks= c(1,2),labels=c("Environment", "Probe2"))+
  geom_errorbar(data=c_tr1, aes(x=Env, ymax = ymax, ymin=ymin),width=0.1, size=1.7)+
  facet_wrap(~Subj,scales="free")+
  xlab('Environment Effort Level')+ylab("Travel Duration (s)")+
  ggtitle("Travel Duration")+
  theme(text=element_text(face="bold", size=30),
        axis.title.y= element_text(margin = margin(t = 0, r = 25, b = 0, l = 10)),
        axis.title.x= element_text(margin = margin(t = 25, r = 0, b = 10, l = 0)),
        plot.title = element_text(size=40, hjust=0.5,face="bold", margin = margin(t = 10, r = 0, b = 25, l = 0)), 
        legend.key.size = unit(0.5,"in"),
        legend.text=element_text( size=35, face="bold"), 
        legend.position = c(0.8,0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


a_mv1 = aggregate(move_scale~Env+Probe+Subj, subj_data1,mean)
b_mv1 = aggregate(move_scale~Env+Probe+Subj, subj_data1, st.err)
c_mv1 = cbind(a_mv1$Env, a_mv1$Probe,  a_mv1$move_scale,a_mv1$Subj, a_mv1$move_scale+b_mv1$move_scale, 
              a_mv1$move_scale-b_mv1$move_scale)
head(c_mv1)
colnames(c_mv1)=c('Env','Probe','move_scale','Subj','ymin','ymax')

c_mv1 <- as.data.frame(c_mv1)

ggplot()+geom_line(data=c_mv1, aes(x=Env, y=move_scale,linetype=factor(Probe)),size=1.7)+
  scale_x_discrete(limits=c("1","2"), labels=c("Low","High"))+
  scale_linetype_discrete(name = "", breaks= c(1,2),labels=c("Environment", "Probe2"))+
  geom_errorbar(data=c_mv1, aes(x=Env, ymax = ymax, ymin=ymin),width=0.1, size=1.7)+
  facet_wrap(~Subj,scales="free")+
  xlab('Environment Effort Level')+ylab("Travel Duration (s)")+
  ggtitle("Travel Duration")+
  theme(text=element_text(face="bold", size=30),
        axis.title.y= element_text(margin = margin(t = 0, r = 25, b = 0, l = 10)),
        axis.title.x= element_text(margin = margin(t = 25, r = 0, b = 10, l = 0)),
        plot.title = element_text(size=40, hjust=0.5,face="bold", margin = margin(t = 10, r = 0, b = 25, l = 0)), 
        # legend.key.size = unit(0.5,"in"),
        # legend.text=element_text( size=35, face="bold"), 
        # legend.position = c(0.8,0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


##################### New plot for probe env interaction with subject grouping ######################

a = aggregate(Pk_Vel~Probe+Env,subj_data1, mean)
b = aggregate(Pk_Vel~Probe+Env,subj_data1, st.err)
cc = cbind(a$Probe,a$Env,a$Pk_Vel, a$Pk_Vel+b$Pk_Vel,a$Pk_Vel-b$Pk_Vel)
cc<-as.data.frame(cc)
head(cc)
colnames(cc) = c('Probe','Env','Pk_Vel','ymax','ymin')
cc<-as.data.frame(cc)

cc = cc[cc$Probe==2,]

a1 = aggregate(Pk_Vel~Probe+Env+Subj,subj_data1, mean)
b1 = aggregate(Pk_Vel~Probe+Env+Subj,subj_data1, st.err)
cc1 = cbind(a1$Probe,a1$Env,a1$Pk_Vel, a1$Subj, a1$Pk_Vel+b1$Pk_Vel,a1$Pk_Vel-b1$Pk_Vel)
cc1<-as.data.frame(cc1)
head(cc1)
colnames(cc1) = c('Probe','Env', 'Pk_Vel', 'Subj','ymax','ymin')
cc1<-as.data.frame(cc1)

cc1=cc1[cc1$Probe==2,]


  ggplot()+
  geom_bar(data=cc, aes(x=Env,y=Pk_Vel),size=1,colour=c( "#F19998", "#D3FFB2"), fill='#9A97C8', stat = "identity")+
  geom_line(data=cc1, aes(x=Env,y=Pk_Vel,group = factor(Subj)),size=0.1)+
  scale_x_discrete(limits=c("1","2"),labels=c("Low","High"))+
  #scale_color_discrete(name="", breaks = c("1","2"),labels=c("Environment","Probe") )+
  geom_errorbar(data=cc, aes(x=Env, ymax = ymax, ymin=ymin),width=0.1, size=1.0)+
  geom_errorbar(data=cc1, aes(x=Env, ymax = ymax, ymin=ymin, group = factor(Subj)),width=0.05, size=.5)+
  xlab('Environment Effort')+ylab('Peak Velocity (m/s)')+scale_y_continuous(limits=c(0,1))+
  ggtitle('Peak Velocity')+theme_classic()+
  theme(text=element_text( face="bold", size=35),
        axis.title.y= element_text(margin = margin(t = 0, r = 25, b = 0, l = 10)),
        axis.title.x= element_text(margin = margin(t = 25, r = 0, b = 10, l = 0)),
        plot.title = element_text(size=40, hjust=0.5,face="bold", margin = margin(t = 10, r = 0, b = 25, l = 0)),
        # legend.key.size = unit(0.5,"in"),
        # legend.text=element_text( size=30, face="bold"),
        # legend.position = c(0.71,0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


  pattern.color=c("#F19998", "#D3FFB2")
  background.color=c('#9A97C8','#9A97C8')
  pattern.type<- c('crosshatch', 'crosshatch')
  density <- c(20,20)
m1 = patternbar(cc1,as.factor(cc1$Env), cc1$Pk_Vel , group=NULL,ylab='Low, High', pattern.type=pattern.type,
             pattern.color=pattern.color, background.color=background.color,
             pattern.line.size=0.5,frame.color=pattern.color)
  
m1 = m1+
  ggplot()+
  geom_line(data=cc1, aes(x=Env,y=Pk_Vel,group = factor(Subj)),size=0.1)+
  scale_x_discrete(limits=c("1","2"),labels=c("Low","High"))+
  #scale_color_discrete(name="", breaks = c("1","2"),labels=c("Environment","Probe") )+
  geom_errorbar(data=cc, aes(x=Env, ymax = ymax, ymin=ymin),width=0.1, size=1.0)+
  geom_errorbar(data=cc1, aes(x=Env, ymax = ymax, ymin=ymin, group = factor(Subj)),width=0.05, size=.5)+
  xlab('Environment Effort')+ylab('Peak Velocity (m/s)')+ scale_y_continuous(limits=c(0,1))+
  ggtitle('Peak Velocity')+theme_classic()+
  theme(text=element_text( face="bold", size=35),
        axis.title.y= element_text(margin = margin(t = 0, r = 25, b = 0, l = 10)),
        axis.title.x= element_text(margin = margin(t = 25, r = 0, b = 10, l = 0)),
        plot.title = element_text(size=40, hjust=0.5,face="bold", margin = margin(t = 10, r = 0, b = 25, l = 0)),
        # legend.key.size = unit(0.5,"in"),
        # legend.text=element_text( size=30, face="bold"),
        # legend.position = c(0.71,0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
ggexport(plot = m1,filename="/Volumes/Shruthi/CU/Research/Alaa_lab/foraging_data/facet1_opp.pdf" , width = 5, height = 5, units = "cm")

####################### New plot for probe env interaction with subject grouping  #########################

a_tr = aggregate(move_scale~Env+Probe, subj_data1,mean)
b_tr = aggregate(move_scale~Env+Probe, subj_data1, st.err)
c_tr = cbind(a_tr$Env, a_tr$Probe,  a_tr$move_scale,a_tr$move_scale+b_tr$move_scale, 
             a_tr$move_scale-b_tr$move_scale)
head(c_tr)
colnames(c_tr)=c('Env','Probe','move_scale','ymin','ymax')

c_tr <- as.data.frame(c_tr)

c_tr= c_tr[c_tr$Probe==2,]

a_tr1 = aggregate(move_scale~Env+Probe+Subj, subj_data1,mean)
b_tr1 = aggregate(move_scale~Env+Probe+Subj, subj_data1, st.err)
c_tr1 = cbind(a_tr1$Env, a_tr1$Probe,  a_tr1$move_scale,a_tr1$Subj, a_tr1$move_scale+b_tr1$move_scale, 
             a_tr1$move_scale-b_tr1$move_scale)
head(c_tr1)
colnames(c_tr1)=c('Env','Probe','move_scale','Subj','ymin','ymax')

c_tr1 <- as.data.frame(c_tr1)
c_tr1 = c_tr1[c_tr1$Probe==2,]

g2 = 
  ggplot()+geom_bar(data=c_tr, aes(x=Env, y=move_scale),size=1.2,colour=c( "#F19998", "#D3FFB2"),fill='#9A97C8',stat="identity")+
  geom_line(data = c_tr1, aes(x=Env, y=move_scale, group=factor(Subj)),size=0.1)+
  scale_x_discrete(limits=c("1","2"), labels=c("Low","High"))+
  #scale_linetype_discrete(name = "", breaks= c(1,2),labels=c("Environment", "Probe"))+
  geom_errorbar(data=c_tr, aes(x=Env, ymax = ymax, ymin=ymin),width=0.1, size=1.0)+
  geom_errorbar(data=c_tr1, aes(x=Env, ymax = ymax, ymin=ymin, group=factor(Subj)),width=0.05, size=.5)+
  xlab('Environment Effort Level')+ylab("Travel Duration (s)")+
  ggtitle("Travel Duration")+theme_classic()+
  theme(text=element_text(face="bold", size=30),
        axis.title.y= element_text(margin = margin(t = 0, r = 25, b = 0, l = 10)),
        axis.title.x= element_text(margin = margin(t = 25, r = 0, b = 10, l = 0)),
        plot.title = element_text(size=40, hjust=0.5,face="bold", margin = margin(t = 10, r = 0, b = 25, l = 0)), 
        legend.key.size = unit(0.5,"in"),
        legend.text=element_text( size=35, face="bold"), 
        legend.position = c(0.3,0.85),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# cc_tr = cbind(c_tr$Env[c_tr$Probe==1], c_tr$travel_scale[c_tr$Probe==1], a_tr$travel_scale[c_tr$Probe==1]+b_tr$travel_scale[c_tr$Probe==1],
#               a_tr$travel_scale[c_tr$Probe==1]-b_tr$travel_scale[c_tr$Probe==1])
# colnames(cc_tr) = c('Env','travel_scale','ymax','ymin')
# cc_tr <- as.data.frame(cc_tr)
  
  
# head(cc_tr)
ggexport(plot = g2,filename="/Volumes/Shruthi/CU/Research/Alaa_lab/foraging_data/facet2_opp.pdf" , width = 5, height = 5, units = "cm")


#Doing quick ANOVA on the data
a_tr1 = aggregate(move_scale~Env+Probe+Subj, subj_data1,mean)
b_tr1 = aggregate(move_scale~Env+Probe+Subj, subj_data1, st.err)
c_tr1 = cbind(a_tr1$Env, a_tr1$Probe,  a_tr1$move_scale,a_tr1$Subj, a_tr1$move_scale+b_tr1$move_scale, 
              a_tr1$move_scale-b_tr1$move_scale)
head(c_tr1)
colnames(c_tr1)=c('Env','Probe','move_scale','Subj','ymin','ymax')

c_tr1 <- as.data.frame(c_tr1)
aovv <- with(c_tr1, aov(move_scale~ Env*Probe +Error(Subj/(Env))))
summary(aovv)
