#LIBRARY###
library(lmerTest) 
library(emmeans)
library(doBy)
library(car)
library(Hmisc)
library(plyr)
library(xtable)
library(effects)
library(gridExtra)
library(plyr)
library(lubridate)
library(plotrix)
library(reshape2)
library(sjPlot)
library(jtools)
library(stringr)
library(stringi)
library(rvg)
library(officer)
library(ggsignif)
library(pmatch)
library(MuMIn)

p.match.v <- function(needle, haystack){
  pos=c()
  for (i in c(1:length(haystack))){
    pos=c(pos, pmatch(needle,haystack[i]))
  }
  which(pos==1)
}

outliersZ.3 <- function(data, zCutOff = 2, replace = NA, values = FALSE, digits = 2) {
  #compute standard deviation (sample version n = n [not n-1])
  stdev <- sqrt(sum((data - mean(data, na.rm = T))^2, na.rm = T) / sum(!is.na(data)))
  #compute absolute z values for each value
  absZ <- abs(data - mean(data, na.rm = T)) / stdev
  #subset data that has absZ greater than the zCutOff and replace them with replace
  #can also replace with other values (such as max/mean of data)
  data[absZ > zCutOff] <- replace 
  
  if (values == TRUE) {
    return(round(absZ, digits)) #if values == TRUE, return z score for each value
  } else {
    return(round(data, digits)) #otherwise, return values with outliers replaced
  }
}

dati="SET WORKING DIRECTORY"
pathOut="SET PATH FOR OUTPUT"

# Releveling: 
dati$age.group=factor(dati$age.group,levels(dati$age.group)[c(3,2,1)])
dati$source=factor(dati$source,levels(dati$source)[c(3,1,2)])
dati$Direction=relevel(dati$Direction, ref="Selfish")
dati$bin.cont.f=relevel(dati$bin.cont.f, ref="Not influenced")
dati$rt.don2[dati$rt.don2<250]="NA"

#Check
options("contrasts") #Should be 'treatment'
levels(dati$age.group) #Should be Young Adolescents, Mid Adolescents and Adults
levels(dati$source) #Should be Teenagers, Adults and Computer
levels(dati$Direction) #Should be Prosocial and Selfish
levels(dati$bin.cont.f) #Should be Influenced and Not Influenced

###############################################################################################
##########################    FIRST DONATIONS CATEGORICAL    ##################################
###############################################################################################

don1=lmer(don1~age.group + (1|sub), dati)
mod=data.frame(emmeans(don1, pairwise ~ age.group,adjust="Bonferroni", type="response")$emmeans)
table = summaryBy(don1 ~ sub + age.group, dati, FUN=c(mean))

ggplot(table, aes(x=age.group, y=don1.mean, fill=age.group)) +
  geom_violin(alpha=0.5)+theme_bw()+#geom_boxplot(width=0.08, fill='white')+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+
  ylab("First donation")+xlab("Age group")+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"), 
        panel.grid.minor = element_blank())+
  guides(fill=FALSE) +geom_point(data=mod, aes(age.group, emmean), shape=15, size=3)+
  geom_errorbar(data=mod, aes(age.group, emmean, ymin = asymp.LCL, ymax = asymp.UCL),width = 0.05,size  = 0.5)

ggsave(file.path(pathOut,paste0("don1.pdf")), width=12, height=8)

###############################################################################################
###########################    FIRST DONATIONS CONTINUOUS    ##################################
###############################################################################################

don1c.red=lmer(don1~poly(age,3) + (1|sub) , dati)
cont.bin.1=don1c.red
m.newdat <- expand.grid(
  age=unique(dati$age), 
  don1= 0
) 

m.mm <- model.matrix(terms(cont.bin.1),m.newdat) 
m.newdat$don1 <- m.mm %*% fixef(cont.bin.1) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(cont.bin.1),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = m.newdat$don1-2*sqrt(m.pvar1)
  , m.phi = m.newdat$don1+2*sqrt(m.pvar1)
)

names=names(fixef(cont.bin.1))
lin=-c(p.match.v("poly(age, 3)2",names),p.match.v("poly(age, 3)3",names))
# quad=-c(p.match.v("poly(age, 3)1",names),p.match.v("poly(age, 3)3",names))
cub=-c(p.match.v("poly(age, 3)1",names),p.match.v("poly(age, 3)2",names))
m.newdat$don1.1 <- m.mm[,lin] %*% fixef(cont.bin.1)[lin]
# m.newdat$rt.don2.2 <- exp(m.mm[,quad] %*% fixef(cont.bin.1)[quad])
m.newdat$don1.3 <- m.mm[,cub] %*% fixef(cont.bin.1)[cub]

table2 <- summaryBy(don1~ sub + age, dati, FUN=c(mean))
table2 <- summaryBy(don1.mean ~ age, table2, FUN=c(mean,length))

ggplot(table2, aes(age, don1.mean.mean))+
  geom_point(aes(size=don1.mean.length), alpha=0.5)+
  scale_size_area(name = "Number of participants")+
  xlab("Age")+ylab("")+geom_line(data=m.newdat,aes(age, don1), size=1) + 
  geom_ribbon(data=m.newdat, aes(age,don1, ymin=m.plo, ymax=m.phi),alpha=0.2)+theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=22,face="bold"),
        panel.grid.minor = element_blank())+ylab("First donation")+
  geom_line(data=m.newdat, aes(age, don1.1, color="linear"), lty=2)+
  geom_line(data=m.newdat, aes(age, don1.3, color="cubic"), lty=2, size=0.4)

ggsave(file.path(pathOut,paste0("don1.cont.pdf")), width=12, height=8)

###############################################################################################
########################   INFLUENCE PROBABILITY CATEGORICAL   ################################
###############################################################################################

cont.bin.1.red=glmer(bin.cont~age.group*Direction+source + (Direction + source|sub),
                     control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
                     family="binomial", dati)
#Source----
table = summaryBy(bin.cont ~ sub + age.group + source, dati, FUN=c(mean))
mod=data.frame(emmeans(cont.bin.1.red, pairwise ~ source,adjust="Bonferroni", type="response")$emmeans)

ggplot(table, aes(x=source, y=bin.cont.mean, fill=source)) +
  geom_violin(alpha=0.5)+theme_bw()+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+
  ylab("Influence probability")+xlab("Source")+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"), 
        panel.grid.minor = element_blank())+
  guides(fill=FALSE) + geom_point(data=mod, aes(source, prob), shape=15, size=3)+
  geom_errorbar(data=mod, aes(source, prob, ymin = asymp.LCL, ymax = asymp.UCL),width = 0.05,size  = 0.5)

ggsave(file.path(pathOut,paste0("Source.pdf")), width=12, height=8)

# Age (Main Effect)----
table = summaryBy(bin.cont ~ sub + age.group, dati, FUN=c(mean))
mod=data.frame(emmeans(cont.bin.1.red, pairwise ~ age.group,adjust="Bonferroni", type="response")$emmeans)

ggplot(table, aes(x=age.group, y=bin.cont.mean, fill=age.group)) +
  geom_violin(alpha=0.5)+theme_bw()+#geom_boxplot(width=0.08, fill='white')+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+
  ylab("Influence probability")+xlab("Age group")+
  theme(axis.text=element_text(size=17),
        axis.title=element_text(size=22,face="bold"), 
        panel.grid.minor = element_blank())+
  guides(fill=FALSE) + geom_point(data=mod, aes(age.group, prob), shape=15, size=3)+
  geom_errorbar(data=mod, aes(age.group, prob, ymin = asymp.LCL, ymax = asymp.UCL),width = 0.05,size  = 0.5)

ggsave(file.path(pathOut,paste0("AgeMain.pdf")), width=12, height=8)

#Interaction: Age x Direction----
table = summaryBy(bin.cont ~ sub + age.group + Direction, dati, FUN=c(mean))
mod=data.frame(emmeans(cont.bin.1.red, pairwise ~ age.group | Direction,adjust="Bonferroni", type="response")$emmeans)

ggplot(table, aes(x=age.group, y=bin.cont.mean, fill=age.group)) +
  geom_violin(alpha=0.5)+theme_bw() +#geom_boxplot(width=0.08, fill='white')+
  geom_point(position=position_jitter(width=0.1, height=0), size=1, alpha=0.2)+
  ylab("Influence probability")+xlab("Age group")+
  theme(axis.text=element_text(size=17),
        axis.title=element_text(size=22,face="bold"), 
        panel.grid.minor = element_blank())+
  guides(fill=FALSE) + geom_point(data=mod, aes(age.group, prob), shape=15, size=3)+
  geom_errorbar(data=mod, aes(age.group, prob, ymin = asymp.LCL, ymax = asymp.UCL),width = 0.05,size  = 0.5)+
  facet_wrap(~Direction)

ggsave(file.path(pathOut,paste0("AgexDirection.pdf")), width=12, height=8)

###############################################################################################
#########################   INFLUENCE PROBABILITY CONTINUOUS   ################################
###############################################################################################

IPc.1red=glmer(bin.cont~I(age^-1)*Direction+source+ (Direction + source|sub),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               family="binomial", dati)
#Age Trend----
cont.bin.1=IPc.1red
m.newdat <- expand.grid(
  source=levels(dati$source),
  Direction = levels(dati$Direction),
  age=unique(dati$age), 
  bin.cont= 0
) 

m.mm <- model.matrix(terms(cont.bin.1),m.newdat) 
m.newdat$bin.cont <- m.mm %*% fixef(cont.bin.1) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(cont.bin.1),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = plogis(m.newdat$bin.cont-2*sqrt(m.pvar1))
  , m.phi = plogis(m.newdat$bin.cont+2*sqrt(m.pvar1))
)
m.newdat$bin.cont=plogis(m.newdat$bin.cont)
m.newdat2=aggregate(.~age, m.newdat[,-c(1,2)], mean)

colnames(m.newdat2)[c(2)]=c("bin.cont")

table <- summaryBy(bin.cont ~ sub + age + Direction, dati, FUN=c(mean))
table <- summaryBy(bin.cont.mean ~ age, table, FUN=c(mean,length))

ggplot(table, aes(age, bin.cont.mean.mean))+
  geom_point(aes(size=bin.cont.mean.length), alpha=0.6)+
  scale_size_area(name = "Number of participants")+
  xlab("Age")+ylab("")+geom_line(data=m.newdat2, aes(age, bin.cont), size=0.7) + 
  geom_ribbon(data=m.newdat2, aes(age,bin.cont, ymin=m.plo, ymax=m.phi),alpha=0.2)+theme_bw()+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"),
        panel.grid.minor = element_blank())+ylab("Influence probability")

ggsave(file.path(pathOut,paste0("Age.cont.pdf")), width=12, height=8)

#Interaction: Age Trend x Direction----
cont.bin.1=IPc.1red
m.newdat <- expand.grid(
  source=levels(dati$source),
  Direction = levels(dati$Direction),
  age=unique(dati$age), 
  bin.cont= 0
) 

m.mm <- model.matrix(terms(cont.bin.1),m.newdat) 
m.newdat$bin.cont <- m.mm %*% fixef(cont.bin.1) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(cont.bin.1),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = plogis(m.newdat$bin.cont-2*sqrt(m.pvar1))
  , m.phi = plogis(m.newdat$bin.cont+2*sqrt(m.pvar1))
)
m.newdat$bin.cont=plogis(m.newdat$bin.cont)

m.newdat2=aggregate(.~age+Direction, m.newdat[,-c(1)], mean)

colnames(m.newdat2)[c(3)]=c("bin.cont")

table <- summaryBy(bin.cont ~ sub + age + Direction, dati, FUN=c(mean))
table <- summaryBy(bin.cont.mean ~  Direction + age, table, FUN=c(mean,length))

ggplot(table, aes(age, bin.cont.mean.mean, group=Direction, fill=Direction))+
  geom_point(aes(size=bin.cont.mean.length, color=Direction), alpha=0.6)+
  scale_size_area(name = "Number of participants")+
  xlab("Age")+ylab("")+geom_line(data=m.newdat2, aes(age, bin.cont, color=Direction), size=0.7) + 
  geom_ribbon(data=m.newdat2, aes(age,bin.cont, ymin=m.plo, ymax=m.phi),alpha=0.2)+theme_bw()+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"),
        panel.grid.minor = element_blank())+ylab("Influence probability")

ggsave(file.path(pathOut,paste0(Sys.Date(),"AgexDirection.cont.pdf")), width=12, height=8)


###############################################################################################
##########################     REACTION TIME CATEGORICAL     ##################################
###############################################################################################

rt3.red.2=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ (Direction + source + bin.cont.f|sub),
               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati) 

mod=data.frame(emmeans(rt3.red.2, pairwise ~ age.group|Direction*bin.cont.f,adjust="none", type="response")$emmeans)
ggplot(mod,  aes(age.group, response, group=interaction(age.group:bin.cont.f), color=bin.cont.f)) +
  geom_point(shape=15, size=3, position = position_dodge(width=0.50))+
  geom_errorbar(aes(min = asymp.LCL, ymax = asymp.UCL), width = 0.05,size  = 0.5,position = position_dodge(width=0.50))+
  theme_bw()+ylab("RT (ms) second donation")+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"), 
        axis.text.x=element_text(angle=45, hjust=1),legend.text=element_text(size=18),
        legend.title=element_text(size=18),strip.text = element_text(size=12, face="bold"),
        panel.grid.minor = element_blank())+xlab("")+facet_wrap(~Direction)

ggsave(file.path(pathOut,paste0("RT.pdf")), width=12, height=8)

###############################################################################################
###########################     REACTION TIME CONTINUOUS     ##################################
###############################################################################################

RTc.3red=lmer(log(rt.don2)~(poly(age,2)+bin.cont.f+Direction)^3+
                +source + (Direction+source+bin.cont.f|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)

cont.bin.1=RTc.3red
m.newdat <- expand.grid(
  source=levels(dati$source),
  Direction = levels(dati$Direction),
  age=unique(dati$age), 
  bin.cont.f=levels(dati$bin.cont.f),
  rt.don2= 0
) 

m.mm <- model.matrix(terms(cont.bin.1),m.newdat) 
m.newdat$rt.don2 <- m.mm %*% fixef(cont.bin.1) 
m.pvar1 <- diag(m.mm %*% tcrossprod(vcov(cont.bin.1),m.mm)) 
m.newdat <- data.frame(
  m.newdat
  , m.plo = exp(m.newdat$rt.don2-2*sqrt(m.pvar1))
  , m.phi = exp(m.newdat$rt.don2+2*sqrt(m.pvar1))
  #, m.tlo = m.newdat$respC-2*sqrt(m.tvar1)
  #, m.thi = m.newdat$respC+2*sqrt(m.tvar1)
)
m.newdat$rt.don2=exp(m.newdat$rt.don2)

names=names(fixef(cont.bin.1))
lin=-c(p.match.v("poly(age, 2)2",names),p.match.v("poly(age, 2)3",names))
quad=-c(p.match.v("poly(age, 2)1",names),p.match.v("poly(age, 2)3",names))
# cub=-c(p.match.v("poly(age, 3)1",names),p.match.v("poly(age, 3)2",names))
m.newdat$rt.don2.1 <- exp(m.mm[,lin] %*% fixef(cont.bin.1)[lin])
m.newdat$rt.don2.2 <- exp(m.mm[,quad] %*% fixef(cont.bin.1)[quad])
# m.newdat$rt.don2.3 <- exp(m.mm[,cub] %*% fixef(cont.bin.1)[cub])

m.newdat3=aggregate(.~age+Direction + bin.cont.f, m.newdat[,-1], mean)
colnames(m.newdat3)[c(4,7,8)]=c("rt","rt.1","rt.2")

table2 <- summaryBy(rt.don2 ~ sub + age + Direction + bin.cont.f, dati, FUN=c(median))
table2 <- summaryBy(rt.don2.median ~ Direction + age + bin.cont.f, table2, FUN=c(median,length))

# Suppressing linear and quadratic estimates for trends that did not interact with the influence term. 
m.newdat3$rt.1[!(m.newdat3$Direction=="Prosocial" & m.newdat3$bin.cont.f=="Not influenced")]=NA
m.newdat3$rt.2[!(m.newdat3$Direction=="Selfish" & m.newdat3$bin.cont.f=="Influenced")]=NA

ggplot(table2, aes(age, rt.don2.median.median, group=bin.cont.f,fill=bin.cont.f))+
  geom_point(aes(size=rt.don2.median.length, color=bin.cont.f), alpha=0.5)+
  scale_size_area(name = "Number of participants")+
  xlab("Age")+ylab("")+geom_line(data=m.newdat3,aes(age, rt, color=bin.cont.f), size=1) + 
  geom_ribbon(data=m.newdat3, aes(age,rt, ymin=m.plo, ymax=m.phi, fill=bin.cont.f),alpha=0.2)+theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=22,face="bold"),
        panel.grid.minor = element_blank())+ylab("RT (ms)")+facet_wrap(~Direction)+
  geom_line(data=m.newdat3, aes(age, rt.1,color=bin.cont.f), lty=2)+geom_line(data=m.newdat3, aes(age, rt.2,color=bin.cont.f), lty=2)

ggsave(file.path(pathOut,paste0("RT.cont.pdf")), width=12, height=8)

###############################################################################################
########################    INFLUENCE MAGNITUDE CATEGORICAL    ################################
###############################################################################################

dati$cont.2=dati$cont
dati$cont.2[dati$cont.2 < 0] = 0 

IM.1red=lmer(cont.2~ age.group*Direction*source +scale(delta.abs)*Direction*age.group + (Direction+source+scale(delta.abs)|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)

slopes=data.frame(emtrends(IM.1red, pairwise ~ Direction | age.group, var= "delta.abs", type="response", adjust='Bonferroni')$emtrends)

ggplot(slopes,  aes(age.group, delta.abs.trend, group=Direction, fill=Direction, color=Direction)) +
  geom_point(shape=15, size=3,position = position_dodge(width=0.40))+
  geom_errorbar(aes(min = asymp.LCL, ymax = asymp.UCL), width = 0.05,size  = 0.5,position = position_dodge(width=0.40))+
  theme_bw()+ylab("Slope for change in donations")+geom_bar(stat="identity", alpha=0.6, width=0.4, position=position_dodge())+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22,face="bold"), 
        axis.text.x=element_text(angle=45, hjust=1),legend.text=element_text(size=18),
        legend.title=element_text(size=18),strip.text = element_text(size=12, face="bold"),
        panel.grid.minor = element_blank())+xlab("")

ggsave(file.path(pathOut,paste0("IM.pdf")), width=12, height=8)
