
#"Prosocial influence and opportunistic conformity in adolescents and young adults". Psychological Science (in press)
#By Gabriele Chierchia, Blanca Piera Pi-Sunyer and Sarah-Jayne Blakemore

#Institute of Cognitive Neuroscience, University College London

# Data is available in the open science framework at: 
# https://osf.io/3e9s6/

#LOAD PACKAGES ####
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
library(dplyr)
library(data.table)
library(formattable)
library(tidyr)
library(htmltools)
library(webshot)    
library(outliers)
library(sessioninfo)

# Session Information####

sessionInfo()

# USER DEFINED FUNCTIONS####

# To later create tables
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

num_lev <- function(x) {
  length(unique(x))
}

# Replaces values higher/lower than 4*SD from the mean (we later re-run models
# without such extreme values)
outliersZ.4 <- function(data, zCutOff = 4, replace = NA, values = FALSE, digits = 2) {
  stdev <- sqrt(sum((data - mean(data, na.rm = T))^2, na.rm = T) / sum(!is.na(data)))
  absZ <- abs(data - mean(data, na.rm = T)) / stdev
  data[absZ > zCutOff] <- replace 
  
  if (values == TRUE) {
    return(round(absZ, digits))
  } else {
    return(round(data, digits))
  }
}

# Counts number of observations per level of factor
obs.per.lev <- function(x) {
  paste0(levels(x)[1],"=", length(which(x==levels(x)[1]))," ", levels(x)[2],"=",length(which(x==levels(x)[2])))
}

#PREPROCESSING ####
# Load data:
dati=read.csv("INSERT HERE DATA FILE IN DIRECTORY")

#Setting Contrasts
options(contrasts=c('contr.treatment','contr.poly')) # (Also the R default)

# Adding 8 column transformation to the dataset:

# # We relabel the factor levels for the influence factor and store this in
# "bin.cont.f"
dati$bin.cont.f=revalue(as.factor(dati$bin.cont), c("1"="Influenced", "0"="Not influenced"))

# We store separate components for the linear, quadratic and cubic components of
# the "age" term (this is only needed to obtain p-values of the omnibus tests
# for the separate components using the Anova function)
dati$age.1=poly(dati$age,3)[,1]
dati$age.2=poly(dati$age,3)[,2]
dati$age.3=poly(dati$age,3)[,3]

# Adding the delta converted to monetary value
dati$delta.conv= dati$delta.abs*dati$guess.conversion.rate

# Adding 4 columns related to computer trials only, 2 for influence magnitude
# and 2 for influence probability: on vector in each pair aggregates over all
# computer trials, the second breaks these down by direction.
# Aggregate computer contagion at the sub x Direction level:   
comp.dir <- summaryBy(bin.cont + cont ~ sub + Direction, droplevels(dati[dati$source=="Computer",]), FUN=c(mean))
colnames(comp.dir)[c(3,4)]=c("comp.bin.cont.dir", "comp.cont.dir")
# Aggregate computer contagion at the sub level:    
comp <- summaryBy(bin.cont + cont ~ sub, droplevels(dati[dati$source=="Computer",]), FUN=c(mean))
colnames(comp)[c(2,3)]=c("comp.bin.cont", "comp.cont")

dati=merge(dati, comp.dir, by=c("sub", "Direction"), all.x=T)
dati=merge(dati, comp, by=c("sub"), all.x=T)

#Setting reference levels of each factor
dati$age.group=factor(dati$age.group,levels(dati$age.group)[c(3,2,1)])
dati$source=factor(dati$source,levels(dati$source)[c(3,1,2)])
dati$Direction=relevel(dati$Direction, ref="Prosocial")
dati$bin.cont.f=relevel(dati$bin.cont.f, ref="Influenced")

dati$bin.cont.anti=NA
dati$bin.cont.anti[dati$cont >= 0]=0
dati$bin.cont.anti[dati$cont < 0]=1

dati$change=0
dati$change[dati$cont!=0]=1

dati$c.ac[dati$cont>0]=1
dati$c.ac[dati$cont==0]=0
dati$c.ac[dati$cont<0]=-1

var<- summaryBy(c.ac~ sub, dati, FUN=c(var))
dati=merge(dati,var, all.x=T)

dati$conf='no.change'
dati$conf[dati$cont>0]='conf'
dati$conf[dati$cont<0]='anti.conf'
dati$conf=as.factor(dati$conf)

#Check
options("contrasts") #Should be 'treatment'
levels(dati$age.group) #Should be Young Adolescents, Mid Adolescents and Adults
levels(dati$source) #Should be Teenagers, Adults and Computer
levels(dati$Direction) #Should be Prosocial and Selfish
levels(dati$bin.cont.f) #Should be Influenced and Not Influenced


##########################    FIRST DONATIONS CATEGORICAL    ##################################

don1=lmer(don1~age.group+ (1|sub), dati)
Anova(don1, type=3)

#Contrasts --> Age
cons=as.data.frame(emmeans(don1, pairwise ~ age.group,adjust="none")$contrasts)
colnames(cons)[6]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#RM#
#Suspicion
don1.susp=lmer(don1~age.group  + (1|sub) , droplevels(dati[dati$evaluation =="Not Problematic",]))
Anova(don1.susp, type=3)

#Special Needs
don1.sn=lmer(don1~age.group  + (1|sub) , droplevels(dati[dati$sn =="none",]))
Anova(don1.sn, type=3)

#No ceiling no floor
don1.ncnf=lmer(don1~age.group  + (1|sub) , droplevels(dati[!dati$cat %in% c("floor", "ceil" ,"no.var"),]))
Anova(don1.ncnf, type=3)

#Defect
don1.defect=lmer(don1~age.group  + (1|sub) , droplevels(dati[dati$eval =="ok",]))
Anova(don1.defect, type=3)

#Only Adolescents
don1.adol=lmer(don1~age.group + (1|sub), droplevels(dati[dati$age.group %in% c("Young adolescents","Mid adolescents"),]))
Anova(don1.adol, type=3)

#Extreme Values
dati$don1.no.out=outliersZ.4(dati$don1)
don1.n.no.out=lmer(don1.no.out~age.group  + (1|sub) , dati)
Anova(don1.n.no.out, type=3)

#CM#
#Gender
don1.gen=lmer(don1~age.group*gender  + (1|sub) , dati)
Anova(don1.gen, type=3)
emmeans(don1.gen, pairwise ~ gender, adjust="Bonferroni")

#GCR
don1.gcr=lmer(don1~age.group + scale(guess.conversion.rate) + (1|sub) ,dati)
Anova(don1.gcr, type=3)

#ART
don1.art=lmer(don1~age.group +scale(propCorr) + (1|sub) , dati)
Anova(don1.art, type=3)

#Group Size
don1.gs=lmer(don1~age.group + scale(group.size) + (1|sub) ,dati)
Anova(don1.gs,type=3)

# Do Deltas vary with Age?
deltas=lmer(delta.abs~Direction*age.group + (Direction|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(deltas, type=3)

cons=as.data.frame(emmeans(deltas, pairwise ~ age.group|Direction,adjust="none")$contrasts)
colnames(cons)[7]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]
cons$p.bonf[cons$p.bonf>1]=1

cons=as.data.frame(emmeans(deltas, pairwise ~ Direction|age.group,adjust="none")$contrasts)
colnames(cons)[7]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]
cons$p.bonf[cons$p.bonf>1]=1

##########################    FIRST DONATIONS CONTINUOUS     ##################################

# Run models with different curve fitting techniques: 
mm.1=lmer(don1~poly(age,1) + (1|sub), dati)
mm.1.exp=lmer(don1~exp(age) + (1|sub), dati)
mm.1.log=lmer(don1~log(age) + (1|sub), dati)
mm.1.inv=lmer(don1~I(age^-1) + (1|sub), dati)
mm.2=lmer(don1~I(age^2) + (1|sub), dati)
mm.3=lmer(don1~I(age^3) + (1|sub), dati)
mm.1.2=lmer(don1~poly(age,2) + (1|sub), dati)
mm.1.3=lmer(don1~age+I(age^3) + (1|sub), dati)
mm.1.2.3=lmer(don1~poly(age,3) + (1|sub), dati)

# Compare model and sort them according to AIC
model.list=mget(grep("mm.", ls(),value=T))
aics<-lapply(model.list,function(x)AIC(x))
sort(unlist(aics))

#Basic Model
don1c=lmer(don1~poly(age,3) + (1|sub), dati)
Anova(don1c,type=3)

don1c.red=lmer(don1~age.1+age.2+age.3 + (1|sub) , dati)
Anova(don1c.red,type=3)

#RM#
#Suspicion
don1c.susp=lmer(don1~age.1+age.2+age.3+(1|sub), droplevels(dati[dati$evaluation =="Not Problematic",]))
Anova(don1c.susp, type=3)

#Special Needs
don1c.sn=lmer(don1~age.1+age.2+age.3+ (1|sub), droplevels(dati[dati$sn =="none",]))
Anova(don1c.sn, type=3)

#No ceiling no floor
don1c.ncnf=lmer(don1~age.1+age.2+age.3 +(1|sub), droplevels(dati[!dati$cat %in% c("floor", "ceil" ,"no.var"),]))
Anova(don1c.ncnf,type=3)

#Defect
don1c.defect=lmer(don1~age.1+age.2+age.3 +(1|sub), droplevels(dati[dati$eval=="ok",]))
Anova(don1c.defect,type=3)

#Only Adolescents
don1c.adol=lmer(don1~age.1+age.2+age.3 + (1|sub) , droplevels(dati[dati$age.group %in% c("Young adolescents","Mid adolescents"),]))
Anova(don1c.adol,type=3)

#Extreme Values
dati$don1c.no.out=outliersZ.4(dati$don1)
don1c.n.no.out=lmer(don1c.no.out~age.1+age.2+age.3  + (1|sub), dati)
Anova(don1c.n.no.out, type=3)

#CM#
#Gender
don1c.gen=lmer(don1~age.1+age.2+age.3 + age.1*gender + age.2*gender + age.3*gender  + (1|sub), dati)
Anova(don1c.gen, type=3)

#GCR
don1c.gcr=lmer(don1~age.1+age.2+age.3 +scale(guess.conversion.rate)+ (1|sub),dati)
Anova(don1c.gcr, type=3)

#ART
don1c.art=lmer(don1~age.1+age.2+age.3 +scale(propCorr) +(1|sub), dati)
Anova(don1c.art, type=3)

#Group Size
don1c.gs=lmer(don1~age.1+age.2+age.3 + scale(group.size) +(1|sub), dati)
Anova(don1c.gs, type=3)

##############################    MANIPULATION CHECKS    ######################################

# Influence frequency manipulation check----
#Adding change column: Whether participants adjusted their donations from Phase
#1 to Phase 2 or not.
dati$change[dati$cont==0]=0
dati$change[dati$cont!=0]=1

# What % of trials did participants change their responses
round(mean(dati$change),2)

# What % of that % did participants conform as opposed to anti-conform?
table <- summaryBy(bin.cont ~ sub, droplevels(dati[dati$change==1,]), FUN=c(mean))
100*round(mean(dati$bin.cont[dati$change==1]),2)

# Is this proportion different from 0.5? 
with(droplevels(dati[dati$change==1,]), 
     binom.test(length(which(bin.cont==1)), length(bin.cont),
                p = 0.5))

# And for each condition separately? 
m <- expand.grid(
  age.group=levels(dati$age.group), 
  source = levels(dati$source),
  Direction = levels(dati$Direction)
) 

for (i in c(1:nrow(m))){
  b=with(droplevels(dati[dati$change==1 & dati$age.group==m$age.group[i] & dati$Direction == m$Direction[i] & dati$source==m$source[i],]), 
         binom.test(length(which(bin.cont==1)), length(bin.cont), p = 0.5))
  m$p.unc[i]=b$p.value
}

m$p.bonf=round(m$p.unc*nrow(m),3)
m$p.unc=round(m$p.unc,3)
m$p.bonf[m$p.bonf>1]=1


# Influence magnitude manipulation check----
table <- summaryBy(cont ~ sub + age.group, dati, FUN=c(mean))
t=t.test(table$cont.mean) 
round(t$estimate,2)
c(round(t$conf.int[1],2), round(t$conf.int[2],2))

table <- summaryBy(cont ~ sub + age.group + source + Direction, dati, FUN=c(mean))
m <- expand.grid(
  age.group=levels(dati$age.group), 
  source = levels(dati$source),
  Direction = levels(dati$Direction)
) 

for (i in c(1:nrow(m))){
  t=t.test(table$cont.mean[table$age.group == m$age.group[i] & table$Direction == m$Direction[i] & table$source ==m$source[i]])
  m$mean[i]=round(t$estimate,3)
  m$CI.lower[i]=as.numeric(round(t$conf.int,3))[1]
  m$CI.upper[i]=as.numeric(round(t$conf.int,3))[2]
  m$t[i]=t$statistic
  m$Df[i]=t$parameter
  m$p.unc[i]=t$p.value
}
m$p.bonf=round(m$p.unc*nrow(m),3)
m$p.unc=round(m$p.unc,3)
m$t=round(m$t,2)
m$p.bonf[m$p.bonf>1]=1

m

# Are age effects of conformity and anti-conformity robust to controlling for response variance ----
var =merge(var, unique(dati[,c('sub', 'age')]), all.x=T)
var.age = lm(c.ac.var ~ I(age^-1), var)
Anova(var.age, type=3)
summary(var.age) # Variance increases with the inverse of age. 

anti<- summaryBy(bin.cont.anti~ sub, dati, FUN=c(mean)) # 
dati=merge(dati,anti, all.x=T)


IP.c.base=glmer(bin.cont~I(age^-1) + (1|sub),
                control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
                family="binomial",dati)
Anova(IP.c.base, type=3) # conforming probability decreases with age. 
summary(IP.c.base)

IP.c.var.base=glmer(bin.cont~I(age^-1) + c.ac.var+ (1|sub),
                    control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
                    family="binomial",dati)
Anova(IP.c.var.base, type=3) # this is robust to controlling for response variance
summary(IP.c.var.base)

IP.anti.c=glmer(bin.cont.anti~I(age^-1) +  (1|sub),
                control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
                family="binomial",dati)
Anova(IP.anti.c, type=3)
summary(IP.anti.c) # Anticonformting probability also increases with the inverse of age. 

IP.anti.c.var=glmer(bin.cont.anti~I(age^-1) + c.ac.var + (1|sub),
                    control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
                    family="binomial",dati)
Anova(IP.anti.c.var, type=3) #Completely wiped out by response variance. 
summary(IP.anti.c.var) 

# Non-linear deltas?----
dati$delta.abs.s=scale(dati$delta.abs)
# Full dataset
dati.t = dati
# If I don't scale model doesn't converge
mm.1=lmer(cont~delta.abs.s + (delta.abs.s|sub),
          control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati.t)
summary(mm.1) # Shows significant linear effect of delta on cont
mm.1.gl=glmer(bin.cont~delta.abs.s + (1|sub),family = 'binomial',
              control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati.t)
summary(mm.1.gl) # same for bin.cont

df = data.frame(var = c(rep('IM',6), rep('IP',6)), direction = expand.grid(levels(dati$Direction), levels(dati$age.group))[,1], 
                age.group = expand.grid(levels(dati$Direction), levels(dati$age.group))[,2], best.fit = NA, lin.slope = NA, p = NA)

for (i in c(1:nrow(df))){
  dati.t=dati[dati$age.group == df$age.group[i] & dati$Direction == df$direction[i],]
  # IM
  mm.im.1=lmer(cont~delta.abs.s + (1|sub),
               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati.t)
  mm.im.2=lmer(cont~I(delta.abs.s^2) + (1|sub),
               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati.t)
  mm.im.3=lmer(cont~I(delta.abs.s^3) + (1|sub),
               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati.t)
  mm.im.4=lmer(cont~I(delta.abs.s^4) + (1|sub),
               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati.t)
  # IP
  mm.ip.1.gl=glmer(bin.cont~delta.abs.s + (1|sub),family = 'binomial',
                   control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati.t)
  mm.ip.2.gl=glmer(bin.cont~I(delta.abs.s^2) + (1|sub),family = 'binomial',
                   control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati.t)
  mm.ip.3.gl=glmer(bin.cont~I(delta.abs.s^3) + (1|sub),family = 'binomial',
                   control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati.t)
  mm.ip.4.gl=glmer(bin.cont~I(delta.abs.s^4) + (1|sub),family = 'binomial',
                   control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati.t)
  
  model.list.im=mget(grep("mm.im", ls(),value=T))
  aics.im<-names(sort(unlist(lapply(model.list.im,function(x)AIC(x))))[1])
  model.list.ip=mget(grep("mm.ip", ls(),value=T))
  aics.ip<-names(sort(unlist(lapply(model.list.ip,function(x)AIC(x))))[1])
  
  if (df$var[i]=="IM"){
    df$best.fit[i] = aics.im
    df$lin.slope[i]=round(fixef(mm.im.1)[2],3)
    df$p[i]=round(summary(mm.im.1)$coefficients[2,5],3)
  }
  if (df$var[i]=="IP"){
    df$best.fit[i] = aics.ip
    df$lin.slope[i]=round(fixef(mm.ip.1.gl)[2],3)
    df$p[i]=round(summary(mm.ip.1.gl)$coefficients[2,4],3)
  }
}
df 

########################   INFLUENCE PROBABILITY CATEGORICAL   ################################

IP.1=glmer(bin.cont~age.group*Direction*source + (Direction + source|sub),
           control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
           family="binomial", dati)
Anova(IP.1,type=3)

#Model reduction
IP.1red=glmer(bin.cont~age.group*Direction+source + (Direction + source|sub),
              control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
              family="binomial", dati)
anova(IP.1,IP.1red) # Justied to drop additional terms. 
Anova(IP.1red,type=3)  

#Contrasts --> Source
cons=as.data.frame(emmeans(IP.1red, pairwise ~ source,adjust="none")$contrasts)
colnames(cons)[6]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#Contrasts --> Age
cons=as.data.frame(emmeans(IP.1red, pairwise ~ age.group,adjust="none")$contrasts)
colnames(cons)[6]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#Contrasts --> Age x Direction
cons=as.data.frame(emmeans(IP.1red, pairwise ~ age.group|Direction,adjust="none")$contrasts)
colnames(cons)[7]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

cons=as.data.frame(emmeans(IP.1red, pairwise ~ Direction|age.group,adjust="none")$contrasts)
colnames(cons)[7]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#RM#
#Suspicion
IP.susp=glmer(bin.cont~age.group*Direction+source + (Direction + source|sub),
              control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
              family="binomial", droplevels(dati[dati$evaluation =="Not Problematic",]))
Anova(IP.susp, type=3)

#Special Needs
IP.sn=glmer(bin.cont~age.group*Direction+source + (Direction + source|sub),
            control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
            family="binomial", droplevels(dati[dati$sn =="none",]))
Anova(IP.sn, type=3)

#No ceiling no floor
IP.ncnf=glmer(bin.cont~age.group*Direction+source + (Direction + source|sub),
              control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
              family="binomial", droplevels(dati[!dati$cat %in% c("floor", "ceil" ,"no.var"),]))
Anova(IP.ncnf, type=3)

#Defect
IP.defect=glmer(bin.cont~age.group*Direction+source + (Direction + source|sub),
             control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
             family="binomial", droplevels(dati[dati$eval=="ok",]))
Anova(IP.defect, type=3)

#Only Adolescents
IP.adol=glmer(bin.cont~age.group*Direction + source + (Direction|sub),
              control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
              family="binomial", droplevels(dati[dati$age.group %in% c("Young adolescents","Mid adolescents"),]))
Anova(IP.adol, type=3)

#Extreme Values
dati$IP.no.out=outliersZ.4(dati$bin.cont)
IP.n.no.out=glmer(IP.no.out~age.group*Direction+source + (Direction + source|sub),
                  control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
                  family="binomial", dati)
Anova(IP.n.no.out, type=3)

#CM#
#Gender
IP.gen=glmer(bin.cont~age.group*Direction+source +gender*age.group + (Direction + source|sub),
             control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
             family="binomial", dati)
Anova(IP.gen, type=3)

#GCR
IP.gcr=glmer(bin.cont~age.group*Direction+source + scale(guess.conversion.rate) + (Direction + source|sub),
             control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
             family="binomial", dati)
Anova(IP.gcr, type=3)

#ART
IP.art=glmer(bin.cont~age.group*Direction+source + scale(propCorr)+ (Direction + source|sub),
             control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
             family="binomial", dati)
Anova(IP.art,type=3) 

#Don1
IP.don1=glmer(bin.cont~age.group*Direction+source + scale(don1) + (Direction + source|sub),
              control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
              family="binomial", dati)
Anova(IP.don1,type=3)

cons=as.data.frame(emmeans(IP.don1, pairwise ~ age.group|Direction,adjust="none")$contrasts)
colnames(cons)[7]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#Delta
IP.delta=glmer(bin.cont~age.group*Direction+source +scale(delta.abs)+ (Direction+source|sub),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
               family="binomial", dati)
Anova(IP.delta,type=3)

#Conversion
IP.conv=glmer(bin.cont~age.group*Direction+source +scale(delta.conv)+ (Direction+source|sub),
                    control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
                    family="binomial", dati)
Anova(IP.conv,type=3)

#Block
IP.block=glmer(bin.cont~age.group*Direction+source +block.f+ (Direction+source|sub),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
               family="binomial", dati)
Anova(IP.block,type=3)

#Group Size
IP.gs=glmer(bin.cont~age.group*Direction+source + scale(group.size)+ (Direction+source|sub),
            control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
            family="binomial", dati)
Anova(IP.gs,type=3)

#Variance in Conforming Trials
IP.conf=glmer(bin.cont~age.group*Direction+source + c.ac.var + (Direction + source|sub),
              control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
              family="binomial", dati)
Anova(IP.conf,type=3)

# Computer Influence
IP.comp.bin=glmer(bin.cont~age.group*Direction+source + scale(comp.bin.cont)+ (Direction+source|sub),
                  control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
                  family="binomial", droplevels(dati[dati$source!="Computer",]))
Anova(IP.comp.bin,type=3)

##########################  INFLUENCE PROBABILITY CONTINUOUS  #################################

# Run models with different curve fitting techniques: 
mm.1=glmer(bin.cont~poly(age,1) + (1|sub), family="binomial",dati)
mm.1.exp=glmer(bin.cont~exp(poly(age,1)) + (1|sub),family="binomial", dati)
mm.1.log=glmer(bin.cont~log(age) + (1|sub),family="binomial", dati)
mm.1.inv=glmer(bin.cont~I(age^-1) + (1|sub),family="binomial", dati)
mm.2=glmer(bin.cont~I(age.2) + (1|sub),family="binomial", dati)
mm.3=glmer(bin.cont~I(age.3) + (1|sub),family="binomial", dati)
mm.1.2=glmer(bin.cont~poly(age,2) + (1|sub),family="binomial", dati)
mm.1.3=glmer(bin.cont~age.1+age.3 + (1|sub),family="binomial", dati)
mm.1.2.3=glmer(bin.cont~poly(age,3) + (1|sub),family="binomial", dati)

model.list=mget(grep("mm.", ls(),value=T))
aics<-lapply(model.list,function(x)AIC(x))
sort(unlist(aics))

#Basic Model
IPc.1=glmer(bin.cont~I(age^-1)*Direction*source+ (Direction + source|sub),
            control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
            family="binomial", dati)
Anova(IPc.1, type=3)

#Model Reduction
IPc.1red=glmer(bin.cont~I(age^-1)*Direction+source+ (Direction + source|sub),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               family="binomial", dati)
anova(IPc.1,IPc.1red)
Anova(IPc.1red,type=3)
summary(IPc.1red)

#Contrasts --> Direction x Age
cons=as.data.frame(emtrends(IPc.1red, pairwise ~ Direction, var='I(age^-1)',adjust="none")$contrasts)
colnames(cons)[6]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

emt=emtrends(IPc.1red, pairwise ~ Direction, var="I(age^-1)",adjust="Bonferroni")
summary(emt, infer =c(T,T))

#RM#
#Suspicion
IPc.susp=glmer(bin.cont~I(age^-1)*Direction+source+ (Direction + source|sub),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               family="binomial", droplevels(dati[dati$evaluation =="Not Problematic",]))
Anova(IPc.susp, type=3)

#Special Needs
IPc.sn=glmer(bin.cont~I(age^-1)*Direction+source+ (Direction + source|sub),
             control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
             family="binomial", droplevels(dati[dati$sn =="none",]))
Anova(IPc.sn, type=3)

#No ceiling no floor
IPc.ncnf=glmer(bin.cont~I(age^-1)*Direction+source+ (Direction + source|sub),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               family="binomial", droplevels(dati[!dati$cat %in% c("floor", "ceil" ,"no.var"),]))
Anova(IPc.ncnf,type=3)

#Defect
IPc.defect=glmer(bin.cont~I(age^-1)*Direction+source+ (Direction + source|sub),
                 control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
                 family="binomial", droplevels(dati[dati$eval=="ok",]))
Anova(IPc.defect,type=3)

#Only Adolescents
IPc.adol=glmer(bin.cont~I(age^-1)*Direction+source+ (Direction + source|sub),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               family="binomial", droplevels(dati[dati$age.group %in% c("Young adolescents","Mid adolescents"),]))
Anova(IPc.adol, type=3)

#Extreme Values
dati$IPc.no.out=outliersZ.4(dati$bin.cont)
IPc.n.no.out=glmer(IPc.no.out~I(age^-1)*Direction+source+ (Direction + source|sub),
                   control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
                   family="binomial", dati)
Anova(IPc.n.no.out, type=3)

#CM#
#Gender
IPc.gen=glmer(bin.cont~I(age^-1)*Direction+source+ gender*I(age^-1)+(Direction+source|sub),
              control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
              family="binomial",dati)
Anova(IPc.gen, type=3)

#GCR
IPc.gcr=glmer(bin.cont~I(age^-1)*Direction+source+ scale(guess.conversion.rate)+(Direction+source|sub),
              control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
              family="binomial",dati)
Anova(IPc.gcr, type=3)

#ART
IPc.art=glmer(bin.cont~I(age^-1)*Direction+source+ scale(propCorr)+(Direction+source|sub),
              control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
              family="binomial", dati)
Anova(IPc.art,type=3)

#First Donation
IPc.don1=glmer(bin.cont~I(age^-1)*Direction+source+ scale(don1)+(Direction+source|sub),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               family="binomial", dati)
Anova(IPc.don1,type=3)

#Delta
IPc.delta=glmer(bin.cont~I(age^-1)*Direction+source+ scale(delta.abs)+(Direction+source|sub),
                control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
                family="binomial", dati)
Anova(IPc.delta,type=3)

#Conversion
IPc.conv=glmer(bin.cont~I(age^-1)*Direction+source+ scale(delta.conv)+(Direction+source|sub),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               family="binomial", dati)
Anova(IPc.conv,type=3)

#Block
IPc.block=glmer(bin.cont~I(age^-1)*Direction+source+ block.f+(Direction+source|sub),
                control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
                family="binomial",dati)
Anova(IPc.block,type=3)

#Group Size
IPc.gs=glmer(bin.cont~I(age^-1)*Direction+source+ scale(group.size)+(Direction+source|sub),
             control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
             family="binomial",dati)
Anova(IPc.gs,type=3)

#Variance in Conforming Trials
IPc.conf=glmer(bin.cont~I(age^-1)*Direction+source+ c.ac.var + (Direction + source|sub),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               family="binomial", dati)
Anova(IPc.conf,type=3)

# Computer Influence
IPc.comp=glmer(bin.cont~I(age^-1)*Direction+source+ scale(comp.bin.cont)+(Direction+source|sub),
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
               family="binomial",droplevels(dati[dati$source !="Computer",]))
Anova(IPc.comp,type=3)
summary(IPc.comp)

##########################     REACTION TIME CATEGORICAL     ##################################

# Replacing RTs lower than 250 ms with NA
dati$rt.don2[dati$rt.don2<250]=NA

RT=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction+source)^3 + (Direction +source + bin.cont.f|sub),
        control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RT,type=3)

#Model Reduction
RT.1=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + age.group*Direction*source+ (Direction + source + bin.cont.f|sub),
          control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
anova(RT,RT.1) # Justified to drop additional terms.

RT.2=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ (Direction + source + bin.cont.f|sub),
          control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati) 
anova(RT.1,RT.2) # Justified to drop additional terms. 
Anova(RT.2,type=3)

#Contrasts --> Influenced or not
cons=as.data.frame(emmeans(RT.2, pairwise ~ bin.cont.f,adjust="none")$contrasts)
colnames(cons)[6]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#Contrasts --> Age x Influenced or not x Direction
cons=as.data.frame(emmeans(RT.2, pairwise ~ age.group|bin.cont.f|Direction,adjust="none")$contrasts)
colnames(cons)[8]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#Contrasts --> Opportunistic Conformity (Within Age)
cons=as.data.frame(emmeans(RT.2, pairwise ~ Direction|bin.cont.f|age.group,adjust="none")$contrasts)
colnames(cons)[8]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons$p.bonf[cons$p.bonf>1]=1
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#RM#
#Suspicion
RT.susp=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ (Direction + source + bin.cont.f|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),droplevels(dati[dati$evaluation =="Not Problematic",]))
Anova(RT.susp, type=3)

#Special Needs
RT.sn=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ (Direction + source + bin.cont.f|sub),
           control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
           droplevels(dati[dati$sn =="none",]))
Anova(RT.sn, type=3)

#No ceiling no floor
RT.ncnf=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ (Direction + source + bin.cont.f|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
             droplevels(dati[!dati$cat %in% c("floor", "ceil" ,"no.var"),]))
Anova(RT.ncnf, type=3)

#Defect
RT.defect=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ (Direction + bin.cont.f|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
            droplevels(dati[dati$eval=="ok",])) 
Anova(RT.defect, type=3)

#Only Adolescents
RT.adol=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ (Direction + source + bin.cont.f|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
             droplevels(dati[dati$age.group %in% c("Young adolescents","Mid adolescents"),]))
Anova(RT.adol, type=3)

#Extreme Values
dati$RT.no.out=outliersZ.4(dati$rt.don2)
RT.n.no.out=lmer(RT.no.out~(age.group+bin.cont.f+Direction)^3 + source+ (Direction + source + bin.cont.f|sub),
                 control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RT.n.no.out, type=3)

#CM#
#Gender
RT.gen=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ gender*age.group+(Direction + bin.cont.f|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RT.gen, type=3)

#GCR
RT.gcr=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+  scale(guess.conversion.rate) + (Direction + source + bin.cont.f|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RT.gcr,type=3)

#ART
RT.art=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+scale(propCorr) + (Direction + bin.cont.f|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RT.art,type=3) 

#Don1
RT.don1=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ scale(don1) + (Direction + source + bin.cont.f|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RT.don1,type=3)

#Delta 
RT.delta=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source + scale(delta.abs)  + (Direction + source + bin.cont.f|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RT.delta,type=3)

#Covnersion

RT.conv=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source + scale(delta.conv)  + (Direction + source + bin.cont.f|sub),
                   control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RT.conv,type=3)

#Block
RT.block=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+  block.f + (Direction + source + bin.cont.f|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RT.block,type=3)

#Group Size
RT.gs=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ scale(group.size) + (Direction + source + bin.cont.f|sub),
           control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RT.gs,type=3)

#Variance in Conforming Trials
RT.conf=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ c.ac.var + (Direction+bin.cont.f|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati) 
Anova(RT.conf,type=3)

# RT in Donation 1
RT.rt1=lmer(log(rt.don2)~(age.group+bin.cont.f+Direction)^3 + source+ log(rt.don1) + (Direction + source + bin.cont.f|sub),
           control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RT.rt1,type=3)

#############################    REACTION TIME CONTINUOUS    ##################################

# Run models with different curve fitting techniques: 
mm.1=lmer(log(rt.don2)~poly(age,1)*bin.cont.f + (bin.cont.f|sub), dati)
mm.1.exp=lmer(log(rt.don2)~exp(poly(age,1))*bin.cont.f + (bin.cont.f|sub), dati)
mm.1.log=lmer(log(rt.don2)~log(age)*bin.cont.f + (bin.cont.f|sub), dati)
mm.1.inv=lmer(log(rt.don2)~I(age^-1)+bin.cont.f + (bin.cont.f|sub), dati)
mm.2=lmer(log(rt.don2)~age.2*bin.cont.f + (bin.cont.f|sub), dati)
mm.3=lmer(log(rt.don2)~age.3*bin.cont.f + (bin.cont.f|sub), dati)
mm.1.2=lmer(log(rt.don2)~poly(age,2)*bin.cont.f + (bin.cont.f|sub), dati)
mm.1.3=lmer(log(rt.don2)~age.1*bin.cont.f + age.3*bin.cont.f+ (bin.cont.f|sub), dati) 
mm.1.2.3=lmer(log(rt.don2)~poly(age,3)*bin.cont.f+ (bin.cont.f|sub), dati)

model.list=mget(grep("mm.", ls(),value=T))
aics<-lapply(model.list,function(x)AIC(x))
sort(unlist(aics))

#Basic Model
RTc.1=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction+source)^3+(age.2+bin.cont.f+Direction+source)^3 
           +(age.3+bin.cont.f+Direction+source)^3 + (Direction+source+bin.cont.f|sub),
           control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
           dati)
Anova(RTc.1,type=3)

#Model Reduction
RTc.1red=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 
              +(age.3+bin.cont.f+Direction)^3 +source + (Direction+source+bin.cont.f|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
anova(RTc.1,RTc.1red)

RTc.2red=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 
              +source +age.3 + (Direction+source+bin.cont.f|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
anova(RTc.1red,RTc.2red)

RTc.3red=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 
              +source + (Direction+source+bin.cont.f|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
anova(RTc.2red,RTc.3red)
Anova(RTc.3red,type=3)

#Contrasts --> Influenced or not x Direction x Age
cons1=data.frame(emtrends(RTc.3red, pairwise~ bin.cont.f|Direction, var = "age.1", type="response", adjust='none')$contrasts)
cons1$trend="linear"
cons2=data.frame(emtrends(RTc.3red, pairwise~ bin.cont.f|Direction, var = "age.2", type="response", adjust="none")$contrasts)
cons2$trend="quadratic"
cons=rbind(cons1,cons2)
cons$p.unc=round(cons$p.value,3)
cons$p.value=round(cons$p.value*nrow(cons),3)
cons=cons[,-which(colnames(cons) %in% c("z.ratio", "df"))]
cons=cons[,c(6,1:5,7)]
cons[,c(4,5)]=round(cons[,c(4,5)],2)

emt1=emtrends(RTc.3red, pairwise~ bin.cont.f|Direction, var = "age.1", type="response", adjust='Bonferroni')
summary(emt1, infer =c(T,T))

emt2=emtrends(RTc.3red, pairwise~ bin.cont.f|Direction, var = "age.2", type="response", adjust='Bonferroni')
summary(emt2, infer =c(T,T))

#RM#
#Suspicion
RTc.susp=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 
              +source + (Direction+source+bin.cont.f|sub),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[dati$evaluation =="Not Problematic",]))
Anova(RTc.susp, type=3)

#Special Needs
RTc.sn=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 
            +source + (Direction+source+bin.cont.f|sub),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[dati$sn =="none",]))
Anova(RTc.sn, type=3)

#No ceiling no floor
RTc.ncnf=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 
              +source + (Direction+source+bin.cont.f|sub),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[!dati$cat %in% c("floor", "ceil" ,"no.var"),]))
Anova(RTc.ncnf,type=3)

#Defect
RTc.defect=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 
             +source + (Direction+source+bin.cont.f|sub),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[dati$eval=="ok",]))
Anova(RTc.defect,type=3)

#Only Adolescents
RTc.adol=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 
              +source + (Direction+source+bin.cont.f|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), 
              droplevels(dati[dati$age.group %in% c("Young adolescents","Mid adolescents"),]))
Anova(RTc.adol, type=3)

#Extreme Values
dati$RTc.no.out=outliersZ.4(dati$rt.don2)
RTc.n.no.out=lmer(log(RTc.no.out)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 
                  +source + (Direction+source+bin.cont.f|sub),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(RTc.n.no.out, type=3)

#CM#
#Gender
RTc.gen=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 +source + gender*age.1 + gender*age.2 + (Direction+source+bin.cont.f|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(RTc.gen, type=3)

#GCR
RTc.gcr=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 +source + scale(guess.conversion.rate) + (Direction+source+bin.cont.f|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(RTc.gcr,type=3)

#ART
RTc.art=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 +source + scale(propCorr) + (Direction+source+bin.cont.f|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(RTc.art,type=3) 

#Don1
RTc.don1=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 +source + scale(don1)+ (Direction+source+bin.cont.f|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(RTc.don1,type=3)

#Delta
RTc.delta=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 +source + scale(delta.abs) + (Direction+source+bin.cont.f|sub),
               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(RTc.delta,type=3)

#Conversion
RTc.conv=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 +source + scale(delta.conv) + (Direction+source+bin.cont.f|sub),
               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(RTc.conv,type=3)

#Block
RTc.block=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 +source + block.f + (Direction+bin.cont.f|sub),
               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(RTc.block,type=3)

#Group Size
RTc.gs=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3 + (age.2+bin.cont.f+Direction)^3 +source + scale(group.size) + (Direction+source+bin.cont.f|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(RTc.gs,type=3)

#Variance in Conforming Trials
RTc.conf=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3+(age.2+bin.cont.f+Direction)^3 
              +source + c.ac.var +(Direction+source+bin.cont.f|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(RTc.conf,type=3)

# RT in Donation 1
RTc.rt1=lmer(log(rt.don2)~(age.1+bin.cont.f+Direction)^3 + (age.2+bin.cont.f+Direction)^3 +source + log(rt.don1) + (Direction+source+bin.cont.f|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(RTc.rt1,type=3)
emtrends(RTc.rt1, pairwise~ bin.cont.f|Direction, var = "age.1", type="response", adjust='none')

########################    INFLUENCE MAGNITUDE CATEGORICAL    ################################

IM=lmer(cont~(age.group+Direction+source+scale(delta.abs))^3+(Direction+source+scale(delta.abs)|sub),
        control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(IM,type=3)

#Model Reduction
IM.1=lmer(cont~age.group*Direction*source +scale(delta.abs)*Direction*age.group + source*scale(delta.abs) + (Direction+source+scale(delta.abs)|sub),
          control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
anova(IM,IM.1)
Anova(IM.1,type=3)

IM.1red=lmer(cont~ age.group*Direction*source +scale(delta.abs)*Direction*age.group + (Direction+source+scale(delta.abs)|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
anova(IM.1,IM.1red)
Anova(IM.1red,type=3)

#Contrasts --> Source
cons=as.data.frame(emmeans(IM.1red, pairwise ~ source,adjust="none")$contrasts)
colnames(cons)[6]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#Contrasts --> Age x Direction x Delta
cons=as.data.frame(emtrends(IM.1red, pairwise ~ age.group|Direction,var="delta.abs",adjust="none")$contrasts)
colnames(cons)[7]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#Contrasts --> Opportunistic Conformity (Within age)
cons=as.data.frame(emtrends(IM.1red, pairwise ~ Direction|age.group,var="delta.abs",adjust="none")$contrasts)
colnames(cons)[7]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#Contrasts --> Age x Direction x Source
cons=as.data.frame(emmeans(IM.1red, pairwise ~ age.group|Direction|source,adjust="none")$contrasts)
colnames(cons)[8]="p.unc"
cons$estimate=round(cons$estimate,2)
cons$SE=round(cons$SE,2)
cons$p.bonf=cons$p.unc*nrow(cons)
cons[,c("p.unc", "p.bonf")]=round(cons[,c("p.unc", "p.bonf")],3)
cons=cons[,-which(colnames(cons)%in% c("z.ratio", "df"))]

#RM#
#Suspicion
IM.susp=lmer(cont~ age.group*Direction*source + scale(delta.abs)*Direction*age.group + (Direction + scale(delta.abs)|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[dati$evaluation =="Not Problematic",]))

#Special Needs
IM.sn=lmer(cont~ age.group*Direction*source + scale(delta.abs)*Direction*age.group + (Direction+source+ scale(delta.abs)|sub),
           control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[dati$sn =="none",]))
Anova(IM.sn, type=3)

#No ceiling no floor
IM.ncnf=lmer(cont~ age.group*Direction*source + scale(delta.abs)*Direction*age.group + (Direction+source+ scale(delta.abs)|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[!dati$cat %in% c("floor", "ceil" ,"no.var"),]))
Anova(IM.ncnf, type=3)

#Defect
IM.defect=lmer(cont~ age.group*Direction*source + scale(delta.abs)*Direction*age.group + (Direction+source+ scale(delta.abs)|sub),
               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[dati$eval=="ok",]))
Anova(IM.defect, type=3)

#Only Adolescents
IM.adol=lmer(cont~ age.group*Direction*source + age.group*Direction*scale(delta.abs) + (Direction+scale(delta.abs)|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
             droplevels(dati[dati$age.group %in% c("Young adolescents","Mid adolescents"),]))
Anova(IM.adol, type=3)

#Extreme Values
dati$IM.no.out=outliersZ.4(dati$cont)
IM.n.no.out=lmer(IM.no.out~ age.group*Direction*source + scale(delta.abs)*Direction*age.group + (Direction+source+ scale(delta.abs)|sub), 
                 control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(IM.n.no.out, type=3)

#CM#
#Gender
IM.gen=lmer(cont~age.group*Direction*source + scale(delta.abs)*Direction*age.group + gender*age.group + (Direction+source+scale(delta.abs)|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IM.gen, type=3)

#GCR
IM.gcr=lmer(cont~age.group*Direction*source + scale(delta.abs)*Direction*age.group + scale(guess.conversion.rate) + (Direction + source +scale(delta.abs)|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IM.gcr, type=3)

#ART
IM.art=lmer(cont~age.group*Direction*source + scale(delta.abs)*Direction*age.group + scale(propCorr)+ (Direction + source +scale(delta.abs)|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IM.art, type=3)

#Don1
IM.don1=lmer(cont~age.group*Direction*source + scale(delta.abs)*Direction*age.group + scale(don1) + (Direction + source+scale(delta.abs)|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IM.don1, type=3)

#Conversion
dati$cont.conv=scale(dati$cont*dati$guess.conversion.rate)
IM.conv=lmer(cont.conv~age.group*Direction*source + scale(delta.conv)*Direction*age.group + (Direction+scale(delta.conv)|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IM.conv,type=3) 

#Block
IM.block=lmer(cont~age.group*Direction*source + scale(delta.abs)*Direction*age.group + block.f + (Direction + source+scale(delta.abs)|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IM.block,type=3)

#Group Size
IM.gs=lmer(cont~age.group*Direction*source + scale(delta.abs)*Direction*age.group + scale(group.size) + (Direction + source+scale(delta.abs)|sub),
           control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IM.gs,type=3)

#Variance in Conforming Trials
IM.conf=lmer(cont~ age.group*Direction*source + age.group*Direction*scale(delta.abs) + c.ac.var + (Direction+source+scale(delta.abs)|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IM.conf,type=3)

# Computer Influence
IM.comp.mag=lmer(cont~age.group*Direction*source + scale(delta.abs)*Direction*age.group + scale(comp.cont) + (Direction + source+scale(delta.abs)|sub),
                 control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),droplevels(dati[dati$source!= "Computer",]))
Anova(IM.comp.mag,type=3)

# Capped Magnitude
to.cap=which(dati$cont > 0 & (dati$cont > dati$delta.abs))
dati$cont.capped=dati$cont
dati$cont.capped[to.cap]=dati$delta.abs[to.cap]
IM.capped=lmer(cont.capped~age.group*Direction*source + scale(delta.abs)*Direction*age.group  + (Direction +scale(delta.abs)|sub),
               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IM.capped,type=3)

###########################   INFLUENCE MAGNITUDE CONTINUOUS   ################################

# Run models with different curve fitting techniques: 
mm.1=lmer(cont~poly(age,1)+ (1|sub),dati)
mm.1.exp=lmer(cont~exp(poly(age,1)) + (1|sub), dati)
mm.1.log=lmer(cont~log(age)+ (1|sub),dati) 
mm.1.inv=lmer(cont~I(age^-1) + (1|sub), dati)
mm.2=lmer(cont~age.2+ (1|sub),dati) 
mm.3=lmer(cont~age.3 + (1|sub), dati)
mm.1.2=lmer(cont~poly(age,2) + (1|sub), dati)
mm.1.3=lmer(cont~age.1+age.3+(1|sub), dati)
mm.1.2.3=lmer(cont~poly(age,3) + (1|sub), dati)

model.list=mget(grep("mm.", ls(),value=T))
aics<-lapply(model.list,function(x)AIC(x))
sort(unlist(aics))

#Basic Model
IMc.1=lmer(cont~(poly(age,3)+Direction+scale(delta.abs)+source)^3+ (Direction + source + scale(delta.abs)|sub),
           control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(IMc.1,type=3)

#Are these higher level interactions due to outliers?
dati$IMc.no.out=outliersZ.4(dati$cont)
IMc.1out=lmer(IMc.no.out~(poly(age,3)+Direction+scale(delta.abs)+source)^3+ (Direction + source + scale(delta.abs)|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(IMc.1out,type=3) #Most higher level are explained by outliers - Revise - Reduce model to lower level interactions.

#Revised Model
IMc.1rev=lmer(cont~poly(age,3)*Direction*scale(delta.abs)+source*scale(delta.abs)+ (Direction + source + scale(delta.abs)|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
anova(IMc.1,IMc.1rev)
Anova(IMc.1rev,type=3)

IMc.1sep=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
           +age.3*Direction*scale(delta.abs)+source*scale(delta.abs)+(Direction+source+scale(delta.abs)|sub),
           control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IMc.1sep,type=3)

#Model Reduction
IMc.1red=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
              +age.3+source*scale(delta.abs)+(Direction+source+scale(delta.abs)|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
anova(IMc.1sep,IMc.1red)

IMc.2red=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
              +source*scale(delta.abs)+(Direction+source+scale(delta.abs)|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
anova(IMc.1red,IMc.2red)
Anova(IMc.2red,type=3)
summary(IMc.2red)

#Contrasts --> Source
emmeans(IMc.2red, pairwise ~ source, adjust="Bonferroni")

#Contrasts --> Age x Direction x Delta
#Prosocial
IMc.prosocial=lmer(cont~age.1*scale(delta.abs)+age.2*scale(delta.abs)+source*scale(delta.abs)
                   +(source+scale(delta.abs)|sub),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
                   droplevels(dati[dati$Direction=="Prosocial",]))
Anova(IMc.prosocial, type=3)
summary(IMc.prosocial)

#Selfish
IMc.selfish=lmer(cont~ age.1*scale(delta.abs)+age.2*scale(delta.abs)+source*scale(delta.abs)
                      +(source+scale(delta.abs)|sub),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
                      droplevels(dati[dati$Direction=="Selfish",]))
Anova(IMc.selfish, type=3)
summary(IMc.selfish)

#RM#
#Suspicion
IMc.susp=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
              +source*scale(delta.abs)+(Direction+scale(delta.abs)|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[dati$evaluation =="Not Problematic",]))
Anova(IMc.susp, type=3)

#Special Needs
IMc.sn=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
            +source*scale(delta.abs)+(Direction+source+scale(delta.abs)|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[dati$sn =="none",]))
Anova(IMc.sn, type=3)

#No ceiling no floor
IMc.ncnf=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
              +source*scale(delta.abs)+(Direction+source+scale(delta.abs)|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[!dati$cat %in% c("floor", "ceil" ,"no.var"),]))
Anova(IMc.ncnf,type=3)

#Defect
IMc.defect=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
                +source*scale(delta.abs)+(Direction+source+scale(delta.abs)|sub),
                control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[dati$eval=="ok",]))
Anova(IMc.defect,type=3)

#Only Adolescents
IMc.adol=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
              +source*scale(delta.abs)+(Direction+scale(delta.abs)|sub),
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
              droplevels(dati[dati$age.group %in% c("Young adolescents","Mid adolescents"),]))
Anova(IMc.adol,type=3)

#Extreme Values
dati$IMc.no.out=outliersZ.4(dati$cont)
IMc.n.no.out=lmer(IMc.no.out~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
                  +source*scale(delta.abs)+(Direction+source+scale(delta.abs)|sub),
                  control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IMc.n.no.out, type=3)

#CM#
#Gender
IMc.gen=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
             +source*scale(delta.abs)+gender*age.1 + gender*age.2 + (Direction+source+scale(delta.abs)|sub) ,
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IMc.gen, type=3)

#GCR
IMc.gcr=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
             +source*scale(delta.abs)+ scale(guess.conversion.rate) + (Direction+source+scale(delta.abs)|sub),
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IMc.gcr, type=3)

#ART
IMc.art=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
             +source*scale(delta.abs)+ scale(propCorr) + (Direction+source+scale(delta.abs)|sub) ,
             control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IMc.art, type=3)

#Don1
IMc.don1=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
              +source*scale(delta.abs)+ scale(don1) + (Direction+source+scale(delta.abs)|sub) ,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IMc.don1, type=3)

#Conversion (NOTE: dati$cont.conv created in IM section)
IMc.conv=lmer(cont~age.1*Direction*scale(delta.conv)+age.2*Direction*scale(delta.conv)
              +source*scale(delta.conv)+ (Direction+scale(delta.conv)|sub) ,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IMc.conv,type=3)

#Block
IMc.block=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
               +source*scale(delta.abs)+ block.f + (Direction+source+scale(delta.abs)|sub) ,
               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IMc.block, type=3)

#Group Size
IMc.gs=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
            +source*scale(delta.abs)+ scale(group.size) + (Direction+source+scale(delta.abs)|sub) ,
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(IMc.gs, type=3)
summary(IMc.gs)

#Variance in Conforming Trials
IMc.conf=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
              +source*scale(delta.abs)+ c.ac.var + (Direction+source+scale(delta.abs)|sub) ,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),dati)
Anova(IMc.conf,type=3)

# Computer Influenced
IMc.comp.mag=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
                  +source*scale(delta.abs)+ scale(comp.cont) + (Direction+source+scale(delta.abs)|sub),
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), droplevels(dati[dati$source != "Computer",]))
Anova(IMc.comp.mag, type=3)

# Capped Magnitude
to.cap=which(dati$cont > 0 & (dati$cont > dati$delta.abs))
dati$cont.capped=dati$cont
dati$cont.capped[to.cap]=dati$delta.abs[to.cap]
IMc.capped=lmer(cont~age.1*Direction*scale(delta.abs)+age.2*Direction*scale(delta.abs)
                +source*scale(delta.abs) + (Direction+scale(delta.abs)|sub) ,
            control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)), dati)
Anova(IMc.capped, type=3)
