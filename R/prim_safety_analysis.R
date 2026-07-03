## Primary safety endpoint (blinded)

# Descriptive statistics
# Main analysis (adjusted with MI)
# Sensitivity analyses

#Reviewed 30JUN2026 Inge Christoffer Olsen
# Comment: The treatment should be included as an explanatory variable in the 
# imputation model because it is a non-inferiority assessment. Yes, this has been added now (CC).

###############################

source("R/external/functions.R")

library(tidyverse)
#library(glmmTMB)
#library(DHARMa)
library(marginaleffects)
library(broom) # for nice model tables


adsl <- read_rds("data/ad/adsl.rds") # with shamrand
baseline <- read_rds("data/td/baseline_td.rds")
#extra <- read_rds("data/td/baseline_extra_td.rds")
saf <- read_rds("data/td/cso_td.rds")

adsl <- adsl %>% select(-c(site,ran_date))
baseline <- baseline %>% select(-site)
#extra <- extra %>% select(subjectid,euroscore=score)

saf <- saf %>% left_join(adsl,by="subjectid")

saf_base <- saf %>% left_join(baseline,by="subjectid")
#saf_base <- saf_base %>% left_join(extra,by="subjectid")

saf_base <- saf_base %>% mutate(smoke2 = 
                                  fct_collapse(smoke,
                                               smoker = c("On a daily basis", "Smokes occasionally")))

# Subset for modelling
saf_base_mod <- saf_base %>% select(safety,ran_trt,age,sex,cad,diabetes,hypertension,prev_stroke,GFR) #,frailty_status,bmi,euroscore
p <- ncol(saf_base_mod)-1

idc <- which(apply(saf_base_mod[, 3:(p+1)], 1, function(x) any(is.na(x))))
saf_base_mod[idc,]

# Do something about covariate missing values (replace by median or most common class)
for (i in 3:(p+1)){
  vv <- unlist(saf_base_mod[,i])
  if (sum(is.na(vv))==0) next
  if (!is.factor(vv)){
    saf_base_mod[is.na(vv),i] <- median(vv,na.rm=T)
  }else{
    tt <- table(vv)
    saf_base_mod[is.na(vv),i] <- names(sort(tt,T)[1])
  }
}
saf_base_mod[idc,]


# Timing of safety event
time_to_event <- saf$event_date-saf$ran_date
summary(time_to_event)

# Checking all missing safety indicator
saf %>% filter(is.na(safety)) %>% select(subjectid,stroke,mi,death,varc,pp)
# the missing are early discontinuators

##### Descriptive statistics #####
n_notmissing <- sum(!is.na(saf$safety))
n_saf <- length(which(saf$safety==1))
# 9% event rate

safsum0 <- saf %>% group_by(ran_trt,safety,.drop=FALSE) %>% summarise(n=n(),.groups = "drop_last") %>%
  pivot_wider(names_from=ran_trt,values_from=n)
n_saf_doac <- as.numeric(safsum0 %>% filter(safety==1)%>%select(DOAC))
n_saf_asa <- as.numeric(safsum0 %>% filter(safety==1)%>%select(ASA))
n_tot_doac <- sum(safsum0 %>% filter(!is.na(safety))%>%select(DOAC))
n_tot_asa <- sum(safsum0 %>% filter(!is.na(safety))%>%select(ASA))
pct_doac <- 100*n_saf_doac/n_tot_doac
pct_asa <- 100*n_saf_asa/n_tot_asa
pct_all <- 100*(n_saf_doac+n_saf_asa)/(n_tot_doac+n_tot_asa)

safsum <- saf %>% group_by(ran_trt,event_type,.drop=FALSE) %>% summarise(n=n(),.groups = "drop_last") %>%
  pivot_wider(names_from=ran_trt,values_from=n)
safsum <- safsum %>% filter(!is.na(event_type))
colnames(safsum0)[1] <- "event_type"
safsum <- bind_rows(safsum0[which(safsum0$event_type==0),],safsum,
                    safsum0[which(is.na(safsum0$event_type)),])
safsum <- safsum %>% mutate(event_type = 
                              fct_recode(event_type,`No safety events`="0"))
safsum <- safsum %>% mutate(event_type = fct_na_value_to_level(event_type, level = "Missing"))
safsum$ASA[is.na(safsum$ASA)] <- 0
safsum <- safsum %>% mutate(Overall=rowSums(safsum[,2:3]))
safsum <- safsum %>% mutate(DOAC=as.character(DOAC),ASA=as.character(ASA),
                            Overall=as.character(Overall))
safsum <- safsum %>% 
  add_row(event_type = "Safety events", DOAC =paste0(n_saf_doac," (",round(pct_doac,1),"%)"),
          ASA =paste0(n_saf_asa," (",round(pct_asa,1),"%)"),
          Overall=paste0(n_saf_doac+n_saf_asa," (",round(pct_all,1),"%)"))
safsum <- safsum[c(6,1,7,2:5),]


##### Main analysis #####
# Adjusted analysis - multiple imputation with imputation model - containing TREATMENT
library(mice)

dat_mice <- mice(saf_base_mod,m=20,printFlag = F,seed=11254)
saf_mod_mice <- with(dat_mice,glm(safety~ran_trt+age+sex+cad+prev_stroke+
                                    diabetes+hypertension+GFR,family="binomial"))
safety_main_diff <- avg_comparisons(saf_mod_mice,variables=list(ran_trt = c("ASA", "DOAC")),
                            equivalence=c(NA,0.119))


safety_main_ratio <- avg_comparisons(saf_mod_mice,variables=list(ran_trt = c("ASA", "DOAC")),
                             comparison="lnratioavg",transform=exp)


# Diagnostics for imputation
plot(dat_mice) # should not have any trends
diagMI <- data.frame(dataset=rep(NA,21),pDOAC=NA,pASA=NA)
diagMI$dataset <- c("complete",paste(1:20))

for (i in 1:21){
  datai <- complete(dat_mice,action=(i-1)) 
  datai <- datai[!is.na(datai$safety),]
  diagMI$pDOAC[i] <- sum(datai$ran_trt=="DOAC"&datai$safety=="1")/sum(datai$ran_trt=="DOAC" )
  diagMI$pASA[i] <- sum(datai$ran_trt=="ASA"&datai$safety=="1")/sum(datai$ran_trt=="ASA")
}
diagMI_long <- pivot_longer(diagMI,cols=2:3)
diagMI_long$value <- as.numeric(diagMI_long$value)
ggplot(diagMI_long,aes(x=dataset,y=value,color=name,shape=name))+geom_point()
# Imputation model imputes a bit more events in the DOAC group 
# but not very much more than in the complete data


# Summary figure
safetyplot <- ggplot(safety_main_diff, aes(y=1, x = 100*estimate)) +
  geom_point(colour="#008837",size=3) +
  geom_errorbar(
    aes(xmin = 100*conf.low, xmax = 100*conf.high),
    linewidth = 0.3,width=0.1,colour="#008837")+
  ylim(0,1.7) +xlim(-10,15)+
  geom_vline(xintercept=11.9,linetype=2)+
  geom_vline(xintercept=0,linetype=1)+
  geom_hline(yintercept=0,linetype=1)+
  labs(y=" ",
       x="Absolute risk difference, DOAC minus ASA (percentage points)",
       title="Safety composite")+
  geom_text(
    aes(label = paste0("95% CI ",round(100*conf.low,1)," to ",
                       round(100*conf.high,1), " pp"),y=1),
    size  = 6,vjust=2,colour="#008837"
  ) +
  geom_text(
    aes(label = paste0(round(100*estimate,1)," pp"),y=1),
    size  = 6,vjust=-1.5,colour="#008837"
  ) +
  geom_text(
    aes(label = "No difference (0)",y=1.7,x=0),
    size  = 4,vjust=-2,
  ) +
  geom_text(
    aes(label = "NI margin (+11.9 pp)",y=1.7,x=11.9),
    size  = 4,vjust=-2,
  ) +
  coord_cartesian(clip = "off") + 
  theme_bw() +
  theme(legend.position = "none", # no legend
        panel.border    = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust=5, face = "bold"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y   = element_blank(),
        plot.margin        = margin(20, 5.5, 5.5, 5.5),
        text = element_text(size = 15))  



######## Sensitivity analyses ########################

### Unadjusted analysis - multiple imputation
saf_mod_mice_sens1 <- with(dat_mice,glm(safety~ran_trt,family="binomial"))
safety_sens1_diff <- avg_comparisons(saf_mod_mice_sens1,variables=list(ran_trt = c("ASA", "DOAC")),
                                    equivalence=c(NA,0.119))


safety_sens1_ratio <- avg_comparisons(saf_mod_mice_sens1,variables=list(ran_trt = c("ASA", "DOAC")),
                                     comparison="lnratioavg",transform=exp)


### Adjusted analysis - complete case
saf_mod_sens2 <- glm(safety~ran_trt+age+sex+cad+prev_stroke+
                       diabetes+hypertension+GFR,family="binomial",data=saf_base_mod)
safety_sens2_diff <- avg_comparisons(saf_mod_sens2,variables=list(ran_trt = c("ASA", "DOAC")),
                                     equivalence=c(NA,0.119))


safety_sens2_ratio <- avg_comparisons(saf_mod_sens2,variables=list(ran_trt = c("ASA", "DOAC")),
                                      comparison="lnratioavg",transform=exp)


### Adjusted analysis - worst-case imputation
idm <- which(is.na(saf_base_mod$safety))
saf_base_mod2 <- saf_base_mod
saf_base_mod2$safety[idm] <- 1
saf_mod_sens3 <- glm(safety~ran_trt+age+sex+cad+prev_stroke+
                       diabetes+hypertension+GFR,family="binomial",data=saf_base_mod2)
safety_sens3_diff <- avg_comparisons(saf_mod_sens3,variables=list(ran_trt = c("ASA", "DOAC")),
                                     equivalence=c(NA,0.119))


safety_sens3_ratio <- avg_comparisons(saf_mod_sens3,variables=list(ran_trt = c("ASA", "DOAC")),
                                      comparison="lnratioavg",transform=exp)



### Saving stuff

save(safsum,safetyplot,safety_main_diff,safety_main_ratio,safety_sens1_diff,
     safety_sens1_ratio,safety_sens2_diff,safety_sens2_ratio,safety_sens3_diff,
     safety_sens3_ratio,file="data/res/prim_saf_tab.RData")
