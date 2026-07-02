## Primary efficacy endpoint (blinded)

# Descriptive statistics
# Main analysis (unadjusted, with MI)
# Sensitivity analyses

# Reviewed and checked primary analysis 30JUN2026 Inge Christoffer Olsen

###############################

source("R/external/functions.R")

library(tidyverse)
#library(glmmTMB)
#library(DHARMa)
library(marginaleffects)
library(broom)


adsl <- read_rds("data/ad/adsl.rds") # with shamrand
baseline <- read_rds("data/td/baseline_td.rds")
eff <- read_rds("data/td/cto_td.rds")

adsl <- adsl %>% select(-site,-ran_date)
baseline <- baseline %>% select(-site)

eff <- eff %>% left_join(adsl,by="subjectid")
eff <- eff %>% mutate(halt2=as.factor(ifelse(halt=="yes" | death=="yes","yes","no")))

eff_base <- eff %>% left_join(baseline,by="subjectid")

eff_base <- eff_base %>% mutate(smoke2 = 
                                  fct_collapse(smoke,
                                               smoker = c("On a daily basis", "Smokes occasionally")))


# Subset for modelling
eff_base_mod <- eff_base %>% select(halt2,ran_trt,bmi,sodium,potassium,ALT,triglycerides,albumin,bilirubin) 
p <- ncol(eff_base_mod)-1

idc <- which(apply(eff_base_mod[, 3:(p+1)], 1, function(x) any(is.na(x))))
eff_base_mod[idc,]

# Do something about covariate missing values (replace by median or most common class)
for (i in 3:(p+1)){
  vv <- unlist(eff_base_mod[,i])
  if (sum(is.na(vv))==0) next
  if (!is.factor(vv)){
    eff_base_mod[is.na(vv),i] <- median(vv,na.rm=T)
  }else{
    tt <- table(vv)
    eff_base_mod[is.na(vv),i] <- names(sort(tt,T)[1])
  }
}
eff_base_mod[idc,]


# Timing of CT
time_to_ct <- eff$ct_date-eff$ran_date
summary(time_to_ct)
which(time_to_ct>400)

##### Descriptive statistics #####
n_notmissing <- sum(!is.na(eff$halt))
n_halt <- length(which(eff$halt=="yes"))
# 24% event rate

n_notmissing <- sum(!is.na(eff$halt2))
n_halt2 <- length(which(eff$halt2=="yes"))
# 28% event rate

effsum0 <- eff %>% group_by(ran_trt,halt) %>% summarise(n=n(),.groups = "drop_last") %>%
  pivot_wider(names_from=ran_trt,values_from=n)
n_halt_doac <- as.numeric(effsum0 %>% filter(halt=="yes")%>%select(DOAC))
n_halt_asa <- as.numeric(effsum0 %>% filter(halt=="yes")%>%select(ASA))
n_tot_doac <- sum(effsum0 %>% filter(!is.na(halt))%>%select(DOAC))
n_tot_asa <- sum(effsum0 %>% filter(!is.na(halt))%>%select(ASA))
pct_doac <- 100*n_halt_doac/n_tot_doac
pct_asa <- 100*n_halt_asa/n_tot_asa
pct_all <- 100*(n_halt_doac+n_halt_asa)/(n_tot_doac+n_tot_asa)

effsum <- eff %>% group_by(ran_trt,cthalt) %>% summarise(n=n(),.groups = "drop_last") %>%
  pivot_wider(names_from=ran_trt,values_from=n)
effsum <- effsum %>% mutate(Overall=rowSums(effsum[,2:3]))
effsum$cthalt <- as.character(effsum$cthalt)
effsum$cthalt[is.na(effsum$cthalt)] <- "Missing HALT"
effsum <- effsum %>% mutate(DOAC=as.character(DOAC),ASA=as.character(ASA),
                            Overall=as.character(Overall))
#effsum <- effsum %>% 
#  add_row(cthalt = "Total HALT", DOAC =paste0(n_halt_doac," (",round(pct_doac,1),"%)"),
#          ASA =paste0(n_halt_asa," (",round(pct_asa,1),"%)"))
effsum



effsum02 <- eff %>% group_by(ran_trt,halt2) %>% summarise(n=n(),.groups = "drop_last") %>%
  pivot_wider(names_from=ran_trt,values_from=n)
n_halt2_doac <- as.numeric(effsum02 %>% filter(halt2=="yes")%>%select(DOAC))
n_halt2_asa <- as.numeric(effsum02 %>% filter(halt2=="yes")%>%select(ASA))
n_tot2_doac <- sum(effsum02 %>% filter(!is.na(halt2))%>%select(DOAC))
n_tot2_asa <- sum(effsum02 %>% filter(!is.na(halt2))%>%select(ASA))
pct2_doac <- 100*n_halt2_doac/n_tot2_doac
pct2_asa <- 100*n_halt2_asa/n_tot2_asa
pct2_all <- 100*(n_halt2_doac+n_halt2_asa)/(n_tot2_doac+n_tot2_asa)

effsum_d <- eff %>% group_by(ran_trt,death) %>% summarise(n=n(),.groups = "drop_last") %>%
  pivot_wider(names_from=ran_trt,values_from=n)
n_death_doac <- as.numeric(effsum_d %>% filter(death=="yes")%>%select(DOAC))
n_death_asa <- as.numeric(effsum_d %>% filter(death=="yes")%>%select(ASA))

effsum2 <- effsum %>% 
  add_row(cthalt = "Death", DOAC =paste0(n_death_doac),ASA =paste0(n_death_asa),
          Overall=paste0(n_death_doac+n_death_asa)) %>% 
  add_row(cthalt = "HALT", DOAC =paste0(n_halt_doac," (",round(pct_doac,1),"%)"),
          ASA =paste0(n_halt_asa," (",round(pct_asa,1),"%)"),
          Overall=paste0(n_halt_doac+n_halt_asa," (",round(pct_all,1),"%)") ) %>% 
  add_row(cthalt = "HALT & Death", DOAC =paste0(n_halt2_doac," (",round(pct2_doac,1),"%)"),
          ASA =paste0(n_halt2_asa," (",round(pct2_asa,1),"%)"),
          Overall=paste0(n_halt2_doac+n_halt2_asa," (",round(pct2_all,1),"%)")) 
effsum2 <- effsum2[c(6,1,8:9,2:5,7),]




##### Main analysis ##### 
library(mice)

dat_mice0 <- mice(eff_base_mod,maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable

dat_mice <- mice(eff_base_mod,m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice <- with(dat_mice,glm(halt2~ran_trt,family="binomial"))
halt_preds <- avg_predictions(mod_mice,by="ran_trt")
halt_main_ratio <- avg_comparisons(mod_mice,variables=list(ran_trt = c("ASA", "DOAC")),
                            comparison="lnratioavg",transform=exp)

halt_main_diff <- avg_comparisons(mod_mice,variables=list(ran_trt = c("ASA", "DOAC")))

# Diagnostics for imputation
plot(dat_mice) # should not have any trends
diagMI <- data.frame(dataset=rep(NA,21),pDOAC=NA,pASA=NA)
diagMI$dataset <- c("complete",paste(1:20))

for (i in 1:21){
  datai <- complete(dat_mice,action=(i-1)) 
  datai <- datai[!is.na(datai$halt2),]
  diagMI$pDOAC[i] <- sum(datai$ran_trt=="DOAC"&datai$halt2=="yes")/sum(datai$ran_trt=="DOAC" )
  diagMI$pASA[i] <- sum(datai$ran_trt=="ASA"&datai$halt2=="yes")/sum(datai$ran_trt=="ASA")
}
diagMI_long <- pivot_longer(diagMI,cols=2:3)
diagMI_long$value <- as.numeric(diagMI_long$value)
ggplot(diagMI_long,aes(x=dataset,y=value,color=name,shape=name))+geom_point()
# Imputation model sometimes imputes more events in the DOAC group 
# but not very much more than in the complete data

# Summary figure
haltplot <- ggplot(halt_preds, aes(x = ran_trt, y = 100*estimate, fill = ran_trt)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  geom_errorbar(
    aes(ymin = 100*conf.low, ymax = 100*conf.high),
    position = position_dodge(width = 0.8),
    linewidth = 0.1,width=0.1
  )+
  geom_vline(xintercept=0.3,linetype=1)+
  geom_hline(yintercept=0,linetype=1)+
  ylim(0,40) +
  labs(x="Error bars: 95% CI",
         y="Patients with HALT (%)",
         title="HALT")+
  geom_text(
    aes(label = paste0(round(100*estimate),"%"),y=100*conf.high+2),
    size  = 4
  ) +
  theme_bw() +
  theme(legend.position = "none", # no legend
        panel.border    = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust=1, face = "bold"))  +
  scale_fill_manual(values = c(
    DOAC = "#0571b0",
    ASA = "#d95f02"
  ))

######## Sensitivity analyses ########################

### Adjusted analysis with MI
halt_mod_mice_sens1 <- with(dat_mice,glm(halt2~ran_trt+bmi+sodium+potassium+
                                           ALT+triglycerides+albumin+bilirubin,family="binomial"))
halt_sens1_ratio <- avg_comparisons(halt_mod_mice_sens1,variables=list(ran_trt = c("ASA", "DOAC")),
                                   comparison="lnratioavg",transform=exp)

halt_sens1_diff <- avg_comparisons(halt_mod_mice_sens1,variables=list(ran_trt = c("ASA", "DOAC")))

### Unadjusted analysis - complete case
halt_mod_sens2 <- glm(halt2~ran_trt,family="binomial",data=eff_base_mod)
halt_sens2_ratio <- avg_comparisons(halt_mod_sens2,variables=list(ran_trt = c("ASA", "DOAC")),
                                    comparison="lnratioavg",transform=exp)

halt_sens2_diff <- avg_comparisons(halt_mod_sens2,variables=list(ran_trt = c("ASA", "DOAC")))

### Unadjusted analysis - worst case imputation
idm <- which(is.na(eff_base_mod$halt2))
eff_base_mod2 <- eff_base_mod
eff_base_mod2$halt2[idm] <- "yes"
halt_mod_sens3 <- glm(halt2~ran_trt,family="binomial",data=eff_base_mod2)
halt_sens3_ratio <- avg_comparisons(halt_mod_sens3,variables=list(ran_trt = c("ASA", "DOAC")),
                                    comparison="lnratioavg",transform=exp)

halt_sens3_diff <- avg_comparisons(halt_mod_sens3,variables=list(ran_trt = c("ASA", "DOAC")))

### Saving stuff
save(effsum2,haltplot, halt_main_ratio,halt_main_diff,halt_sens1_diff,
     halt_sens1_ratio,halt_sens2_diff,halt_sens2_ratio,halt_sens3_diff,
     halt_sens3_ratio,file="data/res/prim_eff_tab.RData")


##########################################################################

# Checking a model
# m00 <- glmmTMB(halt2~ran_trt,family=binomial(),data=eff_base)
# 
# m0 <- glmmTMB(halt2~ran_trt+sex*age+smoke2+bilirubin,family=binomial(),data=eff_base)
# res <- simulateResiduals(m0)
# plot(res)
# 
# summary(m0)
# 
# # Brief look at baseline covariates
# eff_base %>% 
#   filter(!is.na(halt)) %>% 
#   group_by(sex, halt) %>% 
#   summarise(n = n()) %>% 
#   group_by(sex) %>% 
#   mutate(tot = sum(n),
#          pct = round(n/tot*100,digits = 1))
# # women seem to get more HALT (also when including death)
# 
# eff_base %>% 
#   filter(!is.na(halt)) %>% 
#   group_by(smoke, halt) %>% 
#   summarise(n = n()) %>% 
#   group_by(smoke) %>% 
#   mutate(tot = sum(n),
#          pct = round(n/tot*100,digits = 1))
# # smoking may have some effect (previous smokers seem a bit protected from HALT?
# # less so when death is included)
# 
# eff_base %>% 
#   filter(!is.na(halt)) %>% 
#   group_by(valve_type, halt) %>% 
#   summarise(n = n()) %>% 
#   group_by(valve_type) %>% 
#   mutate(tot = sum(n),
#          pct = round(n/tot*100,digits = 1))
# # valve-type has no effect (maybe when including death)
# 
# eff_base %>% 
#   filter(!is.na(halt)) %>% 
#   group_by(tavi_post_dila, halt) %>% 
#   summarise(n = n()) %>% 
#   group_by(tavi_post_dila) %>% 
#   mutate(tot = sum(n),
#          pct = round(n/tot*100,digits = 1))
# # TAVI post dilatation has no effect
# 
# eff_base %>% 
#   filter(!is.na(halt)) %>% 
#   group_by(tavi_supra_pos, halt) %>% 
#   summarise(n = n()) %>% 
#   group_by(tavi_supra_pos) %>% 
#   mutate(tot = sum(n),
#          pct = round(n/tot*100,digits = 1))
# # TAVI supra pos has no effect (maybe when including death)
# 
# eff_base %>% 
#   filter(!is.na(halt)) %>% 
#   group_by(hypertension, halt) %>% 
#   summarise(n = n()) %>% 
#   group_by(hypertension) %>% 
#   mutate(tot = sum(n),
#          pct = round(n/tot*100,digits = 1))
# # hypertension has no effect (maybe when including death)
# 
# eff_base %>% 
#   filter(!is.na(halt)) %>% 
#   group_by(cad, halt) %>% 
#   summarise(n = n()) %>% 
#   group_by(cad) %>% 
#   mutate(tot = sum(n),
#          pct = round(n/tot*100,digits = 1))
# # cad has no effect (maybe when including death)
# 
# 
# eff_base %>% 
#   filter(!is.na(halt)) %>% 
#   group_by(diabetes, halt) %>% 
#   summarise(n = n()) %>% 
#   group_by(diabetes) %>% 
#   mutate(tot = sum(n),
#          pct = round(n/tot*100,digits = 1))
# # diabetes has no efefct
# 
# eff_base %>% 
#   filter(!is.na(halt)) %>% 
#   group_by(ch_obs_pulm, halt) %>% 
#   summarise(n = n()) %>% 
#   group_by(ch_obs_pulm) %>% 
#   mutate(tot = sum(n),
#          pct = round(n/tot*100,digits = 1))
# # pulmonary obstructive has no effect (maybe when inclduing death)
# 
# eff_base %>% 
#   filter(!is.na(halt)) %>% 
#   group_by(prev_stroke, halt) %>% 
#   summarise(n = n()) %>% 
#   group_by(prev_stroke) %>% 
#   mutate(tot = sum(n),
#          pct = round(n/tot*100,digits = 1))
# # previous stroke may have some effect (also when including death)
# 
# 
# ggplot(eff_base, aes(x = age, y = as.numeric(halt == "yes"))) +
#   geom_point(alpha = 0.2, position = position_jitter(height = 0.02)) +
#   geom_smooth(method = "glm",
#               method.args = list(family = "binomial"),
#               se = TRUE) +
#   scale_y_continuous(labels = scales::percent_format(),
#                      limits = c(0, 1)) +
#   labs(y = "P(halt == 'yes')") #+ facet_grid(~sex)
# # no general relationship with age (but maybe within women)
# 
# ggplot(eff_base, aes(x = bmi, y = as.numeric(halt == "yes"))) +
#   geom_point(alpha = 0.2, position = position_jitter(height = 0.02)) +
#   geom_smooth(method = "glm",
#               method.args = list(family = "binomial"),
#               se = TRUE) +
#   scale_y_continuous(labels = scales::percent_format(),
#                      limits = c(0, 1)) +
#   labs(y = "P(halt == 'yes')")
# higher BMI seems protective (but there are few patients at the extremes)

# little relationship with: BP, GFR, platelets, INR, cardiactroponingT, NTProBNP, totchol, LDLchol
# some relationship with: sodium, maybe potassium, hemoglobin, maybe asc_aorta_diam, maybe leukocytes
# ALT, creatinine, glucose, glycaHaemoglobin, triglycerides, HDLchol, albumin, maybe CRP, bilirubin

# OBS there's probably an error in on INR measurment - multiplied by 100?


# ggplot(all_comp1,aes(predicted_hi,predicted_lo)) + 
#   geom_point() +
#   geom_abline(slope=1,intercept=0,linetype=3) +
#   coord_fixed() + 
#   labs(x="ASA",y="DOAC")
# 
# hist1 <- ggplot(all_comp1,aes(estimate)) +
#   geom_histogram(bins=100,fill="grey") + theme_classic() +
#   geom_vline(xintercept=mean(all_comp1$estimate),color="orange",linewidth=1) +
#   geom_vline(xintercept=median(all_comp1$estimate),color="darkgreen",linewidth=1)+
#   labs(x="ASA - DOAC")
