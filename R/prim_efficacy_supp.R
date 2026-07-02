## Primary efficacy endpoint: supplmentary analysis (hypothetical for death)


###############################

source("R/external/functions.R")

library(tidyverse)
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
eff_base_mod <- eff_base %>% select(halt,ran_trt,bmi,sodium,potassium,ALT,triglycerides,albumin,bilirubin) 
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


sum(is.na(eff_base_mod$halt))

### Unadjusted analysis - multiple imputation with imputation model - not containing TREATMENT
library(mice)

dat_mice0 <- mice(eff_base_mod,maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable

dat_mice <- mice(eff_base_mod,m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice <- with(dat_mice,glm(halt~ran_trt,family="binomial"))
halt_supp_ratio <- avg_comparisons(mod_mice,variables=list(ran_trt = c("ASA", "DOAC")),
                                   comparison="lnratioavg",transform=exp)

halt_supp_diff <- avg_comparisons(mod_mice,variables=list(ran_trt = c("ASA", "DOAC")))

# Diagnostics for imputation
plot(dat_mice) # should not have any trends
diagMI <- data.frame(dataset=rep(NA,21),pDOAC=NA,pASA=NA)
diagMI$dataset <- c("complete",paste(1:20))

for (i in 1:21){
  datai <- complete(dat_mice,action=(i-1)) 
  datai <- datai[!is.na(datai$halt),]
  diagMI$pDOAC[i] <- sum(datai$ran_trt=="DOAC"&datai$halt=="yes")/sum(datai$ran_trt=="DOAC" )
  diagMI$pASA[i] <- sum(datai$ran_trt=="ASA"&datai$halt=="yes")/sum(datai$ran_trt=="ASA")
}
diagMI_long <- pivot_longer(diagMI,cols=2:3)
diagMI_long$value <- as.numeric(diagMI_long$value)
ggplot(diagMI_long,aes(x=dataset,y=value,color=name,shape=name))+geom_point()
# Imputation model sometimes imputes more events in the DOAC group 
# but not very much more than in the complete data

### Saving stuff
save(halt_supp_ratio,halt_supp_diff,file="data/res/prim_eff_supp_tab.RData")
