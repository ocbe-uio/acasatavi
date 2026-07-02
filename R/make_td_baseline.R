#######################
##    ACASA-TAVI     ##
##    C. Cunen       ##
#######################

# Make the tabulation datasets (without treatment group): baseline characteristics

# Demographics
# Lab results
# Physical examination
# TAVI procedure characteristics
# Medical history
# Echography parameters (preTAVI and Baseline (postTAVI))

# Computation of frailty index based on some variables from physical examination and lab


####################################

library(tidyverse)
library(lubridate)

source("R/external/functions.R")
raw <- read_rds("data/raw/raw.rds")

# randomisation date 
ran <- pick(raw,"ran")

# demographics (age, sex, site)
dm <- pick(raw,"dm")

# Physical examination
phys <- pick(raw,"pe")

# Lab results
lab <- pick(raw,"lb")

# TAVI
tp <- pick(raw,"tp")

# Medical history
mh <- pick(raw,"mh")

# Baseline medication? Not included yet - not sure what is important
cm <- pick(raw,"cm")

# CTT (for aorta diameter): these are all at screening
ctt <- pick(raw,"ctt")

# Echocardiography at baseline
eco <- pick(raw,"eco")

# checks on ids
id_list <- ran %>% select(subjectid)

# check no duplicates
id_list <- as.character(id_list$subjectid)
sum(duplicated(id_list))

nn <- length(id_list)
# 360 patients were randomized


############################

# only include subjects that are randomized
# (should make no difference)
dm <- dm %>% filter(subjectid %in% id_list)
phys <- phys %>% filter(subjectid %in% id_list)
lab <- lab %>% filter(subjectid %in% id_list)
tp <- tp %>% filter(subjectid %in% id_list)
mh <- mh %>% filter(subjectid %in% id_list)
cm <- cm %>% filter(subjectid %in% id_list)
ctt <- ctt %>% filter(subjectid %in% id_list)
eco <- eco %>% filter(subjectid %in% id_list)

##### Demographics ##### 
dm0 <- dm %>% select(subjectid,site=sitename, age=dmage, sex)

##### Lab results ##### 
lab0 <- lab %>% filter(eventname == 'screening') %>% 
  select(subjectid, hemoglobin=lbhhbres, platelets=lbhplres, leukocytes=lbhlcres,
         INR=lbinrres, sodium=lbsodres, potassium=lbpotres, ALT=lbaltres,
         cardiacTroponin=lbcctres, NTProBNP=lbbnpres, creatinine=lbcreres,
         GFR=lbgfrres, glucose=lbglures, glycaHaemoglobin=lbhghre, 
         triglycerides=lbtrires, totchol=lbtchres, HDLchol=lbhdlres, LDLchol=lbldlres,
         albumin=lbalbres, CRP=lbcrpres, bilirubin=lbbilres)

##### Physical examination ##### 
phys0 <- phys %>% filter(eventname == 'screening')%>% 
  select(subjectid,ageVis=agevis, weight=peweight, height=peheight, 
         bmi=pebmi,systolicBP=pesysbp,diastolicBP=pediabp,smoke=pesmoke,
         body_temp=petemp,pulse_rate=pepulse, resp_rate=peresp,
         mini_cog=pecog,five_sts=pefive,abdomen=peares,abd_find=peadesc,
         abd_sign=peaclsig,#jug_vein=pecres,jug_find=pecdesc,jug_sign=pecclsig,no findings
         resp=perres,resp_find=perdesc,resp_sign=perclsig,circulation=pesres)
         #other=peores,oth_find=peodesc,oth_sign=peoclsig) no findings
phys0 <- phys0 %>% mutate(abnormal_find=ifelse(phys0$abdomen=="Abnormal"|phys0$resp=="Abnormal"|
                                                 phys0$circulation=="Not well circulated",1,0))
phys0$abnormal_find <- as.factor(phys0$abnormal_find)
attributes(phys0$abnormal_find)$label <- "Subject had abnormal results for abdomen, respiration or circulation."
phys0 <- phys0 %>% select(-c(abdomen,abd_find,abd_sign,resp,resp_find,resp_sign,circulation))
phys0$five_sts[phys0$five_sts==0] <- NA

#####  TAVI procedure ##### 
tp0 <- tp %>% select(subjectid, valve_type=tpvalve, tavi_post_dila=tppostd,
                     tavi_supra_pos=tp1)

#####  Medical history ##### 
mhb <- mh %>% mutate(hypertension1 = if_else(mhcat == "Hypertension", 1L, 0L),
                     cad1=if_else(mhcat == "Coronary artery disease", 1L, 0L),
                     diabetes1=if_else(mhcat == "Diabetes mellitus", 1L, 0L),
                     ch_obs_pulm1=if_else(mhcat == "Chronic obstructive pulmonary disease", 1L, 0L),
                     prev_stroke1=if_else(mhcat == "Previous stroke", 1L, 0L),
                     prev_pacemaker=if_else(mhcat == "Permanent pacemaker",1,0))
mh0 <- mhb %>% group_by(subjectid) %>% summarise(hypertension=as.numeric(sum(hypertension1)>0),
                                                 cad=as.numeric(sum(cad1)>0),
                                                 diabetes=as.numeric(sum(diabetes1)>0),
                                                 ch_obs_pulm=as.numeric(sum(ch_obs_pulm1)>0),
                                                 prev_stroke=as.numeric(sum(prev_stroke1)>0),
                                                 prev_pacemaker=as.numeric(sum(prev_pacemaker)>0))
mh0$hypertension <- as.factor(mh0$hypertension)
mh0$cad <- as.factor(mh0$cad)
mh0$diabetes <- as.factor(mh0$diabetes)
mh0$ch_obs_pulm <- as.factor(mh0$ch_obs_pulm)
mh0$prev_stroke <- as.factor(mh0$prev_stroke)
mh0$prev_pacemaker <- as.factor(mh0$prev_pacemaker)
ctt0 <- ctt %>% select(subjectid, asc_aorta_diam=cttaad)

#####  Echocardiography ##### 
# Screening data is PRE-TAVI data - save separately
# Baseline data is at baseline 

# Fixing some wrong values (0s are not possible)
# Hardcoding is unfortunate, but due to time limitation it was not possible to fix these errors in Viedoc.
eco$eclvgls[eco$eclvgls=="."] <- NA # I've asked about this
eco$eclvgls[eco$eclvgls=="0"] <- NA # I've asked about this
eco$eclvgls <- as.numeric(gsub(",",".",eco$eclvgls))
eco$eclvgls[which(eco$eclvgls>0)] <- -eco$eclvgls[which(eco$eclvgls>0)]
eco$eclvef[which(eco$eclvef==0)] <- NA
eco$ecstvo[which(eco$ecstvo==0)] <- NA
eco$ecstvoin[which(eco$ecstvoin==0)] <- NA
eco$ecava[which(eco$ecava==0)] <- NA

eco_scr <- eco %>% filter(eventname=="screening") %>% 
  select(subjectid,ecodate=ececdat,avmg=ecavmg,lvot=eclvot,
         lflg=eclflg,ntm=ecntm,lvef=eclvef,lvgls=eclvgls,
         stvo=ecstvo,svi=ecstvoin,ava=ecava,avr=ecoavr,
         avrvp=ecavrvp,avrm=ecoavrm,dvi=eco4)
eco_bas <- eco %>% filter(eventname=="baseline")  %>% 
  select(subjectid,ecodate=ececdat,avmg=ecavmg,lvot=eclvot,
         lflg=eclflg,ntm=ecntm,lvef=eclvef,lvgls=eclvgls,
         stvo=ecstvo,svi=ecstvoin,ava=ecava,avr=ecoavr,
         avrvp=ecavrvp,avrm=ecoavrm,dvi=eco4)


# Joining
baseline_td <- left_join(dm0,phys0,by="subjectid")
baseline_td <- left_join(baseline_td,lab0,by="subjectid")
baseline_td <- left_join(baseline_td,tp0,by="subjectid")
baseline_td <- left_join(baseline_td,mh0,,by="subjectid")
baseline_td <- left_join(baseline_td,ctt0,,by="subjectid")
baseline_td <- left_join(baseline_td,eco_bas,by="subjectid")

##### Frailty index ##### 
baseline_td$mini_cog[baseline_td$mini_cog=="Not done"] <- NA
baseline_td <- baseline_td %>% mutate(hemoglobin_low=(ifelse(sex=="Male" & round(hemoglobin,1)<13,1,0)+
                                                         ifelse(sex=="Female" & round(hemoglobin,1)<12,1,0)))
baseline_td <- baseline_td %>% mutate(cog_imp=(ifelse(mini_cog=="2" | mini_cog=="1" | mini_cog=="0",1,0)))
baseline_td <- baseline_td %>% mutate(sts=ifelse(round(five_sts,1)>=15,1,0))
baseline_td <- baseline_td %>% mutate(alb=ifelse(round(albumin)<35,1,0))
baseline_td <- baseline_td %>% mutate(tmp_sum = 
                                        rowSums(across(c(sts, cog_imp, 
                                                         hemoglobin_low, alb)), na.rm = TRUE),
    any_na = if_any(c(sts, cog_imp, hemoglobin_low, alb), is.na),
    
    eft = case_when(
      tmp_sum >= 3               ~ 3L,            # cap at 3; ignore remaining NAs
      tmp_sum <  3 & any_na      ~ NA_integer_,   # below 3 AND some NA → NA
      TRUE                       ~ tmp_sum        # below 3 and no NAs → actual sum
    )
  ) %>%
  select(-tmp_sum, -any_na,-hemoglobin_low,-alb)       # optional: clean up
attributes(baseline_td$eft)$label <- "Essential Frailty Toolset"
baseline_td$sts <- as.factor(baseline_td$sts)
baseline_td$cog_imp <- as.factor(baseline_td$cog_imp)
baseline_td <- baseline_td %>% mutate(frailty_status = case_when(eft == 0 ~ "robust",
                                                                 eft > 0 & eft < 3 ~ "pre-frail",
                                                                 eft >= 3  ~ "frail",
                                                                 TRUE ~ "missing"),
                                      frailty_status = factor(frailty_status,
                    levels = c("robust", "pre-frail", "frail","missing")))

# Save
readr::write_rds(eco_scr, "data/td/preTAVI_eco_td.rds")
readr::write_rds(baseline_td, "data/td/baseline_td.rds")
