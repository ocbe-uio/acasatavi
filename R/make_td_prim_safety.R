#######################
##    ACASA-TAVI     ##
##    C. Cunen       ##
#######################

# Make the tabulation datasets (without treatment group): Safety co-primary endpoint
# And secondary safety endpoints related to bleeding

# Checking and combining information from:
# Clinical outcome (adjudicator): the classification of VARC bleeding events
# Clinical outcome (KU): stroke, MI, bleeding, deaths (also TIA which is used in secondary endpoint)
# End-of-study: deaths
# AEs: also deaths (just checking)

# "Safety composite" is defined in the end


####################################

library(tidyverse)
library(lubridate)

source("R/external/functions.R")
raw <- read_rds("data/raw/raw.rds")

# randomisation  
ran <- pick(raw,"ran")
ran <- ran %>% select(subjectid,  site=sitename, ran_date = randat)
ran$ran_dateTime <- ran$ran_date
ran$ran_date <- date(ran$ran_date)

# Clinical outcome (adjudicator)
so <- pick(raw,"so")

# Clinical outcome
ku <- pick(raw,"ku")

# End-of-study
eos <- pick(raw,"eos")

# AEs
ae_td <- readr::read_rds("data/td/ae_td.rds")

# checks on ids
id_list <- ran %>% select(subjectid)

# check no duplicates
id_list <- as.character(id_list$subjectid)
sum(duplicated(id_list))

nn <- length(id_list)
# 360 patients were randomized

table(so$eventname) # OBS some at 12 months and some at EOS


##### VARC events data ##### 
so <- so %>% select(subjectid, sovarc,so_eventdate=eventdate)
tablec(so$sovarc)
which(is.na(so$sovarc))
length(unique(so$subjectid))
# there are 352 patients here
# three patients are registered twice with "no bleeding"

so$varcyn <- ifelse(so$sovarc=="No bleeding","no","yes")

sov <- so %>% group_by(subjectid) %>%
  summarise(varc=ifelse(sum(varcyn=="yes",na.rm=T)>0,"yes","no"),
            varc1=ifelse(sum(sovarc=="Type 1",na.rm=T)>0,"yes","no"),
            varc2=ifelse(sum(sovarc=="Type 2",na.rm=T)>0,"yes","no"),
            varc3=ifelse(sum(sovarc=="Type 3",na.rm=T)>0,"yes","no"),
            varc4=ifelse(sum(sovarc=="Type 4",na.rm=T)>0,"yes","no"))
tablec(sov$varc)

##### Clinical Outcome ##### 
kut <- ku %>% select(subjectid,ku_dth=ku1_0,ku_dtdate=ku1_0dat,
                     strokeyn=ku1_2,strokedate=ku1_2dat,miyn=ku1_3,midate=ku1_3dat,
                     bleedingyn=ku2_0,bldate=ku2_0dat,tiayn=ku1_1,tiadate=ku1_1dat,
                     endocyn=ku3_4,endocdate=ku3_4dat)
suppressWarnings(
  kutt <- kut %>% group_by(subjectid) %>% 
  summarise(ku_death=(sum(ku_dth=="Ja",na.rm=T)>0),
            dtdate=min(ku_dtdate,na.rm=T),
            stroke=(sum(strokeyn=="Ja",na.rm=T)>0),
            stroke_date=min(strokedate,na.rm=T),
            mi=(sum(miyn=="Ja",na.rm=T)>0),
            mi_date=min(midate,na.rm=T),
            bleeding=(sum(bleedingyn=="Ja",na.rm=T)>0),
            bl_date=min(bldate,na.rm=T),
            tia=(sum(tiayn=="Ja",na.rm=T)>0),
            tia_date=min(tiadate,na.rm=T),
            endoc=(sum(endocyn=="Ja",na.rm=T)>0), # Endocarditis - used for HVD stage 1 definition
            endoc_date=min(endocdate,na.rm=T)) )
kutt$dtdate[is.infinite(kutt$dtdate)] <- NA
kutt$stroke_date[is.infinite(kutt$stroke_date)] <- NA
kutt$mi_date[is.infinite(kutt$mi_date)] <- NA
kutt$bl_date[is.infinite(kutt$bl_date)] <- NA
kutt$tia_date[is.infinite(kutt$tia_date)] <- NA
kutt$endoc_date[is.infinite(kutt$endoc_date)] <- NA

##### Find other deaths in EOS and AE ##### 
eosdt <- eos %>% select(subjectid,yn=eosyn,reas=eosreas,eosdtdate=eosdtdat)
eosdt$eos_death <- ifelse(eosdt$reas=="Death" & eosdt$yn=="No",1,0)
# eos_death=0 means that the subject was alive at EOS (not necessarily at 12 months)

aesdt <- ae_td %>% select(subjectid,aesdth,aestdat,aeendat)
aesdt <- aesdt %>% filter(aesdth=="Yes")
aesdt$aedate1 <- aesdt$aeendat
aesdt$aedate1[is.na(aesdt$aedate1)] <- aesdt$aestdat[is.na(aesdt$aedate1)]
aesdt <- aesdt %>% group_by(subjectid) %>%
  summarise(aedt=(sum(aesdth=="Yes")>0),aedate=max(aedate1))

alldt <- left_join(ran,kutt,by="subjectid")
alldt <- left_join(alldt,eosdt,by="subjectid")
alldt <- left_join(alldt,aesdt,by="subjectid")
alldt <- alldt %>% select(subjectid,ku_death,dtdate,eos_death,eosdtdate,aedt,aedate) 
alldt %>% filter(ku_death==TRUE | eos_death==1 | aedt==TRUE)
# EOS and AE agree

alldt$death <- ifelse(alldt$ku_death==TRUE | alldt$eos_death==1 | alldt$aedt==TRUE,"yes",NA)
alldt$death[is.na(alldt$death)] <- alldt$eos_death[is.na(alldt$death)] # those that are alive by EOS (not necessarily at 12 months)
alldt$death[alldt$death==0] <- "no"
alldt$death_date <- alldt$eosdtdate
alldt$death_date[is.na(alldt$death_date)] <- alldt$aedate[is.na(alldt$death_date)]
alldt$death_date[is.na(alldt$death_date)] <- alldt$dtdate[is.na(alldt$death_date)]

alldt <- alldt %>% select(subjectid,death,death_date)

##### All bleeding (checking) #####
kuBl <- kutt %>% select(subjectid,bleeding,bl_date)
allbl <- left_join(ran,sov,by="subjectid")
allbl <- left_join(allbl,kutt,by="subjectid")
allbl %>% filter(varc=="yes") %>% select(subjectid,varc,bleeding,bl_date)
# All VARC bleeding events are registered as bleedings in KU (except one - a patient that died, from bleeding (probably))
allbl %>% filter(bleeding)
# But in KU there are many more bleeding events (are the additional ones deemed as non-procedure-related?)

sovd <- allbl %>% filter(!is.na(varc))
sovd <- sovd %>% select(subjectid, varc,bl_date)
sovd$bl_date[which(sovd$varc!="yes")] <- NA

# Join everything
kutt <- kutt %>% select(subjectid,stroke,stroke_date,mi,mi_date,tia,tia_date,endoc)
which(is.na(kutt$stroke))
which(is.na(kutt$mi))

kutt$stroke <-ifelse(kutt$stroke,"yes","no")
kutt$mi <- ifelse(kutt$mi,"yes","no")

cso <- left_join(ran,kutt,by="subjectid")
cso <- left_join(cso,alldt,by="subjectid")
cso <- left_join(cso,sovd,by="subjectid")

length(unique(cso$subjectid))

cso$bl_date <- gsub("nk-nk","02-15",cso$bl_date) #replace missing day-in-year by 15th of February (we know if was between jan and march 2023)
cso$bl_date <- gsub("nk","15",cso$bl_date) #replace missing day-in-month by 15
cso$bl_date <- as.Date(cso$bl_date)

##### Defining safety composite #####

# Safety=1 if there is at least one "yes" (NAs don't matter)
# Safety=0 if all are "no" (NA if not all of them are "no")
cso$safety <- ifelse(cso$stroke=="yes" | cso$mi=="yes" | cso$death=="yes" | cso$varc=="yes",1,
                  ifelse((cso$stroke=="no" & cso$mi=="no" & cso$death=="no" & cso$varc=="no"),0,NA))
tablec(cso$safety)
# 6 with missing NA/MI and/or VARC:
#cso %>% filter(is.na(safety)) %>% select(subjectid,stroke,mi,death,varc) %>% print(n=17)
#cso %>% filter(safety==0) %>% select(subjectid,stroke,mi,death,varc) %>% print(n=32)

##### Safety2 also includes TIA (for secondary endpoint) #####
cso$safety2 <- ifelse(cso$safety==1 | cso$tia,1,
                     ifelse((cso$safety==0 & cso$tia==FALSE),0,NA))
tablec(cso$safety2)

# Count safety events:
cso$n_safety <- apply(cbind((cso$stroke=="yes"),(cso$mi=="yes"),(cso$death=="yes"),(cso$varc=="yes")),1,sum,na.rm=T)
tablec(cso$n_safety)
if(max(cso$n_safety,na.rm=T)>1){
  print("there is at least one patient with more than one safety events")
}

# Record event type and event date
cso$event_date <- NA
cc <- cso %>% filter(safety==1) %>% mutate(event_date=pmin(pmin(
  pmin(stroke_date,death_date,na.rm=T),mi_date,na.rm=T),bl_date,na.rm=T))
cc <- cc %>% mutate(event_type = case_when(
  event_date == stroke_date ~ "stroke",
  event_date == mi_date     ~ "mi",
  event_date == death_date  ~ "death",
  event_date == bl_date  ~ "bleeding",
      TRUE    ~ NA_character_ ))

cso <- cso %>% mutate(event_date=pmin(pmin(
  pmin(stroke_date,death_date,na.rm=T),mi_date,na.rm=T),bl_date,na.rm=T))
cso <- cso %>% mutate(event_type = case_when(
  event_date == stroke_date ~ "stroke",
  event_date == mi_date     ~ "mi",
  event_date == death_date  ~ "death",
  event_date == bl_date  ~ "bleeding",
  TRUE    ~ NA_character_ ))

# Factorize:
cso <- cso %>% mutate(across(c(stroke,mi,death,varc,safety,event_type),as.factor))

# Save
varc <- allbl %>% select(subjectid,ran_date,varc1,varc2,varc3,varc4,bl_date)
varc$bl_date <- gsub("nk","15",varc$bl_date) #replace missing day-in-month by 15
varc$bl_date <- as.Date(varc$bl_date)

readr::write_rds(cso, "data/td/cso_td.rds")
readr::write_rds(varc, "data/td/varc_td.rds")

# Check date of safety event (did it happen later than 12 mon?)
time <- cc$event_date-cc$ran_date
summary(time)


