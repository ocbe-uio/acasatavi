#######################
##    ACASA-TAVI     ##
##    C. Cunen       ##
#######################

# Make the tabulation datasets (without treatment group): eligibility and adherence

# Checking and combining information from:
# End-of-study (discontinuations)
# Protocol deviations
# Clinical outcome (switching between drugs)
# Medication (trial drug, which and also dosage)

# "Treatment deviation" is defined in the end


####################################


library(tidyverse)
library(lubridate)

source("R/external/functions.R")
raw <- read_rds("data/raw/raw.rds")

# inclusion exclusion (screening dates etc)
ie <- pick(raw,"ie")

# Medication (trial drug)
md <- pick(raw,"md")

# Prior and concomitant medication:
cm <- pick(raw,"cm")

# randomisation 
ran <- pick(raw,"ran")

# end of study
#eos <- pick(raw,"eos")
eos <- read_rds("data/td/dates_end_td.rds")

# protocol deviations
dv <- pick(raw,"dv")

# Safety questions # not used here
#sq <- pick(raw,"sq")

# Clinical outcome
ku <- pick(raw,"ku")

# Comments
co <- pick(raw,"co") # Maybe check? Some explanations for discontinuation

# checks on ids
id_list <- ran %>% select(subjectid)

# check no duplicates
id_list <- as.character(id_list$subjectid)
sum(duplicated(id_list))

nn <- length(id_list)
# 360 patients were randomized

length(which(ran$subjectid %in% ie$subjectid))
length(which(ie$subjectid %in% ran$subjectid))
# all screened patients were deemed elligble and randomized

# only include subjects that are randomized
# (should make no difference)
dv <- dv %>% filter(subjectid %in% id_list)
eos <- eos %>% filter(subjectid %in% id_list)
ie <- ie %>% filter(subjectid %in% id_list)
md <- md %>% filter(subjectid %in% id_list)
cm <- cm %>% filter(subjectid %in% id_list)
ku <- ku %>% filter(subjectid %in% id_list)

##### End-of-study ##### 
tablec(eos$eos_yn) # did the subject complete the study according to protocol?
tablec(eos$eos_reas)
table(eos$eos_yn,eos$eos_reas)

nonCompl <- eos %>% filter(eos_yn=="No")
length(unique(nonCompl$subjectid))

eos <- eos %>% mutate(non_completers=ifelse(eos_yn=="No",1,0))
attributes(eos$non_completers)$label <- "Subject did not complete according to protocol (including death)"
eosa <- eos %>% select(subjectid,ran_date,eos_date,eos_yn,non_completers,eos_reas) 

all <- ran %>% 
  select(subjectid,  site=sitename)
all <- left_join(all,eosa,by="subjectid")
all <- all %>% mutate(time_in_study=eos_date-ran_date)

#####  Protocol deviations ##### 
length(unique(dv$subjectid))
tablec(dv$dvcat)
tablec(dv$dvplclas)
# all protocol deviations are deemed "not important"

idc <- which(nonCompl$subjectid %in% dv$subjectid)
nonCompl[idc,]

idcc <- which(dv$subjectid %in% nonCompl$subjectid)
dv[idcc,]

dv <- dv %>% mutate(pd_not_important1=ifelse(dvplclas=="Not important" | is.na(dvplclas),1,0), #OBS adding NAs here
                    pd_important1=ifelse(dvplclas=="Important",1,0))
dvp <- dv %>% group_by(subjectid) %>% summarise(pd_not_important=(sum(pd_not_important1)>0),
                                                pd_important=(sum(pd_important1,na.rm=T)>0))
all <- left_join(all,dvp,by="subjectid")
attributes(all$pd_important)$label <- "Subject experienced important Protocol Deviation"
attributes(all$pd_not_important)$label <- "Subject experienced not important Protocol Deviation"

#####  Clinical outcome - did patient stop taking study drug? (including switches) ##### 
kut <- ku %>% select(subjectid,changeTRT=ku4_2,change_date=ku4_2dat,change_reas=ku4_2_1,
                     change_dose=ku4_2_2,dose_reas=ku4_2_2b)
suppressWarnings(
  kutt <- kut %>% group_by(subjectid) %>%
  summarize(changeTRT_ku=(sum(changeTRT=="Ja",na.rm=T)>0), # or just sum(changeTRT=="Ja",na.rm=T)
            firstchange_date=min(change_date,na.rm=T),
            change_dose_ku=(sum(change_dose=="Nei",na.rm=T)>0))  )
kutt$firstchange_date[is.infinite(kutt$firstchange_date)] <- NA


all <- left_join(all,kutt,by="subjectid")
attributes(all$changeTRT_ku)$label <- "Did patient stop or switch study drug?" 

##### Medication (also adherence to treatment) #####  
mds <- md %>% select(subjectid,eventname,eventdate,mdmed,mdprev,mdchan)
md_wide <- mds %>% pivot_wider(names_from=eventname,values_from=c(mdmed,mdprev,mdchan,eventdate))
md_wide <- md_wide %>% select(-c("mdmed_baseline","mdmed_12 months","mdprev_screening","mdprev_3 months",
                                 "mdprev_6 months","mdprev_9 months","mdprev_12 months","mdchan_12 months",
                                 "mdprev_baseline","mdchan_screening","mdchan_baseline","eventdate_baseline"))
colnames(md_wide) <- gsub(" months","mo",colnames(md_wide))
md_wide <- md_wide %>% mutate(changeTRT=((mdchan_3mo=="Yes") + (mdchan_6mo=="Yes") + (mdchan_9mo=="Yes")))
attributes(md_wide$changeTRT)$label <- "Number of times the subject changed treatment groups?"
md_wide <- md_wide %>% mutate(change1_date=as.Date(ifelse(mdchan_3mo=="Yes",eventdate_3mo,
                                                   ifelse(mdchan_6mo=="Yes",eventdate_6mo,
                                                          ifelse(mdchan_9mo=="Yes",eventdate_9mo,NA)))))
attributes(md_wide$change1_date)$label <- "Date of first change of treatment group"

md_wide <- md_wide %>% mutate(time_trt=(rowSums(
  across(mdmed_3mo:mdmed_9mo, ~ .x == mdmed_screening),na.rm=T)+1)*3)
attributes(md_wide$time_trt)$label <- "Number of months on rand. treatment" # for those that drop out the number is a "maximum"

md_wide <- md_wide %>% select(-c("eventdate_screening","eventdate_3mo","eventdate_6mo",
                                 "eventdate_9mo","eventdate_12mo"))
  
# Look at dosing and frequency of ASA
table(md$mdasau)
table(md$mdasado)
table(md$mdasafr)
# no-one changes dose or frequency of ASA


# Look at dosing and frequency of DOAC and changes of DOAC type
table(md$mddoac)

table(md$mdapido) #change
table(md$mdapiu) #no change
table(md$mdapifr) # one change

table(md$mdedodo) #change
table(md$mdedou) # no change
table(md$mdedofr) #no change

table(md$mdrivdo) #change
table(md$mdrivu) #no change
table(md$mdrivfr) #no change

mdd <- md %>% select(subjectid,eventname,mddoac,mdapido,mdapifr,mdedodo,mdrivdo)
mdd <- mdd %>% mutate(dose=coalesce(mdapido,mdedodo,mdrivdo)) %>% select(-c(mdapido,mdedodo,mdrivdo))
mdd_wide <- mdd %>% pivot_wider(names_from=eventname,values_from=c(mddoac,dose,mdapifr))
mdd_wide <- mdd_wide %>% select(-c("mddoac_baseline","dose_baseline","mdapifr_baseline",
                                   "mddoac_12 months","dose_12 months","mdapifr_12 months"))
count_f <- function(vec){
  vec <- vec[!is.na(vec)]
  if (length(vec)==0){
    ll <- NA
  }else{
    ll <- length(unique(vec))
  }
  return(ll>1)
}
mdd_wide <- mdd_wide %>% rowwise()%>% 
  mutate(change_doac=count_f(c_across("mddoac_screening":"mddoac_9 months")))
attributes(mdd_wide$change_doac)$label <- "Change of DOAC type?"  

mdd_wide <- mdd_wide %>% rowwise()%>% 
  mutate(change_dose=ifelse(!change_doac,
                            count_f(c_across("dose_screening":"dose_9 months")) ,
                            FALSE),
         change_freq=ifelse(!change_doac,
                            count_f(c_across("mdapifr_screening":"mdapifr_9 months")) ,
                            FALSE))
attributes(mdd_wide$change_dose)$label <- "Change of dose?" 
attributes(mdd_wide$change_freq)$label <- "Change of drug frequency?"  

# Join
md_wide <- left_join(md_wide,mdd_wide[,c("subjectid","change_doac","change_dose","change_freq")],
                     by="subjectid")
md_wideB <- md_wide %>% select(subjectid,changeTRT:change_freq)
allUB <- left_join(all,md_wide,by="subjectid")
allB <- left_join(all,md_wideB,by="subjectid")

##### Define participants with deviations from randomized treatment (switching, stopping and early discontinuation) ##### 
allB$changeTRT[which(is.na(allB$changeTRT) & !is.na(allB$eos_date))] <- 0
#allB$non_completer_notDead <- ifelse(allB$eos_reas=="Voluntary discontinuation by patient" &
#                                      allB$eos_reas!="Death",1,0)
allB$non_completer_notDead <- ifelse(allB$eos_reas!="Death",1,0)
allB$changeTRT_all <- ifelse(allB$changeTRT>0 | allB$changeTRT_ku==T,1,0)
allB$changeTRT_all[which(is.na(allB$changeTRT_all) & allB$eos_reas=="Death")] <- 0
allB <- allB %>% mutate(switch_discont=ifelse(non_completer_notDead==1 | changeTRT_all==1,1,0))
allB <- allB %>% mutate(reason=case_when(changeTRT_all==1 ~ "switching/stopping", 
                                       (changeTRT_all==0 | is.na(changeTRT_all)) & eos_reas=="Voluntary discontinuation by patient" ~ "voluntary discontinuation",
                                       (changeTRT_all==0 | is.na(changeTRT_all)) & eos_reas=="Incorrect randomisation" ~ "incorrect randomisation") )
tablec(allB$switch_discont)
tablec(allB$reason)
#allB %>% filter(is.na(switch_discont))
# many NAs, but we can assume that they have not switched or discontinued
allB$switch_discont[is.na(allB$switch_discont)] <- 0

# Fix missing change dates
#allB %>% filter(switch_discont==1) %>% select(firstchange_date,change1_date) %>% print(n=55)
allB <- allB %>% mutate(change_date=firstchange_date)
allB$change_date[which(is.na(allB$firstchange_date) & !is.na(allB$change1_date))] <- allB$change1_date[which(is.na(allB$firstchange_date) & !is.na(allB$change1_date))]


readr::write_rds(allUB, "data/td/adherence_unblind_td.rds")
readr::write_rds(allB, "data/td/adherence_blind_td.rds")


#################################################
# Checking adherence OLD
# tablec(md$eventname)
# 
# md_scr <- md %>% filter(eventname=="screening") %>%
#   select(subjectid,  site=sitecode, date=eventdate,trt=mdmed,change=mdchan)
# table(md_scr$trt)
# table(md_scr$change) # Will the subject change treatment group?
# 
# # md_bl 
# # nothing entered in the baseline table
# 
# md_3mo <- md %>% filter(eventname=="3 months") %>%
#   select(subjectid,  site=sitecode, date=eventdate,trt=mdmed,prev_trt=mdprev,
#          change=mdchan) 
# sum(is.na(md_3mo$date)) # no NAs but 12 patients are not present here
# tablec(md_3mo$prev_trt)
# tablec(md_3mo$change) # will the subject change treatment group? 24
# table(md_3mo$change,md_3mo$trt) # some switch in both groups
# table(md_3mo$change,md_3mo$prev_trt)
# tablec(md_3mo$trt)
# 
# md_6mo <- md %>% filter(eventname=="6 months") %>%
#   select(subjectid,  site=sitecode, date=eventdate,trt=mdmed,prev_trt=mdprev,
#          change=mdchan) 
# sum(is.na(md_6mo$date)) # no NAs but 13 patients are not present here
# tablec(md_6mo$prev_trt)
# table(md_6mo$change) # will the subject change treatment group? 9
# table(md_6mo$change,md_6mo$trt) # some switch in both groups
# tablec(md_6mo$trt)
# 
# md_9mo <- md %>% filter(eventname=="9 months") %>%
#   select(subjectid,  site=sitecode, date=eventdate,trt=mdmed,prev_trt=mdprev,
#          change=mdchan) 
# sum(is.na(md_9mo$date)) # no NAs but 16 patients are not present here
# table(md_9mo$prev_trt)
# table(md_9mo$change) # will the subject change treatment group? 5
# table(md_9mo$change,md_9mo$trt) # some switch in both groups
# table(md_9mo$trt)
# 
# md_12mo <- md %>% filter(eventname=="12 months") %>%
#   select(subjectid,  site=sitecode, date=eventdate,trt=mdmed,prev_trt=mdprev,
#          change=mdchan) 
# sum(is.na(md_12mo$date)) 
# table(md_12mo$prev_trt)
# table(md_12mo$change) # very many change here - but it's the end of study so...
# table(md_12mo$change,md_12mo$trt) # some switch in both groups
# table(md_12mo$trt)


