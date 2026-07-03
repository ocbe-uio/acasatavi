#######################
##    ACASA-TAVI     ##
##    C. Cunen       ##
#######################

##############################
# Make the randomisation analysis dataset (with shamrand)
# Input: raw 
# Output: adran
###############################

set.seed(42)

library(tidyverse)
library(lubridate)
library(forcats)
library(labelled)

source("R/external/functions.R")
raw <- read_rds("data/raw/raw.rds")


ran <- raw %>% pick("ran")


adran <- ran %>% 
  select(subjectid,  site=sitename, ran_date = randat, 
         ran_trt,ran_trtcd)
adran$ran_dateTime <- adran$ran_date
adran$ran_date <- date(adran$ran_date)
table(adran$ran_trt)

# #############################
# # Introduce pseudorandomisation
# # Remove when running final analysis
# ###############################
# nn <- nrow(adran)
# ran_order <- sample(1:nn,nn,replace=F)
# adran$ran_trt <- adran$ran_trt[ran_order]
# adran$ran_trtcd <- adran$ran_trtcd[ran_order]
# table(adran$ran_trt,adran$ran_trtcd)

# Remove until here
###########################################################


readr::write_rds(adran, "data/ad/adran.rds")


