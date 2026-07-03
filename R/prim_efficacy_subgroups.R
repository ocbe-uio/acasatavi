#######################
##    ACASA-TAVI     ##
##    C. Cunen       ##
#######################

## Primary efficacy endpoint - subgroup analysis

# Descriptive statistics
# Unadjusted, with imputation

###############################

source("R/external/functions.R")

library(tidyverse)
library(marginaleffects)
#library(broom)
#library(forcats)
library(ggplot2)
#library(rlang)
#library(purrr)


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

# Define subgroup factors
eff_base <- eff_base %>%
  mutate(ageF = factor(if_else(age < 75, "<75", ">= 75"),
                       levels = c("<75", ">= 75")))
eff_base <- eff_base %>%
  mutate(renal_function = factor(if_else(GFR < 30, "impaired", "preserved"),
                                 levels = c("impaired", "preserved")))
eff_base <- eff_base %>%
  mutate(ascF = factor(if_else(asc_aorta_diam < 30, "<30 mm", ">= 30 mm"),
                       levels = c("<30 mm", ">= 30 mm")))

eff_base$diabetesF <- factor(ifelse(eff_base$diabetes=="1","yes","no"))
eff_base$hypertensionF <- factor(ifelse(eff_base$hypertension=="1","yes","no"))
eff_base$pdF <- factor(ifelse(eff_base$tavi_post_dila=="Yes","yes","no"))

# Subset for modelling
eff_base_mod <- eff_base %>% select(halt2,ran_trt,bmi,sodium,potassium,ALT,triglycerides,albumin,bilirubin,
                                    ageF,sex,diabetesF,hypertensionF,renal_function,valve_type,
                                    pdF,ascF,frailty_status) 
p <- 8

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

### Unadjusted analysis - multiple imputation with imputation model - not containing TREATMENT
library(mice)

dat_mice0 <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                  "ALT","triglycerides","albumin","bilirubin")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable

dat_mice <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                 "ALT","triglycerides","albumin","bilirubin")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice <- with(dat_mice,glm(halt2~ran_trt,family="binomial"))
halt_main_ratio <- avg_comparisons(mod_mice,variables=list(ran_trt = c("ASA", "DOAC")),
                                   comparison="lnratioavg",transform=exp)


# age
dat_mice0 <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                  "ALT","triglycerides","albumin","bilirubin","ageF")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable
pred[, "ageF"] <- 0 # to ensure that imputation model does not use subgroup 

dat_mice <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                 "ALT","triglycerides","albumin","bilirubin","ageF")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice1 <- with(dat_mice,glm(halt2~ran_trt*ageF,family="binomial"))
halt_sub1 <- avg_comparisons(mod_mice1,variables=list(ran_trt = c("ASA", "DOAC")),by=c("ageF"),
                         comparison="lnratioavg",transform=exp)


# sex
dat_mice0 <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                  "ALT","triglycerides","albumin","bilirubin","sex")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable
pred[, "sex"] <- 0 # to ensure that imputation model does not use subgroup

dat_mice <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                 "ALT","triglycerides","albumin","bilirubin","sex")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice2 <- with(dat_mice,glm(halt2~ran_trt*sex,family="binomial"))
halt_sub2 <- avg_comparisons(mod_mice2,variables=list(ran_trt = c("ASA", "DOAC")),by=c("sex"),
                             comparison="lnratioavg",transform=exp)

# diabetes
idm <- which(is.na(eff_base_mod$diabetesF))
dat_mice0 <- mice(eff_base_mod[-idm,c("halt2","ran_trt","bmi","sodium","potassium",
                                  "ALT","triglycerides","albumin","bilirubin","diabetesF")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable
pred[, "diabetesF"] <- 0 # to ensure that imputation model does not use subgroup 

dat_mice <- mice(eff_base_mod[-idm,c("halt2","ran_trt","bmi","sodium","potassium",
                                 "ALT","triglycerides","albumin","bilirubin","diabetesF")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice3 <- with(dat_mice,glm(halt2~ran_trt*diabetesF,family="binomial"))
halt_sub3 <- avg_comparisons(mod_mice3,variables=list(ran_trt = c("ASA", "DOAC")),by=c("diabetesF"),
                             comparison="lnratioavg",transform=exp)
#plot(dat_mice)

# hypertension
idm <- which(is.na(eff_base_mod$hypertensionF))
dat_mice0 <- mice(eff_base_mod[-idm,c("halt2","ran_trt","bmi","sodium","potassium",
                                  "ALT","triglycerides","albumin","bilirubin","hypertensionF")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable
pred[, "hypertensionF"] <- 0 # to ensure that imputation model does not use subgroup 

dat_mice <- mice(eff_base_mod[-idm,c("halt2","ran_trt","bmi","sodium","potassium",
                                 "ALT","triglycerides","albumin","bilirubin","hypertensionF")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice4 <- with(dat_mice,glm(halt2~ran_trt*hypertensionF,family="binomial"))
halt_sub4 <- avg_comparisons(mod_mice4,variables=list(ran_trt = c("ASA", "DOAC")),by=c("hypertensionF"),
                             comparison="lnratioavg",transform=exp)

# renal function
#idm <- which(is.na(eff_base_mod$renal_function))
dat_mice0 <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                  "ALT","triglycerides","albumin","bilirubin","renal_function")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable
pred[, "renal_function"] <- 0 # to ensure that imputation model does not subgroup 

dat_mice <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                 "ALT","triglycerides","albumin","bilirubin","renal_function")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice5 <- with(dat_mice,glm(halt2~ran_trt*renal_function,family="binomial"))
halt_sub5 <- avg_comparisons(mod_mice5,variables=list(ran_trt = c("ASA", "DOAC")),by=c("renal_function"),
                             comparison="lnratioavg",transform=exp)

# valve_type
idm <- which(is.na(eff_base_mod$valve_type))
dat_mice0 <- mice(eff_base_mod[-idm,c("halt2","ran_trt","bmi","sodium","potassium",
                                  "ALT","triglycerides","albumin","bilirubin","valve_type")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable
pred[, "valve_type"] <- 0 # to ensure that imputation model does not subgroup 

dat_mice <- mice(eff_base_mod[-idm,c("halt2","ran_trt","bmi","sodium","potassium",
                                 "ALT","triglycerides","albumin","bilirubin","valve_type")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice6 <- with(dat_mice,glm(halt2~ran_trt*valve_type,family="binomial"))
halt_sub6 <- avg_comparisons(mod_mice6,variables=list(ran_trt = c("ASA", "DOAC")),by=c("valve_type"),
                             comparison="lnratioavg",transform=exp)

# post-dilatation
idm <- which(is.na(eff_base_mod$pdF))
dat_mice0 <- mice(eff_base_mod[-idm,c("halt2","ran_trt","bmi","sodium","potassium",
                                  "ALT","triglycerides","albumin","bilirubin","pdF")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable
pred[, "pdF"] <- 0 # to ensure that imputation model does not subgroup 

dat_mice <- mice(eff_base_mod[-idm,c("halt2","ran_trt","bmi","sodium","potassium",
                                 "ALT","triglycerides","albumin","bilirubin","pdF")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice7 <- with(dat_mice,glm(halt2~ran_trt*pdF,family="binomial"))
halt_sub7 <- avg_comparisons(mod_mice7,variables=list(ran_trt = c("ASA", "DOAC")),by=c("pdF"),
                             comparison="lnratioavg",transform=exp)

# ascending aorta diameter
#idm <- which(is.na(eff_base_mod$ascF))
dat_mice0 <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                  "ALT","triglycerides","albumin","bilirubin","ascF")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable
pred[, "ascF"] <- 0 # to ensure that imputation model does not subgroup 

dat_mice <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                 "ALT","triglycerides","albumin","bilirubin","ascF")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice8 <- with(dat_mice,glm(halt2~ran_trt*ascF,family="binomial"))
halt_sub8 <- avg_comparisons(mod_mice8,variables=list(ran_trt = c("ASA", "DOAC")),by=c("ascF"),
                             comparison="lnratioavg",transform=exp)

# frailty
#idm <- which(is.na(eff_base_mod$frailty_status))
dat_mice0 <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                  "ALT","triglycerides","albumin","bilirubin","frailty_status")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ran_trt"] <- 0 # to ensure that imputation model does not use treatment variable
pred[, "frailty_status"] <- 0 # to ensure that imputation model does not subgroup 

dat_mice <- mice(eff_base_mod[,c("halt2","ran_trt","bmi","sodium","potassium",
                                 "ALT","triglycerides","albumin","bilirubin","frailty_status")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice9 <- with(dat_mice,glm(halt2~ran_trt*frailty_status,family="binomial"))
halt_sub9 <- avg_comparisons(mod_mice9,variables=list(ran_trt = c("ASA", "DOAC")),by=c("frailty_status"),
                             comparison="lnratioavg",transform=exp)

# Nice figure
summarise_subgroup <- function(data, subgroup_var) {
  subgroup_var <- enquo(subgroup_var)
  
  # Get variable name as a string
  var_name <- quo_name(subgroup_var)
  
  # Remove trailing "F" (if present) and capitalise
  subgroup_label <- var_name %>%
    str_replace_all("_"," ") %>%
    str_remove("F$") %>%   # remove "F" at the end of the string
    str_to_sentence()         # capitalise each word; use str_to_sentence() if you prefer
  
  
  data %>%
    # exclude missing subgroup levels
    filter(!is.na(!!subgroup_var)) %>%
    group_by(!!subgroup_var) %>%
    summarise(
      n       = n(),
      n_DOAC  = sum(halt2 == "yes" & ran_trt == "DOAC", na.rm = TRUE),
      n_ASA   = sum(halt2 == "yes" & ran_trt == "ASA",  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(level = !!subgroup_var) %>%
    mutate(
      subgroup = subgroup_label,
      .before = 1
    )
}

summarise_overall <- function(data) {
  data %>%
    summarise(
      n       = n(),
      n_DOAC  = sum(halt2 == "yes" & ran_trt == "DOAC", na.rm = TRUE),
      n_ASA   = sum(halt2 == "yes" & ran_trt == "ASA",  na.rm = TRUE)
    ) %>%
    mutate(
      subgroup = "Overall",
      level    = "Overall",
      .before  = 1
    )
}

make_summary_table <- function(data, subgroup_vars) {
  # subgroup_vars: a list of symbols, e.g. vars(age, sex)
  
  overall_tbl <- summarise_overall(data)
  
  subgroup_tbls <- subgroup_vars %>%
    map(~ summarise_subgroup(data, !!.x)) %>%
    bind_rows()
  
  bind_rows(overall_tbl, subgroup_tbls) %>%
    select(subgroup, level, n, n_DOAC, n_ASA)
}

df_sum <- make_summary_table(eff_base,vars(ageF,sex,diabetesF,hypertensionF,
                                           renal_function,valve_type,
                                           pdF,ascF,frailty_status))
df_sum$subgroup[df_sum$subgroup=="Pd"] <- "Post-dilatation"
df_sum$subgroup[df_sum$subgroup=="Asc"] <- "Ascending aorta diam."
df_sum$subgroup[df_sum$subgroup=="Frailty Status"] <- "Frailty index"

df_sum <- df_sum %>%
  mutate(
    risk_DOAC  = n_DOAC  / n,
    risk_ASA = n_ASA / n
  )
df_sum$rr <- c(halt_main_ratio$estimate,halt_sub1$estimate,halt_sub2$estimate,
               halt_sub3$estimate,halt_sub4$estimate,halt_sub5$estimate,halt_sub6$estimate,
               halt_sub7$estimate,halt_sub8$estimate,halt_sub9$estimate)
df_sum$lo <- c(halt_main_ratio$conf.low,halt_sub1$conf.low,halt_sub2$conf.low,
               halt_sub3$conf.low,halt_sub4$conf.low,halt_sub5$conf.low,halt_sub6$conf.low,
               halt_sub7$conf.low,halt_sub8$conf.low,halt_sub9$conf.low)
df_sum$hi <- c(halt_main_ratio$conf.high,halt_sub1$conf.high,halt_sub2$conf.high,
               halt_sub3$conf.high,halt_sub4$conf.high,halt_sub5$conf.high,halt_sub6$conf.high,
               halt_sub7$conf.high,halt_sub8$conf.high,halt_sub9$conf.high)



## 3. Add “Age” and “Sex” header rows and prepare labels
df_with_headers <- df_sum %>%
  # add explicit header rows
  bind_rows(
    tibble(subgroup = "Age", level = "Age", is_header = TRUE),
    tibble(subgroup = "Sex", level = "Sex", is_header = TRUE),
    tibble(subgroup = "Diabetes", level = "Diabetes", is_header = TRUE),
    tibble(subgroup = "Hypertension", level = "Hypertension", is_header = TRUE),
    tibble(subgroup = "Renal function", level = "Renal function", is_header = TRUE),
    tibble(subgroup = "Valve type", level = "Valve type", is_header = TRUE),
    tibble(subgroup = "Post-dilatation", level = "Post-dilatation", is_header = TRUE),
    tibble(subgroup = "Ascending aorta diam.", level = "Ascending aorta diam.", is_header = TRUE),
    tibble(subgroup = "Frailty index", level = "Frailty index", is_header = TRUE)
  ) %>%
  mutate(
    is_header = if_else(
      is.na(is_header),
      subgroup %in% c("Age", "Sex","Diabetes","Hypertension","Renal function",
                      "Valve type","Post-dilatation","Ascending aorta diam.",
                      "Frailty index")& level == subgroup,
      is_header
    ),
    # first-column text
    row_lab = case_when(
      subgroup == "Overall" ~ "Overall",
      is_header             ~ as.character(subgroup),        # Age, Sex
      TRUE                  ~ paste0("   ", level)           # indent levels
    )
  ) %>%
  arrange(
    factor(subgroup, levels = c("Overall","Age", "Sex","Diabetes","Hypertension","Renal function",
                                "Valve type","Post-dilatation","Ascending aorta diam.",
                                "Frailty index")),
    desc(is_header),    # header first within subgroup
    level
  ) 
df_with_headers$row_factor <- factor(paste0(df_with_headers$subgroup,df_with_headers$level),
                                     levels=rev(paste0(df_with_headers$subgroup,df_with_headers$level)))

## 4. Text for counts / percentages
df_with_headers <- df_with_headers %>%
  mutate(
    lab_n = case_when(
      subgroup == "Overall" ~ sprintf("%d", n),
      is_header             ~ "",   # no numbers on headers
      TRUE                  ~ sprintf("%d", n)
    ),
    lab_trt = case_when(
      subgroup == "Overall" ~ sprintf("%d (%.1f%%)", n_DOAC, 100 * risk_DOAC),
      is_header             ~ "",
      TRUE                  ~ sprintf("%d (%.1f%%)", n_DOAC, 100 * risk_DOAC)
    ),
    lab_ctl = case_when(
      subgroup == "Overall" ~ sprintf("%d (%.1f%%)", n_ASA, 100 * risk_ASA),
      is_header             ~ "",
      TRUE                  ~ sprintf("%d (%.1f%%)", n_ASA, 100 * risk_ASA)
    )
  )

## 5. Forest plot with headers + side columns
HALTsubgroups <- ggplot() +
  # first column: subgroup labels (Overall, Age, Sex, levels)
  geom_text(
    data = df_with_headers,
    aes(
      y = row_factor, x = 0.01, label = row_lab, hjust=-5, #0.2
      fontface = if_else(is_header | subgroup == "Overall", "bold", "plain")
    ),
    hjust = 0, size = 7
  ) +
  # reference line
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  # CIs and points (exclude header rows)
  geom_errorbar(
    data = df_with_headers %>% filter(!is_header),
    aes(y = row_factor, xmin = lo, xmax = hi),
    height = 0.15, colour = "navy"
  ) +
  geom_point(
    data = df_with_headers %>% filter(!is_header),
    aes(y = row_factor, x = rr),
    shape = 15, size = 3, colour = "navy"
  ) +
  
  # side columns: N, Treat, Control
  geom_text(
    data = df_with_headers,
    aes(y = row_factor, x = 1.5, label = lab_n), #0.4
    hjust = 0, size = 7
  ) +
  geom_text(
    data = df_with_headers,
    aes(y = row_factor, x = 1.7, label = lab_trt), #2.5
    hjust = 0, size = 7
  ) +
  geom_text(
    data = df_with_headers,
    aes(y = row_factor, x = 2, label = lab_ctl),#3.5
    hjust = 0, size = 7
  ) +
  # headers for those text columns
  annotate("text", x = 1.5, y = Inf, label = "N",
           hjust = 0, vjust = -0.1, size = 8, fontface = "bold") +
  annotate("text", x = 1.7, y = Inf, label = "DOAC n(%)",
           hjust = 0, vjust = -0.1, size = 8, fontface = "bold") +
  annotate("text", x = 2, y = Inf, label = "ASA n(%)",
           hjust = 0, vjust = -0.1, size = 8, fontface = "bold") +
  # RR axis
  scale_x_continuous( #scale_x_log10
    breaks = c( 0.5, 1,1.5, 2),
    limits = c(0.01, 2.5), # c(0.2, 4.5)
    name   = "Risk Ratio (DOAC / ASA)"
  ) +
  coord_cartesian(clip = "off") + 
  scale_y_discrete(NULL) +
  theme_minimal(base_size = 11) +
  theme(
    text = element_text(size = 20),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.y       = element_blank(),
    axis.text.y        = element_blank(),  # we drew labels manually
    plot.margin        = margin(20, 5.5, 5.5, 10)
  )

HALTsubgroups

save(HALTsubgroups,halt_sub1,halt_sub2,halt_sub3,halt_sub4,
     halt_sub5,halt_sub6,halt_sub7,halt_sub8,halt_sub9,
     file="data/res/subgroups_prim_eff_tab.RData")




################## Old descriptive figures
# 
# # Descriptive statistics
# overallT <- eff_base%>% group_by(ran_trt) %>% 
#   summarise(n_noHALT=sum(halt=="no",na.rm=T),
#             n_hd=sum(halt2=="yes",na.rm=T), .groups = "drop_last")
# 
# ageT <- eff_base%>% group_by(ageF,ran_trt) %>% 
#   summarise(n_noHALT=sum(halt=="no",na.rm=T),n_HALT=sum(halt=="yes",na.rm=T),
#             n_death=sum(death=="yes",na.rm=T),
#             n_hd=sum(halt2=="yes",na.rm=T), .groups = "drop_last") %>%
#   mutate(n_tot=n_noHALT+n_hd,pct=paste(round(100*n_hd/n_tot,1),"%",sep="")) %>%
#   select(-c(n_hd,n_tot))
# ageT <- head(ageT) %>%  t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "at") 
# ageT$at[3:6] <- c("no HALT","HALT","Deaths","HALT & Deaths")
# 
# 
# sexT <- eff_base%>% group_by(sex,ran_trt) %>% 
#   summarise(n_noHALT=sum(halt=="no",na.rm=T),n_HALT=sum(halt=="yes",na.rm=T),
#             n_death=sum(death=="yes",na.rm=T),
#             n_hd=sum(halt2=="yes",na.rm=T), .groups = "drop_last") %>%
#   mutate(n_tot=n_noHALT+n_hd,pct=paste(round(100*n_hd/n_tot,1),"%",sep="")) %>%
#   select(-c(n_hd,n_tot))
# sexT <- head(sexT) %>%  t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "at") 
# sexT$at[3:6] <- c("no HALT","HALT","Deaths","HALT & Deaths")
# 
# diabT <- eff_base%>% filter (!is.na(diabetes)) %>% group_by(diabetes,ran_trt) %>% 
#   summarise(n_noHALT=sum(halt=="no",na.rm=T),n_HALT=sum(halt=="yes",na.rm=T),
#             n_death=sum(death=="yes",na.rm=T),
#             n_hd=sum(halt2=="yes",na.rm=T), .groups = "drop_last") %>%
#   mutate(n_tot=n_noHALT+n_hd,pct=paste(round(100*n_hd/n_tot,1),"%",sep="")) %>%
#   select(-c(n_hd,n_tot))
# diabT <- head(diabT) %>%  t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "at") 
# diabT$at[3:6] <- c("no HALT","HALT","Deaths","HALT & Deaths")
# 
# 
# hypT <- eff_base%>% filter (!is.na(hypertension)) %>% group_by(hypertension,ran_trt) %>% 
#   summarise(n_noHALT=sum(halt=="no",na.rm=T),n_HALT=sum(halt=="yes",na.rm=T),
#             n_death=sum(death=="yes",na.rm=T),
#             n_hd=sum(halt2=="yes",na.rm=T), .groups = "drop_last") %>%
#   mutate(n_tot=n_noHALT+n_hd,pct=paste(round(100*n_hd/n_tot,1),"%",sep="")) %>%
#   select(-c(n_hd,n_tot))
# hypT <- head(hypT) %>%  t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "at") 
# hypT$at[3:6] <- c("no HALT","HALT","Deaths","HALT & Deaths")
# 
# 
# renT <- eff_base%>% filter (!is.na(impaired_renal_function)) %>% group_by(impaired_renal_function,ran_trt) %>% 
#   summarise(n_noHALT=sum(halt=="no",na.rm=T),n_HALT=sum(halt=="yes",na.rm=T),
#             n_death=sum(death=="yes",na.rm=T),
#             n_hd=sum(halt2=="yes",na.rm=T), .groups = "drop_last") %>%
#   mutate(n_tot=n_noHALT+n_hd,pct=paste(round(100*n_hd/n_tot,1),"%",sep="")) %>%
#   select(-c(n_hd,n_tot))
# renT <- head(renT) %>%  t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "at") 
# renT$at[3:6] <- c("no HALT","HALT","Deaths","HALT & Deaths")
# 
# valT <- eff_base%>% filter (!is.na(valve_type)) %>% group_by(valve_type,ran_trt) %>% 
#   summarise(n_noHALT=sum(halt=="no",na.rm=T),n_HALT=sum(halt=="yes",na.rm=T),
#             n_death=sum(death=="yes",na.rm=T),
#             n_hd=sum(halt2=="yes",na.rm=T), .groups = "drop_last") %>%
#   mutate(n_tot=n_noHALT+n_hd,pct=paste(round(100*n_hd/n_tot,1),"%",sep="")) %>%
#   select(-c(n_hd,n_tot))
# valT <- head(valT) %>%  t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "at") 
# valT$at[3:6] <- c("no HALT","HALT","Deaths","HALT & Deaths")
# 
# posT <- eff_base%>% filter (!is.na(tavi_post_dila)) %>% group_by(tavi_post_dila,ran_trt) %>% 
#   summarise(n_noHALT=sum(halt=="no",na.rm=T),n_HALT=sum(halt=="yes",na.rm=T),
#             n_death=sum(death=="yes",na.rm=T),
#             n_hd=sum(halt2=="yes",na.rm=T), .groups = "drop_last") %>%
#   mutate(n_tot=n_noHALT+n_hd,pct=paste(round(100*n_hd/n_tot,1),"%",sep="")) %>%
#   select(-c(n_hd,n_tot))
# posT <- head(posT) %>%  t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "at") 
# posT$at[3:6] <- c("no HALT","HALT","Deaths","HALT & Deaths")
# 
# supT <- eff_base%>% filter (!is.na(tavi_supra_pos)) %>% group_by(tavi_supra_pos,ran_trt) %>% 
#   summarise(n_noHALT=sum(halt=="no",na.rm=T),n_HALT=sum(halt=="yes",na.rm=T),
#             n_death=sum(death=="yes",na.rm=T),
#             n_hd=sum(halt2=="yes",na.rm=T), .groups = "drop_last") %>%
#   mutate(n_tot=n_noHALT+n_hd,pct=paste(round(100*n_hd/n_tot,1),"%",sep="")) %>%
#   select(-c(n_hd,n_tot))
# supT <- head(supT) %>%  t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "at") 
# supT$at[3:6] <- c("no HALT","HALT","Deaths","HALT & Deaths")
# 
# ascT <- eff_base%>% filter (!is.na(ascF)) %>% group_by(ascF,ran_trt) %>% 
#   summarise(n_noHALT=sum(halt=="no",na.rm=T),n_HALT=sum(halt=="yes",na.rm=T),
#             n_death=sum(death=="yes",na.rm=T),
#             n_hd=sum(halt2=="yes",na.rm=T), .groups = "drop_last") %>%
#   mutate(n_tot=n_noHALT+n_hd,pct=paste(round(100*n_hd/n_tot,1),"%",sep="")) %>%
#   select(-c(n_hd,n_tot))
# ascT <- head(ascT) %>%  t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "at") 
# ascT$at[3:6] <- c("no HALT","HALT","Deaths","HALT & Deaths")
# 
# fraT <- eff_base%>% filter (!is.na(frailty_status)) %>% group_by(frailty_status,ran_trt) %>% 
#   summarise(n_noHALT=sum(halt=="no",na.rm=T),n_HALT=sum(halt=="yes",na.rm=T),
#             n_death=sum(death=="yes",na.rm=T),
#             n_hd=sum(halt2=="yes",na.rm=T), .groups = "drop_last") %>%
#   mutate(n_tot=n_noHALT+n_hd,pct=paste(round(100*n_hd/n_tot,1),"%",sep="")) %>%
#   select(-c(n_hd,n_tot))
# fraT <- head(fraT) %>%  t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "at") 
# fraT$at[3:6] <- c("no HALT","HALT","Deaths","HALT & Deaths")