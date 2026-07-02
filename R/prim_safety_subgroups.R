#######################
##    ACASA-TAVI     ##
##    C. Cunen       ##
#######################

## Primary safety endpoint - subgroup analysis

# Descriptive statistics
# Adjusted, with imputation

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
saf <- read_rds("data/td/cso_td.rds")

adsl <- adsl %>% select(-c(site,ran_date))
baseline <- baseline %>% select(-site)

saf <- saf %>% left_join(adsl,by="subjectid")
saf_base <- saf %>% left_join(baseline,by="subjectid")

saf_base <- saf_base %>% mutate(smoke2 = 
                                  fct_collapse(smoke,
                                               smoker = c("On a daily basis", "Smokes occasionally")))


# Define subgroup factors
saf_base <- saf_base %>%
  mutate(ageF = factor(if_else(age < 75, "<75", ">= 75"),
                       levels = c("<75", ">= 75")))
saf_base <- saf_base %>%
  mutate(renal_function = factor(if_else(GFR < 30, "impaired", "preserved"),
                                 levels = c("impaired", "preserved")))
saf_base <- saf_base %>%
  mutate(ascF = factor(if_else(asc_aorta_diam < 30, "<30 mm", ">= 30 mm"),
                       levels = c("<30 mm", ">= 30 mm")))

saf_base$diabetesF <- factor(ifelse(saf_base$diabetes=="1","yes","no"))
saf_base$hypertensionF <- factor(ifelse(saf_base$hypertension=="1","yes","no"))
saf_base$pdF <- factor(ifelse(saf_base$tavi_post_dila=="Yes","yes","no"))

# Subset for modelling
saf_base_mod <- saf_base %>% select(safety,ran_trt,age,sex,cad,diabetes,hypertension,prev_stroke,GFR,
                                    ageF,sex,diabetesF,hypertensionF,renal_function,valve_type,
                                    pdF,ascF,frailty_status) 
p <- 8

idc <- which(apply(saf_base_mod[, 3:(p+1)], 1, function(x) any(is.na(x))))
saf_base_mod[idc,3:(p+1)]

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
saf_base_mod[idc,3:(p+1)]

# Adjusted analysis - multiple imputation with imputation model - containing TREATMENT
library(mice)

dat_mice <- mice(saf_base_mod[,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                 "diabetes","hypertension","GFR")],m=20,printFlag = F,seed=11254)
saf_mod_mice <- with(dat_mice,glm(safety~ran_trt+age+sex+cad+prev_stroke+
                                    diabetes+hypertension+GFR,family="binomial"))
safety_main_diff <- avg_comparisons(saf_mod_mice,variables=list(ran_trt = c("ASA", "DOAC")),
                                    equivalence=c(NA,0.119))


# age
dat_mice <- mice(saf_base_mod[,c("safety","ran_trt","sex","cad","prev_stroke",
                                 "diabetes","hypertension","GFR","ageF")],m=20,printFlag = F,seed=11254)
mod_mice1 <- with(dat_mice,glm(safety~ran_trt*ageF+sex+cad+prev_stroke+
                                 diabetes+hypertension+GFR,family="binomial"))
saf_sub1 <- avg_comparisons(mod_mice1,variables=list(ran_trt = c("ASA", "DOAC")),by=c("ageF"))
plot(dat_mice)

# sex
dat_mice <- mice(saf_base_mod[,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                 "diabetes","hypertension","GFR")],m=20,printFlag = F,seed=11254)
mod_mice2 <- with(dat_mice,glm(safety~ran_trt*sex+age+sex+cad+prev_stroke+
                                 diabetes+hypertension+GFR,family="binomial"))
saf_sub2 <- avg_comparisons(mod_mice2,variables=list(ran_trt = c("ASA", "DOAC")),by=c("sex"))

# diabetes
idm <- which(is.na(saf_base_mod$diabetesF))
dat_mice <- mice(saf_base_mod[-idm,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                     "hypertension","GFR","diabetesF")],m=20,printFlag = F,seed=11254)
mod_mice3 <- with(dat_mice,glm(safety~ran_trt*diabetesF+age+sex+cad+prev_stroke+hypertension+GFR,family="binomial"))
saf_sub3 <- avg_comparisons(mod_mice3,variables=list(ran_trt = c("ASA", "DOAC")),by=c("diabetesF"))
#plot(dat_mice)

# hypertension
idm <- which(is.na(saf_base_mod$hypertensionF))
dat_mice <- mice(saf_base_mod[-idm,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                 "diabetes","GFR","hypertensionF")],m=20,printFlag = F,seed=11254)
mod_mice4 <- with(dat_mice,glm(safety~ran_trt*hypertensionF+age+sex+cad+prev_stroke+
                                 diabetes+GFR,family="binomial"))
saf_sub4 <- avg_comparisons(mod_mice4,variables=list(ran_trt = c("ASA", "DOAC")),by=c("hypertensionF"))

# renal function
dat_mice <- mice(saf_base_mod[,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                 "diabetes","hypertension","renal_function")],m=20,printFlag = F,seed=11254)
mod_mice5 <- with(dat_mice,glm(safety~ran_trt*renal_function+age+sex+cad+prev_stroke+
                                 diabetes+hypertension,family="binomial"))
saf_sub5 <- avg_comparisons(mod_mice5,variables=list(ran_trt = c("ASA", "DOAC")),by=c("renal_function"))

# valve_type
idm <- which(is.na(saf_base_mod$valve_type))
dat_mice0 <- mice(saf_base_mod[-idm,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                  "diabetes","hypertension","GFR","valve_type")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "valve_type"] <- 0 # to ensure that imputation model does not use subgroup 

dat_mice <- mice(saf_base_mod[-idm,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                 "diabetes","hypertension","GFR","valve_type")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice6 <- with(dat_mice,glm(safety~ran_trt*valve_type+age+sex+cad+prev_stroke+
                                 diabetes+hypertension+GFR,family="binomial"))
saf_sub6 <- avg_comparisons(mod_mice6,variables=list(ran_trt = c("ASA", "DOAC")),by=c("valve_type"))

# post-dilatation
idm <- which(is.na(saf_base_mod$pdF))
dat_mice0 <- mice(saf_base_mod[-idm,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                  "diabetes","hypertension","GFR","pdF")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "pdF"] <- 0 # to ensure that imputation model does not subgroup 

dat_mice <- mice(saf_base_mod[-idm,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                 "diabetes","hypertension","GFR","pdF")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice7 <- with(dat_mice,glm(safety~ran_trt*pdF+age+sex+cad+prev_stroke+
                                 diabetes+hypertension+GFR,family="binomial"))
saf_sub7 <- avg_comparisons(mod_mice7,variables=list(ran_trt = c("ASA", "DOAC")),by=c("pdF"))

# ascending aorta diameter
dat_mice0 <- mice(saf_base_mod[,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                  "diabetes","hypertension","GFR","ascF")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "ascF"] <- 0 # to ensure that imputation model does not use subgroup 

dat_mice <- mice(saf_base_mod[,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                 "diabetes","hypertension","GFR","ascF")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice8 <- with(dat_mice,glm(safety~ran_trt*ascF+age+sex+cad+prev_stroke+
                                 diabetes+hypertension+GFR,family="binomial"))
saf_sub8 <- avg_comparisons(mod_mice8,variables=list(ran_trt = c("ASA", "DOAC")),by=c("ascF"))

# frailty
dat_mice0 <- mice(saf_base_mod[,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                  "diabetes","hypertension","GFR","frailty_status")],maxit=0)
pred <- dat_mice0$predictorMatrix
pred[, "frailty_status"] <- 0 # to ensure that imputation model does not subgroup 

dat_mice <- mice(saf_base_mod[,c("safety","ran_trt","age","sex","cad","prev_stroke",
                                 "diabetes","hypertension","GFR","frailty_status")],m=20,printFlag = F,seed=11254,
                 predictorMatrix = pred)
mod_mice9 <- with(dat_mice,glm(safety~ran_trt*frailty_status+age+sex+cad+prev_stroke+
                                 diabetes+hypertension+GFR,family="binomial"))
saf_sub9 <- avg_comparisons(mod_mice9,variables=list(ran_trt = c("ASA", "DOAC")),by=c("frailty_status"))

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
      n_DOAC  = sum(safety == "1" & ran_trt == "DOAC", na.rm = TRUE),
      n_ASA   = sum(safety == "1" & ran_trt == "ASA",  na.rm = TRUE),
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
      n_DOAC  = sum(safety == "1" & ran_trt == "DOAC", na.rm = TRUE),
      n_ASA   = sum(safety == "1" & ran_trt == "ASA",  na.rm = TRUE)
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

df_sum <- make_summary_table(saf_base,vars(ageF,sex,diabetesF,hypertensionF,
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
df_sum$rd <- c(safety_main_diff$estimate,saf_sub1$estimate,saf_sub2$estimate,
               saf_sub3$estimate,saf_sub4$estimate,saf_sub5$estimate,saf_sub6$estimate,
               saf_sub7$estimate,saf_sub8$estimate,saf_sub9$estimate)
df_sum$lo <- c(safety_main_diff$conf.low,saf_sub1$conf.low,saf_sub2$conf.low,
               saf_sub3$conf.low,saf_sub4$conf.low,saf_sub5$conf.low,saf_sub6$conf.low,
               saf_sub7$conf.low,saf_sub8$conf.low,saf_sub9$conf.low)
df_sum$hi <- c(safety_main_diff$conf.high,saf_sub1$conf.high,saf_sub2$conf.high,
               saf_sub3$conf.high,saf_sub4$conf.high,saf_sub5$conf.high,saf_sub6$conf.high,
               saf_sub7$conf.high,saf_sub8$conf.high,saf_sub9$conf.high)



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
SAFsubgroups <- ggplot() +
  # first column: subgroup labels (Overall, Age, Sex, levels)
  geom_text(
    data = df_with_headers,
    aes(
      y = row_factor, x = -0.5, label = row_lab, #0.2
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
    aes(y = row_factor, x = rd),
    shape = 15, size = 3, colour = "navy"
  ) +
  
  # side columns: N, Treat, Control
  geom_text(
    data = df_with_headers,
    aes(y = row_factor, x = 0.22, label = lab_n), #0.4
    hjust = 0, size = 7
  ) +
  geom_text(
    data = df_with_headers,
    aes(y = row_factor, x = 0.26, label = lab_trt), #2.5
    hjust = 0, size = 7
  ) +
  geom_text(
    data = df_with_headers,
    aes(y = row_factor, x = 0.4, label = lab_ctl),#3.5
    hjust = 0, size = 7
  ) +
  # headers for those text columns
  annotate("text", x = 0.22, y = Inf, label = "N",
           hjust = 0, vjust = -0.1, size = 8, fontface = "bold") +
  annotate("text", x = 0.26, y = Inf, label = "DOAC n(%)",
           hjust = 0, vjust = -0.1, size = 8, fontface = "bold") +
  annotate("text", x = 0.4, y = Inf, label = "ASA n(%)",
           hjust = 0, vjust = -0.1, size = 8, fontface = "bold") +
  # RD axis
  scale_x_continuous( #scale_x_log10
    breaks = c( -0.3,-0.2,-0.1, 0,0.1,0.2,0.3),
    limits = c(-0.5, 0.6), # c(0.2, 4.5)
    name   = "Risk Difference (DOAC - ASA)"
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
    plot.margin        = margin(20, 5.5, 5.5, 5.5)
  )

SAFsubgroups

save(SAFsubgroups,saf_sub1,saf_sub2,saf_sub3,saf_sub4,
     saf_sub5,saf_sub6,saf_sub7,saf_sub8,saf_sub9,
     file="data/res/subgroups_prim_saf_tab.RData")


