#######################
##    ACASA-TAVI     ##
##    C. Cunen       ##
#######################

# Source all these scripts when data is updated



####################################
# 

source("R/make_raw.R")

rm(list = objects())
source("R/make_td_ae.R")
rm(list = objects())
source("R/make_td_baseline.R")
rm(list = objects())
source("R/make_td_dates.R")
rm(list = objects())
source("R/make_td_physlab12.R")
rm(list = objects())
source("R/make_td_adherence.R")
rm(list = objects())
source("R/make_td_prim_efficacy.R")
rm(list = objects())
source("R/make_td_prim_safety.R")
rm(list = objects())
source("R/make_hospitalizations_td.R")
rm(list = objects())
source("R/make_td_kccq.R")
rm(list = objects())
source("R/make_td_key_secondary.R")
rm(list = objects())
source("R/make_td_euroscore.R")
rm(list = objects())

source("R/make_adran.R")
rm(list = objects())
source("R/make_adsl.R")
rm(list = objects())

source("R/make_disp_tabs.R")
rm(list = objects())
source("R/make_baseline_tables.R")
rm(list = objects())
source("R/make_ae_tables.R")
rm(list = objects())

source("R/prim_efficacy_analysis.R")
rm(list = objects())

source("R/prim_safety_analysis.R")
rm(list = objects())

source("R/prim_safety_switching_figure.R")
rm(list = objects())

source("R/prim_efficacy_supp.R")
rm(list = objects())

source("R/prim_safety_supp2.R")
rm(list = objects())


source("R/prim_efficacy_subgroups.R")
rm(list = objects())

source("R/prim_safety_subgroups.R")
rm(list = objects())

source("R/key_secondary.R")
rm(list = objects())
