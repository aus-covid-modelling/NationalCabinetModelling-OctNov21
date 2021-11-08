# Script to produce linelist of arrival spillovers
# 
# David J Price
# 4 October, 2021
# 
library(data.table)
library(tibble)
library(dplyr)
library(readr)

#WA_family_rate <- 18.19/100
#NSW_family_rate <- 17.9/100
percentage_family <- 18.19/100
arrivals_perweek <- 130265

nonfamily_arrivals <- round((1-percentage_family)*arrivals_perweek)
family_arrivals <- round(percentage_family*arrivals_perweek)

# setwd("~/Desktop/COVID_Misc/WP3/")
source("arrivals_functions_mychanges.R")

## Where are the arrivals scenarios files ----
quarantine.output.folder <- "./quarantine_output/"

# Files for each arrivals pathway (for arrivals)#these should be pathway1, pathway2, etc... just repeating for the single file
# 

arrivals.files <- c("7Day_Hotel_Vaccinated/Traveller_breach.csv", 
                    "14Day_Hotel_Vaccinated/Traveller_breach.csv",
                    "7Day_Home_Vaccinated/Traveller_breach.csv",
                    "14Day_Home_Vaccinated/Traveller_breach.csv",
                    "NoQuarantine/Traveller_breach.csv",
                    "14Day_Hotel_Unvaccinated/Traveller_breach.csv",
                    "Family_7Day_Hotel_Vaccinated/Traveller_breach.csv", 
                    "Family_14Day_Hotel_Vaccinated/Traveller_breach.csv",
                    "Family_7Day_Home_Vaccinated/Traveller_breach.csv",
                    "Family_14Day_Home_Vaccinated/Traveller_breach.csv",
                    "Family_NoQuarantine/Traveller_breach.csv",
                    "Family_14Day_Hotel_Unvaccinated/Traveller_breach.csv"
                    )

# Files for each arrivals pathway (for workers)
workers.files <- c("7Day_Hotel_Vaccinated/Worker_breach.csv", 
                   "14Day_Hotel_Vaccinated/Worker_breach.csv",
                   "7Day_Home_Vaccinated/Worker_breach.csv",
                   "14Day_Home_Vaccinated/Worker_breach.csv",
                   "NoQuarantine/Worker_breach.csv",
                   "14Day_Hotel_Unvaccinated/Worker_breach.csv",
                   "Family_7Day_Hotel_Vaccinated/Worker_breach.csv", 
                   "Family_14Day_Hotel_Vaccinated/Worker_breach.csv",
                   "Family_7Day_Home_Vaccinated/Worker_breach.csv",
                   "Family_14Day_Home_Vaccinated/Worker_breach.csv",
                   "Family_NoQuarantine/Worker_breach.csv",
                   "Family_14Day_Hotel_Unvaccinated/Worker_breach.csv"
                    )

# workers.files <- c("Run_7_D7_T2_pV1_Vi36_Vs67_C9_Home_2021_10_08/Worker_breach.csv", 
#                    "Run_17_D14_T3_pV_C1_Hotel_2021_10_08/Worker_breach.csv",
#                    "Run_18_D14_T3_pV1_Vi36_Vs67_C1_Hotel_2021_10_08/Worker_breach.csv",
#                    "Run_1_D14_T3_pV1_Vi36_Vs67_C_Home_2021_10_08/Worker_breach.csv",
#                    "Run_15_D7_T2_pV_C1_Hotel_2021_10_08/Worker_breach.csv"
# )

#### Run a single example
#   generate_one_linelist(arrivals = arrivals.files[1], workers = workers.files[1],
#                         quarantine.folder = quarantine.output.folder, n.arrivals = num.arrivals[1],
#                         t.max = 180)
####


# Create a data.frame with each scenario (the number of arrivals through each pathway) as rows

## NOTE will need to add the prevalence in each pathway into this table and pass to 
## 'create_linelist' to scale appropriately. Baseline prevalence is 0.01 (check) ----

# Number of arrivals through each pathway for each scenario - make sure that this rate is correct. (100 people per 7 Days or 14 days. )
num.arrivals.mat <- tribble(
  ~scenario, ~Hotel7, ~Hotel14, ~Home7, ~Home14, ~NoneV, ~Hotel14Unvax, ~Hotel7Family, ~Hotel14Family, ~Home7Family, ~Home14Family, ~NoneVFamily, ~Hotel14UnvaxFamily,
  "NSW_Arrivals40_Hotel7", nonfamily_arrivals, 0, 0, 0, 0, 0, family_arrivals, 0, 0, 0, 0, 0,
  "NSW_Arrivals40_Hotel14", 0, nonfamily_arrivals, 0, 0, 0, 0, 0, family_arrivals, 0, 0, 0, 0,
  "NSW_Arrivals40_Home7", 0, 0, nonfamily_arrivals, 0, 0, 0, 0, 0, family_arrivals, 0, 0, 0,
  "NSW_Arrivals40_Home14", 0, 0, 0, nonfamily_arrivals, 0, 0, 0, 0, 0, family_arrivals, 0, 0,
  "NSW_Arrivals40_NoneV", 0, 0, 0, nonfamily_arrivals, 0, 0, 0, 0, 0, family_arrivals, 0, 0,
  "NSW_Arrivals40_Hotel14Unvax", 0, 0, 0, 0, nonfamily_arrivals, 0, 0, 0, 0, 0, family_arrivals, 0,
  "NSW_Arrivals40_Home7", 0, 0, 2*nonfamily_arrivals, 0, 0, 0, 0, 0, 2*family_arrivals, 0, 0, 0
)

# Prevalence in arrivals through each pathway for each scenario
path.prev.mat <- tribble(
  ~scenario, ~Hotel7, ~Hotel14, ~Home7, ~Home14, ~NoneV, ~Hotel14Unvax, ~Hotel7Family, ~Hotel14Family, ~Home7Family, ~Home14Family, ~NoneVFamily, ~Hotel14UnvaxFamily,
  1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
  2, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
  3, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
  4, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
  5, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
  6, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
  7, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
) %>% mutate(scenario=num.arrivals.mat$scenario)


path.family.mat  <- tribble(
  ~scenario, ~Hotel7, ~Hotel14, ~Home7, ~Home14, ~NoneV, ~Hotel14Unvax, ~Hotel7Family, ~Hotel14Family, ~Home7Family, ~Home14Family, ~NoneVFamily, ~Hotel14UnvaxFamily,
  1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
  2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
  3, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
  4, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
  5, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
  6, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
  7, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
) %>% mutate(scenario=num.arrivals.mat$scenario)

# Generate n.sims linelists of spillovers for each scenario 
n.sims <- 200
# Default/baseline prevalence used in HQ model
baseline.prevalence <- 0.01
# Might have to think about vaccination for the import prevalence IN ORIGIN COUNTRY. 
# This is explicitly the prevalence in an unvaccinated population. 
# Vaccine efficacy against infection might have to incorporated. 

# Expand number of arrivals and prevalence for each scenario n.sims times
all.sims.arrivals <- num.arrivals.mat[rep(seq_len(nrow(num.arrivals.mat)), each = n.sims), ] %>% 
  mutate(sim = rep(1:n.sims, times = nrow(num.arrivals.mat)))

all.sims.prevs <- path.prev.mat[rep(seq_len(nrow(path.prev.mat)), each = n.sims), ] %>%
  mutate(sim = rep(1:n.sims, times = nrow(path.prev.mat)))

all.sims.family <- path.family.mat[rep(seq_len(nrow(path.family.mat)), each = n.sims), ] %>%
 mutate(sim = rep(1:n.sims, times = nrow(path.family.mat)))


# set.seed(1)
# Produces a list with the file for each scenario, n.sims times
tic <- Sys.time()
all.out <- lapply(X = 1:nrow(all.sims.arrivals), 
                  function(i){create_linelist(num.arrivals = all.sims.arrivals[i,] %>% select(-scenario, -sim) %>% as.numeric(), 
                                              arrivals.files = arrivals.files, workers.files = workers.files, 
                                              pathway.prevalence = all.sims.prevs[i,] %>% select(-scenario, -sim) %>% as.numeric(),t.max = 500,VE_transmission = 0.36, family.pathway = all.sims.family[i,] %>% select(-scenario,-sim) %>% as.numeric())})
(toc <- Sys.time() - tic)



# all.out <- lapply(X = 1:nrow(all.sims.arrivals), 
#                   function(i){create_linelist(num.arrivals = all.sims.arrivals[i,] %>% select(-scenario, -sim) %>% as.numeric(), 
#                                               arrivals.files = arrivals.files, workers.files = workers.files, 
#                                               pathway.prevalence = all.sims.prevs[i,] %>% select(-scenario, -sim) %>% as.numeric(),t.max = 500,VE_transmission = 0.36, 
#                     family.pathway = all.sims.family[i,] %>% select(-scenario,-sim) %>% as.numeric())})
# (toc <- Sys.time() - tic)

# Saves each simulated linelist in afolder/file named according to scenario and sim no.
mapply(FUN = function(out, sc, sim){
  folder.name <- paste0("arrivals_linelist/scenario",sc,"/")
  dir.create(folder.name, showWarnings = FALSE)
  write_csv(out,file = paste0(folder.name,"arrivals_linelist_scenario",sc,"_sim",sim,".csv"))
  },
                       out = all.out, sc = all.sims.arrivals$scenario, sim = all.sims.arrivals$sim)







