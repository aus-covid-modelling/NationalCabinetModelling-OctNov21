#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
folder = args[1]

# setwd("./outputs/Scenario1")
setwd(folder)
filenames <- list.files(pattern = "sim_number*")
library(dplyr)
library(readr)


summary <- tibble()
i <- 1
for(file in filenames){
  sim <- read_csv(file)
  if (nrow(sim)==0)
	  next
  summary <-bind_rows(sim %>% mutate(Date = cut(`Time of symptom onset`,breaks = seq(0,500,by = 1),labels = FALSE),Q= Cluster!=0)%>% group_by(Q, Date ) %>% summarise( n = n()) %>% arrange(Date) %>% mutate(Sim = i),summary)
  i <- i +1

}
outfile <- gsub(pattern="./outputs/", replacement="", x=folder)
write_csv(summary,paste0(outfile, "_CompressedSymptomOnset.csv"))
# 
# sim <- read_csv("sim_number_99.csv")
# sim <- sim_number_9 %>% mutate(Q = Cluster!=0)
# # Group by symptom onset. 
# i <- 1
# 
#   
# # summary <-sim %>% mutate(AgeBracket = cut(Age,breaks = c(seq(0,80,by = 5),Inf),right =FALSE),Date = cut(`Time of symptom onset`,breaks = seq(0,500,by = 1),labels = FALSE),Q= Cluster!=0) %>% 
#   # group_by(`Vaccine at infection`, `Doses at infection`, Q, AgeBracket, Date ) %>% summarise( n = n()) %>% arrange(Date) %>% mutate(Sim = i)
# 
# 
# ggplot(sim,aes(x = `Time of symptom onset`)) + geom_histogram(binwidth= 1) + facet_wrap(~Q)
# 
