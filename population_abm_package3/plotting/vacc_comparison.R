# Script to pull individual figures together for different vaccine coverage thresholds
# 
# 
# David J Price
# 22 October, 2021


## FIGURE 1 ----
the.state <- "NSW"
baseline <- TRUE
the.coverage <- c(70, 80, 90)
arrival.volume <- "40"

scenarios <- paste0("Arrivals",arrival.volume,"_Home7_",the.state,the.coverage)

scenario.list <- data.frame(scenarios) %>% separate(scenarios, into = c("Volume","Quarantine","Vaccine coverage")) %>% 
  mutate(Volume = gsub(pattern = "Arrivals",replacement = "", x = Volume)) %>% 
  mutate(`Vaccine coverage` = gsub(pattern = the.state, replacement = "", x = `Vaccine coverage`)) %>% 
  mutate(Quarantine = case_when(
    Quarantine == "NoneV" ~ "No quarantine (vacc)",
    Quarantine == "Hotel7" ~ "7-day Hotel (vacc)",
    Quarantine == "Hotel14" ~ "14-day Hotel (vacc)",
    Quarantine == "Home7" ~ "7-day Home, 90% Comp. (vacc)",
    Quarantine == "Home14" ~ "14-day Home, 90% Comp. (vacc)",
    Quarantine == "Hotel14Unvax" ~ "14-day Hotel (unvacc)"
  ))


if(isTRUE(baseline)){
  scenarios <- paste0(scenarios, "_baseline")
}

the.folder.names <- paste0("scenario", the.state, "_", scenarios, "_CompressedSymptomOnset")
output.folder <- "./output_25Oct2021/"

the.folders <- paste0(output.folder, the.folder.names,"/")

plot.list.arr <- list()
plot.list.tot <- list()
for (i in 1:length(the.folders)){
  load(file = paste0(the.folders[i],the.folder.names[i], "_plots.Rdata"))
  
  
  the.arr.plot <- arr.inf.plot  
  if(i == 1){
    the.arr.plot <- the.arr.plot + ggtitle(label = "By Cluster", subtitle = paste0("Vol: ",scenario.list$Volume[i],"%; ", 
                                                                                   "Quar: ", scenario.list$Quarantine[i],";\n", 
                                                                                   "Vacc: ", scenario.list$`Vaccine coverage`[i],"%"))
  } else{
    the.arr.plot <- the.arr.plot + ggtitle(label = "", subtitle = paste0("Vol: ",scenario.list$Volume[i],"%; ", 
                                                                         "Quar: ", scenario.list$Quarantine[i],";\n", 
                                                                         "Vacc: ", scenario.list$`Vaccine coverage`[i],"%"))
  }
  
  the.total.plot <- total.inf.plot  
  if(i == 1){
    the.total.plot <- the.total.plot + ggtitle(label = "", subtitle = paste0("Vol: ",scenario.list$Volume[i],"%; ", 
                                                                                       "Quar: ", scenario.list$Quarantine[i],";\n", 
                                                                                       "Vacc: ", scenario.list$`Vaccine coverage`[i],"%"))
  } else{
    the.total.plot <- the.total.plot + ggtitle(label = "", subtitle = paste0("Vol: ",scenario.list$Volume[i],"%; ", 
                                                                             "Quar: ", scenario.list$Quarantine[i],";\n", 
                                                                             "Vacc: ", scenario.list$`Vaccine coverage`[i],"%"))
  }
  
  
  plot.list.arr[[i]] <- the.arr.plot 
  plot.list.tot[[i]] <- the.total.plot 
  
  # assign(paste0(scenarios[i], ".arr.plot"), arr.inf.plot )
  # assign(paste0(scenarios[i], ".total.plot"), total.inf.plot )
}

p1 <- cowplot::plot_grid( plotlist = plot.list.arr, nrow = 1)
p2 <- cowplot::plot_grid( plotlist = plot.list.tot, nrow = 1)


if(isTRUE(baseline)){
  ggsave(paste0("./output_25Oct2021/vacc_coverage_comparison_",the.state,"_baseline_arrivals_plot.png"), plot = p1, height = 6, width = 15, units = "in", bg = "white")
  ggsave(paste0("./output_25Oct2021/vacc_coverage_comparison_",the.state,"_baseline_total_plot.png"), plot = p2, height= 6, width = 15, units = "in", bg = "white")
} else{
  ggsave(paste0("./output_25Oct2021/vacc_coverage_comparison_",the.state,"_arrivals_plot.png"), plot = p1, height = 6, width = 15, units = "in", bg = "white")
  ggsave(paste0("./output_25Oct2021/vacc_coverage_comparison_",the.state,"_total_plot.png"), plot = p2, height= 6, width = 15, units = "in", bg = "white")
}




