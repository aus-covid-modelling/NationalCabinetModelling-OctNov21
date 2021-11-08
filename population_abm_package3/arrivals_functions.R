# Functions for producing linelist of arrival spillovers
# 
# David J Price
# 4 October, 2021
# 


## Function to create linelist for each pathway (and numbers through that pathway) and combine ----
create_linelist <- function(num.arrivals, arrivals.files, workers.files, pathway.prevalence, t.max = 180, base.n.arrivals = 100, VE_transmission= 0.6, baseline.prevalence = 0.01){
  lapply(X = 1:length(num.arrivals), 
         FUN = function(i){generate_one_linelist(arrivals = arrivals.files[i], workers = workers.files[i], 
                                                 quarantine.folder = quarantine.output.folder, n.arrivals = num.arrivals[i], pathway.prevalence = pathway.prevalence[i],
                                                 t.max = t.max, base.n.arrivals = base.n.arrivals, VE_transmission = VE_transmission, baseline.prevalence)}) %>% 
    bind_rows() %>% arrange(total.time)
}


## This function creates a single linelist of spillovers from the corresponding arrival and worker files ----
generate_one_linelist <- function(arrivals, workers = NA, quarantine.folder, n.arrivals, pathway.prevalence, t.max, base.n.arrivals = 100, VE_transmission = 0.66, baseline.prevalence = 0.01){
  # arrivals and workers are the corresponding file names
  # set.seed(1)
  # n.arrivals must be scaled to be the appropriate arrivals PER unit time. 
  if(n.arrivals>0){
    scaling.factor <- n.arrivals / base.n.arrivals * (pathway.prevalence / baseline.prevalence) 
    # scale inter-event times according to numbers of arrivals
    # applies to workers too, assuming that the ratio of arrivals:workers remains constant
    
    # arrivals.dat <- suppressMessages(read_csv(paste0(quarantine.folder,arrivals)))
    arrivals.dat <- suppressMessages(fread(file = paste0(quarantine.folder,arrivals))) %>% as_tibble()
    
    # arrivals.times <-  arrivals.dat %>% 
    #   select(time_discharged, exposure_days, FoI_community, vaccinated) %>% 
    #   mutate(time.diff = diff(c(0,time_discharged))/scaling.factor) %>% 
    #   sample_n(size = 2000, replace = TRUE) %>% 
    #   mutate(total.time = cumsum(time.diff)) %>% 
    #   filter(total.time <= t.max)
    
    arrivals.times <-  arrivals.dat %>% 
      select(time_discharged, exposure_days, FoI_community, vaccinated) %>% 
      mutate(time.diff = diff(c(0,time_discharged))/scaling.factor) %>% 
      sample_n(size = 2000, replace = TRUE) %>% 
      mutate(total.time = cumsum(time.diff))
    
    # while(max(arrivals.times$total.tim)<t.max){
    #     arrivals.times <- arrivals.times %>% bind_rows(.,arrivals.dat %>% 
    #     select(time_discharged, exposure_days, FoI_community, vaccinated) %>% 
    #     mutate(time.diff = diff(c(0,time_discharged))/scaling.factor) %>% 
    #     sample_n(size = 2000, replace = TRUE) %>% 
    #     mutate(total.time = cumsum(time.diff)))
    # }
    # 
    arrivals.times <- arrivals.times %>% filter(total.time <= t.max)
    
    if(!is.na(workers)){
      # workers.dat <- suppressMessages(read_csv(paste0(quarantine.folder,workers))) 
      workers.dat <- suppressMessages(fread(file = paste0(quarantine.folder,workers))) %>% 
        as_tibble() %>% 
        slice(-1) #removing first entry as it looks like this is all 0's
      if(nrow(workers.dat)==0){
        workers.times <- my.empty.df
      }else{
      # workers.times <- workers.dat %>% 
      #   select(time_discharged, exposure_days, FoI_community, vaccinated) %>% 
      #   mutate(time.diff = diff(c(0,time_discharged))/scaling.factor) %>% 
      #   sample_n(size = 20000, replace = TRUE) %>% 
      #   mutate(total.time = cumsum(time.diff)) %>% 
      #   filter(total.time <= t.max)
      #   
      workers.times <- workers.dat %>% 
        select(time_discharged, exposure_days, FoI_community, vaccinated) %>% 
        mutate(time.diff = diff(c(0,time_discharged))/scaling.factor) %>% 
        sample_n(size = 20000, replace = TRUE) %>% 
        mutate(total.time = cumsum(time.diff))
      
        # while(max(workers.times) < t.max){
        #     workers.times <- workers.times %>% bind_rows(.,workers.dat %>% 
        #     select(time_discharged, exposure_days, FoI_community, vaccinated) %>% 
        #     mutate(time.diff = diff(c(0,time_discharged))/scaling.factor) %>% 
        #     sample_n(size = 2000, replace = TRUE) %>% 
        #     mutate(total.time = cumsum(time.diff)))
        # }
      
      workers.times <- workers.times %>% filter(total.time <= t.max)
      
      }
      
    } else{
      workers.times <- my.empty.df
      # return empty data.frame
    }
    
    # Should probably put a while in here to keep simulating until t.max, 
    # but let's gamble that 'size' is big enough that this won't be a problem (it's fast)
    # If considering high prevalence and high n.arrivals, perhaps increase
    # Don't need time_discharged here, but keeping as an identifier so we can look at 
    # entries if needed
    
    all.times <- arrivals.times %>% 
      add_row(workers.times) %>% 
      mutate(modified_FoI = FoI_community * (1 - vaccinated*VE_transmission))  # multiplying by vaccinated status will turn this off for unvaccinated individuals
    # don't bother arranging here, as will arrange for all output
  } else{
    all.times <- my.empty.df
  }
  return(all.times)
  
}


# Define structure of the empty dataframe for inclusion
my.empty.df <- data.frame(time_discharged = Inf, exposure_days = NA, FoI_community = NA, time.diff = Inf, total.time = Inf, vaccinated = NA)[numeric(0),]



