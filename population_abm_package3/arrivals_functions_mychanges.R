# Functions for producing linelist of arrival spillovers
# 
# David J Price
# 4 October, 2021
# 


## Function to create linelist for each pathway (and numbers through that pathway) and combine ----
create_linelist <- function(num.arrivals, arrivals.files, workers.files, pathway.prevalence, t.max = 180, base.n.arrivals = 100, VE_transmission= 0.6, baseline.prevalence = 0.01,family.pathway){
  lapply(X = 1:length(num.arrivals), 
         FUN = function(i){generate_one_linelist(arrivals = arrivals.files[i], workers = workers.files[i], 
                                                 quarantine.folder = quarantine.output.folder, n.arrivals = num.arrivals[i], pathway.prevalence = pathway.prevalence[i],
                                                 t.max = t.max, base.n.arrivals = base.n.arrivals, VE_transmission = VE_transmission, baseline.prevalence,family.pathway = family.pathway[i])}) %>%
  bind_rows() %>% arrange(total.time)
}


## This function creates a single linelist of spillovers from the corresponding arrival and worker files ----
generate_one_linelist <- function(arrivals, workers = NA, quarantine.folder, n.arrivals, pathway.prevalence, t.max, base.n.arrivals = 100, VE_transmission = 0.66, baseline.prevalence = 0.01,family.pathway){
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
      sample_n(size = 20000, replace = TRUE) %>% 
      mutate(total.time = cumsum(time.diff))
    
    while(max(arrivals.times$total.time,na.rm=TRUE)<t.max){
        arrivals.times <- arrivals.times %>% bind_rows(.,arrivals.dat %>% 
        select(time_discharged, exposure_days, FoI_community, vaccinated) %>% 
        mutate(time.diff = diff(c(0,time_discharged))/scaling.factor) %>% 
        sample_n(size = 20000, replace = TRUE))%>% 
          mutate(total.time = cumsum(time.diff))
    }
    
    arrivals.times <- arrivals.times %>% filter(total.time <= t.max) %>% 
      mutate(modified_FoI = FoI_community * (1 - vaccinated*VE_transmission))
    
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
      
        while(max(workers.times$total.time,na.rm=TRUE)<t.max){
            workers.times <- workers.times %>% bind_rows(.,workers.dat %>% 
            select(time_discharged, exposure_days, FoI_community, vaccinated) %>% 
            mutate(time.diff = diff(c(0,time_discharged))/scaling.factor) %>% 
            sample_n(size = 20000, replace = TRUE)) %>% 
              mutate(total.time = cumsum(time.diff))
        }
      
      workers.times <- workers.times %>% filter(total.time <= t.max) %>% 
        mutate(modified_FoI = FoI_community * (1 - vaccinated*VE_transmission))
      
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
      add_row(workers.times)  %>% mutate(family = family.pathway)
    if("child" %in% colnames(arrivals.dat)){
      all.times <- left_join(all.times,arrivals.dat%>% select(time_discharge,FoI_community,child),by = c("time_discharge","FoI_community")) %>% rename(kid = child)
    }else{
      
     all.times<- all.times %>% mutate(kid = as.numeric(family & !vaccinated)) 
      
    }
    
    
    all.times<- all.times %>% mutate(modified_FoI = case_when(kid==1~modified_FoI*0.6,TRUE~modified_FoI)) %>% select(time_discharged,exposure_days,FoI_community,time.diff,vaccinated,total.time,modified_FoI,kid)# multiplying by vaccinated status will turn this off for unvaccinated individuals
    # don't bother arranging here, as will arrange for all output


  } else{
    all.times <- my.empty.df
  }
  
  
  return(all.times)
  
}


# Define structure of the empty dataframe for inclusion
my.empty.df <- data.frame(time_discharged = Inf, exposure_days = NA, FoI_community = NA, time.diff = Inf, vaccinated = NA,total.time = Inf,modified_FoI = NA)[numeric(0),]



