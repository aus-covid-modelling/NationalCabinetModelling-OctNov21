ind_vax <- sim_number_1 %>% mutate(Bracket = cut(Age,breaks = c(0,12,16,seq(20,80,by = 10),Inf),right = FALSE)) %>% group_by(Bracket, `Current Vaccine`,`Current doses`) %>% summarise(n = n())

total_people <- ind_vax %>% ungroup() %>% group_by(Bracket) %>% summarise(total = sum(n))

test<-left_join(ind_vax,total_people) %>% mutate(prop=n/total)

ind_vax <- sim_number_1 %>% mutate(Bracket = cut(Age,breaks = c(0,12,16,Inf),right = FALSE)) %>% group_by(Bracket, `Current doses`) %>% summarise(n = n())

total_people <- ind_vax %>% ungroup() %>% group_by(Bracket) %>% summarise(total = sum(n))
test<-left_join(ind_vax,total_people) %>% mutate(prop=n/total)
