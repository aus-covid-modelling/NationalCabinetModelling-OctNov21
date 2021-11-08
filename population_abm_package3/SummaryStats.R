library(tidyverse)
library(ggplot2)
library(ggdist)
library(MASS)
gi <- read_csv("GI_1_5.csv")
tost <- read_csv("TOST_1_5.csv")
si <- read_csv("SI_1_5.csv")
gi <- read_csv("GI.csv")
tost <- read_csv("TOST.csv")
si <- read_csv("SI.csv")
# x <- seq(0,20,by=0.25)
# gi_vec <- dgamma(x,shape = 6.8,rate =1.23)
# ggplot(gi,aes(x = GI,y = ..density..)) +geom_density(fill="red", alpha = 0.2) + geom_line(data = tibble(x,gi_vec),aes(x = x,y =gi_vec),linetype="dashed") + scale_x_continuous("Generation interval")

discrete_gi <- read_csv("ttiq_gi_distributions.csv")
ggplot(discrete_gi,aes(x = days)) + geom_line(aes(y = none)) + geom_line(aes(y = partial),color = "red")+ geom_line(aes(y = optimal),color="blue")+ scale_x_continuous("Generation interval",limits =c(0,20))


ggplot(gi,aes(x = GI,y = ..density..)) +geom_histogram(binwidth = 1.0,center = 0.5,color = "black",fill = "red", alpha = 0.4) + geom_line(data = discrete_gi,aes(x = days,y =none),linetype="dashed") + scale_x_continuous("Generation interval",limits =c(0,20))

ggplot(gi,aes(x = floor(GI),y = ..density..)) +geom_histogram(binwidth = 1.0,boundary = 0,color = "black",fill = "red", alpha = 0.4) + geom_line(data = discrete_gi,aes(x = days ,y =none),linetype="dashed") + scale_x_continuous("Generation interval",limits =c(0,30))

ggplot(gi,aes(x = GI,y = ..density..)) +geom_density() + geom_line(data = discrete_gi,aes(x = days ,y =none),linetype="dashed") + scale_x_continuous("Generation interval",limits =c(0,20))

ggplot(gi,aes(x = floor(GI),y = ..density..)) +geom_histogram(binwidth =1,center = 0) + geom_line(data = discrete_gi,aes(x = days ,y =none),linetype="dashed") + scale_x_continuous("Generation interval",limits =c(0,20))



ggplot(si,aes(x = SI, y = ..density..)) + geom_histogram(binwidth = 1)
ML_SI <- fitdistr(si$SI,"negative binomial")
ML_SI


x2 <- seq(-10,10,by = 0.25)
tost_vec <- dstudent_t(x2,df = 3.3454,mu = -0.0747,sigma = 1.8567)
ggplot(tost,aes(x = TOST,y = ..density..)) +geom_density(kernel = "rectangular",fill="red", color = "red",alpha = 0.2) + geom_line(data = tibble(x2,tost_vec),aes(x = x2,y =tost_vec),linetype="dashed") + scale_x_continuous("Tost")

ggplot(tost,aes(x = TOST,y = ..density..)) +geom_histogram(binwidth = 1.0,color = "black",fill = "red",alpha = 0.4) + geom_line(data = tibble(x2,tost_vec),aes(x = x2,y =tost_vec),linetype="dashed") + scale_x_continuous("Tost")

length(which(tost$TOST< 0))/nrow(tost)



sim_number_1 <- read_csv("outputs/TestRuns/sim_number_1.csv")
ggplot(sim_number_1,aes(x = `Time of symptom onset`,fill = as.factor(Cluster))) + geom_histogram(binwidth = 1,alpha = 0.4,position = "identity")

+ facet_wrap(~Cluster)

ggplot(sim_number_1,aes(x = `Secondary Infections`)) + geom_histogram() 

+ facet_wrap(~Cluster)

fitdistr(sim_number_1$`Secondary Infections`,"negative binomial")

test <- mutate(gi,fG = floor(GI))
a <- test %>% group_by(fG) %>% summarise(prop = n()/nrow(gi))
ggplot(a,aes(x = fG,y = prop)) +geom_line() + geom_line(data = discrete_gi,aes(x = days ,y =none),linetype="dashed") + scale_x_continuous("Generation interval",limits =c(0,30))
