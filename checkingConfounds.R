library("lme4")

setwd("~/Documents/Scripts/GroupMollyTracking/")

indv <- read.csv("JolleTracksAll_2.csv")
indv$ExpDay <- as.integer(indv$ExpDay)

indv$week <- indv$ExpDay %/% 7 ## Is this integer division
indv$triday <- indv$ExpDay %/% 3

hourly <- read.csv("JolleTracksHourly_2.csv")
hourly$Hour <- as.integer(hourly$Hour)
hourly$ExpDay <- as.integer(hourly$ExpDay)

model.simple <- lm(dist_mean ~ pDist_mean,data=indv)
summary(model.simple)

model.re <- lmer(dist_mean ~ pDist_mean + (1|Pi),data=indv)
summary(model.re)
