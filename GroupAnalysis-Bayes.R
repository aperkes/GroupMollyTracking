
library(nlme)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MCMCglmm)
library(parallel)
library(reshape2)
library(corrplot)
library(MuMIn)
library(gridExtra)
library(ggcorrplot)
library(zoo)
library(egg)
library(patchwork)
library(viridis)
library(formatR)

setwd("~/Documents/Scripts/GroupMollyTracking/")

indv <- read.csv("JolleTracksAll_4.csv")
indv$ExpDay <- as.integer(indv$ExpDay)

indv[indv$Pi == 'pi11',]$ExpDay <- indv[indv$Pi == 'pi11',]$ExpDay - 16
indv[indv$Pi == 'pi12',]$ExpDay <- indv[indv$Pi == 'pi12',]$ExpDay + 15

indv$week <- indv$ExpDay %/% 7 ## Is this integer division
indv$triday <- indv$ExpDay %/% 3

indv.good <- indv[indv$ExpDay >= 0,]
#indv.good <- indv.good[indv.good$Pi != 'pi12',]

indv <- indv.good

hourly <- read.csv("JolleTracksHourly_2.csv")
hourly$Hour <- as.integer(hourly$Hour)
hourly$ExpDay <- as.integer(hourly$ExpDay)

## Some of the pi's fall off over time (e.g., pi41 drops off at day 11)
## Others go on 90 days
piDays <- indv %>%
  group_by(Pi) %>%
  summarize(max = max(ExpDay))

piDays[piDays$max >= 54,]
long_data <- indv %>%
  filter(Pi %in% piDays[piDays$max >= 62,]$Pi)

indv.long <- long_data[long_data$ExpDay <= 63,]

long_data54 <- indv %>%
  filter(Pi %in% piDays[piDays$max >= 54,]$Pi)
indv.long54 <- long_data54[long_data54$ExpDay <= 54,]

## Only look at 4 weeks for fish where we have that
good_data <- indv[indv$Pi != "pi34",]
good_data <- good_data[good_data$Pi != "pi41",]
good_data <- good_data[good_data$ExpDay <= 31,]

indv.com <- good_data

### Same for hourly
good_hourly <- hourly[hourly$Pi != "pi34",]
good_hourly <- good_hourly[good_hourly$Pi != "pi41",]
good_hourly <- good_hourly[good_hourly$ExpDay <= 31,]

hourly.com <- good_hourly

day1 <- hourly %>%
  filter(ExpDay == 0)
day1.com <- day1 %>%
  filter(Pi %in% good_hourly$Pi)

day27 <- hourly %>%
  filter(ExpDay == 27)
day27.com <- day27 %>%
  filter(Pi %in% good_hourly$Pi)


### Prior Specifications (buckle up...copied from Kate's code)

prior.null <- list(R = list(V = 1, nu = 0.002))

prior.id <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1=list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 25^2)))

prior.id.slope <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1=list(V = diag(2), nu = 0.002, alpha.mu = c(0,0),  alpha.V = diag(2)*25^2)))

prior.cov4 <- list(R = list(V = diag(4), nu = 4.002),
                   G = list(G1=list(V = diag(4), nu = 4.002, alpha.mu = rep(0,4), alpha.V = 1000*diag(4))))

prior.cov8 <- list(R = list(V = diag(8), nu = 8.002),
                   G = list(G1=list(V = diag(8), nu = 8.002, alpha.mu = rep(0,8), alpha.V = 1000*diag(8))))

prior.cov9 <- list(R = list(V = diag(9), nu = 9.002),
                   G = list(G1=list(V = diag(9), nu = 9.002, alpha.mu = rep(0,9), alpha.V = 1000*diag(9))))


prior.cov10 <- list(R = list(V = diag(10), nu = 10.002),
                    G = list(G1=list(V = diag(10), nu = 10.002, alpha.mu = rep(0,10), alpha.V = 1000*diag(10))))

### DATA EXPLORATION GOES HERE

## Daily 
pairs(indv.com[,c("vel_mean", "pDist_mean", "dist_mean", "NearN_mean", "pDistC_mean", "velC_mean", "head_mean", "angleC_mean")])

## Hourly 
pairs(hourly.com[,c("vel_mean_", "pDist_mean_", "dist_mean_", "NearN_mean_", "pDistC_mean_", "velC_mean_", "head_mean_", "angleC_mean_")])

### PCA

## Important to scale since data ranges are dramatically different
indv.scaled <- scale(indv.com[,c("vel_mean", "pDist_mean", "dist_mean", "pDistC_mean", "velC_mean", "head_mean", "angleC_mean")],center = TRUE, scale = TRUE)
indv.scaled$Pi <- indv.com$Pi
indv.scaled$ExpDay <- indv.com$ExpDay

## So nothing too dramatic here, other than that vel_mean weights separately 
##  from all the 'cohesion' measures, which is sort of neat. 
pc <- prcomp(indv.scaled[,c("vel_mean", "pDist_mean", "dist_mean", "pDistC_mean", "velC_mean", "head_mean", "angleC_mean")])

pc
## What variable should we use? 


### Looking for Group differences

# pDistC_mean_ is the mean (hourly) correlation in distance from center

lm.test <- lme(pDistC_mean_ ~ Hour, random = ~1|Pi, data = day1.com)

#par(mfrow=c(2,2))
plot(lm.test)
hist(resid(lm.test))
qqnorm(resid(lm.test))
qqline(resid(lm.test))
boxplot(resid(lm.test) ~day1.com$Hour)

## Define some models: 

# Null model -----------
set.seed(58)
dist.day1.0 <-  MCMCglmm(dist_mean_ ~ Hour, 
                          data = day1.com, 
                          family = "gaussian", 
                          prior = prior.null, 
                          nitt=510000, burnin = 10000, thin = 200, 
                          verbose = F)


# Intercepts ID ----------------------
set.seed(3432)
dist.day1.1 <- MCMCglmm(dist_mean_ ~ Hour, 
                         random = ~Pi, 
                         data = day1.com, 
                         family = "gaussian", 
                         prior = prior.id, 
                         nitt=510000, burnin = 10000, thin = 200, 
                         verbose = F)



# Intercepts and Slopes ID -----------------------
set.seed(472)
dist.day1.2 <- MCMCglmm(dist_mean_ ~ Hour, 
                         random = ~us(1 + Hour):Pi, 
                         data = day1.com, 
                         family = "gaussian", 
                         prior = prior.id.slope, 
                         nitt=510000, burnin = 10000, thin = 200, 
                         verbose = F)


DIC(dist.day1.0, dist.day1.1, dist.day1.2)

## Big jump with group ID, only Slight improvement from adding Pi as a random slope, but not dramatic. 
## I guess this means it's important that hour & day are centered, but that already happened. 

summary(dist.day1.2)
posterior.mode(dist.day1.2$Sol)
HPDinterval(dist.day1.2$Sol)

posterior.mode(dist.day1.2$VCV)
HPDinterval(dist.day1.2$VCV)

# only ID repeatability 
rpt.spd1 <- dist.day1.2$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day1.2$VCV[,"(Intercept):(Intercept).Pi"]  + 
                                                                   dist.day1.2$VCV[,"Hour:Hour.Pi"] +
                                                                   dist.day1.2$VCV[,"units"])
posterior.mode(rpt.spd1)
HPDinterval(rpt.spd1)


# to get marginal R2 to explain fixed effects variance
vmVarF<-numeric(2500)

for(i in 1:2500){
  Var<-var(as.vector(dist.day1.2$Sol[i,] %*% t(dist.day1.2$X)))
  vmVarF[i]<-Var}


R2m<-vmVarF/(vmVarF+dist.day1.2$VCV[,1]+ dist.day1.2$VCV[,4] + dist.day1.2$VCV[,5])

posterior.mode(R2m)
HPDinterval(R2m)

R2c<-(vmVarF + dist.day1.2$VCV[,1])/(vmVarF + dist.day1.2$VCV[,1] + dist.day1.2$VCV[,4] + dist.day1.2$VCV[,5])

posterior.mode(R2c)
HPDinterval(R2c)

day1.com %>% ggplot(aes(Hour,dist_mean_,group_by='Pi',color=Pi)) + 
  geom_line() + 
  scale_color_viridis(discrete=TRUE) + 
  theme_classic() + 
  theme(legend.position = "none",
        rect=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 0)) 


### Check last day: 
set.seed(472)
dist.day27.2 <- MCMCglmm(dist_mean_ ~ Hour, 
                        random = ~us(1 + Hour):Pi, 
                        data = day27.com, 
                        family = "gaussian", 
                        prior = prior.id.slope, 
                        nitt=510000, burnin = 10000, thin = 200, 
                        verbose = F)



## Big jump with group ID, only Slight improvement from adding Pi as a random slope, but not dramatic. 
## I guess this means it's important that hour & day are centered, but that already happened. 

summary(dist.day27.2)
posterior.mode(dist.day27.2$Sol)
HPDinterval(dist.day27.2$Sol)

posterior.mode(dist.day27.2$VCV)
HPDinterval(dist.day27.2$VCV)

# only ID repeatability 
rpt.spd1 <- dist.day27.2$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day27.2$VCV[,"(Intercept):(Intercept).Pi"]  + 
                                                              dist.day27.2$VCV[,"Hour:Hour.Pi"] +
                                                              dist.day27.2$VCV[,"units"])
posterior.mode(rpt.spd1)
HPDinterval(rpt.spd1)


# to get marginal R2 to explain fixed effects variance
vmVarF<-numeric(2500)

for(i in 1:2500){
  Var<-var(as.vector(dist.day27.2$Sol[i,] %*% t(dist.day27.2$X)))
  vmVarF[i]<-Var}


R2m<-vmVarF/(vmVarF+dist.day27.2$VCV[,1]+ dist.day27.2$VCV[,4] + dist.day27.2$VCV[,5])

posterior.mode(R2m)
HPDinterval(R2m)

R2c<-(vmVarF + dist.day27.2$VCV[,1])/(vmVarF + dist.day27.2$VCV[,1] + dist.day27.2$VCV[,4] + dist.day27.2$VCV[,5])

posterior.mode(R2c)
HPDinterval(R2c)

day27.com %>% ggplot(aes(Hour,dist_mean_,group_by='Pi',color=Pi)) + 
  geom_line() + 
  scale_color_viridis(discrete=TRUE) + 
  theme_classic() + 
  theme(legend.position = "none",
        rect=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 0)) 


### ADD VARIOUS SENSITIVITY ANALYSES HERE

### OTHER BEHAVIORS

### Week 1 days

## Day 2

day2 <- hourly.com %>%
  filter(ExpDay == 1)


# Null model -----------
set.seed(58)
dist.day2.0 <-  MCMCglmm(pDistC_mean_ ~ Hour, 
                         data = day2, 
                         family = "gaussian", 
                         prior = prior.null, 
                         nitt=510000, burnin = 10000, thin = 200, 
                         verbose = F)


# Intercepts ID ----------------------
set.seed(3432)
dist.day2.1 <- MCMCglmm(pDistC_mean_ ~ Hour, 
                        random = ~Pi, 
                        data = day2, 
                        family = "gaussian", 
                        prior = prior.id, 
                        nitt=510000, burnin = 10000, thin = 200, 
                        verbose = F)



# Intercepts and Slopes ID -----------------------
set.seed(472)
dist.day2.2 <- MCMCglmm(pDistC_mean_ ~ Hour, 
                        random = ~us(1 + Hour):Pi, 
                        data = day2, 
                        family = "gaussian", 
                        prior = prior.id.slope, 
                        nitt=510000, burnin = 10000, thin = 200, 
                        verbose = F)


DIC(dist.day2.0, dist.day2.1, dist.day2.2)


## Big jump with group ID, here the slopes don't help, but I'm including them for consistency. 
## I guess this means it's important that hour & day are centered, but that already happened. 

summary(dist.day2.2)
posterior.mode(dist.day2.2$Sol)
HPDinterval(dist.day2.2$Sol)

posterior.mode(dist.day2.2$VCV)
HPDinterval(dist.day2.2$VCV)

# only ID repeatability 
rpt.spd1 <- dist.day2.2$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day2.2$VCV[,"(Intercept):(Intercept).Pi"]  + 
                                                              dist.day2.2$VCV[,"Hour:Hour.Pi"] +
                                                              dist.day2.2$VCV[,"units"])
posterior.mode(rpt.spd1)
HPDinterval(rpt.spd1)


# to get marginal R2 to explain fixed effects variance
vmVarF<-numeric(2500)

for(i in 1:2500){
  Var<-var(as.vector(dist.day2.2$Sol[i,] %*% t(dist.day2.2$X)))
  vmVarF[i]<-Var}


R2m<-vmVarF/(vmVarF+dist.day2.2$VCV[,1]+ dist.day2.2$VCV[,4] + dist.day2.2$VCV[,5])

posterior.mode(R2m)
HPDinterval(R2m)

R2c<-(vmVarF + dist.day2.2$VCV[,1])/(vmVarF + dist.day2.2$VCV[,1] + dist.day2.2$VCV[,4] + dist.day2.2$VCV[,5])

posterior.mode(R2c)
HPDinterval(R2c)

## ADD : Do this for all the days 1-7

### ENTIRE WEEK 1

set.seed(58)
# Null model ----------------------
dist.week1.0 <-  MCMCglmm(pDistC_mean_ ~ Hour, 
                           data = hourly.com[which(hourly.com$ExpDay < 7),], 
                           family = "gaussian", 
                           prior = prior.null, 
                           nitt=510000, burnin = 10000, thin = 200, 
                           verbose = F)

# Intercepts ID ----------------------
dist.week1.1 <-  MCMCglmm(pDistC_mean_ ~ ExpDay, 
                          random = ~Pi,
                          data = hourly.com[which(hourly.com$ExpDay < 7),], 
                          family = "gaussian", 
                          prior = prior.id, 
                          nitt=510000, burnin = 10000, thin = 200, 
                          verbose = F)

# Intercepts and slopes ID ----------------------
dist.week1.2 <-  MCMCglmm(pDistC_mean_ ~ ExpDay, 
                          random = ~us(1 + ExpDay):Pi, 
                          data = hourly.com[which(hourly.com$ExpDay < 7),], 
                          family = "gaussian", 
                          prior = prior.id.slope, 
                          nitt=510000, burnin = 10000, thin = 200, 
                          verbose = F)


DIC(dist.week1.0, dist.week1.1, dist.week1.2)

## MOdel to report: 

summary(dist.week1.2)
posterior.mode(dist.week1.2$Sol)
HPDinterval(dist.week1.2$Sol)

posterior.mode(dist.week1.2$VCV)
HPDinterval(dist.week1.2$VCV)

# only ID repeatability 
rpt.spd1 <- dist.week1.2$VCV[,"(Intercept):(Intercept).Pi"]/(dist.week1.2$VCV[,"(Intercept):(Intercept).Pi"]  + 
                                                              dist.week1.2$VCV[,"ExpDay:ExpDay.Pi"] +
                                                              dist.week1.2$VCV[,"units"])
posterior.mode(rpt.spd1)
HPDinterval(rpt.spd1)


# to get marginal R2 to explain fixed effects variance
vmVarF<-numeric(2500)

for(i in 1:2500){
  Var<-var(as.vector(dist.week1.2$Sol[i,] %*% t(dist.week1.2$X)))
  vmVarF[i]<-Var}


R2m<-vmVarF/(vmVarF+dist.week1.2$VCV[,1]+ dist.week1.2$VCV[,4] + dist.week1.2$VCV[,5])

posterior.mode(R2m)
HPDinterval(R2m)

R2c<-(vmVarF + dist.week1.2$VCV[,1])/(vmVarF + dist.week1.2$VCV[,1] + dist.week1.2$VCV[,4] + dist.week1.2$VCV[,5])

posterior.mode(R2c)
HPDinterval(R2c)

### ADD OTHER WEEK1 BEHAVIORS

### Don't have any body size I don't think, so that skips TONS of code in kate's notebook
### 1920 -> 2500

### On to Part 2, change over time

## First, model comparison stuff, as before

dist.0 <-  MCMCglmm(pDistC_mean ~ ExpDay, 
                          data = indv.com, 
                          family = "gaussian", 
                          prior = prior.null, 
                          nitt=510000, burnin = 10000, thin = 200, 
                          verbose = F)

# Intercepts ID ----------------------
dist.1 <-  MCMCglmm(pDistC_mean ~ ExpDay, 
                          random = ~Pi,
                          data = indv.com, 
                          family = "gaussian", 
                          prior = prior.id, 
                          nitt=510000, burnin = 10000, thin = 200, 
                          verbose = F)

# Intercepts and slopes ID ----------------------
dist.2 <-  MCMCglmm(pDistC_mean ~ ExpDay, 
                    random = ~us(1 + ExpDay):Pi,
                    data = indv.com, 
                    family = "gaussian", 
                    prior = prior.id.slope, 
                    nitt=510000, burnin = 10000, thin = 200, 
                    verbose = F)

DIC(dist.0, dist.1, dist.2)

### Report Model
summary(dist.2)

### ADD: More model Validation 

### Ok, on to the slicing

# select only those IDs in that vector & only keep up to obs 70 & make obs 1 = 0
indv.30days.cen <- indv.com %>%
  mutate(Day0 = ExpDay, 
         Day3 = ExpDay-3,
         Day6 = ExpDay-6,
         Day9 = ExpDay-9,
         Day12 = ExpDay-12,
         Day15 = ExpDay-15,
         Day18 = ExpDay-18,
         Day21 = ExpDay-21,
         Day24 = ExpDay-24,
         Day27 = ExpDay-27,
         Day30 = ExpDay-30)

set.seed(403)
dist.day0 <- MCMCglmm(pDistC_mean ~ Day0, 
                       random = ~us(1 + Day0):Pi, 
                       data = indv.30days.cen, 
                       family = "gaussian",
                       pr = T,
                       prior = prior.id.slope, 
                       nitt=310000, burnin = 10000, thin = 200, 
                       verbose = F)

set.seed(765)
dist.day3 <- MCMCglmm(pDistC_mean ~ Day3, 
                      random = ~us(1 + Day3):Pi, 
                      data = indv.30days.cen, 
                      family = "gaussian",
                      pr = T,
                      prior = prior.id.slope, 
                      nitt=310000, burnin = 10000, thin = 200, 
                      verbose = F)

dist.day6 <- MCMCglmm(pDistC_mean ~ Day6, 
                      random = ~us(1 + Day6):Pi, 
                      data = indv.30days.cen, 
                      family = "gaussian",
                      pr = T,
                      prior = prior.id.slope, 
                      nitt=310000, burnin = 10000, thin = 200, 
                      verbose = F)

dist.day9 <- MCMCglmm(pDistC_mean ~ Day9, 
                      random = ~us(1 + Day9):Pi, 
                      data = indv.30days.cen, 
                      family = "gaussian",
                      pr = T,
                      prior = prior.id.slope, 
                      nitt=310000, burnin = 10000, thin = 200, 
                      verbose = F)

dist.day12 <- MCMCglmm(pDistC_mean ~ Day12, 
                      random = ~us(1 + Day12):Pi, 
                      data = indv.30days.cen, 
                      family = "gaussian",
                      pr = T,
                      prior = prior.id.slope, 
                      nitt=310000, burnin = 10000, thin = 200, 
                      verbose = F)

dist.day15 <- MCMCglmm(pDistC_mean ~ Day15, 
                       random = ~us(1 + Day15):Pi, 
                       data = indv.30days.cen, 
                       family = "gaussian",
                       pr = T,
                       prior = prior.id.slope, 
                       nitt=310000, burnin = 10000, thin = 200, 
                       verbose = F)

dist.day18 <- MCMCglmm(pDistC_mean ~ Day18, 
                       random = ~us(1 + Day18):Pi, 
                       data = indv.30days.cen, 
                       family = "gaussian",
                       pr = T,
                       prior = prior.id.slope, 
                       nitt=310000, burnin = 10000, thin = 200, 
                       verbose = F)

dist.day21 <- MCMCglmm(pDistC_mean ~ Day21, 
                       random = ~us(1 + Day21):Pi, 
                       data = indv.30days.cen, 
                       family = "gaussian",
                       pr = T,
                       prior = prior.id.slope, 
                       nitt=310000, burnin = 10000, thin = 200, 
                       verbose = F)

dist.day24 <- MCMCglmm(pDistC_mean ~ Day24, 
                       random = ~us(1 + Day24):Pi, 
                       data = indv.30days.cen, 
                       family = "gaussian",
                       pr = T,
                       prior = prior.id.slope, 
                       nitt=310000, burnin = 10000, thin = 200, 
                       verbose = F)

dist.day27 <- MCMCglmm(pDistC_mean ~ Day27, 
                       random = ~us(1 + Day27):Pi, 
                       data = indv.30days.cen, 
                       family = "gaussian",
                       pr = T,
                       prior = prior.id.slope, 
                       nitt=310000, burnin = 10000, thin = 200, 
                       verbose = F)

dist.day30 <- MCMCglmm(pDistC_mean ~ Day30, 
                       random = ~us(1 + Day30):Pi, 
                       data = indv.30days.cen, 
                       family = "gaussian",
                       pr = T,
                       prior = prior.id.slope, 
                       nitt=310000, burnin = 10000, thin = 200, 
                       verbose = F)

plot(dist.day0$VCV)

## Copy pasting so much, getting day: variance
date <- c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30) ## "It's a magic number"

rpt0 <- dist.day0$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day0$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(dist.day0$VCV[,"units"]))
rpt3 <- dist.day3$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day3$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(dist.day3$VCV[,"units"]))
rpt6 <- dist.day6$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day6$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(dist.day6$VCV[,"units"]))
rpt9 <- dist.day9$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day9$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(dist.day9$VCV[,"units"]))
rpt12 <- dist.day12$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day12$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(dist.day12$VCV[,"units"]))
rpt15 <- dist.day15$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day15$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(dist.day15$VCV[,"units"]))
rpt18 <- dist.day18$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day18$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(dist.day18$VCV[,"units"]))
rpt21 <- dist.day21$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day21$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(dist.day21$VCV[,"units"]))
rpt24 <- dist.day24$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day24$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(dist.day24$VCV[,"units"]))
rpt27 <- dist.day27$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day27$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(dist.day27$VCV[,"units"]))
rpt30 <- dist.day30$VCV[,"(Intercept):(Intercept).Pi"]/(dist.day30$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(dist.day30$VCV[,"units"]))

rpt <- c(posterior.mode(rpt0), 
         posterior.mode(rpt3),
         posterior.mode(rpt6),
         posterior.mode(rpt9),
         posterior.mode(rpt12),
         posterior.mode(rpt15),
         posterior.mode(rpt18),
         posterior.mode(rpt21),
         posterior.mode(rpt24),
         posterior.mode(rpt27),
         posterior.mode(rpt30))


ci.rpt <- c(HPDinterval(rpt0)[1:2],
            HPDinterval(rpt3)[1:2],
            HPDinterval(rpt6)[1:2],
            HPDinterval(rpt9)[1:2],
            HPDinterval(rpt12)[1:2],
            HPDinterval(rpt15)[1:2],
            HPDinterval(rpt18)[1:2],
            HPDinterval(rpt21)[1:2],
            HPDinterval(rpt24)[1:2],
            HPDinterval(rpt27)[1:2],
            HPDinterval(rpt30)[1:2])
ci.rpt <- matrix(ci.rpt, nrow = 11, byrow = T)

post.id <- c(posterior.mode(dist.day0$VCV[,"(Intercept):(Intercept).Pi"]), 
             posterior.mode(dist.day3$VCV[,"(Intercept):(Intercept).Pi"]),
             posterior.mode(dist.day6$VCV[,"(Intercept):(Intercept).Pi"]),
             posterior.mode(dist.day9$VCV[,"(Intercept):(Intercept).Pi"]),
             posterior.mode(dist.day12$VCV[,"(Intercept):(Intercept).Pi"]),
             posterior.mode(dist.day15$VCV[,"(Intercept):(Intercept).Pi"]),
             posterior.mode(dist.day18$VCV[,"(Intercept):(Intercept).Pi"]),
             posterior.mode(dist.day21$VCV[,"(Intercept):(Intercept).Pi"]),
             posterior.mode(dist.day24$VCV[,"(Intercept):(Intercept).Pi"]),
             posterior.mode(dist.day27$VCV[,"(Intercept):(Intercept).Pi"]),
             posterior.mode(dist.day30$VCV[,"(Intercept):(Intercept).Pi"]))

ci.id <- c(HPDinterval(dist.day0$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
           HPDinterval(dist.day3$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
           HPDinterval(dist.day6$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
           HPDinterval(dist.day9$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
           HPDinterval(dist.day12$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
           HPDinterval(dist.day15$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
           HPDinterval(dist.day18$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
           HPDinterval(dist.day21$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
           HPDinterval(dist.day24$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
           HPDinterval(dist.day27$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
           HPDinterval(dist.day30$VCV[,"(Intercept):(Intercept).Pi"])[1:2])

ci.id <-matrix(ci.id, nrow = 11, byrow = T)

post.w <- c(posterior.mode(dist.day0$VCV[,"units"]),
            posterior.mode(dist.day3$VCV[,"units"]),
            posterior.mode(dist.day6$VCV[,"units"]),
            posterior.mode(dist.day9$VCV[,"units"]),
            posterior.mode(dist.day12$VCV[,"units"]),
            posterior.mode(dist.day15$VCV[,"units"]),
            posterior.mode(dist.day18$VCV[,"units"]),
            posterior.mode(dist.day21$VCV[,"units"]),
            posterior.mode(dist.day24$VCV[,"units"]),
            posterior.mode(dist.day27$VCV[,"units"]),
            posterior.mode(dist.day30$VCV[,"units"]))

ci.w <- c(HPDinterval(dist.day0$VCV[,"units"])[1:2],   
          HPDinterval(dist.day3$VCV[,"units"])[1:2], 
          HPDinterval(dist.day6$VCV[,"units"])[1:2], 
          HPDinterval(dist.day9$VCV[,"units"])[1:2], 
          HPDinterval(dist.day12$VCV[,"units"])[1:2], 
          HPDinterval(dist.day15$VCV[,"units"])[1:2], 
          HPDinterval(dist.day18$VCV[,"units"])[1:2], 
          HPDinterval(dist.day21$VCV[,"units"])[1:2], 
          HPDinterval(dist.day24$VCV[,"units"])[1:2], 
          HPDinterval(dist.day27$VCV[,"units"])[1:2], 
          HPDinterval(dist.day30$VCV[,"units"])[1:2])

ci.w <- matrix(ci.w, nrow = 11, byrow = T)            
plot(post.id ~ date)
plot(post.w ~ date)
plot(rpt ~ date)

rpt.slice.wide <- data.frame(date, rpt, "lower.rpt" = ci.rpt[,1], "upper.rpt" = ci.rpt[,2],
                             post.id, "lower.id" = ci.id[,1], "upper.id" = ci.id[,2], 
                             post.w, "lower.w" = ci.w[,1], "upper.w" = ci.w[,2])

rpt.slice.long <- data.frame("date" = rep(date, 3), 
                             "type" = rep(c("rpt", "id", "within"), each = 11),
                             "variance" = unname(c(rpt, post.id, post.w)),
                             "lower" = unname(c(ci.rpt[,1], ci.id[,1], ci.w[,1])),
                             "upper" = unname(c(ci.rpt[,2], ci.id[,2], ci.w[,2])))


### This is very dataset specific, you may need to check it: 
ids <- colnames(dist.day0$Sol)[3:14] ## Make sure this matches below
ids <- substr(ids, 16, 19)

'this pulls out the individual intercepts and adds in the overall intercepts

so that way these numbers are absolute values, as opposed to differences from overall'

intercepts0 <- unname(posterior.mode(dist.day0$Sol)[3:14] + posterior.mode(dist.day0$Sol)["(Intercept)"])
intercepts3 <- unname(posterior.mode(dist.day3$Sol)[3:14] + posterior.mode(dist.day0$Sol)["(Intercept)"])
intercepts6 <- unname(posterior.mode(dist.day6$Sol)[3:14] + posterior.mode(dist.day6$Sol)["(Intercept)"])
intercepts9 <- unname(posterior.mode(dist.day9$Sol)[3:14] + posterior.mode(dist.day9$Sol)["(Intercept)"])
intercepts12 <- unname(posterior.mode(dist.day12$Sol)[3:14] + posterior.mode(dist.day12$Sol)["(Intercept)"])
intercepts15 <- unname(posterior.mode(dist.day15$Sol)[3:14] + posterior.mode(dist.day15$Sol)["(Intercept)"])
intercepts18 <- unname(posterior.mode(dist.day18$Sol)[3:14] + posterior.mode(dist.day18$Sol)["(Intercept)"])
intercepts21 <- unname(posterior.mode(dist.day21$Sol)[3:14] + posterior.mode(dist.day21$Sol)["(Intercept)"])
intercepts24 <- unname(posterior.mode(dist.day24$Sol)[3:14] + posterior.mode(dist.day24$Sol)["(Intercept)"])
intercepts27 <- unname(posterior.mode(dist.day27$Sol)[3:14] + posterior.mode(dist.day27$Sol)["(Intercept)"])
intercepts30 <- unname(posterior.mode(dist.day30$Sol)[3:14] + posterior.mode(dist.day30$Sol)["(Intercept)"])

date <- rep(c(0,3,6,9,12,15,18,21,24,27,30), each = 12)
picomp <- rep(ids, 11)
blup <- c(intercepts0, intercepts3, intercepts6, intercepts9, intercepts12, intercepts15, 
          intercepts18, intercepts21, intercepts24, intercepts27, intercepts30)

pred.intercepts <- data.frame(date, picomp, blup)

ggplot(pred.intercepts, aes(x=date, y=blup, group=picomp, color=picomp)) + 
  geom_line()
## START AT LINE 3065 in the other sheet.

indv.com$week <- indv.com$ExpDay %/% 7 ## Is this integer division
indv.com$triday <- indv.com$ExpDay %/% 3
indv.com$fiver <- indv.com$ExpDay %/% 5

### This does it by three day blocks, but that makes the sample size really low
indv.wide <- indv.com[indv.com$triday < 10,] %>%
  mutate(
    obs.day = case_when(
      ExpDay %in% seq(0,31,3) ~ 0, 
      ExpDay %in% seq(1,31,3) ~ 1,
      ExpDay %in% seq(2,31,3) ~ 2,
      ExpDay %in% seq(3,31,3) ~ 3,
      ExpDay %in% seq(4,31,3) ~ 4,
      ExpDay %in% seq(5,31,3) ~ 5,
      ExpDay %in% seq(6,31,3) ~ 6,
      ExpDay %in% seq(7,31,3) ~ 7
    )) %>%
  select(Pi, triday, pDistC_mean, obs.day) %>%
  spread(triday, pDistC_mean, sep = "")

behav.triday.id <- MCMCglmm(cbind(triday0, triday1, triday2, triday3, triday4, triday5, triday6, triday7, triday8, triday9) ~ trait - 1, 
                          random = ~us(trait):Pi,
                          rcov = ~us(trait):units,
                          family = c(rep("gaussian", 10)), 
                          prior = prior.cov10, 
                          pr = T, 
                          nitt = 510000, thin = 200, burnin = 10000, 
                          verbose = F,
                          data = indv.wide)

# Model with only ID ----
id.matrix.triday <- matrix(posterior.mode(posterior.cor(behav.triday.id$VCV[,1:100])),10,10, 
                         dimnames = list(c("Triday 1", "Triday 2", "Triday 3", "Triday 4", "Triday 5", "Triday 6", "Triday 7", "Triday 8", "Triday 9", "Triday 10"), 
                                         c("Triday 1", "Triday 2", "Triday 3", "Triday 4", "Triday 5", "Triday 6", "Triday 7", "Triday 8", "Triday 9", "Triday 10")))

# now to extract the CI estimates
ci.triday <- data.frame(HPDinterval(posterior.cor(behav.triday.id$VCV[,1:100])))

# for corrplot need 3 matrices - estimates, lower CI, upper CI
lower.triday <- matrix(ci.triday[,1],10,10)
upper.triday <- matrix(ci.triday[,2],10,10)

test <- melt(lower.triday) %>%
  mutate(p.value = ifelse(value < 0, 1, 0)) %>%
  select(Var1, Var2, p.value)

p.mat <- diag(10)
p.mat[cbind(test$Var1, test$Var2)] <- p.mat[cbind(test$Var2, test$Var1)] <- test$p.value


triday.corr <- ggcorrplot(id.matrix.triday, 
                         type = "lower", 
                         p.mat = p.mat, 
                         insig= "blank",
                         colors = c("slateblue4","gray", "mediumorchid1"))

triday.corr

## Trying by weeks, it's a tradeoff, but this seems better. 

## Doing longer days reduces our pi sample size down to just a few. This is tuff. 
#indv.wide <- indv.com[indv.com$week < 4,] %>%

indv.wide <- indv.com[indv.com$week < 4,] %>%
  mutate(
    obs.day = case_when(
      ExpDay %in% seq(0,27,7) ~ 1, 
      ExpDay %in% seq(1,27,7) ~ 2,
      ExpDay %in% seq(2,27,7) ~ 3,
      ExpDay %in% seq(3,27,7) ~ 4,
      ExpDay %in% seq(4,27,7) ~ 5,
      ExpDay %in% seq(5,27,7) ~ 6,
      ExpDay %in% seq(6,27,7) ~ 7
    )) %>%
  select(Pi, week, dist_mean, obs.day) %>%
  spread(week, dist_mean, sep = "")

behav.week.id <- MCMCglmm(cbind(week0, week1, week2, week3) ~ trait - 1, 
                          random = ~us(trait):Pi,
                          rcov = ~us(trait):units,
                          family = c(rep("gaussian", 4)), 
                          prior = prior.cov4, 
                          pr = T, 
                          nitt = 510000, thin = 200, burnin = 10000, 
                          verbose = F,
                          data = indv.wide)

# Model with only ID ----
id.matrix.week <- matrix(posterior.mode(posterior.cor(behav.week.id$VCV[,1:16])),4,4, 
                         dimnames = list(c("week 1", "week 2", "week 3", "week 4"), 
                                         c("week 1", "week 2", "week 3", "week 4")))

# now to extract the CI estimates
ci.week <- data.frame(HPDinterval(posterior.cor(behav.week.id$VCV[,1:16])))

# for corrplot need 3 matrices - estimates, lower CI, upper CI
lower.week <- matrix(ci.week[,1],4,4)
upper.week <- matrix(ci.week[,2],4,4)

test <- melt(lower.week) %>%
  mutate(p.value = ifelse(value < 0, 1, 0)) %>%
  select(Var1, Var2, p.value)

p.mat <- diag(4)
p.mat[cbind(test$Var1, test$Var2)] <- p.mat[cbind(test$Var2, test$Var1)] <- test$p.value


week.corr <- ggcorrplot(id.matrix.week, 
                        type = "lower", 
                        p.mat = p.mat, 
                        insig= "blank",
                        colors = c("slateblue4","gray", "mediumorchid1"))

week.corr


ggplot(indv.wide,aes(week0,week3)) + geom_point() + geom_smooth(method='lm')
lm.model.weeks <- lm(week3 ~ week0, indv.wide)
summary(lm.model.weeks)

indv.wide.day1 <- indv.wide[indv.wide$obs.day == 1,]
lm.model.simple <- lm(week3 ~ week0, indv.wide.day1)
summary(lm.model.simple)

ggplot(indv.wide.day1,aes(week0,week3)) + 
  geom_point() + geom_smooth(method='lm') + 
  theme_classic() + 
  theme(legend.position = "none",
        rect=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 0)) 

### This uses the 8 pi's that run 9 weeks, but sample size is too small for it to work. 
if(TRUE){ 
  indv.wide <- indv.long %>%
    mutate(
      obs.day = case_when(
        ExpDay %in% seq(0,62,7) ~ 1, 
        ExpDay %in% seq(1,62,7) ~ 2,
        ExpDay %in% seq(2,62,7) ~ 3,
        ExpDay %in% seq(3,62,7) ~ 4,
        ExpDay %in% seq(4,62,7) ~ 5,
        ExpDay %in% seq(5,62,7) ~ 6,
        ExpDay %in% seq(6,62,7) ~ 7
        )) %>%
    select(Pi, week, dist_mean, obs.day) %>%
    spread(week, dist_mean, sep = "")
  
  indv.wide <- indv.wide[indv.wide$Pi != "pi31",] ## Pi 31 is missing day 35
  
  behav.week.id <- MCMCglmm(cbind(week0, week1, week2, week3, week4, week5, week6, week7, week8) ~ trait - 1, 
                              random = ~us(trait):Pi,
                              rcov = ~us(trait):units,
                              family = c(rep("gaussian", 9)), 
                              prior = prior.cov9, 
                              pr = T, 
                              nitt = 510000, thin = 200, burnin = 10000, 
                              verbose = F,
                              data = indv.wide)
  
  # Model with only ID ----
  id.matrix.week <- matrix(posterior.mode(posterior.cor(behav.week.id$VCV[,1:81])),9,9, 
                             dimnames = list(c("week 1", "week 2", "week 3", "week 4","week 5", "week 6", "week 7", "week 8", "week9"), 
                                             c("week 1", "week 2", "week 3", "week 4","week 5", "week 6", "week 7", "week 8", "week9")))
  
  # now to extract the CI estimates
  ci.week <- data.frame(HPDinterval(posterior.cor(behav.week.id$VCV[,1:81])))
  
  # for corrplot need 3 matrices - estimates, lower CI, upper CI
  lower.week <- matrix(ci.week[,1],9,9)
  upper.week <- matrix(ci.week[,2],9,9)
  
  test <- melt(lower.week) %>%
    mutate(p.value = ifelse(value < 0, 1, 0)) %>%
    select(Var1, Var2, p.value)
  
  p.mat <- diag(9)
  p.mat[cbind(test$Var1, test$Var2)] <- p.mat[cbind(test$Var2, test$Var1)] <- test$p.value
  
  
  week.corr <- ggcorrplot(id.matrix.week, 
                            type = "lower", 
                            p.mat = p.mat, 
                            insig= "blank",
                            colors = c("slateblue4","gray", "mediumorchid1"))
  
  week.corr
  }
## Still no pridictive power here...What is the distribution of cohesion, is mean pdist a stupid metric? 
# With mean distance, there is predictive power, that might be a better place to start then. 

## Kate previously noticed the correlation growing stronger, let's test that
## Result: 4 weeks just isn't enough here, I think we just don't have the sample size to do this

id.matrix.week
# first pull out among-indv corr
df <- melt(replace(id.matrix.week, lower.tri(id.matrix.week, T), NA), na.rm = T)
str(df)

df$start.week <- as.numeric(substr(df$Var1, 6, 8))
df$end.week <- as.numeric(substr(df$Var2, 6, 8))
df$diff <- df$end.week - df$start.week

#now pull out ci for each corr
lower.week <- matrix(ci.week[,1],4,4)
upper.week <- matrix(ci.week[,2],4,4)

test.lower <- melt(replace(lower.week, lower.tri(lower.week, T), NA), na.rm = T)
test.upper <- melt(replace(upper.week, lower.tri(upper.week, T), NA), na.rm = T)

ci.long <- left_join(test.lower, test.upper, by = c("Var1", "Var2")) %>%
  rename(start.week = Var1,
         end.week = Var2,
         lower = value.x, 
         upper = value.y) %>%
  arrange(start.week, end.week)

among.corr <- left_join(df, ci.long, by = c("start.week", "end.week"))



df2 <- among.corr %>%
  filter(diff < 6)

set.seed(421)
corr.mcmc <- MCMCglmm(value ~ start.week*diff, 
                      data = df2, 
                      family = "gaussian", 
                      nitt=510000, burnin = 10000, thin = 200, 
                      verbose = F)

summary(corr.mcmc)

### Multivariate behavioral comparison, let's figure that out: 

## So the distance metrics are correlated, 
## but no correlation between the various cohesion metrics, which isn't ideal....

behav.metric.id <- MCMCglmm(cbind(vel_mean, pDist_mean, dist_mean, NearN_mean, pDistC_mean, velC_mean, head_mean, angleC_mean) ~ trait - 1, 
                          random = ~us(trait):Pi,
                          rcov = ~us(trait):units,
                          family = c(rep("gaussian", 8)), 
                          prior = prior.cov8, 
                          pr = T, 
                          nitt = 510000, thin = 200, burnin = 10000, 
                          verbose = F,
                          data = indv.com)

# Model with only ID ----
id.matrix.metric <- matrix(posterior.mode(posterior.cor(behav.metric.id$VCV[,1:64])),8,8, 
                         dimnames = list(c("Vel.mean", "CenterDist.mean", "ijDistance.mean", "nnDistance.mean", "CenterDistCorr.mean","VelCorr.mean","HeadingCorr.mean","AngleCorr.mean"), 
                                         c("Vel.mean", "CenterDist.mean", "ijDistance.mean", "nnDistance.mean", "CenterDistCorr.mean","VelCorr.mean","HeadingCorr.mean","AngleCorr.mean")))

# now to extract the CI estimates
ci.metric <- data.frame(HPDinterval(posterior.cor(behav.metric.id$VCV[,1:64])))

# for corrplot need 3 matrices - estimates, lower CI, upper CI
lower.metric <- matrix(ci.metric[,1],8,8)
upper.metric <- matrix(ci.metric[,2],8,8)

test <- melt(lower.metric) %>%
  mutate(p.value = ifelse(value < 0, 1, 0)) %>%
  select(Var1, Var2, p.value)

test_low <- melt(upper.week) %>%
  mutate(p.value = ifelse(value > 0, 1, 0)) %>%
  select(Var1, Var2, p.value)

p.mat <- diag(8)
p.mat[cbind(test$Var1, test$Var2)] <- p.mat[cbind(test$Var2, test$Var1)] <- test$p.value


metric.corr <- ggcorrplot(id.matrix.metric, 
                        type = "lower", 
                        p.mat = p.mat, 
                        insig= "blank")#,
                        #colors = c("slateblue4","gray", "mediumorchid1"))

metric.corr

## Ok, so far I've been doing with one or two chosen behaviors,
# But I want to quickly check all the behaviors, so I'm going to define some functions: 

## Define function that calculates DIC for the 3 possible models
func.DIC <- function(depVar,indVar,data) {
  # Null model -----------
  set.seed(58)
  fixed <- as.formula(paste(depVar,"~ ",indVar,sep=""))
  #print('Running simple model')
  model.fixed <- MCMCglmm(fixed = fixed, 
                      data = data, 
                      family = "gaussian", 
                      prior = prior.null, 
                      nitt=510000, burnin = 10000, thin = 200, 
                      verbose = F)
  
  
  # Intercepts ID ----------------------
  set.seed(3432)
  #print('Running random intercept')
  model.randInt <- MCMCglmm(fixed = fixed, 
                      random = ~Pi, 
                      data = data, 
                      family = "gaussian", 
                      prior = prior.id, 
                      nitt=510000, burnin = 10000, thin = 200, 
                      verbose = F)
  
  
  
  # Intercepts and Slopes ID -----------------------
  set.seed(472)
  #print('Running random slope')
  
  random <- as.formula(paste("~us(1 + ",indVar,"):Pi",sep=""))
  model.randSlope <- MCMCglmm(fixed = fixed, 
                      random = random, 
                      data = data, 
                      family = "gaussian", 
                      prior = prior.id.slope, 
                      nitt=510000, burnin = 10000, thin = 200, 
                      verbose = F)
  
  
  print(DIC(model.fixed,model.randInt,model.randSlope))
  print(summary(model.randSlope))
}

func.DIC("dist_mean_","Hour",day27.com)
func.DIC("dist_mean","ExpDay",indv.com)

## Define function that calculates the intercepts on day0
func.hourly.intercepts <- function(depVar,data) {

  fixed <- as.formula(paste(depVar,"~ Hour",sep=""))
  # Null model -----------
  set.seed(58)
  #col.day1.0 <-  MCMCglmm(fixed=fixed, 
  #                         data = day1.com, 
  #                         family = "gaussian", 
  #                         prior = prior.null, 
  #                         nitt=510000, burnin = 10000, thin = 200, 
  #                         verbose = F)
  
  
  # Intercepts ID ----------------------
  set.seed(3432)
  #col.day1.1 <- MCMCglmm(fixed=fixed, 
  #                        random = ~Pi, 
  #                        data = day1.com, 
  #                        family = "gaussian", 
  #                        prior = prior.id, 
  #                        nitt=510000, burnin = 10000, thin = 200, 
  #                        verbose = F)
  
  
  
  # Intercepts and Slopes ID -----------------------
  set.seed(472)
  col.day1.2 <- MCMCglmm(fixed = fixed, 
                          random = ~us(1 + Hour):Pi, 
                          data = data, 
                          family = "gaussian", 
                          prior = prior.id.slope, 
                          nitt=510000, burnin = 10000, thin = 200, 
                          verbose = F)
  
  
  #print(DIC(col.day1.0, col.day1.1, col.day1.2))
  
  ## Big jump with group ID, only Slight improvement from adding Pi as a random slope, but not dramatic. 
  ## I guess this means it's important that hour & day are centered, but that already happened. 
  
  print(summary(col.day1.2))
  print(posterior.mode(col.day1.2$Sol))
  print(HPDinterval(col.day1.2$Sol))
  
  print(posterior.mode(col.day1.2$VCV))
  print(HPDinterval(col.day1.2$VCV))
  
  # only ID repeatability 
  rpt.spd1 <- col.day1.2$VCV[,"(Intercept):(Intercept).Pi"]/(col.day1.2$VCV[,"(Intercept):(Intercept).Pi"]  + 
                                                                col.day1.2$VCV[,"Hour:Hour.Pi"] +
                                                                col.day1.2$VCV[,"units"])
  print('Repeatability:')
  print(posterior.mode(rpt.spd1))
  print(HPDinterval(rpt.spd1))
  
  
  # to get marginal R2 to explain fixed effects variance
  vmVarF<-numeric(2500)
  
  for(i in 1:2500){
    Var<-var(as.vector(col.day1.2$Sol[i,] %*% t(col.day1.2$X)))
    vmVarF[i]<-Var}
  
  
  R2m<-vmVarF/(vmVarF+col.day1.2$VCV[,1]+ col.day1.2$VCV[,4] + col.day1.2$VCV[,5])
  
  print('Marginal R2')
  print(posterior.mode(R2m))
  print(HPDinterval(R2m))
  
  R2c<-(vmVarF + col.day1.2$VCV[,1])/(vmVarF + col.day1.2$VCV[,1] + col.day1.2$VCV[,4] + col.day1.2$VCV[,5])
  
  print('Conditional (?) R2')
  print(posterior.mode(R2c))
  print(HPDinterval(R2c))
  
}
func.hourly.intercepts("dist_mean_",day1.com)


## Define function that calculates the intercepts over multiple days
## Uses variable ~ day + (1+day|pi)


func.triday.intercepts <- function(depVar,data) {
  # select only those IDs in that vector & only keep up to obs 70 & make obs 1 = 0
  indv.30days.cen <- data %>%
    mutate(Day0 = ExpDay, 
           Day3 = ExpDay-3,
           Day6 = ExpDay-6,
           Day9 = ExpDay-9,
           Day12 = ExpDay-12,
           Day15 = ExpDay-15,
           Day18 = ExpDay-18,
           Day21 = ExpDay-21,
           Day24 = ExpDay-24,
           Day27 = ExpDay-27,
           Day30 = ExpDay-30)
  
  #indv.30days.cen <- indv.30days.cen[indv.30days.cen$Pi != "pi12",]
  set.seed(403)
  fixed <- as.formula(paste(depVar,"~ Day0",sep=""))
  print('doing models')
  print('day0')
  col.day0 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Day0):Pi, 
                        data = indv.30days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  set.seed(765)
  fixed <- as.formula(paste(depVar,"~ Day3",sep=""))
  col.day3 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Day3):Pi, 
                        data = indv.30days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Day6",sep=""))
  col.day6 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Day6):Pi, 
                        data = indv.30days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  print('day 9')
  fixed <- as.formula(paste(depVar,"~ Day9",sep=""))
  col.day9 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Day9):Pi, 
                        data = indv.30days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Day12",sep=""))
  col.day12 <- MCMCglmm(fixed = fixed, 
                         random = ~us(1 + Day12):Pi, 
                         data = indv.30days.cen, 
                         family = "gaussian",
                         pr = T,
                         prior = prior.id.slope, 
                         nitt=310000, burnin = 10000, thin = 200, 
                         verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Day15",sep=""))
  col.day15 <- MCMCglmm(fixed = fixed, 
                         random = ~us(1 + Day15):Pi, 
                         data = indv.30days.cen, 
                         family = "gaussian",
                         pr = T,
                         prior = prior.id.slope, 
                         nitt=310000, burnin = 10000, thin = 200, 
                         verbose = F)
  
  print('day 18')
  fixed <- as.formula(paste(depVar,"~ Day18",sep=""))
  col.day18 <- MCMCglmm(fixed = fixed, 
                         random = ~us(1 + Day18):Pi, 
                         data = indv.30days.cen, 
                         family = "gaussian",
                         pr = T,
                         prior = prior.id.slope, 
                         nitt=310000, burnin = 10000, thin = 200, 
                         verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Day21",sep=""))
  col.day21 <- MCMCglmm(fixed = fixed, 
                         random = ~us(1 + Day21):Pi, 
                         data = indv.30days.cen, 
                         family = "gaussian",
                         pr = T,
                         prior = prior.id.slope, 
                         nitt=310000, burnin = 10000, thin = 200, 
                         verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Day24",sep=""))
  col.day24 <- MCMCglmm(fixed = fixed, 
                         random = ~us(1 + Day24):Pi, 
                         data = indv.30days.cen, 
                         family = "gaussian",
                         pr = T,
                         prior = prior.id.slope, 
                         nitt=310000, burnin = 10000, thin = 200, 
                         verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Day27",sep=""))
  col.day27 <- MCMCglmm(fixed = fixed, 
                         random = ~us(1 + Day27):Pi, 
                         data = indv.30days.cen, 
                         family = "gaussian",
                         pr = T,
                         prior = prior.id.slope, 
                         nitt=310000, burnin = 10000, thin = 200, 
                         verbose = F)
  
  print('day 30')
  fixed <- as.formula(paste(depVar,"~ Day30",sep=""))
  col.day30 <- MCMCglmm(fixed = fixed, 
                         random = ~us(1 + Day30):Pi, 
                         data = indv.30days.cen, 
                         family = "gaussian",
                         pr = T,
                         prior = prior.id.slope, 
                         nitt=310000, burnin = 10000, thin = 200, 
                         verbose = F)

  #plot(col.day0$VCV,)
  
  ## Copy pasting so much, getting day: variance
  date <- c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30) ## "It's a magic number"
  
  rpt0 <- col.day0$VCV[,"(Intercept):(Intercept).Pi"]/(col.day0$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.day0$VCV[,"units"]))
  rpt3 <- col.day3$VCV[,"(Intercept):(Intercept).Pi"]/(col.day3$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.day3$VCV[,"units"]))
  rpt6 <- col.day6$VCV[,"(Intercept):(Intercept).Pi"]/(col.day6$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.day6$VCV[,"units"]))
  rpt9 <- col.day9$VCV[,"(Intercept):(Intercept).Pi"]/(col.day9$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.day9$VCV[,"units"]))
  rpt12 <- col.day12$VCV[,"(Intercept):(Intercept).Pi"]/(col.day12$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.day12$VCV[,"units"]))
  rpt15 <- col.day15$VCV[,"(Intercept):(Intercept).Pi"]/(col.day15$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.day15$VCV[,"units"]))
  rpt18 <- col.day18$VCV[,"(Intercept):(Intercept).Pi"]/(col.day18$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.day18$VCV[,"units"]))
  rpt21 <- col.day21$VCV[,"(Intercept):(Intercept).Pi"]/(col.day21$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.day21$VCV[,"units"]))
  rpt24 <- col.day24$VCV[,"(Intercept):(Intercept).Pi"]/(col.day24$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.day24$VCV[,"units"]))
  rpt27 <- col.day27$VCV[,"(Intercept):(Intercept).Pi"]/(col.day27$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.day27$VCV[,"units"]))
  rpt30 <- col.day30$VCV[,"(Intercept):(Intercept).Pi"]/(col.day30$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.day30$VCV[,"units"]))
  
  rpt <- c(posterior.mode(rpt0), 
           posterior.mode(rpt3),
           posterior.mode(rpt6),
           posterior.mode(rpt9),
           posterior.mode(rpt12),
           posterior.mode(rpt15),
           posterior.mode(rpt18),
           posterior.mode(rpt21),
           posterior.mode(rpt24),
           posterior.mode(rpt27),
           posterior.mode(rpt30))
  
  
  ci.rpt <- c(HPDinterval(rpt0)[1:2],
              HPDinterval(rpt3)[1:2],
              HPDinterval(rpt6)[1:2],
              HPDinterval(rpt9)[1:2],
              HPDinterval(rpt12)[1:2],
              HPDinterval(rpt15)[1:2],
              HPDinterval(rpt18)[1:2],
              HPDinterval(rpt21)[1:2],
              HPDinterval(rpt24)[1:2],
              HPDinterval(rpt27)[1:2],
              HPDinterval(rpt30)[1:2])
  ci.rpt <- matrix(ci.rpt, nrow = 11, byrow = T)
  
  post.id <- c(posterior.mode(col.day0$VCV[,"(Intercept):(Intercept).Pi"]), 
               posterior.mode(col.day3$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.day6$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.day9$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.day12$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.day15$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.day18$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.day21$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.day24$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.day27$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.day30$VCV[,"(Intercept):(Intercept).Pi"]))
  
  ci.id <- c(HPDinterval(col.day0$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.day3$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.day6$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.day9$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.day12$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.day15$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.day18$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.day21$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.day24$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.day27$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.day30$VCV[,"(Intercept):(Intercept).Pi"])[1:2])
  
  ci.id <-matrix(ci.id, nrow = 11, byrow = T)
  
  post.w <- c(posterior.mode(col.day0$VCV[,"units"]),
              posterior.mode(col.day3$VCV[,"units"]),
              posterior.mode(col.day6$VCV[,"units"]),
              posterior.mode(col.day9$VCV[,"units"]),
              posterior.mode(col.day12$VCV[,"units"]),
              posterior.mode(col.day15$VCV[,"units"]),
              posterior.mode(col.day18$VCV[,"units"]),
              posterior.mode(col.day21$VCV[,"units"]),
              posterior.mode(col.day24$VCV[,"units"]),
              posterior.mode(col.day27$VCV[,"units"]),
              posterior.mode(col.day30$VCV[,"units"]))
  
  ci.w <- c(HPDinterval(col.day0$VCV[,"units"])[1:2],   
            HPDinterval(col.day3$VCV[,"units"])[1:2], 
            HPDinterval(col.day6$VCV[,"units"])[1:2], 
            HPDinterval(col.day9$VCV[,"units"])[1:2], 
            HPDinterval(col.day12$VCV[,"units"])[1:2], 
            HPDinterval(col.day15$VCV[,"units"])[1:2], 
            HPDinterval(col.day18$VCV[,"units"])[1:2], 
            HPDinterval(col.day21$VCV[,"units"])[1:2], 
            HPDinterval(col.day24$VCV[,"units"])[1:2], 
            HPDinterval(col.day27$VCV[,"units"])[1:2], 
            HPDinterval(col.day30$VCV[,"units"])[1:2])
  
  ci.w <- matrix(ci.w, nrow = 11, byrow = T)            
  #plot(post.id ~ date)
  #plot(post.w ~ date)
  #plot(rpt ~ date)
  
  rpt.slice.wide <- data.frame(date, rpt, "lower.rpt" = ci.rpt[,1], "upper.rpt" = ci.rpt[,2],
                               post.id, "lower.id" = ci.id[,1], "upper.id" = ci.id[,2], 
                               post.w, "lower.w" = ci.w[,1], "upper.w" = ci.w[,2])
  
  rpt.slice.long <- data.frame("date" = rep(date, 3), 
                               "type" = rep(c("rpt", "id", "within"), each = 11),
                               "variance" = unname(c(rpt, post.id, post.w)),
                               "lower" = unname(c(ci.rpt[,1], ci.id[,1], ci.w[,1])),
                               "upper" = unname(c(ci.rpt[,2], ci.id[,2], ci.w[,2])))
  
  
  rpt.slice.rpt <- rpt.slice.long[rpt.slice.long$type == 'rpt',]
  
  ### This is very dataset specific, you may need to check it: 
  ids <- colnames(col.day0$Sol)[3:14] ## Make sure this matches below
  ids <- substr(ids, 16, 19)
  
  'this pulls out the individual intercepts and adds in the overall intercepts
  
  so that way these numbers are absolute values, as opposed to differences from overall'
  
  intercepts0 <- unname(posterior.mode(col.day0$Sol)[3:14] + posterior.mode(col.day0$Sol)["(Intercept)"])
  intercepts3 <- unname(posterior.mode(col.day3$Sol)[3:14] + posterior.mode(col.day0$Sol)["(Intercept)"])
  intercepts6 <- unname(posterior.mode(col.day6$Sol)[3:14] + posterior.mode(col.day6$Sol)["(Intercept)"])
  intercepts9 <- unname(posterior.mode(col.day9$Sol)[3:14] + posterior.mode(col.day9$Sol)["(Intercept)"])
  intercepts12 <- unname(posterior.mode(col.day12$Sol)[3:14] + posterior.mode(col.day12$Sol)["(Intercept)"])
  intercepts15 <- unname(posterior.mode(col.day15$Sol)[3:14] + posterior.mode(col.day15$Sol)["(Intercept)"])
  intercepts18 <- unname(posterior.mode(col.day18$Sol)[3:14] + posterior.mode(col.day18$Sol)["(Intercept)"])
  intercepts21 <- unname(posterior.mode(col.day21$Sol)[3:14] + posterior.mode(col.day21$Sol)["(Intercept)"])
  intercepts24 <- unname(posterior.mode(col.day24$Sol)[3:14] + posterior.mode(col.day24$Sol)["(Intercept)"])
  intercepts27 <- unname(posterior.mode(col.day27$Sol)[3:14] + posterior.mode(col.day27$Sol)["(Intercept)"])
  intercepts30 <- unname(posterior.mode(col.day30$Sol)[3:14] + posterior.mode(col.day30$Sol)["(Intercept)"])
  
  date <- rep(c(0,3,6,9,12,15,18,21,24,27,30), each = 12)
  picomp <- rep(ids, 11)
  blup <- c(intercepts0, intercepts3, intercepts6, intercepts9, intercepts12, intercepts15, 
            intercepts18, intercepts21, intercepts24, intercepts27, intercepts30)
  
  pred.intercepts <- data.frame(date, picomp, blup)
  

  plt.repeat <- ggplot(rpt.slice.rpt,aes(x=date,y=variance)) + 
    geom_point(shape=1) + 
    #geom_line() +
    #geom_errorbar(aes(ymin=lower,ymax=upper)) + 
    ggtitle(depVar) +
    #ylim(0,1) + 
    ylab('Repeatability') + 
    xlab('Days since birth') + 
    theme_classic() + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 14)) 
  ylims <- c(0,1)
  if (depVar == 'dist_mean') {
    ylims <- c(0,500)
  }

  plt.day <- ggplot(pred.intercepts, aes(x=date, y=blup, group=picomp, color=picomp)) + 
    geom_line() + 
    ggtitle(depVar) +
    theme_classic() + 
    #ylim(ylims) + 
    ylab(depVar) + 
    xlab('Days since birth') + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 14))
  
  return(list(plt.day,plt.repeat))
}

testplots <- func.triday.intercepts("dist_mean",indv.com)
testplots[1]
testplots[2]
indv.com.no12 <- indv.com[indv.com$Pi != 'pi12',]
testplot <- func.triday.intercepts("dist_mean",indv.com.no12)

testplot

## What if we do the longer pi's? Still decreasing
depVar <- 'dist_mean'
func.weekly.intercepts <- function(depVar,data) {
# select only those IDs in that vector & only keep up to obs 70 & make obs 1 = 0
  indv.60days.cen <- data %>%
    mutate(Week0 = ExpDay, 
           Week1 = ExpDay-7,
           Week2 = ExpDay-14,
           Week3 = ExpDay-21,
           Week4 = ExpDay-28,
           Week5 = ExpDay-35,
           Week6 = ExpDay-42,
           Week7 = ExpDay-49,
           Week8 = ExpDay-56,
           Week9 = ExpDay-63)
  
  #indv.60days.cen <- indv.36days.cen[indv.60days.cen$Pi != "pi12",]
  set.seed(403)
  fixed <- as.formula(paste(depVar,"~ Week0",sep=""))
  print('week0')
  col.week0 <- MCMCglmm(fixed = fixed, 
                       random = ~us(1 + Week0):Pi, 
                       data = indv.60days.cen, 
                       family = "gaussian",
                       pr = T,
                       prior = prior.id.slope, 
                       nitt=310000, burnin = 10000, thin = 200, 
                       verbose = F)
  
  set.seed(765)
  fixed <- as.formula(paste(depVar,"~ Week1",sep=""))
  col.week1 <- MCMCglmm(fixed = fixed, 
                       random = ~us(1 + Week1):Pi, 
                       data = indv.60days.cen, 
                       family = "gaussian",
                       pr = T,
                       prior = prior.id.slope, 
                       nitt=310000, burnin = 10000, thin = 200, 
                       verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Week2",sep=""))
  col.week2 <- MCMCglmm(fixed = fixed, 
                       random = ~us(1 + Week2):Pi, 
                       data = indv.60days.cen, 
                       family = "gaussian",
                       pr = T,
                       prior = prior.id.slope, 
                       nitt=310000, burnin = 10000, thin = 200, 
                       verbose = F)
  
  print('week3')
  fixed <- as.formula(paste(depVar,"~ Week3",sep=""))
  col.week3 <- MCMCglmm(fixed = fixed, 
                       random = ~us(1 + Week3):Pi, 
                       data = indv.60days.cen, 
                       family = "gaussian",
                       pr = T,
                       prior = prior.id.slope, 
                       nitt=310000, burnin = 10000, thin = 200, 
                       verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Week4",sep=""))
  col.week4 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Week4):Pi, 
                        data = indv.60days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Week5",sep=""))
  col.week5 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Week5):Pi, 
                        data = indv.60days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Week6",sep=""))
  col.week6 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Week6):Pi, 
                        data = indv.60days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Week7",sep=""))
  col.week7 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Week7):Pi, 
                        data = indv.60days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Week8",sep=""))
  col.week8 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Week8):Pi, 
                        data = indv.60days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Week9",sep=""))
  col.week9 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Week9):Pi, 
                        data = indv.60days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  print('done with models')
  
  #plot(col.week0$VCV,)
  
  ## Copy pasting so much, getting day: variance
  #date <- c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30) ## "It's a magic number"
  date <- c(0,7,14,21,28,35,42,49,56,63)
  
  rpt0 <- col.week0$VCV[,"(Intercept):(Intercept).Pi"]/(col.week0$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.week0$VCV[,"units"]))
  rpt1 <- col.week1$VCV[,"(Intercept):(Intercept).Pi"]/(col.week1$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.week1$VCV[,"units"]))
  rpt2 <- col.week2$VCV[,"(Intercept):(Intercept).Pi"]/(col.week2$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.week2$VCV[,"units"]))
  rpt3 <- col.week3$VCV[,"(Intercept):(Intercept).Pi"]/(col.week3$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.week3$VCV[,"units"]))
  rpt4 <- col.week4$VCV[,"(Intercept):(Intercept).Pi"]/(col.week4$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.week4$VCV[,"units"]))
  rpt5 <- col.week5$VCV[,"(Intercept):(Intercept).Pi"]/(col.week5$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.week5$VCV[,"units"]))
  rpt6 <- col.week6$VCV[,"(Intercept):(Intercept).Pi"]/(col.week6$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.week6$VCV[,"units"]))
  rpt7 <- col.week7$VCV[,"(Intercept):(Intercept).Pi"]/(col.week7$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.week7$VCV[,"units"]))
  rpt8 <- col.week8$VCV[,"(Intercept):(Intercept).Pi"]/(col.week8$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.week8$VCV[,"units"]))
  rpt9 <- col.week9$VCV[,"(Intercept):(Intercept).Pi"]/(col.week9$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.week9$VCV[,"units"]))
  
  rpt <- c(posterior.mode(rpt0), 
           posterior.mode(rpt1),
           posterior.mode(rpt2),
           posterior.mode(rpt3),
           posterior.mode(rpt4),
           posterior.mode(rpt5),
           posterior.mode(rpt6),
           posterior.mode(rpt7),
           posterior.mode(rpt8),
           posterior.mode(rpt9))
  
  
  ci.rpt <- c(HPDinterval(rpt0)[1:2],
              HPDinterval(rpt1)[1:2],
              HPDinterval(rpt2)[1:2],
              HPDinterval(rpt3)[1:2],
              HPDinterval(rpt4)[1:2],
              HPDinterval(rpt5)[1:2],
              HPDinterval(rpt6)[1:2],
              HPDinterval(rpt7)[1:2],
              HPDinterval(rpt8)[1:2],
              HPDinterval(rpt9)[1:2])
  ci.rpt <- matrix(ci.rpt, nrow = 10, byrow = T)
  
  post.id <- c(posterior.mode(col.week0$VCV[,"(Intercept):(Intercept).Pi"]), 
               posterior.mode(col.week1$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.week2$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.week3$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.week4$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.week5$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.week6$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.week7$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.week8$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.week9$VCV[,"(Intercept):(Intercept).Pi"]))
  
  ci.id <- c(HPDinterval(col.week0$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.week1$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.week2$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.week3$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.week4$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.week5$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.week6$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.week7$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.week8$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.week9$VCV[,"(Intercept):(Intercept).Pi"])[1:2])
  ci.id <-matrix(ci.id, nrow = 10, byrow = T)
  
  post.w <- c(posterior.mode(col.week0$VCV[,"units"]),
              posterior.mode(col.week1$VCV[,"units"]),
              posterior.mode(col.week2$VCV[,"units"]),
              posterior.mode(col.week3$VCV[,"units"]),
              posterior.mode(col.week4$VCV[,"units"]),
              posterior.mode(col.week5$VCV[,"units"]),
              posterior.mode(col.week6$VCV[,"units"]),
              posterior.mode(col.week7$VCV[,"units"]),
              posterior.mode(col.week8$VCV[,"units"]),
              posterior.mode(col.week9$VCV[,"units"]))
  
  ci.w <- c(HPDinterval(col.week0$VCV[,"units"])[1:2],   
            HPDinterval(col.week1$VCV[,"units"])[1:2], 
            HPDinterval(col.week2$VCV[,"units"])[1:2], 
            HPDinterval(col.week3$VCV[,"units"])[1:2], 
            HPDinterval(col.week4$VCV[,"units"])[1:2], 
            HPDinterval(col.week5$VCV[,"units"])[1:2], 
            HPDinterval(col.week6$VCV[,"units"])[1:2], 
            HPDinterval(col.week7$VCV[,"units"])[1:2], 
            HPDinterval(col.week8$VCV[,"units"])[1:2], 
            HPDinterval(col.week9$VCV[,"units"])[1:2]) 
  
  ci.w <- matrix(ci.w, nrow = 10, byrow = T)            
  
  #plot(post.id ~ date)
  #plot(post.w ~ date)
  #plot(rpt ~ date)
  
  
  rpt.slice.wide <- data.frame(date, rpt, "lower.rpt" = ci.rpt[,1], "upper.rpt" = ci.rpt[,2],
                               post.id, "lower.id" = ci.id[,1], "upper.id" = ci.id[,2], 
                               post.w, "lower.w" = ci.w[,1], "upper.w" = ci.w[,2])
  
  rpt.slice.long <- data.frame("date" = rep(date, 3), 
                               "type" = rep(c("rpt", "id", "within"), each = 10),
                               "variance" = unname(c(rpt, post.id, post.w)),
                               "lower" = unname(c(ci.rpt[,1], ci.id[,1], ci.w[,1])),
                               "upper" = unname(c(ci.rpt[,2], ci.id[,2], ci.w[,2])))
  
  ### This is very dataset specific, you may need to check it: 
  ids <- colnames(col.week0$Sol)[3:10] ## Make sure this matches below
  ids <- substr(ids, 16, 19)
  
  'this pulls out the individual intercepts and adds in the overall intercepts
  
  so that way these numbers are absolute values, as opposed to differences from overall'
  
  intercepts0 <- unname(posterior.mode(col.week0$Sol)[3:10] + posterior.mode(col.week0$Sol)["(Intercept)"])
  intercepts1 <- unname(posterior.mode(col.week1$Sol)[3:10] + posterior.mode(col.week1$Sol)["(Intercept)"])
  intercepts2 <- unname(posterior.mode(col.week2$Sol)[3:10] + posterior.mode(col.week2$Sol)["(Intercept)"])
  intercepts3 <- unname(posterior.mode(col.week3$Sol)[3:10] + posterior.mode(col.week3$Sol)["(Intercept)"])
  intercepts4 <- unname(posterior.mode(col.week4$Sol)[3:10] + posterior.mode(col.week4$Sol)["(Intercept)"])
  intercepts5 <- unname(posterior.mode(col.week5$Sol)[3:10] + posterior.mode(col.week5$Sol)["(Intercept)"])
  intercepts6 <- unname(posterior.mode(col.week6$Sol)[3:10] + posterior.mode(col.week6$Sol)["(Intercept)"])
  intercepts7 <- unname(posterior.mode(col.week7$Sol)[3:10] + posterior.mode(col.week7$Sol)["(Intercept)"])
  intercepts8 <- unname(posterior.mode(col.week8$Sol)[3:10] + posterior.mode(col.week8$Sol)["(Intercept)"])
  intercepts9 <- unname(posterior.mode(col.week9$Sol)[3:10] + posterior.mode(col.week9$Sol)["(Intercept)"])
  
  date <- rep(c(0,7,14,21,28,35,42,49,56,63), each = 8)
  picomp <- rep(ids, 10)
  blup <- c(intercepts0, intercepts1, intercepts2, intercepts3, intercepts4, intercepts5, 
            intercepts6, intercepts7, intercepts8, intercepts9)
  
  pred.intercepts <- data.frame(date, picomp, blup)
  
  rpt.slice.rpt <- rpt.slice.long[rpt.slice.long$type == 'rpt',]
  print(dim(pred.intercepts))
  print(dim(ci.rpt))
  print(blup)
  plt.repeat <- ggplot(rpt.slice.rpt,aes(x=date,y=variance)) + 
    geom_point(shape=1,size=4) + 
    #geom_line() +
    #geom_errorbar(aes(ymin=lower,ymax=upper)) + 
    ggtitle(depVar) +
    #ylim(0,1) + 
    ylab('Repeatability') + 
    xlab('Days since birth') + 
    theme_classic() + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 14)) 
  ylims <- c(0,1)
  if (depVar == 'dist_mean') {
    ylims <- c(0,500)
  }
  plt.week <- ggplot(pred.intercepts, aes(x=date, y=blup, group=picomp, color=picomp)) + 
    geom_line() + 
    ggtitle(depVar) +
    theme_classic() + 
    #ylim(ylims) + 
    ylab(depVar) + 
    xlab('Days since birth') + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 14)) 
  return(list(plt.week,plt.repeat))
}


func.fiver.intercepts <- function(depVar,data) {
  # select only those IDs in that vector & only keep up to obs 70 & make obs 1 = 0
  indv.54days.cen <- data %>%
    mutate(Fiver0 = ExpDay, 
           Fiver1 = ExpDay-5,
           Fiver2 = ExpDay-10,
           Fiver3 = ExpDay-15,
           Fiver4 = ExpDay-20,
           Fiver5 = ExpDay-25,
           Fiver6 = ExpDay-30,
           Fiver7 = ExpDay-35,
           Fiver8 = ExpDay-40,
           Fiver9 = ExpDay-45,
           Fiver10 = ExpDay-50)
  
  #indv.60days.cen <- indv.36days.cen[indv.60days.cen$Pi != "pi12",]
  set.seed(403)
  fixed <- as.formula(paste(depVar,"~ Fiver0",sep=""))
  print('fiver0')
  col.fiver0 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Fiver0):Pi, 
                        data = indv.54days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  set.seed(765)
  fixed <- as.formula(paste(depVar,"~ Fiver1",sep=""))
  col.fiver1 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Fiver1):Pi, 
                        data = indv.54days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Fiver2",sep=""))
  col.fiver2 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Fiver2):Pi, 
                        data = indv.54days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  print('Fiver3')
  fixed <- as.formula(paste(depVar,"~ Fiver3",sep=""))
  col.fiver3 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Fiver3):Pi, 
                        data = indv.54days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Fiver4",sep=""))
  col.fiver4 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Fiver4):Pi, 
                        data = indv.54days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Fiver5",sep=""))
  col.fiver5 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Fiver5):Pi, 
                        data = indv.54days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Fiver6",sep=""))
  col.fiver6 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Fiver6):Pi, 
                        data = indv.54days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Fiver7",sep=""))
  col.fiver7 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Fiver7):Pi, 
                        data = indv.54days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Fiver8",sep=""))
  col.fiver8 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Fiver8):Pi, 
                        data = indv.54days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  fixed <- as.formula(paste(depVar,"~ Fiver9",sep=""))
  col.fiver9 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Fiver9):Pi, 
                        data = indv.54days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  fixed <- as.formula(paste(depVar,"~ Fiver10",sep=""))
  col.fiver10 <- MCMCglmm(fixed = fixed, 
                        random = ~us(1 + Fiver10):Pi, 
                        data = indv.54days.cen, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.id.slope, 
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)
  
  print('done with models')
  
  #plot(col.week0$VCV,)
  
  ## Copy pasting so much, getting day: variance
  #date <- c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30) ## "It's a magic number"
  #date <- c(0,7,14,21,28,35,42,49,56,63)
  date <- c(0,5,10,15,20,25,30,35,40,45,50)
  
  rpt0 <- col.fiver0$VCV[,"(Intercept):(Intercept).Pi"]/(col.fiver0$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.fiver0$VCV[,"units"]))
  rpt1 <- col.fiver1$VCV[,"(Intercept):(Intercept).Pi"]/(col.fiver1$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.fiver1$VCV[,"units"]))
  rpt2 <- col.fiver2$VCV[,"(Intercept):(Intercept).Pi"]/(col.fiver2$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.fiver2$VCV[,"units"]))
  rpt3 <- col.fiver3$VCV[,"(Intercept):(Intercept).Pi"]/(col.fiver3$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.fiver3$VCV[,"units"]))
  rpt4 <- col.fiver4$VCV[,"(Intercept):(Intercept).Pi"]/(col.fiver4$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.fiver4$VCV[,"units"]))
  rpt5 <- col.fiver5$VCV[,"(Intercept):(Intercept).Pi"]/(col.fiver5$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.fiver5$VCV[,"units"]))
  rpt6 <- col.fiver6$VCV[,"(Intercept):(Intercept).Pi"]/(col.fiver6$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.fiver6$VCV[,"units"]))
  rpt7 <- col.fiver7$VCV[,"(Intercept):(Intercept).Pi"]/(col.fiver7$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.fiver7$VCV[,"units"]))
  rpt8 <- col.fiver8$VCV[,"(Intercept):(Intercept).Pi"]/(col.fiver8$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.fiver8$VCV[,"units"]))
  rpt9 <- col.fiver9$VCV[,"(Intercept):(Intercept).Pi"]/(col.fiver9$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.fiver9$VCV[,"units"]))
  rpt10 <- col.fiver10$VCV[,"(Intercept):(Intercept).Pi"]/(col.fiver10$VCV[,"(Intercept):(Intercept).Pi"] + posterior.mode(col.fiver10$VCV[,"units"]))
  
  rpt <- c(posterior.mode(rpt0), 
           posterior.mode(rpt1),
           posterior.mode(rpt2),
           posterior.mode(rpt3),
           posterior.mode(rpt4),
           posterior.mode(rpt5),
           posterior.mode(rpt6),
           posterior.mode(rpt7),
           posterior.mode(rpt8),
           posterior.mode(rpt9),
           posterior.mode(rpt10))
  
  
  ci.rpt <- c(HPDinterval(rpt0)[1:2],
              HPDinterval(rpt1)[1:2],
              HPDinterval(rpt2)[1:2],
              HPDinterval(rpt3)[1:2],
              HPDinterval(rpt4)[1:2],
              HPDinterval(rpt5)[1:2],
              HPDinterval(rpt6)[1:2],
              HPDinterval(rpt7)[1:2],
              HPDinterval(rpt8)[1:2],
              HPDinterval(rpt9)[1:2],
              HPDinterval(rpt10)[1:2])
  ci.rpt <- matrix(ci.rpt, nrow = 11, byrow = T)
  
  post.id <- c(posterior.mode(col.fiver0$VCV[,"(Intercept):(Intercept).Pi"]), 
               posterior.mode(col.fiver1$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.fiver2$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.fiver3$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.fiver4$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.fiver5$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.fiver6$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.fiver7$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.fiver8$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.fiver9$VCV[,"(Intercept):(Intercept).Pi"]),
               posterior.mode(col.fiver10$VCV[,"(Intercept):(Intercept).Pi"]))
  
  ci.id <- c(HPDinterval(col.fiver0$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.fiver1$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.fiver2$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.fiver3$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.fiver4$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.fiver5$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.fiver6$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.fiver7$VCV[,"(Intercept):(Intercept).Pi"])[1:2], 
             HPDinterval(col.fiver8$VCV[,"(Intercept):(Intercept).Pi"])[1:2],
             HPDinterval(col.fiver9$VCV[,"(Intercept):(Intercept).Pi"])[1:2],
             HPDinterval(col.fiver10$VCV[,"(Intercept):(Intercept).Pi"])[1:2])
  ci.id <-matrix(ci.id, nrow = 11, byrow = T)
  
  post.w <- c(posterior.mode(col.fiver0$VCV[,"units"]),
              posterior.mode(col.fiver1$VCV[,"units"]),
              posterior.mode(col.fiver2$VCV[,"units"]),
              posterior.mode(col.fiver3$VCV[,"units"]),
              posterior.mode(col.fiver4$VCV[,"units"]),
              posterior.mode(col.fiver5$VCV[,"units"]),
              posterior.mode(col.fiver6$VCV[,"units"]),
              posterior.mode(col.fiver7$VCV[,"units"]),
              posterior.mode(col.fiver8$VCV[,"units"]),
              posterior.mode(col.fiver9$VCV[,"units"]),
              posterior.mode(col.fiver10$VCV[,"units"]))
  
  ci.w <- c(HPDinterval(col.fiver0$VCV[,"units"])[1:2],   
            HPDinterval(col.fiver1$VCV[,"units"])[1:2], 
            HPDinterval(col.fiver2$VCV[,"units"])[1:2], 
            HPDinterval(col.fiver3$VCV[,"units"])[1:2], 
            HPDinterval(col.fiver4$VCV[,"units"])[1:2], 
            HPDinterval(col.fiver5$VCV[,"units"])[1:2], 
            HPDinterval(col.fiver6$VCV[,"units"])[1:2], 
            HPDinterval(col.fiver7$VCV[,"units"])[1:2], 
            HPDinterval(col.fiver8$VCV[,"units"])[1:2], 
            HPDinterval(col.fiver9$VCV[,"units"])[1:2],
            HPDinterval(col.fiver10$VCV[,"units"])[1:2]) 
  
  ci.w <- matrix(ci.w, nrow = 11, byrow = T)            
  
  #plot(post.id ~ date)
  #plot(post.w ~ date)
  #plot(rpt ~ date)
  
  
  rpt.slice.wide <- data.frame(date, rpt, "lower.rpt" = ci.rpt[,1], "upper.rpt" = ci.rpt[,2],
                               post.id, "lower.id" = ci.id[,1], "upper.id" = ci.id[,2], 
                               post.w, "lower.w" = ci.w[,1], "upper.w" = ci.w[,2])
  
  rpt.slice.long <- data.frame("date" = rep(date, 3), 
                               "type" = rep(c("rpt", "id", "within"), each = 11),
                               "variance" = unname(c(rpt, post.id, post.w)),
                               "lower" = unname(c(ci.rpt[,1], ci.id[,1], ci.w[,1])),
                               "upper" = unname(c(ci.rpt[,2], ci.id[,2], ci.w[,2])))
  
  ### This is very dataset specific, you may need to check it: 
  ids <- colnames(col.fiver0$Sol)[3:12] ## Make sure this matches below
  ids <- substr(ids, 16, 19)
  
  'this pulls out the individual intercepts and adds in the overall intercepts
  
  so that way these numbers are absolute values, as opposed to differences from overall'
  
  intercepts0 <- unname(posterior.mode(col.fiver0$Sol)[3:12] + posterior.mode(col.fiver0$Sol)["(Intercept)"])
  intercepts1 <- unname(posterior.mode(col.fiver1$Sol)[3:12] + posterior.mode(col.fiver0$Sol)["(Intercept)"])
  intercepts2 <- unname(posterior.mode(col.fiver2$Sol)[3:12] + posterior.mode(col.fiver2$Sol)["(Intercept)"])
  intercepts3 <- unname(posterior.mode(col.fiver3$Sol)[3:12] + posterior.mode(col.fiver3$Sol)["(Intercept)"])
  intercepts4 <- unname(posterior.mode(col.fiver4$Sol)[3:12] + posterior.mode(col.fiver4$Sol)["(Intercept)"])
  intercepts5 <- unname(posterior.mode(col.fiver5$Sol)[3:12] + posterior.mode(col.fiver5$Sol)["(Intercept)"])
  intercepts6 <- unname(posterior.mode(col.fiver6$Sol)[3:12] + posterior.mode(col.fiver6$Sol)["(Intercept)"])
  intercepts7 <- unname(posterior.mode(col.fiver7$Sol)[3:12] + posterior.mode(col.fiver7$Sol)["(Intercept)"])
  intercepts8 <- unname(posterior.mode(col.fiver8$Sol)[3:12] + posterior.mode(col.fiver8$Sol)["(Intercept)"])
  intercepts9 <- unname(posterior.mode(col.fiver9$Sol)[3:12] + posterior.mode(col.fiver9$Sol)["(Intercept)"])
  intercepts10 <- unname(posterior.mode(col.fiver10$Sol)[3:12] + posterior.mode(col.fiver10$Sol)["(Intercept)"])
  
  #date <- rep(c(0,7,14,21,28,35,42,49,56,63), each = 8)
  date <- rep(c(0,5,10,15,20,25,30,35,40,45,50), each=10)
  picomp <- rep(ids, 11)
  blup <- c(intercepts0, intercepts1, intercepts2, intercepts3, intercepts4, intercepts5, 
            intercepts6, intercepts7, intercepts8, intercepts9,intercepts10)
  
  pred.intercepts <- data.frame(date, picomp, blup)
  
  rpt.slice.rpt <- rpt.slice.long[rpt.slice.long$type == 'rpt',]
  plt.repeat <- ggplot(rpt.slice.rpt,aes(x=date,y=variance)) + 
    geom_point(shape=1,size=4) + 
    #geom_line() +
    #geom_errorbar(aes(ymin=lower,ymax=upper)) + 
    ggtitle(depVar) +
    #ylim(0,1) + 
    ylab('Repeatability') + 
    xlab('Days since birth') + 
    theme_classic() + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 14)) 
  ylims <- c(0,1)
  if (depVar == 'dist_mean') {
    ylims <- c(0,500)
  }
  plt.week <- ggplot(pred.intercepts, aes(x=date, y=blup, group=picomp, color=picomp)) + 
    geom_line() + 
    ggtitle(depVar) +
    theme_classic() + 
    #ylim(ylims) + 
    ylab(depVar) + 
    xlab('Days since birth') + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 14)) 
  return(list(plt.week,plt.repeat))
}

test.slice.long54 <- func.fiver.intercepts("dist_mean",indv.long54)
test.slice.long54[2]
test.slice.long54[1]

test.slice.long54 <- func.fiver.intercepts("velC_mean",indv.long54)
test.slice.long54[2]
test.slice.long54[1]

test.slice.long <- func.weekly.intercepts("dist_mean",indv.com)
test.slice.long[2]

test.slice.rpt <- test.slice.long[test.slice.long$type == 'rpt',]
plt.repeat <- ggplot(test.slice.rpt,aes(x=date,y=variance)) + geom_line() +
  geom_errorbar(aes(ymin=lower,ymax=upper)) + 
  ggtitle("DepVar") +
  ylim(0,1) + 
  ylab('Repeatability') + 
  xlab('Days since birth') + 
    theme_classic() + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 14)) 
plt.repeat
#testplot.week <- func.weekly.intercepts("dist_mean",indv.com)
plots <- func.weekly.intercepts("dist_mean",indv.long)
testplot.week <- plots[[1]]
testplot.repeat <- plots[[2]]
testplot.week
testplot.repeat



testplot.week
func.hourly.repeat <- function(depVar,data) { 
  
  fixed <- as.formula(paste(depVar," ~ Hour",sep=""))
  var.day <- MCMCglmm(fixed = fixed, 
                          random = ~us(1 + Hour):Pi, 
                          data = data, 
                          family = "gaussian", 
                          prior = prior.id.slope, 
                          nitt=510000, burnin = 10000, thin = 200, 
                          verbose = F)
  
  ## Big jump with group ID, only Slight improvement from adding Pi as a random slope, but not dramatic. 
  ## I guess this means it's important that hour & day are centered, but that already happened. 
  
  #summary(var.day)
  posterior.mode(var.day$Sol)
  HPDinterval(var.day$Sol)
  
  posterior.mode(var.day$VCV)
  HPDinterval(var.day$VCV)
  
  # only ID repeatability 
  rpt.spd1 <- var.day$VCV[,"(Intercept):(Intercept).Pi"]/(var.day$VCV[,"(Intercept):(Intercept).Pi"]  + 
                                                                var.day$VCV[,"Hour:Hour.Pi"] +
                                                                var.day$VCV[,"units"])
  print(posterior.mode(rpt.spd1))
  print(HPDinterval(rpt.spd1))
  
  # to get marginal R2 to explain fixed effects variance
  #vmVarF<-numeric(2500)
  
  #for(i in 1:2500){
  #  Var<-var(as.vector(var.day$Sol[i,] %*% t(var.day$X)))
  #  vmVarF[i]<-Var}
  
  
  #R2m<-vmVarF/(vmVarF+var.day$VCV[,1]+ var.day$VCV[,4] + var.day$VCV[,5])
  
  #posterior.mode(R2m)
  #HPDinterval(R2m)
  
  #R2c<-(vmVarF + var.day$VCV[,1])/(vmVarF + var.day$VCV[,1] + var.day$VCV[,4] + var.day$VCV[,5])
  
  #posterior.mode(R2c)
  #HPDinterval(R2c)
}

func.hourly.repeat('dist_mean_',day1.com)


func.overall.repeat <- function(depVar,data) { 
  fixed <- as.formula(paste(depVar," ~ ExpDay",sep=""))
  var.model <-  MCMCglmm(fixed = fixed,
                            random = ~us(1 + ExpDay):Pi, 
                            data = data, 
                            family = "gaussian", 
                            prior = prior.id.slope, 
                            nitt=510000, burnin = 10000, thin = 200, 
                            verbose = F)
  # only ID repeatability 
  rpt.spd1 <- var.model$VCV[,"(Intercept):(Intercept).Pi"]/(var.model$VCV[,"(Intercept):(Intercept).Pi"]  + 
                                                              var.model$VCV[,"ExpDay:ExpDay.Pi"] +
                                                              var.model$VCV[,"units"])
  print(posterior.mode(rpt.spd1))
  print(HPDinterval(rpt.spd1))
  }
## This wants wide data, which you probably calculated above
func.weekly.predict <- function(depVar,data) {
  
  indv.wide <- indv.com[indv.com$week < 4,] %>%
    mutate(
      obs.day = case_when(
        ExpDay %in% seq(0,27,7) ~ 1, 
        ExpDay %in% seq(1,27,7) ~ 2,
        ExpDay %in% seq(2,27,7) ~ 3,
        ExpDay %in% seq(3,27,7) ~ 4,
        ExpDay %in% seq(4,27,7) ~ 5,
        ExpDay %in% seq(5,27,7) ~ 6,
        ExpDay %in% seq(6,27,7) ~ 7
      )) %>%
    select(Pi, week, all_of(depVar), obs.day) %>%
    spread(week, depVar, sep = "")
  
  behav.week.id <- MCMCglmm(cbind(week0, week1, week2, week3) ~ trait - 1, 
                            random = ~us(trait):Pi,
                            rcov = ~us(trait):units,
                            family = c(rep("gaussian", 4)), 
                            prior = prior.cov4, 
                            pr = T, 
                            nitt = 510000, thin = 200, burnin = 10000, 
                            verbose = F,
                            data = indv.wide)
  
  # Model with only ID ----
  id.matrix.week <- matrix(posterior.mode(posterior.cor(behav.week.id$VCV[,1:16])),4,4, 
                           dimnames = list(c("week 1", "week 2", "week 3", "week 4"), 
                                           c("week 1", "week 2", "week 3", "week 4")))
  
  # now to extract the CI estimates
  ci.week <- data.frame(HPDinterval(posterior.cor(behav.week.id$VCV[,1:16])))
  
  # for corrplot need 3 matrices - estimates, lower CI, upper CI
  lower.week <- matrix(ci.week[,1],4,4)
  upper.week <- matrix(ci.week[,2],4,4)
  
  test <- melt(lower.week) %>%
    mutate(p.value = ifelse(value < 0, 1, 0)) %>%
    select(Var1, Var2, p.value)
  
  p.mat <- diag(4)
  p.mat[cbind(test$Var1, test$Var2)] <- p.mat[cbind(test$Var2, test$Var1)] <- test$p.value
  print(p.mat)
  print(id.matrix.week)
}

names(indv.com)

#'dist_mean',
for (depVar in c('dist_mean','angleC_mean','angMC_mean','pDist_mean','pDistC_mean','velC_mean')) {
#for (depVar in c('dist_mean','polarity_mean','angleC_mean','rotation_mean','angMC_mean','pDistC_mean','velC_mean','vel_mean','pDist_mean')) {
#for (depVar in c('angMC_mean','pDistC_mean','velC_mean','vel_mean','pDist_mean')) {
#for (depVar in c('dist_mean','vel_mean','pDist_mean')) {
#for (depVar in c('pDist_mean')) {
  print(depVar)
  hourly_dep <- paste(depVar,'_',sep="")
  print("DIC")
  #func.DIC(depVar,'ExpDay',indv.com)

  print("Day1 Rpt")
  #func.hourly.repeat(hourly_dep,day1.com)
  
  print("Day27 Rpt")
  #func.hourly.repeat(hourly_dep,day27.com)
  
  print("Week:Week predictability")
  #func.weekly.predict(depVar,indv.com)
  
  print("Overall rpt")
  #func.overall.repeat(depVar,indv.com)
  
  print("Intercepts:")
  plots.weekly <- func.weekly.intercepts(depVar,indv.long)
  print(plots.weekly[1])
  print(plots.weekly[2])
  
  #plots.daily <- func.triday.intercepts(depVar,indv.com)
  print(plots.daily[1])
  print(plots.daily[2])

  ggsave(paste("./figs/plot",depVar,'WeeklyBLUP.jpg',sep='.'),plots.weekly[[1]])
  ggsave(paste("./figs/plot",depVar,'WeeklyRPT.jpg',sep='.'),plots.weekly[[2]])
  ggsave(paste("./figs/plot",depVar,'DailyBLUP.jpg',sep='.'),plots.daily[[1]])
  ggsave(paste("./figs/plot",depVar,'DailyRPT.jpg',sep='.'),plots.daily[[2]])
}



plot.weekly
func.DIC(depVar,'ExpDay',indv.com)
#"dist_mean",
## Heads up, this takes forever to run
for (depVar in c('angleC_mean','head_mean','pDistC_mean','velC_mean')) {
  print(depVar)
  testplot.week <- func.weekly.intercepts(depVar,indv.com)
  print(testplot.week)
}

## What if instead of plotting the blup's, we plot the raw data
ggplot(indv.com,aes(x=ExpDay,y=dist_mean,group=Pi,color=Pi)) + geom_line()

ggplot(indv.long,aes(x=ExpDay,y=dist_mean,group=Pi,color=Pi)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "none",
        rect=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 0)) 

indv.slidingMean <- arrange(indv.long54,Pi,ExpDay) %>% dplyr::mutate(rollingDist=rollapply(dist_mean,5,mean,align='right',fill=NA),rollingVelC=rollapply(velC_mean,5,mean,align='right',fill=NA))

indv.slidingMean$rollingDistNorm <- (800 - indv.slidingMean$rollingDist) / 800

ggplot(indv.slidingMean,aes(x=ExpDay,y=rollingDistNorm,group=Pi,color=Pi)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "none",
        rect=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 0)) 

ggplot(indv.slidingMean,aes(x=ExpDay,y=rollingVelC,group=Pi,color=Pi)) + 
  geom_line() + 
  theme_classic() + 
  theme(legend.position = "none",
        rect=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 0)) 

## Make a simple plot to show DIC numbers:
names <- c("Fixed","Random Intercept","Random Slope")
df.dic <- data.frame(
  score = c(4053.861, 3871.423, 3801.22),
  name = factor(names, levels = names),
  y = seq(length(names)) * 0.9
)
ggplot(df.dic) + geom_col(aes(score, name),width = 0.6) + 
  scale_y_discrete(limits=rev) + 
  coord_cartesian(xlim=c(3700,4200)) + 
  theme_classic() + 
  theme(legend.position = "none",
        rect=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, face = "bold", colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 0)) 
#scale_x_continuous(breaks=seq(3500,4200,length=8),limits=c(3700,4200))


