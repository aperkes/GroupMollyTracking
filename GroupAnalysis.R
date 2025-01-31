## Analysis of the tracking data for Jalle Fish
library("lme4")
library("interactions")
library("ggplot2")
library('dplyr')
data <- read.csv("~/Documents/Scripts/Sandbox/JolleTrack/JalleTracksAll.csv")

### Preliminary DataViz ####
colnames(data)
hist(data$ExpDay)
data$Pi[data$ExpDay == 0]
data$Pi[data$ExpDay == 20]
data$ExpDay[data$Pi == "pi41"]
## Some of the pi's fall off over time (e.g., pi41 drops off at day 11)
data %>%
  group_by(Pi) %>%
  summarize(max = max(ExpDay))

## Only look at 4 weeks for fish where we have that
good_data <- data[data$Pi != "pi34",]
good_data <- good_data[good_data$Pi != "pi41",]
good_data <- good_data[good_data$ExpDay <= 31,]

good_data$vel0 <- good_data$vel_mean[good_data$ExpDay == 0]
good_data %>%
  group_by(Pi) %>%
  summarize(max = max(ExpDay))

data <- good_data


# Get day 1 data score for each individual
day_one_data <- data %>% 
  filter(ExpDay == 0) %>%  # Filter data for day 1 only
  group_by(Pi) %>%  # Group by Individual
  summarize(vel0 = first(vel_mean),
            pdistC0 = first(pDistC_mean),
            dist0 = first(dist_mean))  # Get first confidence value (day 1)

# Merge day 1 vel with original data frame
data <- merge(data, day_one_data, by = "Pi", all.x = TRUE)
data1 <- data[data$ExpDay > 0,]

## Make Slices
data$ExpDay0 <- data$ExpDay
data$ExpDay3 <- data$ExpDay - 3
data$ExpDay6 <- data$ExpDay - 6
data$ExpDay9 <- data$ExpDay - 9
data$ExpDay12 <- data$ExpDay - 12
data$ExpDay15 <- data$ExpDay - 15
data$ExpDay18 <- data$ExpDay - 18
data$ExpDay21 <- data$ExpDay - 21
data$ExpDay24 <- data$ExpDay - 24
data$ExpDay27 <- data$ExpDay - 27
data$ExpDay30 <- data$ExpDay - 30
## distance to wall is "Significantly Different" probably, if we treat our sample size 
## As 40,000, but is it meaningfully different across groups? 
ggplot(data[data$ExpDay == 0,]) +
  geom_bar( aes(x=Pi, y=pDist_mean), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=Pi, ymin=pDist_mean-pDist_std, ymax=pDist_mean+pDist_std), width=0.4, colour="orange", alpha=0.9, size=1.3)

ggplot(data[data$ExpDay == 31,]) +
  geom_bar( aes(x=Pi, y=pDist_mean), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=Pi, ymin=pDist_mean-pDist_std, ymax=pDist_mean+pDist_std), width=0.4, colour="orange", alpha=0.9, size=1.3)

### General Behavior ####

## Velocity
#model <- lm(vel_mean ~ ExpDay*Pi,data=data)
model <- lmer(vel_mean ~ ExpDay + (1+ExpDay|Pi),data=data)
summary(model)
plot(fitted(model),resid(model))
abline(0,0)

qqnorm(resid(model))
qqline(resid(model))

fitted_data <- data.frame(Pi = data$Pi, ExpDay = data$ExpDay, pDist_mean = data$pDist_mean, Fitted = fitted(model))
ggplot(fitted_data, aes(ExpDay,Fitted,group=Pi,colour=Pi)) +
  geom_line()

## Distance to Wall
model <- lm(pDist_mean ~ ExpDay*Pi,data=data)
summary(model)
plot(fitted(model),resid(model))
abline(0,0)

qqnorm(resid(model))
qqline(resid(model))

fitted_data <- data.frame(Pi = data$Pi, ExpDay = data$ExpDay, pDist_mean = data$pDist_mean, Fitted = fitted(model))
ggplot(fitted_data, aes(ExpDay,Fitted,group=Pi,colour=Pi)) +
  geom_line()

### Test for behavioral divergence ####

## Distance from eachother 
#model <- lm(dist_mean ~ ExpDay*Pi,data=data)
model <- lmer(dist_mean ~ ExpDay + (1+ExpDay|Pi),data=data)

summary(model)

plot(fitted(model),resid(model))
abline(0,0)

qqnorm(resid(model))
qqline(resid(model))

fitted_data <- data.frame(Pi = data$Pi, ExpDay = data$ExpDay, pDist_mean = data$pDist_mean, Fitted = fitted(model))
ggplot(fitted_data, aes(ExpDay,Fitted,group=Pi,colour=Pi)) +
  geom_line()

## Cohesion of mean distance varies by day*group (generally increasing by day)

model <- lm(pDistC_mean ~ ExpDay*Pi,data=data)
model <- lmer(pDistC_mean ~ ExpDay + (1+ExpDay|Pi),data=data)

summary(model)

plot(fitted(model),resid(model))
abline(0,0)

qqnorm(resid(model))
qqline(resid(model))

fitted_data <- data.frame(Pi = data$Pi, ExpDay = data$ExpDay, pDist_mean = data$pDist_mean, Fitted = fitted(model))
ggplot(fitted_data, aes(ExpDay,Fitted,group=Pi,colour=Pi)) +
  geom_line()

## Cohesion of mean distance varies by day*group (generally increasing by day)
model <- lm(pDistC_mean ~ ExpDay*Pi,data=data)
summary(model)

plot(fitted(model),resid(model))
abline(0,0)

qqnorm(resid(model))
qqline(resid(model))

fitted_data <- data.frame(Pi = data$Pi, ExpDay = data$ExpDay, pDist_mean = data$pDist_mean, Fitted = fitted(model))
ggplot(fitted_data, aes(ExpDay,Fitted,group=Pi,colour=Pi)) +
  geom_line()

interact_plot(model, pred=ExpDay,modx=Pi,modxvals=c("pi12",'pi13','pi14','pi31','pi32','pi33'))

## Velocity cohesion
model <- lm(velC_mean ~ ExpDay*Pi,data=data)
model <- lmer(velC_mean ~ ExpDay + (1+ExpDay|Pi),data=data)
model3 <- lmer(velC_mean ~ ExpDay3 + (1+ExpDay3|Pi),data=data)
model6 <- lmer(velC_mean ~ ExpDay6 + (1+ExpDay6|Pi),data=data)
model9 <- lmer(velC_mean ~ ExpDay9 + (1+ExpDay9|Pi),data=data)

model27 <- lmer(velC_mean ~ ExpDay27 + (1+ExpDay27|Pi),data=data)
summary(model3)

plot(fitted(model),resid(model))
abline(0,0)

qqnorm(resid(model))
qqline(resid(model))

fitted_data <- data.frame(Pi = data$Pi, ExpDay = data$ExpDay3, pDist_mean = data$pDist_mean, Fitted = fitted(model3))
ggplot(fitted_data, aes(ExpDay,Fitted,group=Pi,colour=Pi)) +
  geom_line()

interact_plot(model, pred=ExpDay,modx=Pi,modxvals=c("pi12",'pi13','pi14','pi31','pi32','pi33'))

model <- lm(head_mean ~ ExpDay*Pi,data=data)
summary(model)

plot(fitted(model),resid(model))
abline(0,0)

qqnorm(resid(model))
qqline(resid(model))

fitted_data <- data.frame(Pi = data$Pi, ExpDay = data$ExpDay, pDist_mean = data$pDist_mean, Fitted = fitted(model))
ggplot(fitted_data, aes(ExpDay,Fitted,group=Pi,colour=Pi)) +
  geom_line()

interact_plot(model, pred=ExpDay,modx=Pi,modxvals=c("pi12",'pi13','pi14','pi31','pi32','pi33'))


### Day 0 differences predict later behavior? ####

## Does day0 velocity predict velocity? Probably yes? 

model <- lmer(vel_mean ~ vel0 + ExpDay + (1|Pi),data = data1)
summary(model)

plot(fitted(model),resid(model))
abline(0,0)

qqnorm(resid(model))
qqline(resid(model))

## Does day0 Cohesion predict later cohesion? Probably not
model <- lmer(pDistC_mean ~ pdistC0 + ExpDay + (1|Pi),data = data1)
summary(model)

plot(fitted(model),resid(model))
abline(0,0)

qqnorm(resid(model))
qqline(resid(model))

## Does day0 distance predict later distance? No
model <- lmer(dist_mean ~ dist0 + ExpDay + (1|Pi),data = data1)
summary(model)

plot(fitted(model),resid(model))
abline(0,0)

qqnorm(resid(model))
qqline(resid(model))
