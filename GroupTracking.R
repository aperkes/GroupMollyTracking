### Organized according to the paper: 

library(nlme)
library(ggplot2)
library(dplyr)
library(zoo)
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
library(matrixStats)

setwd("~/Documents/Scripts/GroupMollyTracking/")

indv <- read.csv("JolleTracksAll_4.csv")
indv <- read.csv("JolleTracksAll_5.csv")
indv$ExpDay <- as.integer(indv$ExpDay)

indv[indv$Pi == 'pi11',]$ExpDay <- indv[indv$Pi == 'pi11',]$ExpDay - 16
#indv[indv$Pi == 'pi12',]$ExpDay <- indv[indv$Pi == 'pi12',]$ExpDay + 15

indv$dist_meanScale <- scale(indv$dist_mean)
indv$vel_meanScale <- scale(indv$vel_mean)
indv$pDist_meanScale <- scale(indv$pDist_mean)

indv$velC_meanScale <- scale(indv$velC_mean)
indv$pDistC_meanScale <- scale(indv$pDistC_mean)
indv$angleC_meanScale <- scale(indv$angleC_mean)


indv$week <- indv$ExpDay %/% 7 ## Is this integer division
indv$triday <- indv$ExpDay %/% 3

indv.good <- indv[indv$ExpDay >= 0,]
#indv.good <- indv.good[indv.good$Pi != 'pi12',]

indv <- indv.good



## Some of the pi's fall off over time (e.g., pi41 drops off at day 11)
## Others go on 90 days
piDays <- indv %>%
  group_by(Pi) %>%
  summarize(max = max(ExpDay))

piDays[piDays$max >= 54,]
#long_data <- indv %>%
#  filter(Pi %in% piDays[piDays$max >= 62,]$Pi)

#indv.long <- long_data[long_data$ExpDay <= 63,]

long_data54 <- indv %>%
  filter(Pi %in% piDays[piDays$max >= 54,]$Pi)
indv.long54 <- long_data54[long_data54$ExpDay <= 54,]

## Only look at 4 weeks for fish where we have that
good_data <- indv[indv$Pi != "pi34",] ## This one is missing 
good_data <- good_data[good_data$Pi != "pi41",] ## This oen is zoomed out, maybe impossible
#good_data <- good_data[good_data$ExpDay <= 31,] ## This one is missing one day

indv.com <- good_data

### Same for hourly (not sure I need this code...)

hourly <- read.csv("JolleTracksHourly_2.csv")
hourly$Hour <- as.integer(hourly$Hour)
hourly$ExpDay <- as.integer(hourly$ExpDay)

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

### DATA EXPLORATION GOES HERE


## How do behaviors vary over time? 
res.cohXtime <- lme(dist_mean ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
summary(res.cohXtime)

res.velXtime <- lme(vel_mean ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
summary(res.velXtime)
#vel.resid <- as.double(VarCorr(res.velXtime)["Residual","Variance"])

res.pDistXtime <- lme(pDist_mean ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
summary(res.pDistXtime)
#pDist.resid <- as.double(VarCorr(res.pDistXtime)["Residual","Variance"])

res.pDistCXtime <- lme(pDistC_mean ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
summary(res.pDistCXtime)
#pDistC.resid <- as.double(VarCorr(res.pDistCXtime)["Residual","Variance"])

res.activityXtime <- lme(prop_active ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
summary(res.activityXtime)
#activity.resid <- as.double(VarCorr(res.activityXtime)["Residual","Variance"])

res.uppervelXtime <- lme(upper_vel ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
summary(res.uppervelXtime)
#uppervel.resid <- as.double(VarCorr(res.uppervelXtime)["Residual","Variance"])
### Add max velocity and activity time

### Need to still control for fish size (but it's decreasing, so as a function of fish size will just increase it)

### Velocity cohesion doesn't change, just hangs out around 0.03, low but sig
res.velCXtime <- lme(velC_mean ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
summary(res.velCXtime)
#velC.resid <- as.double(VarCorr(res.velCXtime)["Residual","Variance"])

res.angleCXtime <- lme(angleC_mean ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
summary(res.angleCXtime)
#angleC.resid <- as.double(VarCorr(res.angleCXtime)["Residual","Variance"])


### Plotting group blups over time

# Need to build long df
long_data54 <- indv %>%
  filter(Pi %in% piDays[piDays$max >= 54,]$Pi)
indv.long54 <- long_data54[long_data54$ExpDay <= 54,]

## Define priors used throughout for MCMCglm


prior.id.slope <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1=list(V = diag(2), nu = 0.002, alpha.mu = c(0,0),  alpha.V = diag(2)*25^2)))

### Define prior for heterogeneous residual variance
total_days <- nlevels(as.factor(indv.long54$ExpDay))

## Priors.cov need to be behavior specific because the prior value actually matters
prior.id.slope.cov <- list(R = list(V = diag(total_days),nu = total_days + 0.002),
                           G = list(G1=list(V = diag(2), nu = 0.002, alpha.mu = c(0,0),  alpha.V = diag(2)*25^2)))

prior.best <- list(R = list(V = diag(total_days)*0.5,nu = total_days + 0.002),
                   G = list(G1=list(V = diag(2)*0.5, nu = 2.002, alpha.mu = c(0,0),  alpha.V = diag(2)*25^2)))

### Define function to get repeatability ####

n_days <- 1

func.ndays.intercepts.het <- function(depVar,df,day_bin=1,prior.cov = prior.id.best) {
  # select only those IDs in that vector & only keep up to obs 70 & make obs 1 = 0

  max_days <- max(df$ExpDay) ## It's 0 indexed
  indv.df <- df
  
  n_pis <- length(unique(df$Pi))
  rpt <- list()
  ci.rpt <- list()
  
  post.among <- list()
  ci.among <- list()
  post.within <- list()
  ci.within <- list()
  
  #blup <- list()
  ### Define fixed and random formulas for given dep variable
  fixed <- as.formula(paste(depVar,"~ ExpDay",sep=""))
  random <- as.formula(paste('~us(1 + ExpDay):Pi',sep=''))
  
  ## Need day as a factor for the rcov line
  indv.df$ExpDay_factor <- as.factor(indv.df$ExpDay)
  model.het <- MCMCglmm(fixed = fixed, 
                        random = random, 
                        rcov = ~idh(ExpDay_factor):units, ## use this line for het. residual variance
                        data = indv.df, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.cov, ## Replace this with prior.id.slope for homo residual variance
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)  
  
  sigma.a0 <- model.het$VCV[,"(Intercept):(Intercept).Pi"]
  sigma.a1 <- model.het$VCV[,"ExpDay:ExpDay.Pi"]
  
  rho <- model.het$VCV[,"(Intercept):ExpDay.Pi"] ### whole covariance term, not just rho
  
  ### Calculate conditional rpt at each point x
  for (x in seq(0,max_days,day_bin)) {
    ## Grab residual at that day
    e_col <- paste("ExpDay_factor",x,".units",sep='')
    sigma.e.x <- model.het$VCV[,e_col]
    
    sigma.among.x <- ( sigma.a0 + sigma.a1*(x**2) + 2*rho*x )
    rpt.x <- sigma.among.x / 
      ( sigma.among.x + sigma.e.x)
    
    rpt <- c(rpt,median(rpt.x))
    ci.x <- HPDinterval(rpt.x)[1:2]
    ci.rpt <- c(ci.rpt,ci.x)
    
    post.among.x <- median(sigma.among.x) ## 
    post.among <- c(post.among,post.among.x)
    
    ci.among.x <- HPDinterval(sigma.among.x)[1:2]
    ci.among <- c(ci.among,ci.among.x)
    
    post.within.x <- median(sigma.e.x)
    post.within <- c(post.within,post.within.x)
    
    ci.within.x <- HPDinterval(sigma.e.x)[1:2]
    ci.within <- c(ci.within,ci.within.x)
    
    'this pulls out the individual intercepts and adds in the overall intercepts
    so that way these numbers are absolute values, as opposed to differences from overall'

    #blup <- c(blup,intercepts.n)
  }
  intercepts <- unname(median(model.het$Sol)[3:(3+n_pis-1)] + median(model.het$Sol)["(Intercept)"])
  #blup <- unlist(blup)
  rpt <- unlist(rpt)
  dates <- seq(0,54,day_bin) ## 
  
  ## Switches back to old syntax here
  ci.rpt <- matrix(unlist(ci.rpt), nrow = length(dates), byrow = T)
  ci.id <- matrix(unlist(ci.among), nrow = length(dates), byrow = T)
  ci.w <- matrix(unlist(ci.within), nrow = length(dates), byrow = T)      
  post.id <- unlist(post.among)
  post.w <- unlist(post.within)

  
  rpt.slice.wide <- data.frame(dates, "rpt" = rpt, "lower.rpt" = ci.rpt[,1], "upper.rpt" = ci.rpt[,2],
                               "post.id" = post.id, "lower.id" = ci.id[,1], "upper.id" = ci.id[,2], 
                               "post.w" = post.w, "lower.w" = ci.w[,1], "upper.w" = ci.w[,2])
  
  rpt.slice.long <- data.frame("date" = rep(dates, 3), 
                               "type" = rep(c("rpt", "id", "within"), each = length(dates)),
                               "variance" = unname(c(rpt, post.id, post.w)),
                               "lower" = unname(c(ci.rpt[,1], ci.id[,1], ci.w[,1])),
                               "upper" = unname(c(ci.rpt[,2], ci.id[,2], ci.w[,2])))
  
  
  rpt.slice.rpt <- rpt.slice.long[rpt.slice.long$type == 'rpt',]
  
  #n_pis <- length(colnames(col.day0$Sol)) / 2 - 1
  ### This is very dataset specific, you may need to check it: 
  ids <- colnames(model.het$Sol)[3:(3 + n_pis - 1)] ## Make sure this matches below
  ids <- substr(ids, 16, 19)
  
  dates.rep <- rep(dates, each = n_pis)
  picomp <- rep(ids, length(dates))
  
  pred.intercepts <- data.frame(dates.rep, picomp, intercepts)
  
  plt.repeat <- 0
  
  pred.rank <- pred.intercepts %>%
    filter(dates.rep == 0) %>%
    mutate(ranking = rank(intercepts)) %>%
    select(picomp, ranking)
  
  pred.rank <- left_join(pred.intercepts, pred.rank) %>%
    arrange(ranking)
  plasma_pal <- viridis::plasma(n = 30)
  plasma_pal <- plasma_pal[1:26]
  
  
  blup.plot <- 0
  rpt.plot <- 0
  
  plt.repeat <- 0
  plt.day <- 0
  return(list(plt.day,plt.repeat,rpt.plot,blup.plot,rpt.slice.wide,pred.rank,model.het))
}

### Get repeatability values ####

## These two are in the paper proper 

plots.dist <- func.ndays.intercepts.het('dist_meanScale',indv.long54,n_days,prior.cov=prior.best)
plots.velC <- func.ndays.intercepts.het('velC_meanScale',indv.long54,n_days,prior.cov=prior.best)

## These Four are supplemental
plots.vel <- func.ndays.intercepts.het('vel_meanScale',indv.long54,n_days,prior.cov=prior.best)
plots.wall_dist <- func.ndays.intercepts.het('pDist_meanScale',indv.long54,n_days,prior.cov=prior.best)
#plots.vel_std <- func.ndays.intercepts('vel_std',df,n_days)
#plots.velC_scale <- func.ndays.intercepts('velC_scale',df,n_days)
plots.wall_distC <- func.ndays.intercepts.het('pDistC_meanScale',indv.long54,n_days,prior.cov=prior.best)
plots.angleC <- func.ndays.intercepts.het('angleC_meanScale',indv.long54,n_days,prior.cov=prior.best)
#plots.angle_std <- func.ndays.intercepts('polarity_std',df,n_days)

### Sliding mean plots ####
### Blups are fancy, but sliding average is much easier to explain:
### Also a lot faster
func.slidingMean <- function(depVar,df,n_days) {
  indv.slidingMean <- arrange(df,Pi,ExpDay) %>%
    dplyr::mutate(rollingVar=rollapply(!! rlang::ensym(depVar),n_days,mean,align='left',fill=NA))
  last_day <- max(indv.slidingMean$ExpDay) - 5
  indv.slidingMean.clipped <- indv.slidingMean[indv.slidingMean$ExpDay <= last_day,]
  plot.sliding <- ggplot(indv.slidingMean.clipped,aes(x=ExpDay,y=rollingVar,group=Pi,color=Pi)) + 
    geom_line() + 
    ylab(paste('Rolling',depVar)) +
    theme_classic() + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 0)) 
  plot.raw <- ggplot(indv.slidingMean,aes(x=ExpDay,y=!! rlang::ensym(depVar),group=Pi,color=Pi)) + 
    geom_line() + 
    ylab(paste('Rolling',depVar)) +
    theme_classic() + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 0)) 
  return(list(plot.sliding,indv.slidingMean,plot.raw))
}

sliding.dist <- func.slidingMean('dist_mean',indv.long54,5)
sliding.velC <- func.slidingMean('velC_mean',indv.long54,5)
sliding.velCScale <- func.slidingMean('velC_meanScale',indv.long54,5)

sliding.vel <- func.slidingMean('vel_mean',indv.long54,5)
#sliding.vel_std <- func.slidingMean('vel_std',indv.long54,5)
#sliding.velC_scale <- func.slidingMean('velC_scale',indv.long54,5)
sliding.wall_dist <- func.slidingMean('pDist_mean',indv.long54,5)
sliding.wall_distC <- func.slidingMean('pDistC_mean',indv.long54,5)
sliding.angleC <- func.slidingMean('angleC_mean',indv.long54,5)
#sliding.angle_std <- func.slidingMean('polarity_std',indv.long54,5)

### Make plot of repeatability over time ####

## rpt.data is calculated as part of intercepts, so you can get it as plots.foo[[5]]
func.rpt.plot <- function(rpt.data,n_days=1,components=F) {
  dates <- seq(0,54,n_days)
  print(rpt.data)
  print(dates)
  bar_scale <- max(rpt.data$post.id) + max(rpt.data$post.w)
  #foo <- plots.vel[[5]]
  rpt.data$lower.id_ <- rpt.data$lower.id / bar_scale
  rpt.data$upper.id_ <- rpt.data$upper.id / bar_scale
  rpt.data$post.id_ <- rpt.data$post.id / bar_scale
  rpt.data$lower.w_ <- rpt.data$lower.w / bar_scale
  rpt.data$upper.w_ <- rpt.data$upper.w / bar_scale
  rpt.data$post.w_ <- rpt.data$post.w / bar_scale
  
  n_dates = length(dates)
  if (n_dates > 10) { 
    breaks <- seq(0,max(dates),5)
  }
  else {
    breaks <- dates
  }  
  if (components == T) {
  rpt.plot <- ggplot(rpt.data, aes(x = dates)) +
    geom_line(aes(x = dates-0.5, y = post.id_, color = "#959595")) +
    geom_point(aes(x = dates-0.5, y = post.id_), shape = 21, color = "#000000", fill = "#959595", size = 1) +
    
    #geom_errorbar(aes(x = dates+0.5, ymin = lower.w_, ymax = upper.w_, width = 1, color = "#000000")) +
    geom_line(aes(x = dates+0.5, y = post.w_, color = "#CCCCCC")) +
    geom_point(aes(x = dates + 0.5, y = post.w_), shape = 21, color = "#000000", fill = "#CCCCCC", size = 1) +
  
    geom_point(aes(y = rpt, color = "#000000"), size = 1) +
    geom_line(aes(x=dates, y = rpt, color = "#000000")) +
    geom_errorbar(aes(ymin = lower.rpt, ymax = upper.rpt, width = 1, color = "#000000")) 
  }
  else {
    rpt.plot <- ggplot(rpt.data, aes(x = dates)) +
    geom_point(aes(y = rpt, color = "#000000"), size = 1) +
    geom_line(aes(x=dates, y = rpt, color = "#000000")) +
    geom_errorbar(aes(ymin = lower.rpt, ymax = upper.rpt, width = 1, color = "#000000")) 
  }
    rpt.plot <- rpt.plot + 
    scale_x_continuous(breaks = breaks, labels = breaks) +
    scale_y_continuous(name = "Variance estimate",limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    labs(x = "Day",
         color = "Legend") +
    scale_color_manual(name = "", 
                       values =c("#000000",  "#959595",  "#CCCCCC"),
                       labels = c("Repeatability", "Among-Group", "Within-Group")) +
    theme_classic() +
    theme(legend.position = c(0.2, 0.93),
          legend.key.size = unit(0.3, 'cm'),
          legend.text = element_text(size = 6),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 10),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 10)) #+
  
  return(rpt.plot)
}

### Calculate Repeatability (these plots are down below)

## Again first two are for the main body
#rpt.plot.dist.old <- func.rpt.plot(plots.dist[[5]])
rpt.plot.dist <- func.rpt.plot(plots.dist[[5]])
rpt.plot.dist2 <- func.rpt.plot(plots.dist[[5]],components=T)
rpt.plot.velC <- func.rpt.plot(plots.velC[[5]])
rpt.plot.velC2 <- func.rpt.plot(plots.velC[[5]],components=T)

rpt.plot.vel <- func.rpt.plot(plots.vel[[5]],components=T)
#rpt.plot.vel_std <- func.rpt.plot(plots.vel_std[[5]])
rpt.plot.wall_dist <- func.rpt.plot(plots.wall_dist[[5]],components=T)
rpt.plot.wall_distC <- func.rpt.plot(plots.wall_distC[[5]],components=T)
rpt.plot.angleC <- func.rpt.plot(plots.angleC[[5]],components=T)
#rpt.plot.angle_std <- func.rpt.plot(plots.angle_std[[5]])

#### Run permutation analysis ####

## Function to make a shuffled dataframe
func.shuffle.df <- function(df) {
  
  df.shuffle <- df
  df.shuffle$ShuffledPi <- 0
  
  for (i in seq(0,max(df.shuffle$ExpDay))) { 
    df.shuffle[df.shuffle$ExpDay == i,'ShuffledPi'] <- sample(df.shuffle[df.shuffle$ExpDay == i,'Pi'])
  }
  
  df.shuffle$ShuffledPi
  df.shuffle$Pi <- df.shuffle$ShuffledPi
  
  return(df.shuffle)
}
### Simplified function to just get repeatability, to be run in parallel

### Simplified function for shuffling and calculating rpt. 
func.shuffled.rpt.i <- function(i,depVar="dist_mean",day_bin=1,prior.cov=prior.best) {
  # select only those IDs in that vector & only keep up to obs 70 & make obs 1 = 0
  set.seed(i)
  df.shuffled <- func.shuffle.df(indv.long54)
  #depVar <- "dist_mean"
  #day_bin <- 1
  max_days <- 54
  
  n_pis <- length(unique(df.shuffled$Pi))
  rpt <- list()
  ci.rpt <- list()
  post.id <- list()
  ci.id <- list()
  post.w <- list()
  ci.w <- list()
  blup <- list()
  
  ### Define fixed and random formulas for given dep variable
  fixed <- as.formula(paste(depVar,"~ ExpDay",sep=""))
  random <- as.formula(paste('~us(1 + ExpDay):Pi',sep=''))
  
  ## Need day as a factor for the rcov line
  df.shuffled$ExpDay_factor <- as.factor(df.shuffled$ExpDay)
  model.het <- MCMCglmm(fixed = fixed, 
                        random = random, 
                        rcov = ~idh(ExpDay_factor):units, ## use this line for het. residual variance
                        data = df.shuffled, 
                        family = "gaussian",
                        pr = T,
                        prior = prior.cov,
                        nitt=310000, burnin = 10000, thin = 200, 
                        verbose = F)  
  
  sigma.a0 <- model.het$VCV[,"(Intercept):(Intercept).Pi"]
  
  sigma.a1 <- model.het$VCV[,"ExpDay:ExpDay.Pi"]
  
  ## Using fixed individual variance
  #sigma.e <- model.het$VCV[,"units"]
  
  rho <- model.het$VCV[,"(Intercept):ExpDay.Pi"] ### this is actually the whole covariance term, not just rho
  
  ### Calculate conditional rpt at each point x
  for (x in seq(0,max_days,day_bin)) {
    ## Grab residual at that day
    e_col <- paste("ExpDay_factor",x,".units",sep='')
    sigma.e <- model.het$VCV[,e_col]
    
    rpt.x <- ( sigma.a0 + sigma.a1*(x**2) + 2*rho*x ) / 
      ( sigma.a0 + sigma.a1*(x**2) + 2*rho*x + sigma.e)
    
    rpt <- c(rpt,median(rpt.x)) # median is preferred to mode
  }
  rpt <- unlist(rpt)
  return(rpt) }

func.shuffled.rpt.i(1)
## Run for real: 
### Set interations:  vvv below. Be sure to set cores  vv. 
rpt.all <- mclapply(1:50,func.shuffled.rpt.i,mc.cores=5L)
rpt.mat.all <- do.call("cbind",rpt.all)
rpt.quants <- rowQuantiles(rpt.mat.all,probs=c(.025,.975))

rpt.coh <- mclapply(1:50,func.shuffled.rpt.i,depVar='dist_meanScale',prior.cov=prior.best,mc.cores=5L)
rpt.mat.coh <- do.call("cbind",rpt.coh)
rpt.quants.coh <- rowQuantiles(rpt.mat.coh,probs=c(.025,.975))

rpt.vel <- mclapply(1:50,func.shuffled.rpt.i,depVar='vel_meanScale',prior.cov=prior.best,mc.cores=5L)
rpt.mat.vel <- do.call("cbind",rpt.vel)
rpt.quants.vel <- rowQuantiles(rpt.mat.vel,probs=c(.025,.975))

rpt.velC <- mclapply(1:50,func.shuffled.rpt.i,depVar='velC_meanScale',prior.cov=prior.best,mc.cores=5L)
rpt.mat.velC <- do.call("cbind",rpt.velC)
rpt.quants.velC <- rowQuantiles(rpt.mat.velC,probs=c(.025,.975))

rpt.pdist <- mclapply(1:50,func.shuffled.rpt.i,depVar='pDist_meanScale',prior.cov=prior.best,mc.cores=5L)
rpt.mat.pdist <- do.call("cbind",rpt.pdist)
rpt.quants.pdist <- rowQuantiles(rpt.mat.pdist,probs=c(.025,.975))

rpt.pdistC <- mclapply(1:50,func.shuffled.rpt.i,depVar='pDistC_meanScale',prior.cov=prior.best,mc.cores=5L)
rpt.mat.pdistC <- do.call("cbind",rpt.pdistC)
rpt.quants.pdistC <- rowQuantiles(rpt.mat.pdistC,probs=c(.025,.975))

rpt.angleC <- mclapply(1:50,func.shuffled.rpt.i,depVar='angleC_meanScale',prior.cov=prior.best,mc.cores=5L)
rpt.mat.angleC <- do.call("cbind",rpt.angleC)
rpt.quants.angleC <- rowQuantiles(rpt.mat.angleC,probs=c(.025,.975))


### Make Plots ####

## First two rows of plots: 
#rpt.plot.dist; sliding.dist[[1]]
#rpt.plot.velC; sliding.velC[[1]]

## Main Plots
rpt.plot.dist_ <- rpt.plot.dist + geom_ribbon(aes(ymin = rpt.quants.coh[,1], ymax = rpt.quants.coh[,2]), fill = "lightblue", alpha = 0.5)
sliding.dist[[1]]; rpt.plot.dist_
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/plots.dist_mean.RPT.jpg',rpt.plot.dist_,width = 3,height=3,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/plots.dist_mean.sliding.jpg',sliding.dist[[1]],width = 3,height=3,units="in")

rpt.plot.velC_ <- rpt.plot.velC + geom_ribbon(aes(ymin = rpt.quants.velC[,1], ymax = rpt.quants.velC[,2]), fill = "lightblue", alpha = 0.5)
sliding.velC[[1]]; rpt.plot.velC_

ggsave('~/Documents/Scripts/GroupMollyTracking/figs/plots.velC.RPT.jpg',rpt.plot.velC_,width = 3,height=3,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/plots.velC.sliding.jpg',sliding.velC[[1]],width = 3,height=3,units="in")

### Supplemental Plots
rpt.plot.dist2_ <- rpt.plot.dist2 + geom_ribbon(aes(ymin = rpt.quants.coh[,1], ymax = rpt.quants.coh[,2]), fill = "lightblue", alpha = 0.5)
rpt.plot.velC2_ <- rpt.plot.velC2 + geom_ribbon(aes(ymin = rpt.quants.velC[,1], ymax = rpt.quants.velC[,2]), fill = "lightblue", alpha = 0.5)

ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.dist_mean.RPT.jpg',rpt.plot.dist2_,width = 3,height=3,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.velC.RPT.jpg',rpt.plot.velC2_,width = 3,height=3,units="in")

rpt.plot.vel_ <- rpt.plot.vel + geom_ribbon(aes(ymin = rpt.quants.vel[,1], ymax = rpt.quants.vel[,2]), fill = "lightblue", alpha = 0.5)
sliding.vel[[1]];rpt.plot.vel_
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.vel.RPT.jpg',rpt.plot.vel_,width = 3,height=3,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.vel.sliding.jpg',sliding.vel[[1]],width = 3,height=3,units="in")

rpt.plot.wall_dist_ <- rpt.plot.wall_dist + geom_ribbon(aes(ymin = rpt.quants.pdist[,1], ymax = rpt.quants.pdist[,2]), fill = "lightblue", alpha = 0.5)
sliding.wall_dist[[1]];rpt.plot.wall_dist_
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.wall_dist.RPT.jpg',rpt.plot.wall_dist_,width = 3,height=3,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.wall_dist.sliding.jpg',sliding.wall_dist[[1]],width = 3,height=3,units="in")

rpt.plot.wall_distC_ <- rpt.plot.wall_distC + geom_ribbon(aes(ymin = rpt.quants.pdistC[,1], ymax = rpt.quants.pdistC[,2]), fill = "lightblue", alpha = 0.5)
sliding.wall_distC[[1]];rpt.plot.wall_distC_
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.wall_distC.RPT.jpg',rpt.plot.wall_distC_,width = 3,height=3,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.wall_distC.sliding.jpg',sliding.wall_distC[[1]],width = 3,height=3,units="in")


rpt.plot.angleC_ <- rpt.plot.angleC + geom_ribbon(aes(ymin = rpt.quants.angleC[,1], ymax = rpt.quants.angleC[,2]), fill = "lightblue", alpha = 0.5)
sliding.angleC[[1]]; rpt.plot.angleC_
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.angleC.RPT.jpg',rpt.plot.angleC_,width = 3,height=3,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.angleC.sliding.jpg',sliding.angleC[[1]],width = 3,height=3,units="in")

## There are also some plots generated in python. 

## Supplemental plots: 
#rpt.plot.vel; sliding.vel[[1]]
#rpt.plot.wall_dist; sliding.wall_dist[[1]]
#rpt.plot.wall_distC; sliding.wall_distC[[1]]
#rpt.plot.angleC; sliding.angleC[[1]]

### Difference vs random (row 3 in fig 2) is calculated using shuffleTracks.py

func.ndays.predict <- function(depVar,df,n_days=5) {
  #df <- indv.long54
  #df$velC_scale <- df$velC_mean * 100
  
  #df
  ## build bin column
  df$bin <- df$ExpDay %/% n_days
  print(df)
  ## Can I get away with not dropping values? 
  #df <- df[df$Pi != 'pi31',] ## This one has missing values. 
  
  maxes <- df %>%
    group_by(Pi) %>%
    summarise(max = max(bin, na.rm=TRUE))
  
  max_days = max(unique(df$ExpDay))
  max_bin = max_days %/% n_days
  ## Generating wide programatically is tricky...
  
  ## This is sort of a trick, it works for any n_days less than 11
  df.wide <- df %>%
    mutate(
      obs.day = case_when(
        ExpDay %in% seq(0,max_days,n_days) ~ 1 %% n_days + 1, 
        ExpDay %in% seq(1,max_days,n_days) ~ 2 %% n_days + 1, 
        ExpDay %in% seq(2,max_days,n_days) ~ 3 %% n_days + 1, 
        ExpDay %in% seq(3,max_days,n_days) ~ 4 %% n_days + 1, 
        ExpDay %in% seq(4,max_days,n_days) ~ 5 %% n_days + 1, 
        ExpDay %in% seq(5,max_days,n_days) ~ 6 %% n_days + 1, 
        ExpDay %in% seq(6,max_days,n_days) ~ 7 %% n_days + 1, 
        ExpDay %in% seq(7,max_days,n_days) ~ 8 %% n_days + 1, 
        ExpDay %in% seq(8,max_days,n_days) ~ 9 %% n_days + 1, 
        ExpDay %in% seq(9,max_days,n_days) ~ 10 %% n_days + 1, 
        ExpDay %in% seq(10,max_days,n_days) ~ 11 %% n_days + 1, 
      )) %>%
    select(Pi, bin, all_of(depVar), obs.day) %>%
    spread(bin, depVar, sep = "")
  
  #df.clean <- df.wide[, colSums(is.na(df.wide)) == 0]
  df.clean <- df.wide[rowSums(is.na(df.wide)) == 0,]
  
  df.wide <- df.clean ## Is this ok? 
  ## Not sure how do handle the next trick...need to cbind an arbitrary bin0
  
  bins.list <- (names(df.clean))[3:length(names(df.clean))]
  bins.str <- paste(bins.list,collapse = ',')
  
  form.bins <- as.formula(paste("cbind(",bins.str,") ~ trait - 1",sep=""))
  print(bins.list)
  n_bins <- length(bins.list)
  
  prior.b <- list(R = list(V = diag(n_bins), nu = n_bins + .002),
                  G = list(G1=list(V = diag(n_bins), nu = n_bins + .002, alpha.mu = rep(0,n_bins), alpha.V = 1000*diag(n_bins))))

  print(form.bins)
  print(bins.list)

  for (b in bins.list) {
    #col_name <- paste(b,"Scale",sep=".")
    df.wide[,b] <- scale(df.wide[,b])

  }

  print(form.bins)
  behav.bin.id <- MCMCglmm(form.bins, 
                           random = ~us(trait):Pi,
                           rcov = ~idh(trait):units,
                           family = c(rep("gaussian", n_bins)), 
                           prior = prior.b, 
                           pr = T, 
                           nitt = 510000, thin = 200, burnin = 10000, 
                           verbose = F,
                           data = df.wide)
  
  # Model with only ID ----
  
  id.matrix.bin <- matrix(colMedians(posterior.cor(behav.bin.id$VCV[,1:(n_bins * n_bins)])),n_bins,n_bins, 
                          dimnames = list(bins.list, 
                                          bins.list))

  # now to extract the CI estimates
  ci.bin <- data.frame(HPDinterval(posterior.cor(behav.bin.id$VCV[,1:(n_bins*n_bins)])))
  print(ci.bin)
  # for corrplot need 3 matrices - estimates, lower CI, upper CI
  lower.bin <- matrix(ci.bin[,1],length(bins.list),length(bins.list))
  upper.bin <- matrix(ci.bin[,2],length(bins.list),length(bins.list))
  
  test <- melt(lower.bin) %>%
    mutate(p.value = ifelse(value < 0, 1, 0)) %>%
    select(Var1, Var2, p.value)
  
  df.corrs <- melt(replace(id.matrix.bin, lower.tri(id.matrix.bin, T), NA), na.rm = T)
  str(df.corrs)
  
  df.corrs$start.bin <- as.numeric(substr(df.corrs$Var1, 4, 5))
  df.corrs$end.bin <- as.numeric(substr(df.corrs$Var2, 4, 5))
  df.corrs$diff <- df.corrs$end.bin - df.corrs$start.bin
  
  #now pull out ci for each corr
  lower.bin <- matrix(ci.bin[,1],n_bins,n_bins)
  upper.bin <- matrix(ci.bin[,2],n_bins,n_bins)
  
  test.lower <- melt(replace(lower.bin, lower.tri(lower.bin, T), NA), na.rm = T)
  test.upper <- melt(replace(upper.bin, lower.tri(upper.bin, T), NA), na.rm = T)
  
  ci.long <- left_join(test.lower, test.upper, by = c("Var1", "Var2")) %>%
    rename(start.bin = Var1,
           end.bin = Var2,
           lower = value.x, 
           upper = value.y) %>%
    arrange(start.bin, end.bin)
  
  ci.long$start.bin <- ci.long$start.bin -1
  ci.long$end.bin <- ci.long$end.bin -1
  among.corr <- left_join(df.corrs, ci.long, by = c("start.bin", "end.bin"))  
  print(among.corr)
  
  p.mat <- diag(n_bins)
  p.mat[cbind(test$Var1, test$Var2)] <- p.mat[cbind(test$Var2, test$Var1)] <- test$p.value
  print(p.mat)
  print(id.matrix.bin)
  
  plt.week.corr <- ggcorrplot(id.matrix.bin, 
                              type = "lower", 
                              p.mat = p.mat, 
                              insig= "blank",
                              colors = c("slateblue4","gray", "mediumorchid1"))
  
  plt.week.corr
  
  s=1
  
  pred.df <-data.frame(Pred=double(),upper=double(),lower=double(),bin=integer(),Pstep=factor())
  
  r.0 <- unlist(id.matrix.bin[1,])
  ci_low.0 <- unlist(lower.bin[1,])
  ci_up.0 <- unlist(upper.bin[1,])
  
  r.last <- unlist(id.matrix.bin[,n_bins])
  ci_low.last <- unlist(lower.bin[,n_bins])
  ci_up.last <- unlist(upper.bin[,n_bins])
  
  pred.df2 <- data.frame(bin=seq(n_bins),bin0_corr=r.0,bin0_lower=ci_low.0,bin0_upper=ci_up.0,
                         binLast_corr=r.last,binLast_lower=ci_low.last,binLast_upper=ci_up.last)
  
  plt.predict.one <- ggplot(pred.df2,aes(x=bin*n_days)) + 
    #geom_point(shape=1,size=4) + 
    geom_line(aes(y=bin0_corr,colour='First Bin Correlation')) +
    geom_line(aes(y=binLast_corr,colour='Last Bin Correlation')) +
    #geom_errorbar(aes(ymin=bin0_lower,ymax=bin0_upper,colour='First Bin Correlation')) + 
    #geom_errorbar(aes(ymin=binLast_lower,ymax=binLast_upper,colour='Last Bin Correlation')) + 
    #ggtitle(depVar) +
    #ylim(0,1) + 
    ylab('Group-level Correlation') + 
    xlab('Days since birth') + 
    theme_classic() + 
    scale_color_discrete(name="Step") +
    guides(colour = guide_legend(reverse = F),) + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 14)) 
  
  plt.predict.one
  
  
  for (s in 1:(n_bins -2)){
    pred.list <- list()
    lower.list <- list()
    upper.list <- list()
    for(i in 1:(n_bins-s)){
      r <- id.matrix.bin[i,i+s]
      ci_low <- lower.bin[i,i+s]
      ci_up <- upper.bin[i,i+s]
      pred.list[i] <- r
      lower.list[i] <- ci_low
      upper.list[i] <- ci_up
    }
    pred.list <- unlist(pred.list)
    lower.list <- unlist(lower.list)
    upper.list <- unlist(upper.list)
    pred.list
    upper.list
    lower.list
    
    pred.sub <- data.frame(pred.list,lower.list,upper.list,seq(1,length(upper.list)),rep(s,length(upper.list)))
    names(pred.sub) <- c("Pred","upper","lower","bin","Pstep")
    pred.df <-rbind(pred.df,pred.sub)
    
  }
  pred.df
  
  colourslist <- scales::hue_pal()(length(unique(pred.df$Pstep)))
  # Name your list of colors
  names(colourslist) <- unique(pred.df$Pstep)
  
  plt.predict.steps <- ggplot(pred.df,aes(x=bin*n_days,y=Pred,group=Pstep,color=as.factor(Pstep))) + 
    #geom_point(shape=1,size=4) + 
    geom_line() +
    #geom_errorbar(aes(ymin=lower,ymax=upper)) + 
    ggtitle(depVar) +
    ylim(0,1) + 
    ylab('Predictability') + 
    xlab('Days since birth') + 
    theme_classic() + 
    scale_color_discrete(name="Step") +
    guides(colour = guide_legend(reverse = F),) + 
    theme(legend.position = "right",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 14)) 
  plt.predict.steps
  pred.df
  
  pred.onestep <- pred.df[pred.df$Pstep == 1,]
  plt.onestep <- ggplot(pred.onestep,aes(x=bin*n_days,y=Pred)) + 
    #geom_point(shape=1,size=4) + 
    geom_line() +
    #geom_errorbar(aes(ymin=lower,ymax=upper)) + 
    ggtitle(s) +
    #ylim(0,1) + 
    ylab('Predictability') + 
    xlab('Days since birth') + 
    theme_classic() + 
    #scale_color_discrete(name="Step") +
    #guides(colour = guide_legend(reverse = F),) + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 14)) 
  plt.onestep
  #return(list(plt.week.corr, plt.predict.one,behav.bin.id,n_bins,id.matrix.bin,ci.bin))
  return(list(plt.week.corr, plt.predict.one,behav.bin.id,n_bins,among.corr,df.wide))
}

plots.predict.dist <- func.ndays.predict('dist_mean',indv.long54,5)


### Binwise correlation
# There is some weirdness in scaling with MCMCglmm
# As a sanity check, we ran this with standard regression, 
# and the values correlate closely. 

bin_names <- colnames(plots.predict.dist[[6]])[3:13]
among.corr.d <- plots.predict.dist[[5]]

for (r in seq(dim(among.corr.d)[1])) { 
  #print(among.corr.100[r,1:2])
  bar <- as.formula(paste(among.corr.d[r,1],'~',among.corr.d[r,2]))

  res.bar <- lme(bar,random=~1|obs.day,data=plots.predict.dist[[6]])
  simp.corr <- res.bar$coefficients$fixed[2]
  among.corr.d[r,'lme.corr'] <- simp.corr
}
cor(among.corr.d$value,among.corr.d$lme.corr)

plots.predict.velC <- func.ndays.predict('velC_mean',indv.long54,5)
plots.predict.vel <- func.ndays.predict('vel_mean',indv.long54,5)

plots.predict.wall_dist <- func.ndays.predict('pDist_meanScale',indv.long54,5)
plots.predict.wall_distC <- func.ndays.predict('pDistC_meanScale',indv.long54,5)
plots.predict.angleC <- func.ndays.predict('angleC_meanScale',indv.long54,5)

#plots.predict.dist.new[[2]]; plots.predict.dist.new[[1]]
#plots.predict.angleC[[2]]; plots.predict.angleC[[1]]
### Building the whole matrix fig is quite complicated

### Define function to plot massive correlation plot
func.megafig <- function(plots.predict) {
  
  
  behav.bin.id <- plots.predict[[3]]
  n_bins <- plots.predict[[4]]
  among.corr <- plots.predict[[5]]
  
  #id.matrix.bin <- plots.predict.dist[[5]]
  #ci.bin <- plots.predict.dist[[6]]
  
  ## Need to have weekly blups?
  # n_bins:len(Sol)
  binwise.blups <- data_frame(Trait = colnames(behav.bin.id$Sol)[(n_bins + 1):ncol(behav.bin.id$Sol)],
                              Value = unname(colMedians(behav.bin.id$Sol)[(n_bins + 1):ncol(behav.bin.id$Sol)])) %>%
    separate(Trait, into = c("bin", "pi", "comp")) %>%
    mutate(picomp = paste(pi, comp, sep = "_"),
           bin = substr(bin, nchar(bin)-3, nchar(bin))) %>%
    select(-pi, -comp) %>%
    spread(bin, Value) 
  
  #binwise.blups
  
  
  ranking.df <- plots.dist[[6]][c("picomp","ranking")] %>% 
    arrange(picomp) %>% 
    unique()
  
  #ranking.df <- ranking.df[ranking.df$picomp != 'pi31',]
  
  ranking <- ranking.df$ranking
  
  
  bin_names <- colnames(binwise.blups)[2:length(colnames(binwise.blups))]
  bin_names.full <- bin_names
  bin_names.full[[11]] <- "bin10"
  n_elements <- (n_bins * (n_bins - 1)) / 2
  
  plots.list <- vector("list",n_elements)
  n_count <- 0
  
  matrix.fig <- plots.predict[[2]]
  print(binwise.blups)
  x_min <- min(binwise.blups[,seq(2,11)])
  print(x_min)
  x_max <- max(binwise.blups[,seq(2,11)])
  print(x_max)
  for (i in seq(n_bins)) {
    for (j in seq(n_bins)) {
      if(i >= j) next
      x_min <- min(binwise.blups[,j+1])
      x_max <- max(binwise.blups[,j+1])
      y_min <- min(binwise.blups[,i+1])
      y_max <- max(binwise.blups[,i+1])
      n_count <- n_count + 1
      slope <- among.corr[among.corr$Var1 == bin_names.full[i] & among.corr$Var2 == bin_names.full[j],]$value

      single.plot <- ggplot(binwise.blups, aes(x = .data[[bin_names[j]]], y = .data[[bin_names[[i]]]])) +
        geom_abline(intercept = 0, slope = slope) +
        geom_point(aes(color = ranking), size = 2) +
        #geom_point(aes(x=bin8,y=bin1)) + 
        scale_color_viridis_c(option = "plasma") +
        #scale_x_continuous(limits = c(-2, 2)) + #, breaks = c(-0.75, 0, 0.75, 1.5, 2.25)) +
        #scale_y_continuous(limits = c(-2, 2)) + #, breaks = c(-0.75, 0, 0.75, 1.5, 2.25)) +
        scale_x_continuous(limits = c(x_min, x_max)) + 
        scale_y_continuous(limits = c(x_min, x_max)) +
        theme_classic() +
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.text = element_blank(),
              plot.margin = margin(0,0,0,0))
      plots.list[[n_count]] <- single.plot
      print(c(n_count,i,j))
      matrix.fig <- matrix.fig + single.plot
    }
  }
  
  ## woof. There might be a better way to do this...
  layout <- "
  123456789a
  #AbBcCdDeE
  ##fFgGhHiI
  ###jJkKlLm
  ####MnNoOp
  #####PqQrR
  0000##sStT
  0000###uUv
  0000####Vw
  0000#####W
  "
  
  matrix.fig <- matrix.fig + plot_layout(design = layout)
  
  matrix.fig
  return(matrix.fig)
}

### Main fgure

megafig.dist <- func.megafig(plots.predict.dist)
megafig.dist

ggsave('~/Documents/Scripts/GroupMollyTracking/figs/plot.F3.megafig.dist.jpg',megafig.dist,width = 6.5,height=6.5,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/plot.F3.megafig.dist.svg',megafig.dist,width = 6.5,height=6.5,units="in")

## Supplemental Figures
megafig.vel <- func.megafig(plots.predict.vel)
megafig.vel
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.S2.megafig.vel.jpg',megafig.vel,width = 6.5,height=6.5,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.S2.megafig.vel.svg',megafig.vel,width = 6.5,height=6.5,units="in")

megafig.velC <- func.megafig(plots.predict.velC)

megafig.velC


ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.S2.megafig.velC.jpg',megafig.vel,width = 6.5,height=6.5,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.S2.megafig.velC.svg',megafig.vel,width = 6.5,height=6.5,units="in")

megafig.pDist <- func.megafig(plots.predict.wall_dist)
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.S2.megafig.pDist.jpg',megafig.vel,width = 6.5,height=6.5,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.S2.megafig.pDist.svg',megafig.vel,width = 6.5,height=6.5,units="in")

megafig.pDistC <- func.megafig(plots.predict.wall_distC)
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.S2.megafig.pDistC.jpg',megafig.pDistC,width = 6.5,height=6.5,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.S2.megafig.pDistC.svg',megafig.pDistC,width = 6.5,height=6.5,units="in")

megafig.angleC <- func.megafig(plots.predict.angleC)
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.S2.megafig.angleC.jpg',megafig.angleC,width = 6.5,height=6.5,units="in")
ggsave('~/Documents/Scripts/GroupMollyTracking/figs/SupPlots.S2.megafig.angleC.svg',megafig.angleC,width = 6.5,height=6.5,units="in")


### I think the last thing I need is to check correlation among behaviors. 
# This more fine scale I can get the better. Hourly might be good enough to start. 
prior.cov6 <- list(R = list(V = diag(6), nu = 6.002),
                   G = list(G1=list(V = diag(6), nu = 6.002, alpha.mu = rep(0,6), alpha.V = 1000*diag(6))))

behav.all.cor <- MCMCglmm(cbind(dist_meanScale, vel_meanScale, velC_meanScale, pDist_meanScale,pDistC_meanScale,angleC_meanScale) ~ trait - 1, 
                          random = ~us(trait):Pi, 
                          rcov = ~us(trait):units,
                          family = c(rep("gaussian", 6)), 
                          prior = prior.cov6, 
                          nitt = 510000, thin = 200, burnin = 10000, 
                          verbose = F,
                          data = indv.com)

# Model for entire observation period  ----
behav.matrix.all <- matrix(posterior.mode(posterior.cor(behav.all.cor$VCV[,1:36])),6,6, 
                           dimnames = list(c("dist", "vel", "velC", "pDist","pDistC","angleC"), 
                                           c("dist", "vel", "velC", "pDist","pDistC","angleC")))


# now to extract the CI estimates
ci.all <- data.frame(HPDinterval(posterior.cor(behav.all.cor$VCV[,1:36])))

# for corrplot need 3 matrices - estimates, lower CI, upper CI
lower.all <- matrix(ci.all[,1],6,6)
upper.all <- matrix(ci.all[,2],6,6)

test <- melt(lower.all) %>%
  mutate(p.value = ifelse(value < 0, 1, 0)) %>%
  select(Var1, Var2, p.value)

p.mat <- diag(6)
p.mat[cbind(test$Var1, test$Var2)] <- p.mat[cbind(test$Var2, test$Var1)] <- test$p.value

colnames(p.mat) <- c("dist", "vel", "velC", "pDist","pDistC","angleC")
row.names(p.mat) <- c("dist", "vel", "velC", "pDist","pDistC","angleC")

corrplot(behav.matrix.all, type = "upper", method = "ellipse", p.mat = p.mat, insig = "blank")



hourly.com$dist_meanScale <- scale(hourly.com$dist_mean)
hourly.com$vel_meanScale <- scale(hourly.com$vel_mean)
hourly.com$pDist_meanScale <- scale(hourly.com$pDist_mean)

hourly.com$velC_meanScale <- scale(hourly.com$velC_mean)
hourly.com$pDistC_meanScale <- scale(hourly.com$pDistC_mean)
hourly.com$angleC_meanScale <- scale(hourly.com$angleC_mean)

behav.hourly.cor <- MCMCglmm(cbind(dist_meanScale, vel_meanScale, velC_meanScale, pDist_meanScale,pDistC_meanScale,angleC_meanScale) ~ trait - 1, 
                          random = ~us(trait):Pi, 
                          rcov = ~us(trait):units,
                          family = c(rep("gaussian", 6)), 
                          prior = prior.cov6, 
                          nitt = 510000, thin = 200, burnin = 10000, 
                          verbose = T,
                          data = hourly.com)


behav.matrix.hourly <- matrix(posterior.mode(posterior.cor(behav.hourly.cor$VCV[,1:36])),6,6, 
                           dimnames = list(c("dist", "vel", "velC", "pDist","pDistC","angleC"), 
                                           c("dist", "vel", "velC", "pDist","pDistC","angleC")))


# now to extract the CI estimates
ci.hourly <- data.frame(HPDinterval(posterior.cor(behav.hourly.cor$VCV[,1:36])))

# for corrplot need 3 matrices - estimates, lower CI, upper CI
lower.hourly <- matrix(ci.hourly[,1],6,6)
upper.hourly <- matrix(ci.hourly[,2],6,6)

test.hourly <- melt(lower.hourly) %>%
  mutate(p.value = ifelse(value < 0, 1, 0)) %>%
  select(Var1, Var2, p.value)

p.mat <- diag(6)
p.mat[cbind(test.hourly$Var1, test.hourly$Var2)] <- p.mat[cbind(test.hourly$Var2, test.hourly$Var1)] <- test.hourly$p.value

colnames(p.mat) <- c("dist", "vel", "velC", "pDist","pDistC","angleC")
row.names(p.mat) <- c("dist", "vel", "velC", "pDist","pDistC","angleC")

corrplot(behav.matrix.hourly, type = "upper", method = "ellipse", p.mat = p.mat, insig = "blank")

corrplot(behav.matrix.hourly, type = "upper", method = "ellipse", insig = "blank")
