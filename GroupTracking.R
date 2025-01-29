### Organized according to the paper: 

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

#indv <- read.csv("JolleTracksAll_4.csv")
indv <- read.csv("JolleTracksAll_5.csv")
indv$ExpDay <- as.integer(indv$ExpDay)

indv[indv$Pi == 'pi11',]$ExpDay <- indv[indv$Pi == 'pi11',]$ExpDay - 16
#indv[indv$Pi == 'pi12',]$ExpDay <- indv[indv$Pi == 'pi12',]$ExpDay + 15

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

#long_data54 <- indv %>%
#  filter(Pi %in% piDays[piDays$max >= 54,]$Pi)
#indv.long54 <- long_data54[long_data54$ExpDay <= 54,]

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
res.cohXtime <- lme(dist_mean ~ ExpDay, random = ~1|Pi, data = good_data)
summary(res.time)

res.velXtime <- lme(vel_mean ~ ExpDay, random = ~1|Pi, data = good_data)
summary(res.velXtime)

res.pdistXtime <- lme(pDist_mean ~ ExpDay, random = ~1|Pi, data = good_data)
summary(res.pdistXtime)

res.pdistCXtime <- lme(pDistC_mean ~ ExpDay, random = ~1|Pi, data = good_data)
summary(res.pdistCXtime)

res.activityXtime <- lme(prop_active ~ ExpDay, random = ~1|Pi, data = good_data)
summary(res.activityXtime)

res.uppervelXtime <- lme(upper_vel ~ ExpDay, random = ~1|Pi, data = good_data)
summary(res.uppervelXtime)
### Add max velocity and activity time

### Control for fish size (but it's decreasing, so as a function of fish size will just increase it)

### Velocity cohesion doesn't change, just hangs out around 0.03, low but sig
res.velCXtime <- lme(velC_mean ~ ExpDay, random = ~1|Pi, data = good_data)
summary(res.velCXtime)

### Plotting group blups over time

# Need to build long df
long_data54 <- indv %>%
  filter(Pi %in% piDays[piDays$max >= 54,]$Pi)
indv.long54 <- long_data54[long_data54$ExpDay <= 54,]

## Define priors used throughout for MCMCglm

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

### Define function to get blups ###

func.ndays.intercepts <- function(depVar,df,n_days) {
  # select only those IDs in that vector & only keep up to obs 70 & make obs 1 = 0
  
  indv.ndays.cen <- df
  
  n_pis <- length(unique(df$Pi))
  rpt <- list()
  ci.rpt <- list()
  post.id <- list()
  ci.id <- list()
  post.w <- list()
  ci.w <- list()
  blup <- list()
  
  for (i in seq(0,54,n_days)) {
    
    col_name <- paste('Day',i,sep='')
    print(col_name)
    indv.ndays.cen[[col_name]] <- with(df,ExpDay-i)
    fixed <- as.formula(paste(depVar,"~ ",col_name,sep=""))
    random <- as.formula(paste('~us(1 + ',col_name,'):Pi',sep=''))
    model.blup <- MCMCglmm(fixed = fixed, 
                           random = random, 
                           data = indv.ndays.cen, 
                           family = "gaussian",
                           pr = T,
                           prior = prior.id.slope, 
                           nitt=310000, burnin = 10000, thin = 200, 
                           verbose = F)  
    model_name <- paste('col.day',i,sep='')
    assign(model_name,model.blup)
    rpt_name <- paste('rpt',i,sep='')
    
    ## More proper conditional repeatability approach. 
    sigma.a0 <- model.blup$VCV[,"(Intercept):(Intercept).Pi"]
    a1_col <- paste('Day',i,":Day",i,".Pi",sep='')
    sigma.a1 <- model.blup$VCV[,a1_col]
    rho_col <- paste("(Intercept):Day",i,".Pi",sep='')
    
    rho <- model.blup$VCV[,rho_col] ### this is actually the whole covariance term, not just rho
    sigma.e <- model.blup$VCV[,"units"]
    rpt.n <- ( sigma.a0 + sigma.a1*i + 2*rho*i ) / 
      ( sigma.a0 + sigma.a1*i + 2*rho*i + sigma.e)
    
    #rpt.n <- model.blup$VCV[,"(Intercept):(Intercept).Pi"]/(model.blup$VCV[,"(Intercept):(Intercept).Pi"] + model.blup$VCV[,"units"])
    assign(rpt_name,rpt.n)
    rpt <- c(rpt,median(rpt.n))
    ci.n <- HPDinterval(rpt.n)[1:2]
    ci.rpt <- c(ci.rpt,ci.n)
    
    post.n <- median(model.blup$VCV[,"(Intercept):(Intercept).Pi"])
    post.id <- c(post.id,post.n)
    
    ci.id.n <- HPDinterval(model.blup$VCV[,"(Intercept):(Intercept).Pi"])[1:2]
    ci.id <- c(ci.id,ci.id.n)
    
    post.u <- median(model.blup$VCV[,"units"])
    post.w <- c(post.w,post.u)
    
    ci.w.u <- HPDinterval(model.blup$VCV[,"units"])[1:2]
    ci.w <- c(ci.w,ci.w.u)
    
    'this pulls out the individual intercepts and adds in the overall intercepts
    so that way these numbers are absolute values, as opposed to differences from overall'
    intercepts.n <- unname(median(model.blup$Sol)[3:(3+n_pis-1)] + median(model.blup$Sol)["(Intercept)"])
    blup <- c(blup,intercepts.n)
  }
  blup <- unlist(blup)
  rpt <- unlist(rpt)
  dates <- seq(0,54,n_days) ## 
  ci.rpt <- matrix(unlist(ci.rpt), nrow = length(dates), byrow = T)
  ci.id <-matrix(unlist(ci.id), nrow = length(dates), byrow = T)
  ci.w <- matrix(unlist(ci.w), nrow = length(dates), byrow = T)      
  post.id <- unlist(post.id)
  post.w <- unlist(post.w)
  #plot(post.id ~ date)
  #plot(post.w ~ date)
  #plot(rpt ~ date)
  
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
  ids <- colnames(col.day0$Sol)[3:(3 + n_pis - 1)] ## Make sure this matches below
  ids <- substr(ids, 16, 19)
  
  dates.rep <- rep(dates, each = n_pis)
  picomp <- rep(ids, length(dates))
  
  pred.intercepts <- data.frame(dates.rep, picomp, blup)
  
  plt.repeat <- ggplot(rpt.slice.rpt,aes(x=dates,y=variance)) + 
    geom_point(shape=1) + 
    #geom_line() +
    geom_errorbar(aes(ymin=lower,ymax=upper)) + 
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
  
  pred.rank <- pred.intercepts %>%
    filter(dates.rep == 0) %>%
    mutate(ranking = rank(blup)) %>%
    select(picomp, ranking)
  
  pred.rank <- left_join(pred.intercepts, pred.rank) %>%
    arrange(ranking)
  plasma_pal <- viridis::plasma(n = 30)
  plasma_pal <- plasma_pal[1:26]
  
  blup.plot <- ggplot(pred.rank, aes(x = dates.rep, y = blup, group = picomp)) +
    geom_line(size = 0.8, aes(color = factor(ranking))) + 
    xlab("") +
    scale_y_continuous(name = "Predicted swimming speed (cm/s)") +
    ylab("Predicted swimming speed") +
    scale_x_continuous(breaks = dates) +
    scale_color_manual(values = plasma_pal) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 10),
          legend.position = "none") +
    annotate("text", x = -1, y = 2.0, label = "A", size = 5)
  
  blup.plot
  rpt.plot <- ggplot(rpt.slice.wide, aes(x = dates)) +
    geom_point(aes(y = rpt, color = "#000000"), size = 3) +
    geom_line(aes(x=dates, y = rpt, color = "#000000")) +
    geom_errorbar(aes(ymin = lower.rpt, ymax = upper.rpt, width = 1, color = "#000000")) +
    
    geom_errorbar(aes(x = dates-0.5, ymin = lower.id, ymax = upper.id, width = 1, color = "#000000")) +
    geom_line(aes(x = dates-0.5, y = post.id, color = "#959595")) +
    geom_point(aes(x = dates-0.5, y = post.id), shape = 21, color = "#000000", fill = "#959595", size = 3) +
    
    geom_errorbar(aes(x = dates+0.5, ymin = lower.w, ymax = upper.w, width = 1, color = "#000000")) +
    geom_line(aes(x = dates+0.5, y = post.w, color = "#CCCCCC")) +
    geom_point(aes(x = dates + 0.5, y = post.w), shape = 21, color = "#000000", fill = "#CCCCCC", size = 3) +
    
    scale_x_continuous(breaks = dates, labels = dates) +
    scale_y_continuous(name = "Variance estimate") + #, limits = c(0, 2.5), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    labs(x = "Day",
         color = "Legend") +
    scale_color_manual(name = "", 
                       values =c("#000000",  "#959595",  "#CCCCCC"),
                       labels = c("Repeatability", "Among-Group", "Within-Group")) +
    theme_classic()+
    theme(legend.position = 'none', #c(0.3, 0.90),
          legend.key.size = unit(0.3, 'cm'),
          legend.text = element_text(size = 6),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 10),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 10)) #+
  #annotate("text", label = "B", size = 5, x = -1, y = 0.8)
  
  rpt.plot
  plt.repeat
  plt.day <- ggplot(pred.intercepts, aes(x=dates.rep, y=blup, group=picomp, color=picomp)) + 
    geom_line() + 
    ggtitle(depVar) +
    theme_classic() + 
    #ylim(ylims) + 
    ylab('blup') + 
    xlab('Days since birth') + 
    theme(legend.position = "none",
          rect=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(size = 12, colour = "black"),
          axis.title = element_text(size = 14, face = "bold", colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 14))
  
  plt.day
  return(list(plt.day,plt.repeat,rpt.plot,blup.plot,rpt.slice.wide,pred.rank))
}
  
## These two are in the paper proper 
plots.dist <- func.ndays.intercepts('dist_mean',indv.long54,n_days)
plots.velC <- func.ndays.intercepts('velC_mean',indv.long54,n_days)

## These Four are supplemental
plots.vel <- func.ndays.intercepts('vel_mean',indv.long54,n_days)
plots.wall_dist <- func.ndays.intercepts('pDist_mean',indv.long54,n_days)
#plots.vel_std <- func.ndays.intercepts('vel_std',df,n_days)
#plots.velC_scale <- func.ndays.intercepts('velC_scale',df,n_days)
plots.wall_distC <- func.ndays.intercepts('pDistC_mean',indv.long54,n_days)
plots.angleC <- func.ndays.intercepts('angleC_mean',indv.long54,n_days)
#plots.angle_std <- func.ndays.intercepts('polarity_std',df,n_days)

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

sliding.vel <- func.slidingMean('vel_mean',indv.long54,5)
#sliding.vel_std <- func.slidingMean('vel_std',indv.long54,5)
#sliding.velC_scale <- func.slidingMean('velC_scale',indv.long54,5)
sliding.wall_dist <- func.slidingMean('pDist_mean',indv.long54,5)
sliding.wall_distC <- func.slidingMean('pDistC_mean',indv.long54,5)
sliding.angleC <- func.slidingMean('angleC_mean',indv.long54,5)
#sliding.angle_std <- func.slidingMean('polarity_std',indv.long54,5)

### Calculate repeatability over time ###

## rpt.data is calculated as part of intercepts, so you can get it as plots.foo[[5]]
func.rpt.plot <- function(rpt.data) {
  bar_scale <- max(rpt.data$post.id) + max(rpt.data$post.w)
  #foo <- plots.vel[[5]]
  rpt.data$lower.id_ <- rpt.data$lower.id / bar_scale
  rpt.data$upper.id_ <- rpt.data$upper.id / bar_scale
  rpt.data$post.id_ <- rpt.data$post.id / bar_scale
  rpt.data$lower.w_ <- rpt.data$lower.w / bar_scale
  rpt.data$upper.w_ <- rpt.data$upper.w / bar_scale
  rpt.data$post.w_ <- rpt.data$post.w / bar_scale
  
  rpt.plot <- ggplot(rpt.data, aes(x = dates)) +
    geom_point(aes(y = rpt, color = "#000000"), size = 3) +
    geom_line(aes(x=dates, y = rpt, color = "#000000")) +
    geom_errorbar(aes(ymin = lower.rpt, ymax = upper.rpt, width = 1, color = "#000000")) +
    
    #geom_errorbar(aes(x = dates-0.5, ymin = lower.id_, ymax = upper.id_, width = 1, color = "#000000")) +
    geom_line(aes(x = dates-0.5, y = post.id_, color = "#959595")) +
    geom_point(aes(x = dates-0.5, y = post.id_), shape = 21, color = "#000000", fill = "#959595", size = 3) +
    
    #geom_errorbar(aes(x = dates+0.5, ymin = lower.w_, ymax = upper.w_, width = 1, color = "#000000")) +
    geom_line(aes(x = dates+0.5, y = post.w_, color = "#CCCCCC")) +
    geom_point(aes(x = dates + 0.5, y = post.w_), shape = 21, color = "#000000", fill = "#CCCCCC", size = 3) +
    
    scale_x_continuous(breaks = dates, labels = dates) +
    scale_y_continuous(name = "Variance estimate",limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    labs(x = "Day",
         color = "Legend") +
    scale_color_manual(name = "", 
                       values =c("#000000",  "#959595",  "#CCCCCC"),
                       labels = c("Repeatability", "Among-Group", "Within-Group")) +
    theme_classic()+
    theme(legend.position = 'none', #c(0.3, 0.90),
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
rpt.plot.dist <- func.rpt.plot(plots.dist[[5]])
rpt.plot.velC <- func.rpt.plot(plots.velC[[5]])

rpt.plot.vel <- func.rpt.plot(plots.vel[[5]])
#rpt.plot.vel_std <- func.rpt.plot(plots.vel_std[[5]])
rpt.plot.wall_dist <- func.rpt.plot(plots.wall_dist[[5]])
rpt.plot.wall_distC <- func.rpt.plot(plots.wall_distC[[5]])
rpt.plot.angleC <- func.rpt.plot(plots.angleC[[5]])
#rpt.plot.angle_std <- func.rpt.plot(plots.angle_std[[5]])

## First two rows of plots: 
rpt.plot.dist; sliding.dist[[1]]
rpt.plot.velC; sliding.velC[[1]]

## Supplemental plots: 
rpt.plot.vel; sliding.vel[[1]]
rpt.plot.wall_dist; sliding.wall_dist[[1]]
rpt.plot.wall_distC; sliding.wall_distC[[1]]
rpt.plot.angleC; sliding.angleC[[1]]

### Difference vs random (row 3 in fig 2) is calculated using shuffleTracks.py

## Define function to get predictability over time

func.ndays.predict <- function(depVar,df,n_days) {
  #df <- indv.long54
  df$velC_scale <- df$velC_mean * 100
  
  df
  ## build bin column
  df$bin <- df$ExpDay %/% n_days
  
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
  
  behav.bin.id <- MCMCglmm(form.bins, 
                           random = ~us(trait):Pi,
                           rcov = ~us(trait):units,
                           family = c(rep("gaussian", n_bins)), 
                           prior = prior.b, 
                           pr = T, 
                           nitt = 510000, thin = 200, burnin = 10000, 
                           verbose = F,
                           data = df.wide)
  
  # Model with only ID ----
  
  id.matrix.bin <- matrix(median(posterior.cor(behav.bin.id$VCV[,1:(n_bins * n_bins)])),n_bins,n_bins, 
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
  
  ## Need these for separate plots
  #test.lower <- melt(replace(lower.bin, lower.tri(lower.bin, T), NA), na.rm = T)
  #test.upper <- melt(replace(upper.bin, lower.tri(upper.bin, T), NA), na.rm = T)
  
  #ci.long <- left_join(test.lower, test.upper, by = c("Var1", "Var2")) %>%
  #  rename(start.bin = Var1,
  #         end.bin = Var2,
  #         lower = value.x, 
  #         upper = value.y) %>%
  #  arrange(start.bin, end.bin)
  
  #df.corrs <- melt(replace(id.matrix.bin, lower.tri(id.matrix.bin, T), NA), na.rm = T)
  
  #df.corrs$start.bin <- as.numeric(substr(df.corrs$Var1, 4, 5))
  #df.corrs$end.bin <- as.numeric(substr(df.corrs$Var2, 4, 5))
  #df.corrs$diff <- df.corrs$end.bin - df.corrs$start.bin
  
  ## Need to shift to have it 0 indexed
  #ci.long$start.bin <- ci.long$start.bin -1
  #ci.long$end.bin <- ci.long$end.bin -1
  #among.corr <- left_join(df.corrs, ci.long, by = c("start.bin", "end.bin"))

  
  ## So I need id.matrix.week and ci.bin
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
    #ylim(0,1) + 
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
  return(list(plt.week.corr, plt.predict.one,behav.bin.id,n_bins,among.corr))
}

plots.predict.dist <- func.ndays.predict('dist_mean',indv.long54,5)
plots.predict.velC <- func.ndays.predict('velC_mean',indv.long54,5)

plots.predict.vel <- func.ndays.predict('vel_mean',indv.long54,5)
plots.predict.wall_dist <- func.ndays.predict('pDist_mean',indv.long54,5)
plots.predict.wall_distC <- func.ndays.predict('pDistC_mean',indv.long54,5)
plots.predict.angleC <- func.ndays.predict('angleC_mean',indv.long54,5)


plots.predict.dist[[2]]; plots.predict.dist[[1]]

### Building the whole matrix fig is quite complicated

### Define function to plot massive correlation plot
func.megafig <- function(plots.predict) {
  
  
  behav.bin.id <- plots.predict[[3]]
  n_bins <- plots.predict[[4]]
  among.corr <- plots.predict[[5]]
  
  #id.matrix.bin <- plots.predict.dist[[5]]
  #ci.bin <- plots.predict.dist[[6]]
  
  ## Need to have weekly blups
  # n_bins:len(Sol)
  binwise.blups <- data_frame(Trait = colnames(behav.bin.id$Sol)[(n_bins + 1):ncol(behav.bin.id$Sol)],
                              Value = unname(median(behav.bin.id$Sol)[(n_bins + 1):ncol(behav.bin.id$Sol)])) %>%
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
  
  for (i in seq(n_bins)) {
    for (j in seq(n_bins)) {
      if(i >= j) next
      n_count <- n_count + 1
      slope <- among.corr[among.corr$Var1 == bin_names.full[i] & among.corr$Var2 == bin_names.full[j],]$value
      
      single.plot <- ggplot(binwise.blups, aes(x = .data[[bin_names[j]]], y = .data[[bin_names[[i]]]])) +
        geom_abline(intercept = 0, slope = slope) +
        geom_point(aes(color = ranking), size = 2) +
        #geom_point(aes(x=bin8,y=bin1)) + 
        scale_color_viridis_c(option = "plasma") +
        scale_x_continuous(limits = c(-116, 84)) + #, breaks = c(-0.75, 0, 0.75, 1.5, 2.25)) +
        scale_y_continuous(limits = c(-116, 84)) + #, breaks = c(-0.75, 0, 0.75, 1.5, 2.25)) +
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

megafig.dist <- func.megafig(plots.predict.dist)
megafig.vel <- func.megafig(plots.predict.vel)

megafig.vel
