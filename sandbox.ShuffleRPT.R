## Code to permute a dataset in order to get null rpt estimates

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

### Read in all the data
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

### Add max velocity and activity time

### Control for fish size (but it's decreasing, so as a function of fish size will just increase it)

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

## most libraries and indv.long54 is defined as part of the bigger notebook, 
## Evetually that will need to be defined here. 

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
func.shuffled.rpt.i <- function(i) {
  # select only those IDs in that vector & only keep up to obs 70 & make obs 1 = 0
  set.seed(i)
  df <- func.shuffle.df(indv.long54)
  depVar <- "dist_mean"
  n_days <- 5
  
  indv.ndays.cen <- df
  
  n_pis <- length(unique(df$Pi))
  rpt <- list()
  ci.rpt <- list()
  post.id <- list()
  ci.id <- list()
  post.w <- list()
  ci.w <- list()
  blup <- list()
  

    
    ### Do i need to do all these different ones? 
    ### I guess I just do it, export it as a list, and then check
    ### Seems like it is different based on time

  fixed <- as.formula(paste(depVar,"~ ExpDay",sep=""))
  random <- as.formula(paste('~us(1 + ExpDay):Pi',sep=''))
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
  
  sigma.a0 <- model.blup$VCV[,"(Intercept):(Intercept).Pi"]

  sigma.a1 <- model.blup$VCV[,"ExpDay:ExpDay.Pi"]

  rho <- model.blup$VCV[,"(Intercept):ExpDay.Pi"] ### this is actually the whole covariance term, not just rho
  sigma.e <- model.blup$VCV[,"units"]
  for (i in seq(0,54,n_days)) {
      rpt.n <- ( sigma.a0 + sigma.a1*(i**2) + 2*rho*i ) / 
      ( sigma.a0 + sigma.a1*(i**2) + 2*rho*i + sigma.e)
    
    #rpt.n <- model.blup$VCV[,"(Intercept):(Intercept).Pi"]/(model.blup$VCV[,"(Intercept):(Intercept).Pi"] + model.blup$VCV[,"units"])
    assign(rpt_name,rpt.n)
    rpt <- c(rpt,median(rpt.n))
  }
  rpt <- unlist(rpt)
  
  return(rpt) }

### Testing parallel permutations (it works!)
system.time(rpt.i <- func.shuffled.rpt.i(1))
system.time(rpt.is <- mclapply(1:5,func.shuffled.rpt.i,mc.cores=5L))




### Old functions: 
func.shuffled.rpt.i.slice <- function(i) {
  # select only those IDs in that vector & only keep up to obs 70 & make obs 1 = 0
  set.seed(i)
  df <- func.shuffle.df(indv.long54)
  depVar <- "dist_mean"
  n_days <- 5
  
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
    
    ### Do i need to do all these different ones? 
    ### I guess I just do it, export it as a list, and then check
    ### Seems like it is different based on time
    col_name <- paste('Day',i,sep='')
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
  }
  rpt <- unlist(rpt)
  
  return(rpt) }

func.shuffled.intercepts <- function(depVar,df,n_days) {
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
  for (i in seq(0,54,n_days)) {
    
    col_name <- paste('Day',i,sep='')
    print(col_name)
    indv.ndays.cen[[col_name]] <- with(df,ExpDay-i)

    model_name <- paste('col.day',i,sep='')
    assign(model_name,model.blup)
    rpt_name <- paste('rpt',i,sep='')
    
    rpt.n <- model.blup$VCV[,"(Intercept):(Intercept).Pi"]/(model.blup$VCV[,"(Intercept):(Intercept).Pi"] + model.blup$VCV[,"units"])
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

## Original function that gets repeatability but also intercepts and plots.
func.shuffled.intercepts.slice <- function(depVar,df,n_days) {
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
    
    rpt.n <- model.blup$VCV[,"(Intercept):(Intercept).Pi"]/(model.blup$VCV[,"(Intercept):(Intercept).Pi"] + model.blup$VCV[,"units"])
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

plots.dist.shuff <- func.shuffled.intercepts('dist_mean',df.shuffle,n_days)

### Same as above, but broken out so I could debug individual lines more easily
#func.shuffled.intercepts <- function(depVar,df,n_days) {
  # select only those IDs in that vector & only keep up to obs 70 & make obs 1 = 0
  df <- indv.long54
  n_days <- 5
  depVar <- 'dist_mean'
  indv.ndays.cen <- df
  
  n_pis <- length(unique(df$Pi))
  rpt <- list()
  ci.rpt <- list()
  post.id <- list()
  ci.id <- list()
  post.w <- list()
  ci.w <- list()
  blup <- list()
  
  col_name = "Day0"
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
  
  for (i in seq(0,54,n_days)) {
    #i<-5  
    #col_name <- paste('Day',i,sep='')
    #print(col_name)
    #indv.ndays.cen[[col_name]] <- with(df,ExpDay-i)
    #fixed <- as.formula(paste(depVar,"~ ",col_name,sep=""))
    #random <- as.formula(paste('~us(1 + ',col_name,'):Pi',sep=''))
    #model.blup <- MCMCglmm(fixed = fixed, 
    #                       random = random, 
    #                       data = indv.ndays.cen, 
    #                       family = "gaussian",
    #                       pr = T,
    #                       prior = prior.id.slope, 
    #                       nitt=310000, burnin = 10000, thin = 200, 
    #                       verbose = F)  
    #model_name <- paste('col.day',i,sep='')
    #assign(model_name,model.blup)
    rpt_name <- paste('rpt',i,sep='')
    
    #rpt.n <- model.blup$VCV[,"(Intercept):(Intercept).Pi"]/(model.blup$VCV[,"(Intercept):(Intercept).Pi"] + model.blup$VCV[,"units"])
    
    sigma.a0 <- model.blup$VCV[,"(Intercept):(Intercept).Pi"]
    #a1_col <- paste('Day',i,":Day",i,".Pi",sep='')
    a1_col <- "Day0:Day0.Pi"
    sigma.a1 <- model.blup$VCV[,"Day0:Day0.Pi"]
    rho_col <- "(Intercept):Day0.Pi"
    rho_col <- "Day0:(Intercept).Pi" ##What is the difference here?
    
    rho <- model.blup$VCV[,rho_col] ### this is actually the whole covariance term, not just rho
    sigma.e <- model.blup$VCV[,"units"]
    rpt.n <- ( sigma.a0 + sigma.a1*(i**2) + 2*rho*i ) / 
      ( sigma.a0 + sigma.a1*(i**2) + 2*rho*i + sigma.e)
    print(min(sigma.e))
    #print(rpt.n)
    print(min(rpt.n))
    print(max(rpt.n))
    #rpt.new <- ( sigma.a0 + sigma.a1*i + 2*rho*sqrt(sigma.a0)*sqrt(sigma.a1)*i ) /   
    #  ( sigma.a0 + sigma.a1*i + 2*rho*sqrt(sigma.a0)*sqrt(sigma.a1)*i + sigma.e)
     

    #(model.blup$VCV[,"(Intercept):(Intercept).Pi"] + model.blup$VCV[,"Day0:Day0.Pi"]*i + 2*model.blup$VCV[,"(Intercept):(Intercept).Pi"]  /(model.blup$VCV[,"(Intercept):(Intercept).Pi"] + model.blup$VCV[,"units"])
    
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
  
  rpt.slice.long <- data.frame("dates" = rep(dates, 3), 
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
#  return(list(plt.day,plt.repeat,rpt.plot,blup.plot,rpt.slice.wide,pred.rank))
#}
