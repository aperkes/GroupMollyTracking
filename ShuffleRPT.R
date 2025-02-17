## Code to permute a dataset in order to get null rpt estimates

library(nlme)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MCMCglmm)

### Need these for permutation
library(parallel)
library(matrixStats)

### Read in all the data
setwd("~/Documents/Scripts/GroupMollyTracking/")

indv <- read.csv("JolleTracksAll_5.csv")
indv$ExpDay <- as.integer(indv$ExpDay)

#### Do various data cleaning ####
indv[indv$Pi == 'pi11',]$ExpDay <- indv[indv$Pi == 'pi11',]$ExpDay - 16 ## There were some videos before fish went in. 

indv.good <- indv[indv$ExpDay >= 0,]
#indv.good <- indv.good[indv.good$Pi != 'pi12',]

indv <- indv.good

## Some of the pi's fall off over time (e.g., pi41 drops off at day 11)
## Others go on 90 days
piDays <- indv %>%
  group_by(Pi) %>%
  summarize(max = max(ExpDay))

piDays[piDays$max >= 54,]

## Only look at 4 weeks for fish where we have that
good_data <- indv[indv$Pi != "pi34",] ## This one is missing 
good_data <- good_data[good_data$Pi != "pi41",] ## This oen is zoomed out, maybe impossible
#good_data <- good_data[good_data$ExpDay <= 31,] ## This one is missing one day

indv.com <- good_data

### Add max velocity and activity time

### Control for fish size (but it's decreasing, so as a function of fish size will just increase it)

### Plotting group blups over time


#### Need to build long df ####
long_data54 <- indv %>%
  filter(Pi %in% piDays[piDays$max >= 54,]$Pi)
indv.long54 <- long_data54[long_data54$ExpDay <= 54,]

#### Define priors used throughout for MCMCglm ####

## This is the standard prior uses for slope/intercept model
prior.id.slope <- list(R = list(V = 1, nu = 0.002),
                       G = list(G1=list(V = diag(2), nu = 0.002, alpha.mu = c(0,0),  alpha.V = diag(2)*25^2)))


## We need a prior with matrix to calculate the residual for each day
total_days <- nlevels(as.factor(indv.long54$ExpDay))
prior.id.slope.cov <- list(R = list(V = diag(total_days),nu = total_days + 0.002),
                       G = list(G1=list(V = diag(2), nu = 0.002, alpha.mu = c(0,0),  alpha.V = diag(2)*25^2)))

### We also need to have specified priors for each one, see other one

#### Function to calculate conditional repeatabilty, using het residual variance ####
func.rpt <- function(df,depVar='dist_mean',day_bin=1) {
    # select only those IDs in that vector & only keep up to obs 70 & make obs 1 = 0
    max_days <- nlevels(as.factor(df$ExpDay)) - 1 ## It's 0 indexed
    
    ### I'm superstitious about global variables...
    indv.df <- df
    
    n_pis <- length(unique(df$Pi))
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
    indv.df$ExpDay_factor <- as.factor(indv.df$ExpDay)
    model.het <- MCMCglmm(fixed = fixed, 
                           random = random, 
                           rcov = ~idh(ExpDay_factor):units, ## use this line for het. residual variance
                           data = indv.df, 
                           family = "gaussian",
                           pr = T,
                           prior = prior.id.slope.cov, ## Replace this with prior.id.slope for homo residual variance
                           nitt=310000, burnin = 10000, thin = 200, 
                           verbose = F)  

    sigma.a0 <- model.het$VCV[,"(Intercept):(Intercept).Pi"]
    sigma.a1 <- model.het$VCV[,"ExpDay:ExpDay.Pi"]
    
    ## If using fixed individual variance
    #sigma.e <- model.het$VCV[,"units"]
    
    rho <- model.het$VCV[,"(Intercept):ExpDay.Pi"] ### whole covariance term, not just rho
    
    ### Calculate conditional rpt at each point x
    for (x in seq(0,max_days,day_bin)) {
      ## Grab residual at that day
      e_col <- paste("ExpDay_factor",x,".units",sep='')
      sigma.e.x <- model.het$VCV[,e_col]
      
      rpt.x <- ( sigma.a0 + sigma.a1*(x**2) + 2*rho*x ) / 
        ( sigma.a0 + sigma.a1*(x**2) + 2*rho*x + sigma.e.x)

      rpt <- c(rpt,median(rpt.x)) # median is preferred to mode
    }
    rpt <- unlist(rpt)
    
    return(rpt)  
  
}

## Run conditional repeatability (for mean_velocity )
rpt.cond.dist <- func.rpt(indv.long54)
print(rpt.cond.dist)

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
func.shuffled.rpt.i <- function(i,depVar="dist_mean",day_bin=1,prior.cov=prior.id.slope.coh) {
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

### Testing parallel permutations (it works!)
system.time(rpt.i <- func.shuffled.rpt.i(1))
system.time(rpt.is <- mclapply(1:5,func.shuffled.rpt.i,mc.cores=5L))

### convert the list of lists to a matrix and get your CI by row
print(rpt.is)
rpt.matrix <- do.call("cbind",rpt.is)
rowQuantiles(rpt.matrix,probs=c(.05,.95))

## Run for real: 
### Set interations:  vvv below. Be sure to set cores  vv. 
rpt.all <- mclapply(1:50,func.shuffled.rpt.i,mc.cores=5L)
rpt.mat.all <- do.call("cbind",rpt.all)
rpt.quants <- rowQuantiles(rpt.mat.all,probs=c(.025,.975))

rpt.coh <- mclapply(1:50,func.shuffled.rpt.i,depVar='dist_mean',prior.cov=prior.id.slope.coh,mc.cores=5L)
rpt.mat.coh <- do.call("cbind",rpt.coh)
rpt.quants.coh <- rowQuantiles(rpt.mat.coh,probs=c(.025,.975))

rpt.vel <- mclapply(1:50,func.shuffled.rpt.i,depVar='vel_mean',prior.cov=prior.id.slope.vel,mc.cores=5L)
rpt.mat.vel <- do.call("cbind",rpt.vel)
rpt.quants.vel <- rowQuantiles(rpt.mat.vel,probs=c(.025,.975))

rpt.velC <- mclapply(1:50,func.shuffled.rpt.i,depVar='velC_mean',prior.cov=prior.id.slope.velC,mc.cores=5L)
rpt.mat.velC <- do.call("cbind",rpt.velC)
rpt.quants.velC <- rowQuantiles(rpt.mat.velC,probs=c(.025,.975))

rpt.pdist <- mclapply(1:50,func.shuffled.rpt.i,depVar='pDist_mean',prior.cov=prior.id.slope.pDist,mc.cores=5L)
rpt.mat.pdist <- do.call("cbind",rpt.pdist)
rpt.quants.pdist <- rowQuantiles(rpt.mat.pdist,probs=c(.025,.975))

rpt.pdistC <- mclapply(1:50,func.shuffled.rpt.i,depVar='pDistC_mean',prior.cov=prior.id.slope.pDistC,mc.cores=5L)
rpt.mat.pdistC <- do.call("cbind",rpt.pdistC)
rpt.quants.pdistC <- rowQuantiles(rpt.mat.pdistC,probs=c(.025,.975))

rpt.angleC <- mclapply(1:50,func.shuffled.rpt.i,depVar='angleC_mean',prior.cov=prior.id.slope.angleC,mc.cores=5L)
rpt.mat.angleC <- do.call("cbind",rpt.angleC)
rpt.quants.angleC <- rowQuantiles(rpt.mat.angleC,probs=c(.025,.975))


### Dumping ground for all sorts o fold code: 

res.cohScale <- lme(dist_meanScale ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
res.velCScale <- lme(velC_meanScale ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
res.velScale <- lme(vel_meanScale ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
res.pDistScale <- lme(pDist_meanScale ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
res.pDistCScale <- lme(pDistC_meanScale ~ ExpDay, random = ~ExpDay|Pi, data = good_data)
res.angleCScale <- lme(angleC_meanScale ~ ExpDay, random = ~ExpDay|Pi, data = good_data)


###
func.build_prior <- function(res.slopeInt,total_days = 55) {
  
  res.resid <- as.double(VarCorr(res.slopeInt)["Residual","Variance"])
  res.corr <- as.double(VarCorr(res.slopeInt)["ExpDay",'Corr'])
  res.slope <- as.double(VarCorr(res.slopeInt)["ExpDay",'Variance'])
  res.int <- as.double(VarCorr(res.slopeInt)["(Intercept)",'Variance'])
  res.matrix <- matrix(c(res.int,0,0,res.slope),nrow=2,ncol=2,byrow=TRUE)
  
  prior <- list(R = list(V = diag(total_days)*res.resid,nu = total_days + 0.002),
                G = list(G1=list(V = diag(2), nu = 0.002, alpha.mu = c(res.int,res.slope),  alpha.V = diag(2)*25^2)))
  print(res.resid)
  print(res.matrix)
  return(prior)
}

prior.coh <- func.build_prior(res.cohXtime)
prior.velC <- func.build_prior(res.velCXtime)
prior.vel <- func.build_prior(res.velXtime)
prior.pDist <- func.build_prior(res.pDistXtime)
prior.pDistC <- func.build_prior(res.pDistCXtime)
prior.angleC <- func.build_prior(res.angleCXtime)


### Scaled Priors
prior.cohScale <- func.build_prior(res.cohScale)
prior.velCScale <- func.build_prior(res.velCScale)
prior.velScale <- func.build_prior(res.velScale)
prior.pDistScale <- func.build_prior(res.pDistScale)
prior.pDistCScale <- func.build_prior(res.pDistCScale)
prior.angleCScale <- func.build_prior(res.angleCScale)


prior.null <- list(R = list(V = 1, nu = 0.002))

prior.id <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1=list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 25^2)))

prior.id.slope.cov2 <- list(R = list(V = diag(total_days),nu = total_days + 0.002),
                            G = list(G1=list(V = diag(2), nu = 2.002, alpha.mu = c(0,0),  alpha.V = diag(2)*25^2)))

prior.id.slope.cov1 <- list(R = list(V = diag(total_days),nu = total_days + 0.002),
                            G = list(G1=list(V = diag(2), nu = 1.002, alpha.mu = c(0,0),  alpha.V = diag(2)*25^2)))


prior.good <- list(R = list(V = diag(total_days),nu = total_days + 0.002),
                   G = list(G1=list(V = diag(2), nu = 0.002, alpha.mu = c(0.5,0.5),  alpha.V = diag(2)*25^2)))

prior.better <- list(R = list(V = diag(total_days),nu = total_days + 0.002),
                     G = list(G1=list(V = diag(2)*0.5, nu = 0.002, alpha.mu = c(0,0),  alpha.V = diag(2)*25^2)))



plots.coh2.scaled.prior <- func.ndays.intercepts.het('dist_meanScale',indv.long54,n_days,prior.cov = prior.id.slope.cov2)
rpt.plot.coh2.scaled <- func.rpt.plot(plots.coh2.scaled.prior[[5]])
rpt.plot.coh2.scaled

plots.coh1.scaled.prior <- func.ndays.intercepts.het('dist_meanScale',indv.long54,n_days,prior.cov = prior.id.slope.cov1)
rpt.plot.coh1.scaled <- func.rpt.plot(plots.coh1.scaled.prior[[5]])
rpt.plot.coh1.scaled

plots.coh0.scaled.prior <- func.ndays.intercepts.het('dist_meanScale',indv.long54,n_days,prior.cov = prior.id.slope.cov)
rpt.plot.coh0.scaled <- func.rpt.plot(plots.coh0.scaled.prior[[5]])
rpt.plot.coh0.scaled

## Using the good prior
plots.coh0.scaled.good <- func.ndays.intercepts.het('dist_meanScale',indv.long54,n_days,prior.cov = prior.good)
rpt.plot.good <- func.rpt.plot(plots.coh0.scaled.good[[5]])
rpt.plot.good

plots.zeros <- func.ndays.intercepts.het('dist_meanScale',indv.long54,n_days,prior.cov = prior.zeros)
rpt.zeros <- func.rpt.plot(plots.zeros[[5]])
rpt.zeros

plots.better <- func.ndays.intercepts.het('dist_meanScale',indv.long54,n_days,prior.cov = prior.better)
rpt.better <- func.rpt.plot(plots.better[[5]])
rpt.better

plots.best <- func.ndays.intercepts.het('dist_meanScale',indv.long54,n_days,prior.cov = prior.best)
rpt.best <- func.rpt.plot(plots.best[[5]])
rpt.best

plots.cohWI.prior <- func.ndays.intercepts.het('dist_mean',indv.long54,n_days,prior.cov = prior.coh)
rpt.plot.cohWI <- func.rpt.plot(plots.cohWI.prior[[5]])
rpt.plot.cohWI

plots.cohWI.scaled.prior <- func.ndays.intercepts.het('dist_meanScale',indv.long54,n_days,prior.cov = prior.cohScale)
rpt.plot.cohWI.scaled <- func.rpt.plot(plots.cohWI.scaled.prior[[5]])
rpt.plot.cohWI.scaled

plots.cohWI.priorPE <- func.ndays.intercepts.het('dist_mean',indv.long54,n_days,prior.cov = prior.coh)
rpt.plot.cohWIPE <- func.rpt.plot(plots.cohWI.priorPE[[5]])
rpt.plot.cohWIPE

plots.cohWI.scaled.priorPE <- func.ndays.intercepts.het('dist_meanScale',indv.long54,n_days,prior.cov = prior.cohScale)
rpt.plot.cohWIPE.scaled <- func.rpt.plot(plots.cohWI.scaled.priorPE[[5]])
rpt.plot.cohWIPE.scaled

plots.velC.scaled <- func.ndays.intercepts.het('velC_meanScale',indv.long54,n_days)
plots.velC.scaled.prior <- func.ndays.intercepts.het('velC_meanScale',indv.long54,n_days,prior.cov = prior.velCScale)


sliding.velC.scaled <- func.slidingMean('velC_meanScale',indv.long54,5)
rpt.plot.velC.scaled <- func.rpt.plot(plots.velC.scaled[[5]])
rpt.plot.velC.scaled.prior <- func.rpt.plot(plots.velC.scaled.prior[[5]])

rpt.plot.velC.scaled; sliding.velC.scaled[[1]]
rpt.plot.velC; rpt.plot.velC.scaled.prior; rpt.plot.velC.scaled
