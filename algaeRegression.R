
library('lme4')
library('MuMIn')
library('DHARMa')
library('ggplot2')

algaeDF <- read.csv('~/Documents/Scripts/GroupMollyTracking/algae_df2.csv',header=FALSE)

names(algaeDF) <- c('pixel','x','y','pi','day','algae','count')

algaeDF.clean <- algaeDF[algaeDF$count < 500,]
algaeDF.clean <- algaeDF.clean[algaeDF.clean$count > 0,]

algaeDF.clean$logCount <- log(algaeDF.clean$count)
ggplot(data = algaeDF.clean,aes(count)) + geom_histogram()
ggplot(data = algaeDF.clean,aes(logCount)) + geom_histogram()
ggplot(data = algaeDF,aes(algae)) + geom_histogram()


model.log <- lmer(formula = logCount ~ algae + (1|pi) + (1|pixel),data = algaeDF.clean)
summary(model.log)
MuMIn::r.squaredGLMM(model.log)

simulationOutput <- simulateResiduals(fittedModel = model.log, plot = F)
plot(simulationOutput)

## I could do poisson instead of log transformed, but the above seems good. 
#model.poisson <- glmer(formula = count ~ algae + (1|pi) + (1|pixel),family = poisson, data = algaeDF)

