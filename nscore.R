## https://www.msu.edu/~ashton/research/code/nscore.R
# nscore.R
# Experimental R functions for performing normal score transforms 
# and back transforms on data. This is useful when 
# conducting geostatistical Gaussian simulation.
# 
# Kriging with normal scores may be necessary though not 
# sufficient for ensuring that a problem reflects the 
# multiGaussian assumptions of sequential Gaussian simulation.
# Checks are available; see Goovaerts (1999) chapter 7.
#
# Functions were developed that are loosely based on fortran 
# routines in GSLIB v.2; see Deutsch & Journel, 1998. Errors, 
# bugs, and deviations are all due to my coding and 
# implementation decisions, however.
#
# The back transform implemented here interpolates linearly 
# between data. Extrapolation is also linear. See backtr() 
# code for details and options. This back transform implementation
# is problematic for two reasons:
# 1. How to evaluate ties (it's a rank transform)
# 2. Extrapolating for small and large values.
#
# Valuable discussions on AI-Geostats archives as well as
# in the formal literature consider these issues, and possibly
# much better ideas are out there (e.g. Saito & Goovaerts (2000))
#
# That said, extrapolation decisions in this code are largely 
# theory-free and the user is warned to treat such values with suspicion.
#
# These functions are experimental and are not written to be 
# robust. In particular, function inputs are not tested for 
# validity, there's no error handling, etc. etc. It is provided
# in hopes that it will be stimulating.
#
# For examples, try:
# source('nscore.R')   # loads the functions, runs nothing
# example1.nscore()
# example2.nscore()
#
# References
# Deutsch, C.V. and Journel, A.G. (1998) GSLIB: Geostatistical
# Software Library and User's Guide. New York:Oxford.
#
# Goovaerts, P. (1997) Geostatistics for Natural Resources
# Evaluation. New York:Oxford.
#
# Dubois, G. (2001) AI-GEOSTATS: SUMMARY: Nscore transform & 
# kriging of log normal data sets.
# http://www.mail-archive.com/ai-geostats@jrc.it/msg00152.html
#
# Saito, H. & Goovaerts, P. (2000) Geostatistical interpolation of 
# positively skewed and censored data in a dioxin-contaminated site.
# Environmental Science & Technology 34(19): 4228-4235.
#
#
# written by Ashton Shortridge, May/June, 2008.

nscore <- function(x) {
   # Takes a vector of values x and calculates their normal scores. Returns 
   # a list with the scores and an ordered table of original values and
   # scores, which is useful as a back-transform table. See backtr().
   nscore <- qqnorm(x, plot.it = FALSE)$x  # normal score 
   trn.table <- data.frame(x=sort(x),nscore=sort(nscore))
   
   return (list(nscore=nscore, trn.table=trn.table))
}

backtr <- function(scores, nscore, tails='none', draw=TRUE) {
   # Given a vector of normal scores and a normal score object 
   # (from nscore), the function returns a vector of back-transformed 
   # values. One major issue is how to extrapolate to the tails. Options 
   # other than none may result in dramatically incorrect tail estimates!
   # tails options:
   # 'none' : No extrapolation; more extreme score values will revert 
   # to the original min and max values. 
   # 'equal' : Calculate magnitude in std deviations of the scores about 
   # initial data mean. Extrapolation is linear to these deviations. 
   # will be based upon deviations from the mean of the original 
   # hard data - possibly quite dangerous!
   # 'separate' :  This calculates a separate sd for values 
   # above and below the mean.
   
   if(tails=='separate') { 
      mean.x <- mean(nscore$trn.table$x)
      small.x <- nscore$trn.table$x < mean.x
      large.x <- nscore$trn.table$x > mean.x
      small.sd <- sqrt(sum((nscore$trn.table$x[small.x]-mean.x)^2)/
                       (length(nscore$trn.table$x[small.x])-1))
      large.sd <- sqrt(sum((nscore$trn.table$x[large.x]-mean.x)^2)/
                       (length(nscore$trn.table$x[large.x])-1))
      min.x <- mean(nscore$trn.table$x) + (min(scores) * small.sd)
      max.x <- mean(nscore$trn.table$x) + (max(scores) * large.sd)
      # check to see if these values are LESS extreme than the
      # initial data - if so, use the initial data.
      #print(paste('lg.sd is:',large.sd,'max.x is:',max.x,'max nsc.x is:',max(nscore$trn.table$x)))
      if(min.x > min(nscore$trn.table$x)) {min.x <- min(nscore$trn.table$x)}
      if(max.x < max(nscore$trn.table$x)) {max.x <- max(nscore$trn.table$x)}
   }
   if(tails=='equal') { # assumes symmetric distribution around the mean
      mean.x <- mean(nscore$trn.table$x)
      sd.x <- sd(nscore$trn.table$x)
      min.x <- mean(nscore$trn.table$x) + (min(scores) * sd.x)
      max.x <- mean(nscore$trn.table$x) + (max(scores) * sd.x)
      # check to see if these values are LESS extreme than the
      # initial data - if so, use the initial data.
      if(min.x > min(nscore$trn.table$x)) {min.x <- min(nscore$trn.table$x)}
      if(max.x < max(nscore$trn.table$x)) {max.x <- max(nscore$trn.table$x)}
   }
   if(tails=='none') {   # No extrapolation
      min.x <- min(nscore$trn.table$x)
      max.x <- max(nscore$trn.table$x)
   }
   min.sc <- min(scores)
   max.sc <- max(scores)
   x <- c(min.x, nscore$trn.table$x, max.x)
   nsc <- c(min.sc, nscore$trn.table$nscore, max.sc)
   
   if(draw) {plot(nsc,x, main='Transform Function')}
   back.xf <- approxfun(nsc,x) # Develop the back transform function
   val <- back.xf(scores)
   
   return(val)
}

modelSoftHard <- function(hard, soft) {
   # Given a hard (primary) dataset and a soft (secondary) vector 
   # of colocated values, return a list object with max, min, 
   # correlation, and linear model coefficients for prediction.
   soft.model <- lm(hard ~ soft)
   min.max <- predict(soft.model, data.frame(soft=c(min(soft), max(soft))))
   model.object <- list(min=min.max[1], max=min.max[2],
                   coefficients=soft.model$coeff, 
                   r=cor(hard,soft), summary=summary(soft.model))
   return(model.object)
}

write.data.frame.geoeas <- function(dat, outfile) {
   # This function writes a data frame to a GEO-EAS text file. 
   # Useful for interacting with GSLIB.
   
   # Write the header
   cat('GSLIB file created in R\n', file=outfile)
   cat(length(names(dat)), file=outfile, append=TRUE)
   cat('\n', file=outfile, append=TRUE)
   write(cbind(names(dat)), file=outfile, append=TRUE)
   write.table(dat,file=outfile, append=TRUE, sep='\t', col.names=FALSE, row.names=FALSE)
}

##########################################################
## Example functions testing these functions below here ##
##########################################################
example1.nscore  <- function() {
   # This function illustrates the use of nscore and backtr on a nonspatial synthetic example.
   par(ask=TRUE)
   rain <- runif(1000,0,7)
   rain<-(rain)^2
   rain[rain < 0] <- 0   # A rainfall-like distribution - right skewed
   hist(rain, main='Hypothetical population of precipitation estimates')   # Yep, right skewed.
   rsam <- sample(rain,100)  # Extract a sample from rain
   hist(rsam)               # Eerily similar, yet different

   rsam.sc <- nscore(rsam)  # normalize the sample
   plot(rsam.sc$trn.table$nscore, rsam.sc$trn.table$x, main='Sample rainfall normal scores')
   lines(approx(rsam.sc$trn.table$nscore, rsam.sc$trn.table$x)) # approx is a linear interp

   rsam.bak <- backtr(rsam.sc$nscore, rsam.sc, tails='separate') # Back transform the nscored transform
   cor(rsam,rsam.bak)  # Nice.

   # Nscore and backtransform the rain population data using the sample transform table
   rain.sc <- nscore(rain)
   rain.back <- backtr(rain.sc$nscore, rsam.sc, tails='none')

   cor(rain,rain.back) # about 0.99; pretty close. Note the issues in the tails.
   summary(cbind(rain,rsam,rain.back))
   par(mfrow=c(2,2))
   hist(rain, main='Hypothetical population\nof precipitation estimates')
   hist(rain.sc$nscore, main='Population scores')
   hist(rain.back, main='Back-transformed precip estimates\nusing sample transform table')
   qqplot(rain, rain.back, cex=0.6, main='Distribution Comparison')
   
   # Now, can we draw 1000 standard normal values and use rsam.sc to reproduce rain?!
   rain.sim.sc <- rnorm(1000,0,1)
   rain.sim <- backtr(rain.sim.sc, rsam.sc, tails='equal') # use the sample backtransform table
   hist(rain, main='Hypothetical population\nof precipitation estimates')
   hist(rain, main='Back-transformed simulated normal scores\nusing the sample transform table ')
   qqplot(rain, rain.sim, cex=0.6, main='Distribution Comparison')
   print(summary(cbind(rain,rsam,rain.sim)))  # summary stats....
   par(mfrow=c(1,1))
   par(ask=FALSE)
}

example2.nscore <- function() {
   # An example using gstat library simulation routines and the Meuse dataset.
   par(ask=TRUE)
   library(gstat)
   data(meuse)
   coordinates(meuse) = ~x+y
   hist(log(meuse$zinc),main='ln(Zinc): Not too Gaussian')
   v <- variogram(log(zinc)~1, meuse)
   m <- fit.variogram(v, vgm(1, "Sph", 300, 1))
   print(plot(v, model = m, main='Variogram model on ln(Zn)'))
   set.seed(131)
   data(meuse.grid)
   gridded(meuse.grid) = ~x+y
   sim.log <- krige(formula = log(zinc)~1, meuse, meuse.grid, model = m, 
           nmax = 15, beta = 5.9, nsim = 4)
   # show all 4 simulations
   print(spplot(sim.log, main='Conditional SK Simulation on ln(Zn)'))
   
   # Back transform
   sim <- sim.log
   sim$sim1<-exp(sim$sim1)
   sim$sim2<-exp(sim$sim2)
   sim$sim3<-exp(sim$sim3)
   sim$sim4<-exp(sim$sim4)
   print(spplot(sim, main='Back-transformed ln(Zn) simulations'))
   
   # Try a normal score transform instead
   meuse.zinc.sc <- nscore(meuse$zinc)
   meuse$zinc.sc <- meuse.zinc.sc$nscore
   hist(meuse$zinc.sc, main="Zinc normal score")
   v <- variogram(zinc.sc~1, meuse)
   m <- fit.variogram(v, vgm(1, "Sph", 300, 1))
   print(plot(v, model = m, main='Variogram model on normal scores Zn'))
   set.seed(131)
   sim.ns <- krige(formula = zinc.sc~1, meuse, meuse.grid, model = m, 
           nmax = 15, beta = 0, nsim = 4)
   # show all 4 ns simulations
   print(spplot(sim.ns, main='Conditional SK Simulation on Zn scores'))
   hist(sim.ns$sim1, main='Histogram of simulation 1 - pretty normal')
   sim.ns$sim1 <- backtr(sim.ns$sim1,meuse.zinc.sc, tails='separate')
   sim.ns$sim2 <- backtr(sim.ns$sim2,meuse.zinc.sc, tails='separate')
   sim.ns$sim3 <- backtr(sim.ns$sim3,meuse.zinc.sc, tails='separate')
   sim.ns$sim4 <- backtr(sim.ns$sim4,meuse.zinc.sc, tails='separate')
   print(spplot(sim.ns, main='Back-transformed Zn nscore simulations'))
   print('Original Zn data')
   print(summary(meuse$zinc))
   print('Back transformed Zn normal score simulation')
   print(summary(sim.ns$sim1))
   print('Back transformed ln(Zn) simulation')
   print(summary(sim$sim1))
   par(ask=FALSE)
}
