library(here)
library(tidyverse)
library(nimble)
library(parallel)
library(coda)

setwd(here::here("data"))
load("whiplashOccModelDataConstants.RData")

#### DATA ####
# the list named "data" contains 6 matrices 

str(data)

# the first is called "covariate"
# these are the static landscape predictors, measured within a 1100-m buffer of camera locations
# rows are sites, columns are the different predictors
# they are, in order: aspect, canopy, canopy^2, elevSD, and shelter
# aspect is the proportion of south-facing slopes, calculated from a digital elevation model
# canopy is the proportion of canopy cover, mean value for the landscape
# canopy^2 is just canopy squared (for quadratic effect)
# elevSD is topographic relief, calculated as the standard deviation of elevation
# shelter is the proportion of shelter habitat, defined as coniferous, mixed forest and wooded wetlands
# all predictors are already scaled

# inspect distribution of predictors
data$covariate %>% 
  as_tibble() %>% 
  setNames(., c("aspect", "canopy", "canopy2", "elevSD", "shelter")) %>% 
  mutate(id = row_number()) %>% 
  pivot_longer(aspect:shelter, names_to = "predictor", values_to = "value") %>% 
  ggplot(., aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(~predictor)

# the second matrix contains daily minimum temperature (tmin) and cumulative winter severity index (cwsi)
# rows are sites, columns days, matrix slices are tmin/cwsi
# columns 1-86 are the winter of 2017-2018, 87-172 are 2018-2019
# values are already scaled

# example: this shows tmin for one site over the course of the first winter
plot(data$weather[1,1:86,1]) 

# here's average tmin of all sites over the course of the first winter
plot(apply(data$weather[,,1], 2, mean)[1:86])

# here's average cwsi of all sites over the course of the winter - notice how it accumulates
plot(apply(data$weather[,,2], 2, mean)[1:86])


# the third matrix, "year", indicates weather a column (date) is in 2017-2018 or 2018-2019
# 0 = 2017-2018, 1 = 2018-2019
plot(data$year[1,])

# the fourth matrix, "zmat" is for the spatial random effect. 
# it provides, for each camera location (row), the distance to 20 knot locations (columns)
# values are already scaled

# the fifth matrix, "period", indicates for each camera location (row), day(column),
# whether a replicate survey period is 00:00-08:00 hrs (0), 
# 08:00-16:00 hrs (1), or 16:00-00:00 hrs(0)
str(data$period)
# for one site on one day
plot(data$period[1,1,])

# finally, the sixth matrix, "y", is the deer detection record
# rows are sites, columns days, matrix slices are replicate survey periods
# values can be 0 (deer not detected), 1 (deer detected), or NA (camera not active)

str(data$y[1,,])

# here's a detection history for a single site. this camera was active only in the first winter
data$y[1,,]

#### CONSTANTS ####
# the object "constants" is a list of constants in the model
# these are used for indexing in the model
str(constants)
# nsite - number of camera locations
# nday - total number of days
# nsurvey - number of replicate surveys in one day
# nknots - number of knots for spatial random effect
# nweather  - number of weather predictors
# wstart - index for where weather predictors start in linear predictor of model
# wend - index for where weather predictors end in linear predictor of model
# tw - index for tmin*cwsi interaction in linear predictor
# tintstart - index for where tmin*landscape interactions start in linear predictor
# tintend - index for where tmin*landscape interactions end in linear predictor
# wintstart - index for where cwsi*landscape interactions start in linear predictor
# wintend - index for where cwsi*landscape interactions end in linear predictor
# npred - total number of predictors in linear predictor

# .......................................................................
#### MODEL CODE ####
# .......................................................................
  
myCode <- nimbleCode({
  
  # .............................................................
  # PRIORS
  # .............................................................
  
  # priors for the predictor coefficients
  # a logistic prior is just a normal distribution with fatter tails
  # and more robust to logistic transformations
  # a logistic(0,1) prior is weakly informative
  # with mean = 0, if a predictor isn't related to the response, the coefficient will shrink to 0
  for(i in 1:npred){
    b[i] ~ dlogis(0, 1)
  }
  
  # priors for the spatial random effect penalized spline
  # gamma is just a normal linear predictor
  # sd_bs is the precicsion
  for(i in 1:nknots){
    gamma[i] ~ dnorm(0, sd_bs[i])
    sd_bs[i] ~ dunif(0, 2)
  }
  
  # prior for the period predictor in the observation submodel
  a1 ~ dlogis(0, 1)
  
  # prior for random effect in observation submodel 
  for(i in 1:nsite){
    for(j in 1:nday){
      eps_p[i, j] ~ dnorm(0, sd_eps_p)
    }
  }
  
  # and its precision
  sd_eps_p ~ dunif(0, 2)
  
  # .............................................................
  # LIKELIHOOD
  # .............................................................
  
  # state model - ecological process
  for(i in 1:nsite){
    for(j in 1:nday){
      
      z[i, j] ~ dbern(psi[i, j])
      
      logit(psi[i, j]) <- inprod(covariate[i,1:nland], b[1:nland]) + # landscape patterns
        inprod(weather[i, j, 1:nweather], b[wstart:wend]) + # tmin & cwsi
        b[tw]*weather[i,j,1]*weather[i, j, 2] + # tmin*cwsi interaction
        inprod(covariate[i, 1:nland]*weather[i, j, 1], b[tintstart:tintend]) + # landscape*tmin interactions
        inprod(covariate[i, 1:nland]*weather[i, j, 2], b[wintstart:wintend]) + # landscape*cwsi interactions
        b[npred]*year[i, j] + # separate intercept for year
        inprod(zmat[i, 1:nknots], gamma[1:nknots]) # SRE
      
      # detection submodel - observation process
      for(k in 1:nsurvey){
        
        y[i, j, k] ~ dbern(z[i, j]*p[i, j, k])
        logit(p[i, j, k]) <- a1*period[i, j, k] + eps_p[i, j]
        
      } # nsurvey
    } # nday
  } # nsite
}) # model

# initial values for running model
inits <- function() {list(p = array(0.5, dim = c(constants$nsite,
                                                 constants$nday,
                                                 constants$nsurvey)),
                          z = array(1, dim = c(constants$nsite, constants$nday)),
                          y = array(1, dim = c(constants$nsite, constants$nday, constants$nsurvey)),
                          b = runif(constants$npred, -1, 1),
                          a1 = runif(1, -1, 1), 
                          gamma = runif(constants$nknots, -1, 1),
                          sd_bs = runif(constants$nknots, 0, 2),
                          eps_p = array(rnorm(constants$nsite*constants$nday, 0, 2),
                                        dim = c(constants$nsite, constants$nday)),
                          sd_eps_p = runif(1, 0, 2))}


# parameters to monitor
keepers <- c("b", "a1", "gamma")


# OPTION 1 - NOT RECOMMENDED - RUN CHAINS ONE AFTER ANOTHER
# this is a big model and takes a long time to run - better to 
# do option #2 and run in parallel 

# model <- nimbleModel(code = myCode,
#                      constants = constants,
#                      data = data,
#                      inits = inits())

# model$initializeInfo()

# Cmodel <- compileNimble(model)

# modelConf <- configureMCMC(model)

# modelConf$addMonitors(keepers)

# modelMCMC <- buildMCMC(modelConf)

# CmodelMCMC <- compileNimble(modelMCMC, project = model)
# 
# out1 <- runMCMC(CmodelMCMC,
#                 nburnin = 9000,
#                 niter = 10000,
#                 nchains = 3)



# OPTION 2 (RECOMMENDED) - RUN CHAINS IN PARALLEL
# again, this is a big model, and especially if
# running chains in parallel, make sure you're equipped for the task
# I'm running this on a lab server
# Intel 4116 CPU @ 2.10 GHz, 768 GB RAM, 48 cores
# With that, it takes roughly 20 hours to run 3 chains for 10,000 iterations
# but, a laptop will probably choke trying to run this sucker

nc <- 3
nb <- 9000
ni <- nb + 1000

start <- Sys.time()

# makes a cluster with a node for each chain to run
cl <- makeCluster(nc)

parallel::clusterExport(cl, c("myCode",
                              "inits",
                              "data",
                              "constants",
                              "keepers",
                              "nb",
                              "ni"))

for(j in seq_along(cl)){
  set.seed(j)
  init <- inits()
  clusterExport(cl[j], "init")
}

out <- clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  model <- nimbleModel(code = myCode,
                       name = "myCode",
                       constants = constants,
                       data = data,
                       inits = init)
  
  Cmodel <- compileNimble(model)
  modelConf <- configureMCMC(model)
  modelConf$addMonitors(keepers)
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  out1 <- runMCMC(CmodelMCMC,
                  nburnin = nb,
                  niter = ni)
  return(as.mcmc(out1))
})

# you'll want to shut down the cluster once you're done, but only do so 
# once you're sure model is coverged and you don't need to extend chains
# stopCluster(cl)
end <- Sys.time()
end - start

# key for parameter names
key <- tibble(
  parameter = c("a1", paste0("b", "[", 1:constants$npred, "]"),
                paste0("gamma[", 1:constants$nknots, "]")),
  name = c("period", "aspect", "canopy", "canopy2", "elevSD", "shelter",
           "tmin", "wsi", "tmin*wsi",
           "aspect*tmin", "canopy*tmin", "canopy2*tmin", "elevSD*tmin", "shelter*tmin",
           "aspect*wsi", "canopy*wsi", "canopy2*wsi",  "elevSD*wsi", "shelter*wsi",
           "year",
           paste0("spline", 1:constants$nknots)))

outmcmc <- as.mcmc(out)

# calculate rhat, should be < 1.1 to be considered converged
rhat <- as_tibble(coda::gelman.diag(outmcmc[,1:40])[[1]],
                  rownames = "parameter") %>% 
  setNames(., c("parameter", "rhat", "upper")) %>% 
  dplyr::select(-upper)

samps <- do.call(rbind, out)

# results bundled up nicely
res <- as.data.frame(samps[, 1:40]) %>% 
  pivot_longer(`a1`:`gamma[20]`, names_to = "parameter", values_to = "value") %>% 
  group_by(parameter) %>% 
  mutate(mean = mean(value),
         lower = quantile(value, c(0.025)),
         upper = quantile(value, c(0.975))) %>% 
  dplyr::select(-value) %>% 
  distinct(.) %>% 
  full_join(key) %>% 
  dplyr::select(parameter, name, mean:upper) %>% 
  full_join(rhat)

# If model is not converged, you can extend the chains

start <- Sys.time()
out2 <- clusterEvalQ(cl, {
  outt2 <- runMCMC(CmodelMCMC,
                   niter = 10000,
                   nburnin = 9000)
  return(as.mcmc(outt2))
})
end <- Sys.time()
end - start
# stop that cluster once you're done...
# stopCluster(cl)