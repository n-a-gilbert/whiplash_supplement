library(here)
library(tidyverse)
library(nimble)
library(parallel)
library(coda)

setwd(here::here("data"))
load("whiplashOccModelDataConstants_v02.RData")

#### DATA ####
# the list named "data" contains 5 matrices 

str(data)

# the first is called "covariate"
# these are the static landscape predictors, measured within a 1100-m buffer of camera locations
# rows are sites, columns are the different predictors
# they are, in order: aspect, deciduous, deciduous^2, elevSD, open, open^2, and shelter
# aspect is the proportion of south-facing slopes, calculated from a digital elevation model
# deciduous is the proportion of deciduous forest in the landscape
# deciduous^2 is just deciduous squared (for quadratic effect)
# elevSD is topographic relief, calculated as the standard deviation of elevation
# open is the proportion of "open" (cropland and grassland) habitat
# open^2 is just open squared
# shelter is the proportion of shelter habitat, defined as coniferous, mixed forest and wooded wetlands
# all predictors are already scaled

# inspect distribution of predictors
data$covariate %>% 
  as_tibble() %>% 
  setNames(., c("relief", "deciduous", "deciduous2", "open", "open2", "conifer", "conifer2")) %>% 
  mutate(id = row_number()) %>% 
  pivot_longer(relief:conifer2, names_to = "predictor", values_to = "value") %>% 
  ggplot(., aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(~predictor, scales = "free")

# the second matrix contains daily minimum temperature (tmin) and cumulative winter severity index (cwsi)
# rows are sites, columns days, matrix slices are tmin/cwsi
# columns 1-86 are the winter of 2017-2018, 87-172 are 2018-2019
# values are already scaled
# example: this shows tmin for one site over the course of the first winter
plot(data$weather[1,1:86,1,1]) 

# here's average tmin of all sites over the course of the first winter
plot(apply(data$weather[,,1,1], 2, mean)[1:86])

# here's average cwsi of all sites over the course of the winter - notice how it accumulates
plot(apply(data$weather[,,1,2], 2, mean)[1:86])

# the fourth matrix, "zmat" is for the spatial random effect. 
# it provides, for each camera location (row), the distance to 20 knot locations (columns)
# values are already scaled

# the fifth matrix, "period", indicates for each camera location (row), day(column),
# whether a replicate survey period is 00:00-08:00 hrs (0), 
# 08:00-16:00 hrs (1), or 16:00-00:00 hrs(0)
str(data$period)
# for one site on one day
plot(data$period[1,1,,1])

# finally, the sixth matrix, "y", is the deer detection record
# rows are sites, columns days, 3rd matrix slices are replicate survey periods, 4th matrix slice is year
# values can be 0 (deer not detected), 1 (deer detected), or NA (camera not active)
str(data$y)

# here's a detection history for a single site. this camera was active only in the first winter
data$y[1,,,]

#### CONSTANTS ####
# the object "constants" is a list of constants in the model
# these are used for indexing in the model
str(constants)
# nsite - number of camera locations
# nday - total number of days
# nsurvey - number of replicate surveys in one day
# nyear - number of years
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
  # sd_bs is the precision
  for(i in 1:nknots){
    gamma[i] ~ dnorm(0, sd_bs[i])
    sd_bs[i] ~ dgamma(1, 2)
  }
  
  # prior for the period predictor in the observation submodel
  a1 ~ dlogis(0, 1)
  
  # prior for random intercept for year
  for(k in 1:nyear){
    eps_psi[k] ~ dnorm(0, sd_eps_psi)
  }
  
  # and precision
  sd_eps_psi ~ dgamma(1, 2)
  
  # prior for site/survey random effect
  for(i in 1:nsite){
    for(j in 1:nday){
      for(k in 1:nyear){
        eps_p[i, j, k] ~ dnorm(0, sd_eps_p)
      }
    }
  }
  
  # and its precision
  sd_eps_p ~ dgamma(1, 2)
  
  # .............................................................
  # LIKELIHOOD
  # .............................................................
  
  # state model - ecological process
  for(i in 1:nsite){
    for(j in 1:nday){
      for(k in 1:nyear){
        
        z[i, j, k] ~ dbern(psi[i, j, k])
        
        logit(psi[i, j, k]) <- inprod(covariate[i,1:nland], b[1:nland]) + # landscape patterns
          inprod(weather[i, j, k, 1:nweather], b[wstart:wend]) + # tmin & wsi
          b[tw]*weather[i,j, k, 1]*weather[i, j, k, 2] + # tmin*wsi interaction
          inprod(covariate[i, 1:nland]*weather[i, j, k, 1], b[tintstart:tintend]) + # landscape*tmin interactions
          inprod(covariate[i, 1:nland]*weather[i, j, k, 2], b[wintstart:wintend]) + # landscape*wsi interactions
          inprod(zmat[i, 1:nknots], gamma[1:nknots]) + # spatial random effect
          eps_psi[k] # random intercept for year
        
        # detection submodel - observation process
        for(m in 1:nsurvey){
          
          y[i, j, m, k] ~ dbern(z[i, j, k]*p[i, j, m, k])
          logit(p[i, j, m, k]) <- a1*period[i, j, m, k] + eps_p[i, j, k]
        } # nsurvey 
      } # nyear
    } # nday
  } # nsite
}) # model


# initial values for running model
inits <- function() {list(p = array(0.5, dim = c(constants$nsite,
                                                 constants$nday,
                                                 constants$nsurvey,
                                                 constants$nyear)),
                          z = array(1, dim = c(constants$nsite, constants$nday, constants$nyear)),
                          y = array(1, dim = c(constants$nsite, constants$nday, constants$nsurvey, constants$nyear)),
                          b = runif(constants$npred, -1, 1),
                          a1 = runif(1, -1, 1), 
                          gamma = runif(constants$nknots, -1, 1),
                          sd_bs = runif(constants$nknots, 0, 2),
                          eps_p = array(rnorm(constants$nsite*constants$nday*constants$nyear, 0, 2),
                                        dim = c(constants$nsite, constants$nday, constants$nyear)),
                          sd_eps_p = runif(1, 0, 2),
                          eps_psi = rnorm(constants$nyear, 0, 2),
                          sd_eps_psi = runif(1, 0, 2))}

# parameters to monitor
keepers <- c("b", "a1", "gamma", "eps_psi")

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
  parameter = c("a1",
                paste0("b", "[", 1:constants$npred, "]"),
                paste0("gamma[", 1:constants$nknots, "]")),
  name = c("period",
           "relief", "deciduous", "deciduous2", "open", "open2", "conifer", "conifer2",
           "tmin", "wsi", "tmin*wsi",
           "relief*tmin", "deciduous*tmin", "deciduous2*tmin", "open*tmin", "open2*tmin", "conifer*tmin", "conifer2*tmin",
           "relief*wsi", "deciduous*wsi", "deciduous2*wsi", "open*wsi", "open2*wsi", "conifer*wsi", "conifer2*wsi",
           paste0("spline", 1:constants$nknots)))

outmcmc <- as.mcmc(out)

# calculate rhat, should be < 1.1 to be considered converged
rhat <- as_tibble(coda::gelman.diag(outmcmc[,1:47])[[1]],
                  rownames = "parameter") %>% 
  setNames(., c("parameter", "rhat", "upper")) %>% 
  dplyr::select(-upper)

samps <- do.call(rbind, out)

# results bundled up nicely
res <- as.data.frame(samps[, 1:47]) %>% 
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