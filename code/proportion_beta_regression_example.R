# code to run beta regression analysis looking at how proportion of deer activity
# during night, dawn, day, & dusk periods varies as a function of weather predictors

# These models are written in NIMBLE, a language derived from BUGS.
# Unlike previous programs (WinBUGS, JAGS), NIMBLE models are programmable objects
# in R, but are compiled in C++ for speed. Therefore, prior to installing NIMBLE,
# you must have Rtools installed so the code can be compiled.
# Please see the NIMBLE website for instructions: https://r-nimble.org/download

library(here)
library(tidyverse)
library(nimble)
library(coda)

setwd(here::here("data"))
  
# season: indicates which winter (1  = 2017-2018, 2 = 2018-2019)
# date: date
# pdawn: proportion of activity occuring during dawn period
# pdusk: proportion of activity occuring during dusk period
# pdiurnal: proportion of activity occurring during day period
# pnocturnal: proportion of activity occuring during night period
# CWSI: cumulative winter severity index, statewide average
# TMIN: daily minimum temperature, statewide average
# SNOW: daily snow depth, statewide average
dat <- readr::read_csv("activity_proportions_predictors_v01.csv")

# inspect
glimpse(dat)

# reformat data so each year of data has a separate column
act <- dat %>% 
  dplyr::mutate_at(vars(CWSI:SNOW), function(x) as.vector(scale(x))) %>% 
  group_by(season) %>% 
  mutate(id = row_number()) %>%
  dplyr::select(-date) %>% 
  pivot_wider(names_from = season, values_from = pdawn:SNOW)

# following?
glimpse(act)
  
# this is a key that will help us track through the different models
key <- tibble(
  index = 1:20, 
  response = sort(rep(c("pnocturnal", "pdawn", "pdiurnal", "pdusk"), 5)), # which period 
  name = rep(paste0("mod", 1:5), 4),
  npreds = rep(c(2, 2, 1, 1, 1), 4), # how many predictors
  pred1 = rep(c("TMIN", "TMIN", "TMIN", "SNOW", "CWSI"), 4), # which predictor?
  pred2 = rep(c("CWSI", "SNOW", NA, NA, NA), 4)) # which predictor? if there's more than 1 in model

# list where we'll send the results
result_list <- list()
for(i in 1:nrow(key)){
  
  # slightly separate model structure for 1 versus 2 predictors
  if(key[[i, "npreds"]] == 1){
    
    # package up data for nimble model
    data <- list(response = unname(as.matrix(dplyr::select(act,
                                                           starts_with(key[[i, "response"]])))),
                 pred1 = unname(as.matrix(dplyr::select(act,
                                                        starts_with(key[[i, "pred1"]])))))
    
    # package up constants; for indexing purposes
    constants <- list(pdays = nrow(data$response),
                      nseason = ncol(data$response))
    
    # here's the code for the nimble model
    modelCode <- nimble::nimbleCode( {
      
      # likelihood
      for(i in 1:pdays){
        for(j in 1:nseason){
          response[i, j] ~ dbeta(alpha[i,j], beta[i,j])
          # link function
          alpha[i,j] <- mu[i,j]*phi
          beta[i,j] <- (1 - mu[i,j])*phi
          # linear predictor with random effect of year
          logit(mu[i,j]) <- b1*pred1[i,j] + eps[j]
        }}
      
      # priors
      
      # random effect of year; gives each season a different intercept
      for(j in 1:nseason){
        eps[j] ~ dnorm(0, tau)
      }
      
      tau ~ dgamma(0.1, 0.1)
      phi ~ dgamma(0.1, 0.1)
      b1 ~ dlogis(0, 1)
      
    }) # end of model code
    
    # parameters we're gonna monitor
    keepers <- c("b1",
                 "tau",
                 "eps")
    
    # initial values
    inits <- function(){
      list(phi = rgamma(1, 0.1, 0.1),
           b1 = runif(1, -1, 1),
           tau = rgamma(1, 0.1, 0.1),
           eps = rnorm(2, 0, 1))
    }
    
    # getting things ready to go
    mod <- nimble::nimbleModel(code = modelCode,
                               constants = constants,
                               data = data,
                               inits = inits())
    
    conf <- configureMCMC(mod,
                          enableWAIC = TRUE)
    
    conf$addMonitors(keepers)
    
    rmcmc <- buildMCMC(conf)
    
    cmodel <- compileNimble(mod)
    
    cmcmc <- compileNimble(rmcmc, project = cmodel)
    
    # run it
    out <- runMCMC(cmcmc,
                   niter = 5000,
                   nburnin = 4000,
                   nchains = 3,
                   samples = TRUE,
                   WAIC = TRUE)
    
    # put results from all chains into one array
    samples <- do.call(rbind, out$samples)
    
    # convert to mcmc for calculating rhat
    samples_mcmc <- as.mcmc(lapply(out$samples, mcmc))
    
    # assess convergence with rhat
    rhat <- as_tibble(coda::gelman.diag(samples_mcmc)$psrf,
                      rownames = "parameter") %>%
      rename(rhat = `Point est.`) %>%
      dplyr::select(parameter, rhat) %>%
      dplyr::filter(grepl("^b", parameter)) %>%
      add_column(pred = key[[i, "pred1"]])
    
    # package up the result as a dataframe 
    result_list[[i]] <- as_tibble(samples) %>% 
      dplyr::select(., starts_with("b")) %>%
      pivot_longer(cols = 1,
                   names_to = "parameter",
                   values_to = "values") %>%
      group_by(parameter) %>%
      summarise(mean = mean(values),
                sd = sd(values),
                lower = quantile(values, c(0.025)),
                upper = quantile(values, c(0.975))) %>%
      full_join(rhat) %>%
      add_column(waic = out$WAIC,
                 response = key[[i, "response"]],
                 index = key[[i, "index"]]) %>%
      dplyr::select(index, response, parameter, pred, mean:rhat, waic)
    
    
    # same thing over again, but for mdoels with multiple predictors
  } else {
    
    data <- list(response = unname(as.matrix(dplyr::select(act,
                                                           starts_with(key[[i, "response"]])))),
                 pred1 = unname(as.matrix(dplyr::select(act,
                                                        starts_with(key[[i, "pred1"]])))),
                 pred2 = unname(as.matrix(dplyr::select(act,
                                                        starts_with(key[[i, "pred2"]])))))
    
    constants <- list(pdays = nrow(data$response),
                      nseason = ncol(data$response))
    
    modelCode <- nimble::nimbleCode( {
      
      for(i in 1:pdays){
        for(j in 1:nseason){
          response[i, j] ~ dbeta(alpha[i,j], beta[i,j])
          alpha[i,j] <- mu[i,j]*phi
          beta[i,j] <- (1 - mu[i,j])*phi
          logit(mu[i,j]) <- b1*pred1[i,j] + b2*pred2[i, j] + b3*pred1[i,j]*pred2[i, j] + eps[j]
        }}
      
      for(j in 1:nseason){
        eps[j] ~ dnorm(0, tau)
      }
      
      tau ~ dgamma(0.1, 0.1)
      
      phi ~ dgamma(0.1, 0.1)
      b1 ~ dlogis(0, 1)
      b2 ~ dlogis(0, 1)
      b3 ~ dlogis(0, 1)
      
    })
    
    keepers <- c("b1",
                 "b2",
                 "b3",
                 "tau",
                 "eps")
    
    inits <- function(){
      list(phi = rgamma(1, 0.1, 0.1),
           b1 = runif(1, -1, 1),
           b2 = runif(1, -1, 1),
           b3 = runif(1, -1, 1),
           tau = rgamma(1, 0.1, 0.1),
           eps = rnorm(2, 0, 1))
    }
    
    mod <- nimble::nimbleModel(code = modelCode,
                               constants = constants,
                               data = data,
                               inits = inits())
    
    conf <- configureMCMC(mod,
                          enableWAIC = TRUE)
    
    conf$addMonitors(keepers)
    
    rmcmc <- buildMCMC(conf)
    
    cmodel <- compileNimble(mod)
    
    cmcmc <- compileNimble(rmcmc, project = cmodel)
    
    out <- runMCMC(cmcmc,
                   niter = 5000,
                   nburnin = 4000,
                   nchains = 3,
                   samples = TRUE,
                   WAIC = TRUE)
    
    samples <- do.call(rbind, out$samples)
    
    samples_mcmc <- as.mcmc(lapply(out$samples, mcmc))
    
    rhat <- as_tibble(coda::gelman.diag(samples_mcmc)$psrf,
                      rownames = "parameter") %>%
      rename(rhat = `Point est.`) %>%
      dplyr::select(parameter, rhat) %>%
      dplyr::filter(grepl("^b", parameter)) %>%
      add_column(pred = c(key[[i, "pred1"]],
                          key[[i, "pred2"]],
                          paste(key[[i, "pred1"]], key[[i, "pred2"]], sep = "*")))
    
    result_list[[i]] <- as_tibble(samples) %>%
      dplyr::select(., starts_with("b")) %>%
      pivot_longer(cols = 1:3,
                   names_to = "parameter",
                   values_to = "values") %>%
      group_by(parameter) %>%
      summarise(mean = mean(values),
                sd = sd(values),
                lower = quantile(values, c(0.025)),
                upper = quantile(values, c(0.975))) %>%
      full_join(rhat) %>%
      add_column(waic = out$WAIC,
                 response = key[[i, "response"]],
                 index = key[[i, "index"]]) %>%
      dplyr::select(index, response, parameter, pred, mean:rhat, waic)
    
  }
  
}

# one sleek data frame of results
result_df <- bind_rows(result_list)


# save results if you want
# setwd(here::here("results"))
# write_csv(result_df, "beta_regression_result_summary_all.csv")
