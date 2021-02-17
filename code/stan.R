rm(list=ls())
setwd("~/R/BAYESIAN STATISTICS/PROJECT/dati")
# R & STAN are friends!
library(rstan)
library(coda)

# for plots
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)

# this file can be patiently generated from the core.R code
# or any other wstot?.RData
load('wstot8.RData')
n <- dim(x_it)[1] 
n
# T <- 1600
# time <- 1:T
T <-100
time <- 1:T

# N.B. You have to pass to the stan function, a list containing a list
# with exactly the same element of MODEL block of the stan file

r=2;
epsilon <- x_it - MU[,r]

name = "stan/niw.stan" # "precision.stan" # "cholesky.stan" "gp_cholesky.stan" "gp.stan" "niw.stan"

if(name=='stan/niw.stan') data<-list(TT=T, n=n, epsilon=epsilon[,time], time=matrix(c(1:T),T), a0=a0, b0=b0, muL=0.5, sigmaL=1, a1 = 4)
if(name!="stan/niw.stan") data<-list(TT=T, n=n, epsilon=epsilon[,time], time=matrix(c(1:T),T), a0=a0, b0=b0, muL=0.5, sigmaL=1)

# FIT THE MODEL ------------------------------------------------------------
# Fit a model defined in the Stan modeling language
# and return the fitted result in an instance of class stanfit
# Steps performed by the function:
# (i) translate the model into c++ code
# (ii) compile the C++ code into a binary shared object
# (iii) sampling 

?stan
fit1 <- stan(file = name, 
             data = data, 
             iter = 3000,
             thin = 50,
             chains = 2, warmup = 500, 
             algorithm = 'NUTS', #adaptive mc
             diagnostic_file = 'diagnostic.txt', 
             verbose = TRUE, #it gives us something back
             seed = 111) #select a seed to initialize the algo, useful to run it again

# Arguments:
# chains = number of chains that are randomly initialized, default is 4
# (iter, thin, warmup are setted for each chain)
# algorithm = 'NUTS', 'HMC'
# sample_file, diagnostic_file = write files with samples or diagnostic data 
# save_dso = whether save the dso object or not
# control = set the parameters of the algorithm and the adaptation parameters

is(fit1)

# We can split the 3 steps into 3 different functions:
# (i) stanc: translate stan model into C++ code

step1 <- stanc(file = name, model_name = 'model_Stan')
names(step1)

# (ii) stan_model: construct a stan model tha can be used to draw samples from the model
# The C++ code is compiled into a Dynamic Shared Object.
# DSO = A dynamic shared object (DSO) is an object file that is meant to be
# used simultaneously (or shared) by multiple applications while they are executing.

step2 <- stan_model(stanc_ret = step1, verbose = FALSE)
is(step2)

# (iii) sampling: draw samples from a stan-model
step3 <- sampling(step2, 
                  data = data, chains = 2, iter = 2000, 
                  warmup = 500, thin = 1, seed = 42)
?sampling
is(step3)

# POST PROCESSING -----------------------------------------------------------

print(fit1, 
      probs = c(0.025, 0.5, 0.975), 
      par = c('sigma2', 'L'))

# Output:
# Info about the iterations and the used sampler
# Summary about parameters includes: the mean, the stardard error of the mean 
# (st_dev/sqrt(n)),
# standard deviation and quantiles, the effective sample size and the split Rhat
# All of them are computed without the warmup samples.

# The split rate R:
# One way to monitor whether a chain has converged to the equilibrium distribution 
# is to compare its behavior to other randomly initialized chains.
# The R statistic measures the ratio of the average variance of samples within each 
# chain to the variance of the pooled samples across chains; if all chains are at 
# equilibrium, these will be the same and R will be one. 
# If the chains have not converged to a common distribution, the R stat wll be greater than 1.

plot(fit1, ask = T, pars = "sigma2", ci_level = 0.95, fill_color = "blue")

# Draw the traceplot corresponding to one or more Markov chains, 
# providing a visual way to inspect
# sampling behavior and assess mixing across chains and convergence.

rstan::traceplot(fit1, pars='L', inc_warmup = TRUE)
rstan::traceplot(fit1, pars='sigma2', inc_warmup = TRUE)

# Extract the values of the chains, discarding the warmup (burnin)
# permuted = T is mixing together different chains

par <- rstan::extract(fit1, pars = c("sigma2", "L"), permuted = TRUE)

is(par)
names(par)

l <- par$L
sig2 <- par$sigma2
length(l)



# Other useful methods:
post_mean <- get_posterior_mean(fit1)
post_mean

# plot the posterior distributions 

plot_post <- fit1 %>% 
  rstan::extract(c("L", "sigma2")) %>% 
  as_tibble() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() + 
  theme(legend.position="bottom")

pairs(fit1, pars = c('L', 'sigma2'), condition = 0.5)

# CODA ------------------------------------------------------
# Thus, you can then get standard rjags/coda-style summary tables and plots.

coda_chain <- As.mcmc.list(fit1)
summary(coda_chain)

# Compute the gelman and rubin's convergence diagnostic
gelman.diag(coda_chain, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

# potential scale reduction factor is calculated for each variable in x, 
# together with upper and lower confidence limits. Approximate convergence is 
# diagnosed when the upper limit is close to 1. For multivariate chains, a multivariate 
# value is calculated that bounds above the potential scale reduction factor for 
# any linear combination of the (possibly transformed) variables

# This plot shows the evolution of Gelman and Rubin's shrink factor as the number 
# of iterations increases

gelman.plot(coda_chain[,1], confidence = 0.95)
gelman.plot(coda_chain[,10], confidence = 0.95)

# Geweke's convergence diagnostic

geweke.diag(coda_chain, frac1=0.1, frac2=0.5)

# OPTIMIZATION: obtain a point estimate by maximing the joint posterior
# from the model defined by class stanmodel

is(step2)
opt <- optimizing(step2, 
                  data = data, 
                  algorithm = 'Newton')
opt

# Arguments:
# algorithm = 'LBFGS', 'BFGS', 'Newton'
