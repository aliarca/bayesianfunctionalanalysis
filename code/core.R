# Clear plots
if(!is.null(dev.list())) dev.off()
# Clean workspace
rm(list=ls())
# Clear console
cat("\014") 

needed_packages  <- c("fda.usc","fda","fields","LaplacesDemon","matlab","mvtnorm", "Rmpfr", "LearnBayes")
new_packages  <- needed_packages[!(needed_packages %in%installed.packages ()[, "Package"])]
if (length(new_packages))  install.packages(new_packages)
lapply(needed_packages , require , character.only = TRUE)


load("data/Codice_analisi_funzionale.RData")

load('samplers/niw.RData')
load('samplers/nor.RData')
load('samplers/gib.RData')
load("samplers/update.RData")

load("functions/calcolo.loglike.RData")
load("functions/bieffe.RData")
load("functions/bieffe.power.RData")


###########################################################################
######################### PRELIMINARY OPERATIONS ##########################

# Building the 1-0  matrix to keep track of patinets and signals
a = rep(1,8)             #the patient has all the signals
# order: ausxSL  audxSL  frsxSL  frdxSL ausxML  audxML  frsxML  frdxML
a6  = c(1,1,1,1,0,1,0,1) #missing au&fr sxML
a10 = c(1,1,1,0,1,1,1,1) #missing frsxML
a11 = c(1,1,1,1,1,1,1,0) #missing frdxML
a14 = c(1,1,1,1,1,0,1,0) #missing au&fr dxML 
a15 = c(1,1,1,1,0,1,0,1) #missing au&fr sxML
a22 = c(1,1,1,1,0,0,0,0) #missing all ML
a23 = c(1,1,1,1,0,0,0,0) #missing all ML
a25 = c(1,0,1,0,1,1,1,1) #missing au&fr dxSL
check_signal = rbind(a,a,a,a,a,a6,a,a,a,a10,a11,a,a,a14,a15,a,a,a,a,a,a,a22,a23,a,a25,a)

# lists in which we will save the results
BF=list()
log2BF=list()
BF.power=list()
log2BF.power=list()
set.seed(111)
# for everey signal we do the entire analysis
for (n_sign in 1:8){
  print(n_sign)
  ###########################################################################
  ################################ THE DATA #################################
  
  # number of patients with the specified signal 
  n_tot <- sum(check_signal[,n_sign])
  indexes <- which(check_signal[,n_sign]==1)
  # functional object containing the signal
  fda_signal <- f.data[[n_sign]] 
  
  ###########################################################################
  ############################## THE BEGINNING ##############################
  
  # Disclaimer, here we directly apply the results of preprocessing 
  # (splines representation and PCA)
  # You can find all the discussion in files "Functional.R" and "fPCA.R"
  #n_tot <- 26 # number of patients
  L <- 15 # number of PCA components
  T <- 1600
  
  ###########################################################################
  ########## F-PCA WITH L COMPONENT, NOT CENTRED (study the mean!) ##########
  
  # functional representation of the mean vector, to ease sampling
  
  n <- dim(fda_signal$data)[1] # number of patients, usually between 22 and 26
  signal <- fda_signal$data    # pure data [n x 1600]
  x_it <- signal
  signal <- t(signal)
  time <- 1:T                  # abscissa
  m_opt <- 377                 # optimal n. basis
  basis_opt <- create.bspline.basis(rangeval=c(0,1600), nbasis=m_opt)
  data_L.fd <- Data2fd(y = signal, argvals = time,basisobj = basis_opt)
  pca_L <- pca.fd(data_L.fd,nharm=L,centerfns=FALSE)
  
  
  # betas and b_t #################################################
  # we set our bases b_t and our betas (from the not scaled fPCA)
  b_t <- pca_L$harmonics
  betas <-  pca_L$scores
  
  colMeans(pca_L$scores)
  beta_mean <- colMeans(betas)
  
  
  # > ALL PATIENTS < < < < < < < < < < < < < < < < < < ####
  ###########################################################################
  ########################### Model for the MEAN ############################
  
  # '- - - NIW - - -' ####
  # we have mu in R^L, where L is the number of bases
  # (mu,Sigma) ~ NIW (m0, k0, nu0, Lambda0) 
  # <=> 
  # mu|Sigma ~ N (m0, Sigma/k0) where m0 in R^L, k0>0
  # Sigma ~ IW (nu0, Lambda0) where nu0>0 and nu>L-1, Lambda0 in R^(LxL)
  
  print("Sample total model")
  
  niter <- 5000 # number of Mc samples I run from my posterior, which equals the niter of the Gibbs
  print("MU")
  niw <- sampling_mean_niw(niter, betas, beta_mean, b_t, data_L.fd, L, T, n)
  MU = niw$MU
  k0 <- niw$k0
  
  # # '- - - Normal - - -' ####
  # # A DIFERENT POSSIBLE STRATEGY OF MODELLING THE MEAN TERM:
  # # NORMAL NORMAL instead of NIW
  # 
  # # Due to the difficulty in extimating the prior parameters for the NIW distribution,
  # # We can assume a Normal Normal distribution, so that we can leave a general loose variance
  # # for the prior distribution of mu, without making useless risky assumptions
  # # on the distribution of this variance
  
  # print("MU_norm")
  # niter <- 5000 # number of Mc samples I run from my posterior, which equals the niter of the Gibbs 
  
  # MU_norm <- sampling_mean_normal(niter, mun, Sigman, mu_t, b_t, L, T, n)
  
  # sigma_2 and phi_t #############################################
  # prior sigma_2 is IG (a0, b0)
  
  # how do I choose my initial parameters? I can start reasoning about the mean of the two distributions
  # E[sigma_2] = b0 / (a0-1), if a0>1 => I certainly want a0>1
  # Now, this sigma_2 represents the common factor of the variance...
  # So I compute the variance in each point, I take its mean, and I tune a0 and b0 so that I obtain it back
  # I may start with a0=3 and b0 around twice the mean variance value => E[sigma_2] = b0
  # Then I hope that by updating using data, I can reach a good modelling for it
  
  # NOTE THAT:
  # phi(t) indep of phi(v) => I'll call phi(t) as phi_t for every (t,v) in (1:1600)x(1:1600) and t!=v
  # sigma_2|{phi_t}_t,{x_it}_i_t,{mu(t)}_t ~ IG( a,b )
  
  # we run a Gibbs Sampler to simulate form the joint posterior of sigma_2 and phi
  # parameters: sigma_2 and phi_t for t = 1:T => #param = n+1 = 1601
  
  # burnin
  burnin <- 2000
  niter <- niter
  print("sigma^2 and phi")
  
  # '- - - NIW - - -' ####
  gib <- gibbs(niter, burnin, T, n, x_it, MU)
  sigma_2_result <- gib$sigma
  phi_result <- gib$phi
  # initialization of the two parameters of sigma^2 distribution; 
  a0 <- gib$a0
  b0 <- gib$b0
  
  # # '- - - Normal - - -' ####
  # gib_norm <- gibbs(niter, burnin, T, a0, b0, c0, d0, n, x_it, MU_norm)
  # sigma_2_result_norm <- gib_norm$sigma
  # phi_result_norm <- gib_norm$phi
  
  
  # some plots ####
  
  matplot(signal, type='l', main = paste('all patients, singal',n_sign))
  
  matplot(MU, type='l', main = paste('MU GLOBAL MODEL, signal',n_sign))
  
  plot(sigma_2_result, main = paste('sigma^2 ALL, signal',n_sign))
  matplot(sigma_2_result,type='l', title=title(c('sigma_2_with_Burn_ALL traceplot',paste('signal',n_sign))))
  hist(sigma_2_result,probability = T)
  
  indici = c(1,3,333,73,65,22,768,53,1121,984,1600,823,1233,1453,653,12,90,25,679,123,124,125)
  for (ind in indici)
    matplot(phi_result[,ind], type='l',title=title(c(paste('phi_ALL(',ind,')',sep=''),paste('signal',n_sign))))
  
  
  # # CREDIBLE INTERVALS FOR THE MEAN -------------- ####
  # #tapply(norm_MU, 1,quantile(x,probs = c(0.025, 0.5, 0.975)))
  MU_quantiles <- t(apply(MU,1,function(x) quantile(x, probs = c(0.025, 0.5, 0.975))))
  #norm_MU_quantiles <- t(apply(MU_norm, 1, function(x) quantile(x,probs = c(0.025, 0.5, 0.975))))
  
  matplot(signal,type='l', col=rep('Dark Grey',n), main = paste('all patients CREDIBLE BANDS for MU, singal',n_sign))
  matlines(MU_quantiles,lty =1, lwd=2, col='magenta') 
  # matlines(norm_MU_quantiles, lty=2, lwd=1, col='blue')
  # legend('top', legend=c("MU form NIW", "MU from N"),
  #        col=c('magenta','blue'), lty=1:2, cex=0.6,lwd=c(2,2))
  # #MU_quantiles
  # #norm_MU_quantiles
  # 
  # #we compare the 2 sets of credible intervals
  # # we observe that thay are almost equivalent
  # # we are satisfied with our results, as our sampled curves for mu
  # # have similar shapes but have a variance that is large enough.
  # # ------------------------------------------------------------- #
  
  # plot phi(t) for every t, simulated at iteration r
  # r <- 1000
  #matplot(phi_result[r,], type='l')          
  #matplot(t(phi_result))
  
  
  # Likelihood computation ####
  logTOT = calcolo.loglike(niter, burnin, T, MU, n, sigma_2_result, phi_result, signal)
  
  #### > > GOSE 1 < < < < < < < < < < < < < < < < < < ####
  ###########################################################################
  ############################## THE BEGINNING ##############################
  GOSE <- data[indexes,"GOSE"]
  
  n1 <- sum(GOSE==1) # number of patients with GOSE==1
  x_it1 <- fda_signal$data[GOSE==1,] 
  
  # betas ###################################################################
  #dim(cbind(ones(n1,1), pca_L$scores))
  betas1 <- pca_L$scores[GOSE==1,]
  
  dim(betas1)
  colMeans(pca_L$scores)
  beta_mean1 <- colMeans(betas1)
  
  ###########################################################################
  ########################### Model for the MEAN ############################
  
  print("Sample GOSE=1 model")
  
  # '- - - NIW - - -' ####
  print("MU1")
  niw1 <- sampling_mean_niw(niter, betas1, beta_mean1, b_t, data_L.fd, L, T, n1)
  k01 <- niw1$k0
  MU1 <- niw1$MU
  
  # # '- - - Normal - - -' ####
  # print("MU1_norm")
  # MU1_norm <- sampling_mean_normal(niter, betas1, beta_mean1, b_t, data_L.fd, L, T, n1)
  
  ###########################################################################
  ######################## Model for the Covariance #########################  
  
  niter <- niter
  burnin <- burnin
  print("sigma^2 and phi 1")
  
  # '- - - NIW - - -' ####
  gib1 <- gibbs(niter, burnin, T, n1, x_it1, MU1)
  sigma_2_result1 <- gib1$sigma
  phi_result1 <- gib1$phi

  a01 <- gib1$a0
  b01 <- gib1$b0
  
  # # '- - - Normal - - -' ####
  # gib1_norm <- gibbs(niter, burnin, T, x_it1, MU1_norm)
  # sigma_2_result1_norm <- gib1_norm$sigma
  # phi_result1_norm <- gib1_norm$phi
  
  
  # some plots ####
  matplot(signal[GOSE==1,], type='l', main = paste('patients with GOSE1, signal',n_sign))
  
  matplot(MU1, type='l', main = paste('MU GOSE1,',n_sign))
  
  
  # # ------------- CREDIBLE INTERVALS FOR THE MEAN -------------- #
  # #tapply(norm_MU, 1,quantile(x,probs = c(0.025, 0.5, 0.975)))
  # MU1_quantiles <- t(apply(MU1,1,function(x) quantile(x, probs = c(0.025, 0.5, 0.975))))
  # norm_MU1_quantiles <- t(apply(MU1_norm, 1, function(x) quantile(x,probs = c(0.025, 0.5, 0.975))))
  # 
  # matplot(t(x_it1),type='l', col=rep('Dark Grey',n1), main = paste('GOSE1 CREDIBLE BANDS for MU, singal',n_sign))
  # #matlines(MU, type='l', lwd=2) #a little slow (5000 curves are added!)
  # matlines(MU1_quantiles,lty =1, lwd=2, col='magenta') 
  # matlines(norm_MU1_quantiles, lty=2, lwd=1, col='blue')
  # legend('top', legend=c("MU1 form NIW", "MU1 from N"),
  #        col=c('magenta','blue'), lty=1:2, cex=0.6,lwd=c(2,2))
  # #MU1_quantiles
  # #norm_MU1_quantiles
  # # ------------------------------------------------------------- #
  
  
  plot(sigma_2_result1, main = paste('sigma^2 GOSE1, signal',n_sign))
  matplot(sigma_2_result1,type='l', title=title(c('sigma_2_with_Burn_GOSE1 traceplot',paste('signal',n_sign))))
  hist(sigma_2_result1,probability = T)
  
  indici = c(1,3,333,73,65,22,768,53,1121,984,1600,823,1233,1453,653,12,90,25,679,123,124,125)
  for (ind in indici)
    matplot(phi_result1[,ind], xlab='iteration', type='l', title=title(c(paste('phi_GOSE1(',ind,')',sep=''),paste('signal',n_sign)))) #expression(traceplot~phi(125)~GOSE1)))
  #matplot(t(phi_result1), type='l')
  
  
  # likelihood computation ####
  logGOSE1 = calcolo.loglike(niter, burnin,T,MU1,n1,sigma_2_result1, phi_result1, t(x_it1))
  
  #### > > > GOSE 2 < < < < < < < < < < < < < < < < < < ####
  ###########################################################################
  ############################## THE BEGINNING ##############################
  
  n2 <- sum(GOSE==2) # number of patients with GOSE==1
  x_it2 <- fda_signal$data[GOSE==2,]
  
  # betas  ##################################################################
  betas2 <- pca_L$scores[GOSE==2,]
  
  dim(betas2)
  colMeans(pca_L$scores)
  beta_mean2 <- colMeans(betas2)
  
  
  
  ###########################################################################
  ########################### Model for the MEAN ############################
  
  print("Sample GOSE=2 model")
  
  # '- - - NIW - - -' ####
  print("MU2")
  niw2 <- sampling_mean_niw(niter, betas2, beta_mean2, b_t, data_L.fd, L, T, n2)
  k02 <- niw2$k0
  MU2 <- niw2$MU
  
  # # '- - - Normal - - -' ####
  # print("MU2_norm")
  # MU2_norm <- sampling_mean_normal(niter, betas2, beta_mean2, b_t, data_L.fd, L, T, n2)
  
  ###########################################################################
  ######################## Model for the Covariance #########################  
  
  niter <- niter
  burnin <- burnin
  print("sigma^2 and phi 2")
  
  # '- - - NIW - - -' ####
  gib2 <- gibbs(niter, burnin, T, n2, x_it2, MU2)
  sigma_2_result2 <- gib2$sigma
  phi_result2 <- gib2$phi
  
  a02 <- gib2$a0
  b02 <- gib2$b0
  
  # # '- - - Normal - - -' ####
  # gib2_norm <- gibbs(niter, burnin, T, x_it2, MU2_norm)
  # sigma_2_result2_norm <- gib2_norm$sigma
  # phi_result2_norm <- gib2_norm$phi
  
  
  
  # some plots ####
  matplot(signal[GOSE==2,],type='l', main = paste('patients with GOSE2, signal',n_sign))
  
  matplot(MU2, type='l', main = paste('MU GOSE2, signal',n_sign))
  
  
  # # ------------- CREDIBLE INTERVALS FOR THE MEAN -------------- #
  # #tapply(norm_MU, 1,quantile(x,probs = c(0.025, 0.5, 0.975)))
  # MU2_quantiles <- t(apply(MU2,1,function(x) quantile(x, probs = c(0.025, 0.5, 0.975))))
  # norm_MU2_quantiles <- t(apply(MU2_norm, 1, function(x) quantile(x,probs = c(0.025, 0.5, 0.975))))
  # 
  # matplot(t(x_it2),type='l', col=rep('Dark Grey',n2),main = paste('GOSE2 CREDIBLE BANDS for MU, singal',n_sign))
  # #matlines(MU, type='l', lwd=2) #a little slow (5000 curves are added!)
  # matlines(MU2_quantiles,lty =1, lwd=2, col='magenta') 
  # matlines(norm_MU2_quantiles, lty=2, lwd=1, col='blue')
  # legend('top', legend=c("MU2 form NIW", "MU2 from N"),
  #        col=c('magenta','blue'), lty=1:2, cex=0.6,lwd=c(2,2))
  # #MU2_quantiles
  # #norm_MU2_quantiles
  # # ------------------------------------------------------------- #
  
  
  
  plot(sigma_2_result2, main = paste('sigma^2 GOSE2, signal',n_sign))
  matplot(sigma_2_result2,type='l', title=title(c('sigma_2_with_Burn_GOSE2 traceplot', paste('signal',n_sign))))
  hist(sigma_2_result2,probability = T)
  
  indici = c(1,3,333,73,65,22,768,53,1121,984,1600,823,1233,1453,653,12,90,25,679,123,124,125)
  for (ind in indici)
    matplot(phi_result2[,ind], xlab = 'iteration', type='l', title=title(c(paste('phi_GOSE2(',ind,')',sep=''),paste('signal',n_sign))))
  #matplot(t(phi_result2))
  
  
  # likelihood computation ####
  logGOSE2 = calcolo.loglike(niter, burnin,T,MU2,n2,sigma_2_result2, phi_result2, t(x_it2))
  
  ###########################################################################
  ############################## BAYES FACTOR ###############################
  
  print("Computing Bayes Factor")
  
  BF[[n_sign]] = bieffe(logTOT,logGOSE1,logGOSE2,70)
  log2BF[[n_sign]]=2*log(BF[[n_sign]])
  
  alpha <- 100
  BF.power[[n_sign]] = bieffe.power(logTOT,logGOSE1,logGOSE2,70,alpha,n_tot,T)
  log2BF.power[[n_sign]]=2*log(BF.power[[n_sign]])
  
  save.image(paste("wstot",n_sign,".RData",sep=''))
}