load("Codice_analisi_funzionale.RData")
a = rep(1,8) #all signals present
# Ordine: ausxSL  audxSL  frsxSL  frdxSL ausxML  audxML  frsxML  frdxML
a6 = c(1,1,1,1,0,1,0,1) #missing au&fr sxML
a10 = c(1,1,1,0,1,1,1,1) #missing frsxML
a11 = c(1,1,1,1,1,1,1,0) #missing frdxML
a14 = c(1,1,1,1,1,0,1,0) #missing au&fr dxML 
a15 = c(1,1,1,1,0,1,0,1) #missing au&fr sxML
a22 = c(1,1,1,1,0,0,0,0) #missing all ML
a23 = c(1,1,1,1,0,0,0,0) #missing all ML
a25 = c(1,0,1,0,1,1,1,1) #problems with au&fr dxSL
check_signal = rbind(a,a,a,a,a,a6,a,a,a,a10,a11,a,a,a14,a15,a,a,a,a,a,a,a22,a23,a,a25,a)

for (n_sign in 1:8){
  print(n_sign)
  
  n_tot <- sum(check_signal[,n_sign])
  indexes <- which(check_signal[,n_sign]==1)
  # functional object containing the signal
  fda_signal <- f.data[[n_sign]] 
  n <- dim(fda_signal$data)[1] # number of patients, usually between 22 and 26
  signal <- fda_signal$data 
  signal <- t(signal)
  
  GOSE <- data[indexes,"GOSE"]
  matplot(signal, type='l', main = paste('all patients signal',n_sign), ylim=c(min(signal),max(signal)), col=indexes)
  matplot(signal[,GOSE==1], type='l', main = paste('GOSE1 patients signal',n_sign),ylim=c(min(signal),max(signal)), col=indexes[GOSE==1])
  matplot(signal[,GOSE==2], type='l', main = paste('GOSE2 patients signal',n_sign), ylim=c(min(signal),max(signal)), col = indexes[GOSE==2])
  }