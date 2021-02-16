library(fda.usc)
library(fda)
library(fields)


load("data/Codice_analisi_funzionale.RData")
#load("WS_fPCA.Rdata")

########### patients ok note: only 7 over 26 
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


# functional object containing the signal

n_sign = 3
n_tot <- sum(check_signal[,n_sign])
indexes <- which(check_signal[,n_sign]==1)
plot.fdata(f.data[[n_sign]], main = "Left frontal lobe, short latency")


p_ok <- data$GOSE-1
p_ok <- p_ok[indexes] #0 not ok, 1 ok

i_ok <- c()
i_not_ok <- c()
for (i in 1:length(indexes)){
  if (p_ok[i]){
    i_ok <- c(i_ok,i)
  }
  else{
    i_not_ok <- c(i_not_ok,i)
  }
}



#### FUNCTIONAL#####
fda_signal = f.data[[n_sign]] #oggetto funzionale
signal = fda_signal$data #soli dati
signal = t(signal)
time <- 1:1600 #ascissa

m_opt = 377
basis_opt <- create.bspline.basis(rangeval=c(0,1600), nbasis=m_opt)
data_L.fd <- Data2fd(y = signal,argvals = time,basisobj = basis_opt)

plot.fd(data_L.fd, main="B-splines")
pca_L <- pca.fd(data_L.fd,nharm=5,centerfns=TRUE)

# scree plot
plot(pca_L$values,xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_L$values)/sum(pca_L$values),xlab='j',ylab='CPV',ylim=c(0.8,1))
#with 3 components we explain 90% of the variance

# First three FPCs
layout(cbind(1,2,3))
plot(pca_L$harmonics[1,],col=1,ylab='FPC1')
plot(pca_L$harmonics[2,],col=2,ylab='FPC2')
plot(pca_L$harmonics[3,],col=3,ylab='FPC3')

# plot of the principal components as perturbation of the mean
media <- mean.fd(data_L.fd)

plot(media,lwd=2,main='FPC1', ylim=c(-200,200))
lines(media+pca_L$harmonic[1,]*sqrt(pca_L$values[1]), col=2)
lines(media-pca_L$harmonic[1,]*sqrt(pca_L$values[1]), col=3)
# variation in amplitude in the centre

plot(media,lwd=2, ylim=c(-200,200),main='FPC2')
lines(media+pca_L$harmonic[2,]*sqrt(pca_L$values[2]), col=2)
lines(media-pca_L$harmonic[2,]*sqrt(pca_L$values[2]), col=3)
# horizontal translation - difference in timing

plot(media,lwd=2, ylim=c(-200,200),main='FPC3')
lines(media+pca_L$harmonic[3,]*sqrt(pca_L$values[3]), col=2)
lines(media-pca_L$harmonic[3,]*sqrt(pca_L$values[3]), col=3)
# variation in amplitude at the boundaries of the domain

# Command of the library fda that automatically does these plots
par(mfrow=c(1,3))
plot.pca.fd(pca_L, nx=100, pointplot=TRUE, harm=c(1,2,3), expand=0, cycle=FALSE)
# Scores -> check patients ok/not ok
layout(cbind(1,2,3))
plot(pca_L$scores[,1],pca_L$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2, main=n_sign)
points(pca_L$scores[i_ok,1],pca_L$scores[i_ok,2],col="green", lwd=2,pch=8)
#points(pca_L$scores[i_not_ok,1],pca_L$scores[i_not_ok,2],col="red", lwd=2)
plot(pca_L$scores[,1],pca_L$scores[,3],xlab="Scores FPC1",ylab="Scores FPC3",lwd=2, main=n_sign)
points(pca_L$scores[i_ok,1],pca_L$scores[i_ok,3],col="green", lwd=2,pch=8)
#points(pca_L$scores[i_not_ok,1],pca_L$scores[i_not_ok,3],col="red", lwd=4)
plot(pca_L$scores[,2],pca_L$scores[,3],xlab="Scores FPC2",ylab="Scores FPC3",lwd=2, main=n_sign)
points(pca_L$scores[i_ok,2],pca_L$scores[i_ok,3], col="green", lwd=2, pch=8)
#points(pca_L$scores[i_not_ok,2],pca_L$scores[i_not_ok,3], col="red")

layout(cbind(1,2,3))
plot(pca_L$scores[,1],pca_L$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2, col="white", main=n_sign)
text(pca_L$scores[i_ok,1],pca_L$scores[i_ok,2],col="green", lwd=2,pch=8)
text(pca_L$scores[i_not_ok,1],pca_L$scores[i_not_ok,2],col="red", lwd=2)
plot(pca_L$scores[,1],pca_L$scores[,3],xlab="Scores FPC1",ylab="Scores FPC3",lwd=2, col="white", main=n_sign)
text(pca_L$scores[i_ok,1],pca_L$scores[i_ok,3],col="green", lwd=2,pch=8)
text(pca_L$scores[i_not_ok,1],pca_L$scores[i_not_ok,3],col="red", lwd=4)
plot(pca_L$scores[,2],pca_L$scores[,3],xlab="Scores FPC2",ylab="Scores FPC3",lwd=2, col="white", main=n_sign)
text(pca_L$scores[i_ok,2],pca_L$scores[i_ok,3], col="green", lwd=2, pch=8)
text(pca_L$scores[i_not_ok,2],pca_L$scores[i_not_ok,3], col="red")



#signal <- eval.fd(time,data_L.fd)

layout(1)
matplot(signal,type='l')
lines(signal[,24],lwd=4, col=3) #24 patient is ok
lines(signal[,1],lwd=4, col="red") #not ok


layout(1)
matplot(signal,type='l')
lines(rowMeans(signal[,i_ok]),lwd=3, col=3) #24 patient is ok
lines(rowMeans(signal[,i_not_ok]),lwd=3, col="red") #24 patient is ok


graphics.off()
#save.image("WS_fPCA.Rdata")
