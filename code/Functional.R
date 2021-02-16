needed_packages  <- c("fda.usc","fda","fields")
new_packages  <- needed_packages[!(needed_packages %in%installed.packages ()[, "Package"])]
if (length(new_packages))  install.packages(new_packages)
lapply(needed_packages , require , character.only = TRUE)

load("data/Codice_analisi_funzionale.RData")

#Results are related to signal 3, but they are very similar to each other
#just change n_sign to try the others
n_sign = 3
plot.fdata(f.data[[n_sign]], main = n_sign)

#### FUNCTIONAL#####
#I study the signal identified with the progressive number n_sign
fda_signal = f.data[[n_sign]]#functional object
signal = fda_signal$data #pure data
signal = t(signal)
dim(signal) #1600x26
#1600 are the time instants
#26 are the 26 different patients
#
#Plot funzionale (same as the plot before) <- fda_signal is a functional object
plot.fdata(fda_signal, main = "n_sign")

#Another kind of plot <- signal is just a table
#x11()
matplot((signal),type='l',main='Frontale sinistro short latency',xlab='Time',ylab='X(t)')

time = 1:1600 #abscissa


#### FOURIER ####
#Fourier is for periodic signals.. not here, it would screw up the behaviour at the boundaries

# Choice 1: we set a high dimensional basis (almost interpolating)
# Pros: no loss of information
# Cons: possible overfitting 
# basis.1 <- create.fourier.basis(rangeval=c(0, 1600),nbasis=1600)
# data_signal.fd.1 <- Data2fd(y = signal,argvals = time,basisobj = basis.1)
# 
# # Choice 2: reduced dimensionality (we set a low dimensional basis)
# # Pros: the data are much smoother and the measurement error is filtered
# # Cons: I could have lost important information
# basis.2 <- create.fourier.basis(rangeval=c(0,1600),nbasis=20)
# data_signal.fd.2 <- Data2fd(y = signal,argvals = time,basisobj = basis.2)
# 
# # Choice 3: compromise between 1 and 2
# basis.3 <- create.fourier.basis(rangeval=c(0,1600),nbasis=500)
# data_signal.fd.3 <- Data2fd(y = signal,argvals = time,basisobj = basis.3)

# x11()
# par(mfrow=c(4,1))
# plot.fd(data_signal.fd.1, main = "1600 basi")
# plot.fd(data_signal.fd.2, main = "20 basi")
# plot.fd(data_signal.fd.3, main = "500 basi")
# plot.fdata(fda_signal, main = "Frontale sinistro short latency")



#### BSPLINE ####
#Better to use a b-spline basis
#n. basis
basis <- create.bspline.basis(rangeval=c(0,1600), nbasis=50) # just a try, to see how it goes
data_signal.fd <- Data2fd(y = signal,argvals = time,basisobj = basis)

x11()
par(mfrow=c(2,1))
plot.fd(data_signal.fd, main="B-splines")
plot.fdata(fda_signal, main = "dati veri")

# Estimate of the mean and of the covariance kernel
layout(cbind(1,2))
plot.fd(data_signal.fd,xaxs='i')
lines(mean(data_signal.fd),lwd=2)
eval <- eval.fd(time,data_signal.fd)
image.plot(time, time, (cov(t(eval))[1:1600,])) #ok? come la dovrei leggere? bo'


# CUBIC SPLINES
m <- 4        # spline order: cubiche 
degree <- m-1    # spline degree 
nbasis <- 50
NT <- length(time)
Xobs0 = signal[,1] # just pick up ine among the 26 signals, to test it

# Create the basi
basis <- create.bspline.basis(rangeval=c(0,1600), nbasis=nbasis, norder=m)
# If breaks are not provided, equally spaced knots are created; rangeval is the abscissa range
names(basis)
plot(basis) 


# Evaluate the basis on the grid of abscissa
basismat<- eval.basis(time, basis) #matrix used to evaluate our functions
dim(basismat) #1600x50 here
head(basismat)

# Fit via LS, least square
lsfit(basismat, Xobs0, intercept=FALSE)$coef #coefficient for the 50 basis

Xsp0 <- basismat %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef #interpolation
head(Xsp0)
dim(Xsp0)


par(mfrow=c(1,1))
plot(time,Xobs0,xlab="t",ylab="observed data")
points(time,Xsp0,type="l",col="blue",lwd=2)



#Differences plot of the datum, to see how much noise it has
rappincX1 <- (Xobs0[3:NT]-Xobs0[1:(NT-2)])/(time[3:NT]-time[1:(NT-2)])
rappincX2 <- ((Xobs0[3:NT]-Xobs0[2:(NT-1)])/(time[3:NT]-time[2:(NT-1)])-(Xobs0[2:(NT-1)]-Xobs0[1:(NT-2)])/(time[2:(NT-1)]-time[1:(NT-2)]))*2/(time[3:(NT)]-time[1:(NT-2)])


#1st der. (argument Lfdobj=1)
basismat1<- eval.basis(time, basis, Lfdobj=1)
Xsp1 <- basismat1 %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef

#2nd der. (argument Lfdobj=2)
basismat2<- eval.basis(time, basis, Lfdobj=2)
Xsp2 <- basismat2 %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef

#less noise
par(mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(time,Xobs0,xlab="t",ylab="observed data")
points(time,Xsp0 ,type="l",col="blue",lwd=2)
plot(time[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(time,Xsp1 ,type="l",col="blue",lwd=2)
plot(time[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(time,Xsp2 ,type="l",col="blue",lwd=2)
plot(basis)

x11() #just a zoom of the second differences
plot(time[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(time,Xsp2 ,type="l",col="blue",lwd=2)

#### CROSS VALIDATION QUADRATIC BSPLINES ####
# Optimal n. bases?
# Generalized cross-validation
nbasis <- 400:450
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(0,1600), nbasis[i], 3)
  gcv[i] <- mean(smooth.basis(time, signal, basis)$gcv)
  #gcv[i] <- smooth.basis(time, signal, basis)$gcv[1]
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]
#424 -> 0.1195

#### CROSS VALIDATION with CUBIC SPLINE####
# optimal n. bases?
# Generalized cross-validation
nbasis <- 4:400
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(0,1600), nbasis[i], 4)
  gcv[i] <- mean(smooth.basis(time, signal, basis)$gcv)
  #gcv[i] <- smooth.basis(time, signal, basis)$gcv[1]
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
plot(nbasis,log(gcv))
nbasis[which.min(gcv)] 
#gcv  377 with 0.1138

nbasis <- 100
basis <- create.bspline.basis(c(0,1600), nbasis, 4)
mean(smooth.basis(time, signal, basis)$gcv)


#### try: 93 BASES #### attention, create.bspline.basis has default norder=4 (CUBIC spline)
m_opt = 93
basis_opt <- create.bspline.basis(rangeval=c(0,1600), nbasis=m_opt)
data_signal.fd_opt <- Data2fd(y = signal,argvals = time,basisobj = basis_opt)

par(mfrow=c(2,1))
plot.fd(data_signal.fd_opt, main="B-splines with optimal number")
plot.fdata(fda_signal, main = "Dati veri")

# Estimate of the mean and of the covariance kernel
layout(cbind(1,2))
plot.fd(data_signal.fd_opt,xaxs='i')
lines(mean(data_signal.fd_opt),lwd=2)

eval <- eval.fd(time,data_signal.fd_opt)
image.plot(time, time, (cov(t(eval))[1:1600,])) 

#### try: 377 BASES #### best for CUBIC splines
m_opt = 377
basis_opt <- create.bspline.basis(rangeval=c(0,1600), nbasis=m_opt)
data_signal.fd_opt <- Data2fd(y = signal,argvals = time,basisobj = basis_opt)

par(mfrow=c(2,1))
plot.fd(data_signal.fd_opt, main="B-splines with optimal number")
plot.fdata(fda_signal, main = "Dati veri")

# Estimate of the mean and of the covariance kernel
layout(cbind(1,2))
plot.fd(data_signal.fd_opt,xaxs='i')
lines(mean(data_signal.fd_opt),lwd=2)

eval <- eval.fd(time,data_signal.fd_opt)
image.plot(time, time, (cov(t(eval))[1:1600,])) 



