
##########################################################################
# Jan 15 2010
# Jeff Goldsmith
#
# this file is an example of the PFR method described in Goldsmith, Feder, 
# Crainiceanu, Caffo, and Reich (2010). several key parameters may be 
# changed to evaluate performance.
##########################################################################


library(nlme)
library(SemiPar)

seed.start <- 1 

## set a few constants
by = 0.1				# density grid on which functions are observed
t=seq(0,10, by=by)		# grid on which functions are observed
n=200 					# sample size
NREPS=100				# number of simulation iterations
N=5						# number of knots in covariance smoothing

## Parameters we want to vary
vareps= 1				# variance of outcome from true mean
varX = 1				# measurement error variance
Kz= 30					# note: kb is taken as min(kz, 35)
kb=min(Kz, 35)


## define true beta function

trueBeta=sin(t*pi/5)
#trueBeta=sqrt(t)
#trueBeta=-dnorm(t, mean=2, sd=.3)+3*dnorm(t, mean=5, sd=.4)+dnorm(t, mean=7.5, sd=.5)

## generate the true functions 

true.funcs <- matrix(0, nrow=n, ncol=length(t))
for(i in 1:n){
	true.funcs[i,]=true.funcs[i,]+runif(1, 0, 5)
	true.funcs[i,]=true.funcs[i,]+rnorm(1, 1, 0.2)*t
 	for(j in 1:10){
		e=rnorm(2, 0, 1/j^(2))
		true.funcs[i,]=true.funcs[i,]+e[1]*sin((2*pi/10)*t*j)
		true.funcs[i,]=true.funcs[i,]+e[2]*cos((2*pi/10)*t*j) 
	}
}

## generate outcomes	
outcomes <- sapply(1:n, function(u) sum(true.funcs[u,]*trueBeta)*by)+rnorm(n, 0, sqrt(vareps))


## add measurement error to obtain the observed functions
funcs=true.funcs + matrix(rnorm(length(t)*n, 0, sqrt(1)), nrow=n, ncol=length(t))


## plot the true and observed functions for several subjects
par(mfrow=c(2,4))
for(i in 1:8){
	plot(t, true.funcs[i,], type='l', ylab=expression(X[i](t)))
	points(t, funcs[i,], type='l', ylab=expression(), lty=2)
}


## obtain the eigen decomposition of the observed functions

varFuncs <- var(funcs)

## Smooth the covariance matrix

## G2 is a M*M matrix for the raw covariance function
G2 <- varFuncs
M <- length(t)
diag(G2)= rep(NA, M)
g2 <- as.vector(G2)

## define a N*N knots for bivariate smoothing
x1 <- rep(seq(0,10,length=N), N)
x2 <- rep(seq(0,10,length=N), each=N)
myknots <- data.frame(x1=x1, x2=x2)

## bivariate smoothing using the spm function
t1 <- rep(t, M)
t2 <- rep(t, each=M)
fit <- spm(g2 ~ f(t1, t2, knots=myknots), omit.missing=T)
newdata <- data.frame(t1=t1,t2=t2)
pred <- predict(fit,newdata)

##  make the fitted covariance function symmetric
temp.g2 <- matrix(pred, M)
fit.g2 <- (temp.g2 + t(temp.g2) )/2
		
eigenDecomp <- eigen(fit.g2)
eigenDecomp$vectors <- eigenDecomp$vectors/sqrt(by)


## plot the first few eigenfunctions
par(mfrow=c(1,4))
for(i in 1:4){
	plot(t, eigenDecomp$vectors[,i], type='l', ylab=expression(X[i](t)))
}


##################################################
## PFR - this functions implements the method from the paper.
##################################################

fglmJJ <- function(outcomes, funcs, eigenDecomp=NULL) {
	require(nlme)
	
	## get the eigen decomposition of the smoothed variance matrix
	
	if(is.null(eigenDecomp)) {
		varFuncs <- var(funcs)

		eigenDecomp <- eigen(varFuncs)
		eigenDecomp$vectors <- eigenDecomp$vectors/sqrt(by)
	}

	# set the basis to be used for beta(t)
	kz <- Kz ## since the notation hasn't been consistent
	num=kb-2
	qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
	knots <- quantile(t, qtiles)
##	Basis = bs(t, kb)
	Basis = cbind(1, t, t^2, sapply(knots, function(k) ((t - k > 0) * (t - k)) ^ 2))
##	Basis = cbind(1, t, sapply(knots, function(k) ((t - k > 0) * (t - k))))

	## set up the matrices to be fit in the function lme()
	C <- by* funcs %*% eigenDecomp$vectors[ ,1:kz ]
	J <- by * t(eigenDecomp$vectors[,1:kz]) %*% Basis
	CJ <- C %*% J

	w <- cbind(1, CJ)
	X=as.matrix(w[,1:4])
	Z=as.matrix(w[,5:dim(w)[2]])
	
	group<-rep(1, n)
	##fit=lme(outcomes~-1+X, random=list(group=pdIdent(~-1+Z)))
	## in the cases where we would be getting NAs, use the "optim" function
	fit=try(lme(outcomes~-1+X, random=list(group=pdIdent(~-1+Z))))
	if(class(fit)=="try-error") {
		fit=lme(outcomes~-1+X, random=list(group=pdIdent(~-1+Z)), control=	lmeControl(opt="optim"))
		print("Used the 'optim' method")
	}
	
	## get the coefficient and betaHat estimates
	coefs <- c(fit$coef$fixed,unlist(fit$coef$random))
	fitted <- as.matrix(w[,1:length(coefs)]) %*% coefs
	betaHat <- Basis %*% coefs[-1]
	
	D=diag(c(rep(0,ncol(X)), rep(1, ncol(Z))))
	lambda <- (fit$sigma)^2/as.numeric(VarCorr(fit)[1,1])
	sigmaEpsHat=(fit$sigma)^2
	varBeta=(sigmaEpsHat * solve(t(w)%*%w + lambda*D) %*% t(w) %*% w %*% solve(t(w)%*%w + lambda*D))[-1,-1]
	varBetaHat=Basis%*%varBeta%*%t(Basis)
	
	ret <- list(fit, coefs, fitted, betaHat, w, X, Z, Basis, varBetaHat)
	names(ret) <- c("fit", "coefs", "fitted", "betaHat", "w", "X", "Z", "Basis", "varBetaHat")
	ret
}


## fit the model using the PFR method (the fglmJJ function
## uses a mixed model to estimate beta(t))
fitJJ = try( fglmJJ(outcomes, funcs, eigenDecomp=eigenDecomp))


## plot the true and estimated coefficient function
plot(t, trueBeta, type='l', lwd=2)
points(t, fitJJ$betaHat, type='l', lwd=2, col="red")


## calculate the MSE
MSE=mean((trueBeta - fitJJ$betaHat)^2)
MSE









