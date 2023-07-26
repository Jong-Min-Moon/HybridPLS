##########################################################################
# Nov 13 2009
# Jeff Goldsmith
#
# this file implements a simulation study to evaluate the multivariate
# functional regression technique in Goldsmith, Feder, Crainiceanu, Caffo,  
# and Reich (2010).
##########################################################################

seed.start = 1

library(SemiPar)
library(nlme)
library(MASS)
library(boot)
library(plotrix)


by = 0.1			# density of grid on which functions are observed
t=seq(0,10, by=by)	# grid on which functions are observed
n=200 				# sample size

NREPS=1000			# number of iterations in the simulation

vareps= .5 			# variance of outcome from true mean (.5, 1)
VarX = 0			# measurement error variance         (0, 1)
kz= 20					# note: kb is taken as min(kz, 35)
kb=min(kz, 35)
smoothCov=FALSE


## define true beta functions
trueBeta1=4*sin(t*pi/5)
trueBeta2=(t/5)^2


# keep track of the MSEs of the estimated functions
MSE_MULT=matrix(0, ncol=2, nrow=NREPS)
MSE_PVE=matrix(0, ncol=2, nrow=NREPS)

# keep track of the seed for reproducibility
seeds = matrix(seed.start + 0:(NREPS-1), nrow=NREPS, ncol=1)


for(i in 1:NREPS){
	
	set.seed(seeds[i])
	
	# generate subject-specific functions
	funcs1.true <- matrix(0, nrow=n, ncol=length(t))
	for(i2 in 1:n){
		funcs1.true[i2,]=funcs1.true[i2,]+runif(1, 0, 5)
		funcs1.true[i2,]=funcs1.true[i2,]+rnorm(1, 1, 0.2)*t
 		for(j2 in 1:10){
			e=rnorm(2, 0, 1/j2^(2))
			funcs1.true[i2,]=funcs1.true[i2,]+e[1]*cos((2*pi/10)*t*j2) 
			funcs1.true[i2,]=funcs1.true[i2,]+e[2]*sin((2*pi/10)*t*j2) 
		}
	}	
	
	funcs2.true <- matrix(0, nrow=n, ncol=length(t))
	for(i2 in 1:n){
		funcs2.true[i2,]=funcs2.true[i2,]+runif(1, 0, 5)
		funcs2.true[i2,]=funcs2.true[i2,]+rnorm(1, 1, 0.2)*t
 		for(j2 in 1:10){
			e=rnorm(2, 0, 1/j2^(2))
			funcs2.true[i2,]=funcs2.true[i2,]+e[1]*cos((2*pi/10)*t*j2) 
			funcs2.true[i2,]=funcs2.true[i2,]+e[2]*sin((2*pi/10)*t*j2) 
		}
	}	
	
	# add measurement error
	funcs1=funcs1.true+rnorm(length(funcs1.true), mean=0, sd=sqrt(VarX))
	funcs2=funcs2.true+rnorm(length(funcs2.true), mean=0, sd=sqrt(VarX))

	# construct and smooth covariance matrices
	covFuncs1=cov(funcs1)

	# smooth the first covariance matrix
	m=dim(covFuncs1)[2]
	gw.temp <- covFuncs1
	diag(gw.temp) <- rep(NA, m)
	gw <- as.vector(gw.temp)
	x1 <- rep(seq(0,1,length=m), m)
	x2 <- rep(seq(0,1,length=m), each=m)   
	# myknots <- data.frame(x1=rep(seq(0,1,length=7)[2:6],5), x2=rep(seq(0,1,length=7)[2:6],each=5) )
	# NOTE it looks like when sigma_epsilon is small, sometimes the above knots give error code and
	# the knots below work without any problem (other times is the other way around)  
	myknots <- data.frame(x1=rep(seq(0,1,length=5),5), x2=rep(seq(0,1,length=5),each=5)   ) 
	fit.w <- spm(gw ~ f(x1, x2, knots = myknots),omit.missing=TRUE)         
	newdata <- data.frame(x1=x1,x2=x2) 
	pred.w <- predict(fit.w,newdata)
	s.gw <-matrix(pred.w,m) 

	smoothCov1 <- (s.gw + t(s.gw) )/2 
	eigenDecomp1 = eigen(smoothCov1)
	eigenDecomp1$vectors = eigenDecomp1$vectors/sqrt(by)


	covFuncs2=cov(funcs2)

	# smooth the first covariance matrix
	m=dim(covFuncs2)[2]
	gw.temp <- covFuncs2
	diag(gw.temp) <- rep(NA, m)
	gw <- as.vector(gw.temp)
	x1 <- rep(seq(0,1,length=m), m)
	x2 <- rep(seq(0,1,length=m), each=m)   
	# myknots <- data.frame(x1=rep(seq(0,1,length=7)[2:6],5), x2=rep(seq(0,1,length=7)[2:6],each=5) )
	# NOTE it looks like when sigma_epsilon is small, sometimes the above knots give error code and
	# the knots below work without any problem (other times is the other way around)  
	myknots <- data.frame(x1=rep(seq(0,1,length=5),5), x2=rep(seq(0,1,length=5),each=5)   ) 
	fit.w <- spm(gw ~ f(x1, x2, knots = myknots),omit.missing=TRUE)         
	newdata <- data.frame(x1=x1,x2=x2) 
	pred.w <- predict(fit.w,newdata)
	s.gw <-matrix(pred.w,m) 

	smoothCov2 <- (s.gw + t(s.gw) )/2 
	eigenDecomp2 = eigen(smoothCov2)
	eigenDecomp2$vectors = eigenDecomp2$vectors/sqrt(by)


	# generate outcomes
	errors=rnorm(n, 0, sqrt(vareps))
	outcomes <- sapply(1:n, function(u) {sum(funcs1.true[u,]*trueBeta1)*by+sum(funcs2.true[u,]*trueBeta2)*by})+errors


	# de-mean the functions
	meanFunc1=apply(funcs1, 2, mean)
	funcs1=sapply(1:length(t), function(u) funcs1[,u]-mean(funcs1[,u]))

	meanFunc2=apply(funcs2, 2, mean)
	funcs2=sapply(1:length(t), function(u) funcs2[,u]-mean(funcs2[,u]))


	# set up the matrices to be fit in the function lme()

	# first, calculate the PC loadings
	C1 <- by * funcs1 %*% eigenDecomp1$vectors[ ,1:kz ] 
	C2 <- by * funcs2 %*% eigenDecomp2$vectors[ ,1:kz ] 

	# set the basis to be used for beta(t)
	num=kb-2
	qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
	knots <- quantile(t, qtiles)
	Basis = cbind(rep(1, length(t)), t, sapply(knots, function(k) ((t - k > 0) * (t - k)) ))

	# evaluate the J matrix
	J1 <- by * t(eigenDecomp1$vectors[,1:kz]) %*% Basis
	J2 <- by * t(eigenDecomp2$vectors[,1:kz]) %*% Basis

	# and the CJ matrix
	CJ1 <- C1[,1:kz] %*% J1
	CJ2 <- C2[,1:kz] %*% J2

	# the following code partitions the design matrices so that
	# the lme function performs the correct estimation

	w <- cbind(1, CJ1, CJ2)
	X=as.matrix(w[,c(1:3, (kb+2):(kb+3))])
	Z1=as.matrix(w[,c(4:(kb+1))])
	Z2=as.matrix(w[,c((kb+4):dim(w)[2])])

	Z=cbind(Z1, Z2)

	re.block.inds=list(1:(kb-2), (kb+1-2):(kb-2+kb-2))
	Z.block=list(length=2)
	for (i2 in 1:length(re.block.inds)) 
		Z.block[[i2]] <- as.formula(paste("~Z[,c(",paste( re.block.inds[[i2]],collapse=","),")]-1"))

	colnames(X)=c("int","1.1", "t.1", "1.2", "t.2")
	#colnames(Z)=c(paste("1.",3:kb, sep=""), paste("2.",3:kb,sep="" ))

	group=rep(1, n)
	grouped=data.frame(X,outcomes,Z)
	model.data=groupedData(outcomes~X|group, data=grouped)

	
	# do the regression using a mixed model

	fit.smooth=try( lme(outcomes~-1+X, random=list(group=pdBlocked(Z.block, pdClass=rep("pdIdent",length(Z.block))))) )
	if(class(fit.smooth) == "try-error") {
		fit.smooth=lme(outcomes~-1+X, random=list(group=pdBlocked(Z.block, pdClass=rep("pdIdent",length(Z.block)))) 
			, control=lmeControl(opt="optim"))
	} else {
	

		# get the coefficient and betaHat estimates
		coefs <- c(fit.smooth$coef$fixed,unlist(fit.smooth$coef$random))
		fitted <- as.matrix(w[,1:length(coefs)]) %*% coefs
		betaHat1.MULT <- Basis %*% c(coefs[c(2:3, 6:(kb+1+2))])
		betaHat2.MULT <- Basis %*% c(coefs[c(4:5, (kb+2+2):length(coefs))])

		#calculate MSEs
		MSE_MULT[i,1]=mean((trueBeta1 - betaHat1.MULT)^2)
		MSE_MULT[i,2]=mean((trueBeta2 - betaHat2.MULT)^2)
	}

	##############################################################
	# additionally, we use an FPCR method choosing the number of 
	# components via PVE. this is an intuitive approach to the
	# multivariate FR problem.
	##############################################################

	# determine the number of loadings to use via PVE
	ro1=cumsum(eigenDecomp1$values)/sum(eigenDecomp1$values)
	if(VarX==0){
		K1 <- min(which( ro1 >= 0.99 ))
	} else {
		K1 <- min(which( ro1 >= 0.99 ))
	}

	ro2=cumsum(eigenDecomp2$values)/sum(eigenDecomp2$values)
	if(VarX==0){
		K2 <- min(which( ro2 >= 0.99 ))
	} else {
		K2 <- min(which( ro2 >= 0.99 ))
	}

	# fit the model
	fit.PVE=lm(outcomes~-1+C1[,1:K1]+C2[,1:K2])
	
	# estimate the coefficient functions
	coefs <- coef(fit.PVE)
	betaHat1.PVE <- eigenDecomp1$vectors[ ,1:K1 ] %*% coefs[1:K1]
	betaHat2.PVE <- eigenDecomp2$vectors[ ,1:K2 ] %*% coefs[(K1+1):(K1+K2)]

	#calculate MSEs
	MSE_PVE[i,1]=mean((trueBeta1 - betaHat1.PVE)^2)
	MSE_PVE[i,2]=mean((trueBeta2 - betaHat2.PVE)^2)
	
	print(i)

}


save(MSE_MULT, MSE_PVE, VarX, vareps, seeds, file="Mult_Sim_Results1.rda")


