
#############################################################
# Dec 4 2009
# Jeff Goldsmith with code from Jen Feder and Chongzhi Di
#
# the following code implements a penalized spline approach to
# a case of multilevel functional regression. the model is that
# examined in Goldsmith, Feder, Crainiceanu, Caffo, and Reich 
# (2010).
#############################################################

library(nlme)
library(SemiPar)
library(splines)
library(R2WinBUGS)

setwd("C:/FDA")			# location of WinBUGS file "scores-projection.txt"

seed.start = 1371
NREPS=100				# number of iterations in the simulation

n=M=200  				# sample size
J = 3					# number of visits
VarX=0					# measurement error variance
vareps=.5				# variance on the outcome
tlength=10

## set the standard parameters
by = 0.1				# density of grid on which functions are observed
t=seq(0,10, by=by)		# grid on which functions are observed
N = length(t)
kz1=10					# number of principal components used to estimate
kz2=3					# X_{i} and U_{ij}

## select true beta functions
trueBeta1=2*sin(t*pi/5)
trueBeta2=(t/2.5)^2


# keep track of the MSEs of the estimated functions
MSE_MULTLEV=matrix(0, ncol=4, nrow=NREPS)
PVE_MULTLEV=matrix(0, ncol=4, nrow=NREPS)

# keep track of the seed for reproducibility
seeds = matrix(seed.start + 0:(NREPS-1), nrow=NREPS, ncol=1)

for(i in 1:NREPS){
	
	set.seed(seeds[i])
	
	## generate true subject-specific mean functions
	funcs <- matrix(0, nrow=n, ncol=length(t))
	for(i2 in 1:n){
		funcs[i2,]=funcs[i2,]+runif(1, 0, 5)
		funcs[i2,]=funcs[i2,]+rnorm(1, 1, 0.2)*t
 		for(j2 in 1:10){
			e=rnorm(2, 0, 1/j2^(2))
			funcs[i2,]=funcs[i2,]+e[1]*sin((2*pi/10)*t*j2)
			funcs[i2,]=funcs[i2,]+e[2]*cos((2*pi/10)*t*j2) 
		}
	}


	###     generate the true eigenvalues for the level 2 functions
	lambda2 <-  2*0.5^(0:(kz2-1))

	###     The true level 2 eigenfunctions are four muturally orthogonal 
	###     polynomial functions
	tt <- t/tlength
	f2 <- matrix(0, nrow=kz2, ncol=N)
	f2[1,] <- rep(1, N)*sqrt(1/tlength)
	f2[2,] <- sqrt(3/tlength) * (2*tt - 1)
	f2[3,] <- sqrt(5/tlength) * (6*tt^2 - 6 * tt + 1)


	###     Generate M random coefficients si's from normal distributions 
	###     to form the subject-visit specific deviations
	si2 <- matrix(0, nrow=M*J, ncol=kz2)
	for(k in 1:kz2) {
    	si2[,k] <- rnorm( M * J ,sd=sqrt(lambda2[k]) )
	}

	# Calculate the observedfunctions as the sum of subject's function
	# the level 2 eigenfunctions (noise is also added)
	y <- matrix(0,nrow=M*J,ncol=N)
	for(m in 1:M) {
    	temp <- funcs[m,]
	    for(j in 1:J) {
    	    y[ (m-1)*J + j ,] <- temp + apply( ( si2[(m-1)*J + j ,] %*% t(rep(1,N)) ) * f2, 2,sum) + rnorm(N, sd=VarX)
	    }
	}



###########################################################################################
####    The following code (from C. Di) performs Multilevel FPCA as described in Di et al 
####    
###########################################################################################

	mu <- apply(y, 2,mean)
	eta <- matrix(0, J, N)
	for(j in 1:J) {
    	eta[j,] <- apply(y[ (0:(M-1)*J) + j,  ], 2, mean) - mu
	}

	###    Calculate residuals by subtracting visit-specific mean from original functions 
	resd <- matrix(0, nrow=M*J, ncol=N) 
	for(j in 1:J) 
    	resd[ (0:(M-1)*J) + j,  ] <- t ( t( y[ (0:(M-1)*J) + j,  ] ) - (mu + eta[j,]) ) 
	resd.rev <- matrix(0, nrow=M*J, ncol=N) 
	for(m in 1:M) {
    	resd.rev[ (m-1)*J + 1:2 ,] <- resd[(m-1)*J + 2:1, ]
	}

	###     Calculate the pairwise difference between different visits within the same subject.
	###     There are J(J-1)/2 possible pairs within the same subject.
   	resd.pair.diff <- matrix(0, M*J*(J-1)/2, ncol=N)
	index <- matrix(0, J*(J-1)/2, 2)
    cur.index <- 1
    for(k in 1:(J-1)) 
        for(l in (k+1):J) {
            index[cur.index,] <- c(k,l)
            cur.index <- cur.index + 1
        }
    for(m in 1:M) {
        for(k in 1: (J*(J-1)/2)) {
            resd.pair.diff[ (m-1)*J*(J-1)/2 + k, ] <- resd[(m-1)*J + index[k,1], ] - resd[(m-1)*J + index[k,2], ]
        }
    }


	###     Estimate the three covariance functions: overall covariance G, 
	###     between covariance Gb and within covariance Gw
    G <- matrix(0, N, N)
    Gb <- matrix(0, N, N)
    Gw <- matrix(0, N, N)
    for(i2 in 1:N) 
        for(j in i2:N) {
            G[i2,j] <- cov(resd[,i2],resd[,j])
            G[j,i2] <- G[i2,j]
            Gw[i2,j] <- cov(resd.pair.diff[,i2],resd.pair.diff[,j]) /2
            Gw[j,i2] <- Gw[i2,j]
        }
        
    Gb = G - Gw    

	##### smooth the kernels in the presence of noise
    gw.temp <- Gw
    diag(gw.temp) <- rep(NA, N)
    gw <- as.vector(gw.temp)
    gb <- as.vector(Gb)
    x1 <- rep(seq(0,1,length=N), N)
    x2 <- rep(seq(0,1,length=N), each=N)

    myknots <- data.frame(x1=rep(seq(0,1,length=10),each=10), x2=rep(seq(0,1,length=10),10))
    fit<- spm(gw ~ f(x1, x2,knots=myknots),omit.missing=T)
    newdata <- data.frame(x1=x1,x2=x2)
    pred <- predict(fit,newdata)
    fit1<- spm(gb ~ f(x1, x2,knots=myknots))
    pred1 <- predict(fit1,newdata)


    sw <- mean( diag(Gw) - diag(matrix(pred,N)) )
    s.gw <-matrix(pred, N) 
    Gw <- (s.gw + t(s.gw) )/2 
    s.gb <- matrix(pred1, N)
    Gb <- (s.gb + t(s.gb))/2

	####    estimate eigen values and eigen funcitons at two levels
	fpca1.value <- eigen(Gb)$values[1:(kz1)]*by
	fpca2.value <- eigen(Gw)$values[1:(kz2)]*by

	fpca1.vectors <- eigen(Gb)$vectors[, 1:kz1]/sqrt(by)
	fpca2.vectors <- eigen(Gw)$vectors[, 1:kz2]/sqrt(by)


	for(i2 in 1:kz2) {
    	v1 <- f2[i2,]
	    v2 <- fpca2.vectors[,i2]
    	tempsign <- sum((v1-v2)^2) - sum((v1+v2)^2)
	    fpca2.vectors[,i2] <- ifelse(tempsign>0, -1,1) * v2
	}





	##################################################################################
	####    Estimate the principal component scores by Bayesian MCMC
	##################################################################################


	# First, project each function (per subject per visit) onto the functional space spanned 
	# by the level 1 eigenfunctions and the space spanned by the level 2 eigenfunctions
	mc <- matrix(0, kz1, kz2)
	a <- matrix(0, M*J, kz1)
	b <- matrix(0, M*J, kz2) 

	####    calculate the inner product between level 1 and level 2 eigenfunctions
    for(i2 in 1:kz1)
        for(j in 1:kz2) 
            mc[i2,j] <- sum(fpca1.vectors[,i2]* fpca2.vectors[,j]) * by
 
	####    calculate the projections of each function onto the level 1 and level 2 eigen functional space
    for(i2 in 1:(M*J))   {
        for(j in 1:kz1)  {
            a[ i2 ,j] <- sum( resd[i2,] * fpca1.vectors[,j] ) * by
        }
        for(j in 1:kz2) {
            b[ i2 ,j] <- sum( resd[i2,] * fpca2.vectors[,j] ) * by
        }
    }

	####    Set eigenvalues and eigenfunctions as estimated previousely
	phi1 <- t(fpca1.vectors)
	phi2 <- t(fpca2.vectors)
	tau1 <- 1/fpca1.value
	tau2 <- 1/fpca2.value


	####    Set up the data set and initial values for MCMC
	data <- list("M","J","kz1","kz2","a","b","mc", "tau1", "tau2")
	inits <- function(){
    	list(tau=1 )
	}  

	####    Call "bugs" functions to run MCMC via WinBUGS
	####    We would use 2500 iterations with 500 burn-ins
	Bayes.fit <- try(bugs(data, inits=inits, model.file = "scores-projection.txt",
    parameters = c("sigma", "s1","s2"),
    n.chains = 1, n.iter = 500, n.burnin=50, n.thin=1,
    bugs.directory = "c:/Program Files/WinBUGS14/",debug=FALSE) )

	if(class(Bayes.fit) == "try-error") {
		MSE_MULTLEV[i,]=c(NA, NA)
		PVE_MULTLEV[i,]=c(NA, NA)
	} else {
		####    Export and summarize results from MCMC
		####        "summ.m3" would contain the posterior mean, SD, "2.5%", "50%", and "97.5%" pencentiles
		####        for every parameter (noise variance and PC scores)
		summ.m3 <- matrix(0, M*(kz1+kz2*J)+1, 5)
		rownames(summ.m3) <- colnames(Bayes.fit$sims.matrix)[1:(M*(kz1+kz2*J)+1)]
		colnames(summ.m3) <- c("mean", "sd", "2.5%", "50%", "97.5%")
		summ.m3[,1] <- apply(Bayes.fit$sims.matrix[,1:(M*(kz1+kz2*J)+1)], 2, mean)
		summ.m3[,2] <- apply(Bayes.fit$sims.matrix[,1:(M*(kz1+kz2*J)+1)], 2, sd)
		summ.m3[,3] <- apply(Bayes.fit$sims.matrix[,1:(M*(kz1+kz2*J)+1)], 2, quantile, prob=0.025)
		summ.m3[,4] <- apply(Bayes.fit$sims.matrix[,1:(M*(kz1+kz2*J)+1)], 2, quantile, prob=0.50)
		summ.m3[,5] <- apply(Bayes.fit$sims.matrix[,1:(M*(kz1+kz2*J)+1)], 2, quantile, prob=0.975)
		round(summ.m3, digits=3)[1:10,]

		####    Transform the posterior mean and SD into matrix form
		####    "si1.est3" contains the posterior mean for level 1 PC scores for every subject. 
		####        Each row corresponds to a subject, and each column corresponds to a component. 
		####    "si2.est3" contains the posterior mean for level 2 PC scores for every subject and visit. 
		####    "si1.sd.est3" contains the posterior SD for level 1 PC scores for every subject. 
		####    "si2.sd.est3" contains the posterior SD for level 2 PC scores for every subject and visit. 
		est.si1 <- matrix(summ.m3[2:(M*kz1+1),1], nrow=M, kz1, byrow=T)
		est.si2 <- matrix(summ.m3[(M*kz1+2):(M*kz1+M*kz2*J+1),1], nrow=M*J, kz2, byrow=T)
		est.si1.sd <- matrix(summ.m3[2:(M*kz1+1),2], nrow=M, kz1, byrow=T)
		est.si2.sd <- matrix(summ.m3[(M*kz1+2):(M*kz1+M*kz2*J+1),2], nrow=M*J, kz2, byrow=T)




		###########################################################################################
		# the following code implements the penalized FR approach using the level 1 PC loadings 
		# estimated via MCMC
		#
		# note that we include two beta functions and two outcome variances in the following.
		###########################################################################################

		## generate outcomes from true Beta function and functions measured without error
		errors=rnorm(n, 0, sqrt(vareps))
		outcomes1.1 <- sapply(1:n, function(u) sum(funcs[u,]*trueBeta1)*by)+errors
		outcomes2.1 <- sapply(1:n, function(u) sum(funcs[u,]*trueBeta2)*by)+errors

		outcomes1.2 <- sapply(1:n, function(u) sum(funcs[u,]*trueBeta1)*by)+sqrt(2)*errors
		outcomes2.2 <- sapply(1:n, function(u) sum(funcs[u,]*trueBeta2)*by)+sqrt(2)*errors
	
		# compare estimated function to the true subject-specific functions; note that the
		# PC loadings estimated via MCMC are for centered functions, so we correct for that.

		phi1=t(phi1)

#		est.si1 = est.si1 + matrix(mu %*% phi1 * by, nrow=M, kz1, byrow=T)

		kz=kz1
		kb=min(kz, 35)

		num=kb-2
		qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
		knots <- quantile(t, qtiles)
		Basis = cbind(1, t, t^2, sapply(knots, function(k) ((t - k > 0) * (t - k)) ^ 2))
		##Basis = cbind(1, t, sapply(knots, function(k) ((t - k > 0) * (t - k))))

		## set up the matrices to be fit in the function lme()
		C <- est.si1
		J.mat <- t(phi1) %*% Basis * (by)
		CJ <- C %*% J.mat

		w <- cbind(1, CJ)
		X=as.matrix(w[,1:4])
		Z=as.matrix(w[,5:dim(w)[2]])
	
		group<-rep(1, n)
		
		# fit the PFR model
		fit.multlev.1.1=try(lme(outcomes1.1~-1+X, random=list(group=pdIdent(~-1+Z))))
		fit.multlev.2.1=try(lme(outcomes2.1~-1+X, random=list(group=pdIdent(~-1+Z))))

		fit.multlev.1.2=try(lme(outcomes1.2~-1+X, random=list(group=pdIdent(~-1+Z))))
		fit.multlev.2.2=try(lme(outcomes2.2~-1+X, random=list(group=pdIdent(~-1+Z))))

	
		if(class(fit.multlev.1.1) == "try-error") {
			MSE_MULTLEV[i,1]=NA
		} else {
			## get the coefficient and betaHat estimates
			coefs1.1 <- c(fit.multlev.1.1$coef$fixed,unlist(fit.multlev.1.1$coef$random))
			betaHatMultLev1.1 <- Basis %*% coefs1.1[-1]

			MSE_MULTLEV[i,1]=mean((trueBeta1 - betaHatMultLev1.1)^2)
		}

		if(class(fit.multlev.2.1) == "try-error") {
			MSE_MULTLEV[i,2]=NA
		} else {
			## get the coefficient and betaHat estimates
			coefs2.1 <- c(fit.multlev.2.1$coef$fixed,unlist(fit.multlev.2.1$coef$random))
			betaHatMultLev2.1 <- Basis %*% coefs2.1[-1]

			MSE_MULTLEV[i,2]=mean((trueBeta2 - betaHatMultLev2.1)^2)
		}
		
		if(class(fit.multlev.1.2) == "try-error") {
			MSE_MULTLEV[i,3]=NA
		} else {
			## get the coefficient and betaHat estimates
			coefs1.2 <- c(fit.multlev.1.2$coef$fixed,unlist(fit.multlev.1.2$coef$random))
			betaHatMultLev1.2 <- Basis %*% coefs1.2[-1]

			MSE_MULTLEV[i,3]=mean((trueBeta1 - betaHatMultLev1.2)^2)
		}
	
		if(class(fit.multlev.2.2) == "try-error") {
			MSE_MULTLEV[i,4]=NA
		} else {
			## get the coefficient and betaHat estimates
			coefs2.2 <- c(fit.multlev.2.2$coef$fixed,unlist(fit.multlev.2.2$coef$random))
			betaHatMultLev2.2 <- Basis %*% coefs2.2[-1]

			MSE_MULTLEV[i,4]=mean((trueBeta2 - betaHatMultLev2.2)^2)
		}


		###########################################################################################
		# the following code implements the PCR-PVE approach using the level 1 PC loadings 
		# estimated via MCMC
		###########################################################################################

		ro1=cumsum(fpca1.value)/sum(fpca1.value)
		K1=min((which(ro1>=.99)))

		C1.PVE=est.si1[,1:K1]
		fit.PVE1.1=lm(outcomes1.1~-1+C1.PVE)
		fit.PVE1.2=lm(outcomes1.2~-1+C1.PVE)
		fit.PVE2.1=lm(outcomes2.1~-1+C1.PVE)
		fit.PVE2.2=lm(outcomes2.2~-1+C1.PVE)

		# estimate the coefficient functions
		coefs.PVE1.1 <- coef(fit.PVE1.1)
		coefs.PVE2.1 <- coef(fit.PVE2.1)
		coefs.PVE1.2 <- coef(fit.PVE1.2)
		coefs.PVE2.2 <- coef(fit.PVE2.2)
		
		betaHatPVE1.1 <- fpca1.vectors[,1:K1] %*% coefs.PVE1.1[1:K1]
		betaHatPVE2.1 <- fpca1.vectors[,1:K1] %*% coefs.PVE2.1[1:K1]
		betaHatPVE1.2 <- fpca1.vectors[,1:K1] %*% coefs.PVE1.2[1:K1]
		betaHatPVE2.2 <- fpca1.vectors[,1:K1] %*% coefs.PVE2.2[1:K1]

		#calculate MSEs
		PVE_MULTLEV[i,1]=mean((trueBeta1 - betaHatPVE1.1)^2)
		PVE_MULTLEV[i,2]=mean((trueBeta2 - betaHatPVE2.1)^2)
		PVE_MULTLEV[i,3]=mean((trueBeta1 - betaHatPVE1.2)^2)
		PVE_MULTLEV[i,4]=mean((trueBeta2 - betaHatPVE2.2)^2)

	}

	print(i)
}





save(MSE_MULTLEV, PVE_MULTLEV, seeds, vareps, VarX, file="Multilevel_Sim_Results1.rda")




