
##########################################################################
# Nov 13 2009
# Jeff Goldsmith
#
# this file implements a simulation study to evaluate the sparse 
# functional regression technique in Goldsmith, Feder, Crainiceanu, Caffo,
# and Reich (2010).
##########################################################################

library(nlme)
library(SemiPar)
library(splines)

seed.start = 500
NREPS=1000					# number of iterations in the simulation

## Parameters we want to vary
PARAMS=cbind(rep(c(0, 1), each=2), rep(c(TRUE, FALSE), 2))

## set the standard parameters
by = 0.1					# density of grid on which functions are observed
t=seq(0,10, by=by)			# grid on which functions are observed
kz=10						# number of PCs used to estimate the X_i(t)
kb=min(kz, 35)
n=200  						# sample size
NGRIDPTS=10					# number of sampled points for each function
N.fit=101 					# number of points in the grid for the estimated functions


## define true beta functions
trueBeta1=2*sin(t*pi/5)
trueBeta2=(t/2.5)^2

# keep track of the MSEs of the estimated functions
MSE_SPARSE=vector("list", length=dim(PARAMS)[1])
for(i in 1:dim(PARAMS)[1]){
	MSE_SPARSE[[i]]=matrix(0, ncol=4, nrow=NREPS)
}
MSE_PVE=MSE_SPARSE

# keep track of the seed for reproducibility
seeds = matrix(seed.start + 0:(NREPS-1), nrow=NREPS, ncol=1)

for(k in 1:dim(PARAMS)[1]){
	
	var.sparse=PARAMS[k,1]
	KNOWN=PARAMS[k,2]

	for(i in 1:NREPS){
	
		set.seed(seeds[i])
	
		## generate random functions (note - here we are not measuring with error)
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


		###############################################
		# the following is the sparse data estimating
		# code, written by chongzhi di. i haven't added any
		# comments to his code. some important variables:
		#
		# NGRIDPTS - number of observed points per function
		# N.fit - number of points in the estimated functions
		# y.hat - matrix of estimated functions
		# phi1 - estimated eigenfunctions
		# est.si1 - estimated PC loadings
		###############################################

		## sparse data stuff

		N=NGRIDPTS		
		M=n
		J=1
		tstart=0
		tlength=10
		N.smooth.single=20
		N.smooth=5


		unique.average <- function(index, y) {
    
    		n <- sum(!is.na(y))
		    y.new <- rep(NA, n)
	    	index.new <- rep(NA, n)
    		count <- rep(0, n)
    
     
		    y.new[1] <- y[1]
    		index.new[1] <- index[1] 
		    count[1] <- 1
	    	current <- 1
    
		    if(n > 1) {
    		    for(i in 2:n) {
        		    if(index[i]==index[i-1]) {
            	    	count[i] <- count[i-1] + 1
                		y.new[current] <- ( y.new[current] * count[i-1] + y[i] )/count[i]
	        	        count[i-1] <- count[i]
                   
        	    }
    	        else {
		                current <- current + 1
    	            	y.new[current] <- y[i] 
        	    	    count[i] <- 1
            		    index.new[current] <- index[i] 
	    	        }
    		    }
	    	}
    	
	    	return(list(index=index.new, y=y.new))
		}


		sampled=t(sapply(1:n, function(u) sort(sample(1:length(t), size=NGRIDPTS  ))))
		tx=t(sapply(1:n, function(u) t[sampled[u,]]))
		y=t(sapply(1:n, function(u) funcs[u,sampled[u,]]))

		y=y+rnorm(length(y), mean=0, sd=sqrt(var.sparse))


		t.knot <- seq(tstart, tstart + tlength, length=N.smooth.single)
		t.fit <- seq(tstart, tstart + tlength, length=N.fit)
    
		N.row <- dim(y)[1]
    
		y1 <- matrix(0, M, N)
		index <- matrix(0, M, N)
		t.fit.matrix <- matrix(rep(t.fit, each=N), nrow=N)


		for(m in 1:M ) {
			row.ind <- m
			index.temp <- apply( (t.fit.matrix - tx[row.ind, ])^2, 1, which.min )
			fit <- unique.average(index.temp, y[row.ind,] )
			y1[row.ind,] <- fit$y
			index[row.ind, ] <- fit$index
		}
    
		n.grid <- apply(!is.na(y1), 1, sum)


		library(SemiPar)
		mu.fit <- rep(0,  N.fit)
    
		mu.obs <- matrix(0, M, N)
    
		###  divide the interval into "N.fit" bins to reduce computational burden for smoothing 
		mu.sum <- rep(0, N.fit)
		mu.count <- rep(0, N.fit) 
		mu.mean <- rep(0, N.fit)
		for (m in 1:M) {
			row.ind <- m
			row.len <- n.grid[row.ind]
			mu.count[ index[row.ind, 1:row.len] ] <- mu.count[ index[row.ind, 1:row.len] ] + 1
			mu.sum[ index[row.ind, 1:row.len] ] <- mu.sum[ index[row.ind, 1:row.len] ] + y1[row.ind, 1:row.len ]
		} 
		mu.mean <- ifelse(mu.count==0, NA,  mu.sum/mu.count)    

    
		temp.data <- data.frame(y2 = mu.mean,  x2 = t.fit )
		attach(temp.data, pos=3)
		fit.mu <- spm( y2~f(x2,knots=t.knot) ,omit.missing=TRUE)
		newx <- data.frame(x2=t.fit)
		mu.fit <- predict(fit.mu,newx)
		detach(pos=3)
		mu.obs <- matrix( mu.fit[ index ], ncol=N)

		resd <- matrix(0, nrow=M, ncol=N) 
		resd <- y - mu.obs

		###  divide the interval into "N.fit" bins to reduce computational burden for smoothing 
		cov.sum <- matrix(0, N.fit, N.fit)
		cov.count <- matrix(0, N.fit, N.fit) 
		cov.mean <-  matrix(0, N.fit, N.fit)
		for (m in 1:M) {
			row.ind <- m
			len.ind <- n.grid[row.ind]
			temp.ind <- index[row.ind, 1:len.ind]
			cov.count[ temp.ind, temp.ind ] <- cov.count[ temp.ind, temp.ind ] + 1
			cov.sum[ temp.ind, temp.ind ] <- cov.sum[ temp.ind, temp.ind ] + resd[row.ind, 1:len.ind ] %*% t( resd[row.ind, 1:len.ind])
		} 
		cov.mean <- ifelse(cov.count==0, NA,  cov.sum/cov.count)    

		cov.mean.nodiag <- cov.mean
		diag(cov.mean.nodiag) <- rep(NA, N.fit)
		resd.pair <- cbind( as.vector(cov.mean.nodiag), rep(t.fit, each=N.fit), rep(t.fit, N.fit) )
		resd.pair.diag <- cbind( diag(cov.mean), t.fit) 



		###     Estimate the three covariance functions: overall covariance G, 
		###     between covariance Gb and within covariance Gw
		G <- matrix(0, N.fit, N.fit)

		### load the "SemiPar" package which contains functions for semiparametric smoothing 
		### using linear mixed models
		t.smooth <- seq(tstart, tstart + tlength ,length=N.smooth)
        
		myknots.t <- data.frame(x1t=rep(t.smooth, each=N.smooth), x2t=rep(t.smooth, N.smooth))         
		data.gt <- data.frame(gt=resd.pair[,1], x1t=resd.pair[,2], x2t=resd.pair[,3])
		attach(data.gt, pos=4)
		fit.t <- spm(gt ~ f(x1t, x2t, knots=myknots.t),omit.missing=TRUE)
		newdata.t <- data.frame(x1t=rep(t.fit, each=N.fit),x2t=rep(t.fit, N.fit))
		pred.t <- predict(fit.t, newdata.t)
		detach(pos=4)

		s.g <- matrix(pred.t, N.fit)
		G <- ( s.g + t(s.g) )/2

		eigenDecomp <- eigen(G)

		#par(mfrow=c(3,5))
		#for(i in 1:15) {
		#	plot(eigenDecomp$vectors[,i], type="l")
		#}


		temp.data <- data.frame(var.y = as.vector(resd^2),  var.x = as.vector(tx) )
		attach(temp.data, pos=6)
		fit.var <- spm( var.y ~ f(var.x,knots=t.knot) ,omit.missing=TRUE)
		newx <- data.frame(var.x=t.fit)
		var.fit <- predict(fit.var,newx)          
		detach(pos=6)

		var.noise <- mean( var.fit - diag(matrix(pred.t,N.fit)) ) 


		###############################################
		# here is a change to the original code. we found
		# that in the case of no measurement error, the
		# regression performed poorly. in those cases, we 
		# add noise with a known variance. here we coerce
		# the noise variance to be the known error variance.
		###############################################
	
		if(KNOWN){
			if(var.sparse==0){
				var.noise=.002
			} else {
				var.noise=var.sparse		
			}
		}
	
    	###############################################
		# that was the change. carry on.
		###############################################
    
    
    
	    fpca1.value <- eigenDecomp$values* tlength / N.fit

    	fpca1.value <- ifelse(fpca1.value>=0, fpca1.value, 0)

    
		K1=min(35, N)

		fpca1.vectors <- eigenDecomp$vectors[, 1:K1]*sqrt(N.fit/tlength)


		###     The eigenfunctions are unique only up to a change of signs.
		###     Select the signs of eigenfunctions so that the integration over the domain 
		###     is non-negative
    	for(i2 in 1:K1) {
	        v2 <- fpca1.vectors[,i2]
        	tempsign <- sum(v2)
    	    fpca1.vectors[,i2] <- ifelse(tempsign<0, -1,1) * v2
	    }


	
		lam1=fpca1.value[1:K1]
		phi1=fpca1.vectors



		###
		###     Estimatimate principal component scores and predict various functions
		###

		library(MASS) 

		index <- matrix(0, M, N)
		mu.hat <- matrix(0, M, N)    
		y.center <- matrix(0, M, N)

		phi1.hat <- array(0, dim=c( M, N, K1) )

		t.fit.matrix <- matrix(rep(t.fit, each=N), nrow=N)

		for(m in 1:M ) {
    	    row.ind <- m
        	index <- apply( (t.fit.matrix - tx[row.ind, ])^2, 1, which.min )
    	    y.center[row.ind, ] <- y[ row.ind ,] - mu.fit[ index ]
	        phi1.hat[row.ind, , ] <- phi1[index, ]
		}



		est.si1 <- matrix(0, nrow=M, ncol=K1)

		est.si1.sd <- matrix(0, nrow=M, ncol=K1)


		y.hat.center <- matrix(0, M, N.fit) 
		y.hat <- matrix(0, M, N.fit)
		y.hat.subject <- matrix(0, M, N.fit)

		y.hat.sd <- matrix(0, M, N.fit)
		y.hat.subject.sd <- matrix(0, M, N.fit)



		for(m in 1:M) {
			cov.y <- matrix(0, N, N)
			cov.xi.y <- matrix(0, K1,  N )

			ind1 <- m
			cov.y [ 1:N , 1:N ] <- phi1.hat[ind1,,] %*% diag( lam1) %*%  t( phi1.hat[ind1,,] )  + diag( rep(var.noise, N) )

			ind1 <- m
			cov.xi.y[1:K1, 1:N ] <- diag( lam1) %*% t( phi1.hat[ind1,,] )


			temp.mat <- cov.xi.y %*% ginv(cov.y)
			score <- temp.mat %*% as.vector( t(y.center[m,]) ) 
			var.score <- diag( c(lam1 ) ) -  temp.mat %*% t(cov.xi.y)

			est.si1[m ,] <- score[1:K1]

			y.hat.subject[m, ] <- mu.fit + phi1 %*% est.si1[m,] 
			y.hat.subject.sd[m,] <-  sqrt( diag( phi1 %*% var.score[1:K1, 1:K1] %*% t( phi1 ) ) )

			row.ind <- m
			y.hat.center[ row.ind ,] <-  phi1 %*% est.si1[m,]
			y.hat[row.ind,] <- y.hat.center[ row.ind ,]  +  mu.fit

			temp1 <- cbind( phi1)
			temp2 <- var.score[ c(1:K1) , c(1:K1) ]
			y.hat.sd[row.ind, ] <- sqrt( diag( temp1 %*% temp2 %*% t(temp1) ) )   
		}


		###############################################
		# end sparse data estimation procedure
		###############################################




		###############################################
		# use sparse data estimates in a functional 
		# regression.
		###############################################

		t=t.fit
		by.sparse=t[2]-t[1]

		## generate outcomes from true Beta function and functions measured without error
		errors1=rnorm(n, 0, sqrt(.5))
		errors2=rnorm(n, 0, sqrt(1))
		outcomes1.5 <- sapply(1:n, function(u) sum(funcs[u,]*trueBeta1)*by)+errors1
		outcomes2.5 <- sapply(1:n, function(u) sum(funcs[u,]*trueBeta2)*by)+errors1
		outcomes11 <- sapply(1:n, function(u) sum(funcs[u,]*trueBeta1)*by)+errors2
		outcomes21 <- sapply(1:n, function(u) sum(funcs[u,]*trueBeta2)*by)+errors2
	
		# note - the following code uses the eigenfunctions and
		# loadings from the sparse data code; the results are basically
		# the same as those from using the fglm function and the 
		# estimated functions

		num=kb-2
		qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
		knots <- quantile(t.fit, qtiles)
		Basis = cbind(1, t.fit, t.fit^2, sapply(knots, function(k) ((t.fit - k > 0) * (t.fit - k)) ^ 2))
		##Basis = cbind(1, t, sapply(knots, function(k) ((t - k > 0) * (t - k))))

		## set up the matrices to be fit in the function lme()
		C <- est.si1
		J <- t(phi1) %*% Basis * (by.sparse)
		CJ <- C %*% J

		w <- cbind(1, CJ)
		X=as.matrix(w[,1:4])
		Z=as.matrix(w[,5:dim(w)[2]])
	
		group<-rep(1, n)
							
		fit.sparse1.5=try(lme(outcomes1.5~-1+X, random=list(group=pdIdent(~-1+Z))))
		fit.sparse2.5=try(lme(outcomes2.5~-1+X, random=list(group=pdIdent(~-1+Z))))
		fit.sparse11=try(lme(outcomes11~-1+X, random=list(group=pdIdent(~-1+Z))))
		fit.sparse21=try(lme(outcomes21~-1+X, random=list(group=pdIdent(~-1+Z))))
	

		# make sure that the true betas we compare to are on the same grid
		# as the estimated functions
		trueBeta1.fit=2*sin(t.fit*pi/5)
		trueBeta2.fit=(t/2.5)^2


		if(class(fit.sparse1.5) == "try-error") {
			MSE_SPARSE[[k]][i,1]=NA
		} else {
			## get the coefficient and betaHat estimates
			coefs <- c(fit.sparse1.5$coef$fixed,unlist(fit.sparse1.5$coef$random))
			betaHatSparse1.5 <- Basis %*% coefs[-1]
		
			## get the MSE
			MSE_SPARSE[[k]][i,1]=mean((trueBeta1.fit - betaHatSparse1.5)^2)
		}
		
		if(class(fit.sparse2.5) == "try-error") {
			MSE_SPARSE[[k]][i,2]=NA
		} else {
			## get the coefficient and betaHat estimates
			coefs <- c(fit.sparse2.5$coef$fixed,unlist(fit.sparse2.5$coef$random))
			betaHatSparse2.5 <- Basis %*% coefs[-1]
		
			## get the MSE
			MSE_SPARSE[[k]][i,2]=mean((trueBeta2.fit - betaHatSparse2.5)^2)
		}
		
		if(class(fit.sparse11) == "try-error") {
			MSE_SPARSE[[k]][i,3]=NA
		} else {
			## get the coefficient and betaHat estimates
			coefs <- c(fit.sparse11$coef$fixed,unlist(fit.sparse11$coef$random))
			betaHatSparse11 <- Basis %*% coefs[-1]
		
			## get the MSE
			MSE_SPARSE[[k]][i,3]=mean((trueBeta1.fit - betaHatSparse11)^2)
		}
		
		if(class(fit.sparse21) == "try-error") {
			MSE_SPARSE[[k]][i,4]=NA
		} else {
			## get the coefficient and betaHat estimates
			coefs <- c(fit.sparse21$coef$fixed,unlist(fit.sparse21$coef$random))
			betaHatSparse21 <- Basis %*% coefs[-1]
		
			## get the MSE
			MSE_SPARSE[[k]][i,4]=mean((trueBeta2.fit - betaHatSparse21)^2)
		}

		# also try a PVE approach
	
		K <- 3

	
		fit.PVE.1.5=lm(outcomes1.5~-1+C[,1:K])
		betaHat.PVE.1.5 <- phi1[ ,1:K ] %*% coef(fit.PVE.1.5)

		fit.PVE.2.5=lm(outcomes2.5~-1+C[,1:K])
		betaHat.PVE.2.5 <- phi1[, 1:K] %*% coef(fit.PVE.2.5)

		fit.PVE.11=lm(outcomes11~-1+C[,1:K])
		betaHat.PVE.11 <- phi1[ ,1:K ] %*% coef(fit.PVE.11)

		fit.PVE.21=lm(outcomes21~-1+C[,1:K])
		betaHat.PVE.21 <- phi1[, 1:K] %*% coef(fit.PVE.21)


		MSE_PVE[[k]][i,1]=mean((trueBeta1.fit - betaHat.PVE.1.5)^2)
		MSE_PVE[[k]][i,2]=mean((trueBeta2.fit - betaHat.PVE.2.5)^2)
		MSE_PVE[[k]][i,3]=mean((trueBeta1.fit - betaHat.PVE.11)^2)
		MSE_PVE[[k]][i,4]=mean((trueBeta2.fit - betaHat.PVE.21)^2)


		# par(mfrow=c(1,2))
		# plot(trueBeta1.fit, type='l')
		# points(betaHatSparse11, type='l', col="red")
		# points(betaHat.PVE.11, type='l', col="blue")
		# plot(trueBeta2.fit, type='l')
		# points(betaHatSparse21, type='l', col="red")
		# points(betaHat.PVE.21, type='l', col="blue")

		if(round(i/10)==i/10){
			print(k)
			print(i)	
		}
	}
}


save(MSE_SPARSE, MSE_PVE, seeds, var.sparse, var.noise, KNOWN, N.fit, NGRIDPTS, file="Sparse_Sim_Results.rda")