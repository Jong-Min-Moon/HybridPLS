##########################################################################
# Nov 30 2009
# Jeff Goldsmith and Jennifer Feder
#
# this file implements a simulation study to evaluate the 
# functional regression technique in Goldsmith, Feder, Crainiceanu, Caffo,
# and Reich (2010). comparisons are made to several competing methods.
##########################################################################


library(nlme)
library(SemiPar)

seed.start <- 1001 
## use a different starting seed each time running the simulation
## we will keep track of the seeds at each iteration so that we 
## can replicate the results of a particular iteration

by = 0.1				# density grid on which functions are observed
t=seq(0,10, by=by)		# grid on which functions are observed
n=200 					# sample size
NREPS=1				# number of simulation iterations
N=15					# number of knots in covariance smoothing

## Parameters we want to vary
VarEps=c(.5, 1)			# variance of outcome from true mean
VarX = c(0, 1)			# measurement error variance
Kz= 30					# note: kb is taken as min(kz, 35)
kb=min(Kz, 35)

## seeds to be used at each generation of random functions (1) and 
## measurment error terms (1) and outcome error (1)
seeds <- matrix(seed.start + 0:(3*NREPS-1), nrow=NREPS, ncol=3)

## construct a data frame of all possible combinations of the
## parameters we vary in our simulation.
Params=data.frame(rep(VarEps, each = length(VarX)*length(Kz)), 
				rep(rep(VarX, each=length(Kz)), length(VarEps)), 
				rep(Kz, length=length(VarEps)*length(VarX)))
colnames(Params)=c("VarpEps", "VarX", "Kz")


## define true beta functions

trueBeta1=sin(t*pi/5)
trueBeta2=(t/2.5)^2
trueBeta3=-dnorm(t, mean=2, sd=.3)+3*dnorm(t, mean=5, sd=.4)+dnorm(t, mean=7.5, sd=.5)
	
trueBeta=cbind(trueBeta1, trueBeta2, trueBeta3)
trueBeta=trueBeta-matrix(colMeans(trueBeta), nrow=length(t), ncol=dim(trueBeta)[2], byrow=TRUE)
	

## create a list with all our random functions.

FUNCS0=vector("list", length=NREPS)
DECOMP0=vector("list", length=NREPS)
DECOMP0smooth=vector("list", length=NREPS)
for(i in 1:NREPS){
	set.seed(seeds[i,1])
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
	
	FUNCS0[[i]]=funcs
	
	varFuncs <- var(funcs)
	eigenDecomp <- eigen(varFuncs)

	# scale the vectors so they are eigen functions
	eigenDecomp$vectors <- eigenDecomp$vectors/sqrt(by)

	DECOMP0[[i]]=eigenDecomp
}


## these are lists for the functions with measurement error (var = 1)
## the PCA decomp takes place AFTER covariance smoothing

FUNCS1=vector("list", length=NREPS)
DECOMP1smooth=vector("list", length=NREPS)
for(i in 1:NREPS){
	set.seed(seeds[i,2])
	funcs = FUNCS0[[i]] + matrix(rnorm(length(t)*n, 0, sqrt(1)), nrow=n, ncol=length(t))
	
	FUNCS1[[i]]=funcs
	
	varFuncs <- var(funcs)

	## Smooth the covariance matrix
		## G2 is a M*M matrix for the raw covariance function
		G2 <- varFuncs
		M <- length(t)
		diag(G2)= rep(NA, M)
		g2 <- as.vector(G2)
		## define a N*N knots for bivariate smoothing
		N <- 15
		x1 <- rep(seq(0,10,length=N), N)
		x2 <- rep(seq(0,10,length=N), each=N)
		myknots <- data.frame(x1=x1, x2=x2)
		## bivariate smoothing using the spm function
		t1 <- rep(t, M)
		t2 <- rep(t, each=M)
		fit <- try(spm(g2 ~ f(t1, t2, knots=myknots), omit.missing=T))
		if(class(fit) == "try-error") {
			DECOMP1smooth[[i]] <- "failed"
		} else {
			newdata <- data.frame(t1=t1,t2=t2)
			pred <- predict(fit,newdata)
			####  make the fitted covariance function symmetric
			temp.g2 <- matrix(pred, M)
			fit.g2 <- (temp.g2 + t(temp.g2) )/2
		
			eigenDecomp <- eigen(fit.g2)
			eigenDecomp$vectors <- eigenDecomp$vectors/sqrt(by)
			DECOMP1smooth[[i]] <- eigenDecomp
		}
	print(i)
}
rm(eigenDecomp, fit.g2, temp.g2, pred, newdata, fit, t1, t2, x1, x2, myknots, g2, G2)


##################################################
## the following collection of functions implement
## the methods compared in our simulations. 
## we include unmodified code from cardot et al 
## and from reiss and ogden.
##################################################



##################################################
## PFR
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
	betaHat <- Basis %*% coefs[-(1)]
	
	D=diag(c(rep(0,ncol(X)), rep(1, ncol(Z))))
	lambda <- (fit$sigma)^2/as.numeric(VarCorr(fit)[1,1])
	sigmaEpsHat=(fit$sigma)^2
	varBeta=solve(t(w)%*% solve(diag((fit$sigma)^2, n)) %*% w + (1/as.numeric(VarCorr(fit)[1,1]))*D) [-(1),-(1)]
	varBetaHat=Basis%*%varBeta%*%t(Basis)
	
	ret <- list(fit, coefs, fitted, betaHat, w, X, Z, Basis, varBetaHat)
	names(ret) <- c("fit", "coefs", "fitted", "betaHat", "w", "X", "Z", "Basis", "varBetaHat")
	ret
}






##################################################
## here is the function that implements 
## the reiss and ogden function
##################################################

##source("~/Dropbox/Ciprian Project/Reiss Code/fpcr_code.r")

fpcr <- function(yy, signals, nc, nint=40, nfold=NULL, givenlam=NULL, s.index=NULL,
                 plotit=FALSE, pentype='deriv2', cv1=F) {

    require(splines)
    require(MASS)
    
    # Center the signals
    signals = scale(signals, scale=F)
       
    if (!exists('basisobj')) {
    	if (is.null(s.index)) stop('Must set s.index argument in order to define B-spline basis')
    	assign('basisobj', make.basis(s.index = s.index, nint), env=.GlobalEnv)
    }
    
    else if (basisobj$nb != nint + 3) stop("Number of intervals differs from that for existing B-spline basis object 'basisobj'; please remove or rename it before proceeding")
    
    bspline.basis = basisobj$bspline.basis
    knots = basisobj$knots
    nbasis = basisobj$nbasis
    s.index = basisobj$s.index
    
    if (!exists('halfpen') || ncol(halfpen)!=nbasis) assign('halfpen', make.halfpen(pentype, s.index = s.index, knots, nbasis), env=.GlobalEnv)

    if (max(nc) > nrow(halfpen)) stop(paste('Maximum number of components for the given basis is', nrow(halfpen)))
    
    if (((length(nc) > 1) || (cv1)) && is.null(nfold)) 
        stop('Must set nfold argument to specify k for k-fold cross-validation')

    if (!is.null(nfold) && max(nc) > 0) {
        vv <- array(0,c(nbasis,max(nc),nfold))
        pvv <- array(0,c(nrow(halfpen),dim(vv)[[2]],nfold))  
        rrr <- array(0,c(dim(vv)[[2]],dim(vv)[[2]],nfold))

        for (fold in 1:nfold) {
            trainers <- (1:length(yy))[(1:length(yy))%%nfold!=(fold-1)]
            sbtrain <- scale(signals[trainers,], scale=F) %*% bspline.basis  
            vv[,,fold] <- svd(sbtrain)$v[,1:max(nc)]
            pvv[,,fold] <- halfpen %*% vv[,,fold]
            rrr[,,fold] <- qr.R(qr(pvv[,,fold]))            
        }
    }

    #################################################
    # Find optimal number of components
    #################################################

    sumsspred <- rep(NA,length(nc))
    names(sumsspred) <- nc

    if (length(nc)==1 && cv1==F) bestnc <- nc
    # but if cv1==T, calculate CV (next section of code) even if length(nc)==1
 
    else {

        for (cc in 1:length(nc)) {
#            cat(nc[cc],'components\n')
            sspred <- rep(NA,nfold)

            for (fold in 1:nfold) {
                trainers <- (1:length(yy))[(1:length(yy))%%nfold!=(fold-1)]
                ytrain <- yy[trainers]
                yval <- yy[-trainers]
                straw <- signals[trainers,]
                strain <- scale(straw, scale=F)
                sbtrain <- strain %*% bspline.basis
                svraw <- signals[-trainers,]
                sval <- svraw + (strain-straw)[1:nrow(svraw),]

                if (nc[cc]==0) {x <- bspline.basis %*% Null(t(halfpen))
                                xx <- cbind(rep(1,length(ytrain)), strain %*% x)
                                z <- bspline.basis %*% ginv(halfpen)}

                else {x <- NULL
                      xx <- matrix(1,length(trainers),1)
                      z <- bspline.basis %*% vv[,1:nc[cc],fold] %*% solve(rrr[1:nc[cc],1:nc[cc],fold])}

                zz <- strain %*% z

                ##########################################

                # Given that number of components = nc[cc], estimate SS.pred by cross-validation
                bestlam <- if (!is.null(givenlam)) givenlam else penmod(ytrain, x, xx, z, zz)$lam                
                mm <- penmod(ytrain, x, xx, z, zz, bestlam)
                sspred[fold] <- crossprod(yval - mm$int - sval %*% mm$fhat)
        
                ##########################################

            }  # fold loop

            # Add up SS.pred for each validation set
            sumsspred[cc] <- sum(sspred)

        }  # cc loop
   
#        cat('SS.pred for each number of components:\n')
#        print(sumsspred)
        bestnc <- nc[which.min(sumsspred)] 
#        cat('Chosen number of components is', bestnc, '\n')  
    } # end of 'else' clause (--> choice of number of components ends here)

    #################################################
    # Fit model given optimal number of components
    #################################################

    # Create design matrices
    if (bestnc==0) {x <- bspline.basis %*% Null(t(halfpen))
                    xx <- cbind(rep(1,length(yy)), signals %*% x)
                    z <- bspline.basis %*% ginv(halfpen)}

    else {
        sb.all <- signals %*% bspline.basis
        vv.all <- as.matrix(svd(sb.all)$v[,1:bestnc])
        pvv.all <- halfpen %*% vv.all
        rrr.all <- qr.R(qr(pvv.all))
        x <- NULL
        xx <- matrix(1,length(yy),1)
        z <- bspline.basis %*% vv.all %*% solve(rrr.all)
    }

    zz <- signals %*% z

    modfit <- penmod(yy, x, xx, z, zz, if (!is.null(givenlam)) givenlam else NULL)                                

    if (plotit==TRUE) {
  	    if (length(nc)==1 && cv1==T) lines(s.index, modfit$fhat, lwd=4)
  	    else plot(s.index, modfit$fhat,type='l')
    }
    ztz <- crossprod(zz)
    df <- try(sum(diag(solve(ztz + diag(as.numeric(modfit$lambda), nrow(ztz)), ztz))))

    list(fhat=modfit$fhat, nc=bestnc, nint=nint, lambda=modfit$lambda, sumsspred=sumsspred, df=df, intercept=modfit$intercept, fitted=modfit$intercept + signals %*% modfit$fhat)
}

make.basis <- function(s.index, nint) {
    smin <- min(s.index)
    smax <- max(s.index)
    ds <- (smax - smin)/nint
    knots <- seq(smin - 3 * ds, smax + 3 * ds, by = ds)
    bspline.basis <- splineDesign(knots, s.index, 4, 0 * s.index)
    nbasis <- ncol(bspline.basis)  
#    cat('Creating B-spline basis\n')
    list(bspline.basis=bspline.basis, knots=knots, nbasis=nbasis, s.index=s.index)
}

make.halfpen <- function(pentype, s.index, knots, nbasis) {
	if (pentype!='deriv2' && substr(pentype,1,4)!='diff') stop("pentype must be 'deriv2' or 'diff1', 'diff2', etc.")
    if (pentype=='deriv2') halfpen <- qr.R(qr(splineDesign(knots, s.index, 4, rep(2, length(s.index)))))[1:(nbasis-2), ]
    else if (substr(pentype,1,4)=='diff') 
        halfpen <- diff(diag(nbasis), diff=as.numeric(substr(pentype, nchar(pentype), nchar(pentype))))
#    cat('Creating penalty matrix\n')
    halfpen
}

#################################################################################
#################################################################################
# penmod
#################################################################################
# Fit a signal regression model with roughness penalty 
#
# If the smoothing parameter lambda is not supplied, it is estimated by REML
# using function lme from package nlme.
# Either way the coefficients are divided into "fixed" and "random" effects.
# This is an internal function called by fpcr.

penmod <- function(y0, x0, xx0, z0, zz0, lambda=NULL) {
    nfix <- ncol(xx0)

    # Find lambda and fhat by REML  
    if (is.null(lambda)) {
        require(nlme)
        lme.gp <- rep(1,length(y0))
        lmefit <- lme(y0~-1+xx0,random=list(lme.gp=pdIdent(~-1+zz0)))
        lambda <- exp(-2*unlist(lmefit$modelStruct))
        fixx <- lmefit$coef$fix
        rndm <- t(lmefit$coef$random$lme.gp) 
        if (nfix > 1) fhat <- x0 %*% lmefit$coef$fix[-1] + z0 %*% t(lmefit$coef$random$lme.gp)
        else fhat <- z0 %*% t(lmefit$coef$random$lme.gp) 
    }

    # Given lambda, find fhat
    else {
        CC <- cbind(xx0,zz0)
        inner <- crossprod(CC)
        diag(inner)[-(1:nfix)] <- diag(inner)[-(1:nfix)] + lambda
        blueblup <- solve(inner, t(CC) %*% y0)
        fixx <- blueblup[(1:nfix), 1]
        rndm <- blueblup[-(1:nfix), 1]
        fixedpart <- 0
        if (nfix > 1) fixedpart <-  x0 %*% as.matrix(fixx)[-1,]
        fhat <- fixedpart + z0 %*% rndm
    }

    list(lambda=lambda, fhat=fhat, intercept=mean(y0), fixx=fixx, rndm=rndm)
}









##################################################
## here is the function that implements 
## pve and cv approaches to FPCR
##################################################

## function to fit the model using PC regression where best K is chosen by CV
fglmPCR <- function(outcomes, funcs, chooseK="cv", eigenDecomp=NULL) {
	## get the eigen decomposition of the smoothed variance matrix
	varFuncs <- var(funcs)

	if(is.null(eigenDecomp)) {
		varFuncs <- var(funcs)

		eigenDecomp <- eigen(varFuncs)
		eigenDecomp$vectors <- eigenDecomp$vectors/sqrt(by)
	}
	
	## choose the optimal K
	if(chooseK=="pve") {
		ro <- cumsum(eigenDecomp$values)/sum(eigenDecomp$values)
		K <- min(which( ro >= 0.99 ))
	}
	if(chooseK=="cv") {
		Kmax <- 30
		cv <- vector("numeric", Kmax)
		for(i in 1:Kmax) {
			C <- by* funcs %*% eigenDecomp$vectors[ ,1:i ]
			fit=lm(outcomes~-1+C)
			hat <- influence(fit)$hat
			cv[i] <- sum(((outcomes - predict(fit))/(1-hat))^2)
		}
		K <- which(cv==min(cv, na.rm=T))
	}
	
	## fit the model
	C <- by* funcs %*% eigenDecomp$vectors[ ,1:K ]
	fit=lm(outcomes~-1+C)
	
	coefs <- coef(fit)
	fitted <- as.matrix(C[,1:length(coefs)]) %*% coefs
	betaHat <- eigenDecomp$vectors[ ,1:K ] %*% coefs
	
	ret <- list(fit, coefs, fitted, betaHat, K)
	names(ret) <- c("fit", "coefs", "fitted", "betaHat", "K")
	ret
}



##################################################
## here is the function that implements 
## Cardot's FPCR with GCV
##################################################



Splinemlfgcv <- function(Y, X, veclambda, veck, order = 4, m = 2)
{
# Functional Linear Model: Spline Estimates
# ESTIMATION DE L'ERREUR PAR Validation Croisee Generalisee
#
# Y: vecteur des n observations de la variable a expliquer
# X: matrice nxp des variables explicatives (p= nbre de points de mesure) 
# veck vecteur de nombre de noeuds
# veclambda: vecteur des parametres de lissage
# order: ordre des B-splines
# m: ordre de derivation utilise pour le calcul de la penalisation
#################################################################
# gcv: matrice des  erreurs GCV selon les valeurs de lambda
#      et du nombre de noeuds
# rss: matrice des sommes des carres des erreurs residuelles 
# df : denominateur de la fonction GCV: [1-1/n * tr(Hatmatrix)]^2
# af : fonction alpha estimee pour les parametres optimaux
# min: minimum de GCV
# yhat: estimation de Y pour les parametres optimaux
# s2hat: estimation de la variance du bruit = RSS/tr(I-Hatmatrix)
#################################################################
	n <- nrow(X)
	p <- ncol(X)
	q <- length(veclambda)
	k <- length(veck)
	vtest <- 10^5
	gcv <- matrix(nrow = k, ncol = q)
	RSS <- matrix(nrow = k, ncol = q)
	DegFree <- matrix(nrow = k, ncol = q)
	nomcv2 <- paste("", veclambda, sep = "")
	nomcv1 <- paste("", veck, sep = "")
	dimnames(gcv) <- list(nomcv1, nomcv2)
	Xc <- sweep(X,2,apply(X,2,mean))
	for(i in 1:k) {
# boucle dimension B-splines
		kk <- veck[i]
		res <- Bspline.ini(t(Xc), kk, order, m)
		A <- rbind(matrix(1,nrow=1,ncol=n),res$A)
		AtA <- A %*% t(A)/n
		Delta <- t(as.matrix(Y)) %*% t(A)/n
PenMat <- matrix(0,ncol=kk+order+1,nrow=kk+order+1)
PenMat[-1,-1] <- res$G
		for(l in 1:q) {
# boucle parametre de lissage
			lambda <- veclambda[l]
			Gamma <- AtA + lambda * PenMat
			Hatn <- t(A) %*% solve(Gamma) %*% A/n
			RSS[i, l] <- mean((Y - Hatn %*% Y)^2)
			DegFree[i, l] <- (1 - sum(diag(Hatn))/n)^2
			gcv[i, l] <- RSS[i, l]/DegFree[i, l]
			if(vtest > gcv[i, l]) {
				alpha0 <- Delta %*% solve(Gamma)
				af <- res$B %*% alpha0[2:(kk+order+1)]
				mu <- alpha0[1]
				vtest <- gcv[i, l]
				yhat <- Hatn %*% Y
				s2hat <- (n * RSS[i, l])/(n - sum(diag(Hatn)))
			}
		}
	}
pos.opt <- order(gcv)[1]
pos.rho.opt <- trunc(pos.opt/k) + 1
pos.k.opt <- pos.opt - (pos.rho.opt - 1)*k
	list(gcv = gcv, rho.opt=veclambda[pos.rho.opt],k.opt=veck[pos.k.opt],rss = RSS, df = DegFree, mu = mu, af = af, min = vtest, yhat = 
		yhat, s2hat = s2hat)
}
Bspline.ini <- function(X, nknot, order, m)
{
###########################
### ARGUMENTS
# X:  Data matrix (size pxn), n=number of curves and 
#         p=number of design points
# nknot:  number of knots
# order:  degree of the B-splines
# m:      order of the derivative used in the penalization
### VALUES
# A:      matrix of the coordinates of the curves
#         in the B-splines basis
# G:      Gram matrix of the B-splines derivatives of order m
# B:      matrix of the B-splines values at the design points
#  
###########################  
require(splines) 
       p <- nrow(X)
        n <- ncol(X)
        A <- matrix(ncol = n, nrow = order + nknot)
        x0 <- seq(0, 1, length = p)
        x <- seq(0, 1, length = 200)
        knot <- quantile(x, (1:nknot)/(nknot + 1))
        delta <- sort(c(rep(range(x), order), knot))
        ## Calcul de la matrice de Gram
        B <- spline.des(delta, x, order)$design
        xdiff <- diff(x, 1)
        DmBj <- spline.des(delta, x, order, derivs = rep(m, 200))$design
        G1 <- t(DmBj[-1,  ]) %*% (DmBj[-1,  ] * xdiff)
        G2 <- t(DmBj[-200,  ]) %*% (DmBj[-200,  ] * xdiff)
        G <- 0.5 * (G1 + G2)    ## Calcul des coordonnees
        B <- spline.des(delta, x0, order)$design
        A <- t(B) %*% X/p
        list(A = A, G = G, B = B)
}






## keep track of the MSEs of each of the competing approaches
	
MSE_OURS=vector("list", length=dim(Params)[1])
for(i in 1:dim(Params)[1]){
	MSE_OURS[[i]]=matrix(0, nrow=NREPS, ncol=dim(trueBeta)[2])
	colnames(MSE_OURS[[i]])=paste("MSE_trueBeta", 1:dim(trueBeta)[2], sep="")
}

MSE_REISS=vector("list", length=dim(Params)[1])
for(i in 1:dim(Params)[1]){
	MSE_REISS[[i]]=matrix(0, nrow=NREPS, ncol=dim(trueBeta)[2])
	colnames(MSE_REISS[[i]])=paste("MSE_trueBeta", 1:dim(trueBeta)[2], sep="")
}

MSE_PVE=vector("list", length=dim(Params)[1])
for(i in 1:dim(Params)[1]){
	MSE_PVE[[i]]=matrix(0, nrow=NREPS, ncol=dim(trueBeta)[2])
	colnames(MSE_PVE[[i]])=paste("MSE_trueBeta", 1:dim(trueBeta)[2], sep="")
}

MSE_CV=vector("list", length=dim(Params)[1])
for(i in 1:dim(Params)[1]){
	MSE_CV[[i]]=matrix(0, nrow=NREPS, ncol=dim(trueBeta)[2])
	colnames(MSE_CV[[i]])=paste("MSE_trueBeta", 1:dim(trueBeta)[2], sep="")
}

MSE_CARDOT=vector("list", length=dim(Params)[1])
for(i in 1:dim(Params)[1]){
	MSE_CARDOT[[i]]=matrix(0, nrow=NREPS, ncol=dim(trueBeta)[2])
	colnames(MSE_CARDOT[[i]])=paste("MSE_trueBeta", 1:dim(trueBeta)[2], sep="")
}

lBound=vector("list", length=dim(Params)[1])
for(k in 1:dim(Params)[1]){
	lBound[[k]]=vector("list", length=dim(trueBeta)[2])
	names(lBound[[k]]) <- paste("MSE_trueBeta", 1:dim(trueBeta)[2], sep="")
	for(j in 1:dim(trueBeta)[2]){
		lBound[[k]][[j]] <- matrix(NA, NREPS, length(t))
	}
}
uBound <- lBound

for(i in 1:NREPS){	
	for(k in 1:dim(Params)[1]){
				
		vareps = Params[k,1]
		varX = Params[k,2]
		kz = Params[k,3]
		
		if(varX==0){
			funcs=FUNCS0[[i]]
			smoothDecomp=NULL
		} else if(varX==1){
			funcs=FUNCS1[[i]]
			smoothDecomp=DECOMP1smooth[[i]]
		} 
		
		if(!is.null(smoothDecomp) & is.character(smoothDecomp)){
			MSE_OURS[[k]][i,]=NA
			MSE_REISS[[k]][i,]=NA
			MSE_PVE[[k]][i,]=NA
			MSE_CV[[k]][i,]=NA
			MSE_CARDOT[[k]][i,]=NA
		} else {		
		
			######################################################
			## Smoothing splines
			######################################################
			## choose kb as necessary (default to 35 REF: Ruppert paper)
			

			mseCurOurs=vector('numeric', dim(trueBeta)[2])
			mseCurReiss=vector('numeric', dim(trueBeta)[2])
			mseCurPVE=vector('numeric', dim(trueBeta)[2])
			mseCurCV=vector('numeric', dim(trueBeta)[2])
			mseCurCardot=vector('numeric', dim(trueBeta)[2])

			## use the same errors on the outcomes from each trueBeta
			set.seed(seeds[i,3])
			errors=rnorm(n, 0, sqrt(vareps))
			#par(mfrow=c(1,3))
			for(j in 1:dim(trueBeta)[2]){
				#################################################################
				## generate outcomes	
				outcomes <- sapply(1:n, function(u) sum(FUNCS0[[i]][u,]*trueBeta[,j])*by)+errors
				#################################################################
			
				
				## our approach:
				
				fitJJ = try( fglmJJ(outcomes, funcs, eigenDecomp=smoothDecomp))
					
				#plot(t, trueBeta[,j], type='l', lwd=2)
				#points(t, fitJJ$betaHat-mean(fitJJ$betaHat), type='l', lwd=2, col="red")
								
				if(class(fitJJ) == "try-error") {
					mseCurOurs[j] <- NA
				} else {
					##mseCurOurs[j]=mean((trueBeta[,j] - (fitJJ$betaHat-mean(fitJJ$betaHat)))^2)
					mseCurOurs[j]=mean((trueBeta[,j] - fitJJ$betaHat)^2)
					
					## confidence bands
					lBound[[k]][[j]][i,] <- with(fitJJ, betaHat - 1.96*sqrt(diag(varBetaHat)))
					uBound[[k]][[j]][i,] <- with(fitJJ, betaHat + 1.96*sqrt(diag(varBetaHat)))

				}
				

				## Reiss and Ogden:
				
				fpcrmod <- fpcr(outcomes, funcs, nc=c(1:10, seq(11,80, 4)), nint=80, nfold=5, s.index=t)

				##mseCurReiss[j]=mean((trueBeta[,j] - (10*fpcrmod$fhat-mean(10*fpcrmod$fhat)))^2)
				mseCurReiss[j]=mean((trueBeta[,j] - 10*fpcrmod$fhat)^2)
				
				#points(t, 10*fpcrmod$fhat-mean(10*fpcrmod$fhat), type='l', lwd=2, col="blue")
				
				rm(basisobj, halfpen)


				## PVE and CV:
				
				fitPCRpve <- fglmPCR(outcomes, funcs, "pve", eigenDecomp=smoothDecomp)
				fitPCRcv <- fglmPCR(outcomes, funcs, "cv", eigenDecomp=smoothDecomp)

				#points(t, fitPCRpve$betaHat-mean(fitPCRpve$betaHat), type='l', lwd=2, col="green")
				#points(t, fitPCRcv$betaHat-mean(fitPCRcv$betaHat), type='l', lwd=2, col="purple")
								
				##mseCurPVE[j]=mean((trueBeta[,j] - (fitPCRpve$betaHat-mean(fitPCRpve$betaHat)))^2)
				mseCurPVE[j]=mean((trueBeta[,j] - fitPCRpve$betaHat)^2)
				##mseCurCV[j]=mean((trueBeta[,j] - (fitPCRcv$betaHat-mean(fitPCRcv$betaHat)))^2)
				mseCurCV[j]=mean((trueBeta[,j] - fitPCRcv$betaHat)^2)
				
				## Cardot:
				## uses same choices of components 'veck' as in Reiss/Ogden approach 'nc'
				## I selected the choices for rho, i.e. 'veclambda' 
				## that seemed to work pretty well for our three true beta functions
				fpcrgcv <- Splinemlfgcv(outcomes, funcs, veclambda=unique(c(seq(1e-8,1e-7,1e-8), 
					seq(1e-7,1e-6,1e-7), seq(1e-6,1e-5,1e-6))), veck=c(1:10, seq(11,80, 4)))
				##mseCurCardot[j]=mean((trueBeta[,j] - (by*fpcrgcv$af-mean(by*fpcrgcv$af)))^2)
				mseCurCardot[j]=mean((trueBeta[,j] - by*fpcrgcv$af)^2)
				
						
			}

			MSE_OURS[[k]][i,]=mseCurOurs
			MSE_REISS[[k]][i,]=mseCurReiss
			MSE_PVE[[k]][i,]=mseCurPVE
			MSE_CV[[k]][i,]=mseCurCV
			MSE_CARDOT[[k]][i,]=mseCurCardot
		}

	}
	print(i)
}


save(MSE_OURS, MSE_REISS, MSE_PVE, MSE_CV, MSE_CARDOT, lBound, uBound, seed.start, file=paste("SmoothSim", seed.start,".rda", sep=""))