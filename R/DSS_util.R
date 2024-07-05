########################################################################
## estimate dispersion for BS-seq data, given means
########################################################################
est.dispersion.BSseq <- function(X, N, estprob, ncores) {
    prior <- est.prior.BSseq.logN(X, N)
    dispersion.shrinkage.BSseq(X, N, prior, estprob, ncores)
}

########################################################################
## A function to estimate dipersion prior for BS-seq, assuming log-normal prior.
## It takes X and N, and only use the sites with big coverages,
## then return the mean and sd of prior distribution.
##
## For single rep data: will use logN(-3,1) as prior.
########################################################################
est.prior.BSseq.logN <- function(X, N) {
    ## rowMeans = DelayedArray::rowMeans
    ## rowSums = DelayedArray::rowSums
    rowVars_ <- function(x){rowVars(x, useNames=TRUE)}

    if(ncol(X) == 1) ## single rep
        return(c(-3, 1))

    ## keep sites with large coverage and no missing data
    ix = rowMeans(N>10)==1 & rowSums(N==0)==0
    if(sum(ix) < 50) {
        warning("The coverages are too low. Cannot get good estimations of prior. Use arbitrary prior N(-3,1).")
        return(c(-3, 1))
    }

    X = X[ix,,drop=FALSE]
    N = N[ix,,drop=FALSE]
    ## compute sample mean/var
    p = X/N
    mm = rowMeans(p)
    mm[mm==0] = 1e-5
    mm[mm==1] = 1-1e-5
    vv = rowVars_(p)
    phi = vv/mm/(1-mm)
    ## exclude those with vv==0. Those are sites with unobservable phis.
    ## But this will over estimate the prior.
    ## What will be the consequences????
    phi = phi[vv>0]
    lphi = log(phi[phi>0])
    prior.mean = median(lphi, na.rm=TRUE)
    prior.sd = IQR(lphi, na.rm=TRUE) /1.39

    ## It seems this over-estimates the truth. Need to use the tricks in
    ## my biostat paper to remove the over-estimation. To be done later.
    c(prior.mean, prior.sd)
}

########################################################################
## Dispersion shrinkage based on log-normal penalized likelihood.
## Takes X, N, estimated mean and prior.
##
## The shrinakge is done in log scale. So data will be shrink to the
## logarithmic means.
########################################################################
dispersion.shrinkage.BSseq <- function(X, N, prior, estprob, ncores) {
    ## penalized likelihood function
    plik.logN <- function(size, X,mu,m0,tau,phi)
        -(sum(dbb(size, X, mu, exp(phi))) + dnorm(phi, mean=m0, sd=tau, log=TRUE))

    ## for CG sites with no coverage, use prior
    shrk.phi=exp(rep(prior[1],nrow(N)))

    ## deal with estprob, make it a matrix if not.
    if(!is.matrix(estprob))
        estprob <- as.matrix(estprob)

    ## skip those without coverage
    ix <- rowSums(N>0) > 0
    X2 <- X[ix, ,drop=FALSE]; N2 <- N[ix,,drop=FALSE]; estprob2 <- estprob[ix,,drop=FALSE]
    shrk.phi2 <- rep(0, nrow(X2))

    ## setup a progress bar
    nCG.pb = round(nrow(X2)/100)
    if(ncores == 1) { ## use single core, no parallelism
        pb <- txtProgressBar(style = 3)
        for(i in 1:nrow(X2)) {
            ## print a progress bar
            if((i %% nCG.pb) == 0)
                setTxtProgressBar(pb, i/nrow(X2))
            ## I can keep the 0's with calculation. They don't make any difference.
            shrk.one=optimize(f=plik.logN, size=N2[i,], X=X2[i,], mu=estprob2[i,], m0=prior[1], tau=prior[2],
            interval=c(-5, log(0.99)),tol=1e-3)
            shrk.phi2[i]=exp(shrk.one$minimum)
        }
        setTxtProgressBar(pb, 1)
        cat("\n")
    } else  { ## use multiple cores.
        foo <- function(i) {
            shrk.one <- optimize(f=plik.logN, size=N2[i,], X=X2[i,],
                                 mu=estprob2[i,], m0=prior[1], tau=prior[2],
                                 interval=c(-5, log(0.99)),tol=1e-3)
            exp(shrk.one$minimum)
        }
        shrk.phi2 <- try(mclapply(1:nrow(X2), foo, mc.cores=ncores))
        if(inherits(shrk.phi2, "try-error")) {
            stop(paste("Shrinkage estimator in parallel computing encountered errors.",
                       "Please switch to single core.\n"))
        } else {
            shrk.phi2 <- unlist(shrk.phi2)
        }
    }

##         shrk.phi2 = foreach ( i=1:nrow(X2), .combine=c )  %dopar% {
##             ## print a progress bar
##             if((i %% nCG.pb) == 0)
##                 setTxtProgressBar(pb, i/nrow(X2))
##             ## I can keep the 0's with calculation. They don't make any difference.
##             shrk.one = optimize(f=plik.logN, size=N2[i,], X=X2[i,], mu=estprob2[i,], m0=prior[1], tau=prior[2],
##             interval = c(-5, log(0.99)),tol=1e-3)
##             exp(shrk.one$minimum)
##         }
##         stopImplicitCluster()
##         setTxtProgressBar(pb, 1)
##         cat("\n")

    shrk.phi[ix] <- shrk.phi2

    return(shrk.phi)
}

#########################################################
## beta-binomial (BB) density function.
## The BB distribution is parametrized by mean and dispersion.
#########################################################
dbb <- function (size, x, mu, phi, log=TRUE)  {
    ## 'size' and/or 'x' could be DelayedArray objects so turn them into
    ## ordinary arrays
    size=as.array(size)
    x=as.array(x)
    ## first convert mu/phi to alpha/beta
    tmp=1/phi-1
    alpha=mu*tmp
    beta=tmp - alpha
    v=lchoose(size,x)-lbeta(beta, alpha)+lbeta(size-x + beta,x+alpha)
    if(!log)
        return(exp(v))
    else return(v)
}
