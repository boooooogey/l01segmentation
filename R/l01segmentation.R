fusedsegmentation <- function(y, lambda2 = NULL, C = NULL, N = NULL, weight = NULL, l = 0, objective = "poisson", format = "compressed"){
    if (!(l %in% c(1,0))){
        stop("l should be either 0 or 1(for L0, L1 penalties respectively).")
    }
    if (!(objective %in% c("gaussian", "poisson", "exponential", "binomial"))){
        stop("objective should be either gaussian or poisson or exponential or binomial(for Gaussian, Poisson, exponential, binomial errors respectively).")
    }
    if (!(format %in% c("compressed", "full", "single"))){
        stop("format should be either compressed or full or single.")
    }
    if( format != "single" ){
        if(is.null(N) && is.null(lambda2)){
            stop("Either N or lambda2 should be given.")
        }
        if(is.null(lambda2) && l == 1){
            stop("N is for segmentation not L1.")
        }
        if(!is.null(N)){
            if(!is.integer(N)){
                stop("N should be an integer.")
            }
            if(N <= 0){
                stop("N should be positive.")
            }
        }
        if(!is.null(lambda2)){
            if (any(lambda2 <= 0)){
                stop("lambda values should be positive.")
            }
            if (!(length(lambda2) %in% c(1, length(y)-1))){
                stop("Wrong number of lambda values given.")
            }
            if (length(lambda2) == 1){
                lambda2 = rep(1,length(y)-1) * lambda2
            }
        }
    }
    if(!is.null(weight)){
        if (any(weight <= 0)){
            stop("weights should be positive.")
        }
        if (!(length(weight) %in% c(1, length(y)))){
            stop("Wrong number of weights given.")
        }
        if (length(weight) == 1){
            weight = rep(1, length(y)) * weight
        }
    }
    else{
        weight = rep(1, length(y))
    }
    if(objective == "gaussian"){
        if(l == 0){
            if (format == "compressed"){
                if(is.null(N)){
                    breakpoints = L0GaussianApproximateCondensed(y, lambda2, weight)
                    breakpoints$start = rev(breakpoints$start)
                    breakpoints$end = rev(breakpoints$end)
                    breakpoints$value = rev(breakpoints$value)
                }
                else{
                    breakpoints = L0GaussianApproximateNCondensed(y, N, weight)
                    breakpoints$start = rev(breakpoints$start)
                    breakpoints$end = rev(breakpoints$end)
                    breakpoints$value = rev(breakpoints$value)
                }
            }
            else{
                if(is.null(N)){
                    signal = L0GaussianApproximate(y, lambda2, weight)
                }
                else{
                    signal = L0GaussianApproximateN(y, N, weight)
                }
            }
        }
        else{
            if (format == "compressed"){
                breakpoints = L1GaussianApproximateCondensed(y,lambda2,weight)
                breakpoints$start = c(0, breakpoints$ii)
                breakpoints$end = c(breakpoints$ii, length(y))
                breakpoints$value = breakpoints$val
            }
            else{
                signal = L1GaussianApproximate(y,lambda2,weight)
            }
        }
    }
    else if(objective == "poisson"){
        if(l == 0){
            if (format == "compressed"){
                if(is.null(N)){
                    breakpoints = L0PoissonApproximateCondensed(y, lambda2, weight)
                    breakpoints$start = rev(breakpoints$start)
                    breakpoints$end = rev(breakpoints$end)
                    breakpoints$value = rev(breakpoints$value)
                }
                else{
                    breakpoints = L0PoissonApproximateNCondensed(y, N, weight)
                    breakpoints$start = rev(breakpoints$start)
                    breakpoints$end = rev(breakpoints$end)
                    breakpoints$value = rev(breakpoints$value)
                }
            } else if (format == "single"){
                    singlepoint = L0PoissonBreakPoint(y, weight)
                    breakpoints = list()
                    breakpoints$start = c(0, singlepoint+1)
                    breakpoints$end = c(singlepointi+1, length(y)) 
                    breakpoints$value = c(mean(y[1:(singlepoint)]), mean(y[(singlepoint+1):length(y)]))
            } else{
                if(is.null(N)){
                    signal = L0PoissonApproximate(y, lambda2, weight)
                }
                else{
                    signal = L0PoissonApproximateN(y, N, weight)
                }
            }
        }
        else{
            if (format == "compressed"){
                breakpoints = L1PoissonApproximateCondensed(y, weight, lambda2[1])
                breakpoints$start = unique(c(0, rev(breakpoints$ii) + 1 ))
                breakpoints$end = unique(c(rev(breakpoints$ii) + 1, length(y)))
                breakpoints$value = exp(rev(breakpoints$val))
            }
            else{
                signal = exp(L1PoissonApproximate(y, weight, lambda2[1]))
            }
        }
    }
    else if(objective == "exponential"){
        if(l == 0){
            if (format == "compressed"){
                if(is.null(N)){
                    breakpoints = L0ExponentialApproximateCondensed(y, lambda2, weight)
                    breakpoints$start = rev(breakpoints$start)
                    breakpoints$end = rev(breakpoints$end)
                    breakpoints$value = rev(breakpoints$value)
                }
                else{
                    breakpoints = L0ExponentialApproximateNCondensed(y, N, weight)
                    breakpoints$start = rev(breakpoints$start)
                    breakpoints$end = rev(breakpoints$end)
                    breakpoints$value = rev(breakpoints$value)
                }
            }
            else{
                if(is.null(N)){
                    signal = L0ExponentialApproximate(y, lambda2, weight)
                }
                else{
                    signal = L0ExponentialApproximateN(y, N, weight)
                }
            }
        }
        else{
            stop("Not implemented.")
        }
    }
    else if(objective == "binomial"){
        if( any( y > C ) ){
            stop("y (methylation) cannot be greater than C (coverage).")
        }
        if( any( y < 0 ) || any( C < 0 ) ){
            stop("y and C must be non negative.")
        }
        if( any(C == 0) ){
            stop("Missing data: C (coverage) is 0 for some indices.")
        }
        if(l == 0){
            y = as.vector(rbind(y,C))
            if (format == "compressed"){
                if(is.null(N)){
                    breakpoints = L0BinomialApproximateCondensed(y, lambda2, weight)
                    breakpoints$start = rev(breakpoints$start)
                    breakpoints$end = rev(breakpoints$end)
                    breakpoints$value = rev(breakpoints$value)
                }
                else{
                    breakpoints = L0BinomialApproximateNCondensed(y, N, weight)
                    breakpoints$start = rev(breakpoints$start)
                    breakpoints$end = rev(breakpoints$end)
                    breakpoints$value = rev(breakpoints$value)
                }
            }
            else{
                if(is.null(N)){
                    signal = L0BinomialApproximate(y, lambda2, weight)
                }
                else{
                    signal = L0BinomialApproximateN(y, N, weight)
                }
            }
        }
        else{
            if (format == "compressed"){
                breakpoints = L1BinomialApproximateCondensed(y, C, lambda2)
                breakpoints$start = unique(c(0, rev(breakpoints$ii) + 1 ))
                breakpoints$end = unique(c(rev(breakpoints$ii) + 1, length(y)))
                breakpoints$value = rev(breakpoints$val)
            }
            else{
                signal = L1BinomialApproximate(y, C, lambda2)
            }
        }
    }
    if( format == "full"){
        return(signal)
    }
    else{
        return(data.frame(start=breakpoints$start+1,end=breakpoints$end,value=breakpoints$value))
    }
}
