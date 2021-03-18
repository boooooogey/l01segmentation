fusedsegmentation <- function(y, lambda2, weight = NULL, l = 0, objective = "squared", format = "compressed", maxSegLength=3000, averageRangeLength=4){
    if (!(l %in% c(1,0))){
        stop("l should be either 0 or 1(for L0, L1 penalties respectively).")
    }
    if (!(objective %in% c("squared", "poisson"))){
        stop("objective should be either squared or poisson(for Squared, Poisson errors respectively).")
    }
    if (!(format %in% c("compressed", "full"))){
        stop("format should be either compressed or full.")
    }
    if(objective == "squared"){
        if(l == 0){
            if (format == "compressed"){
                breakpoints = L0SqrtErrBreakPoints(y,lambda2,weight,maxSegLength,averageRangeLength)
            }
            else{
                signal = L0SqrtErrSeg(y,lambda2,weight,maxSegLength,averageRangeLength)
            }
        }
        else{
            if (format == "compressed"){
                breakpoints = L1SqrtErrBreakPoints(y,lambda2,weight)
            }
            else{
                signal = L1SqrtErrFil(y,lambda2,weight)
            }
        }
    }
    else if(objective == "poisson"){
        if(l == 0){
            if (format == "compressed"){
                breakpoints = L0PoisBreakPoints(y,lambda2,weight,maxSegLength,averageRangeLength)
                breakpoints$val = exp(breakpoints$val)
            }
            else{
                signal = exp(L0PoisErrSeg(y,lambda2,weight,maxSegLength,averageRangeLength))
            }
        }
        else{
            stop("Not implemented.")
        }
    }
    if( format == "compressed"){
        N = length(breakpoints)
        starts = c(1,breakpoints$ii+1)
        ends = c(breakpoints$ii,length(y))
        return(data.frame(start=starts,end=ends,value=breakpoints$val))
    }
    else{
        return(signal)
    }
}

