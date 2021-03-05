fusedsegmentation <- function(y, lambda2, weight = NULL, l = 0, objective = "squared", maxSegLength=3000, averageRangeLength=4){
    if (!(l %in% c(1,0))){
        stop("l should be either 0 or 1(for L0, L1 penalties respectively).")
    }
    if (!(objective %in% c("squared", "poisson"))){
        stop("objective should be either squared or poisson(for Squared, Poisson errors respectively).")
    }
    if(objective == "squared"){
        if(l == 0){
            breakpoints = L0SqrtErrBreakPoints(y,lambda2,weight,maxSegLength,averageRangeLength)
        }
        else{
            breakpoints = L1SqrtErrBreakPoints(y,lambda2,weight)
        }
    }
    else if(objective == "poisson"){
        if(l == 0){
            breakpoints = L0PoisBreakPoints(y,lambda2,weight,maxSegLength,averageRangeLength)
        }
        else{
            stop("Not implemented.")
        }
    }
    N = length(breakpoints)
    starts = c(1,breakpoints$ii+1)
    ends = c(breakpoints$ii,length(y))
    return(data.frame(start=starts,end=ends,value=breakpoints$val))
}

plotsegments <- function(data,segments,label=NULL,title="",ylab="",xlab="",show="fused"){
    if(!(show %in% c("fused","breakpoints"))){
        stop("Unknown data type (Options: \"fused\" - \"breakpoints\").")
    }
    if (is.null(label)){
        plot(data,ylab = ylab,xlab = xlab,main = title)
    }
    else{
        plot(data,col=label+2,ylab = ylab,xlab = xlab,main = title)
    }
    if(show == "fused"){
        for(i in 1:dim(segments)[1]){
            e = segments$end[i]
            s = segments$start[i]
            v = segments$value[i]
            if(s != e){
                lines(c(s,e), c(v,v),col=2)
                points(c(s,e),c(v,v),col=2,pch=6)
            }
            else{
                points(s,v,col=2,pch=6)
            }
            if( i < dim(segments)[1]){
                lines(c(e,segments$start[i+1]),c(v,segments$value[i+1]) ,col=2)
            }
        }
    }
    if(show == "breakpoints"){
        abline(v=segments$start,col=2,lty=2)
        abline(v=segments$end,col=2,lty=2)
    }
}

