l1squaredsegments = function(y, lambda2, weight = NULL, chrom="chr"){
    breakpoints = L1SqrtErrBreakPoints(y,lambda2,weight)
    N = length(breakpoints)
    starts = c(1,breakpoints+1)
    ends = c(breakpoints,length(y))
    return(data.frame(chrom=rep(chrom,N+1),chromStart=starts,chromEnd=ends))
}
