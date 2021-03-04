l0squaredsegments = function(y, lambda2, weight = NULL, maxSegLength=3000, averageRangeLength=4, chrom="chr"){
    breakpoints = L0SqrtErrBreakPoints(y,lambda2,weight,maxSegLength,averageRangeLength)
    N = length(breakpoints)
    starts = c(1,breakpoints+1)
    ends = c(breakpoints,length(y))
    return(data.frame(chrom=rep(chrom,N+1),chromStart=starts,chromEnd=ends))
}
