formatgroupsegments <- function(out, Y, w, n, p){
    nbreak = length(out$ii)
    wa = (c(out$ii[2:nbreak], n) - out$ii)
    
    values = matrix(nrow=p,ncol=nbreak+1)
    starts = rep(0,nbreak+1)
    ends = rep(0,nbreak+1)
    
    starts[1] = 1
    ends [1]= c(out$ii[1])
    values[,1] = rowMeans(Y) - colSums(apply(t(out$B) * w[out$ii],2,cumsum)*wa)/n
    
    for(i in 2:(nbreak+1)){
        starts[i] = ends[i-1] + 1
        if( i != (nbreak+1)){
            ends[i] = out$ii[i]
        }
        else{
            ends[i] = n
        }
        values[,i] = values[,i-1] + out$B[,i-1] * w[out$ii[i-1]]
    }
    
    out = data.frame(start=starts,end=ends)
    for(i in 1:dim(values)[1]){
        out[[paste0("v",i)]] = values[i,]
    }
    
    return(out)
}

groupfusedsegmentation <- function(Y, lambda, w = NULL){
    Yhat = Y - rowMeans(Y)
    p = dim(Y)[1]
    n = dim(Y)[2]
    if (lambda <= 0){
        stop("Lambda should be a positive value.")
    }
    if (is.null(w)){
        w = 1:(n-1)
        w = sqrt(n/(w*(n-w)))
    }
    else{
        if(length(w) != dim(Y)[2]-1){
            stop("The number of weights should be equal to n-1 (n is the number of columns that Y has).")
        }
    }
    out = blockcoordinatedescent(Yhat,lambda,w) 
    return(formatgroupsegments(out,Y,w,n,p))
}
