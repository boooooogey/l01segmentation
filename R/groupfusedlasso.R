groupFusedLasso <- function(Y, lambda, w = NULL){
    Yhat = Y - rowMeans(Y)
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
    B = blockcoordinatedescent(Yhat,lambda,w) 
    BX = dotX(B,w)
    U = BX + rowMeans(Y - BX)
    return(U)
}
