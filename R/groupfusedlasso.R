groupFusedLasso <- function(Y, lambda, w = NULL){
    Yhat = Y - rowMeans(Y)
    if (lambda <= 0){
        stop("Lambda should be a positive value.")
    }
    if (is.null(w)){
        w = rep(1,dim(Y)[2]-1)
    }
    else{
        if(length(w) != dim(Y)[2]-1){
            stop("The number of weights should be equal to n-1 (n is the number of columns that Y has).")
        }
    }
    B = blockcoordinatedescent(Yhat,lambda,w) 
    x = X(w)
    BX = B %*% x 
    U = BX + rowMeans(Y - BX)
    return(U)
}
