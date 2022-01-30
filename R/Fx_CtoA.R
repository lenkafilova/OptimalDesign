Fx_CtoA <- function(Fx, h=NULL, echo=TRUE) {
    # Transformation of Fx for regularized c-optimality
    # For simplicity, the regularization weight is fixed to 0.05 
    
    cl <- match.call()
    verify(cl, Fx = Fx, h = h, echo = echo)
    m <- ncol(Fx)

    if (is.null(h)) {
        if (echo) message("Setting h <- c(0,...,0,1).")
        h <- c(rep(0, m - 1), 1)
    }

    alpha <- 0.05 # This version of the package uses a fixed regularization parameter    
    L <- alpha*diag(m) + (1 - alpha)*m*h %*% t(h)/sum(h^2)
    
    if (rcond(L) < sqrt(.Machine$double.eps))
        warning("The problem of regularized c-optimality may be badly conditioned.")
    
    return(as.matrix(Fx %*% solve(chol(L))))
}
