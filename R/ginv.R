ginv <- function (X, tol = sqrt(.Machine$double.eps)) 
{
    # Based on the ginv function from the package MASS
    if (length(dim(X)) > 2L || !(is.numeric(X))) 
        stop("'X' must be a numeric matrix.")
    if (!is.matrix(X))  X <- as.matrix(X)
    Xsvd <- svd(X)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
}
