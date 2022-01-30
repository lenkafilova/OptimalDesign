od_pool <- function(X, val=NULL, pool.fun="sum", echo=TRUE) {
    # Create unique(X) and pool the vector val according to pool.fun
    # The elements of X are considered equal if they agree to 9 significant digits
    
    cl <- match.call()
    
    verify(cl, X = X, val = val, pool.fun = pool.fun, echo = echo)
    n <- nrow(X); k <- ncol(X)
    if (n > 100000) warning("Pooling large spaces may take a few minutes.")
    if (is.null(val)) val <- rep(1, n)
    X <- as.matrix(signif(X), 9)

    for (i in k:1) {
        ord <- order(X[, i], method = "radix")
        X <- X[ord, ]; val <- val[ord]
    }
    X <- as.matrix(X)
    
    # Future: the pools of the square roots of abs(val) could also be interesting
    # For val=dd, they correspond to standard deviations
    fun <- function(v) {
        if (pool.fun == "sum") return(sum(v))
        if (pool.fun == "max") return(max(v))
        if (pool.fun == "min") return(min(v))
        if (pool.fun == "mean") return(mean(v))
        if (pool.fun == "median") return(median(v))
        if (pool.fun == "0") return(0)
    }

    X.unique <- as.matrix(unique(X)); n.unique <- nrow(X.unique)
    val.pooled <- rep(0, n.unique); val.cut <- c(); j <- 1
    for (i in 1:n) {
        val.cut <- c(val.cut, val[i])
        if (i < n && any(X[i + 1, ] != X[i, ])) {
            val.pooled[j] <- fun(val.cut)
            val.cut <- c(); j <- j + 1
        }
    }
    val.pooled[j] <- fun(val.cut) 
    
    return(list(call = cl, X.unique = X.unique, val.pooled = val.pooled))
}
