od_D_DEL <- function(Fx, w) {
    
    M <- infmat(Fx, w, echo = FALSE)
    if (rcond(M) < sqrt(.Machine$double.eps)) 
        warning("The information matrix may be badly conditioned.")
    n <- nrow(Fx); m <- ncol(Fx)
    M <- infmat(Fx, w, echo = FALSE)
    d.fun <- ((Fx %*% t(chol(solve(M))))^2) %*% rep(1, m)
    eps <- max(d.fun) - m + 1e-12
    thr <- m * (1 + eps/2 - sqrt(eps * (4 + eps - 4/m))/2)
    keep <- (1:n)[d.fun >= thr]
    
    return(keep)
}
