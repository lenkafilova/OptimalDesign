od_A_DEL <- function(Fx, w) {
     
    M <- infmat(Fx, w, echo = FALSE)
    if (rcond(M) < sqrt(.Machine$double.eps)) 
        warning("The information matrix may be badly conditioned.")
    n <- nrow(Fx); m <- ncol(Fx)
    M.inv <- solve(M); t <- sum(diag(M.inv))
    a.fun <- ((Fx %*% solve(M))^2) %*% rep(1, m)
    eps <- max(a.fun) - t + 1e-12
    del <- 1 + eps/t
    alpha <- min(eigen(M.inv, symmetric = TRUE)$values)/t
    gamma <- max(1, 1/del)
    B <- t * min(1, 1/del)
    l.bnd <- sqrt(alpha/gamma)
    u.bnd <- 1/sqrt(gamma)
    nu <- (1 - alpha)^3
    fn <- function(x) {
        alpha/x^2 + nu/(del - alpha * x)^2 - gamma
    }
    omega1 <- uniroot(fn, lower = l.bnd, upper = u.bnd)$root
    keep <- (1:n)[a.fun >= omega1^2 * B]
    
    return(keep)
}
