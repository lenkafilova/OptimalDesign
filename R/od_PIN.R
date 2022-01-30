od_PIN <- function(Fx, alg.PIN="KYM", echo=TRUE) {
    # Creates an efficient saturated exact design to initialize computations
    
    cl <- match.call()
    verify(cl, Fx = Fx, alg.PIN = alg.PIN, echo = echo)
    n <- nrow(Fx); m <- ncol(Fx)
    
    start <- as.numeric(proc.time()[3])
    ind <- rep(0, m)
    if (alg.PIN == "GKM") {
        Fxx <- Fx; v2 <- (Fxx^2) %*% rep(1, m) 
        j <- which.max(v2); ind[1] <- j
        for (i in 2:m) {
            scp <- Fxx %*% Fxx[j, ]
            Fxx <- Fxx - scp %*% t(Fxx[j, ]) / v2[j]
            v2 <- v2 - scp^2 / v2[j]
            j <- which.max(v2); ind[i] <- j
        }
    } else {
        P <- diag(m)
        j <- which.max(abs(Fx %*% rnorm(m))); ind[1] <- j
        for (i in 2:m) {
            fx <- P %*% Fx[j,]
            P <- P - tcrossprod(fx) / sum(fx^2)
            j <- which.max(abs(Fx %*% (P %*% rnorm(m))))
            ind[i] <- j
        }
    }
    
    w <- rep(0, n); w[ind] <- 1
    t.act <- as.numeric(proc.time()[3]) - start
    if (echo)
        print(paste(alg.PIN, " finished in ", round(t.act, 2),
                    " seconds.", sep = ""), quote = FALSE)
    
    return(list(call = cl, w.pin = w, supp = ind, M.pin = infmat(Fx, w, echo = FALSE),
                Phi.D = optcrit(Fx, w, echo = FALSE), t.act = t.act))
}
