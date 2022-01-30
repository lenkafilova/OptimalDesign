od_D_VDM <- function(Fx, w1, eff, it.max, t.max, track) {
    
    start <- as.numeric(proc.time()[3])
    Fx <- as.matrix(Fx); n <- nrow(Fx); m <- ncol(Fx)
    n.iter <- 0; one <- rep(1, m); next.sec <- 0
    
    if (track) {
        info <- paste("Running od_D_VDM for cca", t.max, "seconds")
        info <- paste(info, " starting at ", Sys.time(), ".", sep = "")
        print(info, quote = FALSE)
        info <- paste("The problem size is n=", n, sep = "")
        info <- paste(info, " and m=", m, ".", sep = "")
        print(info, quote = FALSE)
    }
    
    if (is.null(w1)) w1 <- od_PIN(Fx, echo = FALSE)$w.pin/m
    w <- w1; M <- crossprod(sqrt(w) * Fx)
    d.fun <- as.vector(((Fx %*% t(chol(solve(M))))^2) %*% one)
    
    while (TRUE) {
        n.iter <- n.iter + 1; w <- n.iter * w/(n.iter + 1)
        i.max <- which.max(d.fun); w[i.max] <- w[i.max] + 1/(n.iter + 1)
        M <- crossprod(sqrt(w) * Fx)
        d.fun <- as.vector(((Fx %*% t(chol(solve(M))))^2) %*% one)
        eff.act <- m/max(d.fun)
        tm <- as.numeric(proc.time()[3]) - start    
        if (track && tm > next.sec) {
            print(paste("od_D_VDM Time:", round(tm, 1), 
                        "Efficiency:", round(eff.act, 9)))
            next.sec <- ceiling(tm)
        }
        if (eff.act >= eff || n.iter >= it.max || tm > t.max) break
    }
    
    t.act <- round(as.numeric(proc.time()[3]) - start, 2)
    Phi.best <- optcrit(Fx, w, echo = FALSE)
    
    if (track) {
        info <- paste("od_D_VDM finished after", 
            t.act, "seconds at", Sys.time())
        info <- paste(info, "with", n.iter, "iterations.")
        print(info, quote = FALSE)
        print(paste("D-criterion value:", Phi.best), quote = FALSE)
        print(paste("Efficiency at least:", eff.act), quote = FALSE)
    }
    
    return(list(w.best = w, Phi.best = Phi.best, eff.best = eff.act, 
        n.iter = n.iter, t.act = t.act))
}
