od_c_LP <- function(Fx, h, track) {

    start <- as.numeric(proc.time()[3])
    Fx <- as.matrix(Fx); n <- nrow(Fx); m <- ncol(Fx)
    a <- c(1, rep(0, 2*n)); b3 <- c(rep(0, m), 1)
    A3 <- rbind(cbind(-h, t(Fx), -t(Fx)), c(0, rep(1, 2*n)))
    
    res <- lpSolve::lp(direction = "max", objective.in = a, const.mat = A3, 
        const.dir = "=", const.rhs = b3)
    
    wx <- res$solution[2:(2*n + 1)]
    w <- wx[1:n] + wx[(n + 1):(2*n)]
    
    # Sometimes some components of w are very small negative numbers, therefore:
    w <- pmax(0, w); w <- w/sum(w)
    
    c.val <- optcrit(Fx, w, crit = "c", h = h, echo = FALSE)
    
    t.act <- round(as.numeric(proc.time()[3]) - start, 2)
    if (track) {
        info <- paste("Method od_c_LP finished after", t.act, "seconds at", 
            Sys.time(), ".")
        print(info, quote = FALSE)
        print(paste("c-criterion value:", c.val), quote = FALSE)
    }
    supp <- (1:n)[w > 0]; M <- infmat(Fx, w, echo = FALSE) 
    
    return( list(w.best = w, supp = supp, w.supp = w[supp], M.best = M,
                Phi.best = c.val, status = res$status, eff.best = 1,
                t.act = as.numeric(proc.time()[3]) - start))
}
