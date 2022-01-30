od_A_KL_bin <- function(Fx, N, Phi.app, w1, K, L, rest.max, t.max, track) {
    
    n <- nrow(Fx); m <- ncol(Fx); one.m <- rep(1, m)
    next.sec <- 0; n.ex <- 0; n.rest <- 0
    
    if (is.null(Phi.app))
        # We will only compute a very conservative upper bound on the optimal value
        Phi.app <- N*od_REX(Fx, crit = "A", track = FALSE)$Phi.best
    
    start <- as.numeric(proc.time()[3])
    if (track) {
        info <- paste("Running od_A_KL_bin for cca", t.max, "seconds")
        info <- paste(info, " starting at ", Sys.time(), ".", sep = "")
        print(info, quote = FALSE)
        info <- paste("The problem size is n=", n, sep = "")
        info <- paste(info, ", m=", m, ", N=", N, sep = "")
        print(info, quote = FALSE)
    }
    
    if (is.null(K)) 
        K <- max(10, min(ceiling(sqrt(c(N, n)))))
    if (is.null(L)) 
        L <- max(10, ceiling(sqrt(n)))
    
    K <- min(K, N); L <- min(n - N, L) 
    if (track) print(paste("Setting K=", K, ", L=", L, sep = ""), quote = FALSE)
    
    A <- array(0, dim = c(n, m, m))
    for (i in 1:n) A[i, , ] <- tcrossprod(Fx[i, ])
    
    finish.all <- FALSE; trMinv.best <- Inf
    while (!finish.all && n.rest < rest.max) {
        n.rest <- n.rest + 1
        if (is.null(w1)) {
            supp <- od_PIN(Fx, echo = FALSE)$supp
            if (N > m) supp <- c(supp, sample(setdiff(1:n, supp), N - m))
            w <- rep(0, n); w[supp] <- 1
            M <- infmat(Fx, w, echo = FALSE)
        } else {
            w <- w1; M <- crossprod(sqrt(w) * Fx); w1 <- NULL
        }
        
        trMinv <- sum(diag(solve(M)))
        if (trMinv.best > trMinv) {
            w.best <- w; M.best <- M; trMinv.best <- trMinv
        }
        
        finish.all <- finish <- as.numeric(proc.time()[3]) > start + t.max
        while (!finish) {
            tm <- as.numeric(proc.time()[3]) - start
            if ((tm > next.sec) && track) {
                Phi.best <- m/trMinv.best
                info <- paste("od_A_KL_bin Time:", round(tm, 1), "Value:", 
                              round(Phi.best, 6), "Efficiency:", 
                              round(Phi.best/Phi.app, 6))
                print(info, quote = FALSE); next.sec <- ceiling(tm)
            }
            a.fun <- ((Fx %*% solve(M))^2) %*% one.m
            supp <- (1:n)[w > 0.5]
            Kind <- supp[order(a.fun[supp])][1:K]
            n.supp <- (1:n)[w < 0.5]
            Lind <- n.supp[order(a.fun[n.supp])][(n - N - L + 1):(n - N)]
            
            imp <- FALSE
            for (iL in L:1) {
                for (iK in 1:K) {
                    M.temp <- M + A[Lind[iL], , ] - A[Kind[iK], , ]
                    if (rcond(M.temp) >  sqrt(.Machine$double.eps)) {
                        trMinv.temp <- sum(diag(solve(M.temp)))
                    } else {
                        trMinv.temp <- Inf
                    }
                    if (trMinv.temp < trMinv) {
                        w[Lind[iL]] <- w[Lind[iL]] + 1
                        w[Kind[iK]] <- w[Kind[iK]] - 1
                        M <- M.temp; trMinv <- trMinv.temp
                        n.ex <- n.ex + 1; imp <- TRUE
                        break
                    }
                }
                if (imp) break
            }
            if (as.numeric(proc.time()[3]) > start + t.max) finish.all <- TRUE
            if (finish.all || !imp) finish <- TRUE
        }
        M <- crossprod(sqrt(w)*Fx)
        trMinv <- sum(diag(solve(M)))
        if (trMinv < trMinv.best) {
            w.best <- w; M.best <- M; trMinv.best <- trMinv
        }
    }
    
    t.act <- round(as.numeric(proc.time()[3]) - start, 2)

    if (track) {
        info <- paste("od_A_KL_bin finished after", t.act, "seconds at", Sys.time())
        print(info, quote = FALSE)
        info <- paste("with", n.rest, "restarts and", n.ex, "exchanges.")
        print(info, quote = FALSE)
    }
    
    Phi.best <- optcrit(Fx, w.best, crit = "A", echo = FALSE)

    return(list(w.best = w.best, Phi.best = Phi.best, eff.best = Phi.best/Phi.app, 
                n.ex = n.ex, n.rest = n.rest, t.act = t.act))
}
