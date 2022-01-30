od_D_RC <- function(Fx, b, A, w0, Phi.app, w1, rest.max, t.max, track) {
    
    start <- as.numeric(proc.time()[3])
    if (track) {
        info <- paste("Running od_D_RC for cca", t.max, "seconds")
        info <- paste(info, " starting at ", Sys.time(), ".", sep = "")
        print(info, quote = FALSE)
    }
    eps <- sqrt(.Machine$double.eps)
    back.max <- 20; V.max <- 1e+08
    V <- rep(FALSE, V.max + 1)
    n <- nrow(Fx); m <- ncol(Fx)
    if (track) {
        info <- paste("The problem size is n=", n, sep = "")
        info <- paste(info, " and m=", m, ".", sep = "")
        print(info, quote = FALSE)
    }
    next.sec <- 0; n.rest <- 1
    one.n <- rep(1, n); E <- eps * diag(m)
    b <- b + eps; A <- A + 1e-24
    Phi <- function(w) {
        supp <- as.logical(w)
        w.tmp <- w[supp]; Fx.tmp <- Fx[supp, ]
        if (sum(supp) == 1) Fx.tmp <- t(Fx.tmp)
        det(crossprod(sqrt(w.tmp) * Fx.tmp) + E)^(1/m)
    }
    Phi.exact <- function(w) {
        supp <- as.logical(w)
        w.tmp <- w[supp]; Fx.tmp <- Fx[supp, ]
        if (length(supp) == 1) Fx.tmp <- t(Fx.tmp)
        M <- crossprod(sqrt(w.tmp) * Fx.tmp)
        if (rcond(M) > 1e-12) return(max(c(0, det(M)))^(1/m))
        return(0)
    }
    compute.gamma <- function(numer, denom) {
        if (!any(denom > eps)) return(0)
        return(min((numer/denom)[denom > eps]))
    }
    explore.up <- function() {
        wei <- w; val <- -one.n
        for (i in which(d > 0)) {
            wei[i] <- w[i] + 1
            attrei <- floor(V.max * (sin(Phi(wei)) + 1)/2 + 1)
            if (!V[attrei]) {
                dei <- floor(matrixStats::colMins(as.vector(res - A[, i])/A) + eps)
                gammaei <- compute.gamma(res - A[, i], A %*% dei)
                val[i] <- Phi(wei + gammaei * dei)
            }
            wei[i] <- w[i]
        }
        up.index <- which.max(val)
        max.val <- val[up.index]
        if (max.val < 0) up.index <- 0
        up.index
    }
    explore.down <- function() {
        wei <- w; val <- -one.n
        for (i in which(w > w0)) {
            wei[i] <- w[i] - 1
            attrei <- floor(V.max * (sin(Phi(wei)) + 1)/2 + 1)
            if (!V[attrei]) {
                dei <- floor(matrixStats::colMins(as.vector(res + A[, i])/A) + eps)
                gammaei <- compute.gamma(res + A[, i], A %*% dei)
                val[i] <- Phi(wei + gammaei * dei)
            }
            wei[i] <- w[i]
        }
        down.index <- which.max(val)
        max.val <- val[down.index]
        if (max.val < 0) down.index <- 0
        down.index
    }
    random.step <- function() {
        success <- FALSE
        while (!success) {
            dir <- sample(c(FALSE, TRUE), 1)
            if (dir) {
                up.index <- sample(n, 1)
                if (d[up.index] > 0) {
                    w[up.index] <<- w[up.index] + 1
                    res <<- res - A[, up.index]
                    d <<- floor(matrixStats::colMins(as.vector(res)/A) + eps)
                    attrb <<- floor(V.max * (sin(Phi(w)) + 1)/2 + 1)
                    success <- TRUE
                }
            } else {
                down.index <- sample(n, 1)
                if (w[down.index] > w0[down.index]) {
                    w[down.index] <<- w[down.index] - 1
                    res <<- res + A[, down.index]
                    d <<- floor(matrixStats::colMins(as.vector(res)/A) + eps)
                    attrb <<- floor(V.max * (sin(Phi(w)) + 1)/2 + 1)
                    success <- TRUE
                    back.no <<- back.no + 1
                }
            }
        }
    }
    finish <- FALSE; back.no <- 0
    start <- as.numeric(proc.time()[3])
    if (is.null(w0)) w0 <- rep(0, n)
    if (is.null(w1) || sum(w1) == 0) {
        w <- w0; res <- b - A %*% w
        d <- floor(matrixStats::colMins(as.vector(res)/A) + eps)
        maximal <- FALSE
        while (!maximal) {
            val <- (d > 0)
            if (sum(val) == 0) {
                maximal <- TRUE
            }
            else {
                i <- sample(1:n, 1, prob = d)
                w[i] <- w[i] + 1
                res <- b - A %*% w
                d <- floor(matrixStats::colMins(as.vector(res)/A) + eps)
            }
        }
    } else {
        w <- w1; res <- b - A %*% w
        d <- floor(matrixStats::colMins(as.vector(res)/A) + eps)
    }
    attrb <- floor(V.max * (sin(Phi(w)) + 1)/2 + 1)
    w.best <- w; res.best <- res
    d.best <- d; attr.best <- attrb
    Phi.best <- Phi.exact(w)
    while (!finish) {
        tm <- as.numeric(proc.time()[3]) - start
        if ((tm > next.sec) && track) {
            if (!is.null(Phi.app)) {
            info <- paste("od_D_RC Time:", round(tm, 1), "Value / Efficiency:", 
                          round(Phi.exact(w.best), 6), "/",
                          round(Phi.exact(w.best)/Phi.app, 6))
            } else {
                info <- paste("od_D_RC Time:", round(tm, 1), "Value:", 
                              round(Phi.exact(w.best), 6))
            }
            print(info, quote = FALSE)
            next.sec <- ceiling(tm)
        }
        if (!V[attrb]) {
            V[attrb] <- TRUE
            up.index <- explore.up()
            if (up.index) {
                w[up.index] <- w[up.index] + 1
                res <- res - A[, up.index]
                d <- floor(matrixStats::colMins(as.vector(res)/A) + eps)
                attrb <- floor(V.max * (sin(Phi(w)) + 1)/2 + 1)
            } else {
                if (!any(as.logical(d))) {
                    Phi.w.cmp <- Phi.exact(w)
                    if (Phi.w.cmp > Phi.best) {
                        w.best <- w; res.best <- res
                        d.best <- d; attr.best <- attrb
                        Phi.best <- Phi.w.cmp
                        back.no <- 0
                    }
                }
                down.index <- explore.down()
                if (down.index) {
                    w[down.index] <- w[down.index] - 1
                    res <- res + A[, down.index]
                    d <- floor(matrixStats::colMins(as.vector(res)/A) + eps)
                    attrb <- floor(V.max * (sin(Phi(w)) + 1)/2 + 1)
                    back.no <- back.no + 1
                } else {
                    random.step()
                }
            }
        } else {
            down.index <- explore.down()
            if (down.index) {
                w[down.index] <- w[down.index] - 1
                res <- res + A[, down.index]
                d <- floor(matrixStats::colMins(as.vector(res)/A) + eps)
                attrb <- floor(V.max * (sin(Phi(w)) + 1)/2 + 1)
                back.no <- back.no + 1
            }
            else {
                up.index <- explore.up()
                if (up.index) {
                    w[up.index] <- w[up.index] + 1
                    res <- res - A[, up.index]
                    d <- floor(matrixStats::colMins(as.vector(res)/A) + eps)
                    attrb <- floor(V.max * (sin(Phi(w)) + 1)/2 + 1)
                }
                else {
                    random.step()
                }
            }
        }
        if (back.no > back.max) {
            w <- w.best; res <- res.best
            d <- d.best; attrb <- attr.best
            back.no <- 0; n.rest <- n.rest + 1
        }
        ptm <- as.numeric(proc.time()[3]) - start
        if ((ptm > t.max) || (n.rest > rest.max)) finish <- TRUE
    }
    t.act <- round(as.numeric(proc.time()[3]) - start, 2)
    if (track) {
        info <- paste("od_D_RC finished after", t.act, "seconds at", Sys.time())
        print(info, quote = FALSE)
    }
    
    return(list(w.best = w.best, Phi.best = Phi.best, eff.best = Phi.best/Phi.app, 
                n.rest = n.rest, t.act = t.act))
}
