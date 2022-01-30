od_A_REX <- function(Fx, w1, ver, gamma, eff, it.max, t.max, track) {
  
  stepsize <- function(M, w, v)
  { 
    M.inv <- solve(M)
    dv <- Fx[v, ] %*% M.inv %*% t(Fx[v, ]) 
    av <- Fx[v, ] %*% M.inv %*% M.inv %*% t(Fx[v, ])  
    A <- av[2,2] - av[1,1]; C <- dv[2,2] - dv[1,1]; D <- dv[1,1] * dv[2,2] - dv[1,2]^2  
    B <- 2 * dv[1,2] * av[1,2] - dv[2,2] * av[1,1] - dv[1,1] * av[2,2]  
    G <- A * D + B * C; k <- v[1]; l <- v[2]
    
    if ((abs(G) < del.alpha) && (abs(B) > del.alpha)) {
      r <- -A / (2*B)  
      if ((-w[l] <= r) && (r <= w[k])) return(r)  
    }
    if (abs(G) > 0) {  
      r <- -(B + sqrt(B^2 - A * G)) / G
      if ((abs(r) < del.alpha) && (B < 0)) {
        x <- A * G / B^2
        r <- -B * (x^3/16 + x^2/8 + x/2) / G
      }
      if ((-w[l] <= r) && (r <= w[k])) return(r) 
    }
    
    if (A > del.alpha) {
      return(w[k])
    } else if (A < -del.alpha) {
      return(-w[l])
    } else {
      return(0)
    }
  }
  
  start <- as.numeric(proc.time()[3]); del <- 1e-24; del.alpha <- 1e-14
  Fx <- as.matrix(Fx); n <- nrow(Fx); m <- ncol(Fx)
  
  if (track) {
    info <- paste("Running od_A_REX for cca", t.max, "seconds")
    info <- paste(info, " starting at ", Sys.time(), ".", sep = "")
    print(info, quote = FALSE)
    info <- paste("The problem size is n=", n, sep = "")
    info <- paste(info, " and m=", m, ".", sep = "")
    print(info, quote = FALSE)
  }
  
  eff.inv <- 1/eff; n.iter <- 0; L <- min(n, gamma*m)
  lx.vec <- rep(0, L); index <- 1:n; one <- rep(1, m)
  
  if (is.null(w1)) w1 <- od_PIN(Fx, echo = FALSE)$w.pin/m
  supp <- (1:n)[w1 > 0]; K <- length(supp); Fx.supp <- Fx[supp, ] 
  w <- rep(0, n); w[supp] <- w1[supp]/sum(w1[supp]); w.supp <- w[supp]
  M <- crossprod(sqrt(w.supp) * Fx.supp)
  M.inv <- solve(M); a.fun <- ((Fx %*% M.inv)^2) %*% one
  ord <- order(a.fun, decreasing = TRUE)
  lx.vec <- sample(ord[1:L]); kx.vec <- sample(supp)
  
  while (TRUE) {    
    n.iter <- n.iter + 1; ord1 <- which.min(a.fun[supp])
    kb <- supp[ord1]; lb <- ord[1]; v <- c(kb, lb)
    alpha <- stepsize(M, w, v)
    w[kb] <- w[kb] - alpha; w[lb] <- w[lb] + alpha
    M <- M + alpha * (tcrossprod(Fx[lb, ]) - tcrossprod(Fx[kb, ]))
    
    if ((w[kb] < del) && (ver == 1)) {
      # LBE is nullifying and the version is 1
      for (l in 1:L) {                          
        lx <- lx.vec[l]; Alx <- tcrossprod(Fx[lx, ])
        for (k in 1:K) {
          kx <- kx.vec[k]; v <- c(kx, lx)
          alpha <- stepsize(M, w, v)
          wkx.temp <- w[kx] - alpha; wlx.temp <- w[lx] + alpha
          if ((wkx.temp < del) || (wlx.temp < del)) {
            w[kx] <- wkx.temp; w[lx] <- wlx.temp
            M <- M + alpha * (Alx - tcrossprod(Fx[kx, ]))                            
          }
        }
      }
    } else {
      # LBE is non-nullifying or the version is 0
      for (l in 1:L) {                          
        lx <- lx.vec[l]; Alx <- tcrossprod(Fx[lx, ])
        for (k in 1:K) {
          kx <- kx.vec[k]; v <- c(kx, lx)
          alpha <- stepsize(M, w, v)
          w[kx] <- w[kx] - alpha; w[lx] <- w[lx] + alpha
          M <- M + alpha * (Alx - tcrossprod(Fx[kx, ]))                            
        }
      }         
    }
    
    supp <- index[w > del]; K <- length(supp); w.supp <- w[supp]
    M.inv <- solve(M); a.fun <- ((Fx %*% M.inv)^2) %*% one
    ord.ind <- (1:n)[a.fun >= -sort(-a.fun, partial = L)[L]]
    ord <- ord.ind[order(a.fun[ord.ind], decreasing = TRUE)]
    # The two lines above can be replaced by simpler but usually
    # somewhat slower ord <- order(a.fun, decreasing=TRUE)[1:L]
    lx.vec <- sample(ord); kx.vec <- sample(supp)
    
    tm <- as.numeric(proc.time()[3])
    eff.act <-  sum(diag(M.inv)) / a.fun[ord[1]]
    if (track) {
      print(paste("od_A_REX Time:", round(tm - start, 2),
                  "Efficiency:", round(eff.act, 9)), quote = FALSE) 
    }
    if (1/eff.act < eff.inv || n.iter >= it.max || tm > start + t.max) break
  }
  
  t.act <- round(as.numeric(proc.time()[3]) - start, 2)
  Phi.best <- m/sum(diag(M.inv)); eff.best <- sum(diag(M.inv))/a.fun[ord[1]]
  
  if (track) {
    print(paste("od_A_REX finished at", Sys.time()), quote = FALSE)
    print(paste("Computation time:", t.act), quote = FALSE) 
    print(paste("A-criterion value:", Phi.best), quote = FALSE)
    print(paste("Efficiency at least:", eff.best), quote = FALSE)
  }
  
  return(list(w.best = w, Phi.best = Phi.best, eff.best = eff.best,
              n.iter = n.iter, t.act = t.act))
}
