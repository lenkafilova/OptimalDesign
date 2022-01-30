od_D_REX <- function(Fx, w1, ver, gamma, eff, it.max, t.max, track) {
  
  start <- as.numeric(proc.time()[3]); del <- 1e-14; eps <- 1e-24
  Fx <- as.matrix(Fx); n <- nrow(Fx); m <- ncol(Fx)
  
  if (track) {
    info <- paste("Running od_D_REX for cca", t.max, "seconds")
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
  d.fun <- ((Fx %*% t(chol(solve(M))))^2) %*% one / m 
  ord <- order(d.fun, decreasing = TRUE)
  lx.vec <- sample(ord[1:L]); kx.vec <- sample(supp)
  
  while (TRUE) {    
    n.iter <- n.iter + 1; ord1 <- which.min(d.fun[supp])
    kb <- supp[ord1]; lb <- ord[1]; v <- c(kb, lb)
    cv <- Fx[v, ] %*% solve(M, t(Fx[v, ]))
    alpha <- 0.5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + eps)
    alpha <- min(w[kb], alpha)
    w[kb] <- w[kb] - alpha; w[lb] <- w[lb] + alpha
    M <- M + alpha * (tcrossprod(Fx[lb, ]) - tcrossprod(Fx[kb, ]))
    
    if ((w[kb] < del) && (ver == 1)) {
      # LBE is nullifying and the version is 1
      for (l in 1:L) {                          
        lx <- lx.vec[l]; Alx <- tcrossprod(Fx[lx, ])
        for (k in 1:K) {
          kx <- kx.vec[k]; v <- c(kx, lx)
          cv <- Fx[v, ] %*% solve(M, t(Fx[v, ]))
          alpha <- 0.5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + eps)
          alpha <- min(w[kx], max(-w[lx], alpha))
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
          cv <- Fx[v, ] %*% solve(M, t(Fx[v, ]))
          alpha <- 0.5 * (cv[2, 2] - cv[1, 1])/(cv[1, 1] * cv[2, 2] - cv[1, 2]^2 + eps)
          alpha <- min(w[kx], max(-w[lx], alpha))
          w[kx] <- w[kx] - alpha; w[lx] <- w[lx] + alpha
          M <- M + alpha * (Alx - tcrossprod(Fx[kx, ]))                            
        }
      }         
    }
    
    supp <- index[w > del]; K <- length(supp); w.supp <- w[supp]
    d.fun <- ((Fx %*% t(chol(solve(M))))^2) %*% one / m
    ord.ind <- (1:n)[d.fun >= -sort(-d.fun, partial = L)[L]]
    ord <- ord.ind[order(d.fun[ord.ind], decreasing = TRUE)]
    lx.vec <- sample(ord); kx.vec <- sample(supp)    
    tm <- as.numeric(proc.time()[3])
    eff.act <-  1 / d.fun[ord[1]]
    if (track) {
      print(paste("od_D_REX Time:", round(tm - start, 2),
                  "Efficiency:", round(eff.act, 9)), quote = FALSE) 
    }
    if (d.fun[ord[1]] < eff.inv || n.iter >= it.max || tm > start + t.max) break
  }
  
  t.act <- round(as.numeric(proc.time()[3]) - start, 2)
  Phi.best <- det(M)^(1/m); eff.best <- 1/d.fun[ord[1]]
  
  if (track) {
    print(paste("od_D_REX finished at", Sys.time()), quote = FALSE)
    print(paste("Computation time:", t.act), quote = FALSE) 
    print(paste("D-criterion value:", Phi.best), quote = FALSE)
    print(paste("Efficiency at least:", eff.best), quote = FALSE)
  }
  
  return(list(w.best = w, Phi.best = Phi.best, eff.best = eff.best,
              n.iter = n.iter, t.act = t.act))
}
