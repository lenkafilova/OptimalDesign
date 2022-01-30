od_SYM <- function(Fx, w, b1=NULL, A1=NULL, b2=NULL, A2=NULL, b3=NULL, A3=NULL,
                   w0=NULL, crit="D", h=NULL, echo=TRUE) {
  # "Symmetrize" a design without changing its information matrix (much)
  # It is assumed that either the size of the design space is less than 1000
  # or w is almost-optimal with support of size < 1000 (e.g., produced by REX)
  # Note: The design w itself shoud satisfy the constraints
  
  cl <- match.call()
  verify(cl, Fx = Fx, w = w, b1 = b1, A1 = A1, b2 = b2, A2 = A2, b3 = b3, A3 = A3,
         w0 = w0, crit = crit, h = h, echo = echo)
  n <- nrow(Fx); m <- ncol(Fx)
  k1 <- k2 <- k3 <- 0
  eps <- sqrt(.Machine$double.eps)
  
  if (is.null(b1) & is.null(b2) & is.null(b3)) {
    message("b1, b2 and b3 are NULL. Setting b1=sum(w).")
    b1 <- sum(w);  k1 <- 1
  }
  if (!is.null(b1)) {
    k1 <- length(b1)
    if (is.null(A1)) A1 <- matrix(1, nrow = k1, ncol = n)
  } else {
    b1 <- 0; k1 <- 1
    A1 <- matrix(0, nrow = 1, ncol = n)
  }
  if (!is.null(b2)) {
    k2 <- length(b2)
    if (is.null(A2)) A2 <- matrix(1, nrow = k2, ncol = n)
  }else {
    b2 <- 0; k2 <- 1
    A2 <- matrix(0, nrow = 1, ncol = n)
  }
  if (!is.null(b3)) {
    k3 <- length(b3)
    if (is.null(A3)) A3 <- matrix(1, nrow = k3, ncol = n)
  }
  else {
    b3 <- 0; k3 <- 1
    A3 <- matrix(0, nrow = 1, ncol = n)
  }
  if (!is.null(w0) && sum(w0) > 0) {
    A2 <- rbind(A2, diag(n)[w0 > 0,])
    b2 <- c(b2, w0[w0 > 0])
  }
  
  if (crit %in% c("C", "c") && is.null(h)) h <- c(rep(0, m - 1), 1)  
  
  supp <- (1:n)[w > eps]
  
  if (length(supp) <= 1000) {
    
    M <- infmat(Fx, w, echo = FALSE)
    k <- round(0.5*m*(m + 1))
    b.M <- rep(0, k)
    
    if (n <= 1000) {
      A.M <- matrix(0, nrow = k, ncol = n)
      temp <- 1
      for (i in 1:m) {
        for (j in 1:i) {
          A.M[temp, ] <- c(Fx[, i] * Fx[, j])
          b.M[temp] <- M[i, j]; temp <- temp + 1
        }
      }
      A <- rbind(A.M, -A.M, A3, -A3, -A1, A2, diag(n))
      b <- c(b.M, -b.M, b3, -b3, -b1, b2, rep(0, n))
      print(paste("Computing a QP problem, please wait.", Sys.time()))
      res <- quadprog::solve.QP(diag(n), rep(0, n), t(A), b - eps)
      w.sym <- res$solution[1:n]
      return(list(call = cl, w.sym = pmax(w.sym, 0))) 
    }
    
    if (rcond(M) > sqrt(.Machine$double.eps)) {
      var.fun <- varfun(Fx, w, crit, h, echo = FALSE)
      var.fun[supp] <- Inf; ord <- order(var.fun, decreasing = TRUE)
      supp.ex <- (1:n)[ord[1:1000]]; ns <- 1000
      A.M <- matrix(0, nrow = k, ncol = ns)
      Fs <- Fx[supp.ex, ]; temp <- 1
      for (i in 1:m) {
        for (j in 1:i) {
          A.M[temp, ] <- c(Fs[, i] * Fs[, j])
          b.M[temp] <- M[i, j]; temp <- temp + 1
        }
      }
      
      A <- rbind(A.M, -A.M, A3[, supp.ex], -A3[, supp.ex], -A1[, supp.ex], 
                 A2[, supp.ex], diag(ns))
      b <- c(b.M, -b.M, b3, -b3, -b1, b2, rep(0, ns))
      print(paste("Computing a QP problem, please wait.", Sys.time()))
      res <-  quadprog::solve.QP(diag(ns), rep(0, ns), t(A), b - eps)
      w.sym <- rep(0, n); w.sym[supp.ex] <- res$solution[1:ns]
      return(list(call = cl, w.sym = pmax(w.sym, 0)))
    }
  }
  
  warning("The design to be symmetrized is close to singular or its support is too large.")
  return(list(call = cl, w.sym = w))
}
