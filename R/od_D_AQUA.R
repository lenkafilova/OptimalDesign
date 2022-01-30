od_D_AQUA <- function(Fx, b1, A1, b2, A2, b3, A3, bin, M.anchor, ver, conic, t.max) {

  if (!requireNamespace('gurobi', quietly = TRUE)) stop("AQuA requires the gurobi package.")
  
  n <- nrow(Fx); m <- ncol(Fx) 
  k1 <- length(b1); k2 <- length(b2); k3 <- length(b3)
  b <- c(b1, b2, b3); A <- rbind(A1, A2, A3)
  sense <- c(rep("<", k1), rep(">", k2), rep("=", k3))
  
  if (is.null(M.anchor)) {
    if (length(b) == 1 && all(A == matrix(1, ncol = n, nrow = 1))) { 
      M.anchor <- infmat(Fx, b*od_REX(Fx, t.max = t.max/5,
                                      track = FALSE)$w.best, echo = FALSE)
    } else {
      w.app <- od_MISOCP(Fx, b1, A1, b2, A2, b3, A3, bin = bin, t.max = t.max/3,
                         type = "approximate", gap = NULL, echo = FALSE)$w.best
      if (is.null(w.app)) {
        stop("Failure in the computation of M.anchor (and no M.anchor supplied by the user).")
      } else {
        M.anchor <- infmat(Fx, w.app, echo = FALSE)
      }
    }
  }
  
  M1 <- solve(M.anchor)
  
  if (conic) {
    
    s <- m*(m + 1)/2; Gm <- matrixcalc::duplication.matrix(m); eps <- 1e-12
    vec.M1 <- matrixcalc::vec(M1)
    tilde.h <- t(Gm) %*% vec.M1
    
    if (ver == "+") tilde.Q <- 0.5 * t(Gm) %*% (-tcrossprod(vec.M1) / m + M1 %x% M1) %*% Gm
    if (ver == "-") tilde.Q <- (1/6) * t(Gm) %*% (tcrossprod(vec.M1) / m + M1 %x% M1) %*% Gm
    
    if (rcond(tilde.Q) > eps) {
      tilde.C <- t(chol(tilde.Q))
    } else {  
      eig <- eigen(tilde.Q, symmetric = TRUE)
      mx <- max(eig$values); ts <- sum(eig$values/mx > eps)
      tilde.C <- eig$vectors[, 1:ts] %*% diag(sqrt(eig$values[1:ts]))
    }
    diff <- max(abs(tilde.Q - tcrossprod(tilde.C)))
    print(paste("Decomposition instability:", diff))
    
    H <- matrix(0, nrow = n, ncol = s); k <- 0
    for (i1 in 1:m) {
      for (i2 in i1:m) {
        k <- k + 1; H[, k] <- Fx[, i1]*Fx[, i2]
      }
    }
    h <- H %*% tilde.h; S <- H %*% tilde.C; u <- ncol(S)
    
    model <- list(); k <- nrow(A); np <- n + u + 3
    
    Aeq <- matrix(0, nrow = k + 2 + u, ncol = np); Aeq[1:k, 1:n] <- A
    Aeq[k + 1, n + 1] <- 1; Aeq[k + 1, n + 3] <- -1/sqrt(2)
    Aeq[k + 2, n + 2] <- 1; Aeq[k + 2, n + 3] <- 1/sqrt(2)
    Aeq[(k + 3):(k + 2 + u), 1:n] <- t(S)
    Aeq[(k + 3):(k + 2 + u), (n + 4):np] <- -diag(u) 
    model$A <- Matrix::Matrix(Aeq, sparse = TRUE)
    
    beq <- rep(0, k + 2 + u); beq[1:k] <- b
    beq[k + 1] <- 1/2/sqrt(2); beq[k + 2] <- 1/2/sqrt(2)
    model$rhs <- beq
    
    model$sense <- c(rep("<", k1), rep(">", k2), rep("=", k3 + u + 2))
    
    qc1 <- list()
    qc1$Qc <- Matrix::sparseMatrix(i = c(n + 1, n + 2, (n + 4):np),
                                   j = c(n + 1, n + 2, (n + 4):np),
                                   x = c(-1, rep(1, np - n - 2)),
                                   dims = c(np, np), symmetric = TRUE)
    qc1$rhs <- 0; qc1$q <- rep(0, np)
    model$quadcon <- list(qc1)
    
    lb <- rep(-Inf, np); lb[1:(n + 1)] <- 0; model$lb <- lb
    ob <- rep(0, np); ob[1:n] <- -h; ob[n + 3] <- 1
    
    model$obj <- ob; vtypes <- rep("C", np)
    if (bin) {vtypes[1:n] <- "B"} else {vtypes[1:n] <- "I"}
    model$vtypes <- vtypes; model$modelsense <- "min"
    
    params  <-  list(TimeLimit = t.max, MIPFocus = 1)
    result <- gurobi::gurobi(model, params); gc()
    
  } else {
 
    V1 <- Fx %*% M1 %*% t(Fx)
    h <- -diag(V1)
    
    if (ver == "+") Q <- 0.5*V1^2 - 0.5/m*h %*% t(h)
    if (ver == "-") Q <- 1/6*V1^2 + (1/6/m)*h %*% t(h)
    Q <- 0.5*(Q + t(Q)) # Numerical errors may make the original Q slighly non-symmetric
       
    model <- list(); model$obj <- h
    model$Q <- Q; model$A <- Matrix::Matrix(A, sparse = TRUE); model$rhs <- c(b)
    model$sense <- sense; model$lb <- rep(0, n)
    if (bin) {model$vtypes <- "B"} else {model$vtypes <- "I"}
    params <- list(TimeLimit = t.max, MIPFocus = 1)
    result <- gurobi::gurobi(model, params); gc()
  }
  
  res <- result$x
  if (!is.vector(res) || !is.numeric(res) || !all(is.finite(res)) || length(res) < n) {
    w.best <- NULL
  } else {
    w.best <- round(res[1:n])
  }    
  
  return(list(w.best = w.best, status = result$status))
}

