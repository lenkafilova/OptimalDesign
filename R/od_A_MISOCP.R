od_A_MISOCP <- function(Fx, b1, A1, b2, A2, b3, A3, bin, type, gap, t.max) {

  if (!requireNamespace('gurobi', quietly = TRUE)) stop("(MI)SOCP requires package gurobi.")

  n <- nrow(Fx); m <- ncol(Fx)
  k1 <- length(b1); k2 <- length(b2); k3 <- length(b3)
  b <- c(b1, b2, b3); A <- rbind(A1, A2, A3)
  sense <- c(rep("<=", k1), rep(">=", k2), rep("=", k3))
  model  <-  list(); model$modelsense <- "min"
  da <- dim(A); nc <- da[1]; M <- 4*n + m*n
  
  Aeq <- matrix(0, nrow = m*m + nc + 2*n, ncol = M) 
  A2 <- matrix(0, nrow = nc, ncol = M)
  A2[, 1:n] <- A 

  k <- 1
  for (i in 1:(m*m)) {
    for (j in 1:n) {
      u <- i %% m; if (u == 0) u <- m
      Aeq[i, (3 + k)*n + j] = 0.5*Fx[j, u]
    }
    if (i %% m == 0) k <- k + 1
  }
  
  A3 <- matrix(0, nrow = 2*n, ncol = M)
  for (i in 1:n) {
    A3[i, i] <- 1; A3[i, n + i] <- 1
    A3[i, 2*n + i] <- -1; A3[n + i, i] <- 1
    A3[n + i, n + i] <- -1; A3[n + i, 3*n + i] <- 1
  }
  
  beq <- rep(0, m*m + nc + 2*n)
  beq[1:(m^2)] <- as.vector(diag(m))
  
  for (i in 1:nc) {
    Aeq[m*m + i, ] <- A2[i, ]
    beq[m*m + i] <- b[i]
  }
  
  for (i in 1:(2*n)) 
    Aeq[m*m + nc + i, ] <- A3[i, ]
  
  model$A <- Matrix::Matrix(Aeq, sparse = TRUE)
  model$rhs <- beq
  
  sns <- rep("=", m*m + nc + 2*n)
  for (i in (m*m + 1):(m*m + nc)) sns[i] <- sense[i - m^2]
  model$sense <- sns
  
  model$quadcon <- vector("list", n)
  for (i in 1:n) {
    a <- 4:(3 + m)
    qc <- list()
    u <- length(c(2*n + i, a*n + i, 3*n + i))
    qc$Qc <- Matrix::spMatrix(M, M, c(2*n + i, a*n + i, 3*n + i),
                              c(2*n + i, a*n + i, 3*n + i), c(-1, rep(1, u - 1)))
    qc$rhs <- 0
    model$quadcon[[i]] <- qc
  }
  
  lb <- rep(-Inf, M)
  lb[1:(3 * n)] <- 0
  model$lb <- lb
  model$ub <- rep(Inf, M)
  
  c <- rep(0, M)
  c[(n + 1):(2*n)] <- 1
  model$obj <- c
  
  vtypes <- rep("C", M)
  if (type == "exact") {
    vtypes[1:n] <- "I"
    if (bin) vtypes[1:n] <- "B"
  }
  model$vtypes <- vtypes
  
  if (type == "approximate" && bin)
    model$ub[1:n] <- rep(1, n) 
  
  if (!is.null(gap)) {
    params <- list(TimeLimit = t.max, MIPFocus = 1, MIPGap = gap)
  } else {
    params <- list(TimeLimit = t.max, MIPFocus = 1)
  }
  
  result <- gurobi::gurobi(model, params); gc()
  res <- result$x
  
  if (!is.vector(res) || !is.numeric(res) || !all(is.finite(res)) || length(res) < n) {
    w.best <- NULL
  } else {
    w.best <- res[1:n]
    if (type == "exact") w.best <- round(res[1:n])
  }
  
  return(list(w.best = w.best, status = result$status))
}
