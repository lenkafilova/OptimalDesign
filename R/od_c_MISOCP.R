od_c_MISOCP <- function(Fx, b1, A1, b2, A2, b3, A3, bin, h, type, gap, t.max) {
  
  if (!requireNamespace('gurobi', quietly = TRUE)) stop("(MI)SOCP requires package gurobi.")
  
  n <- nrow(Fx); m <- ncol(Fx)
  k1 <- length(b1); k2 <- length(b2); k3 <- length(b3)
  b <- c(b1, b2, b3); A <- rbind(A1, A2, A3)
  sense <- c(rep("<=", k1), rep(">=", k2), rep("=", k3))
  
  model <- list()
  da <- dim(A); nc <- da[1]; M <- 5*n
  Aeq <- matrix(0, nrow = m + nc + 2*n, ncol = M)
  Aeq[1:m, (2 * n + 1):(3 * n)] <- 0.5 * t(Fx)
  Aeq[(m + 1):(m + nc), 1:n] <- A
  Aeq[(m + nc + 1):(m + nc + n), 1:n] <- diag(n)
  Aeq[(m + nc + 1):(m + nc + n), (n + 1):(2*n)] <- diag(n)
  Aeq[(m + nc + 1):(m + nc + n), (3*n + 1):(4*n)] <- -diag(n)
  Aeq[(m + nc + 1 + n):(m + nc + 2*n), 1:n] <- diag(n)
  Aeq[(m + nc + 1 + n):(m + nc + 2*n), (n + 1):(2*n)] <- -diag(n)
  Aeq[(m + nc + 1 + n):(m + nc + 2*n), (4*n + 1):(5*n)] <- -diag(n)
  model$A <- Matrix::Matrix(Aeq, sparse = TRUE)
  beq <- c(0.5*h, b, rep(0, 2*n))
  model$rhs <- beq
  
  sns <- rep("=", m + nc + 2*n)
  for (i in (m + 1):(m + nc)) sns[i] <- sense[i - m]
  model$sense <- sns      
  
  model$quadcon <- vector("list", n)
  for (i in 1:n) {
    qc <- list()
    qc$Qc <- Matrix::spMatrix(M, M, c(3*n + i, 4*n + i, 2*n + i), 
                               c(3*n + i, 4*n + i, 2*n + i), c(-1, 1, 1))
    qc$rhs <- 0
    model$quadcon[[i]] <- qc
  }
  
  lb <- rep(-Inf, M); lb[1:(2*n)] <- 0; lb[(3*n + 1):(4*n)] <- 0
  model$lb <- lb; model$ub <- rep(Inf, M)
  cc <- rep(0, M); cc[(n + 1):(2*n)] <- 1
  model$obj <- cc
  
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
  res <- as.vector(result$x)
  
  if (!is.vector(res) || !is.numeric(res) || !all(is.finite(res)) || length(res) < n) {
    w.best <- NULL
  } else {
    w.best <- res[1:n]
    if (type == "exact") w.best <- round(res[1:n])
  }
  
  return(list(w.best = w.best, status = result$status))
}

