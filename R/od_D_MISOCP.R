od_D_MISOCP <- function(Fx, b1, A1, b2, A2, b3, A3, bin, type, gap, t.max) {

  if (!requireNamespace('gurobi', quietly = TRUE)) stop("(MI)SOCP requires package gurobi.")

  geomeanG <- function(m, n, M) {
    
    cones.index <- NULL
    
    if (m == 2) {
      cones.index <- vector("list", 1)
      cones.index[[1]] <- c(M + 1, M + 2, M + 3)
      newvars <- 3; ncg <- 2; cn <- 1
      Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
      Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
      Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
    }
    
    if (m == 3) {
      cones.index <- vector("list", 3)
      cones.index[[1]] <- c(M + 1, M + 2, M + 3)
      cones.index[[2]] <- c(M + 4, M + 5, M + 6)
      cones.index[[3]] <- c(M + 7, M + 8, M + 9)
      ncg <- 6; newvars <- 9; cn <- 3
      Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
      Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
      Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
      Ag[3, n + 2 * m + 3] <- 1;  Ag[3, M + 9] <- 1/2;  Ag[3, M + 4] <- -1
      Ag[4, n + 2 * m + 3] <- 1;  Ag[4, M + 9] <- -1/2;  Ag[4, M + 5] <- -1
      Ag[5, M + 3] <- 1/2;  Ag[5, M + 6] <- 1/2;  Ag[5, M + 7] <- -1
      Ag[6, M + 3] <- 1/2;  Ag[6, M + 6] <- -1/2;  Ag[6, M + 8] <- -1
    }
    
    if (m == 4) {
      cones.index <- vector("list", 3)
      cones.index[[1]] <- c(M + 1, M + 2, M + 3)
      cones.index[[2]] <- c(M + 4, M + 5, M + 6)
      cones.index[[3]] <- c(M + 7, M + 8, M + 9)
      ncg <- 6; newvars <- 9; cn <- 3
      Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
      Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
      Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
      Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
      Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
      Ag[5, M + 3] <- 1/2;  Ag[5, M + 6] <- 1/2;  Ag[5, M + 7] <- -1
      Ag[6, M + 3] <- 1/2;  Ag[6, M + 6] <- -1/2;  Ag[6, M + 8] <- -1
    }
    
    if (m == 5) {
      cones.index <- vector("list", 6)
      cones.index[[1]] <- c(M + 1, M + 2, M + 3)
      cones.index[[2]] <- c(M + 4, M + 5, M + 6)
      cones.index[[3]] <- c(M + 7, M + 8, M + 9)
      cones.index[[4]] <- c(M + 10, M + 11, M + 12)
      cones.index[[5]] <- c(M + 13, M + 14, M + 15)
      cones.index[[6]] <- c(M + 16, M + 17, M + 18)
      cn <- 6; newvars <- 18; ncg <- 12
      Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
      Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
      Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
      Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
      Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
      Ag[5, n + 4 * m + 5] <- 1;  Ag[5, M + 18] <- 1;  Ag[5, M + 7] <- -1 
      Ag[6, n + 4 * m + 5] <- 1;  Ag[6, M + 18] <- -1;  Ag[6, M + 8] <- -1
      Ag[7, M + 3] <- 1/2;  Ag[7, M + 6] <- 1/2;  Ag[7, M + 10] <- -1
      Ag[8, M + 3] <- 1/2;  Ag[8, M + 6] <- -1/2;  Ag[8, M + 11] <- -1
      Ag[9, M + 9] <- 1/2;  Ag[9, M + 18] <- 1/2;  Ag[9, M + 13] <- -1
      Ag[10, M + 9] <- 1/2;  Ag[10, M + 18] <- -1/2;  Ag[10, M + 14] <- -1
      Ag[11, M + 12] <- 1/2;  Ag[11, M + 15] <- 1/2;  Ag[11, M + 16] <- -1
      Ag[12, M + 12] <- 1/2;  Ag[12, M + 15] <- -1/2;  Ag[12, M + 17] <- -1
    }
    
    if (m == 6) {
      cones.index <- vector("list", 6)
      cones.index[[1]] <- c(M + 1, M + 2, M + 3)
      cones.index[[2]] <- c(M + 4, M + 5, M + 6)
      cones.index[[3]] <- c(M + 7, M + 8, M + 9)
      cones.index[[4]] <- c(M + 10, M + 11, M + 12)
      cones.index[[5]] <- c(M + 13, M + 14, M + 15)
      cones.index[[6]] <- c(M + 16, M + 17, M + 18)
      cn <- 6; newvars <- 18; ncg <- 12
      Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
      Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
      Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
      Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
      Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
      Ag[5, n + 4 * m + 5] <- 1;  Ag[5, n + 5 * m + 6] <- 1;  Ag[5, M + 7] <- -1 
      Ag[6, n + 4 * m + 5] <- 1;  Ag[6, n + 5 * m + 6] <- -1;  Ag[6, M + 8] <- -1
      Ag[7, M + 3] <- 1/2;  Ag[7, M + 6] <- 1/2;  Ag[7, M + 10] <- -1
      Ag[8, M + 3] <- 1/2;  Ag[8, M + 6] <- -1/2;  Ag[8, M + 11] <- -1
      Ag[9, M + 9] <- 1/2;  Ag[9, M + 18] <- 1/2;  Ag[9, M + 13] <- -1
      Ag[10, M + 9] <- 1/2;  Ag[10, M + 18] <- -1/2;  Ag[10, M + 14] <- -1
      Ag[11, M + 12] <- 1/2;  Ag[11, M + 15] <- 1/2;  Ag[11, M + 16] <- -1
      Ag[12, M + 12] <- 1/2;  Ag[12, M + 15] <- -1/2;  Ag[12, M + 17] <- -1
    }
    
    if (m == 7) {
      cones.index <- vector("list", 7)
      cones.index[[1]] <- c(M + 1, M + 2, M + 3)
      cones.index[[2]] <- c(M + 4, M + 5, M + 6)
      cones.index[[3]] <- c(M + 7, M + 8, M + 9)
      cones.index[[4]] <- c(M + 10, M + 11, M + 12)
      cones.index[[5]] <- c(M + 13, M + 14, M + 15)
      cones.index[[6]] <- c(M + 16, M + 17, M + 18)
      cones.index[[7]] <- c(M + 19, M + 20, M + 21)
      cn <- 7; newvars <- 21; ncg <- 14
      Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
      Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
      Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
      Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
      Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
      Ag[5, n + 4 * m + 5] <- 1;  Ag[5, n + 5 * m + 6] <- 1;  Ag[5, M + 7] <- -1 
      Ag[6, n + 4 * m + 5] <- 1;  Ag[6, n + 5 * m + 6] <- -1;  Ag[6, M + 8] <- -1
      Ag[7, n + 6 * m + 7] <- 1;  Ag[7, M + 21] <- 1;  Ag[7, M + 10] <- -1
      Ag[8, n + 6 * m + 7] <- 1;  Ag[8, M + 21] <- -1;  Ag[8, M + 11] <- -1
      Ag[9, M + 3] <- 1/2;  Ag[9, M + 6] <- 1/2;  Ag[9, M + 13] <- -1
      Ag[10, M + 3] <- 1/2;  Ag[10, M + 6] <- -1/2;  Ag[10, M + 14] <- -1
      Ag[11, M + 9] <- 1/2;  Ag[11, M + 12] <- 1/2;  Ag[11, M + 16] <- -1
      Ag[12, M + 9] <- 1/2;  Ag[12, M + 12] <- -1/2;  Ag[12, M + 17] <- -1
      Ag[13, M + 15] <- 1/2;  Ag[13, M + 18] <- 1/2;  Ag[13, M + 19] <- -1
      Ag[14, M + 15] <- 1/2;  Ag[14, M + 18] <- -1/2;  Ag[14, M + 20] <- -1
    }
    
    if (m == 8)  {
      cones.index <- vector("list", 7)
      cones.index[[1]] <- c(M + 1, M + 2, M + 3)
      cones.index[[2]] <- c(M + 4, M + 5, M + 6)
      cones.index[[3]] <- c(M + 7, M + 8, M + 9)
      cones.index[[4]] <- c(M + 10, M + 11, M + 12)
      cones.index[[5]] <- c(M + 13, M + 14, M + 15)
      cones.index[[6]] <- c(M + 16, M + 17, M + 18)
      cones.index[[7]] <- c(M + 19, M + 20, M + 21)
      cn <- 7; newvars <- 21; ncg <- 14
      Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
      Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
      Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
      Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
      Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
      Ag[5, n + 4 * m + 5] <- 1;  Ag[5, n + 5 * m + 6] <- 1;  Ag[5, M + 7] <- -1 
      Ag[6, n + 4 * m + 5] <- 1;  Ag[6, n + 5 * m + 6] <- -1;  Ag[6, M + 8] <- -1
      Ag[7, n + 6 * m + 7] <- 1;  Ag[7, n + 7 * m + 8] <- 1;  Ag[7, M + 10] <- -1
      Ag[8, n + 6 * m + 7] <- 1;  Ag[8, n + 7 * m + 8] <- -1;  Ag[8, M + 11] <- -1
      Ag[9, M + 3] <- 1/2;  Ag[9, M + 6] <- 1/2;  Ag[9, M + 13] <- -1
      Ag[10, M + 3] <- 1/2;  Ag[10, M + 6] <- -1/2;  Ag[10, M + 14] <- -1
      Ag[11, M + 9] <- 1/2;  Ag[11, M + 12] <- 1/2;  Ag[11, M + 16] <- -1
      Ag[12, M + 9] <- 1/2;  Ag[12, M + 12] <- -1/2;  Ag[12, M + 17] <- -1
      Ag[13, M + 15] <- 1/2;  Ag[13, M + 18] <- 1/2;  Ag[13, M + 19] <- -1
      Ag[14, M + 15] <- 1/2;  Ag[14, M + 18] <- -1/2;  Ag[14, M + 20] <- -1
    }
    
    if (m == 9) {
      cones.index <- vector("list", 11)
      cones.index[[1]] <- c(M + 1, M + 2, M + 3)
      cones.index[[2]] <- c(M + 4, M + 5, M + 6)
      cones.index[[3]] <- c(M + 7, M + 8, M + 9)
      cones.index[[4]] <- c(M + 10, M + 11, M + 12)
      cones.index[[5]] <- c(M + 13, M + 14, M + 15)
      cones.index[[6]] <- c(M + 16, M + 17, M + 18)
      cones.index[[7]] <- c(M + 19, M + 20, M + 21)
      cones.index[[8]] <- c(M + 22, M + 23, M + 24)
      cones.index[[9]] <- c(M + 25, M + 26, M + 27)
      cones.index[[10]] <- c(M + 28, M + 29, M + 30)
      cones.index[[11]] <- c(M + 31, M + 32, M + 33)
      cn <- 11; newvars <- 33; ncg <- 22
      Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
      Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
      Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
      Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
      Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
      Ag[5, n + 4 * m + 5] <- 1;  Ag[5, n + 5 * m + 6] <- 1;  Ag[5, M + 7] <- -1 
      Ag[6, n + 4 * m + 5] <- 1;  Ag[6, n + 5 * m + 6] <- -1;  Ag[6, M + 8] <- -1
      Ag[7, n + 6 * m + 7] <- 1;  Ag[7, n + 7 * m + 8] <- 1;  Ag[7, M + 10] <- -1
      Ag[8, n + 6 * m + 7] <- 1;  Ag[8, n + 7 * m + 8] <- -1;  Ag[8, M + 11] <- -1
      Ag[9, n + 8 * m + 9] <- 1;  Ag[9, M + 33] <- 1;  Ag[9, M + 13] <- -1
      Ag[10, n + 8 * m + 9] <- 1;  Ag[10, M + 33] <- -1;  Ag[10, M + 14] <- -1
      Ag[11, M + 3] <- 1/2;  Ag[11, M + 6] <- 1/2;  Ag[11, M + 16] <- -1
      Ag[12, M + 3] <- 1/2;  Ag[12, M + 6] <- -1/2;  Ag[12, M + 17] <- -1
      Ag[13, M + 9] <- 1/2;  Ag[13, M + 12] <- 1/2;  Ag[13, M + 19] <- -1
      Ag[14, M + 9] <- 1/2;  Ag[14, M + 12] <- -1/2;  Ag[14, M + 20] <- -1
      Ag[15, M + 15] <- 1/2;  Ag[15, M + 33] <- 1/2;  Ag[15, M + 22] <- -1
      Ag[16, M + 15] <- 1/2;  Ag[16, M + 33] <- -1/2;  Ag[16, M + 23] <- -1
      Ag[17, M + 18] <- 1/2;  Ag[17, M + 21] <- 1/2;  Ag[17, M + 25] <- -1
      Ag[18, M + 18] <- 1/2;  Ag[18, M + 21] <- -1/2;  Ag[18, M + 26] <- -1
      Ag[19, M + 24] <- 1/2;  Ag[19, M + 33] <- 1/2;  Ag[19, M + 28] <- -1
      Ag[20, M + 24] <- 1/2;  Ag[20, M + 33] <- -1/2;  Ag[20, M + 29] <- -1
      Ag[21, M + 27] <- 1/2;  Ag[21, M + 30] <- 1/2;  Ag[21, M + 31] <- -1
      Ag[22, M + 27] <- 1/2;  Ag[22, M + 30] <- -1/2;  Ag[22, M + 32] <- -1
    }
    
    if (m == 10) {
      cones.index <- vector("list", 11)
      cones.index[[1]] <- c(M + 1, M + 2, M + 3)
      cones.index[[2]] <- c(M + 4, M + 5, M + 6)
      cones.index[[3]] <- c(M + 7, M + 8, M + 9)
      cones.index[[4]] <- c(M + 10, M + 11, M + 12)
      cones.index[[5]] <- c(M + 13, M + 14, M + 15)
      cones.index[[6]] <- c(M + 16, M + 17, M + 18)
      cones.index[[7]] <- c(M + 19, M + 20, M + 21)
      cones.index[[8]] <- c(M + 22, M + 23, M + 24)
      cones.index[[9]] <- c(M + 25, M + 26, M + 27)
      cones.index[[10]] <- c(M + 28, M + 29, M + 30)
      cones.index[[11]] <- c(M + 31, M + 32, M + 33)
      cn <- 11; newvars <- 33; ncg <- 22
      Ag <- matrix(0, nrow = ncg, ncol = M + newvars)
      Ag[1, n + 1] <- 1;  Ag[1, n + m + 2] <- 1;  Ag[1, M + 1] <- -1
      Ag[2, n + 1] <- 1;  Ag[2, n + m + 2] <- -1;  Ag[2, M + 2] <- -1
      Ag[3, n + 2 * m + 3] <- 1;  Ag[3, n + 3 * m + 4] <- 1;  Ag[3, M + 4] <- -1
      Ag[4, n + 2 * m + 3] <- 1;  Ag[4, n + 3 * m + 4] <- -1;  Ag[4, M + 5] <- -1
      Ag[5, n + 4 * m + 5] <- 1;  Ag[5, n + 5 * m + 6] <- 1;  Ag[5, M + 7] <- -1 
      Ag[6, n + 4 * m + 5] <- 1;  Ag[6, n + 5 * m + 6] <- -1;  Ag[6, M + 8] <- -1
      Ag[7, n + 6 * m + 7] <- 1;  Ag[7, n + 7 * m + 8] <- 1;  Ag[7, M + 10] <- -1
      Ag[8, n + 6 * m + 7] <- 1;  Ag[8, n + 7 * m + 8] <- -1;  Ag[8, M + 11] <- -1
      Ag[9, n + 8 * m + 9] <- 1;  Ag[9, n + 9 * m + 10] <- 1;  Ag[9, M + 13] <- -1
      Ag[10, n + 8 * m + 9] <- 1;  Ag[10, n + 9 * m + 10] <- -1;  Ag[10, M + 14] <- -1
      Ag[11, M + 3] <- 1/2;  Ag[11, M + 6] <- 1/2;  Ag[11, M + 16] <- -1
      Ag[12, M + 3] <- 1/2;  Ag[12, M + 6] <- -1/2;  Ag[12, M + 17] <- -1
      Ag[13, M + 9] <- 1/2;  Ag[13, M + 12] <- 1/2;  Ag[13, M + 19] <- -1
      Ag[14, M + 9] <- 1/2;  Ag[14, M + 12] <- -1/2;  Ag[14, M + 20] <- -1
      Ag[15, M + 15] <- 1/2;  Ag[15, M + 33] <- 1/2;  Ag[15, M + 22] <- -1
      Ag[16, M + 15] <- 1/2;  Ag[16, M + 33] <- -1/2;  Ag[16, M + 23] <- -1
      Ag[17, M + 18] <- 1/2;  Ag[17, M + 21] <- 1/2;  Ag[17, M + 25] <- -1
      Ag[18, M + 18] <- 1/2;  Ag[18, M + 21] <- -1/2;  Ag[18, M + 26] <- -1
      Ag[19, M + 24] <- 1/2;  Ag[19, M + 33] <- 1/2;  Ag[19, M + 28] <- -1
      Ag[20, M + 24] <- 1/2;  Ag[20, M + 33] <- -1/2;  Ag[20, M + 29] <- -1
      Ag[21, M + 27] <- 1/2;  Ag[21, M + 30] <- 1/2;  Ag[21, M + 31] <- -1
      Ag[22, M + 27] <- 1/2;  Ag[22, M + 30] <- -1/2;  Ag[22, M + 32] <- -1
    }
    
    return(list(cones.index = cones.index,  newvars = newvars,  ncg = ncg,  Ag = Ag, cn = cn)) 
  }
  
  n <- nrow(Fx); m <- ncol(Fx)
  if (m > 10) stop("This procedure is only implemented for m<=10.")
  k1 <- length(b1); k2 <- length(b2); k3 <- length(b3)
  b <- c(b1, b2, b3); A <- rbind(A1, A2, A3)
  sense <- c(rep("<=", k1), rep(">=", k2), rep("=", k3))
  model  <-  list(); model$modelsense <- "max"
  da <- dim(A); nc <- da[1]; M <- n + m^2 + 4*m*n
  geo <- geomeanG(m, n, M)
  varnum <- geo$newvars + M
  
  nC <- c(m*m, 1/2 * m*(m - 1), m, nc, m*n, m*n, geo$ncg)
  sC <- sum(nC)
  Aeq <- matrix(0, nrow = sum(nC), ncol = varnum)
  for (i in 1:(m*m)) Aeq[i, i + n] <- -1 
  
  for (i in 1:m) {
    for (j in 1:n) {
      for (ii in 1:m) {
        Aeq[(i - 1)*m + ii, (n + m^2 + m*n + (j - 1)*m + ii)] <- 1/2*Fx[j, i] 
      }
    }
  }
  
  ind <- 1
  for (i in 1:(m - 1)) {
    for (j in i:(m - 1)) {
      Aeq[m*m + ind, n + (i - 1)*m + j + 1] <- 1
      ind <- ind + 1
    }
  }
  
  for (i in 1:m) {
    Aeq[m*m + m*(m - 1)/2 + i, n + (i - 1)*m + i] <- -1
    for (j in 1:n) {
      Aeq[m*m + m*(m - 1)/2 + i, n + m*m + (j - 1)*m + i] <- 1
    }
  }
  
  Aeq[(m^2 + 1/2 * m*(m - 1) + m + 1):(m^2 + 0.5*m*(m - 1) + m + nc), 1:n] <- A
  
  S <- sum(nC[1:4])
  for (i in 1:(n*m)) {
    Aeq[S + i, n + m^2 + i] <- 1
    Aeq[S + i, n + m^2 + 2*n*m + i] <- -1
    Aeq[S + i, floor((i - 1)/m) + 1] <- 1
  }
  
  S <- sum(nC[1:5])
  for (i in 1:(n*m)) {
    Aeq[S + i, n + m^2 + i] <- 1
    Aeq[S + i, n + m^2 + 3*n*m + i] <- -1
    Aeq[S + i, floor((i - 1)/m) + 1] <- -1
  }
  
  S <- sum(nC[1:6])
  if (nC[7] > 0) Aeq[(S + 1):sC, ] <- geo$Ag
  model$A <- Matrix::Matrix(Aeq, sparse = TRUE)
  rhs <- rep(0, sC); rhs[(sum(nC[1:3]) + 1):(sum(nC[1:4]))] <- b
  model$rhs <- rhs 
  
  model$sense <- rep("=", sC)
  for (i in (m^2 + 1/2 * m*(m - 1) + 1):(m^2 + 1/2 * m*(m - 1) + m)) 
    model$sense[i] <- "<="
  for (i in (m^2 + 1/2 * m*(m - 1) + m + 1):(m^2 + 1/2 * m*(m - 1) + m + nc)) 
    model$sense[i] <- sense[i - (m^2 + 1/2 * m*(m - 1) + m)]
  
  model$quadcon <- vector("list", n + geo$cn)
  for (i in 1:geo$cn) {
    qc <- list()
    ind <- geo$cones.index[[i]]
    qc$Qc <- Matrix::spMatrix(varnum, varnum, ind, ind, c(-1, 1, 1))
    qc$rhs <- 0
    model$quadcon[[i]] <- qc
  }
  
  for (i in 1:(m*n)) {
    qc <- list()
    qc$Qc <- Matrix::spMatrix(varnum, varnum, c(n + m^2 + 2*m*n + i, n + m^2 + 3*m*n + i, n + m^2 + m*n + i),
                               c(n + m^2 + 2*m*n + i, n + m^2 + 3*m*n + i, n + m^2 + m*n + i), c(-1, 1, 1))
    qc$rhs <- 0
    model$quadcon[[i + geo$cn]] <- qc
  }
   
  blx <- rep(-Inf, varnum); blx[1:n] <- 0
  blx[(n + m^2 + 1):(n + m^2 + n*m)] <- 0
  blx[(n + m^2 + 2*m*n + 1):(n + m^2 + 3*m*n)] <- 0
  
  # 1st indices of all cones in geomean
  if (m == 2) blx[M + 1] <- 0
  if (m == 3) {
    blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
  }
  if (m == 4) {
    blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
  }
  if (m == 5) {
    blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
    blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0
  }
  if (m == 6) {
    blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
    blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0
  }
  if (m == 7) {
    blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
    blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0
    blx[M + 19] <- 0
  }
  if (m == 8) {
    blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0; 
    blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0
    blx[M + 19] <- 0
  }
  if (m == 9) {
    blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0
    blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0
    blx[M + 19] <- 0 
    blx[M + 22] <- 0; blx[M + 25] <- 0; blx[M + 28] <- 0
    blx[M + 31] <- 0
  }
  if (m == 10) {
    blx[M + 1] <- 0; blx[M + 4] <- 0; blx[M + 7] <- 0 
    blx[M + 10] <- 0; blx[M + 13] <- 0; blx[M + 16] <- 0
    blx[M + 19] <- 0
    blx[M + 22] <- 0; blx[M + 25] <- 0; blx[M + 28] <- 0
    blx[M + 31] <- 0
  }
  
  bux <- rep(Inf, varnum)
  model$lb <- blx; model$ub <- bux
  model$lb[1:n] <- 0
  
  cc <- rep(0, varnum); cc[varnum] <- 1
  model$obj <- cc
  
  vtypes <- rep("C", varnum)
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
