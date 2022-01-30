od_MISOCP <- function(Fx, b1=NULL, A1=NULL, b2=NULL, A2=NULL, b3=NULL, A3=NULL,
                      w0=NULL, bin=FALSE, type="exact", crit="D", h=NULL, gap=NULL,
                      t.max=120, echo=TRUE) {
  # Computes optimal exact or approximate design under linear constraints
  
  cl <- match.call()
  verify(cl, Fx = Fx, b1 = b1, A1 = A1, b2 = b2, A2 = A2, b3 = b3, A3 = A3, 
         w0 = w0, bin = bin, type = type, crit = crit, h = h, gap = gap, 
         t.max = t.max, echo = echo)
  
  n <- nrow(Fx); m <- ncol(Fx)
  if (!is.null(b1)) if (is.null(A1)) A1 <- matrix(1, nrow = length(b1), ncol = n)
  if (!is.null(b2)) if (is.null(A2)) A2 <- matrix(1, nrow = length(b2), ncol = n)
  if (!is.null(b3)) if (is.null(A3)) A3 <- matrix(1, nrow = length(b3), ncol = n)
  if (is.null(w0)) w0 <- rep(0, n)

  if (crit %in% c("C", "c") && is.null(h)) h <- c(rep(0, m - 1), 1)  
  
  if (!is.null(w0) && sum(w0) > 0) {
    A2 <- rbind(A2, diag(n)[w0 > 0,])
    b2 <- c(b2, w0[w0 > 0])
  } 
  
  info <- paste("Running od_MISOCP for cca", t.max, "seconds")
  info <- paste(info, " starting at ", Sys.time(), ".", sep = "")
  print(info, quote = FALSE)
  info <- paste("The problem size is n=", n, sep = "")
  info <- paste(info, ", m=", m, sep = "")
  info <- paste(info, ", k=", length(c(b1, b2, b3)), ".", sep = "")
  print(info, quote = FALSE)
  
  start <- as.numeric(proc.time()[3])

  if (crit == "D")
    res <- od_D_MISOCP(Fx, b1, A1, b2, A2, b3, A3, bin, type, gap, t.max)
  if (crit == "A")
    res <- od_A_MISOCP(Fx, b1, A1, b2, A2, b3, A3, bin, type, gap, t.max)
  if (crit == "I")
    res <- od_A_MISOCP(Fx_ItoA(Fx, echo = FALSE), b1, A1, b2, A2, b3, A3, bin,
                       type, gap, t.max)
  if (crit == "C")
    res <- od_A_MISOCP(Fx_CtoA(Fx, h, echo = FALSE), b1, A1, b2, A2, b3, A3,
                     bin, type, gap, t.max)
  if (crit == "c")
    res <- od_c_MISOCP(Fx, b1, A1, b2, A2, b3, A3, bin, h, type, gap, t.max)
  
  w <- res$w.best
  if (!is.null(w)) {
    supp <- (1:n)[w > 1e-5]; w.supp <- w[supp]
    M.best <- infmat(Fx, w, echo = FALSE)
    Phi.best <- optcrit(Fx, w, crit = crit, h = h, echo = FALSE)
    err <- c()
    if (!is.null(b1)) err <- c(err, pmin(A1 %*% w - b1, 0))
    if (!is.null(b2)) err <- c(err, pmin(b2 - A2 %*% w, 0))
    if (!is.null(b3)) err <- c(err, abs(A3 %*% w - b3))
    if (max(err) > 1e-05)
      warning(cat("Some constraints are significantly violated:", round(err, 6)))
  } else {
    supp <- NULL; w.supp <- NULL
    M.best <- NULL; Phi.best <- 0
    warning("MISOCP was not able to find any meaningful solution within the alloted time.")
  }  
  
  t.act <- round(as.numeric(proc.time()[3]) - start, 2)
  
  return(list(call = cl, w.best = w, supp = supp, w.supp = w.supp, M.best = M.best,
              Phi.best = Phi.best, status = res$status, t.act = t.act))
}
