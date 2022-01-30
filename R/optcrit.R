optcrit <- function(Fx, w, crit="D", h=NULL, echo=TRUE) {
   # optimality criterion of a possibly non-normalized design w
   
   cl <- match.call()
   verify(cl, Fx = Fx, w = w, crit = crit, h = h, echo = echo) 
   m <- ncol(Fx)
   
   if (crit %in% c("D", "A", "c"))
      M <- infmat(Fx, w, echo = FALSE)
   if (crit == "I") {
      M <- infmat(Fx_ItoA(Fx, echo = FALSE), w, echo = FALSE)
   } 
   if (crit %in% c("C", "c")) {
      if (is.null(h)) {
         if (echo) message("Setting h <- c(0,...,0,1).")
         h <- c(rep(0, m - 1), 1)
      }
      if (crit == "C")
         M <- infmat(Fx_CtoA(Fx, h, echo = FALSE), w, echo = FALSE)
   }
   
   tol <- 1e-12 #tol <- sqrt(.Machine$double.eps)
   if (rcond(M) < tol) {
      if (!(crit == "c")) { 
         return(0)
      } else {
         M.plus <- ginv(M)
         if (max(abs(M %*% M.plus %*% h - h)) > tol) return(0)
         val <- sum(h^2) / (t(h) %*% M.plus %*% h)[1, 1]
      }
   } else {
      if (crit == "D") val <- (abs(det(M)))^(1/m)
      if (crit %in% c("A", "I", "C")) val <- m * sum(diag(solve(M)))^(-1)
      if (crit == "c") val <- sum(h^2) / (t(h) %*% solve(M) %*% h)[1, 1]
   }
   
   return(as.numeric(val))
}
