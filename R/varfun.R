varfun <- function(Fx, w, crit="D", h=NULL, echo=TRUE) {

   # The "variation function" vector. For crit="D", the i-th component is 
   # sigma^{-2} times the variance of the BLUE of f_i'theta 
   # an increasing linear transformation of varfun is the directional derivative
   
   cl <- match.call() 
   verify(cl, Fx = Fx, w = w, crit = crit, h = h, echo = echo) 
   m <- ncol(Fx)

   if (crit %in% c("C", "c") && is.null(h)) {
      if (echo) message("Setting h <- c(0,...,0,1).")
      h <- c(rep(0, m - 1), 1)
   }

   M <- infmat(Fx, w, echo = FALSE)   
   if (rcond(M) < sqrt(.Machine$double.eps)) 
      warning("The information matrix is badly conditioned and there may be problems computing the variance function.")
   
   if (crit == "D")
      var.fun <- ((Fx %*% t(chol(solve(M))))^2) %*% rep(1, m)
   if (crit == "A")
      var.fun <- ((Fx %*% solve(M))^2) %*% rep(1, m)
   if (crit == "c")
      var.fun <- (Fx %*% ginv(M, tol = 1e-6) %*% h)^2
   if (crit == "I")
      var.fun <- varfun(Fx_ItoA(Fx, echo = FALSE), w, crit = "A", echo = FALSE)
   if (crit == "C")
      var.fun <- varfun(Fx_CtoA(Fx, h, echo = FALSE), w, crit = "A", h = h, echo = FALSE)

   return(as.vector(var.fun))
}
