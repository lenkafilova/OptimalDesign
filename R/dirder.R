dirder <- function(Fx, w, crit="D", h=NULL, echo=TRUE) {

   cl <- match.call() 
   verify(cl, Fx = Fx, w = w, crit = crit, h = h, echo = echo) 
   m <- ncol(Fx); w.norm <- w/sum(w)
   
   if (crit %in% c("C", "c") && is.null(h)) {
      if (echo) message("Setting h <- c(0,...,0,1).")
      h <- c(rep(0, m - 1), 1)
   }
   
   M <- infmat(Fx, w.norm, echo = FALSE)
   if (rcond(M) < sqrt(.Machine$double.eps)) 
      warning("The information matrix is badly conditioned and there may be problems computing the directional derivative.")
   
   if (crit == "D") 
      dir.der <- ((det(M))^(1/m))*(((Fx %*% t(chol(solve(M))))^2) %*% rep(1, m) - m)
   if (crit == "A") {
      trM1 <- sum(diag(solve(M)))
      dir.der <- m*(((Fx %*% solve(M))^2) %*% rep(1, m) - trM1)/(trM1)^2
   }
   if (crit == "c") {
      htM1h <- as.numeric(t(h) %*% ginv(M) %*% h)
      dir.der <- sum(h^2)*((Fx %*% ginv(M, tol = 1e-6) %*% h)^2 - htM1h)*(htM1h)^2
   }
   if (crit == "I")
      dir.der <- dirder(Fx_ItoA(Fx, echo = FALSE), w.norm, crit = "A", echo = FALSE)
   if (crit == "C")
      dir.der <- dirder(Fx_CtoA(Fx, h = h, echo = FALSE), w.norm, crit = "A", h = h, echo = FALSE)
   
   return(as.vector(dir.der))
}

