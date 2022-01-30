effbound <- function(Fx, w, crit="D", h=NULL, echo=TRUE) {
   
   # Efficiency bound on the design w under the size constaint

   cl <- match.call()
   verify(cl, Fx = Fx, w = w, crit = crit, h = h, echo = echo)
   m <- ncol(Fx); w.norm <- w/sum(w)
   
   if (crit %in% c("D", "A", "c"))
      M <- infmat(Fx, w.norm, echo = FALSE)
   if (crit == "I") {
      M <- infmat(Fx_ItoA(Fx, echo = FALSE), w.norm, echo = FALSE)
   } 
   if (crit %in% c("C", "c")) {
      if (is.null(h)) {
         if (echo) message("Setting h <- c(0,...,0,1).")
         h <- c(rep(0, m - 1), 1)
      }
      if (crit == "C")
         M <- infmat(Fx_CtoA(Fx, h, echo = FALSE), w.norm, echo = FALSE)
   }
   
   tol <- 1e-12 #sqrt(.Machine$double.eps)
   if (rcond(M) < tol) return(0)
   
   if (crit == "D") 
      eff.bnd <- m/max(varfun(Fx, w.norm, crit = "D", h, echo = FALSE))
   if (crit == "A")
      eff.bnd <- sum(diag(solve(M)))/max(varfun(Fx, w.norm, crit = "A", h, echo = FALSE))
   if (crit == "I") { 
      FxI <- Fx_ItoA(Fx, echo = FALSE)
      eff.bnd <- sum(diag(solve(M)))/max(varfun(FxI, w.norm, crit = "A", h, echo = FALSE))
   }
   if (crit == "C") { 
      FxC <- Fx_CtoA(Fx, h, echo = FALSE)
      eff.bnd <- sum(diag(solve(M)))/max(varfun(FxC, w.norm, crit = "A", h, echo = FALSE))
   }
   if (crit == "c")
      eff.bnd <- (t(h) %*% solve(M) %*% h) / max(varfun(Fx, w.norm, crit = "c", h = h, echo = FALSE))

   return(as.numeric(eff.bnd))
}
