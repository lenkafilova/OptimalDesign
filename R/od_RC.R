od_RC <- function(Fx, b, A = NULL, w0 = NULL, bin = FALSE, Phi.app = NULL,
                  crit = "D", h=NULL, w1 = NULL, rest.max = Inf,
                  t.max = 120, echo = TRUE, track=TRUE) {
   # Computes an efficient exact design under multiple resource constraints 
   
   cl <- match.call()
   verify(cl, Fx = Fx, b = b, A = A, w0 = w0, bin = bin, Phi.app = Phi.app,
          crit = crit, h = h, w1 = w1, rest.max = rest.max,
          t.max = t.max, echo = echo, track = track) 
   
   n <- nrow(Fx); m <- ncol(Fx)
   if (is.null(A)) A <- matrix(1, nrow = 1, ncol = n)
   if (is.null(w0)) w0 <- rep(0, n)
   if (bin) {
      A <- rbind(A, diag(n))
      b <- c(b, rep(1, n))
   }
   if (any(A %*% w0 > b))
      stop("w0 must satisfy the resource constraints.")
   if (crit == "C" && is.null(h)) h <- c(rep(0, m - 1), 1)  
   if (crit == "c")
      stop("The pure c-optimality is not implemented for RC. Try its regularized version: the C-optimality.")
   if (is.null(w1)) w1 <- w0
   if (is.null(Phi.app))
      message("Phi.app not supplied. Therefore, the lower bound on efficieny will not be reported.")

   if (crit == "D") 
      res <- od_D_RC(Fx = Fx, b = b, A = A, w0 = w0, Phi.app = Phi.app, 
                     w1 = w1, rest.max = rest.max, t.max = t.max, track = track)
   if (crit == "A") 
      res <- od_A_RC(Fx = Fx, b = b, A = A, w0 = w0, Phi.app = Phi.app, 
                     w1 = w1, rest.max = rest.max, t.max = t.max, track = track)
   if (crit == "I") 
      res <- od_A_RC(Fx = Fx_ItoA(Fx, echo = FALSE), b = b, A = A, 
                     w0 = w0, Phi.app = Phi.app, w1 = w1, rest.max = rest.max, 
                     t.max = t.max, track = track)
   if (crit == "C") 
      res <- od_A_RC(Fx = Fx_CtoA(Fx, h, echo = FALSE), b = b, A = A, 
                     w0 = w0, Phi.app = Phi.app, w1 = w1, rest.max = rest.max, 
                     t.max = t.max, track = track)
   
   supp <- (1:n)[res$w.best > 0]
   res <- list(call = cl, w.best = res$w.best, supp = supp, 
               w.supp = res$w.best[supp], M.best = infmat(Fx, res$w.best, echo = FALSE), 
               Phi.best = res$Phi.best, eff.best = res$eff.best, n.rest = res$n.rest, 
               t.act = res$t.act)

   return(res)
}

