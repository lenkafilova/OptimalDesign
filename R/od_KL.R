od_KL <- function(Fx, N, bin=FALSE, Phi.app=NULL, crit="D", h=NULL, w1=NULL,
                  K=NULL, L=NULL, rest.max=Inf, t.max=120, echo=TRUE, track=TRUE) {
  # Computes an efficient exact design under the standard size constraint
  # TODO Add the argument w0
  
  cl <- match.call()
  verify(cl, Fx = Fx, N = N, bin = bin, Phi.app = Phi.app, crit = crit, h = h, w1 = w1,
         K = K, L = L, rest.max = rest.max, t.max = t.max, echo = echo, track = track) 
  
  n <- nrow(Fx); m <- ncol(Fx)
  if (crit == "C" && is.null(h)) h <- c(rep(0, m - 1), 1)
  if (crit == "c")
    stop("The pure c-optimality is not implemented for KL. Try its regularized version: the C-optimality.") 
  if (!is.null(w1) && sum(w1) != N) 
    stop("The initial design is not feasible.")
  
  if (bin && n < N) {
    print("The required size N is:"); print(N)
    print("Fx is:"); print.default(Fx, max = 100)
    stop("There is no replication-free size-N design at n<N design points.")
  }
  if (bin && n == N) {
    w.best <- rep(1, n)
    res <- list(call = cl, w.best = w.best, supp = 1:n, w.supp = rep(1, n),   
                M.best = infmat(Fx, w.best, echo = FALSE),
                Phi.best = optcrit(Fx, w.best, crit = crit, echo = FALSE),
                eff.best = 1, n.rest = 0, n.ex = 0, t.act = 0)
    return(res)
  }
    
  if (bin == FALSE) {
    if (crit == "D")
      res <- od_D_KL(Fx = Fx, N = N, Phi.app = Phi.app, w1 = w1, K = K, L = L,
                     rest.max = rest.max, t.max = t.max, track = track)  
    if (crit == "A")
      res <- od_A_KL(Fx = Fx, N = N, Phi.app = Phi.app, w1 = w1, K = K, L = L,
                     rest.max = rest.max, t.max = t.max, track = track)
    if (crit == "I")
      res <- od_A_KL(Fx = Fx_ItoA(Fx, echo = FALSE), N = N, Phi.app = Phi.app,
                     w1 = w1, K = K, L = L, rest.max = rest.max, t.max = t.max, track = track) 
    if (crit == "C")
      res <- od_A_KL(Fx = Fx_CtoA(Fx, h, echo = FALSE), N = N, Phi.app = Phi.app,
                     w1 = w1, K = K, L = L, rest.max = rest.max, t.max = t.max, track = track)
  } else {
    if (crit == "D")
      res <- od_D_KL_bin(Fx = Fx, N = N, Phi.app = Phi.app, w1 = w1, K = K, L = L,
                         rest.max = rest.max, t.max = t.max, track = track)  
    if (crit == "A")
      res <- od_A_KL_bin(Fx = Fx, N = N, Phi.app = Phi.app, w1 = w1, K = K, L = L,
                         rest.max = rest.max, t.max = t.max, track = track)
    if (crit == "I")
      res <- od_A_KL_bin(Fx = Fx_ItoA(Fx, echo = FALSE), N = N, Phi.app = Phi.app,
                         w1 = w1, K = K, L = L, rest.max = rest.max, t.max = t.max, track = track) 
    if (crit == "C")
      res <- od_A_KL_bin(Fx = Fx_CtoA(Fx, h, echo = FALSE), N = N, Phi.app = Phi.app,
                         w1 = w1, K = K, L = L, rest.max = rest.max, t.max = t.max, track = track)
  }
  
  supp <- (1:n)[res$w.best > 0]
  res <- list(call = cl, w.best = res$w.best, supp = supp, w.supp = res$w.best[supp],   
              M.best = infmat(Fx, res$w.best, echo = FALSE), Phi.best = res$Phi.best,
              eff.best = res$eff.best, n.rest = res$n.rest, n.ex = res$n.ex, t.act = res$t.act)
  
  return(res)
}
