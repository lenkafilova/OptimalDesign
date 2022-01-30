od_REX <- function(Fx, crit="D", h=NULL, w1=NULL, alg.AA="REX",
                   eff=0.999999, it.max=Inf, t.max=60, echo=TRUE, track=TRUE) {
  # The REX procedure for optimal approxmimate design 
  # Includes the implementation of the multiplicative algorithm
  # Also includes the vertex direction method for teaching purposes
  
  cl <- match.call()
  verify(cl, Fx = Fx, crit = crit, h = h, w1 = w1, alg.AA = alg.AA,
         eff = eff, it.max = it.max, t.max = t.max, echo = echo, track = track)
  
  n <- nrow(Fx); m <- ncol(Fx)
  if (crit %in% c("C", "c") && is.null(h)) h <- c(rep(0, m - 1), 1)   
  if (!is.null(w1) && sum(w1) != 1) {
    message("w1 not perfectly normalized; normalizing.") 
    w1 <- w1 / sum(w1)
  }
  
  # Use the appropriate engine
  if (crit == "D") {
    if (alg.AA == "REX") 
      res <- od_D_REX(Fx = Fx, w1 = w1, ver = 1, gamma = 4, eff = eff, 
                      it.max = it.max, t.max = t.max, track = track) 
    if (alg.AA == "MUL") 
      res <- od_D_MUL(Fx = Fx, w1 = w1, eff = eff, it.max = it.max,
                      t.max = t.max, track = track) 
    if (alg.AA == "VDM") 
      res <- od_D_VDM(Fx = Fx, w1 = w1, eff = eff, it.max = it.max,
                      t.max = t.max, track = track)
  } 
  
  if (crit == "A") {
    if (alg.AA == "REX") 
      res <- od_A_REX(Fx = Fx, w1 = w1, ver = 1, gamma = 1, eff = eff, 
                      it.max = it.max, t.max = t.max, track = track) 
    if (alg.AA == "MUL") 
      res <- od_A_MUL(Fx = Fx, w1 = w1, lambda = 0.5, eff = eff,
                      it.max = it.max, t.max = t.max, track = track) 
    if (alg.AA == "VDM") 
      res <- od_A_VDM(Fx = Fx, w1 = w1, eff = eff, it.max = it.max,
                      t.max = t.max, track = track)
  } 
  
  if (crit == "I") {
    if (alg.AA == "REX") 
      res <- od_A_REX(Fx = Fx_ItoA(Fx, echo = FALSE), w1 = w1, ver = 1,
                      gamma = 1, eff = eff, it.max = it.max, t.max = t.max, track = track) 
    if (alg.AA == "MUL") 
      res <- od_A_MUL(Fx = Fx_ItoA(Fx, echo = FALSE), w1 = w1, lambda = 0.5, 
                      eff = eff, it.max = it.max, t.max = t.max, track = track) 
    if (alg.AA == "VDM") 
      res <- od_A_VDM(Fx = Fx_ItoA(Fx, echo = FALSE), w1 = w1, eff = eff, 
                      it.max = it.max, t.max = t.max, track = track)
  } 
  
  if (crit == "C") {
    if (alg.AA == "REX") 
      res <- od_A_REX(Fx = Fx_CtoA(Fx, h, echo = FALSE), w1 = w1, ver = 1,
                      gamma = 1, eff = eff, it.max = it.max, t.max = t.max, track = track) 
    if (alg.AA == "MUL") 
      res <- od_A_MUL(Fx = Fx_CtoA(Fx, h, echo = FALSE), w1 = w1, lambda = 0.5, 
                      eff = eff, it.max = it.max, t.max = t.max, track = track) 
    if (alg.AA == "VDM") 
      res <- od_A_VDM(Fx = Fx_CtoA(Fx, h, echo = FALSE), w1 = w1, eff = eff, 
                      it.max = it.max, t.max = t.max, track = track)
  } 
  
  if (crit == "c") {
    message("Linear programming reformulation will be used for c-optimality.")
    res <- od_c_LP(Fx, h, track)
  }
  
  # Output
  supp <- (1:n)[res$w.best > 0]
  
  res <- list(call = cl, w.best = res$w.best, supp = supp, w.supp = res$w.best[supp],   
              M.best = infmat(Fx, res$w.best, echo = FALSE),
              Phi.best = optcrit(Fx, res$w.best, crit, h, echo = FALSE),
              eff.best = res$eff.best, n.iter = res$n.iter, t.act = res$t.act)
  return(res)
}

