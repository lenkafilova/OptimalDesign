od_DEL <- function(Fx, w, crit="D", h=NULL, echo=TRUE) {
  # Removes the design points that cannot support
  # an approximate optimal design under the standard constraint
  
  cl <- match.call()
  verify(cl, Fx = Fx, w = w, crit = crit, echo = echo)
  n <- nrow(Fx); m <- ncol(Fx)
  
  if (!isTRUE(all.equal(sum(w), 1)))
    warning(paste("sum(w) is ", sum(w),"!=1; w will be normalized."))
  w <- w/sum(w)
  if (crit == "C" && is.null(h)) h <- c(rep(0, m - 1), 1)
  if (crit == "c")
    stop("The pure c-optimality is not implemented for DEL.")
  
  if (crit == "D") keep <- od_D_DEL(Fx = Fx, w = w)
  if (crit == "A") keep <- od_A_DEL(Fx = Fx, w = w)
  if (crit == "I") keep <- od_A_DEL(Fx = Fx_ItoA(Fx, echo = FALSE), w = w)
  if (crit == "C") keep <- od_A_DEL(Fx = Fx_CtoA(Fx, h = h, echo = FALSE), w = w)
  
  return(list(call = cl, keep = keep, w.keep = w[keep]/sum(w[keep]), 
              Fx.keep = Fx[keep, ]))
}
