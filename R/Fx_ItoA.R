Fx_ItoA <- function(Fx, echo=TRUE) {
  # Transformation of Fx for I-optimality  
  
  cl <- match.call()
  verify(cl, Fx = Fx, echo = echo)
  n <- nrow(Fx); m <- ncol(Fx)
  L <- m*infmat(Fx, rep(1, n), echo = FALSE)/sum(Fx^2)
  if (rcond(L) <  sqrt(.Machine$double.eps))
    warning("The problem of I-optimality may be badly conditioned.")
  
  return(as.matrix(Fx %*% solve(chol(L))))
}
