infmat <- function(Fx, w, echo=TRUE) {
  # information matrix of a possibly non-normalized design w
   
  cl <- match.call()
  verify(cl, Fx = Fx, w = w, echo = echo)
  n <- nrow(Fx); m <- ncol(Fx)
   
  supp <- w > 0; s <- sum(supp)
  if (s > 1) {
    if (s*(m + 5) < 1.1*n*(m - 1)) {
      res <- crossprod(sqrt(w[supp])*Fx[supp,])
    } else {
      res <- crossprod(sqrt(w)*Fx)
    }
  } else {
    res <- w[supp]*Fx[supp,] %*% t(Fx[supp,])
  }
  
  return(as.matrix(res))
}

