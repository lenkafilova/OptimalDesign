od_print <- function(Fx, w, X=NULL, h=NULL, echo=TRUE) {
  # Print a compact info about the design w

  cl <- match.call()
  verify(cl, Fx = Fx, w = w, X = X, h = h, echo = echo)
  n <- nrow(Fx)
  
  if (is.null(X)) {
    message("X is NULL; setting X to 1:nrow(Fx)."); X <- 1:n
  }
  X <- as.matrix(X)
  
  ind.supp <- (1:n)[w > sqrt(.Machine$double.eps)]
  n.supp <- length(ind.supp)
  if (n.supp > 1) {
    design <- cbind(X[ind.supp,], weight = round(w[ind.supp], 8))
  } else {
    design <- cbind(t(X[ind.supp,]), weight = 1)
  }
  rownames(design) <- 1:n.supp
    
  D.value <- optcrit(Fx, w, crit = "D", echo = FALSE)
  A.value <- optcrit(Fx, w, crit = "A", echo = FALSE)
  I.value <- optcrit(Fx, w, crit = "I", echo = FALSE)
  C.value <- optcrit(Fx, w, crit = "C", h = h, echo = FALSE)
  c.value <- optcrit(Fx, w, crit = "c", h = h, echo = FALSE)
  
  return(list(call = cl, design = design, M = infmat(Fx, w, echo = FALSE),
              eigenvalues = eigen(infmat(Fx, w, echo = FALSE), symmetric = TRUE)$values,
              D.value = D.value, A.value = A.value, I.value = I.value,
              C.value = C.value, c.value = c.value))
}
