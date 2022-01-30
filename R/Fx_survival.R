Fx_survival <- function(formula, theta0, censor.time, survival.model="phI",
                        lower=NULL, upper=NULL, n.levels=NULL, echo=TRUE) {
  # Generate Fx for a survival model
 
  cl <- match.call()
  verify(cl, formula = formula, survival.model = survival.model, theta0 = theta0,
         censor.time = censor.time, lower = lower, upper = upper,
         n.levels = n.levels, echo = echo)
  
  if (survival.model == "phI")
    u <- function(f) 1 - exp(-censor.time * exp(t(f) %*% theta0))
  if (survival.model == "phrand")
    u <- function(f) 1 - (1 - exp(-censor.time * exp(t(f) %*% theta0))) /
    (censor.time * exp(t(f) %*% theta0))

  F.lin <- Fx_cube(formula, lower, upper, n.levels, echo = FALSE)
  n <- nrow(F.lin); m <- ncol(F.lin)
  Fx <- matrix(0, nrow = n, ncol = m)  
  for (i in 1:n) Fx[i, ] <- sqrt(u(F.lin[i, ])) %*% F.lin[i, ]
  cnms <- rep("", m)
  for (j in 1:m) cnms[j] <- paste("S", j, sep = "")
  colnames(Fx) <- cnms
  
  return(Fx)
}
