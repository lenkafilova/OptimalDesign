Fx_glm <- function(formula, theta0, glm.model="bin-logit",
                   lower=NULL, upper=NULL, n.levels=NULL, echo=TRUE) {
  # Generate Fx for a glm model

  cl <- match.call()
  verify(cl, formula = formula, theta0 = theta0, glm.model = glm.model,
         lower = lower, upper = upper, n.levels = n.levels, echo = echo)

  if (glm.model == "bin-logit")
    u <- function(eta) exp(eta)/(1 + exp(eta))^2
  if (glm.model == "bin-probit")
    u <- function(eta) pnorm(eta)^2/(pnorm(eta)*(1 - pnorm(eta)))
  if (glm.model == "bin-cloglog")
    u <- function(eta) exp(-exp(eta))/(1 - exp(-exp(eta)))*exp(eta)^2
  if (glm.model == "Poisson-log")
    u <- function(eta) exp(eta)

  F.lin <- Fx_cube(formula, lower, upper, n.levels, echo = FALSE)
  n <- nrow(F.lin); m <- ncol(F.lin) 
  Fx <- matrix(0, nrow = n, ncol = m)
  for (i in 1:n) {
    f <- as.vector(F.lin[i, ])
    eta <- as.numeric(t(f) %*% theta0)
    Fx[i, ] <- sqrt(u(eta)) * f
  }
  cnms <- rep("", m)
  for (j in 1:m) cnms[j] <- paste("G", j, sep = "")
  colnames(Fx) <- cnms
    
  return(Fx)
}
