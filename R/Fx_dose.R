Fx_dose <- function(dose.levels, theta0, dose.model="emax", echo=TRUE) {
  # Generate Fx for a dose-response model
 
  cl <- match.call()
  verify(cl, dose.model = dose.model, dose.levels = dose.levels,
         theta0 = theta0, echo = echo)

  if (dose.model == "emax")
    g <- function(x) 
      c(1, x/(x + theta0[3]), -theta0[2]/(x + theta0[3])^2)
  if (dose.model == "loglin")
    g <- function(x) c(1, log(x + theta0[3]), theta0[2]/(x + theta0[3]))
  if (dose.model == "exp")
    g <- function(x) c(1, exp(x/theta0[3]), -theta0[2]*x*exp(x/theta0[3])/theta0[3]^2)

  n <- length(dose.levels); m <- length(theta0)  
  Fx <- matrix(0, nrow = n, ncol = m)
  for (i in 1:n) Fx[i, ] <- g(dose.levels[i])
  cnms <- rep("", m)
  for (j in 1:m) cnms[j] <- paste("D", j, sep = "")
  colnames(Fx) <- cnms

  return(Fx)
}
