\name{od_DEL}
\alias{od_DEL}

\title{
Removal of redundant design points
}
\description{
Removes the design points (or, equivalently, candidate regressors) that cannot support an optimal approximate design.
}
\usage{
od_DEL(Fx, w, crit = "D", h=NULL, echo = TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} (where \code{m>=2}, \code{m<=n}) matrix containing all candidate regressors (as rows), i.e., \code{n} is the number of candidate design points, and \code{m} is the number of parameters
}
\item{w}{
a non-negative vector of length \code{n} representing the design
}
\item{crit}{
the optimality criterion. Possible values are \code{"D"}, \code{"A"}, \code{"I"}, \code{"C"}.
}
\item{h}{
a non-zero vector of length \code{m} corresponding to the coefficients of the linear parameter combination of interest. If \code{crit} is not \code{"C"} nor \code{"c"} then \code{h} is ignored. If \code{crit} is \code{"C"} or \code{"c"} and \code{h=NULL} then \code{h} is assumed to be \code{c(0,...,0,1)}.
}
\item{echo}{
Print the call of the function?
}
}


\value{
  Output is the list with components:
  \item{call}{the call of the function}
  \item{keep}{the indices of \code{w} that have not been removed}
  \item{w.keep}{the approximate design on the reduced space}
  \item{Fx.keep}{the model matrix of the regressors on the reduced space}
}

\references{
Harman R, Pronzato L (2007): Improvements on removing non-optimal support points in D-optimum design algorithms, Statistics & Probability Letters 77, 90-94

Pronzato L (2013): A delimitation of the support of optimal designs for Kiefers Phi_p-class of criteria. Statistics & Probability Letters 83, 2721-2728}

\author{
Radoslav Harman, Lenka Filova
}

\note{
The design vector \code{w} should have a non-singular information matrix. The procedure is valid only for the standard (size) constraint.
}


\examples{
\dontrun{
# Generate a model matrix for the quadratic model
# on a semi-circle with a huge number of design points
form.q <- ~x1 + x2 + I(x1^2) + I(x2^2) + I(x1*x2)
Fx <- Fx_cube(form.q, lower = c(-1, 0), n.levels = c(1001, 501))
remove <- (1:nrow(Fx))[Fx[ ,2]^2 + Fx[ ,3]^2 > 1]
Fx <- Fx[-remove, ]

# Compute an approximate design w with an efficiency of cca 0.999
w <- od_REX(Fx, eff = 0.999)$w.best

# Remove the redundant design points based on w
Fx <- od_DEL(Fx, w)$Fx.keep

# Now an almost perfect design can be computed very rapidly:
w <- od_REX(Fx, eff = 0.9999999999)$w.best

# Plotting of the relevant directional derivative is also faster:
od_plot(Fx, w, Fx[ , 2:3], dd.size = 0.1)
}
}

