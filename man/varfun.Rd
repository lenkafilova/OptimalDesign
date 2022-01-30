\name{varfun}
\alias{varfun}

\title{
Vector of variances
}

\description{
Computes the vector of variances (sensitivities) for a given design \code{w}.
}

\usage{
varfun(Fx, w, crit="D", h=NULL, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} matrix of candidate regressors (as rows), where \code{n} is the number of candidate design points and \code{m} (where \code{m>=2}, \code{m<=n}) is the number of parameters.
}
\item{w}{
a non-negative vector of length \code{n} representing the design.
}
\item{crit}{
the criterion; possible values are \code{"D"}, \code{"A"}, \code{"I"}, \code{"C"} and \code{"c"}.
}
\item{h}{
a non-zero vector of length \code{m} corresponding to the coefficients of the linear parameter combination of interest. If crit is not \code{"C"} nor \code{"c"} then \code{h} is ignored. If \code{crit} is \code{"C"} or \code{"c"} and \code{h=NULL} then \code{h} is assumed to be \code{c(0,...,0,1)}.
}
\item{echo}{
Print the call of the function?
}
}

\details{For D-optimality, the i-th element of the vector of variances is the variance of the best linear unbiased estimator of the mean value of observations under the experimental conditions represented by the i-th design point (where the variance of the observational errors is assumed to be 1). There is a linear transformation relation of the vector of variances and the vector of directional derivatives for the criterion of D-optimality. See the reference paper at \url{http://www.iam.fmph.uniba.sk/design/} for mathematical details.
}

\value{
The vector of variances (sensitivities) for a given design \code{w}.
}


\author{
Radoslav Harman, Lenka Filova
}

\note{
The design \code{w} should have a non-singular information matrix.
}

\seealso{
\code{\link{effbound}, \link{dirder}}
}

\examples{
# The values of the variance function (for crit=D)
# of D-, I-, and C-optimal approximate design

Fx <- Fx_cube(~x1 + I(x1^2), n.levels = 21)
wD <- od_REX(Fx)$w.best
wI <- od_REX(Fx, crit="I")$w.best
wC <- od_REX(Fx, crit="C", h=c(1, 0, 0))$w.best
vD <- varfun(Fx, wD)
vI <- varfun(Fx, wI)
vC <- varfun(Fx, wC)
plot(Fx[, 2], rep(0, nrow(Fx)), ylim = c(0, max(vD, vI, vC)),
     type = "n", xlab = "x", ylab = "var", lwd = 2)
grid()
lines(Fx[, 2], vD, col = "red")
lines(Fx[, 2], vI, col = "blue")
lines(Fx[, 2], vC, col = "green")

# The D-optimal approximate design minimized the maximum
# of the var. function (it is "G-optimal").
# The I-optimal approximate design minimizes the integral of the var. function.
# The C-optimal design with h=f(0) makes the var. function small around 0.
}

