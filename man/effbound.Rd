\name{effbound}
\alias{effbound}

\title{
Lower bound on efficiency
}

\description{
Computes a lower bound on the efficiency of a design \code{w} in the class of all approximate designs of the same size as \code{w}.
}

\usage{
effbound(Fx, w, crit="D", h=NULL, echo=TRUE)
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
a non-zero vector of length \code{m} corresponding to the coefficients of the linear parameter combination of interest. If \code{crit} is not \code{"C"} nor \code{"c"} then \code{h} is ignored. If \code{crit} is \code{"C"} or \code{"c"} and \code{h=NULL} then \code{h} is assumed to be \code{c(0,...,0,1)}.
}
\item{echo}{
Print the call of the function?
}
}

\details{
The lower bounds are based on the standard methods of convex analysis. See the reference paper at \url{http://www.iam.fmph.uniba.sk/design/} for mathematical details.
}

\value{
A lower bound on the D-, A-, I-, c-, or C-efficiency of \code{w} in the class of all approximate designs of the same size as \code{w} at the set of candidate regressors given by \code{Fx}.
}


\author{
Radoslav Harman, Lenka Filova
}

\note{
The design \code{w} should have a non-singular information matrix. Occasionally, the lower bound is very conservative. The exact value of the efficiency of \code{w} is the ratio of the criterion value of \code{w} and the criterion value of the optimal design.
}

\seealso{
\code{\link{varfun}}, \code{\link{dirder}}
}

\examples{
# A lower bound on the D-efficiencies of the uniform designs
# for the quadratic regression on a line grid

Fx <- Fx_cube(~x1 + I(x1^2), n.levels = 101)
effbound(Fx, rep(1/101, 101))

# The precise value of the D-efficiency
# requires computing the D-optimal design:

w.opt <- od_REX(Fx)$w.best
optcrit(Fx, rep(1/101, 101)) / optcrit(Fx, w.opt)

\dontrun{
# Let us do this for polynomial regressions of various degrees:

n <- 101; d.max <- 10; x <- seq(-1, 1, length = n)
effs <- matrix(0, ncol = 2, nrow = d.max)
Fx <- matrix(1, ncol = 1, nrow = n)
for(d in 1:d.max) {
  Fx <- cbind(Fx, x^d)
  effs[d, 1] <- effbound(Fx, rep(1/n, n))
  w.opt <- od_REX(Fx)$w.best
  effs[d, 2] <- optcrit(Fx, rep(1/n, n)) / optcrit(Fx, w.opt)
}
print(effs)

# We see that the lower bound becomes more and more conservative
# compared to the real efficiency which actually increases with d.


# Compute a D-optimal design for the main effects model
# on a random subsample of a 6D cube

n <- 1000000; m <- 6
Fx <- cbind(1, matrix(runif(n*m), ncol = m))
w <- od_REX(Fx, eff = 0.99)$w.best
Fx <- od_DEL(Fx, w)$Fx.keep
w <- od_REX(Fx)$w.best

# Now we will compute a lower bound on efficiency of such design
# on the entire (continuous) cube:
Fx <- rbind(Fx, Fx_cube(~x1 + x2 + x3 + x4 + x5 + x6, lower = rep(0, 6)))
w <- c(w, rep(0, 2^6))
effbound(Fx, w)

# The real D-efficiency of w on the entire cube is
optcrit(Fx, w)/od_REX(Fx)$Phi.best
}
}

