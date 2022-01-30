\name{dirder}
\alias{dirder}

\title{
Vector of directional derivatives
}

\description{
Computes the vector of derivatives at a normalized approximate design \code{w} of length \code{n} in the directions of singular designs e_i, where i ranges from \code{1} to \code{n}.
}

\usage{
dirder(Fx, w, crit="D", h=NULL, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} matrix of candidate regressors (as rows), where \code{n} is the number of candidate design points and \code{m} (where \code{m>=2}, \code{m<=n}) is the number of parameters.
}
\item{w}{
a non-negative vector of length \code{n} representing the design. It is normalized prior to the computation of the directional derivatives.
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

\details{The i-th directional derivative measures the increase of the criterion value provided that we infinitesimally increase the i-th design weight (and decrease other weights by the same proportion). For a concave optimality criterion, an approximate design is optimal in the class of all normalized approximate designs if and only if all its directional derivatives are non-positive. This statement can be rewritten to the form of the so-called equivalence theorem. See the reference paper at \url{http://www.iam.fmph.uniba.sk/design/} for mathematical details.
}

\value{
The vector of directional derivatives of the chosen criterion at \code{w/sum(w)} in the direction of the singular designs e_i, where i ranges from \code{1} to \code{n}.
}


\author{
Radoslav Harman, Lenka Filova
}

\note{
The design \code{w} should have a non-singular information matrix.
}

\seealso{
\code{\link{effbound}, \link{varfun}}
}

\examples{
\dontrun{
# The directional derivatives of the D-optimal approximate design
# for a cubic regression on a square grid.

form.cube <- ~x1 + x2 + I(x1^2) + I(x2^2) + I(x1*x2) +
              I(x1^3) + I(x1^2*x2) + I(x1*x2^2) + I(x2^3)
Fx <- Fx_cube(form.cube, n.levels = c(101, 101))
w <- od_REX(Fx)$w.best

# Because w is optimal approximate, no directional derivative is positive:
boxplot(dirder(Fx, w))

# The yellow values indicate the directional derivative at each design point:
od_plot(Fx, w, Fx[, 2:3])

# An alternative view is a "projection" of the above plot:
od_plot(Fx, w, Fx[, 2], dd.pool = c("max", "min"))
}
}

