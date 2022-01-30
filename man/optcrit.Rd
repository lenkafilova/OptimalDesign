\name{optcrit}
\alias{optcrit}

\title{
Criterion value of a design
}

\description{
Computes the criterion value of a design \code{w} in the model determined by the matrix \code{Fx} of all regressors.
}

\usage{
optcrit(Fx, w, crit="D", h=NULL, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} (where \code{m>=2}, \code{m<=n}) matrix containing all candidate regressors (as rows), i.e., \code{n} is the number of candidate design points, and \code{m} (where \code{m>=2}) is the number of parameters.
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
The package works with optimality criteria as information functions, i.e., the criteria are concave, positive homogeneous and upper semicontinuous on the set of all non-negative definite matrices. The criteria are normalized such that they assign the value of 1 to any design with information matrix equal to the identity matrix.

}

\value{
A non-negative number corresponding to the criterion value.
}


\author{
Radoslav Harman, Lenka Filova
}

\note{
Since the criteria are positive homogeneous, the relative efficiency of two designs is just the ratio of their criterion values.
}

\seealso{
\code{\link{infmat}}
}

\examples{
# The Fx matrix for the spring balance weighing model with 6 weighed items.
Fx <- Fx_cube(~x1 + x2 + x3 + x4 + x5 + x6 - 1, lower = rep(0, 6), n.levels = rep(2, 6))

# Criteria of the design of size 15 that weighs each pair of items exactly once.
w2 <- rep(0, 64); w2[apply(Fx, 1, sum) == 2] <- 1
optcrit(Fx, w2, crit = "D")
optcrit(Fx, w2, crit = "A")
optcrit(Fx, w2, crit = "I")

# Criteria for the design of size 15 that weighs each quadruple of items exactly once.
w4 <- rep(0, 64); w4[apply(Fx, 1, sum) == 4] <- 1
optcrit(Fx, w4, crit = "D")
optcrit(Fx, w4, crit = "A")
optcrit(Fx, w4, crit = "I")
}
