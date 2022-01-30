\name{Fx_CtoA}
\alias{Fx_CtoA}

\title{
Transformation of candidate regressors for regularized c-optimality
}

\description{
Pre-transforms the matrix of all candidate regressors to the form suitable for computing regularized c-optimal designs via A-optimum algorithms.
}

\usage{
Fx_CtoA(Fx, h=NULL, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} matrix of candidate regressors (as rows), where \code{n} is the number of candidate design points and \code{m} (where \code{m>=2}, \code{m<=n}) is the number of parameters.
}
\item{h}{
a non-zero vector of length \code{m} corresponding to the coefficients of the linear parameter combination of interest. If \code{crit} is not \code{"C"} nor \code{"c"} then \code{h} is ignored. If \code{crit} is \code{"C"} or \code{"c"} and \code{h=NULL} then \code{h} is assumed to be \code{c(0,...,0,1)}.
}
\item{echo}{
Print the call of the function?
}
}

\details{
The standard c-optimal designs are often singular, which may render them unsuitable for practical use. The regularized c-optimality, which we call C-optimality, is an approach to computing designs that are non-singular, but still efficient with respect to the criterion of c-optimality. See \url{http://www.iam.fmph.uniba.sk/design/} for more details.
}

\value{
The \code{n} times \code{m} matrix \code{Fx.trans} of all candidate regressors with the following property: The A-optimal design for the problem defined by \code{Fx.trans} is the same as the regularized c-optimal (i.e., C-optimal) design for the problem defined by \code{Fx}.
}


\author{
Radoslav Harman and Lenka Filova
}


\seealso{
\code{\link{Fx_ItoA}}
}

\examples{
# We will compute a C-efficient (regularized c-optimal) design
# for estimating the mean response in x=1 for a quadratic regression
# using and algorithm for A-optimality.

Fx <- Fx_cube(~x1 + I(x1^2), n.levels=101)
Fx.trans <- Fx_CtoA(Fx, h=c(1, 1, 1))
w <- od_REX(Fx.trans, crit="A")$w.best
od_print(Fx, w, h=c(1, 1, 1))

# Compare the design to the (non-regularized) c-optimal design
w.crisp <- od_REX(Fx, crit="c", h=c(1, 1, 1))$w.best
od_print(Fx, w.crisp, h=c(1, 1, 1))

# The c-efficiency of the C-optimal design is about 0.68
# The D-efficiency of the c-optimal design is 0
# The D-efficiency of the C-optimal design is a very decent
optcrit(Fx, w) / od_REX(Fx)$Phi.best
}

