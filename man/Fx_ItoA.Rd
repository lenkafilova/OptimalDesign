\name{Fx_ItoA}
\alias{Fx_ItoA}

\title{
Transformation of candidate regressors for I-optimality
}

\description{
Pre-transforms the matrix of all candidate regressors to the form suitable for computing I-optimal designs via A-optimum algorithms.
}

\usage{
Fx_ItoA(Fx, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} matrix of candidate regressors (as rows), where \code{n} is the number of candidate design points and \code{m} (where \code{m>=2}, \code{m<=n}) is the number of parameters.
}
\item{echo}{
Print the call of the function?
}
}

\details{
It is simple to see that the problem of I-optimality is equivalent to the problem of A-optimality for a transformed matrix of candidate regressors. This function performs the transformation. See \url{http://www.iam.fmph.uniba.sk/design/} for more details.
}

\value{
The \code{n} times \code{m} matrix \code{Fx.trans} of all candidate regressors with the following property: The A-optimal design for the problem defined by \code{Fx.trans} is the same as the I-optimal design for the problem defined by \code{Fx}.
}


\author{
Radoslav Harman, Lenka Filova
}

\note{
It is also simple to transform the \emph{weighted} I-optimality to A-optimality; just multiply the rows of \code{Fx} by the squares roots of weights of individual design points and transform the resulting matrix by \code{Fx_ItoA}.
}

\seealso{
\code{\link{Fx_CtoA}}
}

\examples{
\dontrun{
# Compute an I-efficient exact size 20 design without replications
# for the Scheffe mixture model with 4 components
# using the AQUA heuristic for A-optimality.

Fx <- Fx_simplex(~x1 + x2 + x3 + x4 + I(x1*x2) + I(x1*x3) + I(x1*x4) +
                  I(x2*x3) + I(x2*x4) + I(x3*x4) - 1, 11)

w <- od_AQUA(Fx_ItoA(Fx), b3=24, bin=TRUE, crit="I", conic=FALSE)$w.best
od_plot(Fx, w, Fx[, 2:4])
}
}

