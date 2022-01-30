\name{infmat}
\alias{infmat}

\title{
Information matrix of a design
}

\description{
Computes the information matrix of a design \code{w} in the model determined by the matrix \code{Fx} of candidate regressors.
}

\usage{
infmat(Fx, w, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} matrix of candidate regressors (as rows), where \code{n} is the number of candidate design points and \code{m} (where \code{m>=2}, \code{m<=n}) is the number of parameters.
}
\item{w}{
a non-negative vector of length \code{n} representing the design.
}
\item{echo}{
Print the call of the function?
}
}


\value{
The information matrix of the design \code{w} in the model with all candidate regresors given by the rows of \code{Fx}.
}

\author{
Radoslav Harman, Lenka Filova
}

\note{
The information matrix is standardized, i.e., it assumes that the variance of the errors is 1.
}

\seealso{
\code{\link{optcrit}}
}

\examples{

# Compute its information matrix for the design that is
# uniform on all the points with at most two levels equal to 1
# in the main effects model with 2 factors.

Fx <- Fx_cube(~x1 + x2 + x3 + x4 + x5, lower = rep(0, 5))
w <- rep(0, 2^5)
for (i in 1:(2^5)) if (sum(Fx[i, 2:6]) <= 2) w[i] <- 1
print(M <- infmat(Fx, w))

\dontrun{
# Visualize the correlation matrix of the parameter estimators

V <- solve(M); Y <- diag(1/sqrt(diag(V)))
library(corrplot); corrplot(Y \%*\% V \%*\% Y)
}
}

