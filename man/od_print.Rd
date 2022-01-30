\name{od_print}
\alias{od_print}

\title{
Compact information about a design
}
\description{
Prints various characteristics of an experimental design
}
\usage{
od_print(Fx, w, X=NULL, h=NULL, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} (where \code{m>=2}, \code{m<=n}) matrix containing all candidate regressors (as rows), i.e., \code{n} is the number of candidate design points, and \code{m} is the number of parameters
}
\item{w}{
the vector of non-negative real numbers of length \code{n} representing the design
}
\item{X}{
an \code{n} times \code{k} matrix of coordinates of design points
}
\item{h}{
a non-zero vector of length \code{m} corresponding to the coefficients of the linear parameter combination of interest. If \code{crit} is not \code{"C"} then \code{h} is ignored. If \code{crit} is \code{"C"} and \code{h=NULL} then \code{h} is assumed to be \code{c(0,...,0,1)}
}
\item{echo}{
Print the call of the function?
}
}


\value{
  Output is a list with components
  \item{call}{the call of the function}
  \item{design}{a matrix with the rows of \code{Fx} corresponding to non-zero design weights and the non-weights themselves}
  \item{M}{the information matrix of \code{w}}
  \item{eigenvalues}{the eigenvalues of \code{M}}
  \item{D.value}{the value of the D-optimality criterion for \code{w}}
  \item{A.value}{the value of the A-optimality criterion for \code{w}}
  \item{I.value}{the value of the I-optimality criterion for \code{w}}
  \item{C.value}{the value of the C-optimality criterion for \code{w}}
  \item{c.value}{the value of the c-optimality criterion for \code{w}}
}


\author{
  Radoslav Harman, Lenka Filova
}


\seealso{
\code{\link{od_plot}}, \code{\link{od_pool}}
}

\examples{
Fx <- Fx_cube(~x1 + I(x1^2), n.levels = 11)
w <- 1:11/sum(1:11)
od_print(Fx, w, Fx[, 2])
}


