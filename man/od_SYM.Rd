\name{od_SYM}
\alias{od_SYM}

\title{
Symmetrization of an approximate design
}

\description{
Attempts to "symmetrize" an approximate design \code{w} by minimizing its norm while keeping its information matrix.
}

\usage{
od_SYM(Fx, w, b1=NULL, A1=NULL, b2=NULL, A2=NULL, b3=NULL, A3=NULL, w0=NULL,
      crit="D", h=NULL, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} (where \code{m>=2}, \code{m<=n}) matrix containing all candidate regressors (as rows), i.e., \code{n} is the number of candidate design points, and \code{m} is the number of parameters
}
\item{w}{
a non-negative vector of length \code{n} representing the design
}
\item{b1,A1, b2,A2, b3,A3}{
the real vectors and matrices that define the constraints on permissible designs \code{w} as follows: \code{A1 \%*\% w <= b1}, \code{A2 \%*\% w >= b2}, \code{A3 \%*\% w == b3}. Each of the arguments can be \code{NULL}, but at least one of \code{b1}, \code{b2}, \code{b3} must be non-\code{NULL}. If some \code{bi} is non-\code{NULL} and \code{Ai} is \code{NULL}, then \code{Ai} is set to be \code{matrix(1, nrow =1, ncol = n)}.
}
\item{w0}{
a non-negative vector of length \code{n} representing the design to be augmented (i.e., the function adds the constraint \code{w >= w0} for permissible designs \code{w}). This argument can also be \code{NULL}; in that case, \code{w0} is set to the vector of zeros.
}
\item{crit}{
the optimality criterion. Possible values are \code{"D"}, \code{"A"}, \code{"I"}, \code{"C"}, \code{"c"}.
}
\item{h}{
a non-zero vector of length \code{m} corresponding to the coefficients of the linear parameter combination of interest. If \code{crit} is not \code{"C"} nor \code{"c"} then \code{h} is ignored. If \code{crit} is \code{"C"} or \code{"c"} and \code{h=NULL} then \code{h} is assumed to be \code{c(0,...,0,1)}.
}
\item{echo}{
Print the call of the function?
}
}

\details{
For some models, the optimum approximate design is not unique (although the optimum information matrix usually \emph{is} unique). This function uses one optimal approximate design to produce an optimal approximate design with a minimum Euclidean norm, which is unique and usually more "symmetric".
}

\value{
  A list with the following components:
   \item{call}{The call of the function}
   \item{w.sym}{The resulting "symmetrized" approximate design}
}

\references{
Harman R, Filova L, Richtarik P (2019). A randomized exchange algorithm for computing optimal approximate designs of experiments. Journal of the American Statistical Association, 1-30. (Subsection 5.1)
}

\author{
Radoslav Harman, Lenka Filova
}

\examples{
# Compute a D-optimal approximate design using the randomized method REX.
# Visualize both the design obtained by REX and its symmetrized version.

form.q <- ~x1 + x2 + x3 + I(x1^2) + I(x2^2) + I(x3^2) + I(x1*x2) + I(x1*x3) + I(x2*x3)
Fx <- Fx_cube(form.q, n.levels = c(5, 5, 5))
w.app <- od_REX(Fx)$w.best
od_plot(Fx, w.app, X=Fx[, 2:3])
w.app.sym <- od_SYM(Fx, w.app, b3 = 1)$w.sym
od_plot(Fx, w.app.sym, X=Fx[, 2:3])
}

