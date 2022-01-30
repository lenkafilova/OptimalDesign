\name{od_AQUA}
\alias{od_AQUA}

\title{Efficient exact design using a quadratic approximation}

\description{
Computes an efficient exact design under general linear constraints via a quadratic approximation of the optimality criterion.
}

\usage{
od_AQUA(Fx, b1=NULL, A1=NULL, b2=NULL, A2=NULL, b3=NULL, A3=NULL, w0=NULL,
      bin=FALSE, crit="D", h=NULL, M.anchor=NULL, ver.qa="+", conic=TRUE,
      t.max=120, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} (where \code{m>=2}, \code{m<=n}) matrix containing all candidate regressors (as rows), i.e., \code{n} is the number of design points, and \code{m} (where \code{m>=2}) is the number of parameters
}
\item{b1,A1, b2,A2, b3,A3}{
the real vectors and matrices that define the constraints on permissible designs \code{w} as follows: \code{A1 \%*\% w <= b1}, \code{A2 \%*\% w >= b2}, \code{A3 \%*\% w == b3}. Each of the arguments can be \code{NULL}, but at least one of \code{b1}, \code{b2}, \code{b3} must be non-\code{NULL}. If some \code{bi} is non-\code{NULL} and \code{Ai} is \code{NULL}, then \code{Ai} is set to be \code{matrix(1, nrow =1, ncol = n)}.
}
\item{w0}{
a non-negative vector of length \code{n} representing the design to be augmented (i.e., the function adds the constraint \code{w >= w0} for permissible designs \code{w}). This argument can also be \code{NULL}; in that case, \code{w0} is set to the vector of zeros.
}
\item{bin}{
Should each design point be used at most once?
}
\item{crit}{
The optimality criterion. Possible values are \code{"D"}, \code{"A"}, \code{"I"}, \code{"C"}.
}
\item{h}{
a non-zero vector of length \code{m} corresponding to the coefficients of the linear parameter combination of interest. If \code{crit} is not \code{"C"} then \code{h} is ignored. If \code{crit} is \code{"C"} and \code{h=NULL} then \code{h} is assumed to be \code{c(0,...,0,1)}.
}
\item{M.anchor}{
the \code{m} times \code{m} information matrix of the optimal or nearly-optimal approximate design for the design problem (for the non-normalized version of the problem and including the design constraints). The argument \code{M.anchor} can also be \code{NULL}. In that case the procedure computes \code{M.anchor} using an appropriate approximate design procedure from the package.
}
\item{ver.qa}{
version of the criterion; possible values are \code{"+"} and \cite{"-"}.
}
\item{conic}{
Should the conic reformulation be used?
}
\item{t.max}{
the time limit for the computation.
}
\item{echo}{
Print the call of the function?
}
}

\value{
A list with the following components:
\item{call}{The call of the function.}
\item{w.best}{The permissible design found, or \code{NULL}. The value \code{NULL} indicates a failed computation.}
\item{supp}{The indices of the support of \code{w.best}.}
\item{w.supp}{The weights of \code{w.best} on the support.}
\item{M.best}{The information matrix of \code{w.best} or \code{NULL} if \code{w.best} is \code{NULL}.}
\item{Phi.best}{The value of the criterion of optimality of the design \code{w.best}. If \code{w.best} has a singular information matrix or if the computation fails, the value of \code{Phi.best} is \code{0}.}
\item{status}{The status variable of the gurobi optimization procedure; see the gurobi solver documentation for details.}
\item{t.act}{The actual time of the computation.}
}

\details{
  At least one of \code{b1}, \code{b2}, \code{b3} must be non-\code{NULL}.
  If \code{bi} is non-\code{NULL} and \code{Ai} is \code{NULL} for some \code{i} then \code{Ai} is set to be the vector of ones. If \code{bi} is \code{NULL} for some \code{i} then \code{Ai} is ignored.
}

\references{
  Harman R., Filova L. (2014): Computing efficient exact designs of experiments using integer quadratic programming, Computational Statistics & Data Analysis, Volume 71, pp. 1159-1167

  Filova L., Harman R. (2018). Ascent with Quadratic Assistance for the Construction of Exact Experimental Designs. arXiv preprint arXiv:1801.09124. (Submitted to Computational Statistics)
}

\author{
  Radoslav Harman, Lenka Filova
}

\note{
The function does not support the classical c-optimality, but it includes its regularized version referred to as C-optimality. The computation is generally stable, but it may fail for instance if the model is numerically singular, there is no exact design satisfying the constraints, no permissible exact design was found within the time limit, the set of feasible exact designs is unbounded and so on; see the \code{status} variable for more details. Note, however, that \code{status = "OPTIMAL"} indicates that the auxiliary integer programming problem was completely solved, which for this procedure does not guarantee that the result is a globally optimal design.
}

\seealso{
  \code{\link{od_KL}, \link{od_RC}, \link{od_MISOCP}}
}

\examples{
\dontrun{
# Compute an I-efficient non-replicated exact design of size 51
# for the "special cubic" model with 3 mixture components

# Each factor has 11 levels:
form.sc <- ~x1 + x2 + x3 + I(x1*x2) + I(x1*x3) + I(x2*x3) + I(x1*x2*x3) - 1
Fx <- Fx_simplex(form.sc, 11)
w <- od_AQUA(Fx, b3 = 51, crit = "I", bin = TRUE)$w.best
od_plot(Fx, w, Fx[, 1:3])

# Each factor has 101 levels (memory intensive without the conic trick)
Fx <- Fx_simplex(form.sc, 101)
w <- od_AQUA(Fx, b3 = 51, crit = "I", bin = TRUE, t.max = 10)$w.best
od_plot(Fx, w, Fx[, 1:3])

# Find an A-efficient exact design for the spring balance model
# with 5 items and 10 weighings
Fx <- Fx_cube(~x1 + x2 + x3 + x4 + x5 - 1, lower = rep(0, 5))
w <- od_AQUA(Fx, b3 = 10, crit = "A", t.max = 10)$w.best
od_print(Fx, w)

}
}
