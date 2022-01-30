\name{od_MISOCP}
\alias{od_MISOCP}

\title{
  Optimal exact design using mixed integer second-order cone programming
}

\description{
Computes an optimal or nearly-optimal approximate or exact experimental design using mixed integer second-order cone programming.
}

\usage{
od_MISOCP(Fx, b1=NULL, A1=NULL, b2=NULL, A2=NULL, b3=NULL, A3=NULL, w0=NULL,
          bin=FALSE, type="exact", crit="D", h=NULL, gap=NULL,
          t.max=120, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} (where \code{m>=2}, \code{m<=n}) matrix containing all candidate regressors (as rows), i.e., \code{n} is the number of candidate design points, and \code{m} (where \code{m>=2}) is the number of parameters
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
\item{type}{
the type of the design. Permissible values are \code{"approximate"} and \code{"exact"}.
}
\item{crit}{
the optimality criterion. Possible values are \code{"D"}, \code{"A"}, \code{"I"}, \code{"C"}, \code{"c"}.
}
\item{h}{
a non-zero vector of length \code{m} corresponding to the coefficients of the linear parameter combination of interest. If \code{crit} is not \code{"C"} nor \code{"c"} then \code{h} is ignored. If \code{crit} is \code{"C"} or \code{"c"} and \code{h=NULL} then \code{h} is assumed to be \code{c(0,...,0,1)}.
}
\item{gap}{
the gap for the MISOCP solver to stop the computation. If \code{NULL}, the default gap is used.  Setting \code{gap=0} and \code{t.max=Inf} will ultimately provide the optimal exact design, but the computation may be extremely time consuming.
}
\item{t.max}{
the time limit for the computation.
}
\item{echo}{
Print the call of the function?
}
}

\details{
  At least one of \code{b1}, \code{b2}, \code{b3} must be non-\code{NULL}.
  If \code{bi} is non-\code{NULL} and \code{Ai} is \code{NULL} for some \code{i} then \code{Ai} is set to be the vector of ones. If \code{bi} is \code{NULL} for some \code{i} then \code{Ai} is ignored.
}

\value{
A list with the following components:
\item{call}{the call of the function}
\item{w.best}{the permissible design found, or \code{NULL}. The value \code{NULL} indicates a failed computation}
\item{supp}{the indices of the support of \code{w.best}}
\item{w.supp}{the weights of \code{w.best} on the support}
\item{M.best}{the information matrix of \code{w.best} or \code{NULL} if \code{w.best} is \code{NULL}}
\item{Phi.best}{the value of the criterion of optimality of the design \code{w.best}. If \code{w.best} has a singular information matrix or if the computation fails, the value of \code{Phi.best} is \code{0}}
\item{status}{the status variable of the gurobi optimization procedure; see the gurobi solver documentation for details}
\item{t.act}{the actual time of the computation}
}

\references{
  Sagnol G, Harman R (2015): Computing exact D-optimal designs by mixed integer second order cone programming. The Annals of Statistics, Volume 43, Number 5, pp. 2198-2224.
}

\author{
  Radoslav Harman, Lenka Filova
}


\seealso{
  \code{\link{od_KL}, \link{od_RC}, \link{od_AQUA}}
}

\examples{
\dontrun{
# Compute an A-optimal block size two design
# for 6 treatments and 9 blocks

Fx <- Fx_blocks(6)
w <- od_MISOCP(Fx, b3 = 9, crit = "A", bin = TRUE)$w.best
des <- combn(6, 2)[, as.logical(w)]
print(des)

library(igraph)
grp <- graph_(t(des), from_edgelist(directed = FALSE))
plot(grp, layout=layout_with_graphopt)

# Compute a symmetrized D-optimal approximate design
# for the full quadratic model on a square grid
# with uniform marginal constraints

Fx <- Fx_cube(~x1 + x2 + I(x1^2) + I(x2^2) + I(x1*x2), n.levels = c(21, 21))
A3 <- matrix(0, nrow = 21, ncol = 21^2)
for(i in 1:21) A3[i, (i*21 - 20):(i*21)] <- 1
w <- od_MISOCP(Fx, b3 = rep(1, 21), A3 = A3, crit = "D", type = "approximate")$w.best
w.sym <- od_SYM(Fx, w, b3 = rep(1, 21), A3 = A3)$w.sym
od_plot(Fx, w.sym, Fx[, 2:3], dd.size = 2)
}
}

