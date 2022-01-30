\name{od_RC}
\alias{od_RC}

\title{Efficient exact design using the RC heuristic}

\description{
Computes an efficient exact design under multiple linear resource constraints using the RC heuristic.
}

\usage{
od_RC(Fx, b, A = NULL, w0 = NULL, bin = FALSE, Phi.app = NULL, crit = "D",
      h=NULL, w1 = NULL, rest.max = Inf, t.max = 120,
      echo = TRUE, track=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} (where \code{m>=2}, \code{m<=n}) matrix containing all candidate regressors (as rows), i.e., \code{n} is the number of candidate design points, and \code{m} is the number of parameters.
}
\item{b, A}{
the vector of length \code{k} with positive real components and the \code{k} times \code{n} matrix of non-negative reals numbers. Each column of \code{A} must have at least one strictly positive element. The linear constraints \code{A\%*\%w<=b, w0<=w} define the set of permissible designs \code{w} (where \code{w0} is a described below.) The argument \code{A} can also be \code{NULL}; in that case \code{b} must be a positive number and \code{A} is set to the \code{1} times \code{n} matrix of ones.
   }
\item{w0}{
a non-negative vector of length \code{n} representing the design to be augmented (i.e., the function adds the constraint \code{w >= w0} for permissible designs \code{w}). This argument can also be \code{NULL}; in that case, \code{w0} is set to the vector of zeros.
}
\item{bin}{
Should each design point be used at most once?
}
\item{Phi.app}{
the optimal value of the corresponding approximate (relaxed) problem. If \code{Phi.app = NULL}, a very conservative upper bound on \code{Phi.app} is pre-computed.
}
\item{crit}{
the optimality criterion. Possible values are \code{"D"}, \code{"A"}, \code{"I"}, \code{"C"}.
}
\item{h}{
a non-zero vector of length \code{m} corresponding to the coefficients of the linear parameter combination of interest. If \code{crit} is not \code{"C"} then \code{h} is ignored. If \code{crit} is \code{"C"} and \code{h=NULL} then \code{h} is assumed to be \code{c(0,...,0,1)}.
}
\item{w1}{
an \code{n} times \code{1} nonnegative vector that represents the initial design. The design \code{w1} must satisfy \code{w0<=w1} and \code{A*w1<=b}. The argument \code{w1} can also be \code{NULL}; in that case the procedure sets \code{w1} to be \code{w0}.
}
\item{rest.max}{
the maximum allowed number of restarts of the method.
}
\item{t.max}{
the time limit for the computation.
}
\item{echo}{
Print the call of the function?
}
\item{track}{
Trace the computation?
}
}

\details{
This is an implementation of the algorithm proposed by Harman et al. (2016); see the references. The inequalities \code{A\%*\%w<=b}, \code{w0<=w} with the specific properties mentioned above, form the so-called resource constraints. They encompass many practical restrictions on the design, and lead to a bounded set of feasible solutions.

The information matrix of \code{w1} should preferably have the reciprocal condition number of at least \code{1e-5}. Note that the floor of an optimal approximate design (computed for instance using \code{\link{od_MISOCP}}) is often a good initial design. Alternatively, the initial design can be the result of another optimal design procedure, such as \code{\link{od_AQUA}}. Even if no initial design is provided, the model should be non-singular in the sense that there \emph{exists} an exact design \code{w} with a well conditioned information matrix, satisfying all constraints. If this requirement is not satisfied, the computation may fail, or it may produce a deficient design.

The procedure always returns a permissible design, but in some cases, especially if \code{t.max} is too small, the resulting design can be inefficient. The performance depends on the problem and on the hardware used, but in most cases the function can compute a nearly-optimal exact design for a problem with a few hundreds design points and tens of constraints within minutes of computing time. Because this is a heuristic method, we advise the user to verify the quality of the resulting design by comparing it to the result of an alternative method (such as \code{\link{od_AQUA}} and \code{\link{od_MISOCP}}) and/or by computing its efficiency relative to the corresponding optimal approximate design.

In the very special (but frequently used) case of the single constraint on the experimental size, it is generally more efficient to use the function \code{\link{od_KL}}.
}

\value{
  A list with the following components:
   \item{call}{The call of the function.}
   \item{w.best}{The resulting exact design.}
   \item{supp}{The indices of the support of \code{w.best}.}
   \item{w.supp}{The weights of \code{w.best} on the support.}
   \item{M.best}{The information matrix of \code{w.best}.}
   \item{Phi.best}{The criterion value of \code{w.best}.}
   \item{eff.best}{A lower bound on the efficiency of \code{w.best} with respect to the optimal approximate design.}
   \item{n.rest}{The number of restarts performed.}
   \item{t.act}{The actual time of the computation.}
}

\references{
  Harman R, Bachrata A, Filova L (2016): Heuristic construction of exact experimental designs under multiple resource constraints, Applied Stochastic Models in Business and Industry, Volume 32, pp. 3-17
}

\author{
  Radoslav Harman, Alena Bachrata, Lenka Filova
}

\seealso{
  \code{\link{od_AQUA}, \link{od_MISOCP}, \link{od_KL}}
}

\examples{
\dontrun{
# A D-efficient exact design for a quadratic model with 2 factors
# constrained by the total time and the total cost of the experiment.
# The cost of a single trial in (x1, x2) is 10 + x1 + 2*x2
# The limit on the total cost is 1000
# (we do not know the number of trials in advance)

form.quad <- ~x1 + x2 + I(x1^2) + I(x2^2) + I(x1 * x2)
Fx <- Fx_cube(form.quad, lower = c(0, 0), upper = c(10, 10), n.levels = c(11, 11))
n <- nrow(Fx); A <- matrix(0, nrow = 1, ncol = n)
for(i in 1:n) A[1, i] <- 5 + Fx[i, 2] + 2*Fx[i, 3]
w <- od_RC(Fx, 1000, A, bin = TRUE, t.max = 8)$w.best
od_plot(Fx, w, Fx[, 2:3], dd.size = 3)
}
}
