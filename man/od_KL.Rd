\name{od_KL}
\alias{od_KL}

\title{
The KL exchange algorithm for efficient exact designs
}

\description{
Computes an optimal or near-optimal exact design of experiments under the standard (size) constraint on the size of the experiment.
}

\usage{
od_KL(Fx, N, bin=FALSE, Phi.app=NULL, crit="D", h=NULL, w1=NULL, K=NULL,
      L=NULL, rest.max=Inf, t.max=120, echo=TRUE, track=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} (where \code{m>=2}, \code{m<=n}) matrix containing all candidate regressors (as rows), i.e., \code{n} is the number of candidate design points, and \code{m} (where \code{m>=2}) is the number of parameters.
}
\item{N}{
the size of the experiment (i.e., the required number of trials).
}
\item{bin}{
Should each design point be used at most once?
}
\item{Phi.app}{
the optimal value of the corresponding approximate (relaxed) problem. If \code{Phi.app = NULL}, the value is pre-computed using \code{\link{od_REX}}.
}
\item{crit}{
the optimality criterion. Possible values are \code{"D"}, \code{"A"}, \code{"I"}, \code{"C"}.
}
\item{h}{
a non-zero vector of length \code{m} corresponding to the coefficients of the linear parameter combination of interest. If \code{crit} is not \code{"C"} nor \code{"c"} then \code{h} is ignored. If \code{crit} is \code{"C"} or \code{"c"} and \code{h=NULL} then \code{h} is assumed to be \code{c(0,...,0,1)}.
}
\item{w1}{
the initial design; it must have a non-singular information matrix and the size \code{sum(w1)} of \code{w1} must be \code{N}. The default option \code{w1 = NULL} prompts the algorithm to generate its own initial design using \code{\link{od_PIN}}.
}
\item{K, L}{
integer numbers (or \code{Inf}) representing parameters of the method. Various combinations of \code{K} and \code{L} lead to specific variants of the exchange method. If \code{K = NULL} or \code{L = NULL}, the algorithm automatically chooses appropriate values.
}
\item{rest.max}{
the limit on the number of restarts of the method.
}
\item{t.max}{
the time limit for the computation.
}
\item{echo}{
Print the call of the function?
}
\item{track}{
Display the progress of the computation?
}
}

\details{
This implementation of the KL algorithm is loosely based on the ideas described in Atkinson et al. (2007); see the references.

The tuning parameter \code{K} is the (upper bound on the) number of "least promising" support points of the current design, for which exchanges are attempted. The tuning parameter \code{L} is the (upper bound on the) number of "most promising" candidate design points for which exchanges are attempted.

The implemented method is greedy in the sense that each improving exchange is immediately executed. If the algorithm stops in a local optimum before the allotted time elapsed, the computation is restarted with a random initial design (independent of \code{w1}). The final result is the best design found within all restarts.

The performance of the function depends on the problem, on the chosen parameters, and on the hardware used, but in most cases the function can compute a nearly-optimal exact design for a problem with a ten thousands design points within seconds of computing time. Because this is only a heuristic, we advise the user to verify the quality of the resulting design by comparing it to the result of an alternative method (such as \code{\link{od_RC}}).
}

\value{
  Output is the list with components:
  \item{call}{the call of the function}
  \item{w.best}{the best exact design found by the algorithm}
  \item{supp}{the indices of the support of w.best}
  \item{w.supp}{the weights of w.best on the support}
  \item{M.best}{the information matrix of w.best}
  \item{Phi.best}{the criterion value of w.best}
  \item{eff.best}{a lower bound on the eff of w.best with respect to \code{Phi.app}}
  \item{n.rest}{number of restarts performed}
  \item{n.ex}{the total number of exchanges performed}
  \item{t.act}{the actual time of the computation}
}

\references{
Atkinson AC, Donev AN, Tobias RD (2007): Optimum experimental designs, with SAS. Vol. 34. Oxford: Oxford University Press.
}

\author{
  Radoslav Harman, Lenka Filova
}


\seealso{
\code{\link{od_RC}}, \code{\link{od_AQUA}}, \code{\link{od_MISOCP}}
}

\examples{
\dontrun{
# Compute a D-efficient exact design of size 27 on a unit square
# for the full quadratic model with 2 discretized factors

form.q <- ~x1 + x2 + I(x1^2) + I(x2^2) + I(x1*x2)
Fx <- Fx_cube(form.q, n.levels = c(101, 101))
w <- od_KL(Fx, 13, t.max = 8)$w.best
od_plot(Fx, w, Fx[, 2:3])
od_print(Fx, w)


# Compute an I-efficient exact design of size 100 without replications
# on a discretized L1 ball for the full quadratic model with 3 factors

form.q <- ~x1 + x2 + x3 + I(x1^2) + I(x2^2) + I(x3^2) + I(x1*x2) + I(x1*x3) + I(x2*x3)
Fx <- Fx_cube(form.q, n.levels = c(21, 21, 21))
remove <- (1:nrow(Fx))[apply(abs(Fx[, 2:4]), 1, sum) > 1 + 1e-9]
Fx <- Fx[-remove, ]
w <- od_KL(Fx, 100, bin = TRUE, crit = "I", t.max = 3)$w.best
od_plot(Fx, w, Fx[, 2:4])

# Compute a D-efficient exact design of size 20 on a 4D cube
# for the full quadratic model with 4 continuous factors

# We can begin with a crude discretization and compute
# an initial (already good) exact design using the KL algorithm

form.q <- ~x1 + x2 + x3 + x4 + I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) +
           I(x1*x2) + I(x1*x3) + I(x1*x4) + I(x2*x3) + I(x2*x4) + I(x3*x4)
Fx <- Fx_cube(form.q, n.levels = rep(11, 4))
w <- od_KL(Fx, 20, t.max = 10)$w.best
od_print(Fx, w)$design[, c(2:5, 16)]
print(paste("D-criterion value:", optcrit(Fx, w)))

# Now we can fine-tune the positions of the design points
# using any general-purpose continuous optimization method

F <- Fx[rep(1:nrow(Fx), w), ]
f <- function(x) {c(1, x, x^2, x[1]*x[2], x[1]*x[3], x[1]*x[4],
                               x[2]*x[3], x[2]*x[4], x[3]*x[4])}
obj <- function(x, M.red) {-log(det(M.red + f(x) \%*\% t(f(x))))}
for (i in 1:10)
  for (j in 1:20) {
    F[j, ] <- f(optim(F[j, 2:5], obj, M.red = t(F[-j, ]) \%*\% F[-j, ],
              method = "L-BFGS-B", lower = rep(-1, 3), upper = rep(1, 3))$par)
  }

tune <- od_pool(round(F, 4), rep(1, 20))
Fx.tune <- tune$X.unique; w.tune <- tune$val.pooled
od_print(Fx.tune, w.tune)$design[, c(2:5, 16)]
print(paste("D-criterion value:", optcrit(Fx.tune, w.tune)))
}
}
