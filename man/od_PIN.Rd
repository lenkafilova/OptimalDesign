\name{od_PIN}
\alias{od_PIN}

\title{
Efficient saturated exact design
}

\description{
Use a fast greedy method to compute an efficient saturated subset (saturated exact design).
}

\usage{
od_PIN(Fx, alg.PIN="KYM", echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} (where \code{m>=2}, \code{m<=n}) matrix containing all candidate regressors (as rows), i.e., \code{n} is the number of candidate design points, and \code{m} is the number of parameters.
}
\item{alg.PIN}{
the method used (either \code{"KYM"} or \code{"GKM"}). KYM is randomized, faster but provides somewhat less efficient subsets/designs. GKM is deterministic, slower, but tends to give more efficient subsets/designs.
}
\item{echo}{
Print the call of the function?
}
}

\details{
The function is developed with the criterion of D-optimality in mind, but it also gives reasonably efficient subset/designs with respect to other criteria. The main purpose of \code{od_PIN} is to initialize algorithms for computing optimal approximate and exact designs. It can also be used to verify whether a model, represented by a matrix \code{Fx} of candidate regressors, permits a non-singular design.
}

\value{
  Output is the list with components:
  \item{call}{the call of the function}
  \item{w.pin}{the resulting exact design}
  \item{supp}{the indices of the support of w.pin}
  \item{M.pin}{the information matrix of w.pin}
  \item{Phi.D}{the D-criterion value of w.pin}
  \item{t.act}{the actual time of the computation}
}

\references{
Harman R, Rosa S (2019): On greedy heuristics for computing D-efficient saturated subsets, (submitted to Operations Research Letters), \url{https://arxiv.org/abs/1905.07647}
}

\author{
Radoslav Harman, Samuel Rosa, Lenka Filova
}

\examples{
# Compute a saturated subset of a random Fx
Fx <- matrix(rnorm(10000), ncol = 5)
w.KYM <- od_PIN(Fx)$w.pin
w.GKM <- od_PIN(Fx, alg.PIN = "GKM")$w.pin
w.REX <- 5*od_REX(Fx)$w.best
optcrit(Fx, w.KYM)
optcrit(Fx, w.GKM)
optcrit(Fx, w.REX)

}

