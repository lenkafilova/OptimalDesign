\name{od_PUK}
\alias{od_PUK}

\title{
Efficient rounding of an approximate design
}

\description{
Compute the classical efficient rounding of a non-normalized approximate design \code{w} such that the resulting exact design has size \code{floor(sum(w))}.
}

\usage{
od_PUK(Fx, w, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} (where \code{m>=2}, \code{m<=n}) matrix containing all candidate regressors (as rows), i.e., \code{n} is the number of candidate design points, and \code{m} (where \code{m>=2}) is the number of parameters
}
\item{w}{
the vector of non-negative real numbers of length \code{n} representing the design
}
\item{echo}{
Print the call of the function?
}

}

\value{
The rounded version of \code{w}
}

\references{
Pukelsheim F, Rieder S (1992) Efficient rounding of approximate designs. Biometrika, 79(4), 763--770.
}

\author{
Radoslav Harman and Samuel Rosa
}


\examples{
# Compute a D-optimal approximate design
# Round it using the efficient rounding to various sizes
# Visualize the designs

Fx <- Fx_cube(~x1 + I(x1^2) + I(x1^3), lower = 0, upper = 1, n.levels = 11)
w.app <- od_REX(Fx)$w.best
Phi.app <- optcrit(Fx, w.app)

w.ex10 <- od_PUK(Fx, 10*w.app)$w.round
w.ex20 <- od_PUK(Fx, 20*w.app)$w.round
w.ex30 <- od_PUK(Fx, 30*w.app)$w.round

par(mfrow = c(2, 2))
od_plot(Fx, w.app, main.lab = "Approximate")
od_plot(Fx, w.ex10, main.lab = paste("N=10, Eff:", round(optcrit(Fx, w.ex10)/Phi.app/10, 4)))
od_plot(Fx, w.ex20, main.lab = paste("N=20, Eff:", round(optcrit(Fx, w.ex20)/Phi.app/20, 4)))
od_plot(Fx, w.ex30, main.lab = paste("N=30, Eff:", round(optcrit(Fx, w.ex30)/Phi.app/30, 4)))
par(mfrow = c(1, 1))
}


