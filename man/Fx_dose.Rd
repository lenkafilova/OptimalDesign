\name{Fx_dose}
\alias{Fx_dose}

\title{
Matrix of candidate regressors for a dose-response model
}

\description{
Creates the matrix of all candidate regressors for a linearization of a dose response model.
}

\usage{
Fx_dose(dose.levels, theta0, dose.model="emax", echo=TRUE)
}

\arguments{
\item{dose.levels}{
the \code{n}-dimensional vector of admissible doses.
}
\item{theta0}{
the 3-dimensional vector of values of the unknown parameter in which to linearize the model.
}
\item{dose.model}{
the type of the dose-response model, possible values are \code{"emax"}, \code{"loglin"}, and \code{"exp"}.
}
\item{echo}{
Print the call of the function?
}
}

\details{
For mathematical details, see the referenced paper.
}

\value{
The \code{n} times \code{3} matrix of all candidate regressors of a dose-response model linearized in \code{theta0}.
}

\references{
Dette H, Kiss C, Bevanda M, Bretz F (2010). Optimal designs for the EMAX, log-linear and exponential models. Biometrika, 97(2), 513-518.
}

\author{
Radoslav Harman, Lenka Filova
}

\seealso{
\code{\link{Fx_cube}, \link{Fx_simplex}, \link{Fx_blocks}, \link{Fx_glm}, \link{Fx_survival}}
}

\examples{
# The loglinear model for the doses 1:150
# Localized at the values of theta0=c(0, 0.0797, 1)
Fx <- Fx_dose(1:150, c(0, 0.0797, 1), dose.model="loglin")

# The locally D-optimal approximate design
w_a <- od_REX(Fx)$w.best
od_plot(Fx, w_a, 1:150)

# The locally D-optimal exact design of size 10
w_e <- od_KL(Fx, 10, t.max=3)$w.best
od_plot(Fx, w_e, 1:150)
}


