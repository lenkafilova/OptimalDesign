\name{Fx_survival}
\alias{Fx_survival}

\title{
Matrix of candidate regressors for a survival model
}

\description{
Creates the matrix of all candidate regressors for a linearization of a proportional hazards survival model.
}

\usage{
Fx_survival(formula, theta0, censor.time, survival.model="phI", lower=NULL,
            upper=NULL, n.levels=NULL, echo=TRUE)
}

\arguments{
\item{formula}{
the formula of the linear part of the model. The rules for creating the formula are standard for R but: 1) the formula must not contain the dependent variable (it is one-sided); 2) the \code{d} factors (variables) must be labeled \code{x1,x2,x3,}...
}
\item{theta0}{
the \code{d}-dimensional vector of values of the unknown parameter in which to linearize the model.
}
\item{censor.time}{
the censoring time, a positive constant.
}
\item{survival.model}{
the type of the survival model, can be either proportional hazards with Type I censoring (\code{"phI"}) or with random censoring (\code{"phrand"}). Both models assume a constant baseline hazard.
}
\item{lower}{
the \code{d}-dimensional vector of the smallest values of factors. If \code{lower=NULL}, the program sets \code{lower <- rep(-1, d)}.
}
\item{upper}{
the \code{d}-dimensional vector of the largest values of factors. If \code{upper=NULL}, the program sets \code{upper <- rep(1, d)}.
}
\item{n.levels}{
the \code{d}-dimensional vector of the numbers of levels of each factor. If \code{n.levels=NULL}, the program sets \code{n.levels <- rep(2, d)}.
}
\item{echo}{
Print the call of the function?
}
}

\details{
For mathematical details, see the referenced paper.
}

\value{
The \code{n} times \code{m} matrix of all candidate regressors of a proportional hazards model linearized in \code{theta0}.
}

\references{
Konstantinou M, Biedermann S, Kimber A (2014). Optimal designs for two-parameter nonlinear models with application to survival models. Statistica Sinica, 24(1), 415-428.
}

\author{
Radoslav Harman, Lenka Filova
}


\seealso{
\code{\link{Fx_cube}, \link{Fx_simplex}, \link{Fx_blocks}, \link{Fx_glm}, \link{Fx_dose}}
}

\examples{
# The proportional hazards model with random censoring
# for three binary explanatory variables x1,x2,x3 without intercept
# censoring time 30 and parameter values theta0=c(1,1,1)
Fx <- Fx_survival(~x1 + x2 + x3 - 1, c(1, 1, 1), 30, "phrand",
      lower = c(0, 0, 0), upper = c(1, 1, 1), n.levels = c(2, 2, 2))

# The locally D-optimal approximate design
w <- od_REX(Fx, crit="D")$w.best
od_print(Fx, w, Fx)

\dontrun{
# The proportional hazards model with random censoring
# for explanatory variables x1,x2,x3 in the range [0,1] discretized into 11 points
# censoring time 30 and parameter values theta0=c(1,1,1)
Fx <- Fx_survival(~x1 + x2 + x3 - 1, c(1, 1, 1), 30, "phrand",
      lower = c(0, 0, 0), upper = c(1, 1, 1), n.levels = c(11, 11, 11))

# The locally A-optimal exact design of size 50 without replications
w <- od_KL(Fx, 50, crit="A", bin=TRUE, t.max=5)$w.best
od_plot(Fx, w, Fx)
}
}

