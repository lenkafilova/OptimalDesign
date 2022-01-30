\name{Fx_glm}
\alias{Fx_glm}

\title{
Matrix of candidate regressors for a generalized linear model
}

\description{
Creates the matrix of all candidate regressors for a linearization of a generalized linear model.
}

\usage{
Fx_glm(formula, theta0, glm.model="bin-logit", lower=NULL, upper=NULL,
       n.levels=NULL, echo=TRUE)
}

\arguments{
\item{formula}{
the formula of the linear part of the model. The rules for creating the formula are standard for R but: 1) the formula must not contain the dependent variable (it is one-sided); 2) the \code{d} factors (variables) must be labeled \code{x1,x2,x3,}...
}
\item{theta0}{
the \code{d}-dimensional vector of values of the unknown parameter in which to linearize the model
}
\item{glm.model}{
the type of the generalized linear model. Available models are \code{"bin-logit"}, \code{"bin-probit"},  \code{"bin-cloglog"}, and \code{Poisson-log}.
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
The \code{n} times \code{m} matrix of all candidate regressors of a generalized linear regression model linearized in \code{theta0}.
}

\references{
Atkinson AC, Woods DC (2015). Designs for generalized linear models. Handbook of Design and Analysis of Experiments, 471-514.
}

\author{
Radoslav Harman, Lenka Filova
}


\seealso{
\code{\link{Fx_cube}, \link{Fx_simplex}, \link{Fx_blocks}, \link{Fx_survival}, \link{Fx_dose}}
}

\examples{
# The logistic model with second-order predictors x1, x2 in [-1,1]
# discretized into 21 points and theta0=c(1, 2, 2, -1, -1.5, 1.5)

form.quad <- ~ x1 + x2 + I(x1*x2) + I(x1^2) + I(x2^2)
Fx <- Fx_glm(form.quad, c(1, 2, 2, -1, -1.5, 1.5),
             glm.model="bin-logit", n.levels=c(21,21))

# The locally D-optimal approximate design
w <- od_REX(Fx)$w.best
Fx.lin <- Fx_cube(form.quad, n.levels=c(21,21)) # Just for the plot
od_plot(Fx, w, Fx.lin[, 2:3], dd.size=2)

\dontrun{
#The GLM with Poisson link and 2 linear predictors x1,x2 in [-1,1]
# discretized into 21 points and theta0=c(0,2,2)
Fx <- Fx_glm(~x1+x2, c(0, 2, 2), glm.model="Poisson-log", n.levels=c(21, 21))

# The locally D-optimal exact design of size 50 without replications
w <- od_KL(Fx, 50, bin=TRUE, t.max=5)$w.best
Fx.lin <- Fx_cube(~x1+x2, n.levels=c(21, 21))
od_plot(Fx, w, Fx.lin[, 2:3], w.lim=Inf)
}
}

