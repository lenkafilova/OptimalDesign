\name{Fx_cube}
\alias{Fx_cube}

\title{
Matrix of candidate regressors for a model on a cuboid grid
}

\description{
Creates the matrix of all candidate regressors for a factor regression model on a cuboid grid (up to 9 factors).
}

\usage{
Fx_cube(formula, lower=NULL, upper=NULL, n.levels=NULL, echo=TRUE)
}

\arguments{
\item{formula}{
the formula of the model. The rules for creating the formula are standard for R but: 1) the formula must not contain the dependent variable (it is one-sided); 2) the \code{d} factors (variables) must be labeled \code{x1,x2,x3,}...
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


\value{
The \code{n} times \code{m} matrix of all candidate regressors for a factor regression model on a cuboid grid. The rows of \code{Fx} are the regressors \code{f(x)} for all candidate design points \code{x}.
}

\author{
Radoslav Harman, Lenka Filova
}

\note{
Note that \code{Fx} is \emph{not} the design matrix (which is also sometimes called the regression matrix, or the model matrix). The design matrix depends on \code{Fx} as well as on the exact experimental design \code{w}. For this package, an exact experimental design is formalized as the vector of non-negative integer values corresponding to the replication of trials (observations) in individual design points. Thus, if \code{Fx} is the matrix of all candidate regressors and \code{w} is the exact design then \code{Fx[rep(1:nrow(Fx), w),]} is the actual design matrix for the experiment.
}

\seealso{
\code{\link{Fx_simplex}, \link{Fx_blocks}, \link{Fx_glm}, \link{Fx_survival}, \link{Fx_dose}}
}

\examples{
\dontrun{
# The Fx for the cubic model on a discretized interval
Fx <- Fx_cube(~x1 + I(x1^2) + I(x1^3), lower=0, upper=2, n.levels=101)

# The D-optimal design of size 20
w <- od_KL(Fx, 20, t.max=5)$w.best
od_plot(Fx, w, Fx[, 2])


# The Fx for the full quadratic response surface model on a non-convex region
Fx <- Fx_cube(~x1 + x2 + I(x1^2) + I(x2^2) + I(x1*x2), n.levels=c(51, 51))
keep <- rep(TRUE, nrow(Fx))
for(i in 1:nrow(Fx)) if(prod(abs(Fx[i, 2:3])) > 0.2) keep[i] <- FALSE
Fx <- Fx[keep, ]

# The D-optimal design of size 29 without replications
w <- od_KL(Fx, 29, bin=TRUE, t.max=5)$w.best
od_plot(Fx, w, Fx[, 2:3])


# The Fx for the chemical weighing with 3 items and a bias term
Fx <- Fx_cube(~x1 + x2 + x3, n.levels=c(3, 3, 3))

# The D-optimal design of size 12
w <- od_KL(Fx, 12, t.max=2)$w.best
od_plot(Fx, w, Fx[, 2:4])
}
}
