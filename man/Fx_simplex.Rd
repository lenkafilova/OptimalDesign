\name{Fx_simplex}
\alias{Fx_simplex}

\title{
Matrix of candidate regressors for a regression model on a simplex grid
}

\description{
Creates the matrix of all candidate regressors for a mixture regression model on a regular simplex grid (up to 9 factors).
}

\usage{
Fx_simplex(formula, n.levels.mix=NULL, echo=TRUE)
}

\arguments{
\item{formula}{
the formula of the model. The rules for creating the formula are standard for R but: 1) the formula must not contain the dependent variable (it is one-sided); 2) the \code{d} factors (variables) must be labeled \code{x1,x2,x3,}...
}
\item{n.levels.mix}{
the number of levels of each factor (each factor has the same number of levels). If \code{n.levels=NULL}, the program sets \code{n.levels <- 2*d + 1}.
}
\item{echo}{
Print the call of the function?
}
}


\value{
The \code{n} times \code{m} matrix of all candidate regressors of a mixture regression model on a regular simplex grid.
}


\author{
Radoslav Harman, Lenka Filova
}

\note{
Note that \code{Fx} is \emph{not} the design matrix (which is also sometimes called the regression matrix, or the model matrix). The design matrix depends on \code{Fx} as well as on the exact experimental design \code{w}. For this package, an exact experimental design is formalized as the vector of non-negative integer values corresponding to the replication of trials (observations) in individual design points. Thus, if \code{Fx} is the matrix of all candidate regressors and \code{w} is the exact design then \code{Fx[rep(1:nrow(Fx), w),]} is the actual design matrix for the experiment.
}

\seealso{
\code{\link{Fx_cube}, \link{Fx_glm}, \link{Fx_dose}, \link{Fx_survival}, \link{Fx_blocks}}
}

\examples{
\dontrun{
# The Fx of the Scheffe quadratic mixture model
# with 3 mixture components, each with 21 levels.
Fx <- Fx_simplex(~x1 + x2 + x3 + I(x1*x2) + I(x1*x3) + I(x2*x3) - 1, 21)

# The approximate I-optimal design of size 20
# bound by 1 at each design point
w <- od_MISOCP(Fx, b3=20, bin=TRUE, crit="I", type="approximate")$w.best
od_plot(Fx, w, Fx[, 2:3])

# As above, with constraints on the proportions
r <- c(); for (i in 1:nrow(Fx)) if (max(Fx[i, 2:4]) > 0.7) r <- c(r, i)
w <- od_MISOCP(Fx[-r, ], b3=20, bin=TRUE, crit="I", type="approximate")$w.best
od_plot(Fx[-r, ], w, Fx[-r, 2:3])

# Note that one must be careful when choosing a model for a mixture experiment:
# Let us compute the matrix of regressors of the simple linear mixture model
# with 4 mixture components, each with levels {0, 0.5, 1}.

Fx <- Fx_simplex(~x1 + x2 + x3 + x4, 3)

# The model has only 5 parameters and as many as 10 design points,
# but there is no design that guarantees estimability of the parameters.
# This can be shown by evaluating:
det(infmat(Fx, rep(1, 10)))
}
}
