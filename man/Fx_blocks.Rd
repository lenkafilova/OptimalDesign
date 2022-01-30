\name{Fx_blocks}
\alias{Fx_blocks}

\title{
Matrix of candidate regressors for a block size-two model
}

\description{
Creates the matrix of all candidate regressors of a linear regression model corresponding to the problem of the optimal block size-two design.
}

\usage{
Fx_blocks(n.treats, blocks=NULL, echo=TRUE)
}

\arguments{
\item{n.treats}{
the number of "treatments" in the block experiment.
}
\item{blocks}{
the \code{2} times \code{n} matrix of all permissible blocks (that is, permissible pairings of treatments). If \code{blocks=NULL}, \code{blocks} is set to \code{combn(n.treats, 2)}, which means that all treatment pairings are permissible.
}
\item{echo}{
Print the call of the function?
}
}

\details{
Creates the matrix \code{Fx} of artificial regressors, such that the D- and A-optimal designs for the corresponding artificial LRM are are the same as what is called the D- and A-optimal design in the original block model with blocks of size two.
}

\value{
the \code{n} times \code{m} matrix of all candidate regressors of an auxiliary linear regression model corresponding to the problem of the optimal block size-two design (\code{n} is \code{ncol(blocks)}, \code{m} is \code{n.treats-1}).
}

\references{
Harman R, Filova, L: Computing efficient exact designs of experiments using integer quadratic programming, Computational Statistics and Data Analysis 71 (2014) 1159-1167.

Sagnol G, Harman R: Computing Exact D-optimal designs by mixed integer second-order cone programming, The Annals of Statistics 43 (2015), 2198-2224.
}

\author{
Radoslav Harman, Lenka Filova
}

\note{
This optimal design problem is equivalent to various optimum-subgraph problems, depending on the criterion.
}

\seealso{
\code{\link{Fx_cube}, \link{Fx_simplex}, \link{Fx_glm}, \link{Fx_dose}, \link{Fx_survival}}
}

\examples{
\dontrun{
# Compute a D-efficient block size-two design
# with 15 treatments and 10 blocks of size two

Fx <- Fx_blocks(10)
w <- od_KL(Fx, 15, t.max = 5)$w.best
des <- combn(10, 2)[, as.logical(w)]
print(des)


# We can visualize the design as a graph
library(igraph)
grp <- graph_(t(des), from_edgelist(directed = FALSE))
plot(grp, layout=layout_with_graphopt)
}
}
