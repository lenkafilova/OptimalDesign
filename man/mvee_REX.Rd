\name{mvee_REX}
\alias{mvee_REX}

\title{
Minimum-volume enclosing ellipsoid
}

\description{
Computes the shape matrix \code{H} and the center \code{z} of the minimum-volume ellipsoid enclosing a finite set of data-points.
}

\usage{
mvee_REX(Data, alg.AA="REX", eff=0.999999, it.max=Inf, t.max=60,
         picture=FALSE, echo=TRUE, track=TRUE)
}

\arguments{
\item{Data}{
the \code{n} times \code{d} (where \code{d<n}) matrix containing the \code{d}-dimensional data-vectors as rows.
}
\item{alg.AA}{
the underlying computational method for approximate D-optimal design; possible values are \code{"REX"}, \code{"MUL"} and \code{"VDM"}.
}
\item{eff}{
the minimum required efficiency.
}
\item{it.max}{
a limit on the number of iterations of the underlying D-optimum approximate design algorithm.
}
\item{t.max}{
a limit on the time of computation.
}
\item{picture}{
Should a picture be plotted? (For the picture, the data need to be either two- or three-dimensional.)
}
\item{echo}{
Print the call of the function?
}
\item{track}{
Display the progress of the computation?
}}

\details{
The problem of the minimum-volume data-enclosing ellipsoid (MVEE) is computationally equivalent to the problem of D-optimal approximate design for an artificial problem based on the data. This procedure performs the computation and the proper conversion of the D-optimal approximate design to the MVEE parameters (the center and the shape matrix).
}

\value{
Output is a list with components:
\item{call}{the call of the function}
\item{H}{the shape matrix of the MVEE}
\item{z}{the center of the MVEE}
\item{bpts}{a set containing the boundary points of the MVEE}
\item{vol}{the volume of the MVEE}
\item{eff.best}{the actual precision of the result (1 is the perfect precision)}
\item{t.iter}{the number of iterations of the underlying D-optimum design algorithm}
\item{t.act}{the actual time of the computation}
}

\references{
Harman R, Filova L, Richtarik P (2019). A randomized exchange algorithm for computing optimal approximate designs of experiments. Journal of the American Statistical Association, 1-30.
}

\author{
  Radoslav Harman, Lenka Filova
}

\note{
Note: The affine hull of the rows of \code{X} should be the full space of dimension \code{d}. For the choice of the algorithm, see the comments in \code{\link{od_REX}}.
}

\seealso{
\code{\link{od_REX}}
}

\examples{
# Generate random 1000 points in a 3-dimensional space
# and compute the MVEE

Data <- matrix(rnorm(3000), ncol = 3)
mvee_REX(Data, picture = FALSE)
}


