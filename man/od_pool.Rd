\name{od_pool}
\alias{od_pool}

\title{
Pool of a vector
}

\description{
A function \code{pool.fun} is applied to all the elements of a vector \code{val} that appear within the groups formed by identical rows of a matrix \code{X}. 
}

\usage{
od_pool(X, val=NULL, pool.fun="sum", echo=TRUE)
}

\arguments{
\item{X}{
the \code{n} times \code{k} matrix of real values.
}
\item{val}{
a real vector of length \code{n}. 
}
\item{pool.fun}{
a string denoting the function to be applied to the subgroups of elements of \code{val} corresponding to the identical rows of \code{X}. Possible values are  \code{"sum"}, \code{"min"}, \code{"max"}, \code{"mean"}, \code{"median"} and \code{"0"}. 
}
\item{echo}{
Print the call of the function?
}
}

\details{
This function is useful for plotting (and understanding) of designs of experiments with more factors than the dimension of the plot.  
}

\value{
A list with components:
  \item{call}{the call of the function}
  \item{X.unique}{the matrix of unique rows of \code{X}}
  \item{val.pooled}{the vector of the length \code{nrows(X.unique)} containing the values of \code{val} pooled using \code{pool.fun}} 
}


\author{
Radoslav Harman, Lenka Filova
}

\note{
The function performs a non-trivial operation only if some of the rows of \code{X} are identical.
}

\seealso{
\code{\link{od_plot}}, \code{\link{od_print}}
}

\examples{
v1 <- c(1, 2, 3); v2 <- c(2, 4, 6); v3 <- c(2, 5, 3)
X <- rbind(v1, v1, v1, v1, v2, v3, v2, v3, v3)
val <- c(1, 2, 7, 9, 5, 8, 4, 3, 6)
od_pool(X, val, "sum")

# The result $val.pooled is a vector with components:
# 19 (=1+2+7+9) because the first 4 rows of X are identical
# 9 (=5+4) because the 5th and the 7th rows of X are identical
# 17 (=8+3+6) because the 6th, the 8th and the 9th rows of X are identical
}


