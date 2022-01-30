\name{od_plot}
\alias{od_plot}

\title{
Visualization of a design
}

\description{
Visualizes selected aspects of an experimental design
}

\usage{
od_plot(Fx, w, X=NULL, w.pool=c("sum", "0"), w.color="darkblue",
        w.size=1, w.pch=16, w.cex=0.8, w.lim=0.01, crit="D",
        h=NULL, dd.pool=c("max", "mean"), dd.color="orange",
        dd.size=1.5, dd.pch=15, asp = NA, main.lab="",
        y.lab="", return.pools=FALSE, echo=TRUE)
}

\arguments{
\item{Fx}{
the \code{n} times \code{m} (where \code{m>=2}, \code{m<=n}) matrix containing all candidate regressors (as rows), i.e., \code{n} is the number of candidate design points, and \code{m} is the number of parameters
}
\item{w}{
the vector of non-negative real numbers of length \code{n} representing the design
}
\item{X}{
an \code{n} times \code{k} matrix of coordinates of design points
}
\item{w.pool}{
a vector with components from the set \code{"sum"}, \code{"min"}, \code{"max"}, \code{"mean"}, \code{"median"}, \code{"0"} that determines various "pools" of the design weights along the projections defined by the coordinates provided by \code{X}
}
\item{w.color}{
the color string for plotting the design weight (from the standard list of R colors)
}
\item{w.size}{
the size of the characters/balls that represent the non-zero design weights
}
\item{w.pch}{
the numerical code of the characters used to plot the non-zero design weights
}
\item{w.cex}{
the size of the text labels representing the magnitudes of the design weights
}
\item{w.lim}{
a threshold fraction of the total design weight to plot the labels
}
\item{crit}{
the optimality criterion. Possible values are \code{"D"}, \code{"A"}, \code{"I"}, \code{"C"}
}
\item{h}{
a non-zero vector of length \code{m} corresponding to the coefficients of the linear parameter combination of interest. If \code{crit} is not \code{"C"} then \code{h} is ignored. If \code{crit} is \code{"C"} and \code{h=NULL} then \code{h} is assumed to be \code{c(0,...,0,1)}
}
\item{dd.pool}{
a vector with components from the set \code{"sum"}, \code{"min"}, \code{"max"}, \code{"mean"}, \code{"median"}, \code{"0"} that determines various "pools" of the directional derivatives along the projections defined by the coordinates provided by \code{X}
}
\item{dd.color}{
the color string for plotting the directional derivatives (from the standard list of R colors)
}
\item{dd.size}{
the size of the characters that represent the directional derivatives
}
\item{dd.pch}{
the numerical code of the character used to plot the directional derivatives
}
\item{asp}{
the aspect ration of the axes
}
\item{main.lab}{
the main title of the plot
}
\item{y.lab}{
the label of the y axis
}
\item{return.pools}{
Return the pooled values?
}
\item{echo}{
Print the call of the function?
}
}

\details{
This function performs a simple visualization of some aspects of an experimental design. It visualizes (the selected pools of) the design weights and (the selected pools of) the directional derivative. The type of graph depends on the number of columns in \code{X}.
}

\value{
If \code{return.pool} is set to \code{TRUE}, the procedure returns the data used to plot the figure. The data can be used to plot a different figure according to the user's needs.
}


\author{
Radoslav Harman, Lenka Filova
}

\note{
The labels of the axes correspond to the column names of \code{X}. For a large \code{unique(Fx)}, rendering the plot can take a considerable time. Note also that using RStudio, it may be a good idea to open an external graphical window (using the command \code{windows()}) before running \code{od_plot}.
}

\seealso{
\code{\link{od_pool}}, \code{\link{od_print}}
}

\examples{

# Compute a D-optimal approximate design
# for the 2nd degree Fourier regression on a partial circle
# Use several types of graphs to visualize the design

Fx <- Fx_cube(~I(cos(x1)) + I(sin(x1)) + I(cos(2*x1)) + I(sin(2*x1)),
              lower = -2*pi/3, upper = 2*pi/3, n.levels = 121)
w <- od_REX(Fx)$w.best
par(mfrow = c(2, 2))
od_plot(Fx, w, X = seq(-2*pi/3, 2*pi/3, length = 121), main = "Plot 1")
od_plot(Fx, w, X = Fx[, 2:3], asp = 1, main = "Plot 2")
od_plot(Fx, w, X = Fx[, c(2,5)], asp = 1, main = "Plot 3")
od_plot(Fx, w, X = Fx[, c(3,4)], asp = 1, main = "Plot 4")
par(mfrow = c(1, 1))

\dontrun{
# Compute an I-efficient exact design of size 20 without replications
# for the Scheffe mixture model
# Use several types of graphs to visualize the design

Fx <- Fx_simplex(~x1 + x2 + x3 + I(x1*x2) + I(x1*x3) + I(x2*x3) - 1, 21)
w <- od_AQUA(Fx, b3=20, bin=TRUE, crit="I")$w.best
X <- Fx[, 1:2] %*% matrix(c(1, 1/2, 0, sqrt(3)/2), ncol = 2)
colnames(X) <- c("", "")
od_plot(Fx, w, X, asp = 1, main = "Plot 1")
od_plot(Fx, w, Fx[, 1:3], main = "Plot 2")

# Compute a symmetrized D-optimal approximate design
# for the full quadratic model with 4 factors
# Use several types of graphs to visualize the design

form.q <- ~x1 + x2 + x3 + x4 + I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) +
           I(x1*x2) + I(x1*x3) + I(x1*x4) + I(x2*x3) + I(x2*x4) + I(x3*x4)
Fx <- Fx_cube(form.q, n.levels = rep(11, 4))
w <- od_REX(Fx)$w.best
od_plot(Fx, w, Fx[, 2:3], dd.size=3)
od_plot(Fx, w, Fx[, 2:4], w.lim=Inf)


# A more complex example:

# Compute the D-optimal 17 point exact design
# for the spring-balance weighing model with 4 items

Fx <- Fx_cube(~x1 + x2 + x3 + x4 - 1, lower = rep(0, 4))
w <- od_KL(Fx, 17, t.max = 5)$w.best
od_print(Fx, w)$design
U <- eigen(diag(4) - 0.25 * rep(1, 4) %*% t(rep(1, 4)))$vectors

# A 2D visualization

X <- Fx[, 1:4] %*% U[, 3:4]
X[, 2] <- -2*X[, 2]
colnames(X) <- c("V", "Number of items on the pan")
od_plot(Fx, w+0.001, X)
for(i in 1:16) for(j in 1:16)
  if(sum(abs(Fx[i,1:4]-Fx[j,1:4]))==1)
    lines(X[c(i,j),1], X[c(i,j),2])

# A 3D visualization

X <- Fx[, 1:4] %*% U[, 1:3]
colnames(X) <- c("V1", "V2", "V3")
od_plot(Fx, w+0.001, X)
for(i in 1:16) for(j in 1:16)
  if(sum(abs(Fx[i, 1:4] - Fx[j, 1:4])) == 1)
    rgl::lines3d(X[c(i, j), 1], X[c(i, j), 2], X[c(i, j), 3])
}
}
