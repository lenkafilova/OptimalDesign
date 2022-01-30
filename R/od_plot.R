od_plot <- function(Fx, w, X=NULL, w.pool=c("sum", "0"), w.color="darkblue",
                    w.size=1, w.pch=16, w.cex=0.8, w.lim=0.01, crit="D", h=NULL,
                    dd.pool=c("max", "mean"), dd.color="orange", dd.size=1.5, dd.pch=15,
                    asp = NA, main.lab="", y.lab="", return.pools=FALSE, echo=TRUE) {
  # Draws a plot of (pools of) w
  # Also draws a plot of (pools of) the vector of directional derivatives
  # Pools can help see various aspects of the multidimensional design
  #
  # If d (=ncol(X)) is 1, it is possible to add more pools to w and dd
  # If d is 2, only the first pools of w and dd are relevant
  # If d is 3, only the first pool of w is relevant
  #
  # Notes:
  # The weight numbers are rounded to 2 decimal places
  # The weight numbers of size less than w.lim*sum(w) do no appear in the plot
  # The output of the function is a matrix of the Pools
  # Use this matrix to see the exact pool values and
  # to create a plot according to your needs

  cl <- match.call()

  verify(cl, Fx = Fx, w = w, X = X, w.pool = w.pool, w.color = w.color,
         w.size = w.size, w.pch, w.cex = w.cex, w.lim = w.lim, crit = crit, h = h,
         dd.pool = dd.pool, dd.color = dd.color, dd.size = dd.size, dd.pch,
         asp = asp, main.lab = main.lab, y.lab = y.lab, return.pools = return.pools, echo = echo)
  n <- nrow(Fx); m <- ncol(Fx)
  nwp <- length(w.pool); nddp <- length(dd.pool)

  if (crit %in% c("C", "c") && is.null(h)) h <- c(rep(0, m - 1), 1)
  if (is.null(w) || sum(w) < sqrt(.Machine$double.eps)) {
    print("w is:"); print.default(w, max = 100)
    stop("Cannot plot if w contains only zeros.")
  }
  if (is.null(X)) {
    message("X is NULL; setting X to 1:nrow(Fx)."); X <- 1:n
  }
  X <- as.matrix(X); d <- ncol(X)
  if (d < 1 || d > 3) stop("ncol(X) can only be 1, 2, or 3.")

  res <- od_pool(X, w, w.pool[1], echo = FALSE)
  Pool.X <- res$X.unique; nx <- nrow(Pool.X)
  if (nx == 1) stop("X must contain at least two distinct rows.")
  Pool.w <- matrix(0, nrow = nx, ncol = nwp)
  Pool.w[, 1] <- res$val.pooled
  if (nwp > 1) {
    for (k in 2:nwp) {
      res <- od_pool(X, w, w.pool[k], echo = FALSE)
      Pool.w[, k] <- res$val.pooled
    }
  }

  dd.vals <- dirder(Fx, w, crit = crit, h = h, echo = FALSE)

  if (d == 1) {
    wsc <- as.vector(0.5*col2rgb(w.color)/256 + 0.5*c(1, 1, 1))
    w.stick.color <- rgb(wsc[1], wsc[2], wsc[3])

    Pool.dd <- matrix(0, nrow = nx, ncol = nddp)
    for (k in 1:nddp) {
      res <- od_pool(X, dd.vals, dd.pool[k], echo = FALSE)
      Pool.dd[, k] <- res$val.pooled
    }

    x <- as.vector(Pool.X)
    x.min <- min(Pool.X) - 0.03*(max(Pool.X) - min(Pool.X))
    x.max <- max(Pool.X) + 0.03*(max(Pool.X) - min(Pool.X))
    y.min <- min(0, min(Pool.dd) - 0.001); y.max <- max(0, max(Pool.dd) + 0.001)
    y.max <- max(y.max, (y.max - y.min)/2)

    plot(x, rep(0, nx), xlim = c(x.min, x.max), ylim = c(y.min, y.max),
         type = "n", xlab = colnames(X)[1], ylab = y.lab, main = main.lab)
    grid(col = "darkgray"); lines(c(x.min, x.max), c(0, 0))

    for (k in 1:nddp) {
      points(x, Pool.dd[, k], pch = 15 + k, col = dd.color, cex = dd.size)
    }

    max.w <- max(Pool.w)
    r <- function(w) {w*(y.max - y.min)/4/max.w}

    rMinw <- matrixStats::rowMins(Pool.w)
    rMaxw <- matrixStats::rowMaxs(Pool.w)
    stick.indices <- (1:nx)[rMaxw > 1e-5*sum(w)]
    text.indices <- (1:nx)[rMaxw > w.lim*sum(w)]

    for (j in stick.indices)
      lines(c(x[j], x[j]), c(r(rMinw[j]), r(rMaxw[j])), col = w.stick.color, lwd = 2)
    for (k in 1:nwp) {
      points(x, r(Pool.w[, k]), pch = 15 + k, col = w.color, cex = w.size)
    }
    shift.text <- 0.3*r(max.w)
    for (j in text.indices)
      text(x[j], r(rMaxw[j]) + shift.text, cex = w.cex, font = 2, labels = round(rMaxw[j], 2))
  }

  if (d == 2) {

    res <- od_pool(X, dd.vals, dd.pool[1], echo = FALSE)
    Pool.dd <- matrix(res$val.pooled, nrow = nx, ncol = 1)
    y.min <- min(0, min(Pool.dd) - 0.001); y.max <- max(0, max(Pool.dd) + 0.001)
    cl.alpha <- as.vector(((Pool.dd[, 1]) - y.min)/(y.max - y.min))
    clr <- cl.alpha %*% t(as.vector(col2rgb(dd.color)/256)) + (1 - cl.alpha) %*% t(c(1, 1, 1))
    span1 <- max(Pool.X[, 1]) - min(Pool.X[, 1])
    span2 <- max(Pool.X[, 2]) - min(Pool.X[, 2])
    x1.max <- 0.08*span1 + max(Pool.X[, 1])
    x2.max <- 0.08*span2 + max(Pool.X[, 2])
    plot(Pool.X[, 1], Pool.X[, 2], xlim = c(min(Pool.X[, 1]), x1.max), type = "n",
         ylim = c(min(Pool.X[, 2]), x2.max), xlab = colnames(X)[1], ylab = colnames(X)[2],
         asp = asp, main = main.lab)
    grid(col = "darkgray")
    points(Pool.X[, 1], Pool.X[, 2], pch = dd.pch, cex = 1.5*dd.size,
         col = rgb(clr[, 1], clr[, 2], clr[, 3]))
    points(Pool.X[, 1], Pool.X[, 2], pch = w.pch, col = w.color,
           cex = 5*w.size*sqrt(as.vector(Pool.w[, 1]/sum(Pool.w[, 1]))))

    # TODO: This cycle is extremely inefficient for even not so large nx
    # Similarly in other cycles below.
    rMaxw <- matrixStats::rowMaxs(Pool.w)
    text.indices <- (1:nx)[rMaxw > w.lim*sum(w)]
    for (j in text.indices)
        text(Pool.X[j, 1] + 0.05*span1, Pool.X[j, 2]  + 0.05*span2,
             labels = round(max(Pool.w[j, ]), 2), cex = w.cex, font = 2)
  }

  if (d == 3) {
    span1 <- max(Pool.X[, 1]) - min(Pool.X[, 1])
    span2 <- max(Pool.X[, 2]) - min(Pool.X[, 2])
    span3 <- max(Pool.X[, 3]) - min(Pool.X[, 3])
    enlarge <- 0.1*max(c(span1, span2, span3))*w.size
    rgl::plot3d(Pool.X[, 1], Pool.X[, 2], Pool.X[, 3], xlab = colnames(X)[1],
                ylab = colnames(X)[2], zlab = colnames(X)[3], type = "s", col = w.color,
                radius = enlarge*(pmax(as.vector(Pool.w[, 1]/sum(Pool.w[, 1])), 0.0001))^(1/3),
                box = TRUE, aspect = TRUE)

    rMaxw <- matrixStats::rowMaxs(Pool.w)
    text.indices <- (1:nx)[rMaxw > w.lim*sum(w)]
    for (j in text.indices)
        rgl::text3d(Pool.X[j, 1] + 0.1*span1, Pool.X[j, 2]  + 0.1*span2, Pool.X[j, 3]  + 0.1*span3,
             texts = round(max(Pool.w[j, 1]), 2), cex = w.cex, font = 2)
    rgl::decorate3d(main = main.lab, font = 2, xlab = "", ylab = "", zlab = "")
  }

  if (d < 3) {
    res <- cbind(Pool.X = Pool.X, Pool.w = Pool.w, Pool.dd = Pool.dd)
  } else {
    res <- cbind(Pool.X = Pool.X, Pool.w = Pool.w)
  }
  if (!return.pools) res <- "Pool return suppressed."

  return(list(call = cl, Pools = res))
}

