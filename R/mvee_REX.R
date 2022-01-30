mvee_REX <- function(Data, alg.AA="REX", eff=0.999999, it.max=Inf, t.max=60,
                     picture=FALSE, echo=TRUE, track=TRUE) {
    # Computes the minumum-volume enclosing ellipsoid containing the rows of Data
    # For two- and three-dimensional rows, produces a picture 
    # Important: It is expected that the procedure will be run with eff very close to 1
    # For eff significantly smaller than 1, the result can only be a rough approximation
    
    # TODO: Define the eff at the input/output via the ratio of volumes of ellipsoids
    
    cl <- match.call()
    verify(cl, Data = Data, alg.AA = alg.AA, eff = eff, t.max = t.max,
           picture = picture, echo = echo, track = track)
    n <- nrow(Data); d <- ncol(Data)
    
    Fx <- cbind(rep(1, n), Data)
    if (rcond(infmat(Fx, od_PIN(Fx, echo = FALSE)$w.pin, echo = FALSE)) < sqrt(.Machine$double.eps)) {
        print("Data is:"); print.default(Data, max = 100)
        stop("The points probably lie on a common hyperplane; I cannot compute the MVEE.")
    }
    if (alg.AA == "REX") {
        res <- od_REX(Fx, eff = eff, it.max = it.max, t.max = t.max, track = track)
    } else if (alg.AA == "VDM") {
        res <- od_REX(Fx, eff = eff, it.max = it.max, alg.AA = "VDM", t.max = t.max, track = track)
    } else {
        res <- od_REX(Fx, eff = eff, it.max = it.max, alg.AA = "MUL", t.max = t.max, track = track)
    }
    
    M <- infmat(Fx, res$w.best, echo = FALSE)
    z <- M[2:(d + 1), 1]
    H <- solve(d * (M[2:(d + 1), 2:(d + 1)] - z %*% t(z)))
    vol <- 1/(prod(sqrt(eigen(H, symmetric = TRUE)$values))) * pi^(d/2) * 1/gamma(d/2 + 1)
    bpts <- od_DEL(Fx, res$w.best)$keep
    
    if (picture) {
        if (d > 3)
            warning("The picture for data in dimension more than 3 cannot be drawn.")

        eig <- eigen(solve(H), symmetric = TRUE)
        Lamb <- diag(eig$values); U <- eig$vectors
        if (d == 2) {
            t <- seq(0, 2*pi, length = 1000)
            x <- t(U %*% sqrt(Lamb) %*% rbind(cos(t), sin(t))) +
                matrix(rep(z, 1000), byrow = TRUE, ncol = 2)
            cnm <- colnames(Data); if (is.null(cnm)) cnm <- c("X1", "X2")
            plot(x[, 1], x[, 2], xlab = cnm[1], ylab = cnm[2], asp = 1, type = "n")
            polygon(x[, 1], x[, 2], col = "mediumpurple1", border = "mediumpurple1")
            points(Data, pch = 16, cex = 1)
            points(Data[bpts, ], pch = 19, cex = 1, col = "red")
        }
        if (d == 3) {
            Y <- sqrt(Lamb) %*% t(U) 
            Y <- rbind(Y, -Y) + matrix(rep(z, 6), byrow = TRUE, ncol = 3) 
            cnm <- colnames(Data); if (is.null(cnm)) cnm <- c("X1", "X2", "X3")
            rgl::plot3d(Y[, 1], Y[, 2], Y[, 3], xlab = cnm[1], ylab = cnm[2], zlab = cnm[3], 
                   type = "n", box = TRUE, aspect = FALSE)
            bex <- 0.03*sqrt(max(eig$values))
            if (length(bpts) < n)
                rgl::plot3d(Data[-bpts, 1], Data[-bpts, 2], Data[-bpts, 3], type = "s",
                       radius = bex, add = TRUE)
            rgl::plot3d(Data[bpts, 1], Data[bpts, 2], Data[bpts, 3], type = "s",
                   col = "red", radius = bex, add = TRUE)
            ell <- rgl::ellipse3d(solve(H), centre = z, level = 0.205)
            rgl::plot3d(ell, col = "purple", alpha = 0.2, add = TRUE, box = FALSE)
        }
    }
    
    return(list(call = cl, H = H, z = z, bpts = bpts, vol = vol, 
                eff.best = res$eff.best, n.iter = res$n.iter, t.act = res$t.act))
}

