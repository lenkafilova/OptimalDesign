od_PUK <- function(Fx, w, echo=TRUE) {
    # A variant of the efficient rounding of Pukelsheim and Reider

    cl <- match.call()
    verify(cl, Fx = Fx, w = w, echo = echo)
    n <- nrow(Fx); m <- ncol(Fx)
    N <- floor(sum(w) + sqrt(.Machine$double.eps))
    supp <- w > sqrt(.Machine$double.eps); n.supp <- sum(supp)

    if (n.supp > N)
       stop("od_PUK assumes that N is not smaller than the size of the support of w.") 
    if (n.supp < m)
       warning("The resulting design may be deficient if the support of w is smaller than the number of parameters ncol(Fx).") 

    w.supp <- w[supp] + sqrt(.Machine$double.eps)*runif(n.supp)
    w.temp <- ceiling(w.supp*(1 - n.supp/2/N))

    res <- N - sum(w.temp)
    if (min(res) < 0) { # Greedy removal
         for (i in 1:(-res)) {
             i.max <- which.max((w.temp - 1)/w.supp)
             w.temp[i.max] <- w.temp[i.max] - 1
         }
    }

    if (min(res) > 0) { # Greedy augmentation
        for (i in 1:res) {
             i.min <- which.min(w.temp/w.supp)
             w.temp[i.min] <- w.temp[i.min] + 1
        }
    }

    w.res <- rep(0, n); w.res[(1:n)[supp]] <- w.temp
    
    return(list(call = cl, w.round = w.res))
}
