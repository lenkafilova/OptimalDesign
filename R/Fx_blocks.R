Fx_blocks <- function(n.treats, blocks=NULL, echo=TRUE) {
    # Generate Fx for a block model with size-two blocks
   
    cl <- match.call()
    verify(cl, n.treats = n.treats, blocks = blocks, echo = echo)

    if (is.null(blocks)) blocks <- utils::combn(n.treats, 2)
    n <- ncol(blocks)
    Gx <- matrix(0, nrow = n, ncol = n.treats)
    for (i in 1:n) {   
       Gx[i, blocks[1, i]] <- 1
       Gx[i, blocks[2, i]] <- -1
    }
    
    one.n.treats <- rep(1, n.treats) 
    V <- eigen(diag(n.treats) - one.n.treats %*% t(one.n.treats)/n.treats,
               symmetric = TRUE)$vectors[, 1:(n.treats - 1)]
    Fx <- as.matrix(Gx %*% V)
    cnms <- rep("", n.treats - 1)
    for (j in 1:(n.treats - 1)) cnms[j] <- paste("B", j, sep = "")
    colnames(Fx) <- cnms
    
    return(Fx)
}
