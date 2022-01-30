Fx_simplex <- function(formula, n.levels.mix=NULL, echo=TRUE) {
  # Generate Fx for a linear mixture model (on a discrete simlex)
  
  cl <- match.call()
  verify(cl, formula = formula, n.levels.mix = n.levels.mix, echo = echo)
  
  av <- all.vars(formula); d <- length(av)
  if (is.null(n.levels.mix)) {
    message("n.levels.mix is NULL; setting n.levels.mix to 2*d+1.")
    n.levels.mix <- 2*d + 1
  }
 
  cmb <- t(combn(n.levels.mix + d - 2, d - 1))
  data <- as.data.frame((cbind(cmb, n.levels.mix + d - 1) -
                           cbind(0, cmb) - 1)/(n.levels.mix - 1))
  
  labs <- c()
  for (i in 1:d)
    labs <- c(labs, paste("x", as.character(i), sep = ""))
  
  names(data) <- labs[1:d]
  mf <- model.frame(formula, data)
  Fx <- as.data.frame(model.matrix(attr(mf, "terms"), mf))
  
  return(as.matrix(Fx))
}
