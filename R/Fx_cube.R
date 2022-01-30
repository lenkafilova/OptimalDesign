Fx_cube <- function(formula, lower=NULL, upper=NULL, n.levels=NULL, echo=TRUE) {
  # Generate Fx for a linear factor model (on a cuboid grid)
 
  cl <- match.call()
  verify(cl, formula = formula, lower = lower, upper = upper,
         n.levels = n.levels, echo = echo)
 
  av <- all.vars(formula); d <- length(av)
  if (is.null(lower)) lower <- rep(-1, d)
  if (is.null(upper)) upper <- rep(1, d)
  if (is.null(n.levels)) n.levels <- rep(2, d)

  lst <- labs <- c()
  for (i in 1:d) {
    lst <- c(lst, list(seq(lower[i], upper[i], length = n.levels[i])))
    labs <- c(labs,paste("x", as.character(i), sep = ""))
  }
  data <- expand.grid(lst)
  names(data) <- labs
  mf <- model.frame(formula, data)
  Fx <- as.data.frame(model.matrix(attr(mf, "terms"), mf))
  
  return(as.matrix(Fx))
}
