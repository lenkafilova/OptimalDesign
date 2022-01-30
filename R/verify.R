verify <- function(call.string, ...) {
   # Note: Testing (in)equality should be more liberal
   # Note: It is possible to print the actual value of the failed variable

   l <- list(...); args <- names(l)

   # Boolean echo must always be one of the arguments to be verified.
   echo <- l$echo
   if (!is.vector(echo) || !is.logical(echo) || length(echo) > 1) {
      message(paste("echo is", echo))
      stop("echo must be a single logical value (TRUE or FALSE).", call. = FALSE)
   }
   if (echo) {
      print("Call of the function:", quote = FALSE)
      print(call.string, quote = FALSE)
   }
   
   # nxm, m>=2, n>=m, real matrix Fx, cannot be NULL
   if ("Fx" %in% args) {
      Fx <- l$Fx
      if (!is.matrix(Fx) || !is.numeric(Fx) || !all(is.finite(Fx))) { 
         print.default("Fx is:"); print.default(Fx, max = 100)
         stop("Fx must be a matrix of real numbers.", call. = FALSE)
      }
      n <- nrow(Fx); m <- ncol(Fx)
      if (n < m) {
         print.default("Fx is:"); print.default(Fx, max = 100)
         stop("The number n=nrow(Fx) of design points must be greater than or equal to the number m=ncol(Fx) of model parameters.", call. = FALSE)
      }
      if (m < 2) {
         print.default("Fx is:"); print.default(Fx, max = 100)
         stop("The number m=ncol(Fx) of model parameters must be at least 2. The case of m=1 is an (often trivial) problem of LP or ILP.", call. = FALSE)
      }
   }
   
   # lenght-n real non-negative vector w cannot be NULL
   if ("w" %in% args) {
      w <- l$w
      if (!is.vector(w) || !is.numeric(w) || !all(is.finite(w)) ||
          min(w) < 0 || max(w) <= 0 || length(w) != n) {
         print.default("w is:"); print.default(w, max = 100)
         # For some procedures, the zero vector w is theoretically permissible,
         # but it is probably never a practically interesting as an input argument.
         stop("w must be non-negative, nonzero, finite numeric vector of length nrow(Fx).", call. = FALSE)
      }
   }

   # Character crit cannot be NULL
   if ("crit" %in% args) {
      crit <- l$crit
      if (!is.vector(crit) || !is.character(crit) || (length(crit) != 1) ||
          !(crit %in% c("D", "A", "I", "C", "c"))) {
         print.default("crit is:"); print.default(crit, max = 100)
         stop("crit must be 'D', 'A', 'I', 'C', or 'c'.", call. = FALSE)
      }
   }

   # Real non-zero vector h can be NULL; then it is set to c(0,...,0,1) in the wrapper.
   if ("h" %in% args) {
      h <- l$h
      if (!is.null(h)) { 
         if (!is.vector(h) || !is.numeric(h) || !all(is.finite(h)) ||
             length(h) != m || isTRUE(all.equal(h, rep(0, m)))) {
            print.default("h is:"); print.default(h, max = 100)
            stop("h must be NULL or a nonzero real vector of length ncol(Fx).", call. = FALSE)
            if (crit != "c" && crit != "C")
               warning("h is set to a non-NULL value but it will be ignored because crit is not 'c' nor 'C'.", call. = FALSE)
         }
      }
   }

   # n.treats cannot be NULL
   if ("n.treats" %in% args) {
      n.treats <- l$n.treats
      if (!is.numeric(n.treats) || !is.finite(n.treats) ||
          length(n.treats) != 1 || n.treats < 3 || n.treats != round(n.treats)) {
         print.default("n.treats is:"); print.default(n.treats, max = 100)
         stop("n.treats must be a natural number greater than 2.", call. = FALSE)
      }
   }
            
   # blocks can be NULL; then it is set to the matrix of all pairs in the wrapper
   if ("blocks" %in% args) {
      blocks <- l$blocks
      if (!is.null(blocks)) {
         if (!is.matrix(blocks) || !is.numeric(blocks) || !all(is.finite(blocks)) ||
             nrow(blocks) != 2 || !all(is.element(blocks, 1:n.treats))) {
            print.default("blocks is:"); print.default(blocks, max = 100)
            stop("blocks must be a matrix with two rows and elements in 1:n.treats.", call. = FALSE)
         }
      }
   }

   # formula cannot be NULL
   if ("formula" %in% args) {
      formula <- l$formula
      if (!plyr::is.formula(formula)) {
         print.default("formula is:"); print.default(formula, max = 100)
         stop("The argument 'formula' does not seem to be a valid formula.", call. = FALSE)
      }
      av <- all.vars(formula); d <- length(av)
      pot.vars <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9")
      if (length(intersect(av, pot.vars[1:d])) != d) {
         print.default("formula is:"); print.default(formula, max = 100)
         stop("Variables in the formula must be x1,x2,...,xd, d<=9, without skipping any intermediate variable.", call. = FALSE) 
      }
   }

   # lower can be NULL; in that case it is set to c(-1,...,-1)
   if ("lower" %in% args) {
      lower <- l$lower
      if (!is.null(lower)) {
         if (!is.numeric(lower) || !all(is.finite(lower)) || !is.vector(lower) || 
             length(lower) != d) {
            print.default("lower is:"); print.default(lower, max = 100)
            print.default("formula is:"); print.default(formula, max = 100)
            stop("lower must be NULL, a real number, or a real vector of length d, where d is the number of variables x1,x2,...,xd in the formula.", call. = FALSE)
         }
      }
   }
   
   # upper can be NULL; in that case it is set to c(1,...,1)
   if ("upper" %in% args) {
      upper <- l$upper
      if (!is.null(upper)) {
         if (!is.numeric(upper) || !all(is.finite(upper)) || !is.vector(upper) || 
             length(upper) != d) {
            print.default("upper is:"); print.default(upper, max = 100)
            print.default("formula is:"); print.default(formula, max = 100)
            stop("upper must be NULL, a real number, or a real vector of length d, where d is the number of variables x1,x2,...,xd in the formula.", call. = FALSE)
         }
         if (any(lower >= upper)) {
            print.default("lower is:"); print.default(lower, max = 100)
            print.default("upper is:"); print.default(upper, max = 100)
            stop("lower[i] must be smaller than upper[i] for all i!", call. = FALSE)
         }
      }
   }
   
   # n.levels can be NULL;in that case it is set to c(2,...,2)
   if ("n.levels" %in% args) {
      n.levels <- l$n.levels
      if (!is.null(n.levels)) {
         if (!is.numeric(n.levels) || !all(is.finite(n.levels)) || !is.vector(n.levels) ||
             !all(n.levels == round(n.levels)) || length(n.levels) != d || min(n.levels) < 2) {
            print.default("n.levels is:"); print.default(n.levels, max = 100)
            print.default("formula is:"); print.default(formula, max = 100)
            stop("n.levels must be NULL, a natural number >=2, or a vector of length d consisting of natural numbers >=2, where d is the number of variables x1,x2,...,xd in the formula.", call. = FALSE)
         }
      }
   }
   
   # dose.levels cannot be NULL
   if ("dose.levels" %in% args) {
      dose.levels <- l$dose.levels
      if (!is.vector(dose.levels) || !is.numeric(dose.levels) ||
          !all(is.finite(dose.levels)) || min(dose.levels) < 0) {
         print.default("dose.levels is:"); print.default(dose.levels, max = 100)
         stop("dose.levels must be a vector of non-negative numbers.", call. = FALSE)
      }
   }

   # theta0 cannot be NULL
   if ("theta0" %in% args) {
      theta0 <- l$theta0
      if (!is.vector(theta0) || !is.numeric(theta0) || !all(is.finite(theta0))) {
         print.default("theta0 is:"); print.default(theta0, max = 100)
         stop("theta0 must be a real vector.", call. = FALSE)
      }
      if ("dose.model" %in% args) {
         dose.model <- l$dose.model
         if (length(theta0) != 3 || theta0[2] == 0 || theta0[3] == 0) {
            print.default("theta0 is:"); print.default(theta0, max = 100)
            stop("In the implemented dose models, theta0 must have length 3 and a both theta[2] and theta0[3] must be non-zero", call. = FALSE)
         }
         if (dose.model %in% c("loglin", "exp") && theta0[3] < 0) {
            print.default("theta0 is:"); print.default(theta0, max = 100)
            stop("In the loglin and exp dose models theta0[3] must be positive.", call. = FALSE)
         }
      }
   }

   # dose.model cannot be NULL
   if ("dose.model" %in% args) {
      dose.model <- l$dose.model
      if (!is.vector(dose.model) || !is.character(dose.model) ||
          (length(dose.model) != 1) ||
          !(dose.model %in% c("emax", "loglin", "exp"))) {
         print.default("dose.model is:"); print.default(dose.model, max = 100)
         stop("dose.model must be 'emax', 'loglin', or 'exp'.", call. = FALSE)
      }
   }

   # glm.model cannot be NULL
   if ("glm.model" %in% args) {
      glm.model <- l$glm.model
      if (!is.vector(glm.model) || !is.character(glm.model) || length(glm.model) != 1 ||
          !(glm.model %in% c("bin-logit", "bin-probit", "bin-cloglog", "Poisson-log",
                             "gamma-BoxCox"))) {
         print.default("glm.model is:"); print.default(glm.model, max = 100)
         stop("glm.model must be 'bin-logit', 'bin-probit', 'bin-cloglog', 'Poisson-log', or 'gamma-BoxCox'.", call. = FALSE)
      }
   }

   # natural number n.levels.mix (>=2) cannot be NULL
   if ("n.levels.mix" %in% args) {
      n.levels.mix <- l$n.levels.mix
      if (!is.null(n.levels.mix)) {
         if (!is.vector(n.levels.mix) || !is.numeric(n.levels.mix) ||
             !all(is.finite(n.levels.mix)) || length(n.levels.mix) != 1 ||
             n.levels.mix < 2 || n.levels.mix != round(n.levels.mix)) {
            print.default("n.levels.mix is:"); print.default(n.levels.mix, max = 100)
            stop("n.levels.mix must be a natural number >=2.", call. = FALSE)
         }
      }
   }
   
   # censor.time cannot be NULL
   if ("censor.time" %in% args) {
      censor.time <- l$censor.time
      if (!is.numeric(censor.time) || !is.finite(censor.time) || 
          length(censor.time) != 1 || censor.time <= 0) {
         print.default("censor.time is:"); print.default(censor.time, max = 100)
         stop("censor.time must be a real positive number.", call. = FALSE)
      }
   }
   
   # survival.model cannot be NULL
   if ("survival.model" %in% args) {
      survival.model <- l$survival.model
      if (!is.vector(survival.model) || !is.character(survival.model) || (length(survival.model) != 1) ||
          !(survival.model %in% c("phI", "phrand"))) {
         print.default("survival.model is:"); print.default(survival.model, max = 100)
         stop("survival.model must be 'phI' or 'phrand'.", call. = FALSE)
      }
   }
     
   # nxd real matrix Data cannot be NULL 
   if ("Data" %in% args) {
      # Pozor na to, aby procedury nespadli pre ziadne rozmery n,d matice X
      Data <- l$Data
      if ((!is.vector(Data) && !is.matrix(Data)) || !is.numeric(Data) ||
          !all(is.finite(Data))) {
         print.default("Data is:"); print.default(Data, max = 100)
         stop("Data must be a matrix of real numbers.", call. = FALSE)
      }
      n <- nrow(Data); d <- ncol(Data)
      if (d < 2) {
         print.default("Data is:"); print.default(Data, max = 100)
         stop("Data is 1 dimensional; the MVEE is just the line segment from the minimum to the maximum.", call. = FALSE)
      }
   }
    
   # Character alg.AA cannot be NULL
   if ("alg.AA" %in% args) {
      alg.AA <- l$alg.AA
      if (!is.vector(alg.AA) || !is.character(alg.AA) || length(alg.AA) != 1 ||
          !(alg.AA %in% c("REX", "VDM", "MUL"))) {
         print.default("alg.AA is:"); print.default(alg.AA, max = 100)
         stop("alg.AA must be 'REX', 'VDM', or 'MUL'", call. = FALSE)
      }
   }
           
   # the real number eff must be in [0, 1]
   if ("eff" %in% args) {
      eff <- l$eff
      if (!is.numeric(eff) || !is.finite(eff) || length(eff) != 1 ||
          eff < 0 || eff > 1) {
         print.default("eff is:"); print.default(eff, max = 100)
         stop("eff must be a real number in [0,1].", call. = FALSE)
      }
   }

   # Non-negative real number t.max cannot be NULL (but it can be Inf)
   if ("t.max" %in% args) {
      t.max <- l$t.max
      if (!is.vector(t.max) || !is.numeric(t.max) || length(t.max) != 1 ||
          t.max < 0) {
         print.default("t.max is:"); print.default(t.max, max = 100)
         stop("t.max must be a non-negative number or Inf.", call. = FALSE)
      }
   }
               
   # Boolean argument picture cannot be NULL
   if ("picture" %in% args) {
      picture <- l$picture
      if (!is.vector(picture) || !is.logical(picture) || length(picture) > 1) {
         print.default("picture is:"); print.default(picture, max = 100)
         stop("picture must be a single logical value (TRUE or FALSE).", call. = FALSE)
      }
   }
   
   # Boolean argument track cannot be NULL
   if ("track" %in% args) {
      track <- l$track
      if (!is.vector(track) || !is.logical(track) || length(track) > 1) {
         print.default("track is:"); print.default(track, max = 100)
         stop("track must be a single logical value (TRUE or FALSE).", call. = FALSE)
      }
   }

   ## Verify b1, b2, b3
   if (("b1" %in% args) && ("b2" %in% args) && ("b3" %in% args)) {
      b1 <- l$b1; b2 <- l$b2; b3 <- l$b3
      if (is.null(b1) && is.null(b2) && is.null(b3)) {
         print.default("b1, b2, b3 are:"); print.default(b1, max = 100)
         print.default(b2, max = 100); print.default(b3, max = 100)
         stop("At least one of b1, b2, b3 must be non-NULL.", call. = FALSE)
      }
   }
   
   # Real vector b1 can be NULL and then the real matrix A1 must also be NULL.
   if ("b1" %in% args) {
      b1 <- l$b1; A1 <- l$A1
      if (!is.null(b1)) {
         if (!is.vector(b1) || !is.numeric(b1) || !all(is.finite(b1))) {
            print.default("b1 is:"); print.default(b1, max = 100)
            stop("b1 must be either NULL, or a real number, or a vector of real numbers.", call. = FALSE)
         }
         if (!is.null(A1)) {
            if (!is.matrix(A1) || !is.numeric(A1) || !all(is.finite(A1))) {
               print.default("A1 is:"); print.default(A1, max = 100)
               stop("A1 must be NULL, or a matrix of real numbers.", call. = FALSE)
            }
            if (nrow(A1) != length(b1)) {
               print.default("b1 is:"); print.default(b1, max = 100)
               print.default("A1 is:"); print.default(A1, max = 100)
               stop("A1 must be NULL, or nrow(A1) must be equal to length(b1).", call. = FALSE)
            }
            if (ncol(A1) != n) {
               print.default("A1 is:"); print.default(A1, max = 100)
               print.default("Fx is:"); print.default(Fx, max = 100)
               stop("A1 must be NULL, or ncol(A1) must be equal to nrow(Fx).", call. = FALSE)
            }
         } else {
            if (!length(b1) == 1 || b1 <= 0) { 
               print.default("b1 is:"); print.default(b1, max = 100)
               print.default("A1 is:"); print.default(A1, max = 100)
               stop("If b1 is non-NULL and A1 is NULL then b1 must be a positive number.", call. = FALSE)
            }
         }
      } else {
         if (!is.null(A1)) {
            print.default("b1 is:"); print.default(b1, max = 100)
            print.default("A1 is:"); print.default(A1, max = 100)
            stop("If b1 is NULL then A1 must also be NULL.", call. = FALSE)
         }
      }
   }
   
   # Real vector b2 can be NULL and then the real matrix A2 must also be NULL.
   if ("b2" %in% args) {
      b2 <- l$b2; A2 <- l$A2
      if (!is.null(b2)) {
         if (!is.vector(b2) || !is.numeric(b2) || !all(is.finite(b2))) {
            print.default("b2 is:"); print.default(b2, max = 100)
            stop("b2 must be either NULL, or a real number, or a vector of real numbers.", call. = FALSE)
         }
         if (!is.null(A2)) {
            if (!is.matrix(A2) || !is.numeric(A2) || !all(is.finite(A2))) {
               print.default("A2 is:"); print.default(A2, max = 100)
               stop("A2 must be NULL, or a matrix of real numbers.", call. = FALSE)
            }
            if (nrow(A2) != length(b2)) {
               print.default("b2 is:"); print.default(b2, max = 100)
               print.default("A2 is:"); print.default(A2, max = 100)
               stop("A2 must be NULL, or nrow(A2) must be equal to length(b2).", call. = FALSE)
            }
            if (ncol(A2) != n) {
               print.default("A2 is:"); print.default(A2, max = 100)
               print.default("Fx is:"); print.default(Fx, max = 100)
               stop("A2 must be NULL, or ncol(A2) must be equal to nrow(Fx).", call. = FALSE)
            }
         } else {
            if (!length(b2) == 1 || b2 <= 0) { 
               print.default("b2 is:"); print.default(b2, max = 100)
               print.default("A2 is:"); print.default(A2, max = 100)
               stop("If b2 is non-NULL and A2 is NULL then b1 must be a positive number.", call. = FALSE)
            }
         }
      } else {
         if (!is.null(A2)) {
            print.default("b2 is:"); print.default(b2, max = 100)
            print.default("A2 is:"); print.default(A2, max = 100)
            stop("If b2 is NULL then A2 must also be NULL.", call. = FALSE)
         }
      }
   }
   
   # Real vector b3 can be NULL and then the real matrix A3 must also be NULL.
   if ("b3" %in% args) {
      b3 <- l$b3; A3 <- l$A3
      if (!is.null(b3)) {
         if (!is.vector(b3) || !is.numeric(b3) || !all(is.finite(b3))) {
            print.default("b3 is:"); print.default(b3, max = 100)
            stop("b3 must be either NULL, or a real number, or a vector of real numbers.", call. = FALSE)
         }
         if (!is.null(A3)) {
            if (!is.matrix(A3) || !is.numeric(A3) || !all(is.finite(A3))) {
               print.default("A3 is:"); print.default(A3, max = 100)
               stop("A3 must be NULL, or a matrix of real numbers.", call. = FALSE)
            }
            if (nrow(A3) != length(b3)) {
               print.default("b3 is:"); print.default(b3, max = 100)
               print.default("A3 is:"); print.default(A3, max = 100)
               stop("A3 must be NULL, or nrow(A3) must be equal to length(b3).", call. = FALSE)
            }
            if (ncol(A3) != n) {
               print.default("A3 is:"); print.default(A3, max = 100)
               print.default("Fx is:"); print.default(Fx, max = 100)
               stop("A3 must be NULL, or ncol(A3) must be equal to nrow(Fx).", call. = FALSE)
            }
         } else {
            if (!length(b3) == 1 || b3 <= 0) { 
               print.default("b3 is:"); print.default(b3, max = 100)
               print.default("A3 is:"); print.default(A3, max = 100)
               stop("If b3 is non-NULL and A3 is NULL then b1 must be a positive number.", call. = FALSE)
            }
         }
      } else {
         if (!is.null(A3)) {
            print.default("b3 is:"); print.default(b3, max = 100)
            print.default("A3 is:"); print.default(A3, max = 100)
            stop("If b3 is NULL then A3 must also be NULL.", call. = FALSE)
         }
      }
   }
   
   # Non-negative vector w0 can be NULL; then it is set to c(0,...,0) by the wrapper.
   if ("w0" %in% args) {
      w0 <- l$w0
      if (!is.null(w0)) {
         if (!is.vector(w0) || !is.numeric(w0) || !all(is.finite(w0)) ||
             min(w0) < 0) {
            print.default("w0 is:"); print.default(w0, max = 100)
            stop("w0 must be NULL, or a vector of non-negative real numbers.", call. = FALSE)
         }
         if (length(w0) != n) {
            print.default("w0 is:"); print.default(w0, max = 100)
            print.default("Fx is:"); print.default(Fx, max = 100)
            stop("The length of w0 must be equal to nrow(Fx).", call. = FALSE)
         }
      }
   }

   # Boolean bin cannot be NULL
   if ("bin" %in% args) {
      bin <- l$bin
      if (!is.vector(bin) || !is.logical(bin) || length(bin) > 1) {
         print.default("bin is:"); print.default(bin, max = 100)
         stop("bin must be a single logical value (TRUE or FALSE).", call. = FALSE)
      }
   }

   # The real non-singular mxm matrix M.anchor can be NULL
   if ("M.anchor" %in% args) {
      M.anchor <- l$M.anchor
      if (!is.null(M.anchor)) {
         if (!is.matrix(M.anchor) || !is.numeric(M.anchor) || !all(is.finite(M.anchor))) {
            print.default("M.anchor is:"); print.default(M.anchor, max = 100)
            stop("M.anchor must be a matrix of real numbers.", call. = FALSE)
         }
         m1 <- nrow(M.anchor); m2 <- ncol(M.anchor)
         if (m1 != m || m2 != m) {
            print.default("M.anchor is:"); print.default(M.anchor, max = 100)
            print.default("Fx is:"); print.default(Fx, max = 100)
            stop("M.anchor must be an ncol(Fx) times ncol(Fx) matrix.", call. = FALSE)
         }
         if (rcond(M.anchor) < sqrt(.Machine$double.eps)) 
            warning("M.anchor is badly conditioned and the results of AQUA may be unreliable.", call. = FALSE)
      }
   }
    
   # Character ver.qa cannot be NULL
   if ("ver.qa" %in% args) {
      ver.qa <- l$ver.qa
      if (!is.vector(ver.qa) || !is.character(ver.qa) || (length(ver.qa) != 1) ||
          !(ver.qa %in% c("+", "-"))) {
         print.default("ver.qa is:"); print.default(ver.qa, max = 100)
         stop("ver.qa must be either '+' or '-'.", call. = FALSE)
      }
   }

   # Boolean conic cannot be NULL
   if ("conic" %in% args) {
      conic <- l$conic
      if (!is.vector(conic) || !is.logical(conic) || length(conic) > 1) {
         print.default("conic is:"); print.default(conic, max = 100)
         stop("conic must be a single logical value (TRUE or FALSE).", call. = FALSE)
      }
   }

   # natural number N cannot be NULL
   if ("N" %in% args) {
      N <- l$N
      if (!is.vector(N) || !is.numeric(N) || !all(is.finite(N)) || 
          length(N) != 1 || N <= 0 || N != round(N)) {
         print.default("N is:"); print.default(N, max = 100)
         stop("The required size N of the design must be a natural number.", call. = FALSE)
      }
   }

   # Positive finite real number Phi.app can be NULL 
   if ("Phi.app" %in% args) {
      Phi.app <- l$Phi.app
      if (!is.null(Phi.app)) {   
         if (!is.vector(Phi.app) || !is.numeric(Phi.app) ||
             !is.finite(Phi.app) || length(Phi.app) != 1 || Phi.app <= 0) {
            print.default("Phi.app is:"); print.default(Phi.app, max = 100)
            stop("Phi.app must be NULL or a positive real number.", call. = FALSE)
         }
      }
   }   

   # Real non-negative length-n vector w1 must satisfy constraints or is NULL
   if ("w1" %in% args) {
      w1 <- l$w1
      if (!is.null(w1)) {
         if (!is.vector(w1) || !is.numeric(w1) || !all(is.finite(w1)) || min(w1) < 0) {
            print.default("w1 is:"); print.default(w1, max = 100)
            stop("w1 must be NULL, or a vector of non-negative real numbers.", call. = FALSE)
         }
         if (length(w1) != n) {
            print.default("w1 is:"); print.default(w1, max = 100)
            stop("The length of w1 must be equal to nrow(Fx).", call. = FALSE)
         }
         if ("K" %in% args) {
            # We get here only for od_KL
            if (!all(w1 == round(w1))) {
               print.default("w1 is:"); print.default(w1, max = 100)
               stop("w1 must be an exact (integer) design.", call. = FALSE)
            }
            if (sum(w1) != l$N) {
               print.default("w1 is:"); print.default(w1, max = 100)
               stop("w1 must satisfy the size constraint sum(w1)=N.", call. = FALSE)
            }
         }
         if ("b" %in% args) {
            # We get here only for od_RC
            if (!all(w1 == round(w1))) {
               print.default("w1 is:"); print.default(w1, max = 100)
               stop("w1 must be an exact (integer) design.", call. = FALSE)
            }
         }
      }
   }

   # The natural number (or Inf) K can be NULL
   if ("K" %in% args) {
      K <- l$K
      if (!is.null(K)) {
         if (!is.vector(K) || !is.numeric(K) || length(K) != 1 ||
             K != round(K) || K < 2) {
            print.default("K is:"); print.default(K, max = 100)
            stop("K must be NULL, or a natural number or Inf.", call. = FALSE)
         }
      }
   }
   
   # The natural number (or Inf) L can be NULL
   if ("L" %in% args) {
      L <- l$L
      if (!is.null(L)) {
         if (!is.vector(L) || !is.numeric(L) || length(L) != 1 ||
             L != round(L) || L < 2) {
            print.default("L is:"); print.default(L, max = 100)
            stop("L must be NULL, or a natural number or Inf.", call. = FALSE)
         }
      }
   }

   # Natural number rest.max cannot be NULL (but it can be Inf)
   if ("rest.max" %in% args) {
      rest.max <- l$rest.max
      if (!is.vector(rest.max) || !is.numeric(rest.max) ||
          length(rest.max) != 1 || rest.max != round(rest.max) || rest.max < 2) {
         print.default("rest.max is:"); print.default(rest.max, max = 100)
         stop("rest.max must be a natural number or Inf.", call. = FALSE)
      }
   }
   
   # Character type cannot be NULL
   if ("type" %in% args) {
      type <- l$type
      if (!is.vector(type) || !is.character(type) || (length(type) != 1) ||
          !(type %in% c("exact", "approximate"))) {
         print.default("type is:"); print.default(type, max = 100)
         stop("type must be 'exact', or 'approximate'.", call. = FALSE)
      }
   }
   
   # Non-negative real number gap cannot be NULL
   if ("gap" %in% args) {
      gap <- l$gap
      if (!is.null(gap)) {
         if (!is.vector(gap) || !is.numeric(gap) || length(gap) != 1 ||
             !is.finite(gap) || gap < 0) {
            print.default("gap is:"); print.default(gap, max = 100)
            stop("gap must be a non-negative number.", call. = FALSE)
         }
      }
   }
   
   # Character alg.PIN cannot be NULL
   if ("alg.PIN" %in% args) {
      alg.PIN <- l$alg.PIN
      if (!is.vector(alg.PIN) || !is.character(alg.PIN) || length(alg.PIN) != 1 ||
          !(alg.PIN %in% c("GKM", "KYM"))) {
         print.default("alg.PIN is:"); print.default(alg.PIN, max = 100)
         stop("alg.PIN must be 'GKM', or 'KYM'", call. = FALSE)
      }
   }

   # nxd real matrix X can be NULL (it can also be a vector of length n)
   if ("X" %in% args) {
      X <- l$X
      if (!is.null(X)) {
         if ((!is.vector(X) && !is.matrix(X)) || !is.numeric(X) || !all(is.finite(X))) {
            print.default("X is:"); print.default(X, max = 100)
            stop("X must be a vector or a matrix of real numbers.", call. = FALSE)
         }
         if ("Fx" %in% args) {
            if ((is.matrix(X) && nrow(Fx) != nrow(X)) || (is.vector(X) && nrow(Fx) != length(X))) {
               print.default("Fx is:"); print.default(Fx, max = 100)
               print.default("X is:"); print.default(X, max = 100)
               stop("X must be have the same number of rows as Fx (if Fx is among the arguments).", call. = FALSE)
            }
         }
         X <- as.matrix(X); n <- nrow(X); d <- ncol(X)
      }
   }
      
   # Real vector b can be NULL and then the real matrix A must also be NULL.
   if ("b" %in% args) {
      b <- l$b; A <- l$A
      if (!is.null(b)) {
         if (!is.vector(b) || !is.numeric(b) || !all(is.finite(b)) || min(b) <= 0) {
            print.default("b is:"); print.default(b, max = 100)
            stop("b must be a vector of positive real numbers.", call. = FALSE)
         }
         if (!is.null(A)) {
            if (!is.matrix(A) || !is.numeric(A) || !all(is.finite(A))) {
               print.default("A is:"); print.default(A, max = 100)
               stop("A must be NULL, or a matrix of real numbers.", call. = FALSE)
            }
            if (min(A) < 0 || any(apply(A, 2, max) <= 0)) {
               print.default("A is:"); print.default(A, max = 100)
               stop("The elements of A are not valid coefficients of resource consumption.", call. = FALSE)
            }
            if (nrow(A) != length(b)) {
               print.default("b is:"); print.default(b, max = 100)
               print.default("A is:"); print.default(A, max = 100)
               stop("A must be NULL or nrow(A) must be equal length(b).", call. = FALSE)
            }
            if (ncol(A) != n) {
               print.default("A is:"); print.default(A, max = 100)
               print.default("Fx is:"); print.default(Fx, max = 100)
               stop("A must be NULL or ncol(A) must be equal to nrow(Fx).", call. = FALSE)
            }
         } else {        
            if (length(b) != 1) {
               print.default("b is:"); print.default(b, max = 100)
               print.default("A is:"); print.default(A, max = 100)
               stop("If b is non-NULL and A is NULL then b must be a positive real number.", call. = FALSE)
            }
         }
      } else {
         stop("b cannot be NULL", call. = FALSE)
      }
   }

   # Natural number it.max cannot be NULL (but it can be Inf)
   if ("it.max" %in% args) {
      it.max <- l$it.max
      if (!is.vector(it.max) || !is.numeric(it.max) || length(it.max) != 1 ||
          it.max != round(it.max) || it.max < 2) {
         print.default("it.max is:"); print.default(it.max, max = 100)
         stop("it.max must be a natural number or Inf.", call. = FALSE)
      }
   }

   # w.pool is a vector of strings that contains "sum"
   if ("w.pool" %in% args) {
      w.pool <- l$w.pool
      if (!is.vector(w.pool) || !is.character(w.pool) || 
          !all(w.pool %in% c("sum", "max", "min", "mean", "median", "0")) ||
          !("sum" %in% w.pool)) {
         print.default("w.pool is:"); print.default(w.pool, max = 100)
         stop("w.pool must a vector of strings containing 'sum'. All components of w.pool must be selected from 'sum', 'min', 'mean', 'median', and '0'.", call. = FALSE)
      }
   }

   # Char w.color cannot be NULL. Does not test that w.color is a valid color string.
   if ("w.color" %in% args) {
      w.color <- l$w.color
      if (!is.vector(w.color) || !is.character(w.color) || length(w.color) != 1) {
         print.default("w.color is:"); print.default(w.color, max = 100)
         stop("w.color must be a single character string.", call. = FALSE)
      }
   }
   
   # Non-negative real number w.size cannot be NULL.
   if ("w.size" %in% args) {
      w.size <- l$w.size
      if (!is.vector(w.size) || !is.numeric(w.size) || length(w.size) != 1 ||
          w.size < 0 || !is.finite(w.size)) {
         print.default("w.size is:"); print.default(w.size, max = 100)
         stop("w.size must be a non-negative number.", call. = FALSE)
      }
   }
   
   if ("w.pch" %in% args) {
      w.pch <- l$w.pch
      if (!is.vector(w.pch) || !is.numeric(w.pch) || length(w.pch) != 1 ||
          w.pch < 0 || !is.finite(w.pch) || round(w.pch) != w.pch) {
         print.default("w.pch is:"); print.default(w.pch, max = 100)
         stop("w.pch must be a non-negative integer number.", call. = FALSE)
      }
   }
   
   # Non-negative real number w.cex cannot be NULL.
   if ("w.cex" %in% args) {
      w.cex <- l$w.cex
      if (!is.vector(w.cex) || !is.numeric(w.cex) || length(w.cex) != 1 ||
          w.cex < 0 || !is.finite(w.cex)) {
         print.default("w.cex is:"); print.default(w.cex, max = 100)
         stop("w.cex must be a non-negative number.", call. = FALSE)
      }
   }
   
   # Real number w.lim cannot be NULL (but it can be Inf)
   if ("w.lim" %in% args) {
      w.lim <- l$w.lim
      if (!is.vector(w.lim) || !is.numeric(w.lim) || length(w.lim) != 1) {
         print.default("w.lim is:"); print.default(w.lim, max = 100)
         stop("w.lim must be a real number or Inf.", call. = FALSE)
      }
   }
   
   # dd.pool is a vector of strings that cannot be empty
   if ("dd.pool" %in% args) {
      dd.pool <- l$dd.pool
      if (!is.vector(dd.pool) || !is.character(dd.pool) || length(dd.pool) < 1 ||
          !all(dd.pool %in% c("sum", "max", "min", "mean", "median", "0"))) {
         print.default("dd.pool is:"); print.default(dd.pool, max = 100)
         stop("dd.pool must a vector of strings containing 'sum'. All components of dd.pool must be selected from 'sum', 'min', 'mean', 'median', and '0'.", call. = FALSE)
      }
   }
   
   # Char dd.color cannot be NULL. Does not test that dd.color is a valid color string.
   if ("dd.color" %in% args) {
      dd.color <- l$dd.color
      if (!is.vector(dd.color) || !is.character(dd.color) || length(dd.color) != 1) {
         print.default("dd.color is:"); print.default(dd.color, max = 100)
         stop("dd.color must be a single character string.", call. = FALSE)
      }
   }
   
   # Non-negative real number dd.size cannot be NULL.
   if ("dd.size" %in% args) {
      dd.size <- l$dd.size
      if (!is.vector(dd.size) || !is.numeric(dd.size) || length(dd.size) != 1 ||
          dd.size < 0 || !is.finite(dd.size)) {
         print.default("dd.size is:"); print.default(dd.size, max = 100)
         stop("dd.size must be a non-negative number.", call. = FALSE)
      }
   }
   
   if ("dd.pch" %in% args) {
      dd.pch <- l$dd.pch
      if (!is.vector(dd.pch) || !is.numeric(dd.pch) || length(dd.pch) != 1 ||
          dd.pch < 0 || !is.finite(dd.pch) || round(dd.pch) != dd.pch) {
         print.default("dd.pch is:"); print.default(dd.pch, max = 100)
         stop("dd.pch must be a non-negative integer number.", call. = FALSE)
      }
   }
   
   # Real number asp cannot be NULL (but it can be NA)
   if ("asp" %in% args) {
      asp <- l$asp
      if (!is.na(asp)) {
         if (!is.vector(asp) || !is.numeric(asp) || length(asp) != 1) {
            print.default("asp is:"); print.default(asp, max = 100)
            stop("asp must be NA, or a real number.", call. = FALSE)
         }
      }
   }
   
   # Char main.lab cannot be NULL. 
   if ("main.lab" %in% args) {
      main.lab <- l$main.lab
      if (!is.vector(main.lab) || !is.character(main.lab) || length(main.lab) != 1) {
         print.default("main.lab is:"); print.default(main.lab, max = 100)
         stop("main.lab must be a single character string.", call. = FALSE)
      }
   }   
   
   # lenght-n real vector val can be NULL; in that case it is set to c(1,...,1)
   if ("val" %in% args) {
      val <- l$val
      if (!is.null(val)) {
         if (!is.vector(val) || !is.numeric(val) || !all(is.finite(val)) || length(val) != n) {
            print.default("val is:"); print.default(val, max = 100)
            stop("val must be NULL or a non-negative, nonzero, finite numeric vector of length nrow(X).", call. = FALSE)
         }
      }
   }
   
   # pool.fun is a vector of strings that cannot be empty
   if ("pool.fun" %in% args) {
      pool.fun <- l$pool.fun
      if (!is.vector(pool.fun) || !is.character(pool.fun) || length(pool.fun) < 1 ||
          !all(pool.fun %in% c("sum", "max", "min", "mean", "median", "0"))) {
         print.default("pool.fun is:"); print.default(pool.fun, max = 100)
         stop("pool.fun must a vector of strings containing 'sum'. All components of pool.fun must be selected from 'sum', 'min', 'mean', 'median', and '0'.", call. = FALSE)
      }
   }
   
   # Boolean return.pools cannot be NULL
   if ("return.pools" %in% args) {
      return.pools <- l$return.pools
      if (!is.vector(return.pools) || !is.logical(return.pools) || length(return.pools) > 1) {
         print.default("return.pools is:"); print.default(return.pools, max = 100)
         stop("return.pools must be a single logical value (TRUE or FALSE).", call. = FALSE)
      }
   }
   
}


