addreg.cem <- function(mt, mf, Y, offset, mono, family, start, control, accelerate,
                       control.method, warn) {
  control2 <- control
  control2$trace <- (control$trace - 1)
  
  if (family$family == "poisson") method <- "nnpois"
  else if (substr(family$family,1,7) == "negbin1") method <- "nnnegbin"
  else if (family$family == "binomial") method <- "addbin"
  
  allref <- addreg.allref(mt, mf, "cem", mono, family, start)
  design.numref <- sapply(allref$allref, length)
  
  best.model <- best.param <- NULL
  best.loglik <- -Inf
  allconv <- TRUE
  totaliter <- 0
  
  if (length(allref$allref) == 0) {
    if (control$trace > 0) cat(method, " parameterisation 1/1\n")
    X <- model.matrix(allref$terms, allref$data)
    if (family$family == "poisson")
      best.model <- nnpois(Y, X, standard, offset, allref$start.new, control2,
                           accelerate, list(control.method))
    else if (substr(family$family, 1, 7) == "negbin1")
      best.model <- nnnegbin(Y, X, standard, offset, allref$start.new, control2,
                             accelerate, list(control.method))
    else
      best.model <- addbin(Y, X, allref$start.new, control, allref, model,
                           accelerate, control.accelerate)
    best.loglik <- best.model$loglik
    best.param <- 0
    allconv <- best.model$converged
    totaliter <- totaliter + best.model$iter
    if (control$trace > 0 & control$trace <= 1) {
      if (substr(family$family, 1, 7) == "negbin1")
        cat("Log-likelihood =", best.model$loglik,
            "Iterations -", best.model$iter, "\n")
      else if (method != "addbin")
        cat("Deviance =", best.model$deviance,
            "Iterations -", best.model$iter, "\n")
    }
  } else {
    design.all <- expand.grid(lapply(design.numref, seq_len))
    nparam <- nrow(design.all)
    
    for (param in seq_len(nparam)) {
      if (control$trace > 0) cat(method, " parameterisation ", param, "/", nparam, "\n", sep = "")
      X <- addreg.design(allref$terms, allref$data, "cem", allref$allref, design.all[param,])
      if (family$family == "poisson")
        thismodel <- nnpois(Y, X, standard, offset, if (param == 1) allref$start.new else NULL,
                            control2, accelerate, list(control.method))
      else if (substr(family$family, 1, 7) == "negbin1")
        thismodel <- nnnegbin(Y, X, standard, offset, if (param == 1) allref$start.new else NULL,
                              control2, accelerate, list(control.method))
      else if (family$family == "binomial")
        thismodel <- addbin(Y, X, if (param == 1) allref$start.new else NULL, control, allref,
                            model, accelerate, control.accelerate)
      if (!thismodel$converged) allconv <- FALSE
      
    }
  }
}