addreg.cem <- function(mt, mf, Y, standard, offset, mono, family, start, control, model,
                       accelerate, control.method, warn) {
  control2 <- control
  control2$trace <- (control$trace > 1)
  
  if (family$family == "poisson") method <- "nnpois"
  else if (substr(family$family,1,7) == "negbin1") method <- "nnnegbin"
  else if (substr(family$family,1,6) == "gamma1") method <- "nngamma"
  else if (family$family == "binomial") method <- "addbin"
  
  allref <- addreg.allref(mt, mf, "cem", mono, family, start)
  design.numref <- sapply(allref$allref, length)
  
  best.model <- best.param <- NULL
  best.loglik <- -Inf
  allconv <- TRUE
  totaliter <- 0
  
  if (length(allref$allref) == 0) {
    if (control$trace > 0) cat(method, "parameterisation 1/1\n")
    X <- model.matrix(allref$terms, allref$data)
    if (family$family == "poisson")
      best.model <- nnpois(Y, X, standard, offset, allref$start.new, control2,
                           accelerate, control.method)
    else if (substr(family$family, 1, 7) == "negbin1")
      best.model <- nnnegbin(Y, X, standard, offset, allref$start.new, control2,
                             accelerate, control.method)
    else if (substr(family$family, 1, 6) == "gamma1")
      best.model <- nngamma(Y, X, offset, allref$start.new, control2, accelerate,
                            control.method)
    else
      best.model <- addbin(Y, X, allref$start.new, control, allref, model,
                           "cem", accelerate, control.method)
    best.loglik <- best.model$loglik
    best.param <- 0
    allconv <- best.model$converged
    totaliter <- totaliter + best.model$iter
    if (control$trace > 0 & control$trace <= 1) {
      if (substr(family$family, 1, 7) == "negbin1" | substr(family$family, 1, 6) == "gamma1")
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
      X <- addreg.design(allref$terms, allref$data, "cem", allref$allref,
                         allref$monotonic, design.all[param,])
      if (family$family == "poisson")
        thismodel <- nnpois(Y, X, standard, offset, if (param == 1) allref$start.new else NULL,
                            control2, accelerate, control.method)
      else if (substr(family$family, 1, 7) == "negbin1")
        thismodel <- nnnegbin(Y, X, standard, offset, if (param == 1) allref$start.new else NULL,
                              control2, accelerate, control.method)
      else if (substr(family$family, 1, 6) == "gamma1")
        thismodel <- nngamma(Y, X, offset, if (param == 1) allref$start.new else NULL,
                             control2, accelerate, control.method)
      else if (family$family == "binomial")
        thismodel <- addbin(Y, X, if (param == 1) allref$start.new else NULL, control, allref,
                            model, "cem", accelerate, control.method)
      if (!thismodel$converged) allconv <- FALSE
      totaliter <- totaliter + thismodel$iter
      if (control$trace > 0 & control$trace <= 1)
        if (substr(family$family, 1, 7) == "negbin1" | substr(family$family, 1, 6) == "gamma1")
          cat("Log-likelihood =", thismodel$loglik, "Iterations -", thismodel$iter, "\n")
        else if (method != "addbin")
          cat("Deviance =", thismodel$deviance, "Iterations -", thismodel$iter, "\n")
      if (thismodel$loglik > best.loglik) {
        best.model <- thismodel
        best.loglik <- thismodel$loglik
        best.param <- param
        if (thismodel$converged & !thismodel$boundary) break
      }
    }
  }
  
  if (length(allref$allref) == 0) {
    nn.coefs <- coefs <- coefs.boundary <- best.model$coefficients
    nn.design <- design <- model.matrix(allref$terms, allref$data)
  } else {
    nn.coefs <- best.model$coefficients
    nn.design <- addreg.design(allref$terms, allref$data, "cem", allref$allref, allref$monotonic,
                               design.all[best.param,])
    reparam <- addreg.reparameterise(nn.coefs, mt, mf, "cem", allref$allref, allref$monotonic, 
                                     design.all[best.param,])
    coefs <- reparam$coefs
    design <- reparam$design
    coefs.boundary <- reparam$coefs.boundary
  }
  
  boundary <- any(coefs.boundary < control$bound.tol)
  
  if (warn) {
    if (!best.model$converged |  (!allconv & best.model$boundary))
      if (identical(accelerate, "em"))
        warning(gettextf("%s: algorithm did not converge within %d iterations -- increase 'maxit'.",
                         method, control$maxit), call. = FALSE)
      else
        warning(gettextf("%s(%s): algorithm did not converge within %d iterations -- increase 'maxit' or try with 'accelerate = \"em\"'.",
                         method, accelerate, control$maxit), call. = FALSE)
    if (boundary) {
      if (coefs.boundary[1] < control$bound.tol) {
        if (family$family == "poisson" || substr(family$family, 1, 7) == "negbin1")
          warning(gettextf("%s: fitted rates numerically 0 occurred", method), call. = FALSE)
        else if (family$family == "binomial")
          warning(gettextf("%s: fitted probabilities numerically 0 or 1 occurred", method), call. = FALSE)
      } else warning(gettextf("%s: MLE on boundary of constrained parameter space", method), call. = FALSE)
    }
  }
  
  fit <- list(coefficients = coefs)
  if (substr(family$family, 1, 7) == "negbin1" | substr(family$family, 1, 6) == "gamma1") fit$scale <- best.model$scale
  
  fit2 <- list(residuals = best.model$residuals, fitted.values = best.model$fitted.values,
               rank = best.model$rank, family = best.model$family,
               linear.predictors = best.model$linear.predictors, 
               deviance = best.model$deviance, loglik = best.model$loglik, aic = best.model$aic, 
               aic.c = best.model$aic.c, null.deviance = best.model$null.deviance, 
               iter = c(totaliter, best.model$iter), prior.weights = best.model$prior.weights,
               df.residual = best.model$df.residual, df.null = best.model$df.null,
               y = best.model$y, x = design, standard = standard, offset = offset)
  if (model) {
    fit2$model <- mf
    if (family$family == "binomial") fit2$model.addpois <- best.model$model.addpois
  }
  
  fit3 <- list(converged = best.model$converged, boundary = boundary, nn.coefficients = nn.coefs,
               nn.x = nn.design)
  fit <- c(fit, fit2, fit3)
  fit
}