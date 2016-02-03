addreg.em <- function(mt, mf, Y, standard, offset, mono, family, start, control, model,
                      accelerate, control.method, warn) {
  control2 <- control
  control2$trace <- (control$trace > 1)
  
  if (family$family == "poisson") method <- "nnpois"
  else if (substr(family$family,1,7) == "negbin1") method <- "nnnegbin"
  else if (family$family == "binomial") method <- "addbin"
  
  allref <- addreg.allref(mt, mf, "em", mono, family, start)
  design.numref <- sapply(allref$allref, length)
  
  if (control$trace > 0) cat(method, "parameterisation 1/1\n")
  if (length(allref$allref) == 0)
    X <- model.matrix(allref$terms, allref$data)
  else {
    design.all <- expand.grid(lapply(design.numref, seq_len))
    X <- addreg.design(allref$terms, allref$data, if (family$family == "binomial") "cem" else "em", 
                       allref$allref, allref$monotonic, design.all[1,])
  }
  
  if (family$family == "poisson")
    thismodel <- nnpois(Y, X, standard, offset, allref$start.new, control2, accelerate,
                        list(control.method))
  else if (substr(family$family, 1, 7) == "negbin1")
    thismodel <- nnnegbin(Y, X, standard, offset, allref$start.new, control2, accelerate,
                          list(control.method))
  else if (family$family == "binomial")
    thismodel <- addbin(Y, X, allref$start.new, control, allref, model, "em", accelerate, control.method)
    
  if (control$trace > 0 & control$trace <= 1)
    if (substr(family$family, 1, 7) == "negbin1")
      cat("Log-likelihood =", thismodel$loglik, "Iterations -", thismodel$iter, "\n")
    else if (method != "addbin")
      cat("Deviance =", thismodel$deviance, "Iterations -", thismodel$iter, "\n")

  if (length(allref$allref) == 0) {
    nn.coefs <- coefs <- coefs.boundary <- thismodel$coefficients
    nn.design <- design <- X
  } else {
    nn.coefs <- thismodel$coefficients
    nn.design <- X
    reparam <- addreg.reparameterise(nn.coefs, mt, mf, if (family$family == "binomial") "cem" else "em",
                                     allref$allref, allref$monotonic, design.all[1,])
    coefs <- reparam$coefs
    design <- reparam$design
    coefs.boundary <- reparam$coefs.boundary
  }
  
  nvars <- length(coefs) + as.numeric(substr(family$family, 1, 7) == "negbin1")
  vardiff <- length(nn.coefs) + as.numeric(substr(family$family, 1, 7) == "negbin1") - nvars
  aic.c <- thismodel$aic - 2 * vardiff + 2 * nvars * (nvars + 1) / (NROW(Y) - nvars - 1)
  
  boundary <- any(coefs.boundary < control$bound.tol)
  
  if (warn) {
    if (!thismodel$converged) {
      if (identical(accelerate, "em"))
        warning(gettextf("%s: algorithm did not converge within %d iterations -- increase 'maxit'.",
                         method, control$maxit), call. = FALSE)
      else
        warning(gettextf("%s(%s): algorithm did not converge within %d iterations -- increase 'maxit' or try with 'accelerate = \"em\"'.",
                         method, accelerate, control$maxit), call. = FALSE)
    }
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
  if (substr(family$family, 1, 7) == "negbin1") fit$scale <- thismodel$scale
  
  fit2 <- list(residuals = thismodel$residuals, fitted.values = thismodel$fitted.values,
               rank = nvars, family = thismodel$family, linear.predictors = thismodel$linear.predictors, 
               deviance = thismodel$deviance, loglik = thismodel$loglik, 
               aic = thismodel$aic - 2*vardiff,  aic.c = aic.c, 
               null.deviance = thismodel$null.deviance,  iter = thismodel$iter, 
               prior.weights = thismodel$prior.weights, df.residual = thismodel$df.residual + vardiff,
               df.null = thismodel$df.null, y = thismodel$y, x = design, 
               standard = standard, offset = offset)
  if (model) {
    fit2$model <- mf
    if (family$family == "binomial") fit2$model.addpois <- thismodel$model.addpois
  }
  
  fit3 <- list(converged = thismodel$converged, boundary = boundary, nn.coefficients = nn.coefs,
               nn.x = nn.design)
  fit <- c(fit, fit2, fit3)
  fit
}