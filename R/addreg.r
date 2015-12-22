addreg <- function (formula, mono = NULL, family, data, standard, subset, na.action, 
                    start = NULL, offset, control = list(...), model = TRUE, 
                    method = c("cem"),
                    accelerate = c("em", "squarem", "pem", "qn"), control.method = list(), 
                    warn = TRUE, ...) {
  call <- match.call()
  method <- match.arg(method)
  accelerate <- match.arg(accelerate)
  
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    if (identical(family, negbin1)) family <- family(link = "identity", phi = NA)
    else family <- family(link = "identity")
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
    
  if (family$link!="identity" | !(family$family %in% c("poisson","binomial") | substr(family$family,1,7) == "negbin1"))
    stop("family(link) must be one of: poisson(identity), binomial(identity), negbin1(identity)")
    
  if (missing(data)) data <- environment(formula)
  control <- do.call("addreg.control", control)
  
  outnames <- c("coefficients", "residuals", "fitted.values", "rank", "family",
                "linear.predictors", "deviance", "loglik", "aic", "aic.c",
                "null.deviance", "iter", "prior.weights", "weights",
                "df.residual", "df.null", "y", "x")
  if (model) {
    outcomes <- c(outnames, "model")
    if (family$family == "binomial") outnames <- c(outnames, "model.addpois")
  }
  outnames <- c(outnames, "converged", "boundary", "na.action", "call", "formula",
                "terms", "data", "standard", "offset", "control", "method", "xlevels",
                "xminmax", "nn.coefficients", "nn.x")
  fit <- sapply(outnames, function(x) NULL)
  
  if (method %in% c("glm", "glm2")) {
    stop(paste(method, "is not currently supported by addreg"))
  } else {
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "offset", "standard"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    
    if (is.empty.model(mt)) stop("empty model")
    if (attr(mt,"intercept") != 1) stop("models without intercept are not supported by addreg")
    if (any(attr(mt,"order") > 1)) stop("models with interactions are not supported by addreg")
    if (attr(mt,"response") == 0) stop("missing response")
    
    standard <- as.vector(model.extract(mf,"standard"))
    if(!is.null(standard)) {
      if (length(standard) != NROW(Y))
        stop(gettextf("number of values in 'standard' is %d should equal %d (number of observations)",
                      length(standard), NROW(Y)), domain = NA)
      if (family$family == "binomial")
        warning("'standard' is not used for binomial family", call. = FALSE)
      if (any(standard <= 0))
        stop("standard must be positive")
    }
    
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
      if (length(offset) != NROW(Y))
        stop(gettextf("number of values in 'offset' is %d should equal %d (number of observations)",
                      length(offset), NROW(Y)), domain = NA)
      if (family$family == "binomial")
        warning("'offset' is not supported for binomial family", call. = FALSE)
    }
    
    Y <- model.response(mf, "numeric")
    if (length(dim(Y)) == 1L) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if(!is.null(nm)) names(Y) <- nm
    }
    
    addreg.method <- paste("addreg", method, sep = ".")
    addreg.args <- list(mt = mt, mf = mf, Y = Y, standard = standard, offset = offset, mono = mono,
                        family = family, start = start, control = control)
    if (method %in% c("cem", "em")) addreg.args$accelerate <- accelerate
    addreg.args <- c(addreg.args, list(control.method = control.method, warn = warn))
    res <- do.call(addreg.method, addreg.args)
    mres <- match(outnames, names(res), 0L)
    fit[names(res)[mres]] <- res[mres]
    fit$family <- family
    fit$weights <- rep(1, NROW(Y))
    if (model) fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    fit$terms <- mt
    fit$data <- data
    fit$standard <- standard
    fit$offset <- offset
    fit$xlevels <- .getXlevels(mt, mf)
    fit$xminmax <- .getXminmax(mt, mf)
  }
    
  fit$call <- call
  fit$formula <- formula
  fit$control <- control
  fit$method <- method
  
  class(fit) <- c("addreg","glm","lm")
  fit
}