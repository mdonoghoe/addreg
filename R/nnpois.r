utils::globalVariables("tol")

nnpois <- function(y, x, standard, offset, start, control = addreg.control(),
                   accelerate = c("em", "squarem", "pem", "qn"),
                   control.accelerate = list(list())) 
{
  control <- do.call("addreg.control", control)
  accelerate <- match.arg(accelerate)
  
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if(is.matrix(y))
    rownames(y)
  else names(y)
  
  if(any(x < 0)) stop("x must be non-negative")
  if(any(apply(x, 2, function(col) all(col==0)))) stop("x contains column with all 0")
  
  nvars <- ncol(x)
  nobs <- NROW(y)
  
  fam <- poisson(link = identity)
  eval(fam$initialize)
  
  mu.eta <- fam$mu.eta
  linkinv <- fam$linkinv
  dev.resids <- fam$dev.resids
  aic <- fam$aic
  
  weights <- rep(1, nobs)
  if (is.null(standard)) standard <- rep.int(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)
  if (any(offset < 0))
    stop("offset must be non-negative")
  
  converged <- FALSE
  
  coefold <- if(!is.null(start)) {
    if (length(start) != nvars)
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                nvars, paste(deparse(xnames), collapse = ", ")),
            domain = NA)
    else if(any(start <= control$bound.tol))
        stop("'start' is on our outside the boundary of the parameter space (consider 'bound.tol')", domain = NA)
    else start
  } else {
    simple <- mean(y / standard) / colMeans(x) + 2*control$bound.tol
    trymat <- tryCatch(as.vector(solve(t(x)%*%x) %*% t(x) %*% (y)) + 2*control$bound.tol,
                       error = function(e) NULL)
    if(is.null(trymat)) simple
    else if(any(trymat < control$bound.tol)) simple
    else trymat
  }
  
  fixptfn <- function(p, y, n, x, o, div, fam, bound.tol) {
    eta <- drop(x %*% p) + o
    y.over.fits <- y / fam$linkinv(eta)
    y.over.fits[fam$linkinv(eta) == 0] <- 0
    pnew <- p * colSums(y.over.fits * x) * div
    pnew[pnew <= 0] <- bound.tol / 2
    return(pnew)
  }
  
  objfn <- function(p, y, n, x, o, div, fam, bound.tol) {
    eta <- drop(x %*% p) + o
    mu <- n * fam$linkinv(eta)
    negll <- -sum(dpois(y, mu, log = TRUE))
    return(negll)
  }
  
  validparams <- function(p) return(all(p >= 0))
  
  conv.user <- function(old, new) return(conv.test(old, new, tol))
  
  std.div <- 1 / colSums(standard * x)
  
  res <- turboEM::turboem(par = coefold, fixptfn = fixptfn, objfn = objfn, method = accelerate,
                          pconstr = validparams, y = y, n = standard, x = x, o = offset,
                          div = std.div, fam = fam, bound.tol = control$bound.tol,
                          control.run = list(convtype = "parameter", tol = control$epsilon,
                                             stoptype = "maxiter", maxiter = control$maxit,
                                             convfn.user = conv.user, trace = control$trace),
                          control.method = control.accelerate)
  if (res$fail[1]) stop(res$errors[1])
  coefnew <- res$pars[1,]
  names(coefnew) <- xnames

  eta <- drop(x %*% coefnew) + offset
  mu <- standard * linkinv(eta)
  residuals <- (y - mu) / mu.eta(eta)
  
  names(y) <- names(mu) <- names(eta) <- names(residuals) <- ynames
  
  dev.new <- sum(dev.resids(y, mu, weights))  
  aic.model <- aic(y, nobs, mu, weights, dev.new) + 2 * nvars
  aic.c.model <- aic.model + 2 * nvars * (nvars + 1) / (nobs - nvars - 1)
  
  wtdmu <- standard * sum(weights * y / standard) / sum(weights)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  nulldf <- nobs - 1
  resdf <- nobs - nvars
  
  boundary <- any(coefnew < control$bound.tol)
  
  list(coefficients = coefnew, residuals = residuals, fitted.values = mu, rank = nvars, 
       family = fam, linear.predictors = eta, deviance = dev.new, aic = aic.model, 
       aic.c = aic.c.model, null.deviance = nulldev, iter = res$itr[1], weights = weights,
       prior.weights = weights, standard = standard, df.residual = resdf, df.null = nulldf, 
       y = y, converged = res$convergence[1], boundary = boundary, 
       loglik = -res$value.objfn[1], nn.design = x)
  
}