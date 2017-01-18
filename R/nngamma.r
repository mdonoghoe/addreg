utils::globalVariables("tol")

nngamma <- function(y, x, offset, start, control = addreg.control(),
                    accelerate = c("em", "squarem", "pem", "qn"),
                    control.method = list())
{
  control <- do.call("addreg.control", control)
  accelerate <- match.arg(accelerate)
  
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y))
    rownames(y)
  else names(y)
  
  if (any(x < 0)) stop("x must be non-negative")
  if (any(apply(x, 2, function(col) all(col == 0)))) stop("x contains column with all 0")
  
  nvars <- ncol(x)
  nobs <- NROW(y)
  
  fam <- gamma1(link = identity, phi = NA)
  eval(fam$initialize)
  
  weights <- rep(1, nobs)
  if (!is.null(offset))
    stop("offset is not currently supported for gamma1 family")
  
  converged <- FALSE
  
  if (!is.null(start)) {
    if (length(start) != nvars + 1)
      stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                    nvars + 1, paste(c(deparse(xnames), "scale"), collapse = ", ")),
           domain = NA)
    else if (any(start[-length(start)] <= control$bound.tol))
      stop("'start' is on our outside the boundary of the parameter space (consider 'bound.tol')",
           domain = NA)
    else if (start[length(start)] <= 0)
      stop("'start' value for scale must be > 0 for the gamma1 family", domain = NA)
    else {
      coefold.mu <- start[-length(start)]
      coefold.phi <- start[length(start)]
    }
  } else {
    simple <- mean(y) / colMeans(x) + 2*control$bound.tol
    trymat <- tryCatch(as.vector(solve(t(x)%*%x) %*% t(x) %*% (y)) + 2*control$bound.tol, error = function(e) NULL)
    if (is.null(trymat) || any(trymat < control$bound.tol)) coefold.mu <- simple
    else coefold.mu <- trymat
    coefold.phi <- 1
  }
  
  estep <- function(theta, x, y, a.fit) {
    x.theta <- theta * x
    elogbeta <- x.theta
    elogbeta[x.theta == 0] <- 1
    elogbeta[x.theta != 0] <- digamma(x.theta[x.theta != 0]) - digamma(a.fit[x.theta != 0])
    elogbeta + log(y)
  }
  
  score <- function(theta, x, t.fit, phi) {
    sum(x[x > 0] * (t.fit[x > 0] - log(phi) - digamma(theta * x[x > 0])))
  }
  
  fixptfn <- function(p, y, x, score, bound.tol, epsilon) {
    c.phi <- p[1]
    c.mu <- p[-1]
    c.a <- c.a.old <- c.mu / c.phi
    a.fit <- drop(x %*% c.a)
    for (j in seq_along(c.a)) {
      if (c.a[j] > bound.tol) {
        t.fit <- estep(c.a.old[j], x[,j], y, a.fit)
        pos <- x[,j] > 0
        meanxpos <- mean(x[pos,j])
        sumxpos <- sum(x[pos,j])
        sumtpos <- sum(x[pos,j] * (t.fit[pos] - log(c.phi) - log(x[pos,j])))
        if (is.infinite(exppart <- exp(-sumtpos / sumxpos))) a.max <- bound.tol
        else a.max <- max(1 / (meanxpos * gsl::lambert_W0(1 / meanxpos * exppart)), bound.tol)
        score.0 <- score(bound.tol / 2, x = x[,j], t.fit = t.fit, phi = c.phi)
        score.max <- score(a.max, x = x[,j], t.fit = t.fit, phi = c.phi)
        if (score.0 <= epsilon) c.a[j] <- 0
        else if (score.max >= 0) c.a[j] <- a.max
        else
          c.a[j] <- uniroot(score, interval = c(bound.tol / 2, a.max), x = x[,j],
                            t.fit = t.fit, phi = c.phi, f.lower = score.0,
                            f.upper = score.max, tol = epsilon * 1e-2)$root
      }
    }
    a.fit <- drop(x %*% c.a)
    c.phi <- sum(y) / sum(a.fit)
    c.mu <- c.phi * c.a
    pnew <- c(c.phi, c.mu)
    return(pnew)
  }
  
  objfn <- function(p, y, x, score, bound.tol, epsilon) {
    fam <- gamma1(link = identity, phi = p[1])
    eta <- drop(x %*% p[-1])
    mu <- fam$linkinv(eta)
    nobs <- NROW(y)
    wts <- rep(1, nobs)
    dev <- sum(fam$dev.resids(y, mu, wts))
    negll <- fam$aic(y, nobs, mu, wts, dev) / 2
    return(negll)
  }
  
  validparams <- function(p) return(all(p >= 0))
  
  projfn <- function(p) {
    p[p < 0] <- .Machine$double.eps
    p
  }
  
  conv.user <- function(old, new) return(conv.test(old[1], new[1], tol) &&
                                           conv.test(old[-1], new[-1], tol))
  
  res <- turboEM::turboem(par = c(coefold.phi, coefold.mu), fixptfn = fixptfn, objfn = objfn,
                          method = accelerate, pconstr = validparams, project = projfn, y = y, x = x, score = score,
                          bound.tol = control$bound.tol, epsilon = control$epsilon,
                          control.run = list(convtype = "parameter", tol = control$epsilon,
                                             stoptype = "maxiter", maxiter = control$maxit,
                                             convfn.user = conv.user, trace = control$trace),
                          control.method = list(control.method))
  if (res$fail[1]) stop(res$errors[1])
  coefnew.phi <- res$pars[1,1]
  coefnew.mu <- res$pars[1,-1]
  coefnew.a <- coefnew.mu / coefnew.phi
  
  fam <- gamma1(link = identity, phi = coefnew.phi)
  eta <- drop(x %*% coefnew.mu)
  mu <- fam$linkinv(eta)
  residuals <- (y - mu) / fam$mu.eta(eta)
  
  names(coefnew.mu) <- xnames
  names(coefnew.phi) <- NULL
  names(residuals) <- names(mu) <- names(eta) <- names(y) <- ynames
  
  dev.new <- sum(fam$dev.resids(y, mu, weights))
  aic.model <- fam$aic(y, nobs, mu, weights, dev.new) + 2 * (nvars + 1)
  aic.c.model <- aic.model + 2 * (nvars + 1) * (nvars + 2) / (nobs - nvars - 1)
  
  wtdmu <- rep(sum(weights * y) / sum(weights), nobs)
  nulldev <- sum(fam$dev.resids(y, wtdmu, weights))
  nulldf <- nobs - 1
  resdf <- nobs - nvars - 1
  
  boundary <- any(coefnew.a < control$bound.tol)
  
  list(coefficients = coefnew.mu, scale = coefnew.phi, residuals = residuals,
       fitted.values = mu, rank = nvars + 1, family = fam, linear.predictors = eta, 
       deviance = dev.new, aic = aic.model, aic.c = aic.c.model, null.deviance = nulldev,
       iter = res$itr[1], weights = weights, prior.weights = weights, 
       df.residual = resdf, df.null = nulldf, y = y, converged = res$convergence[1], 
       boundary = boundary, loglik = -res$value.objfn[1], nn.design = x)
}