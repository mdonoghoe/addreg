utils::globalVariables("tol")

nnnegbin <- function(y, x, standard, offset, start, control = addreg.control(),
                     accelerate = c("em", "squarem", "pem", "qn"),
                     control.accelerate = list(list()))
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
  
  fam <- negbin1(link = identity, phi = NA)
  eval(fam$initialize)
  
  weights <- rep(1, nobs)
  if (is.null(standard)) standard <- rep.int(1, nobs)
  if (!is.null(offset))
    stop("offset is not currently supported for negbin1 family")
    
  converged <- FALSE
  
  if (!is.null(start)) {
    if (length(start) != nvars + 1)
      stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
        nvars + 1, paste(c(deparse(xnames), "scale"), collapse = ", ")),
        domain = NA)
    else if (any(start[-length(start)] <= control$bound.tol))
      stop("'start' is on our outside the boundary of the parameter space (consider 'bound.tol')",
        domain = NA)
    else if (start[length(start)] <= 1)
      stop("'start' value for scale must be > 1 for the negbin1 family", domain = NA)
    else {
      coefold.mu <- start[-length(start)]
      coefold.phi <- start[length(start)] - 1
    }
  } else {
    simple <- mean(y / standard) / colMeans(x) + 2*control$bound.tol
    trymat <- tryCatch(as.vector(solve(t(x)%*%x) %*% t(x) %*% (y)) + 2*control$bound.tol, error = function(e) NULL)
    if (is.null(trymat) || any(trymat < control$bound.tol)) coefold.mu <- simple
    else coefold.mu <- trymat
    coefold.phi <- 0.5
  }
  
  dbetabinom <- function(x, size, alpha, beta) {
    exp(lchoose(size, x) + lbeta(alpha + x, size + beta - x) - lbeta(alpha, beta))
  }
  
  digammadiff <- function(y, r) {
    k <- seq_len(y) - 1
    sum(1 / (r + k))
  }
    
  digammaexp <- function(y, r, r.part, r.full) {
    k <- seq_len(y)
    rdiff <- cumsum(1 / (r + k - 1))
    prob <- sapply(k, dbetabinom, size = y, alpha = r.part, beta = r.full - r.part)
    sum(rdiff * prob)
  }
  
  estep <- function(theta, x, y, std, theta.old, r.fits.old) {
    res <- rep(0, length(y))
    r.part.old <- std * theta.old * x
    chk1 <- (y == 0)
    chk2 <- (r.part.old == r.fits.old)
    if(sum(!chk1 & chk2) > 0)
      res[!chk1 & chk2] <- mapply(digammadiff, y = y[!chk1 & chk2],
                    r = (std * theta * x)[!chk1 & chk2])
    if (sum(!chk1 & !chk2) > 0)
      res[!chk1 & !chk2] <- mapply(digammaexp, y = y[!chk1 & !chk2],
                                        r = (std * theta * x)[!chk1 & !chk2],
                                        r.part = r.part.old[!chk1 & !chk2],
                                        r.full = r.fits.old[!chk1 & !chk2])
    res
  }
    
  score <- function(theta, x, y, std, theta.old, r.fits.old, p) {
    sum(std[x > 0] * x[x > 0] * (estep(theta, x[x > 0], y[x > 0], std[x > 0], theta.old,
      r.fits.old[x > 0]) + log(1 - p)))
  }
  
  fixptfn <- function(p, y, n, x, score, bound.tol, epsilon) {
    c.phi <- p[1]
    c.mu <- p[-1]
    c.r <- c.mu / c.phi
    c.p <- c.phi / (c.phi + 1)
    r.fit <- n * drop(x %*% c.r)
    c.max <- drop(c.r * (t(ifelse(r.fit == 0, 0, n * y / r.fit)) %*% x)
                    / (log(1 / (1 - c.p)) * t(n) %*% x))
    for (j in 1L:length(c.r)) {
      if (c.r[j] > bound.tol) {
        score.0 <- score(bound.tol / 2, x = x[,j], y = y, std = n, theta.old = c.r[j],
                         r.fits.old = r.fit, p = c.p)
        score.max <- score(c.max[j], x = x[,j], y = y, std = n, theta.old = c.r[j],
                           r.fits.old = r.fit, p = c.p)
        if (score.0 <= epsilon)
          c.r[j] <- 0
        else if (score.max >= 0)
          c.r[j] <- c.max[j]
        else
          c.r[j] <- uniroot(score, interval = c(bound.tol / 2, c.max[j]), x = x[,j],
                            y = y, std = n, theta.old = c.r[j], r.fits.old = r.fit,
                            p = c.p, f.lower = score.0, f.upper = score.max,
                            tol = epsilon * 1e-2)$root
      }
    }
    r.fit <- n * drop(x %*% c.r)
    c.p <- sum(y) / (sum(y) + sum(r.fit))
    c.phi <- c.p / (1 - c.p)
    c.mu <- c.phi * c.r
    pnew <- c(c.phi, c.mu)
    return(pnew)
  }
  
  objfn <- function(p, y, n, x, score, bound.tol, epsilon) {
    fam <- negbin1(link = identity, phi = p[1])
    eta <- drop(x %*% p[-1])
    mu <- n * fam$linkinv(eta)
    nobs <- NROW(y)
    wts <- rep(1, nobs)
    dev <- sum(fam$dev.resids(y, mu, wts))
    negll <- fam$aic(y, nobs, mu, wts, dev) / 2
    return(negll)
  }
  
  validparams <- function(p) return(all(p >= 0) && p[1] <= 1)
  
  conv.user <- function(old, new) return(conv.test(old[1], new[1], tol) && 
                                         conv.test(old[-1], new[-1], tol))
                                         
  res <- turboEM::turboem(par = c(coefold.phi, coefold.mu), fixptfn = fixptfn, objfn = objfn,
                          method = accelerate, pconstr = validparams, y = y, n = standard,
                          x = x, score = score, bound.tol = control$bound.tol,
                          epsilon = control$epsilon,
                          control.run = list(convtype = "parameter", tol = control$epsilon,
                                             stoptype = "maxiter", maxiter = control$maxit,
                                             convfn.user = conv.user, trace = control$trace),
                          control.method = control.accelerate)
  if (res$fail[1]) stop(res$errors[1])
  coefnew.phi <- res$pars[1,1]
  coefnew.mu <- res$pars[1,-1]
  coefnew.r <- coefnew.mu / coefnew.phi
  coefnew.p <- coefnew.phi / (1 + coefnew.phi)
  
  fam <- negbin1(link = identity, phi = coefnew.phi)
  eta <- drop(x %*% coefnew.mu)
  mu <- standard * fam$linkinv(eta)
  residuals <- (y - mu) / fam$mu.eta(eta)
  
  names(coefnew.mu) <- xnames
  names(residuals) <- names(mu) <- names(eta) <- names(y) <- ynames
  
  dev.new <- sum(fam$dev.resids(y, mu, weights))
  aic.model <- fam$aic(y, nobs, mu, weights, dev.new) + 2 * (nvars + 1)
  aic.c.model <- aic.model + 2 * (nvars + 1) * (nvars + 2) / (nobs - nvars - 1)
  
  wtdmu <- standard * rep(sum(weights * y / standard) / sum(weights), nobs)
  nulldev <- sum(fam$dev.resids(y, wtdmu, weights))
  nulldf <- nobs - 1
  resdf <- nobs - nvars - 1
  
  boundary <- any(coefnew.r < control$bound.tol)
  
  list(coefficients = coefnew.mu, scale = 1 + coefnew.phi, residuals = residuals,
       fitted.values = mu, rank = nvars + 1, family = fam, linear.predictors = eta, 
       deviance = dev.new, aic = aic.model, aic.c = aic.c.model, null.deviance = nulldev,
       iter = res$itr[1], weights = weights, prior.weights = weights, standard = standard, 
       df.residual = resdf, df.null = nulldf, y = y, converged = res$convergence[1], 
       boundary = boundary, loglik = -res$value.objfn[1], nn.design = x)
}