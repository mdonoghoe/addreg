gamma1 <- function (link, phi = stop("'phi' must be given"))
{
  .Phi <- phi
  env <- new.env(parent = .GlobalEnv)
  assign(".Phi", phi, envir = env)
  lnk <- make.link("identity")
  variance <- function(mu) mu * .Phi
  validmu <- function(mu) all(mu >= 0)
  dev.resids <- function(y, mu, wt) {
    sy <- .Phi * distr::igamma(log(y / .Phi))
    2 * ((sy - mu)/.Phi * log(y/.Phi) - lgamma(sy/.Phi) + lgamma(mu/.Phi))
  }
  aic <- function(y, n, mu, wt, dev) -2 * sum(dgamma(y, shape = mu/.Phi, scale = .Phi, log = TRUE))
  initialize <- expression({
    if (any(y <= 0)) stop("non-positive values not allowed for the 'gamma1' family")
    n <- rep.int(1, nobs)
  })
  simfun <- function(object, nsim) {
    ftd <- fitted(object)
    val <- rgamma(nsim * length(ftd), shape = ftd/.Phi, scale = .Phi)
    val
  }
  fname <- paste("gamma1 (phi = ", format(round(phi, 4)), ")", sep = "")
  environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
  structure(list(family = fname, link = "identity", linkfun = lnk$linkfun,
                 linkinv = lnk$linkinv, variance = variance, dev.resids = dev.resids,
                 aic = aic, mu.eta = lnk$mu.eta, initialize = initialize,
                 valid.mu = validmu, valid.eta = lnk$valideta, simulate = simfun),
            class = "family")
}