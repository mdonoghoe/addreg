negbin1 <- function (link, phi = stop("'phi' must be given"))
{
	.Phi <- phi
	env <- new.env(parent = .GlobalEnv)
	assign(".Phi", phi, envir = env)
	lnk <- make.link("identity")
	variance <- function(mu) mu * (1 + .Phi)
	validmu <- function(mu) all(mu >= 0)
	dev.resids <- function(y, mu, wt) {
	  y.un <- unique(y)
	  mu.un <- rep(NA, length(y.un))
	  mu.un[y.un == 0] <- 0
	  mu.un[y.un == 1] <- .Phi / log(1 + .Phi)
	  rootfn <- function(mu, y, phi) digamma(y + mu/phi) - digamma(mu/phi) - log(1 + phi)
	  getroot <- function(y, phi) uniroot(rootfn, c(ifelse(phi <= 1, phi*y, y), phi*y/log(1+phi)), y = y, phi = phi)$root
	  mu.un[y.un > 1] <- sapply(y.un[y.un > 1], getroot, phi = .Phi)
	  mu.s <- mu.un[match(y, y.un)]
		r <- (mu - mu.s)/.Phi * log(1 + .Phi)
		p <- which(y > 0)
		r[p] <- r[p] + lgamma(y[p] + mu.s[p]/.Phi) - lgamma(mu.s[p]/.Phi) - lgamma(y[p] + mu[p]/.Phi) + lgamma(mu[p]/.Phi)
		2 * r
	}
	aic <- function(y, n, mu, wt, dev) -2 * sum(dnbinom(y, size = mu/.Phi, prob = 1/(1+.Phi), log = TRUE))
	initialize <- expression({
		if (any(y < 0)) stop("negative values not allowed for the 'negbin1' family")
		if (any(abs(y - round(y)) > 0.001)) stop("non-integer counts in a negbin model")
		n <- rep.int(1, nobs)
	})
	simfun <- function(object, nsim) {
		ftd <- fitted(object)
		val <- rnbinom(nsim * length(ftd), size = ftd/.Phi, prob = 1/(1+.Phi))
		val
	}
    fname <- paste("negbin1 (phi = ", format(round(phi, 4)), ")", sep = "")
	environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
	structure(list(family = fname, link = "identity", linkfun = lnk$linkfun,
		linkinv = lnk$linkinv, variance = variance, dev.resids = dev.resids,
		aic = aic, mu.eta = lnk$mu.eta, initialize = initialize,
		valid.mu = validmu, valid.eta = lnk$valideta, simulate = simfun),
		class = "family")
}