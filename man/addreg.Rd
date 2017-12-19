\name{addreg}
\alias{addreg}
\title{
Additive Regression for Discrete Data
}
\concept{Binomial regression}
\concept{Poisson regression}
\concept{CEM algorithm}
\concept{Negative binomial regression}
\description{
\code{addreg} fits additive (identity-link) Poisson, negative binomial
 and binomial regression models using a stable combinatorial EM algorithm.
}
\usage{
addreg(formula, mono = NULL, family, data, standard, subset, na.action,
       start = NULL, offset, control = list(...), model = TRUE, 
       method = c("cem", "em"),
       accelerate = c("em", "squarem", "pem", "qn"), 
       control.method = list(), warn = TRUE, ...)
}
\arguments{
  \item{formula}{
	an object of class \code{"\link{formula}"}
  (or one that can be coerced into that class): a symbolic
  description of the model to be fitted. The details of
  model specification are given under "Details".  Note that
  the model must contain an intercept, and 2nd-order terms
  (such as interactions) or above are currently not supported
  --- see "Note".
}
  \item{mono}{
	a vector indicating which terms in
  \code{formula} should be restricted to have a
  monotonically non-decreasing relationship with the
  outcome. May be specified as names or indices of the
  terms.
}
  \item{family}{
	a description of the error distribution to
  be used in the model. This can be a character string
  naming a family function, a family function or the result
  of a call to a family function (see
  \code{\link[stats]{family}} for details of family
  functions), but here it is restricted to be \code{\link{poisson}},
  \code{\link{negbin1}} or \code{\link{binomial}} family with \code{identity} link.
}
  \item{data}{
	an optional data frame, list or environment
  (or object coercible by \code{\link{as.data.frame}} to a
  data frame) containing the variables in the model. If not
  found in \code{data}, the variables are taken from
  \code{\link{environment}(formula)}, typically the environment
  from which \code{addreg} is called.
}
  \item{standard}{
	a numeric vector of length equal to the
  number of cases, where each element is a positive
  constant that (multiplicatively) standardises the fitted
  value of the corresponding element of the response
  vector. Ignored for binomial family (two-column
  specification of response should be used instead).
}
  \item{subset}{
	an optional vector specifying a subset of
  observations to be used in the fitting process.
}
  \item{na.action}{
  a function which indicates what should happen when the data
  contain \code{NA}s. The default is set be the \code{na.action}
  setting of \code{\link{options}}, and is \code{\link{na.fail}}
  if that is unset. The `factory-fresh' default is \code{\link{na.omit}}.
  Another possible value is \code{NULL}, no action. Value
  \code{\link{na.exclude}} can be useful.
}
  \item{start}{
	starting values for the parameters in the
  linear predictor, also with the starting value for
  the \code{scale} as the last element when
  \code{family = \link{negbin1}}.
}
  \item{offset}{
	this can be used to specify an \emph{a
  priori} known component to be included in the linear
  predictor during fitting. This should be \code{NULL} or a
  \emph{non-negative} numeric vector of length equal to the number of cases.
  One or more \code{\link{offset}} terms can be included in
  the formula instead or as well, and if more than one is
  specified their sum is used. See
  \code{\link{model.offset}}.
  
	Ignored for binomial family; not yet implemented for negative binomial models.
}
  \item{control}{
	list of parameters for controlling the
  fitting process, passed to
  \code{\link{addreg.control}}.
}
  \item{model}{
	a logical value indicating whether the \emph{model frame} 
	(and, for binomial models, the equivalent Poisson model)
	should be included as a component of the returned value.
}
  \item{method}{
  a character string that determines which algorithm to use to
  find the MLE: \code{"cem"} for the combinatorial EM algorithm,
  which cycles through a sequence of constrained parameter spaces,
  or \code{"em"} for a single EM algorithm based on an
  overparameterised model.
}
  \item{accelerate}{
    a character string that determines the acceleration
  algorithm to be used, (partially) matching one of \code{"em"} (no acceleration --- the default), 
  \code{"squarem"}, \code{"pem"} or \code{"qn"}. See \code{\link[turboEM]{turboem}}
  for further details. Note that \code{"decme"} is not permitted.
}
  \item{control.method}{
    a list of control parameters for the acceleration algorithm, which are passed to 
    the \code{control.method} argument of \code{\link[turboEM]{turboem}}.
  
    If any items are not specified, the defaults are used.
}
  \item{warn}{
	a logical indicating whether or not warnings should be
  provided for non-convergence or boundary values.
}
  \item{\dots}{
	arguments to be used to form the default
  \code{control} argument if it is not supplied directly.
}
}
\details{
\code{addreg} fits a generalised linear model (GLM) with a
Poisson or binomial error distribution and identity link
function, as well as additive NegBin I models (which are not 
GLMs). Predictors are assumed to be continuous, unless
they are of class \code{\link{factor}}, or are character or
logical (in which case they are converted to
\code{factor}s). Specifying a predictor as monotonic using
the \code{mono} argument means that for continuous terms,
the associated coefficient will be restricted to be
non-negative, and for categorical terms, the coefficients
will be non-decreasing in the order of the factor
\code{levels}. This allows semi-parametric monotonic regression
functions, in the form of unsmoothed step-functions. For 
smooth regression functions see \code{\link{addreg.smooth}}.

As well as allowing monotonicity constraints, the function 
is useful when a standard GLM routine, such as
\code{\link{glm}}, fails to converge with an identity-link Poisson
or binomial model. If \code{\link{glm}} does achieve successful convergence,
and \code{addreg} converges to an interior point, then the two
results will be identical. However, \code{\link{glm}} may still experience convergence 
problems even when \code{addreg} converges to an interior point. 
Note that if \code{addreg} converges to a boundary point, then it 
may differ slightly from \code{\link{glm}} even if \code{\link{glm}} successfully 
converges, because of differences in the definition of the parameter 
space. \code{addreg} produces valid fitted values for covariate 
values within the Cartesian product of the observed range of covariate 
values, whereas \code{\link{glm}} produces valid fitted values just 
for the observed covariate combinations (assuming it successfully 
converges). This issue is only relevant when \code{addreg} 
converges to a boundary point.

The computational method is a combinatorial EM algorithm
(Marschner, 2014), which accommodates the parameter
contraints in the model and is more stable than iteratively
reweighted least squares. A collection of restricted
parameter spaces is defined which covers the full parameter
space, and the EM algorithm is applied within each
restricted parameter space in order to find a collection of
restricted maxima of the log-likelihood function, from
which can be obtained the global maximum over the full
parameter space. See Marschner (2010), Donoghoe and 
Marschner (2014) and Donoghoe and Marschner (2016) for further details.

Acceleration of the EM algorithm can be achieved through the
methods of the \pkg{\link[turboEM:turboem]{turboEM}} package, specified
through the \code{accelerate} argument. However, note that these
methods do not have the guaranteed convergence of the standard
EM algorithm, particularly when the MLE is on the boundary of
its (possibly constrained) parameter space.
}
\value{
\code{addreg} returns an object of class \code{"addreg"},
which inherits from classes \code{"glm"} and \code{"lm"}.
The function \code{\link{summary.addreg}} can be used 
to obtain or print a summary of the results.

The generic accessor functions \code{\link{coefficients}},
\code{\link{fitted.values}} and \code{\link{residuals}} can be used to
extract various useful features of the value returned by
\code{addreg}. Note that \code{\link{effects}} will not work.

An object of class \code{"addreg"} is a list containing the
same components as an object of class \code{"glm"} (see the
"Value" section of \code{\link{glm}}), but without
\code{contrasts}, \code{qr}, \code{R} or \code{effects}
components. It also includes: 
\item{loglik}{the maximised log-likelihood.}
\item{aic.c}{a small-sample corrected
version of Akaike's \emph{An Information Criterion}
(Hurvich, Simonoff and Tsai, 1998). This is used by
\code{\link{addreg.smooth}} to choose the optimal number of
knots for smooth terms.} 
\item{xminmax}{the minimum and maximum observed values
for each of the continuous covariates, to help define
the covariate space of the model.}

As well as, for Poisson and negative binomial models:
\item{nn.coefficients}{estimated coefficients associated
with the non-negative parameterisation corresponding to the
MLE.} 
\item{nn.x}{non-negative model matrix associated with
\code{nn.coefficients}.} 
\item{standard}{the \code{standard} argument.} 

Or, for binomial models:
\item{model.addpois}{if requested, the \code{addreg} object for the associated identity-link
Poisson model.}

The \code{scale} component of the result is fixed at \eqn{1} for
Poisson and binomial models, and is the constant overdispersion parameter
for negative binomial models (that is, \code{scale} = \eqn{1+\phi}) where 
\eqn{Var(\mu) = (1+\phi)\mu}).
}
\references{
Donoghoe, M. W. and I. C. Marschner (2014). Stable computational methods
for additive binomial models with application to adjusted risk differences. 
\emph{Computational Statistics and Data Analysis} 80: 184--196.

Donoghoe, M. W. and I. C. Marschner (2016). Estimation of adjusted rate
differences using additive negative binomial regression. \emph{Statistics
in Medicine} 35(18): 3166--3178.

Hurvich, C. M., J. S. Simonoff and C.-L. Tsai (1998). Smoothing parameter 
selection in nonparametric regression using an improved Akaike 
information criterion.
\emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)} 60(2): 271--293.

Marschner, I. C. (2010). Stable computation of maximum likelihood estimates
in identity link Poisson regression. \emph{Journal of Computational and 
Graphical Statistics} 19(3): 666--683.

Marschner, I. C. (2014). Combinatorial EM algorithms.
\emph{Statistics and Computing} 24(6): 921--940.
}
\author{
Mark W. Donoghoe \email{markdonoghoe@gmail.com}
}
\note{
Due to the way the covariate space is defined in the CEM algorithm,
specifying interactions in the formula is not currently supported 
by \code{addreg}. 2-way interactions between factors can be
included by calculating a new factor term that has levels
corresponding to all possible combinations of the factor
levels. See the Example.
}
\examples{
require(glm2)
data(crabs)

#============================================================================
#  identity-link Poisson model with periodic non-convergence when glm is used
#============================================================================

crabs.boot <- crabs[crabs$Rep1,-c(5:6)]
crabs.boot$width.shifted <- crabs.boot$Width - min(crabs$Width)

fit.glm <- glm(Satellites ~ width.shifted + factor(Dark) + factor(GoodSpine),
  family = poisson(identity), data = crabs.boot, start = rep(1,4), 
  control = glm.control(trace = TRUE))
  
fit.addreg <- addreg(formula(fit.glm), family = poisson, data = crabs.boot,
  trace = 1)

# Speed up convergence by using single EM algorithm
fit.addreg.em <- update(fit.addreg, method = "em")

# Speed up convergence by using acceleration methods
fit.addreg.acc <- update(fit.addreg, accelerate = "squarem")
fit.addreg.em.acc <- update(fit.addreg.em, accelerate = "squarem")
  
# Usual S3 methods work on addreg objects
summary(fit.addreg)
vcov(fit.addreg)
confint(fit.addreg)
summary(predict(fit.addreg), type = "response")

fit.addreg2 <- addreg(update(formula(fit.glm), ~ . - factor(GoodSpine)),
  family = poisson, data = crabs.boot, trace = 1)
anova(fit.addreg2, fit.addreg, test = "LRT")

# Account for overdispersion (use start to speed it up a little)
fit.addreg.od <- addreg(Satellites ~ factor(Dark) + factor(GoodSpine),
  family = negbin1, data = crabs.boot, trace = 1,
  start = c(4.3423675,-2.4059273,-0.4531984,5.969648))
summary(fit.addreg.od)
}
\keyword{ models }
\keyword{ regression }
