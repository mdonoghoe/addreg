\name{addreg.smooth}
\alias{addreg.smooth}
\title{
Smooth Additive Regression for Discrete Data
}
\description{
\code{addreg.smooth} fits additive (identity-link) Poisson, negative binomial
and binomial regression models using a stable EM algorithm. It provides additional
flexibility over \code{\link{addreg}} by allowing for semi-parametric
terms.
}
\usage{
addreg.smooth(formula, mono = NULL, family, data, standard, subset, 
              na.action, offset, control = list(...), model = TRUE, 
              model.addreg = FALSE, ...)
}
\arguments{
  \item{formula}{
	an object of class \code{"\link{formula}"}
  (or one that can be coerced into that class): a symbolic
  description of the model to be fitted. The details of
  model specification are given under "Details".  The model 
  must contain an intercept and at least one semi-parametric
  term, included by using the \code{\link{B}} or \code{\link{Iso}}
  functions. Note that 2nd-order terms (such as interactions) or above 
  are not currently supported (see \code{\link{addreg}}).
}
  \item{mono}{
	a vector indicating which terms in
  \code{formula} should be restricted to have a
  monotonically non-decreasing relationship with the
  outcome. May be specified as names or indices of the
  terms.
  
	\code{\link{Iso}()} terms are always monotonic.
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
  from which \code{addreg.smooth} is called.
}
  \item{standard}{
	a numeric vector of length equal to the
  number of cases, where each element is a positive
  constant that (multiplicatively) standardises the fitted
  value of the corresponding element of the response
  vector. Ignored for binomial family (the two-column
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
  \item{offset}{
	this can be used to specify an \emph{a
  priori} known component to be included in the linear
  predictor during fitting. This should be \code{NULL} or a
  \emph{non-negative} numeric vector of length equal to the number of cases.
  One or more \code{\link{offset}} terms can be included in
  the formula instead or as well, and if more than one is
  specified their sum is used. See
  \code{\link{model.offset}}.
  
	Ignored for binomial family.
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
  \item{model.addreg}{
	a logical value indicating whether the fitted \code{addreg}
	object should be included as a component of the returned value.
}
  \item{\dots}{
	arguments to be used to form the default
  \code{control} argument if it is not supplied directly.
}
}
\details{
\code{addreg.smooth} performs the same fitting process as \code{\link{addreg}},
providing a stable maximum likelihood estimation procedure for identity-link
Poisson, negative binomial or binomial models, with the added flexibility of allowing semi-parametric
\code{\link{B}} and \code{\link{Iso}} terms (note that \code{addreg.smooth} will stop with an 
error if no semi-parametric terms are specified in the right-hand side of the \code{formula}; 
\code{\link{addreg}} should be used instead).

The method partitions the parameter space associated with the semi-parametric part of the
model into a sequence of constrained parameter spaces, and defines a fully parametric 
\code{addreg} model for each. The model with the highest log-likelihood is the MLE for 
the semi-parametric model (see Donoghoe and Marschner, 2015).
}
\value{
An object of class \code{"addreg.smooth"}, which contains the same objects as class
\code{"addreg"} (the same as \code{"\link{glm}"} objects, without \code{contrasts},
\code{qr}, \code{R} or \code{effects} components), as well as:
\item{model.addreg}{if \code{model.addreg} is \code{TRUE}; the \code{addreg} object
for the fully parametric model corresponding to the fitted model.}
\item{xminmax.smooth}{the minimum and maximum observed values for
each of the smooth terms in the model, to help define the 
covariate space.}
\item{full.formula}{the component from \code{\link{interpret.addreg.smooth}(formula)}
that contains the \code{formula} term with any additional arguments to the \code{\link{B}}
function removed.}
\item{knots}{a named list containing the knot vectors for each of the smooth
terms in the model.}
}
\references{
Donoghoe, M. W. and I. C. Marschner (2015). Flexible regression models for rate 
differences, risk differences and relative risks. \emph{International Journal of 
Biostatistics} 11(1): 91--108.

Marschner, I. C. (2014). Combinatorial EM algorithms.
\emph{Statistics and Computing} 24(6): 921--940.
}
\author{
Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}
}
\seealso{
\code{\link{addreg}}
}
\examples{
## Simple example
dat <- data.frame(x1 = c(3.2,3.3,3.4,7.9,3.8,0.7,2.0,5.4,8.4,3.0,1.8,5.6,5.5,9.0,8.2),
  x2 = c(1,0,0,1,0,1,0,0,0,0,1,0,1,1,0),
  n = c(6,7,5,9,10,7,9,6,6,7,7,8,6,8,10),
  y = c(2,1,2,6,3,1,2,2,4,4,1,2,5,7,7))
m1 <- addreg.smooth(cbind(y, n-y) ~ B(x1, knot.range = 1:3) + factor(x2), mono = 1,
  data = dat, family = binomial, trace = 1)
\donttest{
plot(m1, at = data.frame(x2 = 0:1))
points(dat$x1, dat$y / dat$n)
}
}
\keyword{ regression }
\keyword{ smooth }
