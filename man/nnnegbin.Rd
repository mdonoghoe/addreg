\name{nnnegbin}
\alias{nnnegbin}
\title{
ECME Algorithm for Additive Negative Binomial 1 Model
}
\description{
Finds the maximum likelihood estimate of an additive negative binomial (NB1) model
using an ECME algorithm,
where each of the mean coefficients is restricted to be non-negative.
}
\usage{
nnnegbin(y, x, standard, offset, start, control = addreg.control(),
         accelerate = c("em", "squarem", "pem", "qn"),
         control.accelerate = list(list()))
}
\arguments{
  \item{y}{
	non-negative integer response vector.
}
  \item{x}{
	non-negative covariate matrix.
}
  \item{standard}{
	standardising vector, where each element is a positive constant that
	(multiplicatively) standardises the fitted value of the corresponding
	element of the response vector. The default is a vector of ones.
}
  \item{offset}{
	non-negative additive offset vector. The default is a vector of zeros.
}
  \item{start}{
	vector of starting values for the parameter estimates. The last element is
  the starting value of the \code{scale}, and must be \code{> 1}. The remaining
  elements are for the additive mean parameters, and must be
	greater than \code{control$bound.tol}.
}
  \item{control}{
	an \code{\link{addreg.control}} object, which controls the fitting process.
}
  \item{accelerate}{
    a character string that determines the acceleration
  algorithm to be used, (partially) matching one of \code{"em"} (no acceleration -- the default), 
  \code{"squarem"}, \code{"pem"} or \code{"qn"}. See \code{\link[turboEM]{turboem}}
  for further details. Note that \code{"decme"} is not permitted.
}
  \item{control.accelerate}{
	a list of control parameters for the acceleration algorithm. See \code{\link[turboEM]{turboem}}
	for details of the parameters that apply to each algorithm. If not specified, the defaults are used.
}
}
\details{
This is a workhorse function for \code{\link{addreg}}, and runs the ECME algorithm to find the 
constrained non-negative MLE associated with an additive NB1 model.
}
\value{
A list containing the following components
\item{coefficients}{the constrained non-negative maximum likelihood estimate of the mean parameters.}
\item{scale}{the maximum likelihood estimate of the scale parameter.}
\item{residuals}{the residuals at the MLE, that is \code{y - fitted.values}}
\item{fitted.values}{the fitted mean values.}
\item{rank}{the number of parameters in the model (named ``\code{rank}" for compatibility ---
we assume that models have full rank)}
\item{family}{included for compatibility --- will always be \code{\link{negbin1}(identity)}.}
\item{linear.predictors}{included for compatibility --- same as \code{fitted.values} (as this is
an identity-link model).}
\item{deviance}{up to a constant, minus twice the maximised log-likelihood (with respect to
a saturated NB1 model with the same \code{scale}).}
\item{aic}{a version of Akaike's \emph{An Information Criterion},
minus twice the maximised log-likelihood plus twice the 
number of parameters.}
\item{aic.c}{a small-sample corrected
version of Akaike's \emph{An Information Criterion}
(Hurvich, Simonoff and Tsai, 1998).}
\item{null.deviance}{the deviance for the null model, comparable with \code{deviance}.
The null model will include the offset and an intercept.}
\item{iter}{the number of iterations of the EM algorithm used.}
\item{weights}{included for compatibility --- a vector of ones.}
\item{prior.weights}{included for compatibility --- a vector of ones.}
\item{standard}{the \code{standard} vector passed to this function.}
\item{df.residual}{the residual degrees of freedom.}
\item{df.null}{the residual degrees of freedom for the null model.}
\item{y}{the \code{y} vector used.}
\item{converged}{logical. Did the ECME algorithm converge 
(according to \code{\link{conv.test}})?}
\item{boundary}{logical. Is the MLE on the boundary of the parameter
space --- i.e. are any of the \code{coefficients < control$bound.tol}?}
\item{loglik}{the maximised log-likelihood.}
\item{nn.design}{the non-negative \code{x} matrix used.}
}
\references{
Hurvich, C. M., J. S. Simonoff and C.-L. Tsai (1998). Smoothing parameter 
selection in non-parametric regression using an improved Akaike 
information criterion.
\emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)} 60(2): 271--293.
}
\author{
Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}.
}
\keyword{ regression }
