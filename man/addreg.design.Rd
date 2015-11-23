\name{addreg.design}
\alias{addreg.design}
\title{
Construct Design Matrix for addreg Model
}
\description{
\code{addreg.design} constructs the design matrix for an
\code{\link{addreg}} model, given a particular parameterisation.

This is a workhorse function --- it would not normally be called directly.
}
\usage{
addreg.design(terms, data, allref, design.ref)
}
\arguments{
  \item{terms}{
	\code{terms} component of object returned
  from a call to \code{\link{addreg.allref}} for the
  desired model.
}
  \item{data}{
	\code{data} component of object returned from
  a call to \code{\link{addreg.allref}} for the desired
  model.
}
  \item{allref}{
	\code{allref} component of object returned
  from a call to \code{\link{addreg.allref}} for the
  desired model.
}
  \item{design.ref}{
	vector of indices for a particular
  parameterisation in \code{allref}. That is, each element
  corresponds to a term \code{x} in the model, and the
  value of the element indicates which item in the list
  \code{allref[[x]]} is the reference level in this
  parameterisation.
}
}
\details{
In the CEM algorithm
employed by \code{\link{addreg}}, we must consider the
Cartesian product of all possible parameterisations. The
list of these for each term in the model is constructed by
a call to \code{\link{addreg.allref}}, and a list of all
possible combinations created using
\code{\link{expand.grid}}.

For a particular combination of reference levels,
\code{addreg.design} constructs the associated design
matrix by transforming \code{data}. Continuous covariates
are transformed such that their minimum or maximum observed
value corresponds to a transformed value of zero;
categorical covariates are transformed by using either
\code{\link{contr.treatment}} with a specified reference
level or \code{\link{contr.isotonic}} so that the levels
are increasing in the specified order.
}
\value{
A strictly non-negative design matrix to be passed to the
relevant fitting function.
}
\author{
Mark W. Donoghoe \email{Mark.Donoghoe@mq.edu.au}
}
\seealso{
\code{\link{addreg.allref}}, \code{\link{model.matrix}},
\code{\link{contr.treatment}}, \code{\link{contr.isotonic}}
}

\keyword{ models }
