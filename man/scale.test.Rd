\name{scaleTest}
\alias{scaleTest}
\title{
Distribution-free Permutation Scale Tests in a Stratified One-Way Layout
}
\description{
Permutation tests for scale parameters in the
distribution-free stratified K-sample one-way layout for binary, ordinal,
numeric and potentially independently right-censored outcomes.
}
\usage{
scaleTest(object, B = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{free1way}.}
 \item{B}{an integer specifying the number of replicates used in the
    permutation test. The default (\code{B = 0}) performs an asymptotic
    permutation test, \code{B > 0} uses \code{B} Monte Carlo replications to
    approximate the permutation distribution.}
  \item{\dots}{additional parameters, currently ignored.}
}
\details{
The distribution function \eqn{F_1} of the response in the control group
(defined by the first factor level of the grouping factor \code{groups}) is
compared to the distribution functions in the remaining groups \eqn{F_k} for
\eqn{k = 2, \dots, K}. No assumptions are made regarding the form of \eqn{F_1},
and thus inference is distribution-free. However, a
semiparametric model assumes that each distribution function \eqn{F_k} only
depends on \eqn{F_1} and a scalar parameter \eqn{\xi_k}. The link
function \eqn{Q} defines the scale model
\deqn{F_k(y) = Q^{-1}(Q(F_1(y)) / \exp(\xi_k))}
for testing \eqn{H_0^k: \xi_k = 0, k = 2, \dots, K}. Under the null, the
model is fitted using \code{free1way} and permutation tests are reported.
}
\value{
An object of class \code{htest}.
}
\keyword{htest}
