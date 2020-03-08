\name{ramml}
\alias{ramml}
\title{Robust Adaptive Modified Maximum Likelihood}
\description{
Modified Maximum Likelihood (MML) estimators are asymptotically equivalent to the ML estimators but their methodology works under the assumption of a known shape parameter. Robust Adaptive MML estimators weaken this assumption and are robust to vertical outliers as well as leverage points.
}
\usage{
ramml(X,y,p,e)
}
\arguments{
  \item{X}{predictor matrix}
  \item{y}{response variable}
  \item{p}{shape parameter of long-tailed symmetric distribution (considered as robustness tuning constant)}
  \item{e}{parameter for the linearization of the intractable term}
}

\value{
  \item{coef}{vector of coefficients}
  \item{scale}{estimate of sigma}
  \item{fitted.values}{vector with fitted y-values}
  \item{residuals}{vector with y-residuals}
}
\references{S. Acitas, Robust Statistical Estimation Methods for High-Dimensional Data with Applications, tech. rep., TUBITAK 2219, International Post Doctoral Research Fellowship Programme, 2019.}
\author{Sukru Acitas <sacitas@eskisehir.edu.tr>}
\keyword{ multivariate }
