\name{pramml}
\alias{pramml}
\title{Partial Robust Adaptive Modified Maximum Likelihood}
\description{
Robust Adaptive Modified Maximum Likelihood (RAMML) estimators can be used in the context of PLS to obtain scores and loadings in the latent regression model. The corresponding method is called Partial RAMML (PRAMML).
}
\usage{
pramml(X, y, a, reg = "lts", pmml, opt = "l1m", usesvd = FALSE)
}
\arguments{
  \item{X}{predictor matrix}
  \item{y}{response variable}
  \item{a}{number of PLS components}
  \item{reg}{regression procedure to be used to compute initial estimate of parameter for the linearization of the intractable term; choices are LTS regression ("lts") and S regression ("s") }
  \item{pmml}{shape parameter of long-tailed symmetric distribution (considered as robustness tuning constant)}
  \item{opt}{if "l1m" the mean centering is done by the l1-median; otherwise if "median" the coordinate-wise median is taken}
  \item{usesvd}{if TRUE singular value decomposition is performed; logical, default is FALSE}
}

\value{
  \item{coef}{vector with regression coefficients}
  \item{intercept}{coefficient for intercept}
  \item{wy}{vector of length(y) with residual weights}
  \item{wt}{vector of length(y) with weights for leverage}
  \item{w}{overall weights}
  \item{scores}{matrix with PLS X-scores}
  \item{loadings}{matrix with PLS X-loadings}
  \item{fitted.values}{vector with fitted y-values}
  \item{loadings}{column means of X}
  \item{fitted.values}{mean of y}
}
\references{S. Acitas, Robust Statistical Estimation Methods for High-Dimensional Data with Applications, tech. rep., TUBITAK 2219, International Post Doctoral Research Fellowship Programme, 2019.}
\author{Sukru Acitas <sacitas@eskisehir.edu.tr>}
\examples{
U <- c(rep(2,20), rep(5,30))
X <- replicate(6, U+rnorm(50))
beta <- c(rep(1, 3), rep(-1,3))
e <- c(rnorm(45,0,1.5),rnorm(5,-20,1))
y <- X\%*\%beta + e
res <- pramml(X, y, 4,"s", 16.5, opt ="l1m")
}
\keyword{ multivariate }
