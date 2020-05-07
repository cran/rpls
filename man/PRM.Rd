\name{PRM}
\alias{PRM}
\title{Robust PLS}
\description{
Robust PLS by partial robust M-regression.
}
\usage{
PRM(formula,data,a,wfunX,wfunY,center.type,scale.type,usesvd,numit,prec)
}
\arguments{
  \item{formula}{an object of class formula}
  \item{data}{a data frame which contains the variables given in formula}
  \item{a}{number of PLS components}
  \item{wfunX}{weight function to downweight leverage points; predefined weight funcktions "Fair", "Huber", "Tukey" and "Hampel" with respective tuning constants are passed via a list object, e.g. list("Fair",0.95)}
  \item{wfunY}{weight function to downweight residuals; predefined weight funcktions "Fair", "Huber", "Tukey" and "Hampel" with respective tuning constants are passed via a list object, e.g. list("Fair",0.95)}
  \item{center.type}{type of centering of the data in form of a string that matches an R function, e.g. "median"}
  \item{scale.type}{type of scaling for the data in form of a string that matches an R function, e.g. "qn" or alternatively "no" for no scaling}
  \item{numit}{the number of maximal iterations for the convergence of the coefficient estimates}
  \item{prec}{a value for the precision of estimation of the coefficients}
  \item{usesvd}{if TRUE singular value decomposition is performed; logical, default is FALSE}
}
\details{
M regression is used to robustify PLS. Employment of seperate weight functions for leverage points and residuals.
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
}
\references{S. Serneels, C. Croux, P. Filzmoser, and P.J. Van Espen.
  Partial robust M-regression. Chemometrics and Intelligent Laboratory System, Vol. 79(1-2), pp. 55-64, 2005.}
\author{Peter Filzmoser <peter.filzmoser@tuwien.ac.at>}
\examples{
U <- c(rep(2,20), rep(5,30))
X <- replicate(6, U+rnorm(50))
beta <- c(rep(1, 3), rep(-1,3))
e <- c(rnorm(45,0,1.5),rnorm(5,-20,1))
y <- X\%*\%beta + e
d <- as.data.frame(X)
d$y <- y
res <- PRM(y~., data=d, 3, wfunX=list("Fair",0.95),
wfunY=list("Fair",0.95), center.type = "median",
scale.type = "no",usesvd = FALSE,
numit = 100, prec = 0.01)
res$coef
}
\keyword{ multivariate }
