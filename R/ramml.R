ramml <- function(X,y,p,e) {

  X <- as.matrix(X)

  n <- length(y)
  m <- dim(X)[2]
  q <- 2*p-3

  if (m == 1) {
    mx <- median(X)
  }
  else {
    mx <-l1median(X)
  }
  Xmc <- as.matrix(scale(X, center = mx, scale = FALSE))
  wx <- sqrt(apply(Xmc^2, 1, sum))
  du <- wx/median(wx)

  t <- e

  betay <- 1/((1+((1/q)*(t^2)))^2)
  betax <- 1/((1+((1/q)*(du^2)))^2)
  beta <- betay*betax
  alfa <- (1/q*t)/((1+((1/q)*(t^2)))^2)
  alfa <- alfa*betax;

  WB <- diag(as.vector(beta), n,n)
  WA <- diag(as.vector(alfa), n,n)

  K <- solve(t(X)%*%WB%*%X)%*%(t(X)%*%WB%*%y)
  D <- solve(t(X)%*%WB%*%X)%*%(t(X)%*%WA%*%c(rep(1,n)))

  B <- (2*p/q)*t(y-X%*%K)%*%WA%*%c(rep(1,n));
  C <- (2*p/q)*t(y-X%*%K)%*%WB%*%(y-X%*%K);

  sigmamml <- (B+sqrt(B^2+4*n*C))/(2*sqrt(n*(n-m)));
  betamml <- K+D%*%sigmamml;
  yfit  <- X%*%betamml
  residual <- y-yfit

  list(coef = betamml, scale=as.numeric(sigmamml), fitted.values = yfit, residuals = residual)
}
