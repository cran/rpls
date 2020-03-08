pramml <- function (X, y, a, reg = "lts", pmml, opt = "l1m", usesvd = FALSE)
{


  Unisimpls <- function(X, y, a) {
    n <- nrow(X)
    px <- ncol(X)
    if (px > n) {
      dimensions <- 1
      dimension <- px - n
      ressvd <- svd(t(X))
      X <- ressvd$v %*% diag(ressvd$d)
      n <- nrow(X)
      px <- ncol(X)
    }
    else {
      dimensions <- 0
    }
    s <- t(X) %*% y
    U <- matrix(0, nrow = n, ncol = a)
    V <- matrix(0, nrow = px, ncol = a)
    P <- matrix(0, nrow = px, ncol = a)
    R <- matrix(0, nrow = px, ncol = a)
    B <- V
    for (j in 1:a) {
      r <- s
      u <- X %*% r
      u <- u - U[, 1:max(1, j - 1)] %*% (t(U[, 1:max(1,
                                                     j - 1)]) %*% u)
      normu <- drop(sqrt(t(u) %*% u))
      u <- u/normu
      r <- r/normu
      p <- t(X) %*% u
      v <- p - V[, 1:max(1, j - 1)] %*% (t(V[, 1:max(1,
                                                     j - 1)]) %*% p)
      v <- v/drop(sqrt(t(v) %*% v))
      s <- s - v %*% (t(v) %*% s)
      U[, j] <- u
      P[, j] <- p
      R[, j] <- r
      V[, j] <- v
      B[, j] <- R[, 1:j] %*% t(R[, 1:j]) %*% t(X) %*% y
    }
    if (dimensions == 1) {
      B <- ressvd$u %*% B
    }
    list(coefficients = B, scores = U, loadings = P)
  }
  n = nrow(X)
  p = ncol(X)
  if (usesvd == TRUE) {
    if (p > n) {
      dimensions <- 1
      dimension <- p - n
      ressvd <- svd(t(X))
      X <- ressvd$v %*% diag(ressvd$d)
      n = nrow(X)
      p = ncol(X)
    }
    else {
      dimensions <- 0
    }
  }
  if (opt == "l1m") {
    mx <- l1median(X)
  }
  else if (opt == "mean") {
    mx <- apply(X, 2, mean)
  }
  else {
    mx <- apply(X, 2, median)
  }
  if (opt == "mean") {
    my <- mean(y)
  }
  else {
    my <- median(y)
  }
  Xmc <- as.matrix(scale(X, center = mx, scale = FALSE))
  ymc <- as.vector(scale(y, center = my, scale = FALSE))
  wx <- sqrt(apply(Xmc^2, 1, sum))

  wx <- wx/median(wx)
  wx <- 1/((1+((1/(2*pmml-3))*(wx^2)))^2)
  wy <- abs(ymc)
  wy <- wy/median(wy)
  wy <- 1/((1+((1/(2*pmml-3))*(wy^2)))^2)

  w <- wx * wy
  Xw <- Xmc * sqrt(w)
  yw <- ymc * sqrt(w)

  spls1 <- Unisimpls(Xw, yw, a)
  b <- spls1$coef[, a]
  T <- spls1$sco

  if (reg=="lts") {
    res.reg1 <- ltsReg(T, yw, intercept = FALSE)
  }
  else if(reg=="s") {
    res.reg1 <- lmrob.S(T, yw,lmrob.control(nResample = 50,k.max = 800))
  }

  res.ramml1 <- ramml(T,yw,pmml,res.reg1$residuals/res.reg1$scale[1])
  res.ramml2 <- ramml(T,yw,pmml,res.ramml1$residuals/res.ramml1$scale)

  T <- spls1$sco/sqrt(w)
  r <- (ymc-T%*%res.ramml2$coef)/res.ramml2$scale
  wy <- 1/((1+((1/(2*pmml-3))*(r^2)))^2)
  if (opt == "l1m") {
    mt <- l1median(T)
  }
  else {
    mt <- apply(T, 2, median)
  }
  dt <- T - mt
  wt <- sqrt(apply(dt^2, 1, sum))
  wt <- wt/median(wt)
  wt <- 1/((1+((1/(2*pmml-3))*(wt^2)))^2)
  w <- drop(wy * wt)
  w0 <- which(w == 0)
  if (length(w0) != 0) {
    w <- replace(w, list = w0, values = 10^(-6))
  }
  Xw <- Xmc * sqrt(w)
  yw <- ymc * sqrt(w)
  spls <- Unisimpls(Xw, yw, a)
  b <- spls$coef[, a]
  T <- spls$sco/sqrt(w)

  if (usesvd == TRUE) {
    if (dimensions == 1) {
      b <- drop(ressvd$u %*% b)
      if (opt == "mean") {
        b0 <- mean(y - as.matrix(X) %*% t(ressvd$u) %*%
                     b)
      }
      else {
        b0 <- median(y - as.matrix(X) %*% t(ressvd$u) %*%
                       b)
      }
      yfit <- as.matrix(X) %*% t(ressvd$u) %*% b + b0
    }
    else {
      if (opt == "mean") {
        b0 <- mean(y - as.matrix(X) %*% b)
      }
      else {
        b0 <- median(y - as.matrix(X) %*% b)
      }
      yfit <- as.matrix(X) %*% b + b0
    }
  }
  else {
    if (opt == "mean") {
      b0 <- mean(y - as.matrix(X) %*% b)
    }
    else {
      b0 <- median(y - as.matrix(X) %*% b)
    }
    yfit <- as.matrix(X) %*% b + b0
  }
  list(coef = b, intercept = b0, wy = wy, wt = wt, w = w, scores = T,
       loadings = spls$loadings, fitted.values = yfit, mx = mx,
       my = my)
}
