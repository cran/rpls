PRM <- function(formula, data, a, wfunX=list("Fair",0.95), wfunY=list("Fair",0.95), center.type = "median", scale.type = "qn", usesvd = FALSE,
                numit = 30, prec = 0.01)
{

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  y <- as.vector(y)
  X <- model.matrix(mt, mf)
  Xattrib <- attributes(X)
  intercept <- which(Xattrib$assign == 0)
  if (length(intercept)) {
    X <- X[, -intercept, drop = FALSE]
  }
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
    R <- matrix(0, nrow = px, ncol = a)
    P <- matrix(0, nrow = px, ncol = a)
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
      R[, j] <- r
      V[, j] <- v
      P[, j] <- p
      B[, j] <- R[, 1:j] %*% t(R[, 1:j]) %*% t(X) %*% y
    }
    if (dimensions == 1) {
      B <- ressvd$u %*% B
    }
    list(coefficients = B, scores = U, loadings = P)
  }
  scale_matrix <- function(Data){
    if (center.type == "l1m"){
      Data.center <- l1median(Data)
    } 
    else {
      Data.center <- apply(Data,2,center.type)
    }
    Data.scaled <- (Data - matrix(Data.center,nrow=dim(Data)[1],
                                  ncol=dim(Data)[2],byrow=TRUE))
    if(!(scale.type=="no")){
      Data.scale <- apply(Data,2,scale.type)
      if (any(1/Data.scale>1e19) | any(is.nan(Data.scale))){
        Data.scale <- apply(Data,2,sd)
        if (any(1/Data.scale>1e19) | any(is.nan(Data.scale))){
          Data.scale <- rep(1,ncol(Data))
          warning("Routine used scale.type='no' to avoide division by zero or infinity.")
        } else {
          warning("Routine used scale.type='sd' to avoide division by zero or infinity.")
        }
      }
      Data.scaled <- Data.scaled/matrix(Data.scale,nrow=dim(Data)[1],ncol=dim(Data)[2],byrow=TRUE)
    } else {
      Data.scale <- rep(1,ncol(Data))
    }
    return(Data.scaled)
  }

  scale_vector<-function(Data){

    if(center.type=="l1m"){
      Data.center <- l1median(Data)
    } else {
      center_fun <- match.fun(center.type)
      Data.center <- center_fun(Data)
    }
    Data.scaled <- Data - Data.center
    if(!(scale.type=="no")){
      scale_fun <- match.fun(scale.type)
      Data.scale <- scale_fun(Data)

      if (any(1/Data.scale>1e19) | any(is.nan(Data.scale))){
        Data.scale <- sd(Data)
        if (any(1/Data.scale>1e19) | any(is.nan(Data.scale))){
          Data.scale <- rep(1,length(Data))
          warning("Routine used scale.type='no' to avoide division by zero or infinity.")
        } else {
          warning("Routine used scale.type='sd' to avoide division by zero or infinity.")
        }
      }
      Data.scaled <- as.vector(scale(Data,center = Data.center,scale = FALSE))
    }
    return(Data.scaled)
  }

  fairW <- function(z,probct){
    probct<-as.double(probct)
    wy <- z
    wy <- 1/((1 + abs(z/(probct)))^2)
    return(wy)
  }

  huberW <- function(z,probct){
    probct<-as.double(probct)
    wx <- z
    r <- z
    wx[which(abs(r) <= probct)] <- 1
    wx[which(abs(r) > probct)] <- probct/abs(wx[which(abs(r) > probct)])
    return(wx)
  }

  tukeyW <- function(z,probct){
    probct<-as.double(probct)
    wx <- z
    r <- z
    wx[which(abs(r) >= probct)] <- 0
    wx[which(abs(r) < probct)] <- (1-(wx[which(abs(r) < probct)]/probct)^2)^2
    return(wx)
  }

  hampelW <- function(z,probct1,hampelb,hampelr){
    probct<-as.double(probct1)
    hampelb<-as.double(hampelb)
    hampelr <- as.double(hampelr)
    r <-z
    wye <- r
    wye[which(abs(r) <= probct)] <- 1
    wye[which(abs(r) > probct & abs(r) <= hampelb)] <- probct/abs(r[which(abs(r) >
                                                                            probct & abs(r) <= hampelb)])
    wye[which(abs(r) > hampelb & abs(r) <= hampelr)] <- probct *
      (hampelr - abs(r[which(abs(r) > hampelb & abs(r) <= hampelr)]))/(hampelr -
                                                                         hampelb) * 1/abs(r[which(abs(r) > hampelb & abs(r) <= hampelr)])
    wye[which(abs(r) > hampelr)] <- 0
    wy <- wye
    return(wy)
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

  Xmc <- as.matrix(scale_matrix(X))
  ymc <- as.matrix(scale_vector(y))

  wx <- sqrt(apply(Xmc^2, 1, sum))
  wx <- wx/median(wx)

  wy <- abs(ymc)
  wy <- wy/median(wy)

  if(is.function(wfunX)){
    z <- wx
    wx <- wfunX(z)
  }
  else{
    z <- wx
    if(is.list(wfunX)){
      if(wfunX[1] == "Fair"){
        wx <- fairW(z,qnorm(as.double(wfunX[2])))
      }
      else if(wfunX[1] == "Huber"){
        wx <- huberW(z,qnorm(as.double(wfunX[2])))
      }
      else if(wfunX[1] == "Tukey"){
        wx <- tukeyW(z,qnorm(as.double(wfunX[2])))
      }
      else if(wfunX[1] == "Hampel"){
        probct1 <-  qnorm(as.double(wfunX[2]))
        hampelb <- qnorm(as.double(wfunX[3]))
        hampelr <- qnorm(as.double(wfunX[4]))
        wx <- hampelW(z,probct1,hampelb,hampelr)
      }else{
        stop("Wrong format: wfunX")
      }
    }
    else{
      stop("Wrong format: wfunX")
    }
  }

  if(is.function(wfunY)){
    z <- wy
    wy <- wfunY(z)
  }
  else{
    z <- wy
    if(is.list(wfunY)){
      if(wfunY[1] == "Fair"){
        wy <- fairW(z,qnorm(as.double(wfunY[2])))
      }
      else if(wfunY[1] == "Huber"){
        wy <- huberW(z,qnorm(as.double(wfunY[2])))
      }
      else if(wfunY[1] == "Tukey"){
        wy <- tukeyW(z,qnorm(as.double(wfunY[2])))
      }
      else if(wfunY[1] == "Hampel"){
        probct1 <-  qnorm(as.double(wfunY[2]))
        hampelb <- qnorm(as.double(wfunY[3]))
        hampelr <- qnorm(as.double(wfunY[4]))
        wy <- hampelW(z,probct1,hampelb,hampelr)
      } else{
        stop("Wrong format: wfunY")
      }
    }
    else{
      stop("Wrong format: wfunY")
    }
  }

  w <- drop(wx * wy)
  w0 <- which(w == 0)
  if (length(w0) != 0) {
    w <- replace(w, list = w0, values = 10^(-6))
  }
  Xw <- Xmc * sqrt(w)
  yw <- ymc * sqrt(w)

  loops <- 1
  ngamma <- 10^5
  difference <- 1

  while ((difference > prec) && loops < numit) {
    ngammaold <- ngamma
    spls <- Unisimpls(Xw, yw, a)
    b <- spls$coef[, a]
    gamma <- t(t(yw) %*% spls$sco)
    T <- spls$sco/sqrt(w)
    r <- ymc - T %*% gamma
    rc <- r - median(r)
    r <- rc/median(abs(rc))
    dt <- scale_matrix(T)
    wt <- sqrt(apply(dt^2, 1, sum))
    wt <- wt/median(wt)

    if(is.function(wfunX)){
      z <- wt
      wt <- wfunX(z)
    }
    else{
      z <- wt
      if(is.list(wfunX)){
        if(wfunX[1] == "Fair"){
          probct1 <-  qchisq(as.double(wfunX[2]),a)
          wt <- fairW(z,probct1)
        }
        else if(wfunX[1] == "Huber"){
          probct1 <-  qchisq(as.double(wfunX[2]),a)
          wt <- huberW(z,probct1)
        }
        else if(wfunX[1] == "Tukey"){
          probct1 <-  qchisq(as.double(wfunX[2]),a)
          wt <- tukeyW(z,probct1)
        }
        else if(wfunX[1] == "Hampel"){
          probct1 <-  qchisq(as.double(wfunX[2]),a)
          hampelb <- qchisq(as.double(wfunX[3]),a)
          hampelr <- qchisq(as.double(wfunX[4]),a)
          wt <- hampelW(z,probct1,hampelb,hampelr)
          wt <-as.numeric(wt)
        }else{
          stop("Wrong format: wfunX")
        }
      }
      else{
        stop("Wrong format: wfunX")
      }
    }
    if(is.function(wfunY)){
      z <- r
      wy <- wfunY(z)
    }
    else{
      z <- r
      if(is.list(wfunY)){
        if(wfunY[1] == "Fair"){
          wy <- fairW(z,qnorm(as.double(wfunY[2])))
        }
        else if(wfunY[1] == "Huber"){
          wy <- huberW(z,qnorm(as.double(wfunY[2])))
        }
        else if(wfunY[1] == "Tukey"){
          wy <- tukeyW(z,qnorm(as.double(wfunY[2])))
        }
        else if(wfunY[1] == "Hampel"){
          probct1 <-  qnorm(as.double(wfunY[2]))
          hampelb <- qnorm(as.double(wfunY[3]))
          hampelr <- qnorm(as.double(wfunY[4]))
          wy <- hampelW(z,probct1,hampelb,hampelr)
          wy <- as.numeric(wy)
        }else{
          stop("Wrong format: wfunY")
        }
      }
      else{
        stop("Wrong format: wfunY")
      }
    }

    ngamma <- sqrt(sum(gamma^2))
    difference <- abs(ngamma - ngammaold)/ngamma

    w <- drop(wy * wt)
    w0 <- which(w == 0)
    if (length(w0) != 0) {
      w <- replace(w, list = w0, values = 10^(-6))
    }

    Xw <- Xmc * sqrt(w)
    yw <- ymc * sqrt(w)

    loops <- loops + 1
  }

  if (usesvd == TRUE) {
    if (dimensions == 1) {
      b <- drop(ressvd$u %*% b)
      if (center.type == "mean") {
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
      if (center.type == "mean") {
        b0 <- mean(y - as.matrix(X) %*% b)
      }
      else {
        b0 <- median(y - as.matrix(X) %*% b)
      }
      yfit <- as.matrix(X) %*% b + b0
    }
  }
  else {
    if (center.type == "mean") {
      b0 <- mean(y - as.matrix(X) %*% b)
    }
    else {
      b0 <- median(y - as.matrix(X) %*% b)
    }
    yfit <- as.matrix(X) %*% b + b0
  }

  list(coef = b, intercept = b0, wy = wy, wt = wt, w = w, scores = T,
       loadings = spls$loadings, fitted.values = yfit)
}
