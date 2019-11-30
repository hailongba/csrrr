mrbess = function(X, Y, r = NULL, r.max = min(ncol(X), ncol(Y)), k.max = round(min(ncol(X)/2, nrow(X) - 1)),
                  lam = 0, L = NULL, ic.type = c("GIC", "BIC", "AIC","JRRS"),
                  sigma = -1.0, eps = 1e-3, maxit = 100, 
                  inner.maxit = 10, size.max = 1e6, warm.start = TRUE, normalize = TRUE, rho = 1.0){
  n = dim(X)[1]
  p = dim(X)[2]
  q = dim(Y)[2]
  one = rep(1, n)
  if(lam!=0&is.null(L)) stop("Please input a Laplace Matrix!")
  if(is.null(r)){
    r.max = as.integer(r.max)
    if(r.max<1|r.max>min(p, q)) stop("The input parameter r.max is illegal!")
  }else{
    r = as.integer(r)
    if(r<1|r>min(p, q)) stop("The input parameter r is illegal!")
  }   
  k.max = as.integer(k.max)
  if(k.max<1|k.max>=p) stop("The input parameter k.max is illegal!")
  if(k.max>1000) stop("The input parameter k.max should be smaller than 1000!")
  ic.type = match.arg(ic.type)
  if(ic.type == "AIC") type = 'A'
  if(ic.type == "BIC") type = 'B'
  if(ic.type == "GIC") type = 'G'
  if(ic.type == "JRRS") type = 'J'
  if(normalize)
  {
    #center
    meanx = drop(one %*% X)/n
    X = scale(X, meanx, FALSE)
    meany = drop(one %*% Y)/n
    Y = scale(Y, meany, FALSE)
    #normalize
    normx = sqrt(drop(one %*% (X^2)))
    nosignal = normx/sqrt(n) < (.Machine$double.eps)
    if (any(nosignal))  normx[nosignal] = (.Machine$double.eps) * sqrt(n)
    
    names(normx) = NULL
    X = sqrt(n)*scale(X, FALSE, normx)
  }
  if(is.null(r)) {
    A0 = svd(Y)$v[, 1:r.max, drop = FALSE]
    B0 = matrix(0, p, r.max)
    if(lam == 0){
      fit = mrbess_rnotNull(X, Y, A0, B0, matrix(0), matrix(0), matrix(0), r.max, k.max, 0,
                            sigma, eps, maxit, inner.maxit, warm.start, type, size.max, rho)
    }else{
      norm_vec <- function(X) sqrt(sum(X^2))
      X.gen <- t(X)%*%X+lam*L
      if(p>=n) X.gen = X.gen + rho*diag(p)
      X.norm <- diag(X.gen)
      fit = mrbess_rnotNull(X, Y, A0, B0, L, X.gen, X.norm, r.max, k.max, lam, sigma,
                            eps, maxit, inner.maxit, warm.start, type, size.max)
    }
  } else {
    A0 = svd(Y)$v[, 1:r, drop = FALSE]
    B0 = matrix(0, p, r)
    if(lam == 0){
      fit = mrbess_rNull(X, Y, A0, B0, matrix(0), matrix(0), matrix(0), k.max, 0, sigma, r,
                        eps, maxit, inner.maxit, warm.start, type, rho)
    }else
    {
      norm_vec <- function(X) sqrt(sum(X^2))
      X.gen <- t(X)%*%X+lam*L
      if(p>=n) X.gen = X.gen + rho*diag(p)
      X.norm <- diag(X.gen)
      fit = mrbess_rNull(X, Y, A0, B0, L, X.gen, X.norm, k.max, lam, sigma, r,
                        eps, maxit, inner.maxit, warm.start, type)
    }
  }
  if(normalize)
  {
    C = sqrt(n)*fit$C/normx
    C0 = meany-drop(meanx%*%C)
  }else{
    C = fit$C
    C0 = rep(0, q)
  }
  if(!is.null(r)) {
    Coef = array(NA, dim = c(1, k.max, p, q))
    Coef0 = array(NA, dim = c(1, k.max, q))
    if(normalize){
      for(i in 1:k.max){
        Coef[1,i,,] = sqrt(n)*fit$coef[((i-1)*p+1):(i*p),]/normx
        Coef0[1,i,] = meany-drop(meanx%*%Coef[1,i,,])
      }
    }else{
      for(i in 1:k.max){
        Coef[1,i,,] = fit$coef[((i-1)*p+1):(i*p),]
        Coef0[1,i,] = rep(0, q)
      }
    }
  }else if(k.max*r.max*p*q<=size.max){
    Coef = array(NA, dim = c(r.max, k.max, p, q))
    Coef0 = array(NA, dim = c(r.max, k.max, q))
    if(normalize){
      for(i in 1:r.max)
        for(j in 1:k.max){
          kp = k.max*p
          Coef[i,j,,] = sqrt(n)*fit$coefAll[((i-1)*kp+(j-1)*p+1):((i-1)*kp+j*p),]/normx
          Coef0[i,j,] = meany-drop(meanx%*%Coef[i,j,,])
        }
    }else{
      for(i in 1:r.max)
        for(j in 1:k.max){
          kp = k.max*p
          Coef[i,j,,] = fit$coefAll[((i-1)*kp+(j-1)*p+1):((i-1)*kp+j*p),]
          Coef0[i,j,] = rep(0, q)
        }
    }
  }else{
    Coef = NULL
    Coef0 = NULL
  }
  

  out = list(C = C, C0 = C0, coef = list(Coef = Coef, Coef0 = Coef0), rank = fit$rank,
             ms = fit$ms, JRRS = fit$JRRS, AIC = fit$AIC, BIC = fit$BIC, GIC = fit$GIC, ic.type = ic.type)
  class(out) = "mrbess"
  return(out)
}
