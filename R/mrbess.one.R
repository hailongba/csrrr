mrbess.one = function(X, Y, r = NULL, k = 1, lam = 0, L = NULL, eps = 0.001, maxit = 100, inner.maxit = 10, A0 = NULL, B0 = NULL, normalize = TRUE, rho = 1.0){
  ## initialization
  n = dim(X)[1]
  p = dim(X)[2]
  q = dim(Y)[2]
  if(is.null(r)) stop("Please input a rank!")
  if(lam!=0&is.null(L)) stop("Please input a Laplace Matrix!")
  r = as.integer(r)
  k = as.integer(k)
  if(r<1|r>min(p, q)) stop("The input parameter r is illegal!")
  if(k<1|k>=p) stop("The input parameter k is illegal!")
  if(k>1000) stop("The input parameter k should be smaller than 1000!")
  if(is.null(A0)) A0 = svd(Y)$v[,1:r,drop = FALSE]
  if(is.null(B0)) B0 = matrix(0, p, r)
  one = rep(1, n)
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
  if(is.null(L)){
    fit = mrbess_one(X, Y, A0 = A0, B0 = B0, r = r, p0 = k,
                    eps = eps, maxit = maxit, inner_maxit = inner.maxit, rho = rho)
    C = fit$C
  }else{
    norm_vec <- function(X) sqrt(sum(X^2))
    X.gen <- t(X)%*%X+lam*L
    if(p>=n) X.gen = X.gen + rho*diag(p)
    X.norm <- diag(X.gen)
    fit = nmrbess_one(X, Y, A0 = A0, B0 = B0, X.gen, X.norm, r = r, lam = lam, L = L,
                      p0 = k, eps = eps, maxit = maxit, inner_maxit = inner.maxit)
    C = fit$C
  }
  if(normalize)
  {
    C = sqrt(n)*C/normx
    C0 = meany-drop(meanx%*%C)
  }
  SSE = norm(X%*%C-Y, type="F")^2
  AIC = n*log(SSE/n) + 2*k
  BIC = n*log(SSE/n) + log(n)*k
  out = list(C = C, C0 = C0, rank = r, SSE = SSE, AIC = AIC, BIC = BIC)
  class(out) = 'mrbess.one'
  return(out)
}
