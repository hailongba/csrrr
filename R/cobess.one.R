cobess.one <- function(X, Y, r = 3, kx = 10, ky = 3, eps = 1e-3, maxit = 10, inner.maxit = 10, rho = 1.0, A0 = NULL, B0 = NULL, normalize = TRUE){
  n = dim(X)[1]
  p = dim(X)[2]
  q = dim(Y)[2]
  if(is.null(r)) stop("Please input a rank!")
  r = as.integer(r)
  if(r<1|r>min(p, q)) stop("The input parameter r is illegal!")
  
  if(is.null(A0)) A0 <- svd(Y)$v[,1:r]               
  if(is.null(B0)) B0 <- matrix(0, p, r)
  fit = cobess_one(X, Y, A0, B0, r, kx, ky, eps, maxit, inner.maxit, rho, normalize)
  out = list(C = fit$C, C0 = fit$C0, rank = r, obj = fit$obj, A=fit$A, B=fit$B)
  class(out) = 'cobess.one'
  return(out)
}
