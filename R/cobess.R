cobess = function(X, Y, r = NULL, r.max = min(ncol(X), ncol(Y)), kx.max = round(min(ncol(X)/2, nrow(X) - 1)), ky.max = round(min(ncol(Y)/2, nrow(X) - 1)),
                  ic.type = c("BIC", "AIC", "GIC"),
                  sigma = -1.0, eps = 1e-3, maxit = 10, inner.maxit = 10, warm.start = TRUE,
                  normalize = TRUE, rho = 1.0){
  n = dim(X)[1]
  p = dim(X)[2]
  q = dim(Y)[2]
  if(is.null(r)){
    r.max = as.integer(r.max)
    if(r.max<1|r.max>min(p, q)) stop("The input parameter r.max is illegal!")
  }else{
    r = as.integer(r)
    if(r<1|r>min(p, q)) stop("The input parameter r is illegal!")
  }   
  ic.type = match.arg(ic.type)
  if(ic.type == "AIC") type = 'A'
  if(ic.type == "BIC") type = 'B'
  if(ic.type == "GIC") type = 'G'
  if(is.null(r)){
    fit = cobessC(X, Y, r.max, kx.max, ky.max, sigma, eps, maxit, inner.maxit, warm.start, 
                  type, rho, normalize)
    IC = array(NA, dim = c(r.max, kx.max, ky.max))
    for(r in 1:r.max){
      IC[r,,] = matrix(fit$IC[r,], kx.max, ky.max, byrow = TRUE)
    }
    out = list(C = fit$C, C0 = fit$C0, nrank = fit$nrank, J = fit$J, K = fit$K, IC = IC,
               ic.type = ic.type)
  }else{
    fit = cobess_r(X, Y, kx.max, ky.max, r, eps, maxit, inner.maxit, warm.start, 
                   type, rho, normalize)
    IC = array(NA, dim = c(3, kx.max, ky.max))
    for(r in 1:3){
      IC[r,,] = matrix(fit$IC[r,], kx.max, ky.max, byrow = TRUE)
    }
    rownames(IC) = c("AIC", "BIC", "GIC")
    coefC = array(NA, dim = c(kx.max, ky.max, p, q))
    coef0 = array(NA, dim = c(kx.max, ky.max, q))
    for(i in 1:kx.max)
      for(j in 1:ky.max){
        coefC[i,j,,] = fit$coefC[((i-1)*p+1):(i*p), ((j-1)*q+1):(j*q)]
        coef0[i,j,] = fit$coef0[i, ((j-1)*q+1):(j*q)]
      }
    out = list(C = fit$C, C0 = fit$C0, coef = list(Coef = coefC, Coef0 = coef0),
               rank = fit$rank, IC = IC)
  }
  
  
  class(out) = "cobess"
  return(out)
}
