gen.data = function(n = 100, p = 30, q = 10, p0 = 10, q0 = 10, nrank = 3, snr = 1, C = NULL,
                    sigma = NULL, corrX = c("AR","CS"), rhoX = 0.5, corrE = c("I","AR","CS"), rhoE = 0){
  corrX = match.arg(corrX)
  corrE = match.arg(corrE)
  if(corrX=="AR"){
    SigmaX = matrix(0, p, p)
    for(i in 1:p)
      for(j in 1:p)
        SigmaX[i,j] = rhoX^abs(i-j)
  }else{
    SigmaX = matrix(rhoX, p, p) + diag(1-rhoX, p)
  }
  SigmaE = diag(rep(1,q))
  if(corrE=="AR"){
    SigmaE = matrix(0, q, q)
    for(i in 1:p)
      for(j in 1:p)
        SigmaE[i,j] = rhoE^abs(i-j)
  }else{
    SigmaE = matrix(rhoE, q, q) + diag(1-rhoE, q)
  }
  A = B = NULL
  if(is.null(C)){
    A1 = matrix(ncol = nrank, nrow = q0, rnorm(q0 * nrank))
    A0 = matrix(ncol = nrank, nrow = q - q0, 0)
    A = rbind(A1, A0)
    B1 = matrix(ncol = nrank, nrow = p0, rnorm(p0 * nrank))
    B0 = matrix(ncol = nrank, nrow = p - p0, 0)
    B = rbind(B1, B0)
    C = B %*% t(A)
  }

  X = MASS::mvrnorm(n, rep(0, p), SigmaX)
  E = MASS::mvrnorm(n, rep(0, q), SigmaE)

  if(is.null(sigma)){
    sigma = sqrt( sum(diag(t(C)%*%cov(X)%*%C))/ (sum(diag(cov(E)))*snr) )
  }
  Y = X %*% C + sigma*E
  list(Y = Y, X = X, C = C, A = A, B = B,  sigma = sigma,
       SigmaX = SigmaX, SigmaE = SigmaE)

}