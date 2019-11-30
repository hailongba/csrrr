predict.mrbess=function(object, newdata, type = c("opt", "all"),...){
  type <- match.arg(type)
  if(is.null(colnames(newdata))) {
    newx = as.matrix(newdata)
  }else{
    vn = rownames(object$beta)
    if(any(is.na(match(vn, colnames(newdata))))) stop("names of newdata don't match training data!")
    newx = as.matrix(newdata[,vn])
  }
  C = object$C
  C0 = object$C0
  if(type == "opt"|is.null(object$coef[[1]])){
    y = newx%*%C + matrix(C0, nrow(newx), ncol(C), byrow = TRUE)
  } else {
    Coef = object$coef[[1]]
    Coef0 = object$coef[[2]]
    y = array(NA, dim = c(dim(Coef0)[1], dim(Coef0)[2], nrow(newx), ncol(C)))
    for(i in 1:dim(Coef0)[1])
      for(j in 1:dim(Coef0)[2]){
      y[i,j,,] = newx%*%Coef[i,j,,] + matrix(Coef0[i,j,], nrow(newx), ncol(C), byrow = TRUE)
    }
  }
  return(y)
}
predict.cobess=function(object, newdata, type = c("opt", "all"),...){
  type <- match.arg(type)
  if(is.null(colnames(newdata))) {
    newx = as.matrix(newdata)
  }else{
    vn = rownames(object$beta)
    if(any(is.na(match(vn, colnames(newdata))))) stop("names of newdata don't match training data!")
    newx = as.matrix(newdata[,vn])
  }
  C = object$C
  C0 = object$C0
  if(type == "opt"|is.null(object$coef[[1]])){
    y = newx%*%C + matrix(C0, nrow(newx), ncol(C), byrow = TRUE)
  } else {
    Coef = object$coef[[1]]
    Coef0 = object$coef[[2]]
    y = array(NA, dim = c(dim(Coef0)[1], dim(Coef0)[2], nrow(newx), ncol(C)))
    for(i in 1:dim(Coef0)[1])
      for(j in 1:dim(Coef0)[2]){
        y[i,j,,] = newx%*%Coef[i,j,,] + matrix(Coef0[i,j,], nrow(newx), ncol(C), byrow = TRUE)
      }
  }
  return(y)
}

predict.mrbess.one=function(object,newdata, ...){
  if(is.null(colnames(newdata))) {
    newx = as.matrix(newdata)
  }else{
    vn = rownames(object$beta)
    if(any(is.na(match(vn, colnames(newdata))))) stop("names of newdata don't match training data!")
    newx = as.matrix(newdata[,vn])
  }
  C = object$C
  C0 = object$C0
  y = newx%*%C + matrix(C0, nrow(newx), ncol(C), byrow = TRUE)
  return(y)
}
predict.cobess.one=function(object,newdata, ...){
  if(is.null(colnames(newdata))) {
    newx = as.matrix(newdata)
  }else{
    vn = rownames(object$beta)
    if(any(is.na(match(vn, colnames(newdata))))) stop("names of newdata don't match training data!")
    newx = as.matrix(newdata[,vn])
  }
  C = object$C
  C0 = object$C0
  y = newx%*%C + matrix(C0, nrow(newx), ncol(C), byrow = TRUE)
  return(y)
}