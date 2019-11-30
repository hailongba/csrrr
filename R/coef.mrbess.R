coef.mrbess = function(object, type = c("opt", "all"),...){
  type <- match.arg(type)
  if(type == "opt"|is.null(object$coef)){
    C = object$C
  } else {
    C = object$coef[[1]]
  }
  return(C)
}
coef.cobess = function(object, type = c("opt", "all"),...){
  type <- match.arg(type)
  if(type == "opt"|is.null(object$coef)){
    C = object$C
  } else {
    C = object$coef[[1]]
  }
  return(C)
}


coef.mrbess.one = function(object,...){
  C = object$C
  return(C)
}
coef.cobess.one = function(object,...){
  C = object$C
  return(C)
}