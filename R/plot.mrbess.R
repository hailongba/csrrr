plot.mrbess = function(object, ic.type = c("JRRS", "BIC", "AIC"), xvar = c("rank", "ms"), r = NULL,...){
  fit = object
  ic.type = match.arg(ic.type)
  xvar = match.arg(xvar)
  val = switch(ic.type,
                JRRS = fit$JRRS,
                BIC = fit$BIC,
                AIC = fit$AIC)
  if(is.vector(val)) {
    plot(val, type="b", xlab = "Model size", ylab = ic.type,
         main = paste0("Plot of ", ic.type, " against model size with specified rank = ", fit$rank))
  }else if(xvar=="ms"){
    if(is.null(r)){
      matplot(t(val), type="b", xlab = "Model size", ylab = ic.type, main = paste0("Plot of ", ic.type, " against model size"))
    }else{
      if(r>nrow(val)) stop("rank is too large!")
      plot(val[r,], type="b", xlab = "Model size", ylab = ic.type,
           main = paste0("Plot of ", ic.type, " against model size with specified rank = ", r))
    }
  }else{
    matplot(val, type="b", xlab = "Rank", ylab = ic.type, main = paste0("Plot of ", ic.type, " against rank"))
  }
}