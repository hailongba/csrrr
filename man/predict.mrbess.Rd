\name{predict.mrbess}
\alias{predict.mrbess}
\alias{predict.mrbess.one}
\title{make predictions from a "mrbess" object.}
\description{
  Similar to other predict methods, which returns predictions from a fitted "\code{mrbess}" or "\code{mrbess.one}" object.
}
\usage{
  \method{predict}{mrbess}(object, newdata, type = c("opt", "all"),...)
  \method{predict}{mrbess.one}(object,newdata, ...)
}
\arguments{
  \item{object}{Output from the \code{mrbess} or the \code{mrbess.one} function.}
  \item{newdata}{New data used for prediction.}
  \item{type}{Types of prediction required.  \code{ type = "opt"} corresponds to predictors in best model; \code{ type = "ALL"} corresponds to all predictors in the \code{mrbess} object. Default is \code{ type = "opt"}. }
  \item{\dots}{Additional arguments affecting the predictions produced.}
}
\value{The object returned depends on the \code{type}.}

\author{
  Canhong Wen, Shijie Quan and Xueqin Wang\cr
  Maintainer: Canhong Wen \email{wench@ustc.edu.cn}
}
\references{
  Wen, C., Zhang, A., Quan, S. and Wang, X. (2017) \emph{BeSS: an R package for best subset selection in linear, logistic and CoxPH models}. \url{https://arxiv.org/abs/1709.06254}.\cr
  Wen, C., Zhu, J. and Wang, X. (2018) \emph{Best Subset Selection in Reduced Rank Regression}, Technical Report.
}

\seealso{
\code{\link{mrbess}}, \code{\link{mrbess.one}}
}


