\name{aic}
\alias{aic}
\alias{bic}
\alias{JRRS}
\title{Extract the IC from a "mrbess" or "mrbess.one" object.}
\description{
  These functions are used by \code{mrbess} or \code{mrbess.one} to compute the information criteria from a fitted model object.
}
\usage{
  aic(object,...)
  bic(object,...)
  JRRS(object,...)
}
\arguments{
  \item{object}{Output from the \code{mrbess} function or the \code{mrbess.one} function.}
  \item{\dots}{Additional arguments affecting the predictions produced.}
}
\value{The value of Information Criteria extracted from the "\code{mrbess}" or "\code{mrbess.one}" object.}

\author{
  Canhong Wen, Shijie Quan and Xueqin Wang\cr
  Maintainer: Canhong Wen \email{wench@ustc.edu.cn}
}
\references{
  Wen, C., Zhang, A., Quan, S. and Wang, X. (2017) \emph{BeSS: an R package for best subset selection in linear, logistic and CoxPH models}. \url{https://arxiv.org/abs/1709.06254}.\cr
  Wen, C., Zhu, J. and Wang, X. (2018) \emph{Best Subset Selection in Reduced Rank Regression}, Technical Report.
}

\seealso{
\code{\link{mrbess}}, \code{\link{mrbess.one}}.
}

