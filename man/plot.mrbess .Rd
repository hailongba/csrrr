\name{plot.mrbess}
\alias{plot.mrbess}
\title{plot information criteria from a "mrbess" object}
\description{
Produces an information criteria paths for a fitted "mrbess" object.
}
\usage{
  \method{plot}{mrbess}(object, ic.type = c("JRRS", "BIC", "AIC"), xvar = c("rank", "ms"), r = NULL,...)
}
\arguments{
  \item{object}{Output from the \code{mrbess.one} function.}
  \item{ic.type}{Which information criteria path to plot}
  \item{xvar}{ What is on the X-axis.Either size or model size.}
  \item{r}{ If rank = "ms", this parameter will choose the model with rank "r" to plot.}
  \item{\dots}{Other graphical parameters to plot}
}

\author{
  Canhong Wen, Shijie Quan and Xueqin Wang\cr
  Maintainer: Canhong Wen \email{wench@ustc.edu.cn}
}
\references{
  Wen, C., Zhang, A., Quan, S. and Wang, X. (2017) \emph{BeSS: an R package for best subset selection in linear, logistic and CoxPH models}. \url{https://arxiv.org/abs/1709.06254}.\cr
  Wen, C., Zhu, J. and Wang, X. (2018) \emph{Best Subset Selection in Reduced Rank Regression}, Technical Report.
}
\seealso{
\code{\link{mrbess}}}



