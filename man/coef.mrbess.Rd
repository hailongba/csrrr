\name{coef.mrbess}
\alias{coef.mrbess}
\alias{coef.mrbess.one}
\title{Provides estimated coefficient matrix from a fitted "mrbess" or "mrbess.one" object.}
\description{
  Similar to other prediction methods, this function provides estimated coefficient matrix from a fitted "\code{mrbess}" or "\code{mrbess.one}" object.
}
\usage{
    \method{coef}{mrbess}(object, type = c("opt", "all"),...)
    \method{coef}{mrbess.one}(object,...)
}
\arguments{
  \item{object}{A "\code{mrbess}" object or a "\code{mrbess.one}" object.}
  \item{type}{Types of coefficients returned.  \code{ type = "opt"} corresponds to the coefficient in the optimal model; \code{ type = "all"} corresponds to all the coefficient matrices. Default is { type = "opt"}. }
  \item{\dots}{Other arguments.}
}
\value{
The object returned depends on \code{object} and \code{type}.
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
\code{\link{mrbess}}, \code{\link{mrbess.one}}
}

\examples{

}
