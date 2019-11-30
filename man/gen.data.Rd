\name{gen.data}
\alias{gen.data}
\title{Generate simulated data}
\description{
   Generate data for simulations.
}
\usage{
gendata(n = 100, p = 30, q = 10, p0 = 10, q0 = 10, nrank = 3, snr = 1, C = NULL,
        sigma = NULL, corrX = c("AR","CS"), rhoX = 0.5, corrE = c("I","AR","CS"), rhoE = 0)
}
\arguments{
  \item{n}{The number of observations.}
  \item{p}{The number of covariates.}
  \item{q}{The number of response variables.}
}
\value{
A list with the following components: x, y, true coefficients.
 \item{X}{Design matrix of predictors.}
 \item{Y}{Response variable}
 \item{C}{The coefficients used in the underlying regression model.}
}
\references{
  Wen, C., Zhang, A., Quan, S. and Wang, X. (2017) \emph{BeSS: an R package for best subset selection in linear, logistic and CoxPH models}. \url{https://arxiv.org/abs/1709.06254}.\cr
  Wen, C., Zhu, J. and Wang, X. (2018) \emph{Best Subset Selection in Reduced Rank Regression}, Technical Report.
}

\author{
  Canhong Wen, Shijie Quan and Xueqin Wang\cr
  Maintainer: Canhong Wen \email{wench@ustc.edu.cn}
}

