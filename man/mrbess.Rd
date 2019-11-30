\name{mrbess}
\alias{mrbess}
\title{
  Multivariate Response Best subset selection in Reduced Rank Regression
}
\description{
  Implement the multivariate response best subset selection (MrBeSS) method for selecting relevant predictors in reduced rank regression. Provide two alternative ways to tune the parameters of model size and rank: 1) When a rank is prespecified, an optimal model size is determined by BIC; 2) When a rank is not given, an optimal pair of parameters (model size and rank) is determined by the JRRS criterion introduced by Bunea et al. (2012). 
}
\usage{
mrbess(X, Y, r = NULL, r.max = min(ncol(X), ncol(Y)), k.max = round(min(ncol(X)/2, nrow(X) - 1)),
                  lam = 0, L = NULL, ic.type = c("BIC", "AIC", "JRRS"),
                  sigma = -1.0, eps = 1e-3, maxit = 100, 
                  inner.maxit = 10, size.max = 1e6, warm.start = TRUE,
                  normalize = TRUE, rho = 1.0)
}
\arguments{
  \item{X}{
  Covariate matrix, of dimension n * p; each row is an observation vector.
}
  \item{Y}{
  Response matrix, of dimension n * q; each row is an observation vector.
}
  \item{r}{
  Prespecified rank. Default is \code{r = NULL}.
}
  \item{r.max}{
  The maximum rank to be considered, used only when \code{r = NULL} is not given.
}
  \item{k.max}{
  Limit the maximum number of predictors ever to be nonzero; default is \code{min(p/2, n - 1)}.
}
  \item{lam}{
  Regularization parameter of the Laplacian matrix. Default value is \code{lam = 1}.
}
  \item{L}{
  Laplacian matrix of network in \code{x}. The Laplacian matrix of network in \code{x} is defined by \eqn{L = D - A}, where \eqn{D = diag{d_1, d_2,\dots, d_p}} is the degree matrix and \eqn{A} is the adjacency matrix. If there exists a network, then a network based \code{mrbess.one} is applied; for details, see Wen et al.(2018, 2). Default is \code{L = NULL}, which corresponds to no network in \code{x}.
}
 \item{ic.type}{
  The type of information criterion used for tune the parameters.
}
  \item{sigma}{
  The estimated standard deviation for the error term, used for calculating the JRRS criterion introduced by Bunea et al. (2012).
}
  \item{eps}{
  Convergence threshold for the algorithm with given model size. Each loop continues until the absolute change in the objective function after any coefficient matrix update is less than \code{eps}. Defaults value is \code{1e-3}.
}
  \item{maxit}{
  Maximum number of iterations in unpdating the coefficient matrix for all model sizes; default is \code{100}.
}
  \item{inner.maxit}{
  Maximum number of iterations in inner loop for unpdating the left matrix of dimension \code{p * r}. The default value is \code{10}.
}
  \item{size.max}{
  Thresholding value to determine whether to output all the estimated coefficient matrices. If \code{r.max * k.max * p * q > size.max}, then \code{coef} in output stores all coefficient matrices; otherwise, \code{coef} returns a list with null contents. The defaulst value is \code{1e6}.
}
  \item{warm.start}{
  Logical flag for warm start used in solving the multivariate response best subset selection problem with a list of model sizes. If \code{warm.start = TRUE}, then it use the solution from the best subset selection with size \code{k-1} as a initialization value for problem with size \code{k}. Default is \code{warm.start = TRUE}.}
  \item{normalize}{
  Logical flag for \code{x} variable normalization and \code{y} variable centralization, prior to fitting the model sequence. The coefficient matrix is always returned on the original scale. Default is \code{normalize = TRUE}. 
}
  \item{rho}{
  
  }
}
\details{
The best subset selection problem with model size \code{k} in reduced rank regression is
\deqn{min ||y - C_0 - x C||_F^2, s.t. rank(C) \leq r, ||C||_{2,0} = k,}
where \eqn{1 \le r \le min(rank(X), q, k)}. If \code{x} is normalized and \code{y} is centerlized, then the above problem reduces to 
\deqn{min ||\tilde{y} - \tilde{x} C||_F^2, s.t. rank(C) \leq r, ||C||_{2,0} = k.}

By ranging \code{k} from 1 to \eqn{min(n,p)}, one could derive the solution for the best subset selection problem in reduced rank regression. 

To solve the above problem, we express coefficient matrix \code{C} as a product of two matrices, \code{C=BV^T}, with \code{V } being an orthogonal matrix of dimension \code{q*r} and \code{B} being a matrix of dimension \code{p*r}. Then the best subset problem can be rewritten as
\deqn{min ||\tilde{y} - \tilde{x}BA^T||_F^2, s.t. V^T V = I, ||B||_{2,0} = k.}

For a given rank \code{r} and each candidate model size \code{k}, the best subset selection problem for reduced rank regression is solved by a block-wise algorithm between \eqn{B} and \eqn{A}. The optimization problem involving \eqn{A} reduces to an orthogonal Procrustes problem and has explicit solution. The optimization with respect to \eqn{B} is done with the primal dual active set (PDAS) algorithm, see Wen et al.(2017) and Wen et al.(2018) for details. This algorithm utilizes an active set updating strategy via primal and dual variables and fits the sub-model by exploiting the fact that their support set are non-overlap and complementary. 

When the rank \code{r} is pre-specified, an optimal model size \code{k} is determined by BIC. When the rank \code{r} is not specified, \code{r = NULL}, an optimal pair of parameters (model size \code{k} and rank \code{r} ) is determined by the JRRS criterion introduced by Bunea et al. (2012). 
}
\value{
An list with class attribute \code{"mrbess"} and named components:
  \item{C }{The optimal coefficient matrix of dimension \code{p * q}.}
  \item{C0 }{Intercept vector of length \code{q}.}
  \item{coef }{A list contains all coefficient matrices and intercept vectors for all candidate model sizes, i.e., \code{k = 1,\dots, k.max}. The first component is an array of dimension \code{r.max * k.max * p * q} representing the coefficient matrix; the second component is an array of dimension \code{ r.max * k.max * q} representing the intercept vectors.}
  \item{rank }{The optimal rank determined by the JRRS criterion when \code{r = NULL}; the prespecified rank when \code{r} is given.}
  \item{ms }{The optimal model size. Determined by the JRRS criterion for \code{r = NULL} or by BIC for \code{r} being specified.}
    \item{JRRS }{The values of the JRRS criterion. A matrix of dimension \code{r.max * k.max} when \code{r = NULL}; a vector of length \code{k.max} when \code{r} is given. }
  \item{AIC}{The values of AIC. A matrix of dimension \code{r.max * k.max} when \code{r = NULL}; a vector of length \code{k.max} when \code{r} is given. }
  \item{BIC}{The values of BIC. A matrix of dimension \code{r.max * k.max} when \code{r = NULL}; a vector of length \code{k.max} when \code{r} is given. }
   \item{ic.type}{The type of information criterion used for tune the parameters.}
}
\references{
  Wen, C., Zhang, A., Quan, S. and Wang, X. (2017) \emph{BeSS: an R package for best subset selection in linear, logistic and CoxPH models}. \url{https://arxiv.org/abs/1709.06254}.\cr
  Wen, C., Zhu, J. and Wang, X. (2018) \emph{Best Subset Selection in Reduced Rank Regression}, Technical Report.
  
}
\author{
  Canhong Wen, Shijie Quan and Xueqin Wang\cr
  Maintainer: Canhong Wen \email{wench@ustc.edu.cn}
}
\seealso{
  \code{print}, \code{predict}, \code{coef} and \code{plot} methods, and the \code{mrbess.one} function. 
}