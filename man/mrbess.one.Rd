\name{mrbess.one}
\alias{mrbess.one}
\title{
  Multivariate Response Best Subset Selection with a Specified Model Size
}
\description{
  Implement the multivariate response best subset selection method for selecting a given number of relevant predictors in reduced rank regression. The rank in the reduced rank regression is prespecified in argument \code{r}.   
}
\usage{
mrbess.one(X, Y, r = NULL, k = 1, lam = 0, L = NULL, eps = 0.001, maxit = 100,
inner.maxit = 10, A0 = NULL, B0 = NULL, normalize = TRUE, rho = 1.0)
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
  \item{k}{
  Size of the selected model.It controls number of nonzero coefficients to be allowed in the model. Default value is \code{k = 1}.
}
  \item{lam}{
  Regularization parameter of the Laplacian matrix. Default value is \code{lam = 1}.
}
  \item{L}{
  Laplacian matrix of network in \code{x}. The Laplacian matrix of network in \code{x} is defined by \eqn{L = D - A}, where \eqn{D = diag{d_1, d_2,\dots, d_p}} is the degree matrix and \eqn{A} is the adjacency matrix. If there exists a network, then a network based \code{mrbess.one} is applied; for details, see Wen et al.(2018, 2). Default is \code{L = NULL}, which corresponds to no network in \code{x}.
}
  \item{eps}{
  Convergence threshold for the algorithm with given model size. Each loop continues until the absolute change in the objective function after any coefficient matrix update is less than \code{eps}. Defaults value is \code{1e-3}.
}
  \item{maxit}{
  Maximum number of iterations in updating the coefficient matrix for all model sizes; default is \code{100}.
}
  \item{inner.maxit}{
  Maximum number of iterations in inner loop for updating the left matrix of dimension \code{p * r}. The default value is \code{10}.
}
  \item{A0}{
  Initialization value for \eqn{A}. If not specified, matrix \eqn{A} will be initialized as a matrix with all elements being zero. Default is \code{A0 = NULL}.
  }
    \item{B0}{
  Initialization value for \eqn{B}. If not specified, matrix \eqn{B} will be initialized as a matrix whose columns contain the first \code{r} right singular vectors of \code{y}. Default is \code{B0 = NULL}.
  }
  \item{normalize}{
  Logical flag for \code{x} variable normalization and \code{y} variable centralization, prior to fitting the model sequence. The coefficient matrix is always returned on the original scale. Default is \code{normalize = TRUE}. 
}
  \item{rho}{
    
  }
}

\details{
The best subset selection problem with model size \code{k} in reduced rank regression is
\deqn{min ||y - C_0 - x C||_F^2, s.t. rank(C) \leq r, ||C||_{2,0} = k,}
where \eqn{1 \le r \le min(rank(X), q, k)}. If \code{x} is normalized and \code{y} is centralized, then the above problem reduces to 
\deqn{min ||\tilde{y} - \tilde{x} C||_F^2, s.t. rank(C) \leq r, ||C||_{2,0} = k.}

By ranging \code{k} from 1 to \eqn{min(n,p)}, one could derive the solution for the best subset selection problem in reduced rank regression. 

To solve the above problem, we express coefficient matrix \code{C} as a product of two matrices, \code{C=BV^T}, with \code{V } being an orthogonal matrix of dimension \code{q*r} and \code{B} being a matrix of dimension \code{p*r}. Then the best subset problem can be rewritten as
\deqn{min ||\tilde{y} - \tilde{x}BA^T||_F^2, s.t. V^T V = I, ||B||_{2,0} = k.}

For a given rank \code{r} and each candidate model size \code{k}, the best subset selection problem for reduced rank regression is solved by a block-wise algorithm between \eqn{B} and \eqn{A}. The optimization problem involving \eqn{A} reduces to an orthogonal Procrustes problem and has explicit solution. The optimization with respect to \eqn{B} is done with the primal dual active set (PDAS) algorithm, see Wen et al.(2017) and Wen et al.(2018) for details. This algorithm utilizes an active set updating strategy via primal and dual variables and fits the sub-model by exploiting the fact that their support set are non-overlap and complementary. 

}

\value{
An list with class attribute \code{"mrbess"} and named components:
  \item{C }{The estimated coefficient matrix of dimension \code{p * q}.}
  \item{C0 }{Intercept vector of length \code{q}.}
  \item{rank }{The rank \code{r}.}
  \item{SSE}{The value of sum of squared error. }
  \item{AIC}{The value of AIC. }
  \item{BIC}{The value of BIC. }
}
\references{
  Wen, C., Zhang, A., Quan, S. and Wang, X. (2017) \emph{BeSS: an R package for best subset selection in linear, logistic and CoxPH models}. \url{https://arxiv.org/abs/1709.06254}.\cr
  Wen, C., Zhu, J. and Wang, X. (2018) \emph{Best Subset Selection in Reduced Rank Regression}, Technical Report.\cr
  Wen, C., and Wang, X. (2018) \emph{Network-based Best Subset Selection in Reduced Rank Regression}, Technical Report.
  
}
\author{
  Canhong Wen, Shijie Quan and Xueqin Wang\cr
  Maintainer: Canhong Wen \email{wench@ustc.edu.cn}
}
\seealso{
  \code{print}, \code{predict}, \code{coef} and \code{plot} methods, and the \code{mrbess} function. 
}
