// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
using namespace Rcpp;
using namespace std;
void mrbess_one_inner(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, Eigen::MatrixXd& A, Eigen::MatrixXd& B, Eigen::MatrixXd& C, int r=3, int p0=10, double eps=0.001, int maxit=100, int inner_maxit=10, double rho = 1.0){
  int n = X.rows();
  int p = X.cols();
  Eigen::MatrixXd Baset(p0, r);
  Eigen::MatrixXd diset(p-p0, r);
  Eigen::MatrixXd Xaset = Eigen::MatrixXd::Ones(n, p0);
  Eigen::MatrixXd Xiset(n, p-p0);
  vector<int>E(p);
  for(int k=0;k<=p-1;k++) {
    E[k]=k;
  }
  vector<int>iset(p-p0);
  vector<int>aset(p0);
  vector<int>aset0(p0);
  Eigen::VectorXd bd = Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd C0 = B0*A0.adjoint();
  
  double Eps = 0.0;
  for(int j=0;j<maxit;j++){
    for(int k=0;k<p0;k++) {
      aset0[k]=0;
    }
    Eigen::MatrixXd YA = Y*A0;
    Eigen::MatrixXd d = (X.adjoint()*(YA-X*B0))/double(n);
    //cout<<x_norm.replicate(1, r)<<endl;
    B = B0;
    for(int k=0;k<inner_maxit;k++){
      bd = (B+d).rowwise().squaredNorm();
      /*for(int ii=0;ii<p;ii++)
      cout<<bd(ii)<<" ";
      cout<<endl;*/
      for(int mm=0;mm<=p0-1;mm++) {
        bd.maxCoeff(&aset[mm]);
        bd(aset[mm])=0.0;
        //cout<<aset[mm]<<" ";
      }
      //cout<<endl;
      sort (aset.begin(),aset.end());
      set_difference(E.begin(),E.end(), aset.begin(),aset.end(),iset.begin());
      
      B = Eigen::MatrixXd::Zero(p, r);
      d = Eigen::MatrixXd::Zero(p, r);
      
      for(int mm=0;mm<p0;mm++) {
        Xaset.col(mm)=X.col(aset[mm]);
      }
      for(int mm=0;mm<p-p0;mm++) {
        Xiset.col(mm)=X.col(iset[mm]);
      }
      if(p0<n)
        Baset = Xaset.colPivHouseholderQr().solve(YA);
      else{
        Eigen::MatrixXd X_inv = Eigen::MatrixXd::Identity(p0, p0)*rho + Xaset.adjoint()*Xaset;
        X_inv = X_inv.colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(p0, p0));
        Baset = X_inv*X.adjoint()*YA;
      }
      for(int mm=0;mm<p0;mm++) B.row(aset[mm]) = Baset.row(mm);
      diset = Xiset.adjoint()*(YA-Xaset.rightCols(p0)*Baset.bottomRows(p0))/double(n);
      
      for(int mm=0;mm<p-p0;mm++) d.row(iset[mm]) = diset.row(mm);
      if(aset0==aset) break; else  aset0 = aset;
    }
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Y.adjoint()*X*B, Eigen::ComputeFullV | Eigen::ComputeFullU);
    A = svd.matrixU().leftCols(r)*svd.matrixV().adjoint();
    C = B*A.adjoint();
    
    Eps = (Y-X*C).cwiseAbs().colwise().sum().maxCoeff() - (Y-X*C0).cwiseAbs().colwise().sum().maxCoeff();
    //cout<<(Y-X*C).cwiseAbs().colwise().sum().maxCoeff()<<endl;
    //cout<<(Y-X*C0).cwiseAbs().colwise().sum().maxCoeff()<<endl;
    //cout<<Eps<<endl;
    if(abs(Eps)<eps) break; else{
      C0 = C;
      B0 = B;
      A0 = A;
    }
  }
}

// [[Rcpp::export]]
List mrbess_one(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, int r=3, int p0=10, double eps=0.001, int maxit=100, int inner_maxit=10, double rho=1.0){
  int p = X.cols();
  int q = Y.cols();
  Eigen::MatrixXd A(q, r);
  Eigen::MatrixXd B(p ,r);
  Eigen::MatrixXd C(p, q);
  mrbess_one_inner(X, Y, A0, B0, A, B, C, r, p0, eps, maxit, inner_maxit, rho);
  return List::create(Named("A")=A, Named("B")=B, Named("C")=C);
}
void nmrbess_one_inner(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, Eigen::MatrixXd& A, Eigen::MatrixXd& B, Eigen::MatrixXd& C, Eigen::MatrixXd& X_gen, Eigen::VectorXd& x_norm, int r, double lam, Eigen::MatrixXd& L, int p0, double eps=0.001, int maxit=100, int inner_maxit=10){
  int n = X.rows();
  int p = X.cols();
  Eigen::MatrixXd Baset(p0, r);
  Eigen::MatrixXd diset(p-p0, r);
  Eigen::MatrixXd Xaset(n, p0);
  Eigen::MatrixXd Xgenaset(p0, p0);
  vector<int>E(p);
  for(int k=0;k<=p-1;k++) {
    E[k]=k;
  }
  vector<int>iset(p-p0);
  vector<int>aset(p0);
  vector<int>aset0(p0);
  Eigen::VectorXd bd = Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd C0 = B0*A0.adjoint();
  //Eigen::MatrixXd X_gen = X.adjoint()*X + lam*L;
  //Eigen::VectorXd x_norm = X_gen.colwise().squaredNorm().cwiseSqrt();
  double Eps = 0.0;
  for(int j=0;j<maxit;j++){
    for(int k=0;k<p0;k++) {
      aset0[k]=0;
    }
    Eigen::MatrixXd YA = Y*A0;
    Eigen::MatrixXd d = (X.adjoint()*YA-X_gen*B0).cwiseQuotient(x_norm.replicate(1, r));
    //cout<<x_norm.replicate(1, r)<<endl;
    B = B0;
    for(int k=0;k<inner_maxit;k++){
      bd = (B+d).rowwise().squaredNorm();
      /*for(int ii=0;ii<p;ii++)
      cout<<bd(ii)<<" ";
      cout<<endl;*/
      for(int mm=0;mm<=p0-1;mm++) {
        bd.maxCoeff(&aset[mm]);
        bd(aset[mm])=0.0;
        //cout<<aset[mm]<<" ";
      }
      //cout<<endl;
      sort (aset.begin(),aset.end());
      set_difference(E.begin(),E.end(), aset.begin(),aset.end(),iset.begin());
      
      B = Eigen::MatrixXd::Zero(p, r);
      d = Eigen::MatrixXd::Zero(p, r);
      
      for(int mm=0;mm<p0;mm++) {
        Xaset.col(mm)=X.col(aset[mm]);
      }
      for(int mm1=0;mm1<p0;mm1++)
        for(int mm2=0;mm2<p0;mm2++){
          Xgenaset(mm1, mm2)=X_gen(aset[mm1], aset[mm2]);
        }
      Baset = Xgenaset.colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(p0,p0))*Xaset.adjoint()*YA;
      for(int mm=0;mm<p0;mm++) B.row(aset[mm]) = Baset.row(mm);
      
      d = (X.adjoint()*YA-X_gen*B).cwiseQuotient(x_norm.replicate(1, r));
      for(int mm=0;mm<p0;mm++) d.row(aset[mm]) = Eigen::VectorXd::Zero(r);
      if(aset0==aset) break; else  aset0 = aset;
    }
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Y.adjoint()*X*B, Eigen::ComputeFullV | Eigen::ComputeFullU);
    A = svd.matrixU().leftCols(r)*svd.matrixV().adjoint();
    C = B*A.adjoint();
    
    Eps = (Y-X*C).cwiseAbs().colwise().sum().maxCoeff() - (Y-X*C0).cwiseAbs().colwise().sum().maxCoeff();
    //cout<<(Y-X*C).cwiseAbs().colwise().sum().maxCoeff()<<endl;
    //cout<<(Y-X*C0).cwiseAbs().colwise().sum().maxCoeff()<<endl;
    //cout<<Eps<<endl;
    if(abs(Eps)<eps) break; else{
      C0 = C;
      B0 = B;
      A0 = A;
    }
  }
}
// [[Rcpp::export]]
List nmrbess_one(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, Eigen::MatrixXd& X_gen, Eigen::VectorXd& x_norm, int r, double lam, Eigen::MatrixXd& L, int p0, double eps=0.001, int maxit=100, int inner_maxit=10){
  int p = X.cols();
  int q = Y.cols();
  Eigen::MatrixXd A(q, r);
  Eigen::MatrixXd B(p ,r);
  Eigen::MatrixXd C(p, q);
  nmrbess_one_inner(X, Y, A0, B0, A, B, C, X_gen, x_norm, r, lam, L, p0, eps, maxit, inner_maxit);
  return List::create(Named("A")=A, Named("B")=B, Named("C")=C);
}
// [[Rcpp::export]]
List mrbess_rNull(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, Eigen::MatrixXd& L, Eigen::MatrixXd& X_gen, Eigen::VectorXd& x_norm, int p0_max, double lam = 0.0, double sigma = -1.0, int r = 3, double eps = 1e-3, int maxit = 100, int inner_maxit = 10, bool warm_start = true, char opt = 'A', double rho=1.0){
  int n = X.rows();
  int p = X.cols();
  int q = Y.cols();
  double J;
  double pen;
  Eigen::MatrixXd coefC(p0_max*p, q);
  Eigen::VectorXd SSE = Eigen::VectorXd::Zero(p0_max);
  //Eigen::VectorXd MSE = Eigen::VectorXd::Zero(p0_max);
  Eigen::VectorXd AIC = Eigen::VectorXd::Zero(p0_max);
  Eigen::VectorXd BIC = Eigen::VectorXd::Zero(p0_max);
  Eigen::VectorXd GIC = Eigen::VectorXd::Zero(p0_max);
  Eigen::VectorXd JRRS = Eigen::VectorXd::Zero(p0_max);
  Eigen::MatrixXd A = A0;
  Eigen::MatrixXd B = B0;
  Eigen::MatrixXd A_out(q, r);
  Eigen::MatrixXd B_out(p, r);
  Eigen::MatrixXd C_out(p, q);
  Eigen::MatrixXd C0(p, q);
  const double constant = 4.0;
  double tmp = abs(lam-0.0);
  if(sigma == -1.0){
    if(tmp<1e-16)
      mrbess_one_inner(X, Y, A0, B0, A_out, B_out, C_out, r, p0_max, eps, maxit, inner_maxit, rho);
    else
      nmrbess_one_inner(X, Y, A0, B0, A_out, B_out, C_out, X_gen, x_norm, r, lam, L, p0_max, eps, maxit, inner_maxit);
    Eigen::MatrixXd C_full = C_out;
    double S = (Y-X*C_full).squaredNorm()/double(n*q-(q+p0_max-r)*r);
    sigma = sqrt(S);
  }
  for(int i=0;i<p0_max;i++){
    if(warm_start){
      if(tmp<1e-16)
        mrbess_one_inner(X, Y, A, B, A_out, B_out, C_out, r, i+1, eps, maxit, inner_maxit, rho);
      else
        nmrbess_one_inner(X, Y, A, B, A_out, B_out, C_out, X_gen, x_norm, r, lam, L, i+1, eps, maxit, inner_maxit);
    }else{
      if(tmp<1e-16){
        mrbess_one_inner(X, Y, A0, B0, A_out, B_out, C_out, r, i+1, eps, maxit, inner_maxit, rho);
      }else
        nmrbess_one_inner(X, Y, A0, B0, A_out, B_out, C_out, X_gen, x_norm, r, lam, L, i+1, eps, maxit, inner_maxit);
    }
    A = A_out;
    B = B_out;
    C0 = C_out;
    coefC.middleRows(i*p, p) = C0;
    Eigen::VectorXd Csum = C0.rowwise().sum();
    Csum = (Csum.array()!=0).select(1, Csum);
    J = Csum.sum();
    pen = constant*sigma*double(r)*(2*double(q)+(log(2.0)+1.0)*J+J*(1.0+log(double(p)/J)));
    SSE(i) = (X*C0 - Y).squaredNorm();
    //MSE(i) = SSE(i)/double(n);
    AIC(i) = double(n)*log(SSE(i)/double(n)) + 2*(i+1);
    BIC(i) = double(n)*log(SSE(i)/double(n)) +log(n)*(i+1);
    GIC(i) = double(n)*log(SSE(i)/double(n)) + 0.1*double(r)*(log(double(p))*log(log(double(n)))*(i+1) + double(n)/log(double(n)));
    JRRS(i) = SSE(i) + pen;
  }
  int mark = 0;
  if(opt == 'A') AIC.minCoeff(&mark);
  if(opt == 'B') BIC.minCoeff(&mark);
  if(opt == 'G') GIC.minCoeff(&mark);
  if(opt == 'J') JRRS.minCoeff(&mark);
  
  C_out = coefC.middleRows(mark*p, p);
  return List::create(Named("C")=C_out, Named("JRRS")=JRRS, Named("coef")=coefC, Named("rank")=r, Named("SSE")=SSE, Named("AIC")=AIC, Named("BIC")=BIC, Named("GIC")=GIC, Named("ms")=mark+1);
}
// [[Rcpp::export]]
List mrbess_rnotNull(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd& A0, Eigen::MatrixXd& B0, Eigen::MatrixXd& L, Eigen::MatrixXd& X_gen, Eigen::VectorXd& x_norm, int r_max, int p0_max, double lam = 0.0, double sigma = -1.0, double eps = 1e-3, int maxit = 100, int inner_maxit = 10, bool warm_start = true, char opt = 'A', int size_max = 1e6, double rho=1.0){
  int n = X.rows();
  int p = X.cols();
  int q = Y.cols();
  int size = p0_max*r_max*p*q;
  Eigen::MatrixXd coefAll;
  if(size<=size_max) coefAll = Eigen::MatrixXd::Zero(p0_max*r_max*p, q);
  Eigen::VectorXd SSE(p0_max);
  Eigen::VectorXd J(p0_max);
  Eigen::VectorXd pen(p0_max);
  Eigen::MatrixXd C_out(p, q);
  Eigen::MatrixXd C(p, q);
  Eigen::MatrixXd AIC(r_max, p0_max);
  Eigen::MatrixXd BIC(r_max, p0_max);
  Eigen::MatrixXd GIC(r_max, p0_max);
  int rowmin, colmin;
  const double constant = 4.0;
  double tmp = abs(lam-0.0);
  Eigen::MatrixXd JRRS(r_max, p0_max);
  if(sigma == -1.0){
    Eigen::MatrixXd A_out(q, r_max);
    Eigen::MatrixXd B_out(p, r_max);
    if(tmp<1e-16)
      mrbess_one_inner(X, Y, A0, B0, A_out, B_out, C_out, r_max, p0_max, eps, maxit, inner_maxit, rho);
    else
      nmrbess_one_inner(X, Y, A0, B0, A_out, B_out, C_out, X_gen, x_norm, r_max, lam, L, p0_max, eps, maxit, inner_maxit);
    Eigen::MatrixXd C_full = C_out;
    double S = (Y-X*C_full).squaredNorm()/double(n*q-(q+p0_max-r_max)*r_max);
    sigma = sqrt(S);
  }
  for(int i=0;i<r_max;i++){
    Eigen::MatrixXd A = A0.leftCols(i+1);
    Eigen::MatrixXd B = B0.leftCols(i+1);
    List fit = mrbess_rNull(X, Y, A, B, L, X_gen, x_norm, p0_max, lam, sigma, i+1, eps, maxit, inner_maxit, warm_start, opt, rho);
    Eigen::VectorXd SSE0 = fit[4];
    Eigen::VectorXd AIC0 = fit[5];
    Eigen::VectorXd BIC0 = fit[6];
    Eigen::VectorXd GIC0 = fit[7];
    SSE = SSE0;
    AIC.row(i) = AIC0;
    BIC.row(i) = BIC0;
    GIC.row(i) = GIC0;
    Eigen::MatrixXd coefC = fit[2];
    if(size<=size_max) coefAll.middleRows(i*p0_max*p, p0_max*p) = coefC;
    for(int k=0;k<p0_max;k++){
      Eigen::MatrixXd C0 = coefC.middleRows(k*p, p);
      Eigen::VectorXd Csum = C0.rowwise().sum();
      Csum = (Csum.array()!=0).select(1, Csum);
      J(k) = Csum.sum();
      pen(k) = constant*sigma*double(i+1)*(2*double(q)+(log(2.0)+1.0)*J(k)+J(k)*(1.0+log(double(p)/J(k))));
    }
    JRRS.row(i) = SSE + pen;
  }
  if(opt == 'A') AIC.minCoeff(&rowmin, &colmin);
  if(opt == 'B') BIC.minCoeff(&rowmin, &colmin);
  if(opt == 'G') GIC.minCoeff(&rowmin, &colmin);
  if(opt == 'J') JRRS.minCoeff(&rowmin, &colmin);
  Eigen::MatrixXd A = A0.rightCols(rowmin+1);
  Eigen::MatrixXd B = B0.rightCols(rowmin+1);
  Eigen::MatrixXd A_out(q, rowmin+1);
  Eigen::MatrixXd B_out(p, rowmin+1);
  if(tmp<1e-16)
    mrbess_one_inner(X, Y, A, B, A_out, B_out, C_out, rowmin+1, colmin+1, eps, maxit, inner_maxit, rho);
  else
    nmrbess_one_inner(X, Y, A, B, A_out, B_out, C_out, X_gen, x_norm, rowmin+1, lam, L, colmin+1, eps, maxit, inner_maxit);
  return List::create(Named("C")=C_out, Named("JRRS")=JRRS, Named("coefAll")=coefAll, Named("rank")=rowmin+1, Named("ms")=colmin+1, Named("sigma")=sigma, Named("AIC")=AIC, Named("BIC")=BIC, Named("GIC")=GIC);
}
