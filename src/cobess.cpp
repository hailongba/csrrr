// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
using namespace Rcpp;
using namespace std;
void Normalize(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::VectorXd& meanx, Eigen::VectorXd& meany, Eigen::VectorXd& normx){
  int n = X.rows();
  int p = X.cols();
  int q = Y.cols();
  meanx = X.colwise().mean();
  meany = Y.colwise().mean();
  for(int i=0;i<p;i++){
    X.col(i) = X.col(i).array() - meanx(i);
  }
  for(int i=0;i<q;i++){
    Y.col(i) = Y.col(i).array() - meany(i);
  }
  for(int i=0;i<p;i++)
  normx = X.colwise().norm();
  for(int i=0;i<p;i++){
    X.col(i) = sqrt(double(n))*X.col(i)/normx(i);
  }
}
// [[Rcpp::export]]
Eigen::VectorXd cobess_one_inner(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd& A0, Eigen::MatrixXd& B0, Eigen::MatrixXd& A, Eigen::MatrixXd& B, Eigen::MatrixXd& C, int r=3, int kx=10, int ky=3, double eps=0.001, int maxit=10, int inner_maxit=10, double rho = 1.0){
  int n = X.rows();
  int p = X.cols();
  int q = Y.cols();
  Eigen::MatrixXd Baset(kx, r);
  Eigen::MatrixXd Xaset(n, kx);
  Eigen::MatrixXd Ybset(ky, n);
  Eigen::MatrixXd XB(n, r);
  Eigen::MatrixXd bb(r, r);
  Eigen::MatrixXd bb_inv(r, r);
  Eigen::MatrixXd A1(kx, r);
  
  vector<int>aset(kx);
  vector<int>aset0(kx);
  vector<int>bset(ky);
  vector<int>bset0(ky);
  Eigen::VectorXd bd(p);
  Eigen::VectorXd bdy(q);
  Eigen::MatrixXd C0 = B0*A0.adjoint();
  Eigen::MatrixXd d(p, r);
  Eigen::MatrixXd dy(q, r);
  Eigen::MatrixXd pm(n, n);
  Eigen::VectorXd ypm(q);
  Eigen::VectorXd obj(maxit+1);
  
  obj(0) = (Y - X*C0).norm();
  int j = 0;
  for(j=0;j<maxit;j++){
    for(int k=0;k<kx;k++) {
      aset0[k]=0;
    }
    for(int k=0;k<ky;k++) {
      bset0[k]=0;
    }
    Eigen::MatrixXd YA = Y*A0;
    d = (X.adjoint()*(YA-X*B0))/double(n);
    //cout<<x_norm.replicate(1, r)<<endl;
    B = B0;
    for(int k=0;k<inner_maxit;k++){
      bd = (B+d).rowwise().norm();
      /*for(int ii=0;ii<p;ii++)
        cout<<bd(ii)<<" ";
      cout<<endl;*/

      for(int mm=0;mm<=kx-1;mm++) {
        bd.maxCoeff(&aset[mm]);
        bd(aset[mm])=-1.0;
       // cout<<aset[mm]<<" ";
      }
      //cout<<endl;
      B = Eigen::MatrixXd::Zero(p, r);
      
      for(int mm=0;mm<kx;mm++) {
        Xaset.col(mm)=X.col(aset[mm]);
      }
      /*if(kx<n)
        Baset = Xaset.colPivHouseholderQr().solve(YA);
      else{*/
        Eigen::MatrixXd X_inv = Eigen::MatrixXd::Identity(kx, kx)*rho + Xaset.adjoint()*Xaset;
        X_inv = X_inv.colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(kx, kx));
        Baset = X_inv*X.adjoint()*YA;
      //}
      d = X.adjoint()*(YA-Xaset*Baset)/double(n);
      for(int mm=0;mm<kx;mm++){
		    B.row(aset[mm]) = Baset.row(mm);
	    	d.row(aset[mm]) = Eigen::VectorXd::Zero(r);
      }
      if(aset0 == aset) break; else  aset0 = aset;
    }
    
    XB = X*B;
    bb = XB.adjoint()*XB;
    bb_inv = bb.colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(r, r));
    if(bb.squaredNorm()*bb_inv.squaredNorm() > 10e12){
      bb_inv = Eigen::MatrixXd::Zero(r, r);
      for(int mm=0;mm<r;mm++){
        bb_inv(mm, mm) = 1.0/bb(mm, mm);
      }
    }
    
    pm = XB*bb_inv*XB.adjoint();
    for(int k=0;k<q;k++){
      ypm(k) = Y.col(k).adjoint()*pm*Y.col(k);
    }
    for(int mm=0;mm<ky;mm++) {
      ypm.maxCoeff(&bset[mm]);
      ypm(bset[mm])=-1.0;
      //cout<<bset[mm]<<" "; 
    }
    for(int mm=0;mm<ky;mm++) {
      Ybset.row(mm)=Y.col(bset[mm]);
    }
    A = Eigen::MatrixXd::Zero(q, r);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Ybset*X*B, Eigen::ComputeFullV | Eigen::ComputeFullU);   
    if(ky>=r)
      A1 = (svd.matrixU().leftCols(r))*(svd.matrixV().adjoint());
    else
      A1 = (svd.matrixU())*(svd.matrixV().leftCols(ky).adjoint());
    for(int mm=0;mm<ky;mm++) {
      A.row(bset[mm]) = A1.row(mm);
    }
  
    C = B*A.adjoint();
    obj(j+1) = (Y-X*C).norm();
    if(abs(obj(j+1)-obj(j))<eps) break;
    else{
      C0 = C;
      B0 = B;
      A0 = A;
    }
  }
  if(j == maxit) j = j-1;
  return obj.head(j+2);
}
// [[Rcpp::export]]
List cobess_one(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, int r=3, int kx=10, int ky=3, double eps=0.001, int maxit=10, int inner_maxit=10, double rho=1.0, bool normalize=true){
  int n = X.rows();
  int p = X.cols();
  int q = Y.cols();
  Eigen::VectorXd meanx(p);
  Eigen::VectorXd normx(p);
  Eigen::VectorXd meany(q);
  Eigen::VectorXd C0(q);
  if(normalize){
    Normalize(X, Y, meanx, meany, normx);
  }
  Eigen::MatrixXd A(q, r);
  Eigen::MatrixXd B(p ,r);
  Eigen::MatrixXd C(p, q);
  Eigen::VectorXd obj = cobess_one_inner(X, Y, A0, B0, A, B, C, r, kx, ky, eps, maxit, inner_maxit, rho);
  if(normalize){
    for(int i=0;i<q;i++){
      Eigen::VectorXd tmp = C.col(i);
      C.col(i) = sqrt(double(n))*tmp.array()/normx.array();
      C0(i) = meany(i) - C.col(i).dot(meanx);
    }
  }
  return List::create(Named("A")=A, Named("B")=B, Named("C")=C, Named("C0")=C0, Named("rank")=r, Named("obj")=obj);
}
// [[Rcpp::export]]
Eigen::MatrixXd cobess_r_inner(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd& A, Eigen::MatrixXd& C, Eigen::MatrixXd& ICmat, int& rowmin, int& colmin, int kx_max, int ky_max, int r = 3, double eps = 1e-3, int maxit = 10, int inner_maxit = 10, bool warm_start = true, char opt = 'A', double rho = 1.0){
  int n = X.rows();
  int p = X.cols();
  int q = Y.cols();
  rowmin = 0; 
  colmin = 0;
  int IC;
  if(opt == 'A') IC = 0;
  if(opt == 'B') IC = 1;
  if(opt == 'C') IC = 2;
  double df;
  double SSE;
  double tmp = 1.0e20;
  Eigen::MatrixXd coefC(kx_max*p, ky_max*q);  //output
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Y, Eigen::ComputeFullV | Eigen::ComputeFullU);
  Eigen::MatrixXd A00 = svd.matrixV().leftCols(r);
  Eigen::MatrixXd A0 = A00;
  Eigen::MatrixXd A_opt = A00;
  Eigen::MatrixXd B00 = Eigen::MatrixXd::Zero(p, r);
  Eigen::MatrixXd B0 = B00;
  Eigen::MatrixXd B(p, r);
  for(int i=0;i<kx_max;i++)
    for(int j=0;j<ky_max;j++){
      //cout<<"i="<<i<<" j="<<j<<endl;
      if(warm_start){
        cobess_one_inner(X, Y, A0, B0, A, B, C, r, i+1, j+1, eps, maxit, inner_maxit, rho);
        A0 = A;
        B0 = B;
      }
      else
        cobess_one_inner(X, Y, A00, B00, A, B, C, r, i+1, j+1, eps, maxit, inner_maxit, rho);
      
      coefC.block(i*p, j*q, p, q) = C;
      df = double(i*r*min(n,p))/double(p) + j*r - r*r;
      SSE = (Y - X*C).squaredNorm();
      ICmat(0, i*ky_max+j) = log(SSE/double(n*q)) + 2.0*df/double(n*q);
      ICmat(1, i*ky_max+j) = log(SSE/double(n*q)) + df*log(double(n*q))/double(n*q);
      ICmat(2, i*ky_max+j) = log(SSE/double(n*q)) + r*(df*log(log(double(n*q)))*log(double(p*q)) + double(n*q)/(double(n*q)));
      ICmat(3, i*ky_max+j) = SSE;
     
      if(ICmat(IC, i*ky_max+j)<tmp){
        A_opt = A;
        tmp = ICmat(IC, i*ky_max+j);
        rowmin = i;
        colmin = j;
      } 
     // cout<<"rowmin="<<rowmin<<" colmin="<<colmin<<endl;
    }
    C = coefC.block(rowmin*p, colmin*q, p, q);
    A = A_opt;
    return coefC;
}
// [[Rcpp::export]]
List cobess_r(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int kx_max, int ky_max, int r=3, double eps=0.001, int maxit=10, int inner_maxit=10, bool warm_start = true, char opt = 'A', double rho=1.0, bool normalize=true){
  int n = X.rows();
  int p = X.cols();
  int q = Y.cols();
  int rowmin; 
  int colmin;
  Eigen::VectorXd meanx(p);
  Eigen::VectorXd normx(p);
  Eigen::VectorXd meany(q);
  Eigen::VectorXd tmp(p);
  Eigen::MatrixXd coefC(kx_max*p, ky_max*q);
  Eigen::MatrixXd coef0 = Eigen::MatrixXd::Zero(kx_max, ky_max*q);
  Eigen::MatrixXd ICmat(4, kx_max*ky_max);
  if(normalize){
    Normalize(X, Y, meanx, meany, normx);
  }
  Eigen::MatrixXd A(q, r);
  Eigen::MatrixXd B(p ,r);
  Eigen::MatrixXd C(p, q);
  Eigen::VectorXd C0 = Eigen::VectorXd::Zero(q);
  Eigen::MatrixXd Ctmp(p, q);
  coefC = cobess_r_inner(X, Y, A, C, ICmat, rowmin, colmin, kx_max, ky_max, r, eps, maxit, inner_maxit, warm_start, opt, rho);
  if(normalize){
    for(int i=0;i<kx_max;i++)
      for(int j=0;j<ky_max;j++){
        Ctmp = coefC.block(i*p, j*q, p, q);
        for(int k=0;k<q;k++){
          tmp = Ctmp.col(k);
          Ctmp.col(k) = sqrt(double(n))*tmp.array()/normx.array();
          coef0(i,j*q+k) = meany(k) - Ctmp.col(k).dot(meanx);
        }
        coefC.block(i*p, j*q, p, q) = Ctmp;
        //cout<<1<<endl;
      }
    for(int i=0;i<q;i++){
      tmp = C.col(i);
      C.col(i) = sqrt(double(n))*tmp.array()/normx.array();
      C0(i) = meany(i) - C.col(i).dot(meanx);
    }
  }
  return List::create(Named("A")=A, Named("C")=C, Named("C0")=C0, Named("coefC")=coefC, Named("coef0")=coef0, Named("rank")=r, Named("IC")=ICmat);
}
// [[Rcpp::export]]
List cobessC(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int r_max, int kx_max, int ky_max, double sigma = -1.0, double eps=0.001, int maxit=10, int inner_maxit=10, bool warm_start = true, char opt = 'A', double rho=1.0, bool normalize=true){
  int n = X.rows();
  int p = X.cols();
  int q = Y.cols();
  int rowmin;
  int colmin;
  double constant = 4.0;
  Eigen::MatrixXd IC(r_max, kx_max*ky_max);
  Eigen::MatrixXd ICmat(4, kx_max*ky_max);
  Eigen::MatrixXd penmat(r_max, kx_max*ky_max);
  Eigen::VectorXd meanx(p);
  Eigen::VectorXd normx(p);
  Eigen::VectorXd meany(q);
  Eigen::MatrixXd coefC(kx_max*p, ky_max*q);
  Eigen::VectorXd MSE(kx_max*ky_max);
  Eigen::VectorXd JX(kx_max*ky_max);
  Eigen::VectorXd JY(kx_max*ky_max);
  Eigen::VectorXd pen(kx_max*ky_max);
  Eigen::VectorXd tmp(kx_max*ky_max);
  if(normalize){
    Normalize(X, Y, meanx, meany, normx);
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Y, Eigen::ComputeFullV | Eigen::ComputeFullU);
  if(sigma == -1.0){
    Eigen::MatrixXd A0 = svd.matrixV().leftCols(r_max);
    Eigen::MatrixXd B0 = Eigen::MatrixXd::Zero(p, r_max);
    Eigen::MatrixXd A(q, r_max);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(p, r_max);
    Eigen::MatrixXd C(p, q);
    cobess_one_inner(X, Y, A0, B0, A, B, C, r_max, kx_max, ky_max, eps, maxit, inner_maxit, rho);   
    double S = (Y-X*C).squaredNorm()/(double(n*q-(ky_max+kx_max-r_max)*r_max));
    sigma = sqrt(S);
  }
  Eigen::MatrixXd C(p, q);
  Eigen::VectorXd C0(q);
  Eigen::MatrixXd Copt(p, q);
  double mintmp = 1.0e20;
  int ropt = 0;
  int rowopt = 0;
  int colopt = 0;
  for(int r=0;r<r_max;r++){
    Eigen::MatrixXd A(q, r);
    coefC = cobess_r_inner(X, Y, A, C, ICmat, rowmin, colmin, kx_max, ky_max, r+1, eps, maxit, inner_maxit, warm_start, opt, rho);
    MSE = ICmat.row(3)/double(n*q);
    for(int i=0;i<kx_max;i++)
      for(int j=0;j<ky_max;j++){
        JX(i*ky_max+j) = i+1;
        JY(i*ky_max+j) = j+1;
      }
    //cout<<MSE<<"\n"<<endl;
    tmp = 1.0 + log(double(p)) - JX.array().log();
    pen = constant*sigma*(r+1)*(JX+JY) + JX.cwiseProduct(tmp);
    tmp = 1.0 + log(double(q)) - JY.array().log();
    pen = pen + JY.cwiseProduct(tmp);
    penmat.row(r) = pen;
    IC.row(r) = pen + MSE;
    if((IC.row(r)).minCoeff()<mintmp){
      ropt = r;
      Copt = C;
      rowopt = rowmin;
      colopt = colmin;
      mintmp = (IC.row(r)).minCoeff();
    }
  }
  cout<<ropt+1<<" "<<rowopt+1<<" "<<colopt+1<<endl;
  if(normalize){
    for(int i=0;i<q;i++){
      Eigen::VectorXd tmp2 = Copt.col(i);
      Copt.col(i) = sqrt(double(n))*tmp2.array()/normx.array();
      C0(i) = meany(i) - Copt.col(i).dot(meanx);
    }
  }
  return List::create(Named("C")=Copt, Named("C0")=C0, Named("IC")=IC, Named("nrank")=ropt+1, Named("J")=rowopt+1, Named("K")=colopt+1);
}
