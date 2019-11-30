// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/csrrr.h"
#include <RcppEigen.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// cobess_one_inner
Eigen::VectorXd cobess_one_inner(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd& A0, Eigen::MatrixXd& B0, Eigen::MatrixXd& A, Eigen::MatrixXd& B, Eigen::MatrixXd& C, int r, int kx, int ky, double eps, int maxit, int inner_maxit, double rho);
static SEXP _csrrr_cobess_one_inner_try(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP rSEXP, SEXP kxSEXP, SEXP kySEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type B(BSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type C(CSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type kx(kxSEXP);
    Rcpp::traits::input_parameter< int >::type ky(kySEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type inner_maxit(inner_maxitSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(cobess_one_inner(X, Y, A0, B0, A, B, C, r, kx, ky, eps, maxit, inner_maxit, rho));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _csrrr_cobess_one_inner(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP rSEXP, SEXP kxSEXP, SEXP kySEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP rhoSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_csrrr_cobess_one_inner_try(XSEXP, YSEXP, A0SEXP, B0SEXP, ASEXP, BSEXP, CSEXP, rSEXP, kxSEXP, kySEXP, epsSEXP, maxitSEXP, inner_maxitSEXP, rhoSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// cobess_one
List cobess_one(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, int r, int kx, int ky, double eps, int maxit, int inner_maxit, double rho, bool normalize);
static SEXP _csrrr_cobess_one_try(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP rSEXP, SEXP kxSEXP, SEXP kySEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP rhoSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type kx(kxSEXP);
    Rcpp::traits::input_parameter< int >::type ky(kySEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type inner_maxit(inner_maxitSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(cobess_one(X, Y, A0, B0, r, kx, ky, eps, maxit, inner_maxit, rho, normalize));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _csrrr_cobess_one(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP rSEXP, SEXP kxSEXP, SEXP kySEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP rhoSEXP, SEXP normalizeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_csrrr_cobess_one_try(XSEXP, YSEXP, A0SEXP, B0SEXP, rSEXP, kxSEXP, kySEXP, epsSEXP, maxitSEXP, inner_maxitSEXP, rhoSEXP, normalizeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// cobess_r_inner
Eigen::MatrixXd cobess_r_inner(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd& A, Eigen::MatrixXd& C, Eigen::MatrixXd& ICmat, int& rowmin, int& colmin, int kx_max, int ky_max, int r, double eps, int maxit, int inner_maxit, bool warm_start, char opt, double rho);
static SEXP _csrrr_cobess_r_inner_try(SEXP XSEXP, SEXP YSEXP, SEXP ASEXP, SEXP CSEXP, SEXP ICmatSEXP, SEXP rowminSEXP, SEXP colminSEXP, SEXP kx_maxSEXP, SEXP ky_maxSEXP, SEXP rSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP warm_startSEXP, SEXP optSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type C(CSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type ICmat(ICmatSEXP);
    Rcpp::traits::input_parameter< int& >::type rowmin(rowminSEXP);
    Rcpp::traits::input_parameter< int& >::type colmin(colminSEXP);
    Rcpp::traits::input_parameter< int >::type kx_max(kx_maxSEXP);
    Rcpp::traits::input_parameter< int >::type ky_max(ky_maxSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type inner_maxit(inner_maxitSEXP);
    Rcpp::traits::input_parameter< bool >::type warm_start(warm_startSEXP);
    Rcpp::traits::input_parameter< char >::type opt(optSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(cobess_r_inner(X, Y, A, C, ICmat, rowmin, colmin, kx_max, ky_max, r, eps, maxit, inner_maxit, warm_start, opt, rho));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _csrrr_cobess_r_inner(SEXP XSEXP, SEXP YSEXP, SEXP ASEXP, SEXP CSEXP, SEXP ICmatSEXP, SEXP rowminSEXP, SEXP colminSEXP, SEXP kx_maxSEXP, SEXP ky_maxSEXP, SEXP rSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP warm_startSEXP, SEXP optSEXP, SEXP rhoSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_csrrr_cobess_r_inner_try(XSEXP, YSEXP, ASEXP, CSEXP, ICmatSEXP, rowminSEXP, colminSEXP, kx_maxSEXP, ky_maxSEXP, rSEXP, epsSEXP, maxitSEXP, inner_maxitSEXP, warm_startSEXP, optSEXP, rhoSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// cobess_r
List cobess_r(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int kx_max, int ky_max, int r, double eps, int maxit, int inner_maxit, bool warm_start, char opt, double rho, bool normalize);
static SEXP _csrrr_cobess_r_try(SEXP XSEXP, SEXP YSEXP, SEXP kx_maxSEXP, SEXP ky_maxSEXP, SEXP rSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP warm_startSEXP, SEXP optSEXP, SEXP rhoSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type kx_max(kx_maxSEXP);
    Rcpp::traits::input_parameter< int >::type ky_max(ky_maxSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type inner_maxit(inner_maxitSEXP);
    Rcpp::traits::input_parameter< bool >::type warm_start(warm_startSEXP);
    Rcpp::traits::input_parameter< char >::type opt(optSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(cobess_r(X, Y, kx_max, ky_max, r, eps, maxit, inner_maxit, warm_start, opt, rho, normalize));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _csrrr_cobess_r(SEXP XSEXP, SEXP YSEXP, SEXP kx_maxSEXP, SEXP ky_maxSEXP, SEXP rSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP warm_startSEXP, SEXP optSEXP, SEXP rhoSEXP, SEXP normalizeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_csrrr_cobess_r_try(XSEXP, YSEXP, kx_maxSEXP, ky_maxSEXP, rSEXP, epsSEXP, maxitSEXP, inner_maxitSEXP, warm_startSEXP, optSEXP, rhoSEXP, normalizeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// cobessC
List cobessC(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int r_max, int kx_max, int ky_max, double sigma, double eps, int maxit, int inner_maxit, bool warm_start, char opt, double rho, bool normalize);
static SEXP _csrrr_cobessC_try(SEXP XSEXP, SEXP YSEXP, SEXP r_maxSEXP, SEXP kx_maxSEXP, SEXP ky_maxSEXP, SEXP sigmaSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP warm_startSEXP, SEXP optSEXP, SEXP rhoSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type r_max(r_maxSEXP);
    Rcpp::traits::input_parameter< int >::type kx_max(kx_maxSEXP);
    Rcpp::traits::input_parameter< int >::type ky_max(ky_maxSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type inner_maxit(inner_maxitSEXP);
    Rcpp::traits::input_parameter< bool >::type warm_start(warm_startSEXP);
    Rcpp::traits::input_parameter< char >::type opt(optSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(cobessC(X, Y, r_max, kx_max, ky_max, sigma, eps, maxit, inner_maxit, warm_start, opt, rho, normalize));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _csrrr_cobessC(SEXP XSEXP, SEXP YSEXP, SEXP r_maxSEXP, SEXP kx_maxSEXP, SEXP ky_maxSEXP, SEXP sigmaSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP warm_startSEXP, SEXP optSEXP, SEXP rhoSEXP, SEXP normalizeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_csrrr_cobessC_try(XSEXP, YSEXP, r_maxSEXP, kx_maxSEXP, ky_maxSEXP, sigmaSEXP, epsSEXP, maxitSEXP, inner_maxitSEXP, warm_startSEXP, optSEXP, rhoSEXP, normalizeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// mrbess_one
List mrbess_one(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, int r, int p0, double eps, int maxit, int inner_maxit, double rho);
static SEXP _csrrr_mrbess_one_try(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP rSEXP, SEXP p0SEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type inner_maxit(inner_maxitSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(mrbess_one(X, Y, A0, B0, r, p0, eps, maxit, inner_maxit, rho));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _csrrr_mrbess_one(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP rSEXP, SEXP p0SEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP rhoSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_csrrr_mrbess_one_try(XSEXP, YSEXP, A0SEXP, B0SEXP, rSEXP, p0SEXP, epsSEXP, maxitSEXP, inner_maxitSEXP, rhoSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// nmrbess_one
List nmrbess_one(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, Eigen::MatrixXd& X_gen, Eigen::VectorXd& x_norm, int r, double lam, Eigen::MatrixXd& L, int p0, double eps, int maxit, int inner_maxit);
static SEXP _csrrr_nmrbess_one_try(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP X_genSEXP, SEXP x_normSEXP, SEXP rSEXP, SEXP lamSEXP, SEXP LSEXP, SEXP p0SEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X_gen(X_genSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type x_norm(x_normSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type inner_maxit(inner_maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(nmrbess_one(X, Y, A0, B0, X_gen, x_norm, r, lam, L, p0, eps, maxit, inner_maxit));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _csrrr_nmrbess_one(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP X_genSEXP, SEXP x_normSEXP, SEXP rSEXP, SEXP lamSEXP, SEXP LSEXP, SEXP p0SEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_csrrr_nmrbess_one_try(XSEXP, YSEXP, A0SEXP, B0SEXP, X_genSEXP, x_normSEXP, rSEXP, lamSEXP, LSEXP, p0SEXP, epsSEXP, maxitSEXP, inner_maxitSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// mrbess_rNull
List mrbess_rNull(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, Eigen::MatrixXd& L, Eigen::MatrixXd& X_gen, Eigen::VectorXd& x_norm, int p0_max, double lam, double sigma, int r, double eps, int maxit, int inner_maxit, bool warm_start, char opt, double rho);
static SEXP _csrrr_mrbess_rNull_try(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP LSEXP, SEXP X_genSEXP, SEXP x_normSEXP, SEXP p0_maxSEXP, SEXP lamSEXP, SEXP sigmaSEXP, SEXP rSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP warm_startSEXP, SEXP optSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type L(LSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X_gen(X_genSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type x_norm(x_normSEXP);
    Rcpp::traits::input_parameter< int >::type p0_max(p0_maxSEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type inner_maxit(inner_maxitSEXP);
    Rcpp::traits::input_parameter< bool >::type warm_start(warm_startSEXP);
    Rcpp::traits::input_parameter< char >::type opt(optSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(mrbess_rNull(X, Y, A0, B0, L, X_gen, x_norm, p0_max, lam, sigma, r, eps, maxit, inner_maxit, warm_start, opt, rho));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _csrrr_mrbess_rNull(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP LSEXP, SEXP X_genSEXP, SEXP x_normSEXP, SEXP p0_maxSEXP, SEXP lamSEXP, SEXP sigmaSEXP, SEXP rSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP warm_startSEXP, SEXP optSEXP, SEXP rhoSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_csrrr_mrbess_rNull_try(XSEXP, YSEXP, A0SEXP, B0SEXP, LSEXP, X_genSEXP, x_normSEXP, p0_maxSEXP, lamSEXP, sigmaSEXP, rSEXP, epsSEXP, maxitSEXP, inner_maxitSEXP, warm_startSEXP, optSEXP, rhoSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// mrbess_rnotNull
List mrbess_rnotNull(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd& A0, Eigen::MatrixXd& B0, Eigen::MatrixXd& L, Eigen::MatrixXd& X_gen, Eigen::VectorXd& x_norm, int r_max, int p0_max, double lam, double sigma, double eps, int maxit, int inner_maxit, bool warm_start, char opt, int size_max, double rho);
static SEXP _csrrr_mrbess_rnotNull_try(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP LSEXP, SEXP X_genSEXP, SEXP x_normSEXP, SEXP r_maxSEXP, SEXP p0_maxSEXP, SEXP lamSEXP, SEXP sigmaSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP warm_startSEXP, SEXP optSEXP, SEXP size_maxSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type L(LSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type X_gen(X_genSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type x_norm(x_normSEXP);
    Rcpp::traits::input_parameter< int >::type r_max(r_maxSEXP);
    Rcpp::traits::input_parameter< int >::type p0_max(p0_maxSEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type inner_maxit(inner_maxitSEXP);
    Rcpp::traits::input_parameter< bool >::type warm_start(warm_startSEXP);
    Rcpp::traits::input_parameter< char >::type opt(optSEXP);
    Rcpp::traits::input_parameter< int >::type size_max(size_maxSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(mrbess_rnotNull(X, Y, A0, B0, L, X_gen, x_norm, r_max, p0_max, lam, sigma, eps, maxit, inner_maxit, warm_start, opt, size_max, rho));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _csrrr_mrbess_rnotNull(SEXP XSEXP, SEXP YSEXP, SEXP A0SEXP, SEXP B0SEXP, SEXP LSEXP, SEXP X_genSEXP, SEXP x_normSEXP, SEXP r_maxSEXP, SEXP p0_maxSEXP, SEXP lamSEXP, SEXP sigmaSEXP, SEXP epsSEXP, SEXP maxitSEXP, SEXP inner_maxitSEXP, SEXP warm_startSEXP, SEXP optSEXP, SEXP size_maxSEXP, SEXP rhoSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_csrrr_mrbess_rnotNull_try(XSEXP, YSEXP, A0SEXP, B0SEXP, LSEXP, X_genSEXP, x_normSEXP, r_maxSEXP, p0_maxSEXP, lamSEXP, sigmaSEXP, epsSEXP, maxitSEXP, inner_maxitSEXP, warm_startSEXP, optSEXP, size_maxSEXP, rhoSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _csrrr_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Eigen::VectorXd(*cobess_one_inner)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,int,int,int,double,int,int,double)");
        signatures.insert("List(*cobess_one)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd,Eigen::MatrixXd,int,int,int,double,int,int,double,bool)");
        signatures.insert("Eigen::MatrixXd(*cobess_r_inner)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,int&,int&,int,int,int,double,int,int,bool,char,double)");
        signatures.insert("List(*cobess_r)(Eigen::MatrixXd&,Eigen::MatrixXd&,int,int,int,double,int,int,bool,char,double,bool)");
        signatures.insert("List(*cobessC)(Eigen::MatrixXd&,Eigen::MatrixXd&,int,int,int,double,double,int,int,bool,char,double,bool)");
        signatures.insert("List(*mrbess_one)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd,Eigen::MatrixXd,int,int,double,int,int,double)");
        signatures.insert("List(*nmrbess_one)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXd&,Eigen::VectorXd&,int,double,Eigen::MatrixXd&,int,double,int,int)");
        signatures.insert("List(*mrbess_rNull)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::VectorXd&,int,double,double,int,double,int,int,bool,char,double)");
        signatures.insert("List(*mrbess_rnotNull)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::VectorXd&,int,int,double,double,double,int,int,bool,char,int,double)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _csrrr_RcppExport_registerCCallable() { 
    R_RegisterCCallable("csrrr", "_csrrr_cobess_one_inner", (DL_FUNC)_csrrr_cobess_one_inner_try);
    R_RegisterCCallable("csrrr", "_csrrr_cobess_one", (DL_FUNC)_csrrr_cobess_one_try);
    R_RegisterCCallable("csrrr", "_csrrr_cobess_r_inner", (DL_FUNC)_csrrr_cobess_r_inner_try);
    R_RegisterCCallable("csrrr", "_csrrr_cobess_r", (DL_FUNC)_csrrr_cobess_r_try);
    R_RegisterCCallable("csrrr", "_csrrr_cobessC", (DL_FUNC)_csrrr_cobessC_try);
    R_RegisterCCallable("csrrr", "_csrrr_mrbess_one", (DL_FUNC)_csrrr_mrbess_one_try);
    R_RegisterCCallable("csrrr", "_csrrr_nmrbess_one", (DL_FUNC)_csrrr_nmrbess_one_try);
    R_RegisterCCallable("csrrr", "_csrrr_mrbess_rNull", (DL_FUNC)_csrrr_mrbess_rNull_try);
    R_RegisterCCallable("csrrr", "_csrrr_mrbess_rnotNull", (DL_FUNC)_csrrr_mrbess_rnotNull_try);
    R_RegisterCCallable("csrrr", "_csrrr_RcppExport_validate", (DL_FUNC)_csrrr_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_csrrr_cobess_one_inner", (DL_FUNC) &_csrrr_cobess_one_inner, 14},
    {"_csrrr_cobess_one", (DL_FUNC) &_csrrr_cobess_one, 12},
    {"_csrrr_cobess_r_inner", (DL_FUNC) &_csrrr_cobess_r_inner, 16},
    {"_csrrr_cobess_r", (DL_FUNC) &_csrrr_cobess_r, 12},
    {"_csrrr_cobessC", (DL_FUNC) &_csrrr_cobessC, 13},
    {"_csrrr_mrbess_one", (DL_FUNC) &_csrrr_mrbess_one, 10},
    {"_csrrr_nmrbess_one", (DL_FUNC) &_csrrr_nmrbess_one, 13},
    {"_csrrr_mrbess_rNull", (DL_FUNC) &_csrrr_mrbess_rNull, 17},
    {"_csrrr_mrbess_rnotNull", (DL_FUNC) &_csrrr_mrbess_rnotNull, 18},
    {"_csrrr_RcppExport_registerCCallable", (DL_FUNC) &_csrrr_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_csrrr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}