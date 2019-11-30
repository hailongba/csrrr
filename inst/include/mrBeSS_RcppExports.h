// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_MrBeSS_RCPPEXPORTS_H_GEN_
#define RCPP_MrBeSS_RCPPEXPORTS_H_GEN_

#include <RcppEigen.h>
#include <Rcpp.h>

namespace MrBeSS {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("MrBeSS", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("MrBeSS", "_MrBeSS_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in MrBeSS");
            }
        }
    }

    inline Eigen::VectorXd cobess_one_inner(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd& A0, Eigen::MatrixXd& B0, Eigen::MatrixXd& A, Eigen::MatrixXd& B, Eigen::MatrixXd& C, int r = 3, int kx = 10, int ky = 3, double eps = 0.001, int maxit = 10, int inner_maxit = 10, double rho = 1.0) {
        typedef SEXP(*Ptr_cobess_one_inner)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cobess_one_inner p_cobess_one_inner = NULL;
        if (p_cobess_one_inner == NULL) {
            validateSignature("Eigen::VectorXd(*cobess_one_inner)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,int,int,int,double,int,int,double)");
            p_cobess_one_inner = (Ptr_cobess_one_inner)R_GetCCallable("MrBeSS", "_MrBeSS_cobess_one_inner");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cobess_one_inner(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(A0)), Shield<SEXP>(Rcpp::wrap(B0)), Shield<SEXP>(Rcpp::wrap(A)), Shield<SEXP>(Rcpp::wrap(B)), Shield<SEXP>(Rcpp::wrap(C)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(kx)), Shield<SEXP>(Rcpp::wrap(ky)), Shield<SEXP>(Rcpp::wrap(eps)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(inner_maxit)), Shield<SEXP>(Rcpp::wrap(rho)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Eigen::VectorXd >(rcpp_result_gen);
    }

    inline List cobess_one(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, int r = 3, int kx = 10, int ky = 3, double eps = 0.001, int maxit = 10, int inner_maxit = 10, double rho = 1.0, bool normalize = true) {
        typedef SEXP(*Ptr_cobess_one)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cobess_one p_cobess_one = NULL;
        if (p_cobess_one == NULL) {
            validateSignature("List(*cobess_one)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd,Eigen::MatrixXd,int,int,int,double,int,int,double,bool)");
            p_cobess_one = (Ptr_cobess_one)R_GetCCallable("MrBeSS", "_MrBeSS_cobess_one");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cobess_one(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(A0)), Shield<SEXP>(Rcpp::wrap(B0)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(kx)), Shield<SEXP>(Rcpp::wrap(ky)), Shield<SEXP>(Rcpp::wrap(eps)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(inner_maxit)), Shield<SEXP>(Rcpp::wrap(rho)), Shield<SEXP>(Rcpp::wrap(normalize)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline Eigen::MatrixXd cobess_r_inner(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd& A, Eigen::MatrixXd& C, Eigen::MatrixXd& ICmat, int& rowmin, int& colmin, int kx_max, int ky_max, int r = 3, double eps = 1e-3, int maxit = 10, int inner_maxit = 10, bool warm_start = true, char opt = 'A', double rho = 1.0) {
        typedef SEXP(*Ptr_cobess_r_inner)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cobess_r_inner p_cobess_r_inner = NULL;
        if (p_cobess_r_inner == NULL) {
            validateSignature("Eigen::MatrixXd(*cobess_r_inner)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,int&,int&,int,int,int,double,int,int,bool,char,double)");
            p_cobess_r_inner = (Ptr_cobess_r_inner)R_GetCCallable("MrBeSS", "_MrBeSS_cobess_r_inner");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cobess_r_inner(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(A)), Shield<SEXP>(Rcpp::wrap(C)), Shield<SEXP>(Rcpp::wrap(ICmat)), Shield<SEXP>(Rcpp::wrap(rowmin)), Shield<SEXP>(Rcpp::wrap(colmin)), Shield<SEXP>(Rcpp::wrap(kx_max)), Shield<SEXP>(Rcpp::wrap(ky_max)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(eps)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(inner_maxit)), Shield<SEXP>(Rcpp::wrap(warm_start)), Shield<SEXP>(Rcpp::wrap(opt)), Shield<SEXP>(Rcpp::wrap(rho)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Eigen::MatrixXd >(rcpp_result_gen);
    }

    inline List cobess_r(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int kx_max, int ky_max, int r = 3, double eps = 0.001, int maxit = 10, int inner_maxit = 10, bool warm_start = true, char opt = 'A', double rho = 1.0, bool normalize = true) {
        typedef SEXP(*Ptr_cobess_r)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cobess_r p_cobess_r = NULL;
        if (p_cobess_r == NULL) {
            validateSignature("List(*cobess_r)(Eigen::MatrixXd&,Eigen::MatrixXd&,int,int,int,double,int,int,bool,char,double,bool)");
            p_cobess_r = (Ptr_cobess_r)R_GetCCallable("MrBeSS", "_MrBeSS_cobess_r");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cobess_r(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(kx_max)), Shield<SEXP>(Rcpp::wrap(ky_max)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(eps)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(inner_maxit)), Shield<SEXP>(Rcpp::wrap(warm_start)), Shield<SEXP>(Rcpp::wrap(opt)), Shield<SEXP>(Rcpp::wrap(rho)), Shield<SEXP>(Rcpp::wrap(normalize)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List cobessC(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, int r_max, int kx_max, int ky_max, double sigma = -1.0, double eps = 0.001, int maxit = 10, int inner_maxit = 10, bool warm_start = true, char opt = 'A', double rho = 1.0, bool normalize = true) {
        typedef SEXP(*Ptr_cobessC)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_cobessC p_cobessC = NULL;
        if (p_cobessC == NULL) {
            validateSignature("List(*cobessC)(Eigen::MatrixXd&,Eigen::MatrixXd&,int,int,int,double,double,int,int,bool,char,double,bool)");
            p_cobessC = (Ptr_cobessC)R_GetCCallable("MrBeSS", "_MrBeSS_cobessC");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_cobessC(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(r_max)), Shield<SEXP>(Rcpp::wrap(kx_max)), Shield<SEXP>(Rcpp::wrap(ky_max)), Shield<SEXP>(Rcpp::wrap(sigma)), Shield<SEXP>(Rcpp::wrap(eps)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(inner_maxit)), Shield<SEXP>(Rcpp::wrap(warm_start)), Shield<SEXP>(Rcpp::wrap(opt)), Shield<SEXP>(Rcpp::wrap(rho)), Shield<SEXP>(Rcpp::wrap(normalize)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List mrbess_one(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, int r = 3, int p0 = 10, double eps = 0.001, int maxit = 100, int inner_maxit = 10, double rho = 1.0) {
        typedef SEXP(*Ptr_mrbess_one)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_mrbess_one p_mrbess_one = NULL;
        if (p_mrbess_one == NULL) {
            validateSignature("List(*mrbess_one)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd,Eigen::MatrixXd,int,int,double,int,int,double)");
            p_mrbess_one = (Ptr_mrbess_one)R_GetCCallable("MrBeSS", "_MrBeSS_mrbess_one");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mrbess_one(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(A0)), Shield<SEXP>(Rcpp::wrap(B0)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(p0)), Shield<SEXP>(Rcpp::wrap(eps)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(inner_maxit)), Shield<SEXP>(Rcpp::wrap(rho)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List nmrbess_one(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, Eigen::MatrixXd& X_gen, Eigen::VectorXd& x_norm, int r, double lam, Eigen::MatrixXd& L, int p0, double eps = 0.001, int maxit = 100, int inner_maxit = 10) {
        typedef SEXP(*Ptr_nmrbess_one)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_nmrbess_one p_nmrbess_one = NULL;
        if (p_nmrbess_one == NULL) {
            validateSignature("List(*nmrbess_one)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXd&,Eigen::VectorXd&,int,double,Eigen::MatrixXd&,int,double,int,int)");
            p_nmrbess_one = (Ptr_nmrbess_one)R_GetCCallable("MrBeSS", "_MrBeSS_nmrbess_one");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_nmrbess_one(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(A0)), Shield<SEXP>(Rcpp::wrap(B0)), Shield<SEXP>(Rcpp::wrap(X_gen)), Shield<SEXP>(Rcpp::wrap(x_norm)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(lam)), Shield<SEXP>(Rcpp::wrap(L)), Shield<SEXP>(Rcpp::wrap(p0)), Shield<SEXP>(Rcpp::wrap(eps)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(inner_maxit)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List mrbess_rNull(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd A0, Eigen::MatrixXd B0, Eigen::MatrixXd& L, Eigen::MatrixXd& X_gen, Eigen::VectorXd& x_norm, int p0_max, double lam = 0.0, double sigma = -1.0, int r = 3, double eps = 1e-3, int maxit = 100, int inner_maxit = 10, bool warm_start = true, char opt = 'A', double rho = 1.0) {
        typedef SEXP(*Ptr_mrbess_rNull)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_mrbess_rNull p_mrbess_rNull = NULL;
        if (p_mrbess_rNull == NULL) {
            validateSignature("List(*mrbess_rNull)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd,Eigen::MatrixXd,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::VectorXd&,int,double,double,int,double,int,int,bool,char,double)");
            p_mrbess_rNull = (Ptr_mrbess_rNull)R_GetCCallable("MrBeSS", "_MrBeSS_mrbess_rNull");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mrbess_rNull(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(A0)), Shield<SEXP>(Rcpp::wrap(B0)), Shield<SEXP>(Rcpp::wrap(L)), Shield<SEXP>(Rcpp::wrap(X_gen)), Shield<SEXP>(Rcpp::wrap(x_norm)), Shield<SEXP>(Rcpp::wrap(p0_max)), Shield<SEXP>(Rcpp::wrap(lam)), Shield<SEXP>(Rcpp::wrap(sigma)), Shield<SEXP>(Rcpp::wrap(r)), Shield<SEXP>(Rcpp::wrap(eps)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(inner_maxit)), Shield<SEXP>(Rcpp::wrap(warm_start)), Shield<SEXP>(Rcpp::wrap(opt)), Shield<SEXP>(Rcpp::wrap(rho)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

    inline List mrbess_rnotNull(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::MatrixXd& A0, Eigen::MatrixXd& B0, Eigen::MatrixXd& L, Eigen::MatrixXd& X_gen, Eigen::VectorXd& x_norm, int r_max, int p0_max, double lam = 0.0, double sigma = -1.0, double eps = 1e-3, int maxit = 100, int inner_maxit = 10, bool warm_start = true, char opt = 'A', int size_max = 1e6, double rho = 1.0) {
        typedef SEXP(*Ptr_mrbess_rnotNull)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_mrbess_rnotNull p_mrbess_rnotNull = NULL;
        if (p_mrbess_rnotNull == NULL) {
            validateSignature("List(*mrbess_rnotNull)(Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::MatrixXd&,Eigen::VectorXd&,int,int,double,double,double,int,int,bool,char,int,double)");
            p_mrbess_rnotNull = (Ptr_mrbess_rnotNull)R_GetCCallable("MrBeSS", "_MrBeSS_mrbess_rnotNull");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mrbess_rnotNull(Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(A0)), Shield<SEXP>(Rcpp::wrap(B0)), Shield<SEXP>(Rcpp::wrap(L)), Shield<SEXP>(Rcpp::wrap(X_gen)), Shield<SEXP>(Rcpp::wrap(x_norm)), Shield<SEXP>(Rcpp::wrap(r_max)), Shield<SEXP>(Rcpp::wrap(p0_max)), Shield<SEXP>(Rcpp::wrap(lam)), Shield<SEXP>(Rcpp::wrap(sigma)), Shield<SEXP>(Rcpp::wrap(eps)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(inner_maxit)), Shield<SEXP>(Rcpp::wrap(warm_start)), Shield<SEXP>(Rcpp::wrap(opt)), Shield<SEXP>(Rcpp::wrap(size_max)), Shield<SEXP>(Rcpp::wrap(rho)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<List >(rcpp_result_gen);
    }

}

#endif // RCPP_MrBeSS_RCPPEXPORTS_H_GEN_
