#include <Rcpp.h>
using namespace Rcpp;

//' @useDynLib eratosthenes
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::NumericMatrix gibbs_ad_cpp(Rcpp::NumericMatrix a, Rcpp::IntegerVector tpq_idx, Rcpp::IntegerMatrix phi, Rcpp::List phiList, Rcpp::IntegerVector taq_idx, Rcpp::IntegerMatrix psi, Rcpp::List psiList, Rcpp::IntegerVector prc_idx) {

for (int m = 1; m < a.ncol(); m++) {
    a( _ , m ) = a ( _ , m-1);
    for (int i = 0; i < phiList.length(); i++) {
        int idx = tpq_idx[i] - 1;

        int plength = 0;
        for(int ii = 0; ii < phi.ncol(); ii++) {
            if (phi( idx , ii ) == 1 ) {
                plength += 1;
            }
        }

        Rcpp::IntegerVector postea (plength);

        int j = 0;
        for(int ii = 0; ii < phi.ncol(); ii++) {
            if (phi( idx , ii ) == 1 ) {
                postea[j] = ii;
                j += 1;
            }
        }

        Rcpp::NumericVector post (plength);
        for (int ii = 0; ii < plength; ii++) {
            post[ii] = a( postea[ii], m );
        }
        double U = Rcpp::min(post);

        Rcpp::List phi0 = Rcpp::as<List>(phiList[i]);
        Rcpp::NumericVector isamples0 = as<NumericVector>(phi0["samples"]);
        Rcpp::NumericVector isamples = isamples0[isamples0 < U];
        double tpq0 = Rcpp::sample(isamples, 1).at(0);

        a( idx , m ) = tpq0;

    }
    for (int i = 0; i < psiList.length(); i++) {
        int idx = taq_idx[i] - 1;

        int sum = 0;
        for(int ii = 0; ii < psi.ncol(); ii++) {
            if (psi( idx , ii ) == 1 ) {
                sum += 1;
            }
        }

        Rcpp::IntegerVector antea (sum);

        int j = 0;
        for(int ii = 0; ii < psi.ncol(); ii++) {
            if (psi( idx , ii ) == 1 ) {
                antea[j] = ii;
                j += 1;
            }
        }

        Rcpp::NumericVector ante (sum);
        for (int ii = 0; ii < sum; ii++) {
            ante[ii] = a( antea[ii], m );
        }
        double L = Rcpp::max(ante);

        Rcpp::List psi0 = Rcpp::as<List>(psiList[i]);
        Rcpp::NumericVector isamples0 = as<NumericVector>(psi0["samples"]);
        Rcpp::NumericVector isamples = isamples0[isamples0 > L];
        double taq0 = Rcpp::sample(isamples, 1).at(0);

        a( idx , m ) = taq0;

    }
    for (int i = 0; i < prc_idx.size(); i++) {
        int idx = prc_idx[i] - 1;

        int plength = 0;
        for(int ii = 0; ii < phi.ncol(); ii++) {
            if (phi( idx , ii ) == 1 ) {
                plength += 1;
            }
        }

        Rcpp::IntegerVector postea (plength);

        int j = 0;
        for(int ii = 0; ii < phi.ncol(); ii++) {
            if (phi( idx , ii ) == 1 ) {
                postea[j] = ii;
                j += 1;
            }
        }

        Rcpp::NumericVector post (plength);
        for (int ii = 0; ii < plength; ii++) {
            post[ii] = a( postea[ii], m );
        }
        double U = Rcpp::min(post);



        int alength = 0;
        for(int ii = 0; ii < psi.ncol(); ii++) {
            if (psi( idx , ii ) == 1 ) {
                alength += 1;
            }
        }

        Rcpp::IntegerVector antea (alength);

        j = 0;
        for(int ii = 0; ii < psi.ncol(); ii++) {
            if (psi( idx , ii ) == 1 ) {
                antea[j] = ii;
                j += 1;
            }
        }

        Rcpp::NumericVector ante (alength);
        for (int ii = 0; ii < alength; ii++) {
            ante[ii] = a( antea[ii], m );
        }
        double L = Rcpp::max(ante);

        NumericVector s = Rcpp::runif(1, L, U);

        a( idx , m ) = s[0];

    }
}

return a;
}



