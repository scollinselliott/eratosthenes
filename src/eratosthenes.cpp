#include <Rcpp.h>
using namespace Rcpp;

//' @useDynLib eratosthenes
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericVector gibbs_ad_use_cpp(Rcpp::NumericMatrix marginal, Rcpp::List tpq_list, Rcpp::List taq_list) {

    //block 1 - tpq, block 2 - taq, block3 - use
    Rcpp::NumericMatrix prd ( tpq_list.length(), marginal.ncol()  );
    Rcpp::NumericMatrix dep ( tpq_list.length(), marginal.ncol()  );
    Rcpp::NumericMatrix use ( tpq_list.length(), marginal.ncol()  );

    for (int i = 0; i < tpq_list.length(); i++) {
        prd(i, 0) = marginal(i , 0);
        dep(i, 0) = marginal(i + tpq_list.length() , 0);
        use(i, 0) = marginal(i + 2 * tpq_list.length() , 0);
    }

    for (int k = 1; k < marginal.ncol(); k++) {
        for (int i = 0; i < tpq_list.length(); i++) {
            Rcpp::NumericVector tpq0_ = Rcpp::as<NumericVector>(tpq_list[i]);
            Rcpp::NumericVector tpq1_ = tpq0_[tpq0_ < use(i ,k-1)];
            double L = Rcpp::sample(tpq1_, 1).at(0);
            prd(i, k) = L;
            Rcpp::NumericVector taq0_ = Rcpp::as<NumericVector>(taq_list[i]);
            Rcpp::NumericVector taq1_ = taq0_[taq0_ > use(i ,k-1)];
            double U = Rcpp::sample(taq1_, 1).at(0);
            dep(i, k) = U;
            Rcpp::NumericVector s = Rcpp::runif(1, L, U);
            use(i, k) = s[0];
        }
    }

    for (int i = 0; i < tpq_list.length(); i++) {
        marginal(i , _) = prd(i , _);
        marginal(i + tpq_list.length() , _) = dep(i , _);
        marginal(i + 2 * tpq_list.length() , _) = use(i, _);
    }

    return marginal;
}



//' @useDynLib eratosthenes
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
NumericVector gibbs_ad_use_init_cpp(Rcpp::List tpq_list, Rcpp::List taq_list, int n_samples) {
    Rcpp::NumericMatrix marginal ( tpq_list.length()*3, n_samples  );

    //block 1 - tpq, block 2 - taq, block3 - use
    Rcpp::NumericMatrix prd ( tpq_list.length(), n_samples  );
    Rcpp::NumericMatrix dep ( tpq_list.length(), n_samples  );
    Rcpp::NumericMatrix use ( tpq_list.length(), n_samples  );

    for (int i = 0; i < tpq_list.length(); i++) {
        Rcpp::NumericVector tpq0_ = Rcpp::as<NumericVector>(tpq_list[i]);
        double L = min(tpq0_);
        prd(i, 0) = L;
        Rcpp::NumericVector taq0_ = Rcpp::as<NumericVector>(taq_list[i]);
        double U = max(taq0_);
        dep(i, 0) = U;
        Rcpp::NumericVector s = Rcpp::runif(1, L, U);
        use(i, 0) = s[0];
    }
    for (int k = 1; k < n_samples; k++) {
        for (int i = 0; i < tpq_list.length(); i++) {
            Rcpp::NumericVector tpq0_ = Rcpp::as<NumericVector>(tpq_list[i]);
            Rcpp::NumericVector tpq1_ = tpq0_[tpq0_ < use(i ,k-1)];
            double L = Rcpp::sample(tpq1_, 1).at(0);
            prd(i, k) = L;
            Rcpp::NumericVector taq0_ = Rcpp::as<NumericVector>(taq_list[i]);
            Rcpp::NumericVector taq1_ = taq0_[taq0_ > use(i ,k-1)];
            double U = Rcpp::sample(taq1_, 1).at(0);
            dep(i, k) = U;
            Rcpp::NumericVector s = Rcpp::runif(1, L, U);
            use(i, k) = s[0];
        }
    }

    for (int i = 0; i < tpq_list.length(); i++) {
        marginal(i , _) = prd(i , _);
        marginal(i + tpq_list.length() , _) = dep(i , _);
        marginal(i + 2 * tpq_list.length() , _) = use(i, _);
    }

    return marginal;
}



//' @useDynLib eratosthenes
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::IntegerMatrix quae_postea_matrix_cpp(int elem, Rcpp::List obj) {

    Rcpp::IntegerMatrix qp(elem,elem);

    for (int i = 0; i < obj.length(); i++) {
        Rcpp::IntegerVector sequentia = obj[i];
        for (int j = 0; j < (sequentia.length()-1); j++) {
            int elem1 = sequentia[j];
            for (int k = (j+1); k < sequentia.length(); k++) {
                int elem2 = sequentia[k];
                qp((elem1-1), (elem2-1)) = 1;
            }
        }
    }

    int diff = 0;
    for (int i = 0; i < elem; i++) {
        Rcpp::IntegerVector tmp = qp( i , _ );
        Rcpp::IntegerVector checked (elem);
        diff = 0;
        for (int p = 0; p < elem; p++) {
            int diff1 = tmp[p] - checked[p];
            diff += diff1 * diff1;
        }

        while (diff > 0) {
            for (int j = 0; j < elem; j++) {
                int tmpj = tmp[j];
                int checkedj = checked[j];
                if (tmpj == 1) {
                    if (checkedj == 0) {
                        for (int jj = 0; jj < elem; jj++) {
                            int tmp2 = qp(j, jj);
                            if (tmp2 == 1) {
                                qp(i,jj) = 1;
                            }
                        checked[j] = 1;
                        }
                    }
                }
            }

            tmp = qp( i , _ );
            
            diff = 0;
            for (int p = 0; p < elem; p++) {
                int diff1 = tmp[p] - checked[p];
                diff += diff1 * diff1;
            }
        }
    }

return qp;
}



//' @useDynLib eratosthenes
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::IntegerMatrix quae_antea_matrix_cpp(int elem, Rcpp::List obj) {

    Rcpp::IntegerMatrix qa(elem,elem);

    for (int i = 0; i < obj.length(); i++) {
        Rcpp::IntegerVector sequentia = obj[i];
        for (int j = 1; j < (sequentia.length()); j++) {
            int elem1 = sequentia[j];
            for (int k = 0; k < j; k++) {
                int elem2 = sequentia[k];
                qa((elem1-1), (elem2-1)) = 1;
            }
        }
    }

    int diff = 0;
    for (int i = 0; i < elem; i++) {
        Rcpp::IntegerVector tmp = qa( i , _ );
        Rcpp::IntegerVector checked (elem);
        diff = 0;
        for (int p = 0; p < elem; p++) {
            int diff1 = tmp[p] - checked[p];
            diff += diff1 * diff1;
        }

        while (diff > 0) {
            for (int j = 0; j < elem; j++) {
                int tmpj = tmp[j];
                int checkedj = checked[j];
                if (tmpj == 1) {
                    if (checkedj == 0) {
                        for (int jj = 0; jj < elem; jj++) {
                            int tmp2 = qa(j, jj);
                            if (tmp2 == 1) {
                                qa(i,jj) = 1;
                            }
                        checked[j] = 1;
                        }
                    }
                }
            }

            tmp = qa( i , _ );
            
            diff = 0;
            for (int p = 0; p < elem; p++) {
                int diff1 = tmp[p] - checked[p];
                diff += diff1 * diff1;
            }
        }
    }

return qa;
}



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



//' @useDynLib eratosthenes
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::NumericMatrix gibbs_ad_init_cpp(Rcpp::NumericMatrix a, Rcpp::IntegerVector tpq_idx, Rcpp::IntegerMatrix phi, Rcpp::List phiList, Rcpp::IntegerVector taq_idx, Rcpp::IntegerMatrix psi, Rcpp::List psiList, Rcpp::IntegerVector prc_idx) {

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
    


//' @useDynLib eratosthenes
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::NumericVector gibbs_ad_initial_cpp(Rcpp::NumericVector a, Rcpp::IntegerVector tpq_idx, Rcpp::IntegerMatrix phi, Rcpp::List phiList, Rcpp::IntegerVector taq_idx, Rcpp::IntegerMatrix psi, Rcpp::List psiList, Rcpp::IntegerVector prc_idx, int subsample) {

// x is the initial value to assign
// only use samples less than x
// use assigned to keep track of assigned values

Rcpp::IntegerVector assigned (a.length());

// where assigned indicate 1
for (int i = 0; i < tpq_idx.length(); i++) {
    assigned[tpq_idx[i] - 1] = 1;
}
for (int i = 0; i < taq_idx.length(); i++) {
    assigned[taq_idx[i] - 1] = 1;
}
assigned[assigned.length() - 2] = 1;
assigned[assigned.length() - 1] = 1;


// initial value of x

// for each of the relative events
for (int x = 0; x < prc_idx.size(); x++) {
    int idx = prc_idx[x] - 1;

    int plength = 0;
    for (int ii = 0; ii < phi.ncol(); ii++) {
        if (phi( idx , ii ) == 1 ) {
            if (assigned[ii] == 1) {
                plength += 1;
            }
        }
    }

    Rcpp::IntegerVector postea (plength);

    int j = 0;
    for (int ii = 0; ii < phi.ncol(); ii++) {
        if (phi( idx , ii ) == 1 ) {
            if (assigned[ii] == 1) {
                postea[j] = ii;
                j += 1;
            }
        }
    }

    Rcpp::NumericVector post (plength);

    for (int ii = 0; ii < plength; ii++) {
        post[ii] = a( postea[ii] );
    }

    double U = Rcpp::min(post);

    int alength = 0;
    for( int ii = 0; ii < psi.ncol(); ii++) {
        if (psi( idx , ii ) == 1 ) {
            if (assigned[ii] == 1) {
                alength += 1;
            }
        }
    }

    Rcpp::IntegerVector antea (alength);

    j = 0;
    for (int ii = 0; ii < psi.ncol(); ii++) {
        if (psi( idx , ii ) == 1 ) {
            if (assigned[ii] == 1) {
                antea[j] = ii;
                j += 1;
            }
        }
    }

    Rcpp::NumericVector ante (alength);

    for (int ii = 0; ii < alength; ii++) {
        ante[ii] = a( antea[ii] );
    }
    double L = Rcpp::max(ante);

    NumericVector s = Rcpp::runif(1, L, U);

    a(idx) = s[0];
    assigned[idx] = 1;

    for (int m = 0; m < subsample; m++) {
    
        for (int i = 0; i < phiList.length(); i++) {
            int idx = tpq_idx[i] - 1;
    
            int plength = 0;
            for(int ii = 0; ii < phi.ncol(); ii++) {
                if (phi( idx , ii ) == 1 ) {
                    if (assigned[ii] == 1) {
                        plength += 1;
                    }
                }
            }
    
            Rcpp::IntegerVector postea (plength);
    
            int j = 0;
            for(int ii = 0; ii < phi.ncol(); ii++) {
                if (phi( idx , ii ) == 1 ) {
                    if (assigned[ii] == 1) {
                        postea[j] = ii;
                        j += 1;
                    }
                }
            }
    
            Rcpp::NumericVector post (plength);
            for (int ii = 0; ii < plength; ii++) {
                post[ii] = a( postea[ii] );
            }
            double U = Rcpp::min(post);
    
            Rcpp::List phi0 = Rcpp::as<List>(phiList[i]);
            Rcpp::NumericVector isamples0 = as<NumericVector>(phi0["samples"]);
            Rcpp::NumericVector isamples = isamples0[isamples0 < U];
            double tpq0 = Rcpp::sample(isamples, 1).at(0);
    
            a(idx) = tpq0;
    
        }
        for (int i = 0; i < psiList.length(); i++) {
            int idx = taq_idx[i] - 1;
    
            int sum = 0;
            for(int ii = 0; ii < psi.ncol(); ii++) {
                if (psi( idx , ii ) == 1 ) {
                    if (assigned[ii] == 1) {
                        sum += 1;
                    }
                }
            }
    
            Rcpp::IntegerVector antea (sum);
    
            int j = 0;
            for(int ii = 0; ii < psi.ncol(); ii++) {
                if (psi( idx , ii ) == 1 ) {
                    if (assigned[ii] == 1) {
                        antea[j] = ii;
                        j += 1;
                    }
                }
            }
    
            Rcpp::NumericVector ante (sum);
            for (int ii = 0; ii < sum; ii++) {
                ante[ii] = a( antea[ii] );
            }
            double L = Rcpp::max(ante);
    
            Rcpp::List psi0 = Rcpp::as<List>(psiList[i]);
            Rcpp::NumericVector isamples0 = as<NumericVector>(psi0["samples"]);
            Rcpp::NumericVector isamples = isamples0[isamples0 > L];
            double taq0 = Rcpp::sample(isamples, 1).at(0);
    
            a(idx) = taq0;
    
        }

        // relative events
        for (int i = 0; i < prc_idx.size(); i++) {
            int idx = prc_idx[i] - 1;
    
            int plength = 0;
            for(int ii = 0; ii < phi.ncol(); ii++) {
                if (phi( idx , ii ) == 1 ) {
                    if (assigned[ii] == 1) {
                        plength += 1;
                    }
                }
            }
    
            Rcpp::IntegerVector postea (plength);
    
            int j = 0;
            for(int ii = 0; ii < phi.ncol(); ii++) {
                if (phi( idx , ii ) == 1 ) {
                    if (assigned[ii] == 1) {
                        postea[j] = ii;
                        j += 1;
                    }
                }
            }
    
            Rcpp::NumericVector post (plength);
            for (int ii = 0; ii < plength; ii++) {
                post[ii] = a( postea[ii] );
            }
            double U = Rcpp::min(post);
    
            int alength = 0;
            for(int ii = 0; ii < psi.ncol(); ii++) {
                if (psi( idx , ii ) == 1 ) {
                    if (assigned[ii] == 1) {
                        alength += 1;
                    }
                }
            }
    
            Rcpp::IntegerVector antea (alength);
    
            j = 0;
            for(int ii = 0; ii < psi.ncol(); ii++) {
                if (psi( idx , ii ) == 1 ) {
                    if (assigned[ii] == 1) {
                        antea[j] = ii;
                        j += 1;
                    }
                }
            }
    
            Rcpp::NumericVector ante (alength);
            for (int ii = 0; ii < alength; ii++) {
                ante[ii] = a( antea[ii] );
            }
            double L = Rcpp::max(ante);
    
            NumericVector s = Rcpp::runif(1, L, U);
    
            a(idx) = s[0];

        }
    }
}
return a;
}



