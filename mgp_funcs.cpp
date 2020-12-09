#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;


// [[Rcpp::export]]
Rcpp::NumericMatrix eta_lin(arma::mat lambda, arma::vec ps, int k, int n, arma::mat Y) {
    // --- UPDATE eta --- //
    arma::mat Lmsg = lambda.each_col() % ps;
    arma::mat Veta1 = eye<arma::mat>(k, k) + Lmsg.t() * lambda;
    arma::mat S = inv(trimatu(chol(Veta1)));
    arma::mat Veta = S * S.t();
    arma::mat Meta = Y * Lmsg * Veta;
    arma::mat noise = randn<arma::mat>(n, k);
    arma::mat eta = Meta + noise * S.t();
    return Rcpp::wrap(eta);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix lam_lin(arma::mat eta, arma::mat Plam, arma::vec ps, int k, int p, arma::mat Y) {
    // --- UPDATE lambda --- //
    arma::mat lambda(p, k);
    arma::mat eta2 = eta.t() * eta;    // prepare eta crossproduct before the loop
    for (int j = 0; j < p; ++j) {
        arma::mat Llamt = trimatu(chol(diagmat(Plam.row(j)) + ps(j) * eta2));
        arma::mat Llam = trimatl(Llamt.t());
        lambda.row(j) = (solve(Llamt, randn<arma::vec>(k)) +
            solve(Llamt, solve(Llam, ps(j) * eta.t() * Y.col(j)))).t();
    }
    return Rcpp::wrap(lambda);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix psi_mg(arma::mat lambda, arma::vec tauh, arma::vec ps, int k, int p, double df) {
    // --- UPDATE psihj --- //
    arma::mat psijh(p, k);
    arma::mat lambda_sq = square(lambda);
    arma::mat shape = lambda_sq.each_row() % tauh.t();
    for (int l = 0; l < p; l++) {
        for (int j = 0; j < k; j++) {
            psijh(l, j) = randg(distr_param(df / 2 + 0.5, 1 / (df / 2 + shape(l, j) / 2)));
        }
    }
    return Rcpp::wrap(psijh);
}

// [[Rcpp::export]]
Rcpp::NumericVector del_mg(arma::mat lambda, arma::mat psijh, arma::vec tauh, arma::vec delta, int k, int p,
    double ad1, double bd1, double ad2, double bd2) {
    // --- UPDATE delta & tauh --- //
    arma::mat matr = psijh % square(lambda);
    double ad = ad1 + 0.5 * p * k;
    double bd = bd1 + 0.5 / delta(0) * sum(tauh.t() % sum(matr, 0));
    delta(0) = randg(distr_param(ad, 1 / bd));
    tauh = cumprod(delta);

    for (int h = 1; h < k; ++h) {
        ad = ad2 + 0.5 * p * (k - h);
        arma::vec tauh_sub = tauh.subvec(h, k - 1);
        bd = bd2 + 0.5 / delta(h) * sum(tauh_sub.t() % sum(matr.cols(h, k - 1), 0));
        delta(h) = randg(distr_param(ad, 1 / bd));
        tauh = cumprod(delta);
    }
    return Rcpp::wrap(delta);
}

// [[Rcpp::export]]
Rcpp::NumericVector sig_lin(arma::mat lambda, arma::mat eta, int k, int p, int n, arma::mat Y,
    double as, double bs) {
    // --- UPDATE sigma --- //
    arma::mat Ytil = Y - eta * lambda.t();
    rowvec bsvec = bs + 0.5 * sum(square(Ytil));
    auto lam = [as, n](double val) {return randg<double>(distr_param(as + 0.5 * n, 1 / val)); };
    return Rcpp::wrap((bsvec.transform(lam)).t());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix plm_mg(arma::mat psijh, arma::vec tauh) {
    // --- UPDATE Plam --- //
    return Rcpp::wrap(psijh.each_row() % tauh.t());
}
