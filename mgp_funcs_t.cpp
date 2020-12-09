#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;

struct tree_val {
    arma::mat thetaminus;
    arma::mat rmatminus;
    arma::mat thetaplus;
    arma::mat rmatplus;
    arma::mat thetaprime;
    int nprime;
    int sprime;


    tree_val& operator=(const tree_val& z) {
        if (this == &z) {
            return *this;
        }

        thetaminus = z.thetaminus;
        rmatminus = z.rmatminus;
        thetaplus = z.thetaplus;
        rmatplus = z.rmatplus;
        thetaprime = z.thetaprime;
        nprime = z.nprime;
        sprime = z.sprime;

        return *this;
    }

    tree_val& operator<<=(const tree_val& z) {
        if (this == &z) {
            return *this;
        }

        thetaminus = z.thetaminus;
        rmatminus = z.rmatminus;
        thetaprime = z.thetaprime;
        nprime = z.nprime;
        sprime = z.sprime;

        return *this;
    }

    tree_val& operator>>=(const tree_val& z) {
        if (this == &z) {
            return *this;
        }

        thetaplus = z.thetaplus;
        rmatplus = z.rmatplus;
        thetaprime = z.thetaprime;
        nprime = z.nprime;
        sprime = z.sprime;

        return *this;
    }
};


tree_val eta_buildtree_2(arma::mat& eta, arma::mat& rmat, int u, int v, int j, double epss,
    arma::mat& lambda, arma::mat& Y, arma::vec& ps, int n, int p, int k, double nu,
    arma::mat& Ilamtpslam, arma::mat& lamtpsy, arma::vec& treecount) {
    //Rcpp::Rcout << "Tree constructed" << std::endl;
    treecount += 1;
    if (j == 0) {
        tree_val tree_out;

        double vepss = v * epss;

        arma::vec nuvec(n); nuvec.fill(nu);

        //leapfrog
        //1
        arma::mat grad = -eta * Ilamtpslam + lamtpsy;
        arma::vec mahaleta = sum(eta % eta, 1);
        arma::mat ybar = Y - eta * lambda.t();
        arma::mat ybar2 = ybar % ybar;
        ybar2.each_row() %= ps.t();
        arma::vec mahaly = sum(ybar2, 1);
        arma::vec temp = 1 / (nuvec + mahaly + mahaleta);
        temp *= nu + p + k;
        grad.each_col() %= temp;
        arma::mat rmatprime = rmat + (0.5 * vepss) * grad;

        //2
        arma::mat etaprime = eta + vepss * rmatprime;

        //3
        grad = -etaprime * Ilamtpslam + lamtpsy;
        mahaleta = sum(etaprime % etaprime, 1);
        ybar = Y - etaprime * lambda.t();
        ybar2 = ybar % ybar;
        ybar2.each_row() %= ps.t();
        mahaly = sum(ybar2, 1);
        temp = 1 / (nuvec + mahaly + mahaleta);
        temp *= nu + p + k;
        grad.each_col() %= temp;
        rmatprime += (0.5 * vepss) * grad;

        // etaprime target value
        arma::vec targetvec = (-0.5 * (nu + p + k)) * log((1 / nu) * (nuvec + mahaly + mahaleta));
        double targettotal = sum(targetvec);

        double temp7 = arma::accu(rmatprime % rmatprime);
        bool nprimebool = u <= exp(targettotal - 0.5 * temp7);

        tree_out.nprime = int(nprimebool);
        bool sprimebool = targettotal - 0.5 * temp7 > log(u) - 100000;
        //bool sprimebool = 1;
        tree_out.sprime = int(sprimebool);

        tree_out.thetaminus = etaprime;
        tree_out.rmatminus = rmatprime;
        tree_out.thetaplus = etaprime;
        tree_out.rmatplus = rmatprime;
        tree_out.thetaprime = etaprime;

        return tree_out;
    }
    else {
        tree_val tree_out2 = eta_buildtree_2(eta, rmat, u, v, j - 1, epss, lambda, Y, ps, n, p, k, nu, Ilamtpslam, lamtpsy, treecount);

        arma::mat etaprime = tree_out2.thetaprime;
        int nprime = tree_out2.nprime;
        int sprime = tree_out2.sprime;

        if (sprime == 1) {
            if (v == -1) {
                arma::mat etaminus = tree_out2.thetaminus;
                arma::mat rmatminus = tree_out2.rmatminus;
                tree_out2 <<= eta_buildtree_2(etaminus, rmatminus, u, v, j - 1, epss, lambda, Y, ps, n, p, k, nu, Ilamtpslam, lamtpsy, treecount);
            }
            else {
                arma::mat etaplus = tree_out2.thetaplus;
                arma::mat rmatplus = tree_out2.rmatplus;
                tree_out2 >>= eta_buildtree_2(etaplus, rmatplus, u, v, j - 1, epss, lambda, Y, ps, n, p, k, nu, Ilamtpslam, lamtpsy, treecount);
            }
            //arma::mat etatwoprime = tree_out2.thetaprime;
            int ntwoprime = tree_out2.nprime;

            double prob = ntwoprime / (nprime + ntwoprime);
            double uu = randu<double>();
            if (!(uu < prob)) {
                tree_out2.thetaprime = etaprime;
            }

            bool criterion1 = arma::accu((tree_out2.thetaplus - tree_out2.thetaminus) % tree_out2.rmatplus) >= 0;
            bool criterion2 = arma::accu((tree_out2.thetaplus - tree_out2.thetaminus) % tree_out2.rmatminus) >= 0;
            if (!(criterion1 && criterion2)) {
                tree_out2.sprime = 0;
            }
            tree_out2.nprime += nprime;
        }

        return tree_out2;
    }
}


// [[Rcpp::export]]
Rcpp::NumericMatrix eta_NUTS_2(arma::mat eta, arma::mat lambda, arma::vec ps, int k, int n, int p, arma::mat Y, double nu, double epss, arma::vec& treecount) {
    //Rcpp::Rcout << "NUTS update for eta ###########################" << std::endl;
    arma::mat r0 = randn<arma::mat>(n, k);

    // 공통
    arma::mat Lmsg = lambda.each_col() % ps;
    arma::mat Ilamtpslam = eye<arma::mat>(k, k) + (Lmsg.t() * lambda); //공통
    arma::mat lamtpsy = Y * Lmsg; //공통

    // eta target value
    arma::vec mahaleta = sum(eta % eta, 1);
    arma::mat ybar = Y - eta * lambda.t();
    arma::mat ybar2 = ybar % ybar;
    ybar2.each_row() %= ps.t();
    arma::vec mahaly = sum(ybar2, 1);
    arma::vec nuvec(n); nuvec.fill(nu);
    arma::vec targetvec = (-0.5 * (nu + p + k)) * log((1 / nu) * nuvec + mahaly + mahaleta);
    double targettotal = sum(targetvec);

    double temp7 = arma::accu(r0 % r0);
    double upper = exp(targettotal - 0.5 * temp7);
    double u = randu<double>(); u *= upper;

    tree_val tree_vals;

    tree_vals.thetaminus = eta;
    tree_vals.thetaplus = eta;
    tree_vals.thetaprime = eta;
    tree_vals.rmatminus = r0;
    tree_vals.rmatplus = r0;
    int j = 0;
    arma::mat etanew = eta;
    tree_vals.nprime = 1;
    tree_vals.sprime = 1;


    while (tree_vals.sprime == 1) {
        int prevn = tree_vals.nprime;

        double vdouble = randu<double>();
        int v = 1;
        if (vdouble < 0.5) {
            v = -1;
        }

        if (v == -1) {
            arma::mat etaminus = tree_vals.thetaminus;
            arma::mat rmatminus = tree_vals.rmatminus;
            tree_vals <<= eta_buildtree_2(etaminus, rmatminus, u, v, j, epss, lambda, Y, ps, n, p, k, nu, Ilamtpslam, lamtpsy, treecount);
        }
        else {
            arma::mat etaplus = tree_vals.thetaplus;
            arma::mat rmatplus = tree_vals.rmatplus;
            tree_vals >>= eta_buildtree_2(etaplus, rmatplus, u, v, j, epss, lambda, Y, ps, n, p, k, nu, Ilamtpslam, lamtpsy, treecount);
        }

        if (tree_vals.sprime == 1) {
            double updateprob = tree_vals.nprime / prevn;
            double uu = randu<double>();
            if (uu < updateprob) {
                etanew = tree_vals.thetaprime;
            }
        }
        tree_vals.nprime += prevn;

        bool criterion1 = arma::accu((tree_vals.thetaplus - tree_vals.thetaminus) % tree_vals.rmatplus) >= 0;
        bool criterion2 = arma::accu((tree_vals.thetaplus - tree_vals.thetaminus) % tree_vals.rmatminus) >= 0;
        if (!(criterion1 && criterion2)) {
            tree_vals.sprime = 0;
        }
        j += 1;
    }

    return Rcpp::wrap(etanew);
}


// [[Rcpp::export]]
Rcpp::NumericVector gam_gibbs_t(arma::mat lambda, arma::mat eta, arma::vec ps, int n, int k, double p, arma::mat Y, double nu) {
    // --- UPDATE gamma --- //
    arma::mat temp11 = Y - eta * lambda.t();
    arma::mat temp115 = temp11 % temp11;
    arma::mat temp12 = temp115.each_row() % ps.t();
    arma::vec temp13 = sum(temp12, 1);
    rowvec rate = zeros<arma::rowvec>(n);
    for (int i = 0; i < n; i++) {
        double mahal = temp13(i) + arma::dot(eta.row(i).t(), eta.row(i).t());
        rate[i] = 0.5 * (nu + mahal);
    }
    auto lam = [nu, p, k](double val) {return randg<double>(distr_param(0.5 * (nu + p + k), 1 / val)); };
    return Rcpp::wrap((rate.transform(lam)).t());
}


// [[Rcpp::export]]
Rcpp::NumericMatrix lam_lin_t(arma::mat eta, arma::mat Plam, arma::vec ps, int k, int p, arma::mat Y, arma::vec gam) {
    // --- UPDATE lambda --- //
    arma::mat lambda(p, k);
    arma::mat eta22 = eta.each_col() % gam;
    arma::mat eta2 = eta.t() * eta22;    // prepare eta crossproduct before the loop

    for (int j = 0; j < p; j++) {
        arma::mat temp = diagmat(Plam.row(j));
        temp += ps(j) * eta2;

        arma::mat Llamt = trimatu(chol(symmatu(temp)));
        arma::mat Llam = trimatl(Llamt.t());
        lambda.row(j) = (solve(Llamt, randn<arma::vec>(k)) +
            solve(Llamt, solve(Llam, ps(j) * eta.t() * diagmat(gam) * Y.col(j)))).t();
    }
    return Rcpp::wrap(lambda);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix psi_mg_t(arma::mat lambda, arma::vec tauh, arma::vec ps, int k, int p, double df) {
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
Rcpp::NumericVector del_mg_t(arma::mat lambda, arma::mat psijh, arma::vec tauh, arma::vec delta, int k, int p,
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
Rcpp::NumericVector sig_lin_t(arma::mat lambda, arma::mat eta, int k, int p, int n, arma::mat Y,
    double as, double bs, arma::vec gam) {
    // --- UPDATE sigma --- //
    arma::mat Ytil = Y - eta * lambda.t();
    arma::mat Ytilgam = Ytil.each_col() % sqrt(gam);
    rowvec bsvec = bs + 0.5 * sum(square(Ytilgam));
    auto lam = [as, n](double val) {return randg<double>(distr_param(as + 0.5 * n, 1 / val)); };

    return Rcpp::wrap((bsvec.transform(lam)).t());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix plm_mg_t(arma::mat psijh, arma::vec tauh) {
    // --- UPDATE Plam --- //
    return Rcpp::wrap(psijh.each_row() % tauh.t());
}

// [[Rcpp::export]]
double a1_mg_t(double a1, arma::vec& delta, double a1sigma, arma::vec& rejectcount1) {
    double delta1 = delta(0);
    double a1prime = 0;
    while (a1prime <= 2) {
        double increment = randn<double>();
        a1prime = a1 + a1sigma * increment;
    }
    double logR = lgamma(a1) - lgamma(a1prime);
    logR += (a1prime - a1) * log(delta1);
    logR += log(a1prime) - log(a1);
    logR += (a1 - a1prime);
    double R = exp(logR);
    R *= (1 - normcdf((2 - a1) / (a1sigma))) / (1 - normcdf((2 - a1prime) / (a1sigma)));

    double u = randu<double>();
    if (u < R) {
        return a1prime;
    }
    else {
        rejectcount1 += 1;
        return a1;
    }
}

// [[Rcpp::export]]
double a2_mg_t(double a2, arma::vec& delta, int k, double a2sigma, arma::vec& rejectcount2) {
    double a2prime = 0;
    while (a2prime <= 3) {
        double increment = randn<double>();
        a2prime = a2 + a2sigma * increment;
    }
    double logR = (k - 1) * (lgamma(a2) - lgamma(a2prime));
    logR += (a2prime - a2) * (arma::accu(log(delta)) - log(delta(0)));
    logR += log(a2prime) - log(a2);
    logR += a2 - a2prime;
    double R = exp(logR);
    R *= (1 - normcdf((3 - a2) / (a2sigma))) / (1 - normcdf((3 - a2prime) / (a2sigma)));

    double u = randu<double>();
    if (u < R) {
        return a2prime;
    }
    else {
        rejectcount2 += 1;
        return a2;
    }
}

