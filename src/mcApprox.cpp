#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double pfhtrbm_c_db(double t, double x0, double nu, double kappa, double sigma){
  double k, smt, lambda_k, c_k, temp_k, ret;
  if (t <= 0) {
    ret = 0;
    return(ret);
  }

  k = 1.0; smt = 0.0; lambda_k = 0.0; c_k = 0.0; temp_k = 1.0;

  while (std::abs(temp_k) > 0.0) {
    lambda_k = pow((2 * k - 1), 2) * pow(sigma, 2) *
      pow(M_PI, 2) / (8 * pow((kappa - nu), 2));
    c_k = pow(-1, k + 1) * 4 / ((2 * k - 1) * M_PI) *
      cos((2 * k - 1) * M_PI * (kappa - x0) / (2 * (kappa - nu)));
    temp_k = c_k * exp(- lambda_k * t);
    smt = smt + temp_k;
    k++;
  }

  // dealing with numerical error
  if (0 <= (1-smt) & (1-smt) <= 1) {
    ret = 1-smt;
  } else if (1 < (1-smt)) {
    ret = 1;
  } else {
    ret = 0;
  }
  return(ret);
}

// [[Rcpp::export]]
double lklhSubj(NumericVector t,
                NumericVector event,
                double x0,
                double nu,
                double kappa,
                double sigma){
  NumericVector t_vec = clone(t);
  double Lklh = 1, Lklh_temp = 1;
  int n = t_vec.size();
  NumericVector t_vec_end = t_vec + 0.5;
  NumericVector t_vec_start = t_vec - 0.5;

  for (int i = 0; i < n; i++){
    if (event[i] == 1) {
      Lklh_temp = pfhtrbm_c_db(t_vec_end[i], x0, nu, kappa, sigma) -
        pfhtrbm_c_db(t_vec_start[i], x0, nu, kappa, sigma);
    } else {
      Lklh_temp = 1 - pfhtrbm_c_db(t_vec_end[i], x0, nu, kappa, sigma);
    }
    Lklh = Lklh * Lklh_temp;
  }
  return(Lklh);
}


// [[Rcpp::export]]
NumericVector obsLklh_c(NumericVector t,
                      NumericVector event,
                      arma::vec covar_sigma,
                      arma::vec covar_kappa,
                      arma::mat mcmc_sigma,
                      arma::mat mcmc_kappa,
                      arma::mat mcmc_gamma,
                      arma::mat mcmc_theta,
                      arma::vec mc_z1,
                      arma::vec mc_z2,
                      String frailty,
                      double x0,
                      double nu){
  int mcmcsize = mcmc_sigma.n_rows;
  int mcsize = mc_z1.n_elem;
  NumericVector result(mcmcsize);

  // for safe
  NumericVector t_vec = clone(t);
  NumericVector event_vec = clone(event);

  double sigmaa, kappaa;
  NumericVector obsLklh_allMC(mcsize);

  if (frailty == "correlated") {
    arma::mat alphas = mcmc_kappa * covar_kappa;
    arma::mat betas = mcmc_sigma * covar_sigma;

    for(int i = 0; i < mcmcsize; i++){
      arma::vec mc_z1_i_iter = mc_z1;
      mc_z1_i_iter = sqrt(mcmc_theta(i, 0)) * mc_z1_i_iter;
      arma::vec mc_z2_i_iter = mc_z2;
      mc_z2_i_iter = sqrt(mcmc_theta(i, 1)) * mc_z2_i_iter;

      for(int j = 0; j < mcsize; j++){
        kappaa = x0 + exp(alphas(i, 0) + mc_z2_i_iter[j] +
          mcmc_gamma(i, 0) * mc_z1_i_iter[j]);
        sigmaa = exp(betas(i, 0) + mc_z1_i_iter[j]);
        obsLklh_allMC[j] = lklhSubj(t_vec, event_vec, x0, nu, kappaa, sigmaa);
      }
      result[i] = mean(obsLklh_allMC);
    }

  }else if (frailty == "independent"){
    arma::mat alphas = mcmc_kappa * covar_kappa;
    arma::mat betas = mcmc_sigma * covar_sigma;

    for(int i = 0; i < mcmcsize; i++){
      arma::vec mc_z1_i_iter = mc_z1;
      mc_z1_i_iter = sqrt(mcmc_theta(i, 0)) * mc_z1_i_iter;
      arma::vec mc_z2_i_iter = mc_z2;
      mc_z2_i_iter = sqrt(mcmc_theta(i, 1)) * mc_z2_i_iter;

      for(int j = 0; j <mcsize; j++){
        kappaa = x0 + exp(alphas(i, 0) + mc_z2_i_iter[j]);
        sigmaa = exp(betas(i, 0)  + mc_z1_i_iter[j]);
        obsLklh_allMC[j] = lklhSubj(t_vec, event_vec, x0, nu, kappaa, sigmaa);
      }
      result[i] = mean(obsLklh_allMC);
    }

  } else if (frailty == "shared") {
    arma::mat alphas = mcmc_kappa * covar_kappa;
    arma::mat betas = mcmc_sigma * covar_sigma;

    for(int i = 0; i < mcmcsize; i++){
      arma::vec mc_z1_i_iter = mc_z1;
      mc_z1_i_iter = sqrt(mcmc_theta(i, 0)) * mc_z1_i_iter;
      // arma::vec mc_z2_i_iter = mc_z2;
      // mc_z2_i_iter = sqrt(mcmc_theta(i, 1)) * mc_z2_i_iter;

      for(int j = 0; j < mcsize; j++){
        kappaa = x0 + exp(alphas(i, 0) + mcmc_gamma(i, 0) * mc_z1_i_iter[j]);
        sigmaa = exp(betas(i, 0)  + mc_z1_i_iter[j]);
        obsLklh_allMC[j] = lklhSubj(t_vec, event_vec, x0, nu, kappaa, sigmaa);
      }
      result[i] = mean(obsLklh_allMC);
    }

  } else {
    arma::mat alphas = mcmc_kappa * covar_kappa;
    arma::mat betas = mcmc_sigma * covar_sigma;

    for(int i = 0; i < mcmcsize; i++){
      kappaa = x0 + exp(alphas(i, 0));
      sigmaa = exp(betas(i, 0));
      result[i] = lklhSubj(t_vec, event_vec, x0, nu, kappaa, sigmaa);
    }
  }
  return(result);
}


