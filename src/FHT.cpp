#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dFHT_c(NumericVector t, double x0, double nu, double kappa, double sigma){
  int n = t.size();
  NumericVector ret(n);
  double k, smt, lambda_k, c_k, temp_k;
  
  for (int i = 0; i < n; ++i){
    // check if out of support
    if (t[i] <= 0){
      ret[i] = 0;
      continue;
    }
    
    // reset for next element in vector t
    k = 1.0; smt = 0.0; lambda_k = 0.0; c_k = 0.0; temp_k = 1.0;
    
    // do infinite sum of series 
    while (std::abs(temp_k) > 0.0) {
      lambda_k = pow((2 * k - 1), 2) * pow(sigma, 2) * 
        pow(M_PI, 2) / (8 * pow(kappa - nu, 2));
      c_k = pow(-1, k + 1) * 4  / ((2 * k - 1) * M_PI) * 
        cos((2 * k - 1) * M_PI * (kappa - x0) / (2 * (kappa - nu)));
      temp_k = c_k * lambda_k * exp(- lambda_k * t[i]);
      smt = smt + temp_k;
      k++;
    }
    
    // dealing with numerical error 
    // might exist numerical error such that smt < 0
    if (smt < 0.0) smt = 0.0;
    ret[i] = smt;
  }
  return(ret);
}


// [[Rcpp::export]]
NumericVector pFHT_c(NumericVector t, double x0, double nu, double kappa, double sigma){
  int n = t.size();
  NumericVector ret(n);
  double k, smt, lambda_k, c_k, temp_k;
  
  for (int i = 0; i < n; ++i){
    // check if out of support
    if (t[i] <= 0) {
      ret[i] = 0;
      continue;
    }
    // reset for next element in vector t
    k = 1.0; smt = 0.0; lambda_k = 0.0; c_k = 0.0; temp_k = 1.0;

    while (std::abs(temp_k) > 0.0) {
      lambda_k = pow((2 * k - 1), 2) * pow(sigma, 2) * 
        pow(M_PI, 2) / (8 * pow((kappa - nu), 2));
      c_k = pow(-1, k + 1) * 4 / ((2 * k - 1) * M_PI) * 
        cos((2 * k - 1) * M_PI * (kappa - x0) / (2 * (kappa - nu)));
      temp_k = c_k * exp(- lambda_k * t[i]);
      smt = smt + temp_k;
      k++;
    }

    // dealing with numerical error 
    if (0 <= (1-smt) & (1-smt) <= 1) {
      ret[i] = 1-smt;
    } else if (1 < (1-smt)) {
      ret[i] = 1;
    } else {
      ret[i] = 0;
    }
  }
  return(ret);
}

// [[Rcpp::export]]
NumericVector solveQuad(double slope, double intercept, double t_start, double u){
  NumericVector ret(1);
  double a_root = slope / 2;
  double b_root = intercept;
  double c_root = -slope * pow(t_start, 2)/2 - intercept * t_start - u;
  ret[0] = ( -b_root + sqrt( pow(b_root, 2) - (4 * a_root * c_root) ) ) / (2*a_root);
  return(ret);
}

// [[Rcpp::export]]
List rFHT_c(int n, double x0, double nu, double kappa, double sigma, double k, double lambda_1, 
                NumericVector t_point,
                NumericVector slopes,
                NumericVector intercepts,
                NumericVector M,
                NumericVector I,
                NumericVector pt_point){
  NumericVector canddone(n);
  NumericVector U, u_whichpart, u, ratio, u_exp;
  NumericVector cand(1);
  LogicalVector u_logic(k);
  // count for calculate the rejection rate
  double count = 0;
  int part;
  for (int i = 0; i < n; i++){
    u_whichpart = Rcpp::runif(1, 0, 1);
    u_logic = u_whichpart[0] > pt_point;
    part = sum(u_logic);
    if (part < (k-1)){
      while(true){
        U = Rcpp::runif(1, 0, 1);
        u = Rcpp::runif(1, 0, I[part]);
        cand = solveQuad(slopes[part], intercepts[part], t_point[part], u[0]);
        count ++;
        ratio = dFHT_c(cand, x0, nu, kappa, sigma) / (M[part] * (slopes[part] * cand + intercepts[part]) );
        if (U[0] <= ratio[0]){
          canddone[i] = cand[0];
          break;
        }
      }
    }else{
      while(true){
        U = Rcpp::runif(1, 0, 1);
        u = Rcpp::runif(1, 1-I[part], 1);
        cand =  Rcpp::qexp(u, lambda_1);
        count ++;
        ratio = dFHT_c(cand, x0, nu, kappa, sigma)/(M[part] * Rcpp::dexp(cand, lambda_1, false) );
        if (U[0] <= ratio[0]){
          canddone[i] = cand[0];
          break;
        }
      }
    }
  }
  double rate = 1 - (n/count);
  // Rcout << "Count: " << count << "\n";
  List ret = List::create(_["samples"] = canddone,
                          _["rejectRate"] = rate);
  return(ret);
}
