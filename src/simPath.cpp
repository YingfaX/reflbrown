#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List simPath_c(double x0, double nu,
               double kappa, double sigma, NumericVector dt){
  int n = dt.size();
  int i = 1;
  int count_idx = 1;
  double path_temp, Z;
  NumericVector path(n+1);
  NumericVector proc_idx(n+1);
  path[0] = x0;
  proc_idx[0] = count_idx;

  while(i <= n){
    Z = R::rnorm(0, 1);
    path_temp = path[i-1] + sigma * sqrt(dt[i-1]) * Z;
    if (path_temp > kappa){
      // reflection
      path[i] = 2 * kappa - path_temp;
      proc_idx[i] = count_idx;
    } else if (path_temp > nu) {
      path[i] = path_temp;
      proc_idx[i] = count_idx;
    } else {
      // hit boundary
      path[i] = nu;
      proc_idx[i] = count_idx;
      // reset process
      i++;
      count_idx++;
      proc_idx[i] = count_idx;
      path[i] = x0;
    }
    i++;
  }
  List ret = List::create(_["Xt"] = path,
                          _["proc_idx"] = proc_idx);
  return(ret);
}
