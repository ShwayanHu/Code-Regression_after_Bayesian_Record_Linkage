#include <Rcpp.h>
using namespace Rcpp;

// Helper to convert field number to field name
std::string num2field(int num) {
  switch (num) {
  case 1: return "f1";
  case 2: return "f2";
  case 3: return "f3";
  case 4: return "f4";
  default: return "";
  }
}

// Sample from Dirichlet distribution
NumericVector rdirichlet_cpp(NumericVector alpha) {
  int K = alpha.size();
  NumericVector x(K);
  double sum_x = 0.0;
  
  for (int k = 0; k < K; ++k) {
    x[k] = R::rgamma(alpha[k], 1.0);
    sum_x += x[k];
  }
  
  for (int k = 0; k < K; ++k) {
    x[k] /= sum_x;
  }
  
  return x;
}

// [[Rcpp::export]]
List update_u_cpp(List u, List Z, DataFrame gamma, IntegerVector L, List beta) {
  IntegerVector record_1 = gamma["record_1"];
  IntegerVector record_2 = gamma["record_2"];
  int n_gamma = gamma.nrows();
  
  for (int f = 1; f <= 4; ++f) {
    std::string field = num2field(f);
    IntegerVector gamma_f = gamma[field];
    int max_level = L[f - 1];
    
    NumericVector count(max_level + 1);
    
    for (int l = 0; l <= max_level; ++l) {
      for (int i = 0; i < n_gamma; ++i) {
        if (gamma_f[i] == l) {
          int idx_i = record_1[i];
          int idx_j = record_2[i];
          
          IntegerVector Z_j = Z[idx_j - 1];  // adjust for 0-based indexing
          bool found = false;
          
          for (int k = 0; k < Z_j.size(); ++k) {
            if (Z_j[k] == idx_i) {
              found = true;
              break;
            }
          }
          
          if (!found) {
            count[l] += 1;
          }
        }
      }
    }
    
    NumericVector beta_f = beta[field];
    NumericVector dirichlet_param = count + beta_f;
    
    u[field] = rdirichlet_cpp(dirichlet_param);
  }
  
  return u;
}