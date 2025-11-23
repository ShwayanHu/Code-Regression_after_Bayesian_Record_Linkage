#include <Rcpp.h>
#include <unordered_set>
#include <random>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector update_Z_Sadinle_fast_cpp(
    IntegerVector Z_tmp,
    int n1,
    int n2,
    DataFrame gamma,
    List m_tmp,
    List u_tmp,
    double beta_pi,
    double alpha_pi
) {
  int nrow_gamma = gamma.nrows();
  
  // Step 1: Build gamma_array (3D array)
  std::vector<NumericVector> f_list = {
    gamma["f1"], gamma["f2"], gamma["f3"], gamma["f4"]
  };
  IntegerVector record_1 = gamma["record_1"];
  IntegerVector record_2 = gamma["record_2"];
  
  std::vector<std::vector<std::vector<int>>> gamma_array(
      n1, std::vector<std::vector<int>>(n2, std::vector<int>(4, -1))
  );
  
  for (int i = 0; i < nrow_gamma; ++i) {
    int r1 = record_1[i] - 1;
    int r2 = record_2[i] - 1;
    for (int f = 0; f < 4; ++f) {
      gamma_array[r1][r2][f] = f_list[f][i];
    }
  }
  
  // Random engine
  std::random_device rd;
  std::mt19937 gen(rd());
  
  for (int j = 0; j < n2; ++j) {
    // Step 2: Get Z_minus_j and n_12
    std::unordered_set<int> linked_set;
    for (int jj = 0; jj < n2; ++jj) {
      if (jj == j) continue;
      int val = Z_tmp[jj];
      if (val <= n1)
        linked_set.insert(val);
    }
    int n12 = linked_set.size();
    
    // Step 3: Compute weights
    std::vector<int> possible_qs(n1 + 1);
    std::iota(possible_qs.begin(), possible_qs.end() - 1, 1);
    possible_qs[n1] = n1 + j + 1;  // R is 1-indexed
    
    std::vector<double> w_j(possible_qs.size());
    
    for (int idx = 0; idx < possible_qs.size(); ++idx) {
      int q = possible_qs[idx];
      
      if (q <= n1) {
        if (linked_set.find(q) != linked_set.end()) {
          w_j[idx] = 0.0;
        } else {
          double log_ratio = 0.0;
          for (int f = 0; f < 4; ++f) {
            int g_val = gamma_array[q - 1][j][f];
            if (g_val >= 0) {
              NumericVector m_vec = m_tmp[f];
              NumericVector u_vec = u_tmp[f];
              int l = g_val;  // Already zero-indexed
              log_ratio += std::log(m_vec[l] / u_vec[l]);
            }
          }
          w_j[idx] = std::exp(log_ratio);
        }
      } else {
        double prior_term = (n1 - n12) * (n2 - n12 - 1 + beta_pi) / (n12 + alpha_pi);
        w_j[idx] = prior_term;
      }
    }
    
    // Normalize w_j
    double total_w = std::accumulate(w_j.begin(), w_j.end(), 0.0);
    for (double& w : w_j) w /= total_w;
    
    // Sample from discrete distribution
    std::discrete_distribution<> d(w_j.begin(), w_j.end());
    int sampled_idx = d(gen);
    Z_tmp[j] = possible_qs[sampled_idx];
  }
  
  return Z_tmp;
}