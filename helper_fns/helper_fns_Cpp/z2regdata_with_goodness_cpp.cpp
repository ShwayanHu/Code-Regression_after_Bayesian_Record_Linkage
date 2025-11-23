#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List z2regdata_with_goodness_cpp(IntegerVector Z, 
                                 DataFrame file_1, 
                                 DataFrame file_2,
                                 List m, 
                                 List u,
                                 DataFrame gamma) {
  
  int n1 = file_1.nrows();
  int n2 = file_2.nrows();
  int z_length = Z.size();
  
  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> G;
  
  NumericVector file1_x = file_1["x"];
  NumericVector file2_y = file_2["y"];
  
  IntegerVector gamma_record1 = gamma["record_1"];
  IntegerVector gamma_record2 = gamma["record_2"];
  
  // Pre-extract m and u vectors for each field
  NumericVector m1 = m["field1"];
  NumericVector m2 = m["field2"];
  NumericVector m3 = m["field3"];
  NumericVector m4 = m["field4"];
  
  NumericVector u1 = u["field1"];
  NumericVector u2 = u["field2"];
  NumericVector u3 = u["field3"];
  NumericVector u4 = u["field4"];
  
  for(int j = 0; j < z_length; j++) {
    if(Z[j] == n1 + j + 1) continue; // +1 for R to C++ index conversion
    
    int slot_1 = Z[j] - 1; // Convert to 0-based index
    
    X.push_back(file1_x[slot_1]);
    Y.push_back(file2_y[j]);
    
    double log_lik_ratio = 0.0;
    
    // Find matching rows in gamma
    for(int k = 0; k < gamma.nrows(); k++) {
      if(gamma_record1[k] == slot_1 + 1 && gamma_record2[k] == j + 1) {
        // Add 1 to convert to R indices when accessing gamma columns
        int l1 = gamma[k]; // Assuming these are the field values
        int l2 = gamma[k]; // You'll need to adjust these to match your actual gamma structure
        int l3 = gamma[k];
        int l4 = gamma[k];
        
        log_lik_ratio += log(m1[l1]) - log(u1[l1]);
        log_lik_ratio += log(m2[l2]) - log(u2[l2]);
        log_lik_ratio += log(m3[l3]) - log(u3[l3]);
        log_lik_ratio += log(m4[l4]) - log(u4[l4]);
        
        break; // Assuming one match per record pair
      }
    }
    
    G.push_back(log_lik_ratio);
  }
  
  if(X.size() == 0) return R_NilValue;
  
  // Prepare results
  int n = X.size();
  NumericMatrix X_mat(n, 2);
  NumericVector Y_vec(n);
  NumericVector G_vec(n);
  
  for(int i = 0; i < n; i++) {
    X_mat(i, 0) = 1.0;
    X_mat(i, 1) = X[i];
    Y_vec[i] = Y[i];
    G_vec[i] = G[i];
  }
  
  return List::create(
    Named("X") = X_mat,
    Named("Y") = Y_vec,
    Named("G") = G_vec
  );
}