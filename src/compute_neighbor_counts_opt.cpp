#include <Rcpp.h>
#include <algorithm> // for std::lower_bound
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_neighbor_counts_opt(NumericVector x, NumericVector y, NumericVector radii) {
  int n = x.size();
  int nr = radii.size();
  NumericMatrix counts(n, nr);
  
  // Precompute squared radii.
  NumericVector radii2(nr);
  for(int k = 0; k < nr; k++){
    radii2[k] = radii[k] * radii[k];
  }
  
  // Loop over all pairs of points.
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i == j) continue; // skip self-comparison
      double dx = x[j] - x[i];
      double dy = y[j] - y[i];
      double d2 = dx * dx + dy * dy;
      
      // Use binary search (std::lower_bound) on the sorted radii2 vector.
      // This finds the first index where d2 <= radii2[index].
      NumericVector::iterator lb = std::lower_bound(radii2.begin(), radii2.end(), d2);
      int index = lb - radii2.begin();
      if(index < nr) {
        // For all indices from 'index' to end, the neighbor qualifies.
        for (int k = index; k < nr; k++){
          counts(i, k) += 1;
        }
      }
    }
  }
  return counts;
}