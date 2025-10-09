#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_neighbor_counts_grid(NumericVector x,
                                           NumericVector y,
                                           NumericVector radii) {
  const int n  = x.size();
  if (y.size() != n) stop("x and y must have same length");
  const int nr = radii.size();
  if (nr == 0) stop("radii must be non-empty");
  
  // Precompute squared radii, require sorted ascending (to match original semantics)
  std::vector<double> r2(nr);
  for (int k = 0; k < nr; ++k) r2[k] = radii[k] * radii[k];
  if (!std::is_sorted(r2.begin(), r2.end()))
    stop("radii must be sorted ascending (same requirement as the original function).");
  
  const double rmax2 = r2.back();
  const double cell  = std::sqrt(rmax2); // grid cell size = r_max
  
  // Bounding box
  double xmin = x[0], xmax = x[0], ymin = y[0], ymax = y[0];
  for (int i = 1; i < n; ++i) {
    double xi = x[i], yi = y[i];
    if (xi < xmin) xmin = xi;
    if (xi > xmax) xmax = xi;
    if (yi < ymin) ymin = yi;
    if (yi > ymax) ymax = yi;
  }
  
  // Grid dims (ensure at least 1)
  const int gx = std::max(1, (int)std::floor((xmax - xmin) / cell) + 1);
  const int gy = std::max(1, (int)std::floor((ymax - ymin) / cell) + 1);
  
  // Map points into grid buckets
  std::vector< std::vector<int> > grid((size_t)gx * (size_t)gy);
  grid.shrink_to_fit();
  for (int i = 0; i < n; ++i) {
    int cx = std::min(gx - 1, std::max(0, (int)std::floor((x[i] - xmin) / cell)));
    int cy = std::min(gy - 1, std::max(0, (int)std::floor((y[i] - ymin) / cell)));
    grid[(size_t)cx + (size_t)cy * (size_t)gx].push_back(i);
  }
  
  // Output (numeric to match your original)
  NumericMatrix counts(n, nr);
  
  // Work buffers reused per anchor
  std::vector<int> hist(nr);
  
  // For each anchor point, only scan neighbors in its 3x3 neighborhood
  for (int i = 0; i < n; ++i) {
    std::fill(hist.begin(), hist.end(), 0);
    
    int acx = std::min(gx - 1, std::max(0, (int)std::floor((x[i] - xmin) / cell)));
    int acy = std::min(gy - 1, std::max(0, (int)std::floor((y[i] - ymin) / cell)));
    
    for (int dcy = -1; dcy <= 1; ++dcy) {
      int cy = acy + dcy;
      if (cy < 0 || cy >= gy) continue;
      for (int dcx = -1; dcx <= 1; ++dcx) {
        int cx = acx + dcx;
        if (cx < 0 || cx >= gx) continue;
        
        const std::vector<int>& bucket = grid[(size_t)cx + (size_t)cy * (size_t)gx];
        for (int j : bucket) {
          if (j == i) continue; // skip self
          double dx = (double)x[j] - (double)x[i];
          double dy = (double)y[j] - (double)y[i];
          double d2 = dx * dx + dy * dy;
          if (d2 > rmax2) continue;
          
          auto lb = std::lower_bound(r2.begin(), r2.end(), d2);
          int k = (int)std::distance(r2.begin(), lb);
          if (k < nr) hist[k] += 1; // first qualifying radius bin
        }
      }
    }
    
    // Prefix sum to convert first-hit histogram into cumulative counts over radii
    int c = 0;
    for (int k = 0; k < nr; ++k) {
      c += hist[k];
      counts(i, k) = (double)c;
    }
  }
  
  return counts;
}
