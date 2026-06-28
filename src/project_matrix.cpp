#include <Rcpp.h>
#include <queue>
#include <vector>
using namespace Rcpp;

namespace {
inline int matrix_index(int row, int col, int size) {
  return row + col * size;
}
}

// [[Rcpp::export]]
NumericMatrix project_matrix_cpp_impl_(NumericMatrix S, IntegerVector sigma) {
  const int size = S.nrow();

  if (S.ncol() != size) {
    stop("`S` must be square.");
  }
  if (sigma.size() != size) {
    stop("`sigma` must have the same length as the size of `S`.");
  }

  std::vector<int> next(size);
  for (int i = 0; i < size; ++i) {
    const int image = sigma[i] - 1;
    if (image < 0 || image >= size) {
      stop("`sigma` must contain values between 1 and nrow(S).");
    }
    next[i] = image;
  }

  NumericMatrix projected(size, size);
  std::vector<unsigned char> visited(size * size, 0);

  for (int start = 0; start < size * size; ++start) {
    if (visited[start]) {
      continue;
    }

    std::vector<int> orbit;
    std::queue<int> pending;

    visited[start] = 1;
    orbit.push_back(start);
    pending.push(start);

    while (!pending.empty()) {
      const int current = pending.front();
      pending.pop();

      const int row = current % size;
      const int col = current / size;
      const int candidates[2] = {
        matrix_index(next[row], next[col], size),
        matrix_index(col, row, size)
      };

      for (int candidate_index = 0; candidate_index < 2; ++candidate_index) {
        const int candidate = candidates[candidate_index];
        if (!visited[candidate]) {
          visited[candidate] = 1;
          orbit.push_back(candidate);
          pending.push(candidate);
        }
      }
    }

    double sum = 0.0;
    for (std::size_t orbit_index = 0; orbit_index < orbit.size(); ++orbit_index) {
      const int index = orbit[orbit_index];
      const int row = index % size;
      const int col = index / size;
      sum += S(row, col);
    }
    const double mean = sum / orbit.size();

    for (std::size_t orbit_index = 0; orbit_index < orbit.size(); ++orbit_index) {
      const int index = orbit[orbit_index];
      const int row = index % size;
      const int col = index / size;
      projected(row, col) = mean;
    }
  }

  return projected;
}
