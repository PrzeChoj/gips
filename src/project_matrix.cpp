#include <Rcpp.h>
#include <algorithm>
#include <queue>
#include <vector>
using namespace Rcpp;

namespace {
inline int matrix_index(int row, int col, int size) {
  return row + col * size;
}
}

// [[Rcpp::export]]
List project_matrices_cpp_impl_(List matrices, IntegerVector sigma) {
  const int matrix_count = matrices.size();
  if (matrix_count < 1) {
    stop("`matrices` must contain at least one matrix.");
  }

  std::vector<NumericMatrix> inputs;
  std::vector<NumericMatrix> projected;
  inputs.reserve(matrix_count);
  projected.reserve(matrix_count);

  NumericMatrix first_matrix(matrices[0]);
  const int size = first_matrix.nrow();

  if (sigma.size() != size) {
    stop("`sigma` must have the same length as the size of `S`.");
  }
  for (int matrix_index = 0; matrix_index < matrix_count; ++matrix_index) {
    NumericMatrix current_matrix(matrices[matrix_index]);
    if (current_matrix.nrow() != size || current_matrix.ncol() != size) {
      stop("Every matrix in `matrices` must be square and match the length of `sigma`.");
    }
    inputs.push_back(current_matrix);
    projected.push_back(NumericMatrix(size, size));
  }

  std::vector<int> next(size);
  for (int i = 0; i < size; ++i) {
    const int image = sigma[i] - 1;
    if (image < 0 || image >= size) {
      stop("`sigma` must contain values between 1 and nrow(S).");
    }
    next[i] = image;
  }

  std::vector<unsigned char> visited(size * size, 0);
  std::vector<double> sums(matrix_count, 0.0);

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

    std::fill(sums.begin(), sums.end(), 0.0);
    for (std::size_t orbit_index = 0; orbit_index < orbit.size(); ++orbit_index) {
      const int index = orbit[orbit_index];
      const int row = index % size;
      const int col = index / size;
      for (int matrix_index = 0; matrix_index < matrix_count; ++matrix_index) {
        sums[matrix_index] += inputs[matrix_index](row, col);
      }
    }

    for (std::size_t orbit_index = 0; orbit_index < orbit.size(); ++orbit_index) {
      const int index = orbit[orbit_index];
      const int row = index % size;
      const int col = index / size;
      for (int matrix_index = 0; matrix_index < matrix_count; ++matrix_index) {
        projected[matrix_index](row, col) = sums[matrix_index] / orbit.size();
      }
    }
  }

  List out(matrix_count);
  for (int matrix_index = 0; matrix_index < matrix_count; ++matrix_index) {
    out[matrix_index] = projected[matrix_index];
  }

  return out;
}
