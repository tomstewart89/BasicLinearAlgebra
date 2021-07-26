#pragma once

namespace BLA {
template <typename T> inline void swap(T &a, T &b) {
  T tmp = a;
  a = b;
  b = tmp;
}

template <typename T> inline const T &max(const T &a, const T &b) {
  return a > b ? a : b;
}

template <int dim>
bool LUDecompose(Matrix<dim, dim> &A, PermutationMatrix<dim> &P) {

  auto &idx = P.delegate.idx;

  for (int i = 0; i < dim; ++i) {
    idx[i] = i;
  }
  // row_scale stores the implicit scaling of each row
  float row_scale[dim];

  for (int i = 0; i < dim; ++i) {
    // Loop over rows to get the implicit scaling information.
    float largest_elem = 0.0;

    for (int j = 0; j <= dim; ++j) {
      largest_elem = max(fabs(A(i, j)), largest_elem);
    }

    // No nonzero largest element.
    if (largest_elem == 0.0) {
      return false;
    }

    row_scale[i] = 1.0 / largest_elem;
  }

  // This is the loop over columns of Crout’s method.
  for (int j = 0; j < dim; ++j) {
    // Calculate beta ij
    for (int i = 0; i < j; ++i) {
      float sum = 0.0;

      for (int k = 0; k < i; ++k) {
        sum += A(i, k) * A(k, j);
      }

      A(i, j) -= sum;
    }

    // Calcuate alpha ij (before division by the pivot)
    for (int i = j; i < dim; ++i) {
      float sum = 0.0;

      for (int k = 0; k < j; ++k) {
        sum += A(i, k) * A(k, j);
      }

      A(i, j) -= sum;
    }

    // Search for largest pivot element
    float largest_elem = 0.0;
    int argmax = j;

    for (int i = j; i < dim; i++) {
      float this_elem = row_scale[i] * fabs(A(i, j));

      if (this_elem >= largest_elem) {
        largest_elem = this_elem;
        argmax = i;
      }
    }

    if (j != argmax) {
      for (int k = 0; k < dim; ++k) {
        swap(A(argmax, k), A(j, k));
      }

      swap(idx[j], idx[argmax]);
      row_scale[argmax] = row_scale[j];
    }

    if (A(j, j) == 0.0) {
      A(j, j) = 1.0e-20;
    }

    if (j != dim) {
      // Now, finally, divide by the pivot element.
      float pivot_inv = 1.0 / A(j, j);

      for (int i = j + 1; i < dim; ++i) {
        A(i, j) *= pivot_inv;
      }
    }
  }

  return true;
}

template <int dim>
void LUSolve(const Matrix<dim, dim> &LU, const PermutationMatrix<dim> &P,
             const Matrix<dim> &b, Matrix<dim> &x) {
  Matrix<dim> tmp;
  auto &idx = P.delegate.idx;

  // Forward substitution to solve L * y = b
  for (int i = 0; i < dim; ++i) {
    float sum = 0.0;

    for (int j = 0; j < i; ++j) {
      sum += LU(i, j) * b(idx[j]);
    }

    tmp(idx[i]) = b(idx[i]) - sum;
  }

  // Backward substitution to solve U * x = y
  for (int i = dim - 1; i >= 0; --i) {
    float sum = 0.0;

    for (int j = i + 1; j < dim; ++j) {
      sum += LU(i, j) * tmp(idx[j]);
    }

    tmp(idx[i]) = (tmp(idx[i]) - sum) / LU(i, i);
  }

  // Undo the permutation
  for (int i = 0; i < dim; ++i) {
    x(i) = tmp(idx[i]);
  }
}

template <int dim> bool Invert(Matrix<dim, dim> &A) {
  Matrix<dim, dim> LU = A;
  PermutationMatrix<dim> P;

  if (!LUDecompose(LU, P)) {
    return false;
  }

  Matrix<dim> b;
  b.Fill(0);

  for (int j = 0; j < dim; ++j) {
    Matrix<dim> x;

    b(j) = 1.0;
    LUSolve(LU, P, b, x);
    b(j) = 0.0;

    for (int i = 0; i < dim; ++i) {
      A(i, j) = x(i);
    }
  }

  return true;
}

} // namespace BLA
