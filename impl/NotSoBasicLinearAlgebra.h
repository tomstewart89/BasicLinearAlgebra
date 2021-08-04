#pragma once

namespace BLA
{
template <typename T>
inline void _swap(T &a, T &b)
{
    T tmp = a;
    a = b;
    b = tmp;
}

template <typename T>
inline const T &_max(const T &a, const T &b)
{
    return a > b ? a : b;
}

template <int dim, class MemT>
struct LUDecomposition
{
    bool singular;
    Permutation<dim, typename MemT::elem_t> permutation;
    LowerTriangleOnesDiagonal<MemT> lower;
    UpperTriangle<MemT> upper;

    LUDecomposition(Matrix<dim, dim, MemT> &A) : lower(A.storage), upper(A.storage) {}

    PermutationMatrix<dim, typename MemT::elem_t> P()
    {
        return PermutationMatrix<dim, typename MemT::elem_t>(permutation);
    }
    LowerTriangularDiagonalOnesMatrix<dim, dim, MemT> L()
    {
        return LowerTriangularDiagonalOnesMatrix<dim, dim, MemT>(lower);
    }
    UpperTriangularMatrix<dim, dim, MemT> U() { return UpperTriangularMatrix<dim, dim, MemT>(upper); }
};

template <int dim, class MemT>
LUDecomposition<dim, MemT> LUDecompose(Matrix<dim, dim, MemT> &A)
{
    LUDecomposition<dim, MemT> decomp(A);
    auto &idx = decomp.permutation.idx;

    for (int i = 0; i < dim; ++i)
    {
        idx[i] = i;
    }

    // row_scale stores the implicit scaling of each row
    typename MemT::elem_t row_scale[dim];

    for (int i = 0; i < dim; ++i)
    {
        // Loop over rows to get the implicit scaling information.
        typename MemT::elem_t largest_elem = 0.0;

        for (int j = 0; j <= dim; ++j)
        {
            typename MemT::elem_t this_elem = fabs(A(i, j));
            largest_elem = _max(this_elem, largest_elem);
        }

        // No nonzero largest element.
        if (largest_elem == 0.0)
        {
            decomp.singular = true;
            return decomp;
        }

        row_scale[i] = 1.0 / largest_elem;
    }

    // This is the loop over columns of Crout’s method.
    for (int j = 0; j < dim; ++j)
    {
        // Calculate beta ij
        for (int i = 0; i < j; ++i)
        {
            typename MemT::elem_t sum = 0.0;

            for (int k = 0; k < i; ++k)
            {
                sum += A(i, k) * A(k, j);
            }

            A(i, j) -= sum;
        }

        // Calcuate alpha ij (before division by the pivot)
        for (int i = j; i < dim; ++i)
        {
            typename MemT::elem_t sum = 0.0;

            for (int k = 0; k < j; ++k)
            {
                sum += A(i, k) * A(k, j);
            }

            A(i, j) -= sum;
        }

        // Search for largest pivot element
        typename MemT::elem_t largest_elem = 0.0;
        int argmax = j;

        for (int i = j; i < dim; i++)
        {
            typename MemT::elem_t this_elem = row_scale[i] * fabs(A(i, j));

            if (this_elem >= largest_elem)
            {
                largest_elem = this_elem;
                argmax = i;
            }
        }

        if (j != argmax)
        {
            for (int k = 0; k < dim; ++k)
            {
                _swap(A(argmax, k), A(j, k));
            }

            _swap(idx[j], idx[argmax]);
            row_scale[argmax] = row_scale[j];
        }

        if (A(j, j) == 0.0)
        {
            A(j, j) = 1.0e-20;
        }

        if (j != dim)
        {
            // Now, finally, divide by the pivot element.
            typename MemT::elem_t pivot_inv = 1.0 / A(j, j);

            for (int i = j + 1; i < dim; ++i)
            {
                A(i, j) *= pivot_inv;
            }
        }
    }

    decomp.singular = false;
    return decomp;
}

template <int dim, class MemT1, class MemT2>
Matrix<dim> LUSolve(const LUDecomposition<dim, MemT1> &decomp, const Matrix<dim, 1, MemT2> &b)
{
    Matrix<dim, 1, MemT2> x, tmp;

    auto &idx = decomp.permutation.idx;
    auto &LU = decomp.lower.parent;

    // Forward substitution to solve L * y = b
    for (int i = 0; i < dim; ++i)
    {
        typename MemT2::elem_t sum = 0.0;

        for (int j = 0; j < i; ++j)
        {
            sum += LU(i, j) * tmp(idx[j]);
        }

        tmp(idx[i]) = b(idx[i]) - sum;
    }

    // Backward substitution to solve U * x = y
    for (int i = dim - 1; i >= 0; --i)
    {
        typename MemT2::elem_t sum = 0.0;

        for (int j = i + 1; j < dim; ++j)
        {
            sum += LU(i, j) * tmp(idx[j]);
        }

        tmp(idx[i]) = (tmp(idx[i]) - sum) / LU(i, i);
    }

    // Undo the permutation
    for (int i = 0; i < dim; ++i)
    {
        x(i) = tmp(idx[i]);
    }

    return x;
}

template <int dim, class MemT>
bool Invert(Matrix<dim, dim, MemT> &A)
{
    Matrix<dim, dim> A_copy = A;

    auto decomp = LUDecompose(A_copy);

    if (decomp.singular)
    {
        return false;
    }

    Matrix<dim> b = Zeros<dim>();

    for (int j = 0; j < dim; ++j)
    {
        b(j) = 1.0;
        A.Column(j) = LUSolve(decomp, b);
        b(j) = 0.0;
    }

    return true;
}

template <int dim, class MemT>
typename MemT::elem_t Determinant(const Matrix<dim, dim, MemT> &A)
{
    typename MemT::elem_t det = 0;

    // Add the determinants of all the minors
    for (int i = 0; i < dim; i++)
    {
        Minor<MemT> del(A.storage, i, 0);
        Matrix<dim - 1, dim - 1, Minor<MemT>> m(del);

        if (i % 2)
        {
            det -= Determinant(m) * A(i, 0);
        }
        else
        {
            det += Determinant(m) * A(i, 0);
        }
    }

    return det;
}

template <class MemT>
typename MemT::elem_t Determinant(Matrix<2, 2, MemT> &A)
{
    return A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);
}

}  // namespace BLA
