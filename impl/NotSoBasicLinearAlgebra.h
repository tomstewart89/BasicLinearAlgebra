#pragma once

namespace BLA
{
template <typename ParentType, typename Dtype>
void Swap(MatrixBase<ParentType, ParentType::Rows, ParentType::Cols, Dtype> &A,
          MatrixBase<ParentType, ParentType::Rows, ParentType::Cols, Dtype> &B)
{
    Dtype tmp;
    for (int i = 0; i < ParentType::Rows; i++)
    {
        for (int j = 0; j < ParentType::Cols; j++)
        {
            tmp = A(i, j);
            A(i, j) = B(i, j);
            B(i, j) = tmp;
        }
    }
}

template <typename ParentTypeA, typename ParentTypeB, int Cols>
Matrix<3, Cols, typename ParentTypeA::DType> CrossProduct(
    const MatrixBase<ParentTypeA, 3, Cols, typename ParentTypeA::DType> &matA,
    const MatrixBase<ParentTypeB, 3, Cols, typename ParentTypeA::DType> &matB)
{
    Matrix<3, Cols, typename ParentTypeA::DType> ret;
    for (int i = 0; i < Cols; ++i)
    {
        ret(0,i) = matA(1,i) * matB(2,i) - matB(1,i) * matA(2,i);
        ret(1,i) = matA(2,i) * matB(0,i) - matB(2,i) * matA(0,i);
        ret(2,i) = matA(0,i) * matB(1,i) - matB(0,i) * matA(1,i);
    }
    return ret;
}

template <typename ParentType>
struct LUDecomposition
{
    bool singular;
    typename ParentType::DType parity;
    PermutationMatrix<ParentType::Rows, typename ParentType::DType> P;
    LowerUnitriangularMatrix<ParentType> L;
    UpperTriangularMatrix<ParentType> U;

    LUDecomposition(MatrixBase<ParentType, ParentType::Rows, ParentType::Cols, typename ParentType::DType> &A)
        : L(static_cast<ParentType &>(A)), U(static_cast<ParentType &>(A))
    {
        static_assert(ParentType::Rows == ParentType::Cols, "Input matrix must be square");
    }
};

template <typename ParentType>
struct CholeskyDecomposition
{
    bool positive_definite = true;
    LowerTriangularMatrix<ParentType> L;

    CholeskyDecomposition(MatrixBase<ParentType, ParentType::Rows, ParentType::Cols, typename ParentType::DType> &A)
        : L(static_cast<ParentType &>(A))
    {
    }
};

template <typename ParentType, int Dim>
LUDecomposition<ParentType> LUDecompose(MatrixBase<ParentType, Dim, Dim, typename ParentType::DType> &A)
{
    LUDecomposition<ParentType> decomp(A);
    auto &idx = decomp.P.idx;
    decomp.parity = 1.0;

    for (int i = 0; i < Dim; ++i)
    {
        idx[i] = i;
    }

    // row_scale stores the implicit scaling of each row
    typename ParentType::DType row_scale[Dim];

    for (int i = 0; i < Dim; ++i)
    {
        // Loop over rows to get the implicit scaling information.
        typename ParentType::DType largest_elem = 0.0;

        for (int j = 0; j < Dim; ++j)
        {
            typename ParentType::DType this_elem = fabs(A(i, j));
            largest_elem = max(this_elem, largest_elem);
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
    for (int j = 0; j < Dim; ++j)
    {
        // Calculate beta ij
        for (int i = 0; i < j; ++i)
        {
            typename ParentType::DType sum = 0.0;

            for (int k = 0; k < i; ++k)
            {
                sum += A(i, k) * A(k, j);
            }

            A(i, j) -= sum;
        }

        // Calcuate alpha ij (before division by the pivot)
        for (int i = j; i < Dim; ++i)
        {
            typename ParentType::DType sum = 0.0;

            for (int k = 0; k < j; ++k)
            {
                sum += A(i, k) * A(k, j);
            }

            A(i, j) -= sum;
        }

        // Search for largest pivot element
        typename ParentType::DType largest_elem = 0.0;
        int argmax = j;

        for (int i = j; i < Dim; i++)
        {
            typename ParentType::DType this_elem = row_scale[i] * fabs(A(i, j));

            if (this_elem >= largest_elem)
            {
                largest_elem = this_elem;
                argmax = i;
            }
        }

        if (j != argmax)
        {
            auto row_argmax = A.Row(argmax);
            auto row_j = A.Row(j);
            Swap(row_argmax, row_j);

            decomp.parity = -decomp.parity;

            // swap indices
            {
                auto tmp = idx[j];
                idx[j] = idx[argmax];
                idx[argmax] = tmp;
            }

            row_scale[argmax] = row_scale[j];
        }

        if (A(j, j) == 0.0)
        {
            decomp.singular = true;
            return decomp;
        }

        if (j != Dim)
        {
            // Now, finally, divide by the pivot element.
            typename ParentType::DType pivot_inv = 1.0 / A(j, j);

            for (int i = j + 1; i < Dim; ++i)
            {
                A(i, j) *= pivot_inv;
            }
        }
    }

    decomp.singular = false;
    return decomp;
}

template <int Dim, class LUType, class BType>
Matrix<Dim, 1, typename BType::DType> LUSolve(const LUDecomposition<LUType> &decomp,
                                              const MatrixBase<BType, Dim, 1, typename BType::DType> &b)
{
    Matrix<Dim, 1, typename BType::DType> x, tmp;

    auto &idx = decomp.P.idx;
    auto &LU = decomp.L.parent;

    // Forward substitution to solve L * y = b
    for (int i = 0; i < Dim; ++i)
    {
        typename BType::DType sum = 0.0;

        for (int j = 0; j < i; ++j)
        {
            sum += LU(i, j) * tmp(idx[j]);
        }

        tmp(idx[i]) = b(idx[i]) - sum;
    }

    // Backward substitution to solve U * x = y
    for (int i = Dim - 1; i >= 0; --i)
    {
        typename BType::DType sum = 0.0;

        for (int j = i + 1; j < Dim; ++j)
        {
            sum += LU(i, j) * tmp(idx[j]);
        }

        tmp(idx[i]) = (tmp(idx[i]) - sum) / LU(i, i);
    }

    // Undo the permutation
    for (int i = 0; i < Dim; ++i)
    {
        x(i) = tmp(idx[i]);
    }

    return x;
}

template <typename ParentType, int Dim>
CholeskyDecomposition<ParentType> CholeskyDecompose(MatrixBase<ParentType, Dim, Dim, typename ParentType::DType> &A)
{
    CholeskyDecomposition<ParentType> chol(A);

    for (int i = 0; i < Dim; ++i)
    {
        for (int j = i; j < Dim; ++j)
        {
            float sum = A(i, j);

            for (int k = i - 1; k >= 0; --k)
            {
                sum -= A(i, k) * A(j, k);
            }

            if (i == j)
            {
                if (sum <= 0.0)
                {
                    chol.positive_definite = false;
                    return chol;
                }
                A(i, i) = sqrt(sum);
            }
            else
            {
                A(j, i) = sum / A(i, i);
            }
        }
    }

    return chol;
}

template <int Dim, class LUType, class BType>
Matrix<Dim, 1, typename BType::DType> CholeskySolve(const CholeskyDecomposition<LUType> &decomp,
                                                    const MatrixBase<BType, Dim, 1, typename BType::DType> &b)
{
    Matrix<Dim, 1, typename BType::DType> x;
    auto &A = decomp.L.parent;

    for (int i = 0; i < Dim; ++i)
    {
        float sum = b(i);

        for (int k = i - 1; k >= 0; --k)
        {
            sum -= A(i, k) * x(k);
        }

        x(i) = sum / A(i, i);
    }

    for (int i = Dim - 1; i >= 0; --i)
    {
        float sum = x(i);

        for (int k = i + 1; k < Dim; ++k)
        {
            sum -= A(k, i) * x(k);
        }

        x(i) = sum / A(i, i);
    }

    return x;
}

template <int Dim, typename InType, typename OutType, typename DType>
bool Invert(const MatrixBase<InType, Dim, Dim, DType> &A, MatrixBase<OutType, Dim, Dim, DType> &out)
{
    Matrix<Dim, Dim, DType> A_copy = A;

    auto decomp = LUDecompose(A_copy);

    if (decomp.singular)
    {
        return false;
    }

    Matrix<Dim, 1, DType> b = Zeros<Dim, 1, DType>();

    for (int j = 0; j < Dim; ++j)
    {
        b(j) = 1.0;
        out.Column(j) = LUSolve(decomp, b);
        b(j) = 0.0;
    }

    return true;
}

template <int Dim, class ParentType>
bool Invert(MatrixBase<ParentType, Dim, Dim, typename ParentType::DType> &A)
{
    return Invert(A, A);
}

template <int Dim, class ParentType>
Matrix<Dim, Dim, typename ParentType::DType> Inverse(
    const MatrixBase<ParentType, Dim, Dim, typename ParentType::DType> &A)
{
    Matrix<Dim, Dim, typename ParentType::DType> out;
    Invert(A, out);
    return out;
}

// LU-Decomposition only works for floating point numbers. Use Bareiss algorithm for (signed) integer types.
template <typename ParentType, typename Dtype, int Dim>
typename Types::enable_if<Types::is_floating_point<Dtype>::value, Dtype>::type
DeterminantLUDecomposition(const MatrixBase<ParentType, Dim, Dim, Dtype> &A)
{
    Matrix<Dim, Dim, Dtype> A_copy = A;

    auto decomp = LUDecompose(A_copy);

    Dtype det = decomp.parity;

    for (int i = 0; i < Dim; ++i)
    {
        det *= decomp.U(i, i);
    }

    return det;
}

// Bareiss algorithm works for all (signed) types, but for floating-point numbers LU-Decomposition is faster.
template <typename ParentType, typename Dtype, int Dim>
typename Types::enable_if<Types::is_signed<Dtype>::value, Dtype>::type
DeterminantBareissAlgorithm(const MatrixBase<ParentType, Dim, Dim, Dtype> &A)
{
    Matrix<Dim, Dim, Dtype> A_copy = A;

    int sign = 1;
    Dtype prev = 1;

    for (int i = 0; i < Dim; i++)
    {
        if (A_copy(i, i) == 0)
        {
            int j = i + 1;
            for (; j < Dim; j++)
            {
                if (A_copy(j, i) != 0) break;
            }
            if (j == Dim) return 0;
            auto row_i = A_copy.Row(i);
            auto row_j = A_copy.Row(j);
            Swap(row_i, row_j);
            sign = - sign;
        }
        for (int j = i + 1; j < Dim; j++)
        {
            for (int k = i + 1; k < Dim; k++)
            {
                A_copy(j, k) = (A_copy(j, k) * A_copy(i, i) - A_copy(j, i) * A_copy(i, k)) / prev;
            }
        }
        prev = A_copy(i, i);
    }
    return sign * A_copy(Dim - 1, Dim - 1);
}

template <typename ParentType, typename Dtype, int Dim>
typename Types::enable_if<Types::is_floating_point<Dtype>::value, Dtype>::type
Determinant(const MatrixBase<ParentType, Dim, Dim, Dtype> &A) { return DeterminantLUDecomposition(A); }

template <typename ParentType, typename Dtype, int Dim>
typename Types::enable_if<Types::is_signed_integer<Dtype>::value, Dtype>::type
Determinant(const MatrixBase<ParentType, Dim, Dim, Dtype> &A) { return DeterminantBareissAlgorithm(A); }

template <typename DerivedType>
typename DerivedType::DType Norm(const DownCast<DerivedType> &A)
{
    typename DerivedType::DType sum_sq = 0.0;

    for (int i = 0; i < DerivedType::Rows; ++i)
    {
        for (int j = 0; j < DerivedType::Cols; ++j)
        {
            sum_sq += A(i, j) * A(i, j);
        }
    }
    return sqrt(sum_sq);
}

template <class DerivedType>
typename DerivedType::DType Trace(const DownCast<DerivedType> &A)
{
    typename DerivedType::DType sum_diag = 0.0;

    for (int i = 0; i < DerivedType::Rows; ++i)
    {
        sum_diag += A(i, i);
    }
    return sum_diag;
}

template <int Inputs, int Outputs, typename DType>
struct MatrixFunctor
{
    virtual Matrix<Outputs, 1, DType> operator()(const Matrix<Inputs, 1, DType> &x) const = 0;
};

template <int Inputs, int Outputs, typename InType>
Matrix<Outputs, Inputs, typename InType::DType> Jacobian(
    const MatrixFunctor<Inputs, Outputs, typename InType::DType> &f,
    const MatrixBase<InType, Inputs, 1, typename InType::DType> &x, const typename InType::DType h = 1e-4)
{
    using DType = typename InType::DType;

    Matrix<Outputs, Inputs, DType> jacobian;
    Matrix<Outputs> f_x = f(x);
    Matrix<Inputs, 1, DType> h_vec = Zeros<Inputs, 1, DType>();

    for (int i = 0; i < Inputs; i++)
    {
        h_vec(i) = h;
        Matrix<Outputs, 1, DType> f_xh = f(x + h_vec);
        jacobian.Column(i) = (f_xh - f_x) / h;
        h_vec(i) = 0;
    }

    return jacobian;
}

}  // namespace BLA
