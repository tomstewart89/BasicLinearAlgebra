#pragma once

namespace BLA
{
template <typename T>
inline void bla_swap(T &a, T &b)
{
    T tmp = a;
    a = b;
    b = tmp;
}

template <int Dim, class ParentType>
struct LUDecomposition
{
    bool singular;
    typename ParentType::DType parity;
    PermutationMatrix<Dim, typename ParentType::DType> permutation;
    LowerTriangularDiagonalOnesMatrix<ParentType> lower;
    UpperTriangularMatrix<ParentType> upper;

    LUDecomposition(MatrixBase<ParentType, ParentType::Rows, ParentType::Cols, typename ParentType::DType> &A)
        : lower(A), upper(A)
    {
        static_assert(ParentType::Rows == ParentType::Cols);
    }

    PermutationMatrix<ParentType::Rows, typename ParentType::DType> P() { return permutation; }
    LowerTriangularDiagonalOnesMatrix<ParentType> L() { return lower; }
    UpperTriangularMatrix<ParentType> U() { return upper; }
};

template <typename ParentType, int Dim>
LUDecomposition<Dim, ParentType> LUDecompose(MatrixBase<ParentType, Dim, Dim, typename ParentType::DType> &A)
{
    LUDecomposition<Dim, ParentType> decomp(A);
    auto &idx = decomp.permutation.idx;
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
            for (int k = 0; k < Dim; ++k)
            {
                bla_swap(A(argmax, k), A(j, k));
            }

            decomp.parity = -decomp.parity;

            bla_swap(idx[j], idx[argmax]);
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
Matrix<Dim, 1, typename BType::DType> LUSolve(const LUDecomposition<Dim, LUType> &decomp,
                                              const MatrixBase<BType, Dim, 1, typename BType::DType> &b)
{
    Matrix<Dim, 1, typename BType::DType> x, tmp;

    auto &idx = decomp.permutation.idx;
    auto &LU = decomp.lower.parent;

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

template <typename ParentType, int Dim>
typename ParentType::DType Determinant(const MatrixBase<ParentType, Dim, Dim, typename ParentType::DType> &A)
{
    Matrix<Dim, Dim, typename ParentType::DType> A_copy = A;

    auto decomp = LUDecompose(A_copy);

    typename ParentType::DType det = decomp.parity;

    for (int i = 0; i < Dim; ++i)
    {
        det *= decomp.upper(i, i);
    }

    return det;
}

template <typename ParentType>
typename ParentType::DType Norm(
    const MatrixBase<ParentType, ParentType::Rows, ParentType::Cols, typename ParentType::DType> &A)
{
    typename ParentType::DType sum_sq = 0.0;

    for (int i = 0; i < ParentType::Rows; ++i)
    {
        for (int j = 0; j < ParentType::Cols; ++j)
        {
            sum_sq += A(i, j) * A(i, j);
        }
    }
    return sqrt(sum_sq);
}

template <class ParentType>
typename ParentType::DType Trace(
    const MatrixBase<ParentType, ParentType::Rows, ParentType::Cols, typename ParentType::DType> &A)
{
    typename ParentType::DType sum_diag = 0.0;

    for (int i = 0; i < ParentType::Rows; ++i)
    {
        sum_diag += A(i, i);
    }
    return sum_diag;
}


template<int rows, int cols, class MemT = Array<rows, cols, float>> 
class VectorValuedFunction{
    private:
        Matrix<rows, 1, Array<rows, 1, typename MemT::elem_t>> (*_f)(Matrix<cols, 1, Array<cols,1, typename MemT::elem_t>> x);

    public:
        VectorValuedFunction(){}

        VectorValuedFunction(Matrix<rows, 1, Array<rows, 1, typename MemT::elem_t>> (&f)(Matrix<cols, 1, Array<cols,1, typename MemT::elem_t>> x)){_f = f;}
        VectorValuedFunction(const VectorValuedFunction &V){_f = V._f;}
        
        ~VectorValuedFunction(){}

        virtual Matrix<rows, 1, Array<rows, 1, typename MemT::elem_t>> vv_f(Matrix<cols, 1, Array<cols,1, typename MemT::elem_t>> x){
            if(_f == nullptr){
                return Matrix<rows, 1, Array<rows, 1, typename MemT::elem_t>>();
            }else{
                return _f(x);
            }
        }

        Matrix<rows, 1, Array<rows, 1, typename MemT::elem_t>> operator ()(Matrix<cols, 1, Array<cols,1, typename MemT::elem_t>> x){
            return vv_f(x);
        }
};

template<int rows, int cols, class MemT = Array<rows, cols, float>>
using VVF = VectorValuedFunction<rows, cols, MemT>;


/**
 * @brief numerically approximate the jacobian of a vector valued function. Cannot be used at the boundary of the domain
 * 
 * @tparam n inputs
 * @tparam m outputs
 * @tparam MemT 
 * @param f vector valued functor 
 * @param _x Point to take a derivative at
 * @param h step size. Set at 0.0001, but for more precise applications, it can be lowered
 * @return Matrix<n, m, MemT> 
 */
template <int n, int m, class MemT = Array<n, m, float> >
Matrix<n, m, MemT> Jacobian(VectorValuedFunction<n,m,MemT> f, Matrix<m, 1, Array<m,1, typename MemT::elem_t>> _x,typename MemT::elem_t h = 0.0001)
{  
    Matrix<n, 1, Array<n, 1, typename MemT::elem_t>> f_x = f.vv_f(_x);

    Matrix<n, m, MemT> Jacob;
  
    for(int i = 0; i < m; i++){
        Matrix<m, 1, Array<m, 1, typename MemT::elem_t>> h_vec = Zeros<m>();
        h_vec(i) = 1.00;
        h_vec = h_vec * h;

        Matrix<n, 1, Array<n, 1, typename MemT::elem_t>> f_xh = f.vv_f(_x + h_vec);

        Matrix<n, 1, Array<n, 1, typename MemT::elem_t>> df = (f_xh - f_x)/h;

        
        for(int j = 0; j < n; j++){
            Jacob(j,i) = df(j);
        }
    }
    return Jacob;
}

}  // namespace BLA
