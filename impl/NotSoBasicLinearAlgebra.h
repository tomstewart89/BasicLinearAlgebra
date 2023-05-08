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

template <int dim, class MemT>
struct LUDecomposition
{
    bool singular;
    typename MemT::elem_t parity;
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
    decomp.parity = 1.0;

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

        for (int j = 0; j < dim; ++j)
        {
            typename MemT::elem_t this_elem = fabs(A(i, j));
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

template <int dim, class LUType, class BType>
Matrix<dim, 1, > LUSolve(const LUDecomposition<dim, MemT1> &decomp, const Matrix<dim, 1, MemT2> &b)
{
    ArrayMatrix<dim, 1, typename MemT2::elem_t> x, tmp;

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
bool Invert(const Matrix<dim, dim, MemT> &A, Matrix<dim, dim, MemT> &out)
{
    ArrayMatrix<dim, dim, typename MemT::elem_t> A_copy = A;

    auto decomp = LUDecompose(A_copy);

    if (decomp.singular)
    {
        return false;
    }

    ArrayMatrix<dim, 1, typename MemT::elem_t> b = Zeros<dim>();

    for (int j = 0; j < dim; ++j)
    {
        b(j) = 1.0;
        out.Column(j) = LUSolve(decomp, b);
        b(j) = 0.0;
    }

    return true;
}

template <int dim, class MemT>
bool Invert(Matrix<dim, dim, MemT> &A)
{
    return Invert(A, A);
}

template <int dim, class MemT>
Matrix<dim, dim, MemT> Inverse(const Matrix<dim, dim, MemT> &A)
{
    ArrayMatrix<dim, dim, typename MemT::elem_t> out;
    Invert(A, out);
    return out;
}

template <int dim, class MemT>
typename MemT::elem_t Determinant(const Matrix<dim, dim, MemT> &A)
{
    Matrix<dim, dim> A_copy = A;

    auto decomp = LUDecompose(A_copy);

    typename MemT::elem_t det = decomp.parity;

    for (int i = 0; i < dim; ++i)
    {
        det *= decomp.upper(i, i);
    }

    return det;
}

template <int rows, int cols, class MemT>
typename MemT::elem_t Norm(const Matrix<rows, cols, MemT> &A)
{
    typename MemT::elem_t sum_sq = 0.0;

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            sum_sq += A(i, j) * A(i, j);
        }
    }
    return sqrt(sum_sq);
}

template <int rows, int cols, class MemT>
typename MemT::elem_t Trace(const Matrix<rows, cols, MemT> &A)
{
    typename MemT::elem_t sum_diag = 0.0;

    for (int i = 0; i < rows; ++i)
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
