#ifndef BLA_H
#define BLA_H

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Arduino.h"
#include "MemoryDelegate.hpp"

///////////////////////////////////////////////////////////////// Matrix ///////////////////////////////////////////////////////////////////

// Represents a Slice of rows or columns used in the () operator to select a submatrix - I'll replace this with an iterator eventually
template<int start, int end> struct Slice { };

template<int rows, int cols = 1, class ElemT = float, class MemT = Array<rows,cols,ElemT> > class Matrix
{
public:
    MemT delegate;

    // Constructors
    Matrix<rows,cols,ElemT,MemT>() { }
    Matrix<rows,cols,ElemT,MemT>(MemT &d) : delegate(d) { }
    Matrix<rows,cols,ElemT,MemT>(ElemT arr[rows][cols]) { *this = arr; }

    template<class opElemT, class opMemT>
    Matrix<rows,cols,ElemT,MemT>(const Matrix<rows,cols,opElemT,opMemT> &obj) { *this = obj; }

    // Element Access
    ElemT &operator()(int row, int col = 0) const;
    template<int rowStart, int rowEnd, int colStart, int colEnd> Matrix<rowEnd-rowStart,colEnd-colStart,ElemT,Reference<ElemT,MemT> > Submatrix(Slice<rowStart,rowEnd>, Slice<colStart,colEnd>) const;
    Matrix<rows,cols,ElemT,Reference<ElemT,MemT> > Ref() const;

    // Concatenation
    template<int operandCols, class opElemT, class opMemT> Matrix<rows,cols + operandCols,ElemT,HorzCat<cols,ElemT,MemT,opMemT> > operator||(const Matrix<rows,operandCols,opElemT,opMemT> &obj);
    template<int operandRows, class opElemT, class opMemT> Matrix<rows + operandRows,cols,ElemT,VertCat<rows,ElemT,MemT,opMemT> > operator&&(const Matrix<operandRows,cols,opElemT,opMemT> &obj);

    // Assignment
    template<class opElemT, class opMemT> Matrix<rows,cols,ElemT,MemT> &operator=(const Matrix<rows,cols,opElemT,opMemT> &obj);
    Matrix<rows,cols,ElemT,MemT> &operator=(ElemT arr[rows][cols]);
    Matrix<rows,cols,ElemT,MemT> &Fill(const ElemT &val);

    // Addition
    template<class opElemT, class opMemT> Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > operator+(const Matrix<rows,cols,opElemT,opMemT> &obj);
    template<class opElemT, class opMemT> Matrix<rows,cols,ElemT,MemT> &operator+=(const Matrix<rows,cols,opElemT,opMemT> &obj);

    // Subtraction
    template<class opElemT, class opMemT> Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > operator-(const Matrix<rows,cols,opElemT,opMemT> &obj);
    template<class opElemT, class opMemT> Matrix<rows,cols,ElemT,MemT> &operator-=(const Matrix<rows,cols,opElemT,opMemT> &obj);

    // Multiplication
    template <int operandCols, class opElemT, class opMemT> Matrix<rows,operandCols,ElemT,Array<rows,operandCols,ElemT> > operator*(const Matrix<cols,operandCols,opElemT,opMemT> &operand);
    template<class opElemT, class opMemT> Matrix<rows,cols,ElemT,MemT> &operator*=(const Matrix<rows,cols,opElemT,opMemT> &operand);

    // Negation
    Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > operator-();

    // Transposition
    Matrix<cols,rows,ElemT,Trans<ElemT,MemT> > operator~();

    // Scaling
    Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > operator*(ElemT k);
    Matrix<rows,cols,ElemT,MemT> &operator*=(ElemT k);

    // Returns the inverse of this matrix - only supports square matrices
    Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > Inverse(int *res = NULL);
    
    // Returns the determinant of this matrix
    ElemT Det() { return Determinant(*this); }

    int Rows() { return rows; }
    int Cols() { return cols; }
};

//////////////////////////////////////////////////////////////// Element Access ////////////////////////////////////////////////////////////////

template<int rows, int cols, class ElemT, class MemT>
ElemT &Matrix<rows,cols,ElemT,MemT>::operator()(int row, int col) const
{
    return delegate(row,col);
}

// this must have template arguments so that it can return a matrix of he right size, right? can they be inferred somehow>
template<int rows, int cols, class ElemT, class MemT>
template<int rowStart, int rowEnd, int colStart, int colEnd>
Matrix<rowEnd-rowStart,colEnd-colStart,ElemT,Reference<ElemT,MemT> > Matrix<rows,cols,ElemT,MemT>::Submatrix(Slice<rowStart,rowEnd>, Slice<colStart,colEnd>) const
{
    Reference<ElemT,MemT> ref(delegate, rowStart, rowEnd);
    return Matrix<rowEnd-rowStart,colEnd-colStart,ElemT,Reference<ElemT,MemT> >(ref);
}

template<int rows, int cols, class ElemT, class MemT>
Matrix<rows,cols,ElemT,Reference<ElemT,MemT> > Matrix<rows,cols,ElemT,MemT>::Ref() const
{
    Reference<ElemT,MemT> ref(delegate, 0,0);
    return Matrix<rows,cols,ElemT,Reference<ElemT,MemT> >(ref);
}

///////////////////////////////////////////////////////////////// Concatenation ////////////////////////////////////////////////////////////////

template<int rows, int cols, class ElemT, class MemT>
template<int operandCols, class opElemT, class opMemT>
Matrix<rows,cols+operandCols,ElemT,HorzCat<cols,ElemT,MemT,opMemT> > Matrix<rows,cols,ElemT,MemT>::operator||(const Matrix<rows,operandCols,opElemT,opMemT> &obj)
{
    HorzCat<cols,ElemT,MemT,opMemT> ref(delegate, obj.delegate);
    return Matrix<rows,cols+operandCols,ElemT,HorzCat<cols,ElemT,MemT,opMemT> >(ref);
}

template<int rows, int cols, class ElemT, class MemT>
template<int operandRows, class opElemT, class opMemT>
Matrix<rows+operandRows,cols,ElemT,VertCat<rows,ElemT,MemT,opMemT> > Matrix<rows,cols,ElemT,MemT>::operator&&(const Matrix<operandRows,cols,opElemT,opMemT> &obj)
{
    VertCat<rows,ElemT,MemT,opMemT> ref(delegate, obj.delegate);
    return Matrix<rows+operandRows,cols,ElemT,VertCat<rows,ElemT,MemT,opMemT> >(ref);
}



///////////////////////////////////////////////////////////////// Assignment ///////////////////////////////////////////////////////////////////

template<int rows, int cols, class ElemT, class MemT>
template<class opElemT, class opMemT>
Matrix<rows,cols,ElemT,MemT> &Matrix<rows,cols,ElemT,MemT>::operator=(const Matrix<rows,cols,opElemT,opMemT> &obj)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j)  = obj(i,j);

    return *this;
}

template<int rows, int cols, class ElemT, class MemT>
Matrix<rows,cols,ElemT,MemT> &Matrix<rows,cols,ElemT,MemT>::operator=(ElemT arr[rows][cols])
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j)  = arr[i][j];

    return *this;
}

template<int rows, int cols, class ElemT, class MemT>
Matrix<rows,cols,ElemT,MemT> &Matrix<rows,cols,ElemT,MemT>::Fill(const ElemT &val)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j) = val;

    return *this;
}

////////////////////////////////////////////////////////////////// Addition ////////////////////////////////////////////////////////////////////

template<int rows, int cols, class ElemT, class MemT>
template<class opElemT, class opMemT>
Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > Matrix<rows,cols,ElemT,MemT>::operator+(const Matrix<rows,cols,opElemT,opMemT> &obj)
{
    Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > ret;
    Add(*this,obj,ret);

    return ret;
}

template<int rows, int cols, class ElemT, class MemT>
template<class opElemT, class opMemT>
Matrix<rows,cols,ElemT,MemT> &Matrix<rows,cols,ElemT,MemT>::operator+=(const Matrix<rows,cols,opElemT,opMemT> &obj)
{
    Add(*this,obj,*this);

    return *this;
}

template<int rows, int cols, class ElemT, class MemT, class opElemT, class opMemT, class retElemT, class retMemT>
Matrix<rows,cols,retElemT,retMemT> &Add(const Matrix<rows,cols,ElemT,MemT> &A, const Matrix<rows,cols,opElemT,opMemT> &B, Matrix<rows,cols,retElemT,retMemT> &C)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C.delegate(i,j) = A.delegate(i,j) + B.delegate(i,j);

    return C;
}

////////////////////////////////////////////////////////////////// Subtraction /////////////////////////////////////////////////////////////////

template<int rows, int cols, class ElemT, class MemT>
template<class opElemT, class opMemT>
Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > Matrix<rows,cols,ElemT,MemT>::operator-(const Matrix<rows,cols,opElemT,opMemT> &obj)
{
    Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > ret;
    Subtract(*this,obj,ret);

    return ret;
}

template<int rows, int cols, class ElemT, class MemT>
template<class opElemT, class opMemT>
Matrix<rows,cols,ElemT,MemT> &Matrix<rows,cols,ElemT,MemT>::operator-=(const Matrix<rows,cols,opElemT,opMemT> &obj)
{
    Subtract(*this,obj,*this);

    return *this;
}

template<int rows, int cols, class ElemT, class MemT, class opElemT, class opMemT, class retElemT, class retMemT>
Matrix<rows,cols,retElemT,retMemT> &Subtract(const Matrix<rows,cols,ElemT,MemT> &A, const Matrix<rows,cols,opElemT,opMemT> &B, Matrix<rows,cols,retElemT,retMemT> &C)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C(i,j) = A(i,j) - B(i,j);

    return C;
}

//////////////////////////////////////////////////////////////// Multiplication ////////////////////////////////////////////////////////////////

template<int rows, int cols, class ElemT, class MemT>
template <int operandCols, class opElemT, class opMemT>
Matrix<rows,operandCols,ElemT,Array<rows,operandCols,ElemT> > Matrix<rows,cols,ElemT,MemT>::operator*(const Matrix<cols,operandCols,opElemT,opMemT> &operand)
{
    Matrix<rows,operandCols,ElemT,Array<rows,operandCols,ElemT> > ret;
    Multiply(*this,operand,ret);

    return ret;
}

template<int rows, int cols, class ElemT, class MemT>
template <class opElemT, class opMemT>
Matrix<rows,cols,ElemT,MemT> &Matrix<rows,cols,ElemT,MemT>::operator*=(const Matrix<rows,cols,opElemT,opMemT> &operand)
{
    Matrix<rows,cols,ElemT,MemT> ret;
    Multiply(*this,operand,ret);
    *this = ret;

    return *this;
}

// Multiplies two matrices and stores the result in a third matrix C, this is slightly faster than using the operators
template<int rows, int cols, int operandCols, class ElemT, class MemT, class opElemT, class opMemT, class retElemT, class retMemT>
Matrix<rows,operandCols,retElemT,retMemT> &Multiply(const Matrix<rows,cols,ElemT,MemT> &A, const Matrix<cols,operandCols,opElemT,opMemT> &B, Matrix<rows,operandCols,retElemT,retMemT> &C)
{
    int i,j,k;

    for(i = 0; i < rows; i++)
        for(j = 0; j < operandCols; j++)
        {
            if(cols > 0)
                C.delegate(i,j) = A.delegate(i,0) * B.delegate(0,j);

            for(k = 1; k < cols; k++)
                C.delegate(i,j) += A.delegate(i,k) * B.delegate(k,j);
        }

    return C;
}

/////////////////////////////////////////////////////////////////// Negation ///////////////////////////////////////////////////////////////////

template<int rows, int cols, class ElemT, class MemT>
Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > Matrix<rows,cols,ElemT,MemT>::operator-()
{
    Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > ret;

    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            ret(i,j) = -(*this)(i,j);

    return ret;
}

///////////////////////////////////////////////////////////////// Transposition ////////////////////////////////////////////////////////////////

template<int rows, int cols, class ElemT, class MemT>
Matrix<cols,rows,ElemT,Trans<ElemT,MemT> > Matrix<rows,cols,ElemT,MemT>::operator~()
{
    Trans<ElemT,MemT> ref(delegate);
    Matrix<cols,rows,ElemT,Trans<ElemT,MemT> >tmp(ref);

    return tmp;
}

//////////////////////////////////////////////////////////////////// Scaling ///////////////////////////////////////////////////////////////////

template<int rows, int cols, class ElemT, class MemT>
Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > Matrix<rows,cols,ElemT,MemT>::operator*(ElemT k)
{
    Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > ret;
    Scale(*this,k,ret);

    return ret;
}

template<int rows, int cols, class ElemT, class MemT>
Matrix<rows,cols,ElemT,MemT> &Matrix<rows,cols,ElemT,MemT>::operator*=(ElemT k)
{
    Scale(*this,k,*this);

    return *this;
}

// Multiplies two matrices and stores the result in a third matrix C, this is slightly faster than using the operators
template<int rows, int cols, int operandCols, class ElemT, class MemT, class opElemT, class retElemT, class retMemT>
Matrix<rows,operandCols,retElemT,retMemT> &Scale(const Matrix<rows,cols,ElemT,MemT> &A, const opElemT &B, Matrix<rows,operandCols,retElemT,retMemT> &C)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C(i,j) = A(i,j) * B;

    return C;
}

/////////////////////////////////////////////////////////////////// Determinant //////////////////////////////////////////////////////////////////

template<int dim, class ElemT, class MemT>
ElemT Determinant(const Matrix<dim,dim,ElemT,MemT> &A)
{
    ElemT det;

    // Add the determinants of all the minors
    for(int i = 0; i < dim; i++)
    {
        Minor<ElemT,MemT> del(A.delegate, i, 0);
        Matrix<dim-1,dim-1,ElemT, Minor<ElemT,MemT> > m(del);

        if(!i)
            det = Determinant(m);
        else
            det += Determinant(m) * (i % 2? -A(i,0) : A(i,0));
    }

    return det;
}

template<class ElemT, class MemT>
ElemT Determinant(Matrix<2,2,ElemT,MemT> &A)
{
    return A(0,0) * A(1,1) - A(1,0) * A(0,1);
}

/////////////////////////////////////////////////////////////////// Inversion //////////////////////////////////////////////////////////////////

template<int rows, int cols, class ElemT, class MemT>
Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > Matrix<rows,cols,ElemT,MemT>::Inverse(int *res)
{
    Matrix<rows,cols,ElemT,Array<rows,cols,ElemT> > ret(*this);
    return Invert(ret, res);
}

// Matrix Inversion Routine - modified from code written by Charlie Matlack: http://playground.arduino.cc/Code/MatrixMath
// This function inverts a matrix based on the Gauss Jordan method. Specifically, it uses partial pivoting to improve numeric stability.
// The algorithm is drawn from those presented in NUMERICAL RECIPES: The Art of Scientific Computing.
template<int dim, class ElemT, class MemT>
Matrix<dim,dim,ElemT,MemT> &Invert(Matrix<dim,dim,ElemT,MemT> &A, int *res = NULL)
{
    int pivrow, pivrows[dim]; 	// keeps track of current pivot row and row swaps
    int i,j,k;
    ElemT tmp;		// used for finding max value and making column swaps

    for (k = 0; k < dim; k++)
    {
        // find pivot row, the row with biggest entry in current column
        tmp = 0;
        for (i = k; i < dim; i++)
        {
            if(fabs(A(i,k)) >= tmp)
            {
                tmp = fabs(A(i,k));
                pivrow = i;
            }
        }

        // check for singular matrix
        if (A(pivrow,k) == 0.0f)
            if(res)
                *res = 1;

        // Execute pivot (row swap) if needed
        if (pivrow != k)
        {
            // swap row k with pivrow
            for (j = 0; j < dim; j++)
            {
                tmp = A(k,j);
                A(k,j) = A(pivrow,j);
                A(pivrow,j) = tmp;
            }
        }
        pivrows[k] = pivrow;	// record row swap (even if no swap happened)

        tmp = 1.0f / A(k,k);	// invert pivot element
        A(k,k) = 1.0f;		// This element of input matrix becomes result matrix

        // Perform row reduction (divide every element by pivot)
        for (j = 0; j < dim; j++)
            A(k,j) = A(k,j) * tmp;

        // Now eliminate all other entries in this column
        for (i = 0; i < dim; i++)
        {
            if (i != k)
            {
                tmp = A(i,k);
                A(i,k) = 0.0f;  // The other place where in matrix becomes result mat

                for (j = 0; j < dim; j++)
                    A(i,j) = A(i,j) - A(k,j) * tmp;
            }
        }
    }

    // Done, now need to undo pivot row swaps by doing column swaps in reverse order
    for (k = dim-1; k >= 0; k--)
    {
        if (pivrows[k] != k)
        {
            for (i = 0; i < dim; i++)
            {
                tmp = A(i,k);
                A(i,k) = A(i,pivrows[k]);
                A(i,pivrows[k]) = tmp;
            }
        }
    }

    if(res)
        *res = 0;

    return A;
}

////////////////////////////////////////////////////////////////// Insertion ///////////////////////////////////////////////////////////////////

inline Print &operator <<(Print &strm, const int obj)
{
    strm.print(obj); return strm;
}

inline Print &operator <<(Print &strm, const float obj)
{
    strm.print(obj); return strm;
}

inline Print &operator <<(Print &strm, const char *obj)
{
    strm.print(obj); return strm;
}

inline Print &operator <<(Print &strm, const char obj)
{
    strm.print(obj); return strm;
}

// Stream inserter operator for printing to strings or the serial port
template<int rows, int cols, class ElemT, class MemT>
Print &operator<<(Print &strm, const Matrix<rows,cols,ElemT,MemT> &obj)
{
    strm << '{';

    for(int i = 0; i < rows; i++)
    {
        strm << '{';

        for(int j = 0; j < cols; j++)
            strm << obj(i,j) << ((j == cols - 1)? '}' : ',');

        strm << (i == rows - 1? '}' : ',');
    }
    return strm;
}


#endif // BLA_H
