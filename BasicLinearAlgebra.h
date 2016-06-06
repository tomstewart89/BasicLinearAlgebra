#ifndef BLA_H
#define BLA_H

#include "Arduino.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Represents a range of rows or columns used in the () operator to select a submatrix
template<int length> struct Range
{
    int offset;
    Range(int off) : offset(off) {}
};

// Promises that the subclass can be addressed using the () operator
template<class T> struct Indexable
{
    virtual T &operator()(int row, int col = 0, bool *invalidAccess = NULL) const = 0;
};

////////////////////////////////////////////////////////////////// Matrix ///////////////////////////////////////////////////////////////////

template<int rows, int cols = 1, class T = float> class Matrix : public Indexable<T>
{
    mutable T m[rows*cols]; // underlying data storage - access it with ()
public:

    // Constructors
    Matrix<rows,cols,T>() { } // from nothing
    Matrix<rows,cols,T>(T arr[rows][cols]) { (*this) = arr; } // from an array
    template<class opT> Matrix<rows,cols,T>(const Matrix<rows,cols,opT> &obj) { (*this) = obj; } // from another matrix

    // Assignment
    template<class opT> Matrix<rows,cols,T> &operator=(const Matrix<rows,cols,opT> &obj);
    Matrix<rows,cols,T> &operator=(T arr[rows][cols]);
    Matrix<rows,cols,T> &Fill(const T &val);

    // Element Access
    T &operator()(int row, int col, bool *invalidAccess) const;
    template<int height, int width> Matrix<height,width,T&> operator()(Range<height> rowRange, Range<width> colRange) const;

    // Addition
    template<class opT> Matrix<rows,cols,T> operator+(const Matrix<rows,cols,opT> &obj);
    template<class opT> Matrix<rows,cols,T> &operator+=(const Matrix<rows,cols,opT> &obj);

    // Subtraction
    template<class opT> Matrix<rows,cols,T> operator-(const Matrix<rows,cols,opT> &obj);
    template<class opT> Matrix<rows,cols,T> &operator-=(const Matrix<rows,cols,opT> &obj);

    // Multiplication
    template <int operandCols, class opT> Matrix<rows,operandCols,T> operator*(const Matrix<cols,operandCols,opT> &operand);
    template<class opT> Matrix<rows,cols,T> &operator*=(const Matrix<rows,cols,opT> &operand);

    // Negation
    Matrix<rows,cols,T> operator-();

    // Scaling
    Matrix<rows,cols,T> operator*(T k);
    Matrix<rows,cols,T> &operator*=(T k);

    // Returns a transpose of this matrix
    Matrix<cols,rows,T> Transpose();

    // Returns the inverse of this matrix - only supports square matrices
    Matrix<rows,cols,T> Inverse(int *res);

    int Rows() { return rows; }
    int Cols() { return cols; }
};

////////////////////////////////////////////////////////////////// Reference Matrix ///////////////////////////////////////////////////////////////////

template<int rows, int cols, class T> class Matrix<rows,cols,T&> : Indexable<T>
{
    const Indexable<T> &parent; // reference to the parent matrix
    int rowOffset, colOffset;
    template<int, int, class> friend class Matrix; // Matrices are all friends with each other

public:
    // Constructors
    Matrix<rows,cols,T&>(const Indexable<T> &obj, int rowOff, int colOff) : parent(obj), rowOffset(rowOff), colOffset(colOff) { } // from any Indexable type

    template<int parentRows, int parentCols>  // from a reference matrix - this would also work through the base class but this'll save indirection when making refs of refs
    Matrix<rows,cols,T&>(const Matrix<parentRows,parentCols,T&> &obj, int rowOff, int colOff) : parent(obj.parent), rowOffset(rowOff+obj.rowOffset), colOffset(colOff+obj.colOffset) { }

    // Element Access
    T &operator()(int row, int col, bool *invalidAccess) const;
    template<int height, int width> Matrix<height,width,T&> operator()(Range<height> rowRange, Range<width> colRange) const;

    // Assignment
    template<class opT> Matrix<rows,cols,T&> &operator=(const Matrix<rows,cols,opT> &obj);
    Matrix<rows,cols,T&> &operator=(T arr[rows][cols]);
    Matrix<rows,cols,T&> &Fill(const T &val);

    // Addition
    template<class opT> Matrix<rows,cols,T> operator+(const Matrix<rows,cols,opT> &obj);
    template<class opT> Matrix<rows,cols,T&> &operator+=(const Matrix<rows,cols,opT> &obj);

    // Subtraction
    template<class opT> Matrix<rows,cols,T> operator-(const Matrix<rows,cols,opT> &obj);
    template<class opT> Matrix<rows,cols,T&> &operator-=(const Matrix<rows,cols,opT> &obj);

    // Multiplication
    template <int operandCols, class opT> Matrix<rows,operandCols,T> operator*(const Matrix<cols,operandCols,opT> &operand);
    template<class opT> Matrix<rows,cols,T&> &operator*=(const Matrix<rows,cols,opT> &operand);

    // Negation
    Matrix<rows,cols,T> operator-();

    // Scaling
    Matrix<rows,cols,T> operator*(T k);
    Matrix<rows,cols,T&> &operator*=(T k);

    // Returns a transpose of this matrix
    Matrix<cols,rows,T> Transpose();

    // Returns the inverse of this matrix - only supports square matrices
    Matrix<rows,cols,T> Inverse(int *res);

    int Rows() { return rows; }
    int Cols() { return cols; }
};

///////////////////////////////////////////////////////////////////////// Assignment ///////////////////////////////////////////////////////////////////

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T> &Matrix<rows,cols,T>::operator=(const Matrix<rows,cols,opT> &obj)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j)  = obj(i,j);

    return *this;
}

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T&> &Matrix<rows,cols,T&>::operator=(const Matrix<rows,cols,opT> &obj)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j)  = obj(i,j);

    return *this;
}

template<int rows, int cols, class T>
Matrix<rows,cols,T> &Matrix<rows,cols,T>::operator=(T arr[rows][cols])
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j)  = arr[i][j];

    return *this;
}

template<int rows, int cols, class T>
Matrix<rows,cols,T&> &Matrix<rows,cols,T&>::operator=(T arr[rows][cols])
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j)  = arr[i][j];

    return *this;
}

template<int rows, int cols, class T>
Matrix<rows,cols,T> &Matrix<rows,cols,T>::Fill(const T &val)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j)  = val;

    return *this;
}

template<int rows, int cols, class T>
Matrix<rows,cols,T&> &Matrix<rows,cols,T&>::Fill(const T &val)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j)  = val;

    return *this;
}

///////////////////////////////////////////////////////////////////////// Element Access ///////////////////////////////////////////////////////////////////

template<int rows, int cols, class T>
T &Matrix<rows,cols,T>::operator()(int row, int col = 0, bool *invalidAccess = NULL) const
{
    static T dummy;

    if(row > rows || col > cols)
    {
        if(invalidAccess)
            *invalidAccess = true;

        return dummy;
    }
    else
    {
        if(invalidAccess)
            *invalidAccess = false;

        return m[row * cols + col];
    }
}

template<int rows, int cols, class T>
T &Matrix<rows,cols,T&>::operator()(int row, int col = 0, bool *invalidAccess = NULL) const
{
    return parent(row + rowOffset, col + colOffset, invalidAccess);
}

template<int rows, int cols, class T>
template<int height, int width> Matrix<height,width,T&> Matrix<rows,cols,T>::operator()(Range<height> rowRange, Range<width> colRange) const
{
    return Matrix<height,width,T&>(*this, rowRange.offset, colRange.offset);
}

template<int rows, int cols, class T>
template<int height, int width> Matrix<height,width,T&> Matrix<rows,cols,T&>::operator()(Range<height> rowRange, Range<width> colRange) const
{
    rowRange.offset += rowOffset;
    colRange.offset += colOffset;

    return Matrix<height,width,T&>(parent, rowRange.offset, colRange.offset); // so I don't know anything about the parent at this point, which can't be helped. But I can't return
}

//////////////////////////////////////////////////////////////////////////// Addition ///////////////////////////////////////////////////////////////////////

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T> Matrix<rows,cols,T>::operator+(const Matrix<rows,cols,opT> &obj)
{
    Matrix<rows,cols,T> ret;

    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            ret(i,j) = (*this)(i,j) + obj(i,j);

    return ret;
}

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T> Matrix<rows,cols,T&>::operator+(const Matrix<rows,cols,opT> &obj)
{
    Matrix<rows,cols,T> ret;

    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            ret(i,j) = (*this)(i,j) + obj(i,j);

    return ret;
}

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T> &Matrix<rows,cols,T>::operator+=(const Matrix<rows,cols,opT> &obj)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j) += obj(i,j);

    return *this;
}

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T&> &Matrix<rows,cols,T&>::operator+=(const Matrix<rows,cols,opT> &obj)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j) += obj(i,j);

    return *this;
}

template<int rows, int cols, class T, class opT, class retT>  Matrix<rows,cols,retT> &Add(const Matrix<rows,cols,T> &A, const Matrix<rows,cols,opT> &B, Matrix<rows,cols,retT> &C)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C(i,j) = A(i,j) + B(i,j);

    return C;
}

//////////////////////////////////////////////////////////////////////////// Subtraction ///////////////////////////////////////////////////////////////////////

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T> Matrix<rows,cols,T>::operator-(const Matrix<rows,cols,opT> &obj)
{
    Matrix<rows,cols,T> ret;

    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            ret(i,j) = (*this)(i,j) - obj(i,j);

    return ret;
}

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T> Matrix<rows,cols,T&>::operator-(const Matrix<rows,cols,opT> &obj)
{
    Matrix<rows,cols,T> ret;

    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            ret(i,j) = (*this)(i,j) - obj(i,j);

    return ret;
}

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T> &Matrix<rows,cols,T>::operator-=(const Matrix<rows,cols,opT> &obj)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j) -= obj(i,j);

    return *this;
}

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T&> &Matrix<rows,cols,T&>::operator-=(const Matrix<rows,cols,opT> &obj)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j) -= obj(i,j);

    return *this;
}

template<int rows, int cols, class T, class opT, class retT>
Matrix<rows,cols,retT> &Subtract(Matrix<rows,cols,T&> &A, Matrix<rows,cols,opT> &B, Matrix<rows,cols,retT> &C)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C(i,j) = A(i,j) - B(i,j);

    return C;
}

///////////////////////////////////////////////////////////////////////// Multiplication ///////////////////////////////////////////////////////////////////////

template<int rows, int cols, class T>
template <int operandCols, class opT> Matrix<rows,operandCols,T> Matrix<rows,cols,T>::operator*(const Matrix<cols,operandCols,opT> &operand)
{
    Matrix<rows,operandCols,T> ret;
    int i,j,k;

    for(i = 0; i < rows; i++)
        for(j = 0; j < operandCols; j++)
        {
            if(cols > 0)
                ret(i,j) = (*this)(i,0) * operand(0,j);

            for(k = 1; k < cols; k++)
                ret(i,j) += (*this)(i,k) * operand(k,j);
        }

    return ret;
}

template<int rows, int cols, class T>
template <int operandCols, class opT> Matrix<rows,operandCols,T> Matrix<rows,cols,T&>::operator*(const Matrix<cols,operandCols,opT> &operand)
{
    Matrix<rows,operandCols,T> ret;
    int i,j,k;

    for(i = 0; i < rows; i++)
        for(j = 0; j < operandCols; j++)
        {
            if(cols > 0)
                ret(i,j) = (*this)(i,0) * operand(0,j);

            for(k = 1; k < cols; k++)
                ret(i,j) += (*this)(i,k) * operand(k,j);
        }

    return ret;
}

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T> &Matrix<rows,cols,T>::operator*=(const Matrix<rows,cols,opT> &operand)
{
    int i,j,k;
    Matrix<rows,cols,T> ret;

    for(i = 0; i < rows; i++)
        for(j = 0; j < cols; j++)
        {
            if(cols > 0)
                ret(i,j) = (*this)(i,0) * operand(0,j);

            for(k = 1; k < cols; k++)
                ret(i,j) += (*this)(i,k) * operand(k,j);
        }

    *this = ret;
    return *this;
}

template<int rows, int cols, class T>
template<class opT> Matrix<rows,cols,T&> &Matrix<rows,cols,T&>::operator*=(const Matrix<rows,cols,opT> &operand)
{
    int i,j,k;
    Matrix<rows,cols,T> ret;

    for(i = 0; i < rows; i++)
        for(j = 0; j < cols; j++)
        {
            if(cols > 0)
                ret(i,j) = (*this)(i,0) * operand(0,j);

            for(k = 1; k < cols; k++)
                ret(i,j) += (*this)(i,k) * operand(k,j);
        }

    *this = ret;
    return *this;
}

// Multiplies two matrices and stores the result in a third matrix C, this is slightly faster than using the operators
template<int rows, int cols, int operandCols, class T, class opT, class retT>
Matrix<rows,operandCols,retT> &Multiply(const Matrix<rows,cols,T> &A, const Matrix<cols,operandCols,opT> &B, Matrix<rows,operandCols,retT> &C)
{
    int i,j,k;

    for(i = 0; i < rows; i++)
        for(j = 0; j < cols; j++)
        {
            if(cols > 0)
                C(i,j) = A(i,0) * B(0,j);

            for(k = 1; k < cols; k++)
                C(i,j) += A(i,k) * B(k,j);
        }

    return C;
}

///////////////////////////////////////////////////////////////////////// Negation ///////////////////////////////////////////////////////////////////

template<int rows, int cols, class T>
Matrix<rows,cols,T> Matrix<rows,cols,T>::operator-()
{
    Matrix<rows,cols,T> ret;

    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            ret(i,j) = -(*this)(i,j);

    return ret;
}

template<int rows, int cols, class T>
Matrix<rows,cols,T> Matrix<rows,cols,T&>::operator-()
{
    Matrix<rows,cols,T> ret;

    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            ret(i,j) = -(*this)(i,j);

    return ret;
}

///////////////////////////////////////////////////////////////////////// Scaling ///////////////////////////////////////////////////////////////////

template<int rows, int cols, class T>
Matrix<rows,cols,T> Matrix<rows,cols,T>::operator*(T k)
{
    Matrix<rows,cols,T> ret;

    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            ret(i,j) = (*this)(i,j) * k;

    return ret;
}

template<int rows, int cols, class T>
Matrix<rows,cols,T> Matrix<rows,cols,T&>::operator*(T k)
{
    Matrix<rows,cols,T> ret;

    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            ret(i,j) = (*this)(i,j) * k;

    return ret;
}

template<int rows, int cols, class T>
Matrix<rows,cols,T> &Matrix<rows,cols,T>::operator*=(T k)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j) *= k;

    return *this;
}

template<int rows, int cols, class T>
Matrix<rows,cols,T&> &Matrix<rows,cols,T&>::operator*=(T k)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            (*this)(i,j) *= k;

    return *this;
}

///////////////////////////////////////////////////////////////////////// Transposition ///////////////////////////////////////////////////////////////////

template<int rows, int cols, class T>
Matrix<cols,rows,T> Matrix<rows,cols,T>::Transpose()
{   
    Matrix<cols,rows,T> ret;

    for (int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            ret(j,i) = (*this)(i,j);

    return ret;
}

template<int rows, int cols, class T>
Matrix<cols,rows,T> Matrix<rows,cols,T&>::Transpose()
{
    Matrix<cols,rows,T> ret;

    for (int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            ret(j,i) = (*this)(i,j);

    return ret;
}

template<int rows, int cols, class T, class retT>
Matrix<cols,rows,retT> &Transpose(const Matrix<rows,cols,T> &A, Matrix<cols,rows,retT> &C)
{
    for (int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C(j,i) = A(i,j);

    return C;
}

///////////////////////////////////////////////////////////////////////// Inversion ///////////////////////////////////////////////////////////////////

template<int rows, int cols, class T>
Matrix<rows,cols,T> Matrix<rows,cols,T>::Inverse(int *res = NULL)
{
    Matrix<rows,cols,T> ret = *this;
    return Invert(ret, res);
}

template<int rows, int cols, class T>
Matrix<rows,cols,T> Matrix<rows,cols,T&>::Inverse(int *res = NULL)
{
    Matrix<rows,cols,T> ret = *this;
    return Invert(ret, res);
}

// Matrix Inversion Routine - modified from code written by Charlie Matlack: http://playground.arduino.cc/Code/MatrixMath
// This function inverts a matrix based on the Gauss Jordan method. Specifically, it uses partial pivoting to improve numeric stability.
// The algorithm is drawn from those presented in NUMERICAL RECIPES: The Art of Scientific Computing.
template<int dim, class T>
Matrix<dim,dim,T> &Invert(Matrix<dim,dim,T> &A, int *res = NULL)
{
    int pivrow, pivrows[dim]; 	// keeps track of current pivot row and row swaps
    int i,j,k;
    T tmp;		// used for finding max value and making column swaps

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

template<int dim, class T>
Matrix<dim,dim,T&> &Invert(Matrix<dim,dim,T&> &A, int *res = NULL)
{
    int pivrow, pivrows[dim]; 	// keeps track of current pivot row and row swaps
    int i,j,k;
    T tmp;		// used for finding max value and making column swaps

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

///////////////////////////////////////////////////////////////////////// Concatenation ///////////////////////////////////////////////////////////////////

template<int rows, int cols, int operandCols, class T, class opT>
Matrix<rows,cols+operandCols,T> HorzCat(const Matrix<rows,cols,T> &A, const Matrix<rows,operandCols,opT> &B)
{
    Matrix<rows,cols + operandCols,T> ret;
    ret(Range<rows>(0),Range<cols>(0)) = A;
    ret(Range<rows>(0),Range<operandCols>(cols)) = B;

    return ret;
}

template<int rows, int cols, int operandRows, class T, class opT>
Matrix<rows + operandRows,cols,T> VertCat(const Matrix<rows,cols,T> &A, const Matrix<operandRows,cols,opT> &B)
{
    Matrix<rows + operandRows,cols,T> ret;
    ret(Range<rows>(0),Range<cols>(0)) = A;
    ret(Range<operandRows>(rows),Range<cols>(0)) = B;

    return ret;
}

///////////////////////////////////////////////////////////////////////// Insertion ///////////////////////////////////////////////////////////////////

template<class T>
inline Print &operator <<(Print &strm, const T &obj)
{
    strm.print(obj); return strm;
}

// Stream inserter operator for printing to strings or the serial port
template<int rows, int cols, class T>
Print &operator<<(Print &strm, const Matrix<rows,cols,T> &obj)
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
