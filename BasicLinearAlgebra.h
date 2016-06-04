#ifndef BLA_H
#define BLA_H

#include "Arduino.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

template<int rows, int cols = 1, class T = float> class Matrix
{
    static T dummyElement; // dummy variable returned when the () operator is used to access memory outside the matrix. It's a static so it shouldn't take up lots of memory
    mutable T m[rows*cols]; // underlying data storage - it's public but it's better to access the elements using ()
public:

    mutable bool invalidAccessFlag;  // a flag to indicate that an invalid access has occurred via the () operator or in the Set() function

    Matrix<rows,cols,T>() : invalidAccessFlag(false) { /*Clear();*/ }

    // Constructor to allow the matrix to be filled from an appropriately sized & typed 2D array
    Matrix<rows,cols,T>(T arr[rows][cols]) : invalidAccessFlag(false)
    {
        (*this) = arr;
    }

    Matrix<rows,cols,T>(const Matrix<rows,cols,T> &obj) : invalidAccessFlag(false)
    {
        (*this) = obj;
    }

    T &operator()(int row, int col = 0) const
    {
        // Accessing the dummyElement from outside it's class isn't allowed, that means that the matrix must be square for this function to compile
        if(row > rows || col > cols)
        {
            invalidAccessFlag = true;
            return dummyElement;
        }
        else
            return m[row * cols + col];
    }

    // Unary addition
    template<class opT> Matrix<rows,cols,T> &operator+=(const Matrix<rows,cols,opT> &obj)
    {
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                (*this)(i,j) += obj(i,j);

        return *this;
    }

    // Unary subtraction
    template<class opT> Matrix<rows,cols,T> &operator-=(const Matrix<rows,cols,opT> &obj)
    {
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                (*this)(i,j) -= obj(i,j);

        return *this;
    }

    // Binary addition
    template<class opT> Matrix<rows,cols,T> operator+(const Matrix<rows,cols,opT> &obj)
    {
        Matrix<rows,cols,T> ret;

        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                ret(i,j) = (*this)(i,j) + obj(i,j);

        return ret;
    }

    // Binary subtraction
    template<class opT> Matrix<rows,cols,T> operator-(const Matrix<rows,cols,opT> &obj)
    {
        Matrix<rows,cols,T> ret;

        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                ret(i,j) = (*this)(i,j) - obj(i,j);

        return ret;
    }

    // Negation
    Matrix<rows,cols,T> operator-()
    {
        Matrix<rows,cols,T> ret;

        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                ret(i,j) = -(*this)(i,j);

        return ret;
    }

    // Binary matrix multiplication
    template <int operandCols, class opT> Matrix<rows,operandCols,T> operator*(const Matrix<cols,operandCols,opT> &operand)
    {
        Matrix<rows,operandCols,T> ret;
        int i,j,k;

        for(i = 0; i < rows; i++)
            for(j = 0; j < operandCols; j++)
                for(k = 0; k < cols; k++)
                    ret(i,j) += (*this)(i,k) * operand(k,j);

        return ret;
    }

    template<class opT> Matrix<rows,cols,T> &operator*=(const Matrix<rows,cols,opT> &operand)
    {
        int i,j,k;
        Matrix<rows,cols,T> ret;

        for(i = 0; i < rows; i++)
            for(j = 0; j < cols; j++)
                for(k = 0; k < cols; k++)
                    ret(i,j) += (*this)(i,k) * operand(k,j);

        *this = ret;
        return *this;
    }

    // Scaling
    Matrix<rows,cols,T> &operator*=(T k)
    {
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                (*this)(i,j) *= k;

        return *this;
    }

    Matrix<rows,cols,T> operator*(T k)
    {
        Matrix<rows,cols,T> ret;

        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                ret(i,j) = (*this)(i,j) * k;

        return ret;
    }

    // Copies the value of a matrix of equal size
    template<class opT> Matrix<rows,cols,T> &operator=(const Matrix<rows,cols,opT> &obj)
    {
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                (*this)(i,j)  = obj(i,j);

        return *this;
    }

    // Set the contents of the matrix using a 2D array
    Matrix<rows,cols,T> &operator=(T arr[rows][cols])
    {
        for(int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                (*this)(i,j)  = arr[i][j];

        return *this;
    }

    // Set a subsection of the matrix of dimensions (width,height & starting at destRow,destCol) to the value of the input object 'obj' starting at (srcRow, srcCol)
    template<int objRows, int objCols, class opT> void Set(const Matrix<objRows,objCols,opT> &obj, int destRow = 0, int destCol = 0, int srcRow = 0, int srcCol = 0, int height = objRows, int width = objCols)
    {
        // first make sure that the for loops won't try to access memory outside the matrices - this way the functions defaults to taking as much as will fit
        height = min(height, min(rows - destRow, objRows - srcRow));
        width = min(width,min(cols - destCol,objCols - srcCol));

        for(int i = 0; i < height; i++)
            for(int j = 0; j < width; j++)
                (*this)(i+destRow,j+destCol) = obj(i+srcRow,j+srcCol);
    }

    // Set a subsection of the matrix of dimensions (width,height & starting at destRow,destCol) to the value of the input object 'obj' starting at (srcRow, srcCol)
    template<int subRows, int subCols> Matrix<subRows,subCols,T> Submatrix(int srcRow, int srcCol)
    {
        Matrix<subRows,subCols,T> ret;
        ret.Set(*this,0,0,srcRow,srcCol,subRows,subCols);
        return ret;
    }

    // Returns a transpose of this matrix
    Matrix<cols,rows,T> Transpose()
    {
        Matrix<cols,rows,T> ret;

        for (int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                ret(i,j) = (*this)(j,i);

        return ret;
    }

    // Returns the inverse of this matrix - only supports square matrices
    Matrix<rows,cols,T> Inverse(int *res = NULL)
    {
        Matrix<rows,cols,T> ret = *this;
        return Invert(ret, res);
    }

    int Rows() { return rows; }
    int Cols() { return cols; }
};

template <int rows, int cols, class T> T Matrix<rows,cols,T>::dummyElement;

// Multiply two matrices and store the result in a third matrix C, this is slightly faster than using the operator
template<int rows, int cols, int operandCols, class T, class opT, class retT> Matrix<rows,operandCols,retT> &Multiply(const Matrix<rows,cols,T> &A, const Matrix<cols,operandCols,opT> &B, Matrix<rows,operandCols,retT> &C)
{
    int i,j,k;
    C.Clear();

    for(i = 0; i < rows; i++)
        for(j = 0; j < operandCols; j++)
            for(k = 0; k < cols; k++)
                C(i,j) += A(i,k) * B(k,j);

    return C;
}

template<int rows, int cols, class T, class opT, class retT>  Matrix<rows,cols,retT> &Add(const Matrix<rows,cols,T> &A, const Matrix<rows,cols,opT> &B, Matrix<rows,cols,retT> &C)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C(i,j) = A(i,j) + B(i,j);

    return C;
}

template<int rows, int cols, class T, class opT, class retT> Matrix<rows,cols,retT> &Subtract(Matrix<rows,cols,T> &A, Matrix<rows,cols,opT> &B, Matrix<rows,cols,retT> &C)
{
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C(i,j) = A(i,j) - B(i,j);

    return C;
}

// Create a new matrix by horizontally concatenating two input matrices
template<int rows, int cols, int operandCols, class T, class opT> Matrix<rows,cols + operandCols,T> HorzCat(const Matrix<rows,cols,T> &A, const Matrix<rows,operandCols,opT> &B)
{
    Matrix<rows,cols + operandCols,T> ret;
    ret.Set(A,0,0);
    ret.Set(B,0,cols);

    return ret;
}

// Create a new matrix by vertically concatenating two input matrices
template<int rows, int cols, int operandRows, class T, class opT> Matrix<rows + operandRows,cols,T> VertCat(const Matrix<rows,cols,T> &A, const Matrix<operandRows,cols,opT> &B)
{
    Matrix<rows + operandRows,cols,T> ret;
    ret.Set(A,0,0);
    ret.Set(B,rows,0);

    return ret;
}

// Transpose operation - the result of the transpose goes into the input Matrix 'C'
template<int rows, int cols, class T, class retT> Matrix<cols,rows,retT> &Transpose(const Matrix<rows,cols,T> &A, Matrix<cols,rows,retT> &C)
{
    for (int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C(i,j) = A(j,i);

    return C;
}

// Matrix Inversion Routine - modified from code written by Charlie Matlack: http://playground.arduino.cc/Code/MatrixMath
// This function inverts a matrix based on the Gauss Jordan method. Specifically, it uses partial pivoting to improve numeric stability.
// The algorithm is drawn from those presented in NUMERICAL RECIPES: The Art of Scientific Computing.
template<int dim, class T> Matrix<dim,dim,T> &Invert(Matrix<dim,dim,T> &A, int *res = NULL)
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

template<class T> inline Print &operator <<(Print &strm, const T &obj)
{
    strm.print(obj); return strm;
}

// Stream inserter operator for printing to strings or the serial port
template<int rows, int cols, class T> Print &operator<<(Print &strm, const Matrix<rows,cols,T> &obj)
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

