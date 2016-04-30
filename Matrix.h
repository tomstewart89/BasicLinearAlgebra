#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>
#include <string.h>

template<int rows, int cols = 1, class T = float> class Matrix
{
    static T dummyElement; // dummy variable returned when the () operator is used to access memory outside the matrix. It's a static so it shouldn't take up lots of memory

public:
    T m[rows*cols]; // underlying data storage - it's public but it's better to access the elements using ()
    bool invalidAccessFlag;  // a flag to indicate that an invalid access has occurred via the () operator or in the Set() function

    Matrix<rows,cols,T>() : invalidAccessFlag(false) { Clear(); }

    // Constructor to allow the matrix to be filled from an appropriately sized & typed 2D array
    Matrix<rows,cols,T>(T arr[rows][cols]) : invalidAccessFlag(false)
    {
        (*this) = arr;
    }

    Matrix<rows,cols,T>(const Matrix<rows,cols,T> &obj) : invalidAccessFlag(false)
    {
        memcpy(m,(T*)obj.m, rows*cols * sizeof(T));
    }

    // Unary addition
    Matrix<rows,cols,T> &operator+=(const Matrix<rows,cols,T> &obj)
    {
        for(int i = 0; i < rows * cols; i++)
            m[i] += obj.m[i];

        return *this;
    }

    // Unary subtraction
    Matrix<rows,cols,T> &operator-=(const Matrix<rows,cols,T> &obj)
    {
        for(int i = 0; i < rows * cols; i++)
            m[i] -= obj.m[i];

        return *this;
    }

    // Binary addition
    Matrix<rows,cols,T> operator+(const Matrix<rows,cols,T> &obj)
    {
        Matrix<rows,cols,T> tmp;

        for(int i = 0; i < rows * cols; i++)
            tmp.m[i] = m[i] + obj.m[i];

        return tmp;
    }

    // Binary subtraction
    Matrix<rows,cols,T> operator-(const Matrix<rows,cols,T> &obj)
    {
        Matrix<rows,cols,T> tmp;

        for(int i = 0; i < rows * cols; i++)
            tmp.m[i] = m[i] - obj.m[i];

        return tmp;
    }


    // Element-wise multiplication
    Matrix<rows,cols,T> &operator*=(T k)
    {
        for(int i = 0; i < rows * cols; i++)
            m[i] *= k;

        return *this;
    }

    // Copies the value of a matrix of equal size
    Matrix<rows,cols,T> &operator=(const Matrix<rows,cols,T> &obj)
    {
        memcpy(m,obj.m, rows*cols * sizeof(T));
        return *this;
    }

    // Set the contents of the matrix using a 2D array
    Matrix<rows,cols,T> &operator=(T arr[rows][cols])
    {
        memcpy(m,(T*)arr, rows*cols * sizeof(T));
        return *this;
    }

    // Binary matrix multiplication with
    template <int operandCols> Matrix<rows,operandCols,T> operator*(const Matrix<cols,operandCols,T> &operand)
    {
        Matrix<rows,operandCols,T> ret;
        int i,j,k;

        for(i = 0; i < rows; i++)
            for(j = 0; j < operandCols; j++)
                for(k = 0; k < cols; k++)
                    ret.m[i * operandCols + j] += m[i * cols + k] * operand.m[k * operandCols + j];

        return ret;
    }

    Matrix<rows,cols,T> &operator*=(const Matrix<rows,cols,T> &operand)
    {
        int i,j,k;
        Matrix<rows,cols,T> tmp;

        for(i = 0; i < rows; i++)
            for(j = 0; j < cols; j++)
                for(k = 0; k < cols; k++)
                    tmp.m[i * cols + j] += m[i * cols + k] * operand.m[k * cols + j];

        (*this) = tmp;
        return *this;
    }

    Matrix<rows,cols,T> operator*(T k)
    {
        Matrix<rows,cols,T> ret;

        // no need to make a nested loop here
        for(int i = 0; i < rows * cols; i++)
            ret.m[i] = m[i] * k;

        return ret;
    }

    // Set a subsection of the matrix (starting at i,j) to the value of the input object 'obj' if the
    template<int objRows, int objCols> void Set(const Matrix<objRows,objCols,T> &obj, int row, int col)
    {
        if(row + objRows > rows || col + objCols > cols)
            invalidAccessFlag = true;

        for(int i = 0; i < objRows && (i + row) < rows; i++)
            for(int j = 0; j < objCols && (j + col) < cols; j++)
                m[(i+row) * cols + (j+col)] = obj.m[i * objCols + j];
    }

    // Set the value of every element to 0
    void Clear() { memset(m,'\0',rows*cols*sizeof(T)); }

    T &operator()(int row, int col = 0)
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

    // Unary transpose operation - only supports square matrices
    void Transpose()
    {
        // Accessing the dummyElement from outside it's class isn't allowed, that means that the matrix must be square for this function to compile
        Matrix<rows,rows>::dummyElement;

        T tmp;
        for (int i = 0; i < rows; i++)
            for(int j = 0; j < rows; j++)
            {
                tmp = m[i * rows + j];
                m[i * rows + j] = m[j * rows + i];
                m[j * rows + i] = tmp;
            }
    }

    // Matrix Inversion Routine - modified from code written by Charlie Matlack: http://playground.arduino.cc/Code/MatrixMath
    // This function inverts a matrix based on the Gauss Jordan method. Specifically, it uses partial pivoting to improve numeric stability.
    // The algorithm is drawn from those presented in NUMERICAL RECIPES: The Art of Scientific Computing.
    void Invert(int *res = NULL)
    {
        // Accessing the dummyElement from outside it's class isn't allowed, that'll enforce that the Matrix be square
        Matrix<rows,rows>::dummyElement;

        int pivrow, pivrows[rows]; 	// keeps track of current pivot row and row swaps
        int i,j,k;
        T tmp;		// used for finding max value and making column swaps

        for (k = 0; k < rows; k++)
        {
            // find pivot row, the row with biggest entry in current column
            tmp = 0;
            for (i = k; i < rows; i++)
            {
                if (abs(m[i*rows+k]) >= tmp)
                {
                    tmp = abs(m[i*rows+k]);
                    pivrow = i;
                }
            }

            // check for singular matrix
            if (m[pivrow*rows+k] == 0.0f)
                if(res)
                    *res = 1;

            // Execute pivot (row swap) if needed
            if (pivrow != k)
            {
                // swap row k with pivrow
                for (j = 0; j < rows; j++)
                {
                    tmp = m[k*rows+j];
                    m[k*rows+j] = m[pivrow*rows+j];
                    m[pivrow*rows+j] = tmp;
                }
            }
            pivrows[k] = pivrow;	// record row swap (even if no swap happened)

            tmp = 1.0f / m[k*rows+k];	// invert pivot element
            m[k*rows+k] = 1.0f;		// This element of input matrix becomes result matrix

            // Perform row reduction (divide every element by pivot)
            for (j = 0; j < rows; j++)
                m[k*rows+j] = m[k*rows+j]*tmp;

            // Now eliminate all other entries in this column
            for (i = 0; i < rows; i++)
            {
                if (i != k)
                {
                    tmp = m[i*rows+k];
                    m[i*rows+k] = 0.0f;  // The other place where in matrix becomes result mat

                    for (j = 0; j < rows; j++)
                        m[i*rows+j] = m[i*rows+j] - m[k*rows+j]*tmp;
                }
            }
        }

        // Done, now need to undo pivot row swaps by doing column swaps in reverse order
        for (k = rows-1; k >= 0; k--)
        {
            if (pivrows[k] != k)
            {
                for (i = 0; i < rows; i++)
                {
                    tmp = m[i*rows+k];
                    m[i*rows+k] = m[i*rows+pivrows[k]];
                    m[i*rows+pivrows[k]] = tmp;
                }
            }
        }

        if(res)
            *res = 0;
    }

    int Rows() { return rows; }
    int Cols() { return cols; }
};

template <int rows, int cols, class T> T Matrix<rows,cols,T>::dummyElement;

// Create a new matrix by horizontally concatenating two input matrices
template<int rows, int cols, int operandCols, class T> Matrix<rows,cols + operandCols,T> HorzCat(const Matrix<rows,cols,T> &A, const Matrix<rows,operandCols,T> &B)
{
    Matrix<rows,cols + operandCols,T> ret;
    ret.Set(A,0,0);
    ret.Set(B,0,cols);

    return ret;
}

// Create a new matrix by vertically concatenating two input matrices
template<int rows, int cols, int operandRows, class T> Matrix<rows + operandRows,cols,T> VertCat(const Matrix<rows,cols,T> &A, const Matrix<operandRows,cols,T> &B)
{
    Matrix<rows + operandRows,cols,T> ret;
    ret.Set(A,0,0);
    ret.Set(B,rows,0);

    return ret;
}

// Multiply two matrices and store the result in a third matrix C, this is slightly faster than
template<int rows, int cols, int operandCols, class T> Matrix<rows,operandCols,T> &Multiply(const Matrix<rows,cols,T> &A, const Matrix<cols,operandCols,T> &B, Matrix<rows,operandCols,T> &C)
{
    int i,j,k;
    C.Clear();

    for(i = 0; i < rows; i++)
        for(j = 0; j < operandCols; j++)
            for(k = 0; k < cols; k++)
                C.m[i * operandCols + j] += A.m[i * cols + k] * B.m[k * operandCols + j];

    return C;
}

template<int rows, int cols, class T>  Matrix<rows,cols,T> &Add(const Matrix<rows,cols,T> &A, const Matrix<rows,cols,T> &B, Matrix<rows,cols,T> &C)
{
    for(int i = 0; i < rows * cols; i++)
        C.m[i] = A.m[i] + B.m[i];

    return C;
}

template<int rows, int cols, class T> Matrix<rows,cols,T> &Subtract(Matrix<rows,cols,T> &A, Matrix<rows,cols,T> &B, Matrix<rows,cols,T> &C)
{
    for(int i = 0; i < rows * cols; i++)
        C.m[i] = A.m[i] - B.m[i];

    return C;
}

// Binary invert operation - the result of the inversion is overwritten onto the input Matrix
template<int rows, int cols, class T> Matrix<rows,cols,T> &Invert(Matrix<rows,cols,T> &A, Matrix<rows,cols,T> &C, int *res = NULL)
{
    C = A;
    return Invert(C, res);
}

// Unary transpose operation - result is returned via copying
template<int rows, int cols, class T> Matrix<cols,rows,T> Transpose(const Matrix<rows,cols,T> &A)
{
    Matrix<cols,rows,T> C;

    for (int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C.m[i * cols + j] = A.m[j * cols + i];

    return C;
}

// Binary transpose operation - the result of the transpose goes into the input Matrix 'C'
template<int rows, int cols, class T> Matrix<cols,rows,T> &Transpose(const Matrix<rows,cols,T> &A, Matrix<cols,rows,T> &C)
{
    for (int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            C.m[i * cols + j] = A.m[j * cols + i];

    return C;
}

// Unary invert operation - result is returned via copying
template<int dim, class T> Matrix<dim,dim,T> Invert(const Matrix<dim,dim,T> &A, int *res = NULL)
{
    Matrix<dim,dim,T> C = A;
    C.Invert(res);

    return C;
}

// Binary transpose operation - the result of the inversion goes into the input Matrix 'C'. This is just to keep the interface consistent and doesn't save any memory
template<int rows, int cols, class T> Matrix<cols,rows,T> &Invert(const Matrix<rows,cols,T> &A, Matrix<cols,rows,T> &C, int *res = NULL)
{
    C = A;
    C.Invert(res);

    return C;
}

template<class T> inline Print &operator <<(Print &strm, const T &obj) { strm.print(obj); return strm; }

// Stream inserter operator for printing to strings or the serial port
template<int rows, int cols, class T> Print &operator<<(Print &strm, const Matrix<rows,cols,T> &obj)
{
    strm << '{';

    for(int i = 0; i < rows; i++)
    {
        strm << '{';

        for(int j = 0; j < cols; j++)
            strm << obj.m[i * cols + j] << ((j == cols - 1)? '}' : ',');

        strm << (i == rows - 1? '}' : ',');
    }
    return strm;
}

#endif // MATRIX_H

