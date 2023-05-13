#pragma once

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Arduino.h"

namespace BLA
{
template <typename MatAType, typename MatBType, int MatARows, int MatACols, int MatBCols, typename DType>
Matrix<MatARows, MatBCols, DType> operator*(const MatrixBase<MatAType, MatARows, MatACols, DType> &matA,
                                            const MatrixBase<MatBType, MatACols, MatBCols, DType> &matB)
{
    Matrix<MatARows, MatBCols, DType> ret;

    for (int i = 0; i < MatARows; ++i)
    {
        for (int j = 0; j < MatBCols; ++j)
        {
            if (MatACols > 0)
            {
                ret(i, j) = matA(i, 0) * matB(0, j);
            }

            for (int k = 1; k < MatACols; k++)
            {
                ret(i, j) += matA(i, k) * matB(k, j);
            }
        }
    }
    return ret;
}

template <typename MatAType, typename MatBType, int MatARows, int MatACols, int MatBCols, typename DType>
MatrixBase<MatAType, MatARows, MatACols, DType> &operator*=(MatrixBase<MatAType, MatARows, MatACols, DType> &matA,
                                                            const MatrixBase<MatBType, MatACols, MatBCols, DType> &matB)
{
    matA = matA * matB;
    return matA;
}

template <typename MatAType, typename MatBType, int Rows, int Cols, typename DType>
MatrixBase<MatAType, Rows, Cols, DType> &operator+=(MatrixBase<MatAType, Rows, Cols, DType> &matA,
                                                    const MatrixBase<MatBType, Rows, Cols, DType> &matB)
{
    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            matA(i, j) += matB(i, j);
        }
    }

    return matA;
}

template <typename MatAType, typename MatBType, int Rows, int Cols, typename DType>
MatrixBase<MatAType, Rows, Cols, DType> &operator-=(MatrixBase<MatAType, Rows, Cols, DType> &matA,
                                                    const MatrixBase<MatBType, Rows, Cols, DType> &matB)
{
    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            matA(i, j) -= matB(i, j);
        }
    }

    return matA;
}

template <typename MatType, int Rows, int Cols, typename DType>
MatrixBase<MatType, Rows, Cols, DType> &operator*=(MatrixBase<MatType, Rows, Cols, DType> &mat, const DType k)
{
    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            mat(i, j) *= k;
        }
    }
    return mat;
}

template <typename MatType, int Rows, int Cols, typename DType>
MatrixBase<MatType, Rows, Cols, DType> &operator/=(MatrixBase<MatType, Rows, Cols, DType> &mat, const DType k)
{
    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            mat(i, j) /= k;
        }
    }
    return mat;
}
template <typename MatType, int Rows, int Cols, typename DType>
MatrixBase<MatType, Rows, Cols, DType> &operator+=(MatrixBase<MatType, Rows, Cols, DType> &mat, const DType k)
{
    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            mat(i, j) += k;
        }
    }
    return mat;
}
template <typename MatType, int Rows, int Cols, typename DType>
MatrixBase<MatType, Rows, Cols, DType> &operator-=(MatrixBase<MatType, Rows, Cols, DType> &mat, const DType k)
{
    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            mat(i, j) -= k;
        }
    }
    return mat;
}

template <typename MatAType, typename MatBType, int Rows, int Cols, typename DType>
Matrix<Rows, Cols, DType> operator+(const MatrixBase<MatAType, Rows, Cols, DType> &matA,
                                    const MatrixBase<MatBType, Rows, Cols, DType> &matB)
{
    Matrix<Rows, Cols, DType> ret = matA;
    ret += matB;
    return ret;
}

template <typename MatAType, typename MatBType, int Rows, int Cols, typename DType>
Matrix<Rows, Cols, DType> operator-(const MatrixBase<MatAType, Rows, Cols, DType> &matA,
                                    const MatrixBase<MatBType, Rows, Cols, DType> &matB)
{
    Matrix<Rows, Cols, DType> ret = matA;
    ret -= matB;
    return ret;
}

template <typename MatType, int Rows, int Cols, typename DType>
Matrix<Rows, Cols, DType> operator+(const MatrixBase<MatType, Rows, Cols, DType> &mat, const DType k)
{
    Matrix<Rows, Cols, DType> ret = mat;
    ret += k;
    return ret;
}

template <typename MatType, int Rows, int Cols, typename DType>
Matrix<Rows, Cols, DType> operator-(const MatrixBase<MatType, Rows, Cols, DType> &mat, const DType k)
{
    Matrix<Rows, Cols, DType> ret = mat;
    ret -= k;
    return ret;
}

template <typename MatType, int Rows, int Cols, typename DType>
Matrix<Rows, Cols, DType> operator*(const MatrixBase<MatType, Rows, Cols, DType> &mat, const DType k)
{
    Matrix<Rows, Cols, DType> ret = mat;
    ret *= k;
    return ret;
}

template <typename MatType, int Rows, int Cols, typename DType>
Matrix<Rows, Cols, DType> operator/(const MatrixBase<MatType, Rows, Cols, DType> &mat, const DType k)
{
    Matrix<Rows, Cols, DType> ret = mat;
    ret /= k;
    return ret;
}

template <typename DerivedType, typename OperandType, int Rows, typename DType>
HorizontalConcat<DerivedType, OperandType> operator||(
    const MatrixBase<DerivedType, Rows, DerivedType::Cols, DType> &left,
    const MatrixBase<OperandType, Rows, OperandType::Cols, DType> &right)
{
    return HorizontalConcat<DerivedType, OperandType>(static_cast<const DerivedType &>(left),
                                                      static_cast<const OperandType &>(right));
}

template <typename DerivedType, typename OperandType, int Cols, typename DType>
VerticalConcat<DerivedType, OperandType> operator&&(
    const MatrixBase<DerivedType, DerivedType::Rows, Cols, DType> &top,
    const MatrixBase<OperandType, OperandType::Rows, Cols, DType> &bottom)

{
    return VerticalConcat<DerivedType, OperandType>(static_cast<const DerivedType &>(top),
                                                    static_cast<const OperandType &>(bottom));
}

template <typename MatAType, typename MatBType, int Rows, int Cols, typename DType>
Matrix<Rows, Cols, bool> operator==(const MatrixBase<MatAType, Rows, Cols, DType> &matA,
                                    const MatrixBase<MatBType, Rows, Cols, DType> &matB)
{
    Matrix<Rows, Cols, bool> ret;

    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            ret(i, j) = matA(i, j) == matB(i, j);
        }
    }
    return ret;
}

template <typename MatAType, typename MatBType, int Rows, int Cols>
Matrix<Rows, Cols, bool> operator&(const MatrixBase<MatAType, Rows, Cols, bool> &matA,
                                   const MatrixBase<MatBType, Rows, Cols, bool> &matB)
{
    Matrix<Rows, Cols, bool> ret;

    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            ret(i, j) = matA(i, j) & matB(i, j);
        }
    }
    return ret;
}

template <typename MatAType, typename MatBType, int Rows, int Cols>
Matrix<Rows, Cols, bool> operator|(const MatrixBase<MatAType, Rows, Cols, bool> &matA,
                                   const MatrixBase<MatBType, Rows, Cols, bool> &matB)
{
    Matrix<Rows, Cols, bool> ret;

    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            ret(i, j) = matA(i, j) | matB(i, j);
        }
    }
    return ret;
}

template <typename DerivedType>
Matrix<DerivedType::Rows, DerivedType::Cols, bool> operator!(
    const MatrixBase<DerivedType, DerivedType::Rows, DerivedType::Cols, bool> &matA)
{
    Matrix<DerivedType::Rows, DerivedType::Cols, bool> ret;

    for (int i = 0; i < DerivedType::Rows; ++i)
    {
        for (int j = 0; j < DerivedType::Cols; ++j)
        {
            ret(i, j) = !matA(i, j);
        }
    }
    return ret;
}

template <typename MatAType, typename MatBType, int Rows, int Cols, typename DType>
Matrix<Rows, Cols, bool> operator>(const MatrixBase<MatAType, Rows, Cols, DType> &matA,
                                   const MatrixBase<MatBType, Rows, Cols, DType> &matB)
{
    Matrix<Rows, Cols, bool> ret;

    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            ret(i, j) = matA(i, j) > matB(i, j);
        }
    }
    return ret;
}

template <typename MatAType, typename MatBType, int Rows, int Cols, typename DType>
Matrix<Rows, Cols, bool> operator<(const MatrixBase<MatAType, Rows, Cols, DType> &matA,
                                   const MatrixBase<MatBType, Rows, Cols, DType> &matB)
{
    Matrix<Rows, Cols, bool> ret;

    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            ret(i, j) = matA(i, j) < matB(i, j);
        }
    }
    return ret;
}

template <typename MatAType, typename MatBType, int Rows, int Cols, typename DType>
Matrix<Rows, Cols, bool> operator<=(const MatrixBase<MatAType, Rows, Cols, DType> &matA,
                                    const MatrixBase<MatBType, Rows, Cols, DType> &matB)
{
    Matrix<Rows, Cols, bool> ret;

    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            ret(i, j) = matA(i, j) <= matB(i, j);
        }
    }
    return ret;
}

template <typename MatAType, typename MatBType, int Rows, int Cols, typename DType>
Matrix<Rows, Cols, bool> operator>=(const MatrixBase<MatAType, Rows, Cols, DType> &matA,
                                    const MatrixBase<MatBType, Rows, Cols, DType> &matB)
{
    Matrix<Rows, Cols, bool> ret;

    for (int i = 0; i < Rows; ++i)
    {
        for (int j = 0; j < Cols; ++j)
        {
            ret(i, j) = matA(i, j) >= matB(i, j);
        }
    }
    return ret;
}

template <typename DerivedType>
Matrix<DerivedType::Rows, DerivedType::Cols, bool> operator>(const DownCast<DerivedType> &mat,
                                                             const typename DerivedType::DType k)
{
    Matrix<DerivedType::Rows, DerivedType::Cols, bool> ret;

    for (int i = 0; i < DerivedType::Rows; ++i)
    {
        for (int j = 0; j < DerivedType::Cols; ++j)
        {
            ret(i, j) = mat(i, j) > k;
        }
    }
    return ret;
}

template <typename DerivedType>
Matrix<DerivedType::Rows, DerivedType::Cols, bool> operator>=(const DownCast<DerivedType> &mat,
                                                              const typename DerivedType::DType k)
{
    Matrix<DerivedType::Rows, DerivedType::Cols, bool> ret;

    for (int i = 0; i < DerivedType::Rows; ++i)
    {
        for (int j = 0; j < DerivedType::Cols; ++j)
        {
            ret(i, j) = mat(i, j) >= k;
        }
    }
    return ret;
}

template <typename DerivedType>
Matrix<DerivedType::Rows, DerivedType::Cols, bool> operator<(const DownCast<DerivedType> &mat,
                                                             const typename DerivedType::DType k)
{
    Matrix<DerivedType::Rows, DerivedType::Cols, bool> ret;

    for (int i = 0; i < DerivedType::Rows; ++i)
    {
        for (int j = 0; j < DerivedType::Cols; ++j)
        {
            ret(i, j) = mat(i, j) < k;
        }
    }
    return ret;
}

template <typename DerivedType>
Matrix<DerivedType::Rows, DerivedType::Cols, bool> operator<=(const DownCast<DerivedType> &mat,
                                                              const typename DerivedType::DType k)
{
    Matrix<DerivedType::Rows, DerivedType::Cols, bool> ret;

    for (int i = 0; i < DerivedType::Rows; ++i)
    {
        for (int j = 0; j < DerivedType::Cols; ++j)
        {
            ret(i, j) = mat(i, j) <= k;
        }
    }
    return ret;
}

template <typename DerivedType>
bool Any(const MatrixBase<DerivedType, DerivedType::Rows, DerivedType::Cols, bool> &matA)
{
    for (int i = 0; i < DerivedType::Rows; ++i)
    {
        for (int j = 0; j < DerivedType::Cols; ++j)
        {
            if (matA(i, j))
            {
                return true;
            }
        }
    }

    return false;
}

template <typename DerivedType>
bool All(const MatrixBase<DerivedType, DerivedType::Rows, DerivedType::Cols, bool> &matA)
{
    for (int i = 0; i < DerivedType::Rows; ++i)
    {
        for (int j = 0; j < DerivedType::Cols; ++j)
        {
            if (!matA(i, j))
            {
                return false;
            }
        }
    }

    return true;
}

inline Print &operator<<(Print &strm, const int obj)
{
    strm.print(obj);
    return strm;
}

inline Print &operator<<(Print &strm, const float obj)
{
    strm.print(obj);
    return strm;
}

inline Print &operator<<(Print &strm, const char *obj)
{
    strm.print(obj);
    return strm;
}

inline Print &operator<<(Print &strm, const char obj)
{
    strm.print(obj);
    return strm;
}

// Stream inserter operator for printing to strings or the serial port
template <typename DerivedType, int Rows, int Cols, typename DType>
Print &operator<<(Print &strm, const MatrixBase<DerivedType, Rows, Cols, DType> &mat)
{
    strm << '[';

    for (int i = 0; i < Rows; i++)
    {
        strm << '[';

        for (int j = 0; j < Cols; j++) strm << mat(i, j) << ((j == Cols - 1) ? ']' : ',');

        strm << (i == Rows - 1 ? ']' : ',');
    }
    return strm;
}

}  // namespace BLA
