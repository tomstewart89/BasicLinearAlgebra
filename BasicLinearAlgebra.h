#pragma once

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Arduino.h"
#include "Printable.h"

#include "ElementStorage.h"

namespace BLA
{

template <typename DerivedType, int rows, int cols, typename d_type>
struct MatrixBase : public Printable
{
   public:
    constexpr static int Rows = rows;
    constexpr static int Cols = cols;
    using DType = d_type;

    DType &operator()(int i, int j = 0) { return static_cast<DerivedType *>(this)->operator()(i, j); }

    DType operator()(int i, int j = 0) const { return static_cast<const DerivedType *>(this)->operator()(i, j); }

    MatrixBase() = default;

    template <typename MatType>
    MatrixBase(const MatrixBase<MatType, Rows, Cols, DType> &mat)
    {
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                static_cast<DerivedType &>(*this)(i, j) = mat(i, j);
            }
        }
    }

    MatrixBase &operator=(const MatrixBase &mat)
    {
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                static_cast<DerivedType &>(*this)(i, j) = mat(i, j);
            }
        }

        return static_cast<DerivedType &>(*this);
    }

    template <typename MatType>
    MatrixBase &operator=(const MatrixBase<MatType, Rows, Cols, DType> &mat)
    {
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                static_cast<DerivedType &>(*this)(i, j) = mat(i, j);
            }
        }

        return static_cast<DerivedType &>(*this);
    }

    DerivedType &operator=(DType elem)
    {
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                static_cast<DerivedType &>(*this)(i, j) = elem;
            }
        }

        return static_cast<DerivedType &>(*this);
    }

    void Fill(const DType &val) { *this = val; }

    template <typename DestType>
    Matrix<Rows, Cols, DestType> Cast()
    {
        Matrix<Rows, Cols, DestType> ret;

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                ret(i, j) = (DestType)(*this)(i, j);
            }
        }

        return ret;
    }

    template <int SubRows, int SubCols>
    RefMatrix<DerivedType, SubRows, SubCols> Submatrix(int row_start, int col_start)
    {
        return RefMatrix<DerivedType, SubRows, SubCols>(static_cast<DerivedType &>(*this), row_start, col_start);
    }

    template <int SubRows, int SubCols>
    RefMatrix<const DerivedType, SubRows, SubCols> Submatrix(int row_start, int col_start) const
    {
        return RefMatrix<const DerivedType, SubRows, SubCols>(static_cast<const DerivedType &>(*this), row_start,
                                                              col_start);
    }

    RefMatrix<DerivedType, 1, Cols> Row(int row_start)
    {
        return RefMatrix<DerivedType, 1, Cols>(static_cast<DerivedType &>(*this), row_start, 0);
    }

    RefMatrix<const DerivedType, 1, Cols> Row(int row_start) const
    {
        return RefMatrix<const DerivedType, 1, Cols>(static_cast<const DerivedType &>(*this), row_start, 0);
    }

    RefMatrix<DerivedType, Rows, 1> Column(int col_start)
    {
        return RefMatrix<DerivedType, Rows, 1>(static_cast<DerivedType &>(*this), 0, col_start);
    }

    RefMatrix<const DerivedType, Rows, 1> Column(int col_start) const
    {
        return RefMatrix<const DerivedType, Rows, 1>(static_cast<const DerivedType &>(*this), 0, col_start);
    }

    MatrixTranspose<DerivedType> operator~() { return MatrixTranspose<DerivedType>(static_cast<DerivedType &>(*this)); }

    MatrixTranspose<const DerivedType> operator~() const
    {
        return MatrixTranspose<const DerivedType>(static_cast<const DerivedType &>(*this));
    }

    Matrix<Rows, Cols, DType> operator-() const
    {
        Matrix<Rows, Cols, DType> ret;

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                ret(i, j) = -(*this)(i, j);
            }
        }

        return ret;
    }

    size_t printTo(Print& p) const final
    {
        size_t n;
        n = p.print('[');

        for (int i = 0; i < Rows; i++)
        {
            n += p.print('[');

            for (int j = 0; j < Cols; j++)
            {
                n += p.print(static_cast<const DerivedType *>(this)->operator()(i, j));
                n += p.print((j == Cols - 1) ? ']' : ',');
            }

            n += p.print((i == Rows - 1) ? ']' : ',');
        }
        return n;
    }
};

template <typename DerivedType>
using DownCast = MatrixBase<DerivedType, DerivedType::Rows, DerivedType::Cols, typename DerivedType::DType>;

}  // namespace BLA

#include "impl/Types.h"
#include "impl/BasicLinearAlgebra.h"
#include "impl/NotSoBasicLinearAlgebra.h"
