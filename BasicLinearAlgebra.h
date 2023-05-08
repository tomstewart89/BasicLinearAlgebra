#pragma once

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ElementStorage.h"

namespace BLA
{

template <typename DerivedType, int rows, int cols, typename d_type>
struct MatrixBase
{
   public:
    constexpr static int Rows = rows;
    constexpr static int Cols = cols;
    using DType = d_type;

    DType &operator()(int i, int j = 0)
    {
        assert(i < Rows && j < Cols);
        return static_cast<DerivedType *>(this)->operator()(i, j);
    }

    DType operator()(int i, int j = 0) const
    {
        assert(i < Rows && j < Cols);
        return static_cast<const DerivedType *>(this)->operator()(i, j);
    }

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

    template <int SubRows, int SubCols>
    MatrixRef<DerivedType, SubRows, SubCols> Submatrix(int row_start, int col_start)
    {
        return MatrixRef<DerivedType, SubRows, SubCols>(static_cast<DerivedType &>(*this), row_start, col_start);
    }

    template <int SubRows, int SubCols>
    MatrixRef<const DerivedType, SubRows, SubCols> Submatrix(int row_start, int col_start) const
    {
        return MatrixRef<const DerivedType, SubRows, SubCols>(static_cast<const DerivedType &>(*this), row_start,
                                                              col_start);
    }

    MatrixRef<DerivedType, 1, Cols> Row(int row_start)
    {
        return MatrixRef<DerivedType, 1, Cols>(static_cast<DerivedType &>(*this), row_start, 0);
    }

    MatrixRef<const DerivedType, 1, Cols> Row(int row_start) const
    {
        return MatrixRef<const DerivedType, 1, Cols>(static_cast<const DerivedType &>(*this), row_start, 0);
    }

    MatrixRef<DerivedType, Rows, 1> Column(int col_start)
    {
        return MatrixRef<DerivedType, Rows, 1>(static_cast<DerivedType &>(*this), 0, col_start);
    }

    MatrixRef<const DerivedType, Rows, 1> Column(int col_start) const
    {
        return MatrixRef<const DerivedType, Rows, 1>(static_cast<const DerivedType &>(*this), 0, col_start);
    }

    MatrixTranspose<DerivedType> operator~() { return MatrixTranspose<DerivedType>(static_cast<DerivedType &>(*this)); }

    MatrixTranspose<const DerivedType> operator~() const
    {
        return MatrixTranspose<const DerivedType>(static_cast<const DerivedType &>(*this));
    }

    template <typename OperandType, int OperandCols>
    HorizontalConcat<DerivedType, OperandType> operator||(
        const MatrixBase<OperandType, Rows, OperandCols, DType> &obj) const
    {
        return HorizontalConcat<DerivedType, OperandType>(static_cast<const DerivedType &>(*this),
                                                          static_cast<const OperandType &>(obj));
    }

    template <typename OperandType, int OperandRows>
    VerticalConcat<DerivedType, OperandType> operator&&(
        const MatrixBase<OperandType, OperandRows, Cols, DType> &obj) const
    {
        return VerticalConcat<DerivedType, OperandType>(static_cast<const DerivedType &>(*this),
                                                        static_cast<const OperandType &>(obj));
    }
};

}  // namespace BLA

#include "impl/BasicLinearAlgebra.h"
// #include "impl/NotSoBasicLinearAlgebra.h"
