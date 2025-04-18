#include <gtest/gtest.h>

#include "../BasicLinearAlgebra.h"

using namespace BLA;

TEST(LinearAlgebra, CrossProduct)
{
    // Axis 0
    const Matrix<3> A = {1, 2,  3};
    const Matrix<3> B = {5, 7, 11};
    const Matrix<3> AB_expected = {1, 4, -3};
    const auto AB_result = CrossProduct(A, B);

    for (int i = 0; i < AB_result.Rows; i++)
    {
        for (int j = 0; j < AB_result.Cols; j++)
        {
            EXPECT_FLOAT_EQ(AB_result(i, j), AB_expected(i, j));
        }
    }
    // The cross product will always be computed on the second axis.
    // This can be changed by transposing the matricies.
    const Matrix<3, 3> E = { 1,  2,  3,  5,  7, 11, 13, 17, 19};
    const Matrix<3, 3> F = {23, 29, 31, 37, 41, 43, 47, 53, 59};
    const Matrix<3, 3> EF_expected = {-25, 38, -17, -150, 192, -54, -4, 126, -110};
    const auto EF_result = CrossProduct(~E, ~F);

    for (int i = 0; i < EF_expected.Rows; i++)
    {
        for (int j = 0; j < EF_expected.Cols; j++)
        {
            EXPECT_FLOAT_EQ(EF_result(i, j), EF_expected(j, i));
        }
    }
}

TEST(LinearAlgebra, DotProduct)
{
    // Test with Dimension 3, int
    Matrix<3, 1, int> A = {2, 3, 4};
    Matrix<3, 1, int> B = {7, 8, 9};
    EXPECT_FLOAT_EQ(DotProduct(A, B), 74);

    // Test with Dimension 4, float
    Matrix<4> C = {2.5f, -3.4f, 4.3f, 5.2f};
    Matrix<4> D = {6.7f, 7.6f, -8.5f, 9.4f};
    EXPECT_NEAR(DotProduct(C, D), 3.23999f, 1e-5);
}

TEST(LinearAlgebra, LUDecomposition)
{
    Matrix<7, 7> A = {16, 78, 50, 84, 70, 63, 2, 32, 33, 61, 40, 17, 96, 98, 50, 80, 78, 27, 86, 49, 57, 10, 42, 96, 44,
                      87, 60, 67, 16, 59, 53, 8, 64, 97, 41, 90, 56, 22, 48, 32, 12, 4,  45, 78, 43, 11, 7,  8,  12};

    auto A_orig = A;

    auto decomp = LUDecompose(A);

    EXPECT_FALSE(decomp.singular);

    auto A_reconstructed = decomp.P * decomp.L * decomp.U;

    for (int i = 0; i < A.Rows; ++i)
    {
        for (int j = 0; j < A.Cols; ++j)
        {
            EXPECT_FLOAT_EQ(A_reconstructed(i, j), A_orig(i, j));
        }
    }
}

TEST(LinearAlgebra, LUSolution)
{
    Matrix<3, 3> A{2, 5, 8, 0, 8, 6, 6, 7, 5};
    Matrix<3, 1> b{10, 11, 12};
    Matrix<3, 1> x_expected = {0.41826923, 0.97115385, 0.53846154};

    auto decomp = LUDecompose(A);

    auto x = LUSolve(decomp, b);

    for (int i = 0; i < x_expected.Rows; ++i)
    {
        EXPECT_FLOAT_EQ(x_expected(i), x(i));
    }
}

TEST(LinearAlgebra, CholeskyDecomposition)
{
    // clang-format off

    // We could fill in this lower triangle but since A is required to be symmetric it can be (and is) inferred
    // from the upper triangle
    Matrix<4, 4> A = {0.60171582, -0.20854924,  0.52925771,  0.24206045,
                      0.0,         0.33012847, -0.28941531, -0.33854164,
                      0.0,         0.0,         3.54506632,  1.56758518,
                      0.0,         0.0,         0.0,         1.75291733};
    // clang-format on

    auto A_orig = A;

    auto chol = CholeskyDecompose(A);

    EXPECT_TRUE(chol.positive_definite);

    auto A_reconstructed = chol.L * ~chol.L;

    // Compare the recontruction to the upper triangle of A (the lower triangle will be overwritten by decompose)
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            if (i <= j)
            {
                EXPECT_FLOAT_EQ(A_reconstructed(i, j), A_orig(i, j));
            }
        }
    }
}

TEST(LinearAlgebra, CholeskySolution)
{
    Matrix<5, 5> A = {0.78183123,  0.08385324,  0.37172332,  -0.72518705, -1.11317593, 0.08385324, 0.56011595,
                      0.19965695,  -0.17488402, -0.12703805, 0.37172332,  0.19965695,  0.52769031, -0.19284881,
                      -0.45321194, -0.72518705, -0.17488402, -0.19284881, 2.19127456,  2.13045896, -1.11317593,
                      -0.12703805, -0.45321194, 2.13045896,  3.50184434};

    Matrix<5> b = {1.0, 2.0, 3.0, 4.0, 5.0};
    Matrix<5> x_expected = {3.15866835, 2.12529984, 5.23818026, 0.98626514, 2.58690994};

    Matrix<5, 5> A_copy = A;

    auto chol = CholeskyDecompose(A);

    auto x = CholeskySolve(chol, b);

    for (int i = 0; i < 5; ++i)
    {
        EXPECT_NEAR(x_expected(i), x(i), 1e-5);
    }
}

TEST(LinearAlgebra, Inversion)
{
    BLA::Matrix<3, 3> A = {9.79, 9.33, 11.62, 7.77, 14.77, 14.12, 11.33, 15.72, 12.12};

    auto A_inv = A;
    Invert(A_inv);

    auto I = A_inv * A;

    for (int i = 0; i < A.Rows; ++i)
    {
        for (int j = 0; j < A.Cols; ++j)
        {
            if (i == j)
            {
                EXPECT_NEAR(I(i, j), 1.0, 1e-5);
            }
            else
            {
                EXPECT_NEAR(I(i, j), 0.0, 1e-5);
            }
        }
    }

    BLA::Matrix<3, 3> singular = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0};

    EXPECT_FALSE(Invert(singular));
}

TEST(LinearAlgebra, DoublePrecisionInverse)
{
    Matrix<6, 6, double> A = {1. / 48.,  0,          0,         0, 0, 0, 0, 1. / 48.,  0,        0,       0, 0, 0,
                              -1. / 48., 1. / 48.,   0,         0, 0, 0, 0, 0,         1. / 24., 0,       0, 0, 0,
                              0,         -1. / 28.8, 1. / 28.8, 0, 0, 0, 0, -1. / 12., 1. / 24., 1. / 24.};

    auto A_inv = Inverse(A * 1.8);

    EXPECT_DOUBLE_EQ(A_inv(0, 0), 80.0 / 3.0);
    EXPECT_DOUBLE_EQ(A_inv(5, 5), 40.0 / 3.0);
}

TEST(Arithmetic, Determinant)
{
    Matrix<6, 6> B = {0.05508292, 0.82393504, 0.34938018, 0.63818054, 0.18291131, 0.1986636,  0.56799604, 0.81077491,
                      0.71472733, 0.68527613, 0.72759853, 0.25983183, 0.99035713, 0.76096889, 0.26130098, 0.16855372,
                      0.0253581,  0.47907605, 0.58735833, 0.0913456,  0.03221577, 0.5210331,  0.61583369, 0.33233299,
                      0.20578816, 0.356537,   0.70661899, 0.6569476,  0.90074756, 0.59771572, 0.20054716, 0.41290408,
                      0.70679818, 0.321249,   0.81886099, 0.77819212};

    float det_numpy = -0.03919640039505248;

    EXPECT_FLOAT_EQ(Determinant(B), det_numpy);

    BLA::Matrix<3, 3> singular = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0};

    EXPECT_FLOAT_EQ(Determinant(singular), 0.0);

    BLA::Matrix<4, 4, int16_t> C = {8, 5, 5, 8, 3, 1, 3, 2, 1, 1, 3, 0, 3, 3, 5, 9};
    int16_t det_C_expected = -140;

    EXPECT_EQ(Determinant(C), det_C_expected);

    BLA::Matrix<3, 3, int> singular_int = {1, 0, 0, 1, 0, 0, 1, 0, 0};

    EXPECT_EQ(Determinant(singular_int), 0);
}

template <typename SparseMatA, typename SparseMatB, int OutTableSize = 100>
SparseMatrix<SparseMatA::Rows, SparseMatB::Cols, typename SparseMatA::DType, OutTableSize> sparse_mul(
    const SparseMatA &A, const SparseMatB &B)
{
    static_assert(A.Cols == B.Rows, "Incompatible dimensions for sparse matrix multiplication");

    SparseMatrix<A.Rows, B.Cols, typename SparseMatA::DType, OutTableSize> out;

    for (int i = 0; i < SparseMatA::Size; ++i)
    {
        for (int j = 0; j < SparseMatB::Size; ++j)
        {
            const auto &elem_a = A.table[i];
            const auto &elem_b = B.table[j];

            if (elem_a.row >= 0 && elem_b.row >= 0)
            {
                out(elem_a.row, elem_b.col) += elem_a.val * elem_b.val;
            }
        }
    }

    return out;
}

TEST(Examples, SparseMatrix)
{
    SparseMatrix<1, 3000, float, 100> A;
    SparseMatrix<3000, 1, float, 100> B;

    A(0, 1000) = 5.0;
    B(1000, 0) = 5.0;

    auto C = sparse_mul(A, B);

    EXPECT_EQ(C(0, 0), 25.0);
}

class JacobianTestFunctor : public MatrixFunctor<2, 3, float>
{
    Matrix<3, 1, float> operator()(const Matrix<2, 1, float> &x) const override
    {
        Matrix<3> p1 = {cos(x(0)), sin(x(0)), x(0)};
        Matrix<3> p2 = {cos((x(0) + x(1))), sin((x(0) + x(1))), x(1)};

        return p1 + p2;
    }
};

TEST(Examples, NumericJacobian)
{
    // JacobianTestFunctor(x) = [cos(x1) + cos(x1 + x2), sin(x1) + sin(x1 + x2), x1 + x2]
    // jacobian =
    //       [[-sin(x1) - sin(x1 + x2), -sin(x1 + x2)]
    //        [cos(x1) + cos(x1 + x2) ,  cos(x1 + x2)]
    //        [      1                ,          1   ]]

    Matrix<2> x = {0.0f, 0.0f};
    JacobianTestFunctor functor;

    Matrix<3, 2> jacobian = Jacobian<2, 3>(JacobianTestFunctor(), x);

    EXPECT_NEAR(jacobian(0, 0), 0, 1e-4);
    EXPECT_NEAR(jacobian(0, 1), 0, 1e-4);

    EXPECT_NEAR(jacobian(1, 0), 2, 1e-4);
    EXPECT_NEAR(jacobian(1, 1), 1, 1e-4);

    EXPECT_NEAR(jacobian(2, 0), 1, 1e-4);
    EXPECT_NEAR(jacobian(2, 1), 1, 1e-4);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
