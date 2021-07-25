#include <gtest/gtest.h>
#include "../BasicLinearAlgebra.h"

using namespace BLA;

TEST(LinearAlgebra, LUDecomposition)
{
    Matrix<3, 3> A{2, 5, 8, 0, 8, 6, 6, 7, 5};

    auto LU = A;
    ArrayMatrix<3, 1, int> P;

    LUDecompose(LU, P);
    Matrix<3, 3> L, U;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            // Fill in the diagonal
            if (i == j)
            {
                L(i, j) = 1.0;
                U(i, j) = LU(i, j);
            }
            // Upper triangle
            else if (i < j)
            {
                L(i, j) = 0.0;
                U(i, j) = LU(i, j);
            }
            // Lower triangle
            else
            {
                L(i, j) = LU(i, j);
                U(i, j) = 0.0;
            }
        }
    }

    auto A_reconstructed = L * U;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_FLOAT_EQ(A_reconstructed(P(i), j), A(i, j));
        }
    }
}

TEST(LinearAlgebra, LUSolution)
{
    Matrix<3, 3> A{2, 5, 8, 0, 8, 6, 6, 7, 5};
    Matrix<3, 1> b{10, 11, 12};
    Matrix<3, 1> x;
    Matrix<3, 1> x_expected = {0.41826923, 0.97115385, 0.53846154};

    ArrayMatrix<3, 1, int> P;

    LUDecompose(A, P);

    LUSolve(A, P, b, x);

    for (int i = 0; i < 3; ++i)
    {
        EXPECT_FLOAT_EQ(x_expected(i), x(i));
    }
}

TEST(LinearAlgebra, Inversion)
{
    Matrix<2, 2> A{1, 2, 3, 4};
    Matrix<2, 2> A_inv_numpy{-2, 1, 1.5, -0.5};

    auto A_inv = A;
    Invert(A_inv);

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            EXPECT_FLOAT_EQ(A_inv_numpy(i, j), A_inv(i, j));
        }
    }
}

TEST(Arithmetic, Determinant)
{
    Matrix<3, 3> B = {0.13397823, 0.21299481, 0.34636886, 0.75536007, 0.98363649, 0.59470267, 0.71449196, 0.08585212, 0.55099433};

    float det_numpy = -0.15333813;

    EXPECT_FLOAT_EQ(Determinant(B), det_numpy);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
