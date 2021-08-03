#include "../BasicLinearAlgebra.h"
#include <gtest/gtest.h>

using namespace BLA;

TEST(LinearAlgebra, LUDecomposition)
{
    BLA::Matrix<3, 3> A = {9.79, 9.33, 11.62,
                           7.77, 14.77, 14.12,
                           11.33, 15.72, 12.12};

    auto A_orig = A;

    auto decomp = LUDecompose(A);

    PermutationMatrix<3> P(decomp.permutation);
    LowerTriangularDiagonalOnesMatrix<3, 3, Array<3, 3, float>> L(decomp.lower);
    UpperTriangularMatrix<3, 3, Array<3, 3, float>> U(decomp.upper);

    auto A_reconstructed = P * L * U;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
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

    for (int i = 0; i < 3; ++i)
    {
        EXPECT_FLOAT_EQ(x_expected(i), x(i));
    }
}

TEST(LinearAlgebra, Inversion)
{
    BLA::Matrix<3, 3> A = {9.79, 9.33, 11.62,
                           7.77, 14.77, 14.12,
                           11.33, 15.72, 12.12};

    auto A_inv = A;
    Invert(A_inv);

    auto I = A_inv * A;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
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
}

TEST(Arithmetic, Determinant)
{
    Matrix<3, 3> B = {0.13397823, 0.21299481, 0.34636886,
                      0.75536007, 0.98363649, 0.59470267,
                      0.71449196, 0.08585212, 0.55099433};

    float det_numpy = -0.15333813;

    EXPECT_FLOAT_EQ(Determinant(B), det_numpy);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
