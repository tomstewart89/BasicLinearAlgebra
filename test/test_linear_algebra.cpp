#include <gtest/gtest.h>

#include "../BasicLinearAlgebra.h"

using namespace BLA;

TEST(LinearAlgebra, LUDecomposition)
{
    Matrix<7, 7> A = {16, 78, 50, 84, 70, 63, 2, 32, 33, 61, 40, 17, 96, 98, 50, 80, 78, 27, 86, 49, 57, 10, 42, 96, 44,
                      87, 60, 67, 16, 59, 53, 8, 64, 97, 41, 90, 56, 22, 48, 32, 12, 4,  45, 78, 43, 11, 7,  8,  12};

    auto A_orig = A;

    auto decomp = LUDecompose(A);

    EXPECT_FALSE(decomp.singular);

    auto A_reconstructed = decomp.P() * decomp.L() * decomp.U();

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
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
