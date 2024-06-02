#include <gtest/gtest.h>

#include "../BasicLinearAlgebra.h"

using namespace BLA;

TEST(Arithmetic, BraceInitialisation)
{
    Matrix<2, 2> B{0.0, 45.34, 32.98, 1456.1222};

    EXPECT_FLOAT_EQ(B(0, 0), 0);
    EXPECT_FLOAT_EQ(B(0, 1), 45.34);
    EXPECT_FLOAT_EQ(B(1, 0), 32.98);
    EXPECT_FLOAT_EQ(B(1, 1), 1456.1222);

    Matrix<2, 2> C{1.54, 5.98};

    EXPECT_FLOAT_EQ(C(0, 0), 1.54);
    EXPECT_FLOAT_EQ(C(0, 1), 5.98);
    EXPECT_FLOAT_EQ(C(1, 0), 0.0);
    EXPECT_FLOAT_EQ(C(1, 1), 0.0);
}

TEST(Arithmetic, Fill)
{
    Matrix<2, 2> A;
    A.Fill(0.0f);

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            EXPECT_FLOAT_EQ(A(i, j), 0);
        }
    }
}

TEST(Arithmetic, OnesTest)
{
    Matrix<2, 2> A = Zeros<2, 2>();
    Matrix<2, 2> B = Ones<2, 2>();

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            EXPECT_FLOAT_EQ(A(i, j), 0.0f);
            EXPECT_FLOAT_EQ(B(i, j), 1.0f);
        }
    }
}

TEST(Arithmetic, EyeTest)
{
    auto I = BLA::Eye<2, 2>();
    auto Z = BLA::Zeros<2, 2>();
    auto R = I + Z;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            EXPECT_FLOAT_EQ(I(i, j), i == j ? 1.0f : 0.0f);
            EXPECT_FLOAT_EQ(Z(i, j), 0.0f);
            EXPECT_FLOAT_EQ(R(i, j), i == j ? 1.0f : 0.0f);
        }
    }
}

TEST(Arithmetic, AdditionSubtraction)
{
    Matrix<3, 3> A = {3.25, 5.67, 8.67, 4.55, 7.23, 9.00, 2.35, 5.73, 10.56};

    Matrix<3, 3> B = {6.54, 3.66, 2.95, 3.22, 7.54, 5.12, 8.98, 9.99, 1.56};

    auto C = A + B;
    auto D = A - B;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_FLOAT_EQ(C(i, j), A(i, j) + B(i, j));
            EXPECT_FLOAT_EQ(D(i, j), A(i, j) - B(i, j));
        }
    }

    C -= B;
    D += B;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_FLOAT_EQ(C(i, j), D(i, j));
        }
    }
}

TEST(Arithmetic, OtherDTypes)
{
    Matrix<3, 3, int> A = {3, 6, 5, 8, 34, 7, 3, 7, 9};

    Matrix<3, 3, bool> B = {true, false, true, true, false, false};

    auto C = A + 5;
    auto D = B * false;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            EXPECT_EQ(C(i, j), A(i, j) + 5);
            EXPECT_FALSE(D(i, j));
        }
    }
}

TEST(Arithmetic, ElementwiseOperations)
{
    Matrix<3, 3> A = {3.25f, 5.67f, 8.67f, 4.55f, 7.23f, 9.00f, 2.35f, 5.73f, 10.56f};

    auto C = A + 2.5f;
    auto D = A - 3.7f;
    auto E = A * 1.2f;
    auto F = A / 6.7f;
    auto G = -A;
    const auto H = 2.5f + A;
    const auto I = 3.7f - A;
    const auto J = 1.2f * A;
    const auto K = 6.7f / A;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            EXPECT_FLOAT_EQ(C(i, j), A(i, j) + 2.5);
            EXPECT_FLOAT_EQ(D(i, j), A(i, j) - 3.7);
            EXPECT_FLOAT_EQ(E(i, j), A(i, j) * 1.2);
            EXPECT_FLOAT_EQ(F(i, j), A(i, j) / 6.7);
            EXPECT_FLOAT_EQ(G(i, j), -A(i, j));
            EXPECT_FLOAT_EQ(H(i, j), 2.5 + A(i, j));
            EXPECT_FLOAT_EQ(I(i, j), 3.7 - A(i, j));
            EXPECT_FLOAT_EQ(J(i, j), 1.2 * A(i, j));
            EXPECT_FLOAT_EQ(K(i, j), 6.7 / A(i, j));
        }
    }

    C -= 2.5f;
    D += 3.7f;
    E /= 1.2f;
    F *= 6.7f;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            EXPECT_FLOAT_EQ(C(i, j), A(i, j));
            EXPECT_FLOAT_EQ(D(i, j), A(i, j));
            EXPECT_FLOAT_EQ(E(i, j), A(i, j));
            EXPECT_FLOAT_EQ(F(i, j), A(i, j));
        }
    }
}

TEST(Arithmetic, Multiplication)
{
    Matrix<3, 3> A = {3., 5., 8., 4., 7., 9., 2., 5.0, 10.};

    Matrix<3, 3> B = {6., 3., 2., 3., 7., 5., 8., 9., 1.};

    auto C = A * B;

    EXPECT_FLOAT_EQ(C(0, 0), 97.);
    EXPECT_FLOAT_EQ(C(0, 1), 116.);
    EXPECT_FLOAT_EQ(C(0, 2), 39.);
    EXPECT_FLOAT_EQ(C(1, 0), 117.);
    EXPECT_FLOAT_EQ(C(1, 1), 142.);
    EXPECT_FLOAT_EQ(C(1, 2), 52.);
    EXPECT_FLOAT_EQ(C(2, 0), 107);
    EXPECT_FLOAT_EQ(C(2, 1), 131.);
    EXPECT_FLOAT_EQ(C(2, 2), 39.);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_FLOAT_EQ(C(i, j), (A.Row(i) * B.Column(j))(0, 0));
        }
    }

    A *= B;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_FLOAT_EQ(A(i, j), C(i, j));
        }
    }
}

TEST(Arithmetic, Concatenation)
{
    Matrix<3, 3> A = {3.25, 5.67, 8.67, 4.55, 7.23, 9.00, 2.35, 5.73, 10.56};

    Matrix<3, 3> B = {6.54, 3.66, 2.95, 3.22, 7.54, 5.12, 8.98, 9.99, 1.56};

    auto AleftOfB = A || B;
    auto AonTopOfB = A && B;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_FLOAT_EQ(AleftOfB(i, j), A(i, j));
            EXPECT_FLOAT_EQ(AleftOfB(i, j + 3), B(i, j));
            EXPECT_FLOAT_EQ(AonTopOfB(i, j), A(i, j));
            EXPECT_FLOAT_EQ(AonTopOfB(i + 3, j), B(i, j));
        }
    }
}

TEST(Arithmetic, OuterProduct)
{
    Matrix<3> v = {1.0, 2.0, 3.0};

    Matrix<3, 3> A = v * ~v;

    for (int i = 0; i < A.Rows; ++i)
    {
        for (int j = 0; j < A.Cols; ++j)
        {
            EXPECT_FLOAT_EQ(A(i, j), v(i) * v(j));
        }
    }
}

TEST(Arithmetic, Reference)
{
    const Matrix<3, 3> A = {3.25, 5.67, 8.67, 4.55, 7.23, 9.00, 2.35, 5.73, 10.56};
    Matrix<3, 3> B = {2.25, 6.77, 9.67, 14.55, 0.23, 3.21, 5.67, 6.75, 11.56};

    const auto A_ref = A.Submatrix<2, 2>(0, 1);
    auto B_ref = B.Submatrix<2, 2>(0, 1);

    B_ref(0, 0) = A(0, 0) * A_ref(0, 0);

    EXPECT_FLOAT_EQ(B(0, 1), A(0, 0) * A(0, 1));

    B.Submatrix<3, 3>(0, 0) = A;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_FLOAT_EQ(B(i, j), A(i, j));
        }
    }
}

TEST(Arithmetic, Cast)
{
    Matrix<3, 3> A = {3.25, 5.67, 8.67, 4.55, 7.23, 9.00, 2.35, 5.73, 10.56};

    auto bool_A = A.Cast<bool>();

    EXPECT_TRUE(All(bool_A));

    Matrix<3, 3, double> A_double = A.Cast<double>();

    EXPECT_LT(Norm(A * A_double.Cast<float>() - A * A), 1e-5);
}

TEST(Arithmetic, LogicalOperators)
{
    Matrix<3, 3> A = {3.25, 5.67, 8.67, 4.55, 7.23, 9.00, 2.35, 5.73, 10.56};
    Matrix<3, 3> B = {3.25, 6.77, 9.67, 14.55, 0.23, 3.21, 5.67, 6.75, 11.56};

    auto A_less_than_three = A < 3.0;
    auto A_greater_than_or_equal_to_three_and_a_bit = A <= 3.25;
    auto A_greater_than_one_hundred = A > 100.0;
    auto A_greater_than_or_equal_to_ten_point_fiveish = A >= 10.56;
    auto A_less_than_B = A < B;
    auto A_less_than_or_equal_to_B = A <= B;
    auto A_greater_than_B = A > B;
    auto A_greater_than_or_equal_to_B = A >= B;
    auto A_equals_B = A == B;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_TRUE(A_less_than_three(i, j) == A(i, j) < 3.0);
            EXPECT_TRUE(A_greater_than_or_equal_to_three_and_a_bit(i, j) == A(i, j) <= 3.25);
            EXPECT_TRUE(A_greater_than_one_hundred(i, j) == A(i, j) > 100.0);
            EXPECT_TRUE(A_greater_than_or_equal_to_ten_point_fiveish(i, j) == A(i, j) >= 10.56);
            EXPECT_TRUE(A_less_than_B(i, j) == A(i, j) < B(i, j));
            EXPECT_TRUE(A_less_than_or_equal_to_B(i, j) == A(i, j) <= B(i, j));
            EXPECT_TRUE(A_greater_than_B(i, j) == A(i, j) > B(i, j));
            EXPECT_TRUE(A_greater_than_or_equal_to_B(i, j) == A(i, j) >= B(i, j));
        }
    }

    EXPECT_TRUE(!All(A_greater_than_one_hundred));
    EXPECT_FALSE(Any(A_greater_than_one_hundred));
    EXPECT_TRUE(Any(A_less_than_three));
    EXPECT_FALSE(All(A_less_than_three));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
