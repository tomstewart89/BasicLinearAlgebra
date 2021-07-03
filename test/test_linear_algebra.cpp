#include <gtest/gtest.h>
#include "../BasicLinearAlgebra.h"

using namespace BLA;

// TEST(LinearAlgebra, Inversion)
// {
//     Matrix<3, 3> A{0.18246341, 0.94775365, 0.65449397, 0.25632926, 0.29677899, 0.95827993, 0.75542384, 0.40581397, 0.23745302};
//     Matrix<3, 3> A_inv_numpy{-4.58985214, 2.69814355, 4.14465124, 2.67688548, -0.44356033, -3.19670431, 0.14779867, -0.77569282, 5.92154198};

//     auto A_inv = Invert(A);

//     for (int i = 0; i < 3; ++i)
//     {
//         for (int j = 0; j < 3; ++j)
//         {
//             EXPECT_FLOAT_EQ(A_inv_numpy(i, j), A_inv(i, j));
//         }
//     }
// }

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
