#include <gtest/gtest.h>

#include <BasicLinearAlgebra.h>

using namespace BLA;

// namespace References
// {
// #include "../examples/References/References.ino"
// }

// TEST(Examples, References)
// {
//     References::setup();

//     EXPECT_STREQ(Serial.buf.str().c_str(), "still ones: [[1.00,1.00,1.00,1.00],[1.00,1.00,1.00,1.00],[1.00,1.00,1.00,1.00],[1.00,1.00,1.00,1.00]]\n"
//                                            "scaled rows: [[1.00,1.00,1.00,1.00],[2.00,2.00,2.00,2.00],[3.00,3.00,3.00,3.00],[4.00,4.00,4.00,4.00]]");
// }

namespace SolveLinearEquations
{
#include "../examples/SolveLinearEquations/SolveLinearEquations.ino"
}

TEST(Examples, SolveLinearEquations)
{
    SolveLinearEquations::setup();

    EXPECT_STREQ(Serial.buf.str().c_str(), "");
}

namespace Tensor
{
#include "../examples/Tensor/Tensor.ino"
}

TEST(Examples, Tensor)
{
    Tensor::setup();

    EXPECT_STREQ(Serial.buf.str().c_str(), "Hyper B: [[[[6.00,10.00],[10.00,18.00]],[[10.00,14.00],[18.00,26.00]]],[[[10.00,18.00],[14.00,26.00]],[[18.00,26.00],[26.00,38.00]]]]");
}

namespace CustomMatrix
{
#include "../examples/CustomMatrix/CustomMatrix.ino"
}

TEST(Examples, CustomMatrix)
{
    CustomMatrix::setup();

    EXPECT_STREQ(Serial.buf.str().c_str(), "still ones: [[1.00,1.00,1.00,1.00],[1.00,1.00,1.00,1.00],[1.00,1.00,1.00,1.00],[1.00,1.00,1.00,1.00]]\n"
                                           "scaled rows: [[1.00,1.00,1.00,1.00],[2.00,2.00,2.00,2.00],[3.00,3.00,3.00,3.00],[4.00,4.00,4.00,4.00]]");
}

namespace HowToUse
{
#include "../examples/HowToUse/HowToUse.ino"
}

TEST(Examples, HowToUse)
{
    HowToUse::setup();

    EXPECT_STREQ(Serial.buf.str().c_str(), "v(1): 43.67\n"
                                           "B: [[9.79,9.33,11.62],[7.77,14.77,14.12],[11.33,15.72,12.12]]\n"
                                           "identity matrix: [[1.00,-0.00,-0.00],[0.00,1.00,-0.00],[0.00,0.00,1.00]]");
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
