#include <BasicLinearAlgebra.h>

/*
 * This example sketch should show you everything you need to know in order to work with the Matrix library. It's a
 * faily long explantion but try to read through it carefully. Remember that everything is built around templates so the
 * compatibility of matrices used in all operations can be checked at compile time. For that reason though, if you try
 * to multiply matrices of the wrong dimensions etc you'll be met with some very cryptic compile errors so be careful!
 */

// All the functions in BasicLinearAlgebra are wrapped up inside the namespace BLA, so specify that we're using it like
// so:
using namespace BLA;

void setup()
{
    Serial.begin(115200);

    // Let's being by declaring a few matrices. Matrix is a template class so you'll have to specfify the dimensions as
    // <row,column> after Matrix like so:
    BLA::Matrix<3, 3> A;

    // NOTE: Turns out that one of the libraries included when using an Arduino DUE also declares a class called Matrix,
    // so to resolve the ambiguity you'll need to explicitly identify this library's Matrix class when declaring them by
    // prepending the declaration with BLA:: If you're using other boards then just writing using namespace BLA; will
    // suffice.

    // The cols parameters has a default of 1 so to declare a vector of length 3 you can simply write:
    BLA::Matrix<3> v;

    // Just like any other variable, Matrices should be initialised to some value before you use them. To set every
    // element of the Matrix to a certain value you can use the Fill function like so:
    v.Fill(0);

    // The matrix elements can be accessed individually using the brackets operator like so:
    v(2, 0) = 5.3;
    v(1) = 43.67;  // you can also just write v(2) = 5.3; since v has only one column

    // Or you can set the entire matrix like so:
    A = {3.25, 5.67, 8.67, 4.55, 7.23, 9.00, 2.35, 5.73, 10.56};

    // You can also set the entire array on construction like so:
    BLA::Matrix<3, 3> B = {6.54, 3.66, 2.95, 3.22, 7.54, 5.12, 8.98, 9.99, 1.56};

    // Ask the matrix how many rows and columns it has:
    v.Cols;
    v.Rows;

    // Now you can do some matrix math! The Matrix class supports addition and subtraction between matrices of the same
    // size:
    BLA::Matrix<3, 3> C = A + B;

    // You can't do addition between matrices of a different size. Since the class is built around templates, the
    // compiler will tell you if you've made a mistake. Try uncommenting the next line to see what I mean A + v;

    // You can use the Matrix class's unary operators, like so:
    C -= B;

    // Or like so:
    B += A;

    // As well as addition and subtraction, we can also do matrix multiplication. Note that again, the matrices,
    // including the one in which the result will stored must have the appropriate dimensions. The compiler will let you
    // know if they aren't
    BLA::Matrix<3, 1> D = A * v;

    // As well as algebra, Matrix supports a few other matrix related operations including transposition:
    BLA::Matrix<1, 3> D_T = ~D;

    // And concatenation, both horizontally...
    BLA::Matrix<3, 6> AleftOfB = A || B;

    // And vertically
    BLA::Matrix<6, 3> AonTopOfB = A && B;

    // An inverse of a matrix can also be calculated for square matrices via the Invert function:
    BLA::Matrix<3, 3> C_inv = C;
    bool is_nonsingular = Invert(C_inv);

    // If the matrix is singular, the inversion won't work. In those cases Invert will return false

    // If you want to print out the value of any element in the array you can do that like so:
    Serial << "v(1): " << v(1) << '\n';

    // Alternatively, you can write the whole matrix to Serial using a C++ style inserter, like so:
    Serial << "B: " << B << '\n';

    // You can even write some quite complex compound statements very succinctly. For example:
    // Serial << "identity matrix: " << C << C_inv << '\n';
    Serial << "identity matrix: " << AleftOfB * AonTopOfB - (A * A + B * B) + C * C_inv;

    // Or as a more real life example, here's how you might calculate an updated state estimate for a third order state
    // space model:
    BLA::Matrix<3> x;
    BLA::Matrix<2> u;
    BLA::Matrix<3, 2> G;
    BLA::Matrix<3, 3> F;
    float dt;
    x += (F * x + G * u) * dt;

    // Enjoy!
}

void loop() {}
