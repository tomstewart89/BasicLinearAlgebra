#include <Matrix.h>

/*
 * This example sketch should show you everything you need to know in order to work with the Matrix library. It's a faily long explantion
 * but try to read through it carefully. Remember that everything is built around templates so there's no heap memory useage and the compatibility
 * of matrices used in all operations is checked at compile time. For that reason though, if you try to multiply matrices of the wrong dimensions etc
 * you'll be met with some very cryptic compile errors so be careful!
 */

void setup() 
{
  Serial.begin(115200);

  // Let's being by declaring a few matrices. Matrix is a template class so you'll have to specfify the dimensions as well
  // as the underlying datatype as <row,column,type> after Matrix like so:
  Matrix<3,3,float> A;

  // The template parameters have defaults so if you're satisfied with a matrix of type 'float' with only 1 column you can simply write:
  Matrix<3> v;

  // The matrix elements can be accessed individually using the brackets operator like so:
  v(2,0) = 5.3; 
  v(1) = 43.67; // you can also just write v(2) = 5.3; since v has only one column

  // Note that the brackets operator will always return some value. If by mistake you ask for some element outside it's bounds
  // it'll return a dummy value and set a flag to indicate the violation, like so:
  if(A.invalidAccessFlag == false)
    Serial << "so far so good\n";
  
  A(1000,1000) = 65.9;

  if(A.invalidAccessFlag == true)
      Serial << "invalid access\n";

  // Or you can set the entire matrix at once using a c-array of the same dimensions and type, like so:
  float arrayA[3][3] = {{3.25,5.67,8.67},{4.55,7.23,9.00},{2.35,5.73,10.56}};
  A = arrayA;

  // You can also do this when you create the matrix, like this:
  float arrayB[3][3] = {{13.54,3.66,2.95},{3.22,7.54,5.12},{8.98,9.99,1.56}};
  Matrix<3,3> B(arrayB);

  // Now you can do some matrix math! The Matrix class supports addition and subtraction between matrices of the same size:
  Matrix<3,3> C = A + B;

  // You can't do addition between matrices of a different size. Since the class is built around templates, the compiler
  // will tell you if you've made a mistake. Try uncommenting the next line to see what I mean
  // A + v;

  // Next, it may not seem obvious, but that A + B operation actually returns a temporary matrix which is used to set C. 
  // To avoid creating that extra matrix and the associated memory useage you can use the Add functions to store the result directly into C, like so:
  Add(A,B,C);

  // Of course since those operations modify the 'C' parameter as they execute, you shouldn't write Add(C,B,C) for example.
  // Instead, you can use the Matrix class's unary operators, like so:
  C -= B;

  // Or like so:
  B += A;

  // As well as addition and subtraction, we can also do matrix multiplication. Note that again, the matrices, including 
  // the one in which the result will stored must have the appropriate dimensions. The compiler will let you know if they aren't
  Matrix<3,1> D = A * v;

  // Also, like addition / subtraction that operation will create an extra temporary matrix. As before, you can get around that by using Multiply, like so:
  Multiply(A,v,D);

  // As well as algebra, Matrix supports a few other matrix related operations including transposition:
  Matrix<1,3> D_T = Transpose(D);

  // The transpose function will return a temporary copy. If the matrix is square then that memory useage can be avoided using 
  // the Transpose member function. Be aware though, that this will modify the original matrix
  C.Transpose();

  // In situ Transposition of non-square matrices isn't allowed. To test this, try uncommenting the following line:
  // D.Transpose();

  // Matrix also supports inversion on square matrices via the invert function:
  Matrix<3,3> C_inv = Invert(C);

  // If the matrix is singular, the inversion won't work. In those cases Invert will still return a matrix but not a valid inverse. 
  // To check whether that has happened you can pass in pointer like so:
  int res;
  C_inv = Invert(C, &res); // after this call res will be 0 if the inversion was successful and non-zero otherwise.

  // Invert can also be called as a member function which modifies the matrix it is called on. It's more memory efficient so use it unless the original matrix needs unmodified
  C.Invert();

  // In addition to matrix math  there are functions to concatenate two matrices. You can concatenate matrices horizontally like so:
  Matrix<3,6> AleftOfB = HorzCat(A,B);

  // Or vertically, like so:
  Matrix<6,3> AonTopOfB = VertCat(A,B);

  // If you want to print out the value of any element in the array you can do that like so:
  Serial << "v(1): " << v(1) << '\n';

  // Alternatively, you can write the whole matrix to Serial using a C++ style inserter, like so:
  Serial << "B: " << B << '\n';

  // You can even write some quite complex compound statements very succinctly. For example:
  Serial << "identity matrix: " << AleftOfB * AonTopOfB - (A * A + B * B) + C * Invert(C) << '\n';

  // Or as a more real life example, here's how you might calculate an updated state estimate for a third order state space model:
  Matrix<3> x; Matrix<2> u; Matrix<3,2> G; Matrix<3,3> F; float dt;
  x += (F * x + G * u) * dt;

  // Enjoy!
}

void loop() { }
