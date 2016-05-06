#include <BasicLinearAlgebra.h>

/*
 * This example sketch shows how the template parameter can be nested to form high dimensional Matrices called 'Tensors'.
 * I should point out that I'm not really advocating using arduinos for Tensor computation, this is just an interesting
 * way of showing the flexbility and extensibility of the Matrix library
 */

void setup() 
{
  Serial.begin(115200);

  // The default underlying type of the Matrix class is float. If you want to use a different type, say int for example, then just pass it as a third template parameter like so:
  Matrix<3,3,int> intA;

  // From here you'll be able to do everything you'd be able to do with a float Matrix, but with int precision and memory useage.
  int array[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
  intA = array;

  // You can actually pass any datatype you like to the template and it'll do it's best to make a Matrix out of it.
  Matrix<3,3,unsigned char> charA;
  Matrix<3,3,double> doubleA; // etc

  // This includes parameters of type Matrix, meaning that you can declare matrices of more than two dimensions. For example:
  Matrix<4,4,Matrix<4> > cubeA, cubeB; // a 4x4x4 matrix (3rd order tensor)

  // And so on:
  Matrix<2,2,Matrix<2,2> > hyperA; // a 2x2x2x2 dimensional matrix (4th order tensor)
  Matrix<3,3,Matrix<3,3, Matrix<3,3> > > tensorA; // a 3x3x3x3x3x3 dimensional matrix (6th order tensor)

  // You can access the elements of an arbitrary rank tensor with the brackets operator like so:
  cubeA(0,1)(1) = cubeB(2,3)(3) = 56.34;
  hyperA(1,0)(1,1) = 56.34;
  tensorA(2,0)(1,2)(1,1) = 0.056;

  // Addition, subtraction and transposition all work as usual
  cubeA + cubeB;

  // As does concatenation
  Matrix<4,8,Matrix<4> > cubeAleftOfcubeB = HorzCat(cubeA,cubeB);

  // You can also do multiplication on square tensors with an even rank
  for(int i = 0; i < 2; i++)
      for(int j = 0; j < 2; j++)
          for(int k = 0; k < 2; k++)
              for(int l = 0; l < 2; l++)
                  hyperA(i,j)(k,l) = i+j+k+l;

  Matrix<2,2,Matrix<2,2> > hyperB = (hyperA * hyperA);

  // Everything can be printed too
  Serial << "Hyper B: " << hyperB;

  // Inversion doesn't work. If it did, it'd probably take quite a while for arduino to calculate anyway so maybe it's for the best
}

void loop() { }
