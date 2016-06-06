#include <BasicLinearAlgebra.h>
/*
 * This example sketch shows how to use references matrices to isolate a section of a larger matrix and use it in arithmetic with smaller matrices of the
 * appropriate dimensions. It was a bit of a pain to implement, but it's surprisingly useful so it's worth knowing about if you're using this library
 */

void setup()
{
  Serial.begin(115200);

  // If you've been through the HowToUse example you'll know that you can allocate a Matrix and explicitly specify it's type like so:
  Matrix<8,8,float> bigMatrix;

  // And as before it's a good idea to set the matrix before we use it
  bigMatrix.Fill(0);

  // The underlying storage for a Matrix is an array; but, since it's illegal to allocate an array of references, I've hijacked this Matrix type to implement a special matrix type with
  // slightly different behaviour to the regular matrix types.

  // To allocate a reference matrix just use the same syntax as a regular matrix and specify the underlying type as a reference by adding an & to the type.
  // As you would expect from a reference type you'll need to let the matrix know what it is referring to when you allocate it. Do that like so:
  Matrix<4,4,float&> bigMatrixRef(bigMatrix,4,2);

  // So bigMatrixRef now refers to a 4x4 submatrix of bigMatrix whose element (0,0) is positioned at the (4,2) element of bigMatrix
  // That is to say, if we set the (0,0) element of bigMatrixRef, we're effectively setting the (4,2) element of bigMatrix. So let's do that
  bigMatrixRef(0,0) = 45.67434;

  // And we can see that the original matrix has been set accordingly
  Serial << "bigMatrix(4,2): " << bigMatrix(4,2) << "\n";

  // As well as allocating a reference matrix as per above, you can also create them using the () operator of the parent matrix by passing in a pair of Range objects like so:
  Matrix<2,4,float&> anotherRef = bigMatrix(Range<2>(2),Range<4>(1)); // this creates a 2x4 reference matrix starting at element (2,1) of bigMatrix

  // For all other intents and purposes you can use matrix references as regular matrices. You can set them with an array..
  float arr[4][4] = {{23.44,43.23,12.45,6.23},{93.94,27.23,1.44,101.23},{1.23,3.21,4.56,8.76},{12.34,34.56,76.54,21.09}};
  bigMatrixRef = arr;

  // Or fill them
  anotherRef.Fill(5.678);

  // Do arithmetic with them
  anotherRef * bigMatrixRef;

  // Invert them
  Invert(bigMatrixRef);

  // Print them
  Serial << "bigMatrixRef: " << bigMatrixRef << "\n";

  // You can even make a reference to a reference matrix, do arithmetic with that and then print the result
  Serial  << "result of complicated operation: " << (anotherRef += bigMatrixRef(Range<2>(0),Range<4>(0))) << "\n";

  // The only thing that you can't really do is operate on two matrix references whose underlying memory overlaps,
  // particularly when doing matrix multiplication.

  // Lastly, let's look at what became of bigMatrix after all of this, you might be able to make out the values of bigMatrixRef and anotherRef in their respective areas of bigMatrix.
  Serial << "bigMatrix: " << bigMatrix << "\n";
}

void loop() { }

