#include <BasicLinearAlgebra.h>

/*
 * This example sketch shows how to use a reference matrix to isolate a section of a larger matrix and use it in arithmetic with smaller matrices of the
 * appropriate dimensions.
 */

void setup()
{
  Serial.begin(115200);

  // If you've been through the HowToUse example you'll know that you can allocate a Matrix and explicitly specify it's type like so:
  Matrix<8,8,float> bigMatrix;

  // And as before it's a good idea to fill the matrix before we use it
  bigMatrix.Fill(0);

  // Now we'll make a refence matrix which refers to bigMatrix. The reference has no internal memory of it's own, it just goes and gets whatever it needs from bigMatrix when it takes part in
  // an operation. It'll also stores the result of an operation there if needs be too. The handy thing about reference matrices is we can use them to select just a submatrix of the original
  // and use that in operations with smaller matrices.

  // To allocate a reference matrix use the RefMatrix type and pass in the template parameters in the same way you would a regular Matrix. Since we can't have references floating around
  // with no memory to refer to, you must pass the original matrix to the contructor along with a row and column offset.
  RefMatrix<4,4,float> bigMatrixRef(bigMatrix,4,2);

  // So bigMatrixRef now refers to a 4x4 submatrix of bigMatrix whose element (0,0) is positioned at the (4,2) element of bigMatrix
  // That is to say, if we set the (0,0) element of bigMatrixRef, we're effectively setting the (4,2) element of bigMatrix. So let's do that
  bigMatrixRef(0,0) = 45.67434;

  // And we can see that the original matrix has been set accordingly
  Serial << "bigMatrix(4,2): " << bigMatrix(4,2) << "\n";

  // As well as allocating a reference matrix as per above, you can also create them using the Submatrix() function of the parent matrix and passing in a pair of Range objects like so:
  RefMatrix<2,4,float> anotherRef = bigMatrix.Submatrix(Range<2>(2),Range<4>(1)); // this creates a 2x4 reference matrix starting at element (2,1) of bigMatrix

  // For all other intents and purposes you can use matrix references as regular matrices.

  //You can set them with an array...
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
  Serial  << "result of convoluted operation: " << (anotherRef += bigMatrixRef.Submatrix(Range<2>(0),Range<4>(0))) << "\n";

  // The only thing that you can't really do is operate on two matrix references whose underlying memory overlaps,
  // particularly when doing matrix multiplication.

  // Lastly, let's look at what became of bigMatrix after all of this, you might be able to make out the values of bigMatrixRef and anotherRef in their respective areas of bigMatrix.
  Serial << "bigMatrix: " << bigMatrix << "\n";
}

void loop() { }
