
# Basic Linear Algebra
A lightweight library for matrix operations on arduino.

This library is built around the template class Matrix which allows matrices of arbitrary dimensions to be defined and allocated using stack memory. It also allows for compile time checking of the compatibility of matrices used in any of the operations supported by the library.

Matrix supports: addition, subtraction and multiplication through operator overloads as well as inversion, transposition and concatenation through a series of helper template functions.

The matrix Invert function was adapted from code written by Charlie Matlack as part of the MatrixMath library that can be found here: http://playground.arduino.cc/Code/MatrixMath ... thanks!

To get started just go through the HowToUse sketch in the examples folder; that'll explain everything that most users will need to know to work with the library. For more advanced useage on how to take references to matrices check out the References, or for more on how to make a matrix containing elements of your own custom class have a look at the Tensor example.

