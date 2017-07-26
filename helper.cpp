//==============================================================================
//
// Implementation of helper functions
//
//==============================================================================
#include <iostream>
#include "helper.h"

// allocate memory and set pointers
void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
		int const extentOne)
{
	arrayPtr = new double*[extentZero];
	arrayPtr[0] = new double[extentZero * extentOne];
	for (int i = 1; i < extentZero; ++i) {
		arrayPtr[i] = arrayPtr[i-1] + extentOne;
	}

	// initialize
	for (int i = 0; i < extentZero; ++i) {
		for (int j = 0; j < extentOne; ++j) {
			arrayPtr[i][j] = 0.0;
		}
	}
}

// deallocate memory
void Deallocate2DArray(double**& arrayPtr) {
	if (arrayPtr != 0) delete [] arrayPtr[0];
	delete [] arrayPtr;

	// nullify pointer
	arrayPtr = 0;
}



//TODO delete
void print_mat(double** mat, int rows, int cols) {
  for (int m=0; m<rows; m++) {
    for (int n=0; n<cols; n++) {
      std::cout<<mat[m][n]<< " ";
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
}



