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

// allocate memory and set pointers
void AllocateAndInitialize1DArray(double*& arrayPtr, int const extent)
{
  arrayPtr = new double[extent];
  for (int i = 0; i < extent; ++i) {
    arrayPtr[i] = 0.0;
  }
}

// deallocate memory
void Deallocate1DArray(double*& arrayPtr) {
	if (arrayPtr != 0) delete [] arrayPtr;
	// nullify pointer
	arrayPtr = 0;
}


