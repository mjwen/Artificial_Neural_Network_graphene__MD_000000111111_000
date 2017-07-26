#ifndef HELPER_H_
#define HELPER_H_

// helper routine declarations
void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
                                  int const extentOne);
void Deallocate2DArray(double**& arrayPtr);

//TODO delete
void print_mat(double** mat, int rows, int cols);

#endif
