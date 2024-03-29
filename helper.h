#ifndef HELPER_H_
#define HELPER_H_

#include <vector>

// helper routine declarations
void AllocateAndInitialize3DArray(double***& arrayPtr, int const extentZero,
                                  int const extentOne, int const extentTwo);
void Deallocate3DArray(double***& arrayPtr);

void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
                                  int const extentOne);
void Deallocate2DArray(double**& arrayPtr);

void AllocateAndInitialize1DArray(double*& arrayPtr, int const extent);
void Deallocate1DArray(double*& arrayPtr);

void ComputeMeanAndStdev(std::vector<double> const & data, double& mean, double& stdev);

#endif
