//==============================================================================
//
// Implementation of helper functions
//
//==============================================================================
#include <iostream>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <iterator>
#include "helper.h"

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
void AllocateAndInitialize3DArray(double***& arrayPtr, int const extentZero,
		int const extentOne, int const extentTwo)
{
  arrayPtr = new double**[extentZero];
  arrayPtr[0] = new double*[extentZero * extentOne];
  arrayPtr[0][0] = new double[extentZero * extentOne * extentTwo];

  for (int i = 1; i < extentZero; ++i) {
    arrayPtr[i] = arrayPtr[i-1] + extentOne;
    arrayPtr[i][0] = arrayPtr[i-1][0] + extentOne * extentTwo;
  }

  for (int i = 0; i < extentZero; ++i) {
    for (int j = 1; j < extentOne; ++j) {
      arrayPtr[i][j] = arrayPtr[i][j-1] + extentTwo;
    }
  }

  // initialize
  for (int i = 0; i < extentZero; ++i) {
    for (int j = 0; j < extentOne; ++j) {
      for (int k = 0; k < extentTwo; ++k) {
        arrayPtr[i][j][k] = 0.0;
      }
    }
  }
}

// deallocate memory
void Deallocate3DArray(double***& arrayPtr) {
  if (arrayPtr != 0) {
    if (arrayPtr[0] != 0) {
      delete [] arrayPtr[0][0];
    }
    delete [] arrayPtr[0];
  }
  delete [] arrayPtr;

	// nullify pointer
	arrayPtr = 0;
}



void virialTally2(bool isComputeVirial, bool isComputeParticleVirial,
  double de, double r, const double* dx, int i, int j,
  double *const virial, VectorOfSize6*const  particleVirial)
{

  if (isComputeVirial || isComputeParticleVirial) {
    double vir[6];
    double v = de/r;

    vir[0] = v * dx[0] * dx[0];
    vir[1] = v * dx[1] * dx[1];
    vir[2] = v * dx[2] * dx[2];
    vir[3] = v * dx[1] * dx[2];
    vir[4] = v * dx[0] * dx[2];
    vir[5] = v * dx[0] * dx[1];

    if (isComputeVirial) {
      virial[0] += vir[0];
      virial[1] += vir[1];
      virial[2] += vir[2];
      virial[3] += vir[3];
      virial[4] += vir[4];
      virial[5] += vir[5];
    }

    if (isComputeParticleVirial) {
      vir[0] =0.5 * v * dx[0] * dx[0];
      vir[1] =0.5 * v * dx[1] * dx[1];
      vir[2] =0.5 * v * dx[2] * dx[2];
      vir[3] =0.5 * v * dx[1] * dx[2];
      vir[4] =0.5 * v * dx[0] * dx[2];
      vir[5] =0.5 * v * dx[0] * dx[1];
      particleVirial[i][0] += vir[0];
      particleVirial[i][1] += vir[1];
      particleVirial[i][2] += vir[2];
      particleVirial[i][3] += vir[3];
      particleVirial[i][4] += vir[4];
      particleVirial[i][5] += vir[5];

      particleVirial[j][0] += vir[0];
      particleVirial[j][1] += vir[1];
      particleVirial[j][2] += vir[2];
      particleVirial[j][3] += vir[3];
      particleVirial[j][4] += vir[4];
      particleVirial[j][5] += vir[5];
    }
  }

}








// compute the mean and standard deviation of vector
void ComputeMeanAndStdev(std::vector<double> const & v, double& mean, double& stdev)
{
  double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
  mean =  sum / v.size();

  double accum = 0.0;
  std::for_each (std::begin(v), std::end(v), [&](const double d) {
    accum += (d - mean) * (d - mean);
  });

  //stdev = std::sqrt(accum / (v.size()-1));    // corrected version
  stdev = std::sqrt(accum / (v.size()));    // uncorrected

}


