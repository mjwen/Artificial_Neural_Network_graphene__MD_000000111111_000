#ifndef HELPER_H_
#define HELPER_H_

#include <vector>
#include <string>
#include <iostream>

typedef double VectorOfSize6[6];

// helper routine declarations
void AllocateAndInitialize3DArray(double***& arrayPtr, int const extentZero,
                                  int const extentOne, int const extentTwo);
void Deallocate3DArray(double***& arrayPtr);

void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
                                  int const extentOne);
void Deallocate2DArray(double**& arrayPtr);

void AllocateAndInitialize1DArray(double*& arrayPtr, int const extent);
void Deallocate1DArray(double*& arrayPtr);


void virialTally2(bool isComputeVirial, bool isComputeParticleVirial,
  double de, double r, const double* dx, int i, int j,
  double *const virial, VectorOfSize6*const particleVirial);


void ComputeMeanAndStdev(std::vector<double> const & data, double& mean, double& stdev);


template<class T>
void print_vector(std::string name, std::vector<T> const& vec) {

  size_t size = vec.size();
  std::cout<<"vector: "<< name <<std::endl;
  std::cout<<"size: "<< size <<std::endl;
  for (size_t i=0; i<size; i++) {
    std::cout<< vec[i] << "  ";
  }
  std::cout<< std::endl <<std::endl;
}


#endif
