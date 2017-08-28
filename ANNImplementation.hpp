//
// CDDL HEADER START
//
// The contents of this file are subject to the terms of the Common Development
// and Distribution License Version 1.0 (the "License").
//
// You can obtain a copy of the license at
// http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
// specific language governing permissions and limitations under the License.
//
// When distributing Covered Code, include this CDDL HEADER in each file and
// include the License file in a prominent location with the name LICENSE.CDDL.
// If applicable, add the following below this CDDL HEADER, with the fields
// enclosed by brackets "[]" replaced with your own identifying information:
//
// Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
//
// CDDL HEADER END
//

//
// Copyright (c) 2015, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Mingjian Wen
//


#ifndef ANN_IMPLEMENTATION_HPP_
#define ANN_IMPLEMENTATION_HPP_

#include "KIM_API_status.h"
#include "ANN.hpp"
#include "descriptor.h"
#include "network.h"
#include "helper.h"

#define DIM 3
#define ONE 1.0
#define HALF 0.5

#define MAX_PARAMETER_FILES 1


//==============================================================================
//
// Type definitions, and helper function prototypes
//
//==============================================================================

// type declaration for get neighbor functions
typedef int (GetNeighborFunction)(void**, int*, int*, int*, int*, int**,
                                  double**);
// type declaration for vector of constant dimension
typedef double VectorOfSizeDIM[DIM];


//==============================================================================
//
// Helper class definitions
//
//==============================================================================

// Iterator object for Locator mode access to neighbor list
class LocatorIterator
{
 private:
  KIM_API_model* const pkim_;
  GetNeighborFunction* const get_neigh_;
  int const baseconvert_;
  int const cachedNumberContributingParticles_;
  int request_;
  int const mode_;
 public:
  LocatorIterator(KIM_API_model* const pkim,
                  GetNeighborFunction* const get_neigh,
                  int const baseconvert,
                  int const cachedNumberContributingParticles,
                  int* const i,
                  int* const numnei,
                  int** const n1atom,
                  double** const pRij)
      : pkim_(pkim),
        get_neigh_(get_neigh),
        baseconvert_(baseconvert),
        cachedNumberContributingParticles_(cachedNumberContributingParticles),
        request_(-baseconvert_),  // set to first value (test-based indexing)
    mode_(1)  // locator mode
  {
    next(i, numnei, n1atom, pRij);
  }
  bool done() const
  {
    return !(request_ + baseconvert_ <= cachedNumberContributingParticles_);
  }
  int next(int* const i, int* const numnei, int** const n1atom,
           double** const pRij)
  {
    int ier;
    // Allow for request_ to be incremented to one more than contributing
    // without causing an error/warning from the openkim-api
    int req = std::min(request_,
                       cachedNumberContributingParticles_-baseconvert_-1);
    ier = (*get_neigh_)(
        reinterpret_cast<void**>(const_cast<KIM_API_model**>(&pkim_)),
        (int*) &mode_,
        &req,
        (int*) i,
        (int*) numnei,
        (int**) n1atom,
        (double**) pRij);
    *i += baseconvert_;  // adjust index of current particle

    ++request_;
    return ier;
  }
};


//==============================================================================
//
// Declaration of ANNImplementation class
//
//==============================================================================

//******************************************************************************
class ANNImplementation
{
 public:
  ANNImplementation(
      KIM_API_model* const pkim,
      char const* const parameterFileNames,
      int const parameterFileNameLength,
      int const numberParameterFiles,
      int* const ier);
  ~ANNImplementation();  // no explicit Destroy() needed here

  int Reinit(KIM_API_model* pkim);
  int Compute(KIM_API_model* pkim);

 private:
  // Constant values that never change
  //   Set in constructor (via SetConstantValues)
  //
  //
  // KIM API: Conventions
  int baseconvert_;

	//
  // ANNImplementation: constants
  int numberOfSpeciesIndex_;
  int numberOfParticlesIndex_;
  int particleSpeciesIndex_;
  int particleStatusIndex_;
  int coordinatesIndex_;
  int get_neighIndex_;
  int process_dEdrIndex_;
  int process_d2Edr2Index_;
  //
  // KIM API: Model Output indices
  int cutoffIndex_;
  int energyIndex_;
  int forcesIndex_;
  int particleEnergyIndex_;
  //
  // LennardJones612Implementation: constants
  int numberModelSpecies_;
  int numberUniqueSpeciesPairs_;



  // Constant values that are read from the input files and never change
  //   Set in constructor (via functions listed below)
  //
  //
  // KIM API: Model Fixed Parameters
  //   Memory allocated in   AllocateFixedParameterMemory()
  //   Memory deallocated in destructor
  //   Data set in ReadParameterFile routines
  // none
  //
  // KIM API: Model Free Parameters whose (pointer) values never change
  //   Memory allocated in   AllocateFreeParameterMemory() (from constructor)
  //   Memory deallocated in destructor
  //   Data set in ReadParameterFile routines OR by KIM Simulator
  double* cutoffs_;

  // Mutable values that only change when reinit() executes
  //   Set in Reinit (via SetReinitMutableValues)
  //
  //
  // KIM API: Model Fixed Parameters
  // none
  //
  // KIM API: Model Free Parameters
  // none
  //
  // ANNImplementation: values
  double** cutoffsSq2D_;

  // Mutable values that can change with each call to Reinit() and Compute()
  //   Memory may be reallocated on each call
  //
  //
  // ANNImplementation: values that change
  int cachedNumberOfParticles_;
  int cachedNumberContributingParticles_;

	// descriptor;
	Descriptor* descriptor_;
	NeuralNetwork* network_;



	// Helper methods
  //
  //
  // Related to constructor
  int SetConstantValues(KIM_API_model* const pkim);
  void AllocateFreeParameterMemory();
  static int OpenParameterFiles(
      KIM_API_model* const pkim,
      char const* const parameterFileNames,
      int const parameterFileNameLength,
      int const numberParameterFiles,
      FILE* parameterFilePointers[MAX_PARAMETER_FILES]);
  static void CloseParameterFiles(
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
      int const numberParameterFiles);
  int ProcessParameterFiles(
      KIM_API_model* const pkim,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
      int const numberParameterFiles);
  void getNextDataLine(FILE* const filePtr, char* const nextLine,
                       int const maxSize, int* endOfFileFlag);
  int getXdouble(char* linePtr, const int N, double* list);
  int getXint(char* linePtr, const int N, int* list);
  void lowerCase(char* linePtr);
  int ConvertUnits(KIM_API_model* const pkim);
  int RegisterKIMParameters(KIM_API_model* const pkim) const;
  int RegisterKIMFunctions(KIM_API_model* const pkim) const;
  //
  // Related to Reinit()
  int SetReinitMutableValues(KIM_API_model* const pkim);
  //
  // Related to Compute()
  int SetComputeMutableValues(KIM_API_model* const pkim,
                              bool& isComputeProcess_dEdr,
                              bool& isComputeProcess_d2Edr2,
                              bool& isComputeEnergy,
                              bool& isComputeForces,
                              bool& isComputeParticleEnergy,
                              int const*& particleSpecies,
                              GetNeighborFunction *& get_neigh,
                              VectorOfSizeDIM const*& coordinates,
                              double*& energy,
                              double*& particleEnergy,
                              VectorOfSizeDIM*& forces);
  int CheckParticleSpecies(KIM_API_model* const pkim,
                           int const* const particleSpecies) const;
  int GetComputeIndex(const bool& isComputeProcess_dEdr,
                      const bool& isComputeProcess_d2Edr2,
                      const bool& isComputeEnergy,
                      const bool& isComputeForces,
                      const bool& isComputeParticleEnergy) const;

  // compute functions
  template< class Iter,
            bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
            bool isComputeEnergy, bool isComputeForces,
            bool isComputeParticleEnergy>
  int Compute(KIM_API_model* const pkim,
              const int* const particleSpecies,
              GetNeighborFunction* const get_neigh,
              const VectorOfSizeDIM* const coordinates,
              double* const energy,
              VectorOfSizeDIM* const forces,
              double* const particleEnergy);


};

//==============================================================================
//
// Definition of ANNImplementation::Compute functions
//
// NOTE: Here we rely on the compiler optimizations to prune dead code
//       after the template expansions.  This provides high efficiency
//       and easy maintenance.
//
//==============================================================================

template< class Iter,
          bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
          bool isComputeEnergy, bool isComputeForces,
          bool isComputeParticleEnergy>
int ANNImplementation::Compute(
    KIM_API_model* const pkim,
    const int* const particleSpecies,
    GetNeighborFunction* const get_neigh,
    const VectorOfSizeDIM* const coordinates,
    double* const energy,
    VectorOfSizeDIM* const forces,
    double* const particleEnergy)
{
  int ier = KIM_STATUS_OK;
  const int Nparticles = cachedNumberOfParticles_;
  const int Ncontrib = cachedNumberContributingParticles_;

  if ((isComputeEnergy == false) &&
      (isComputeParticleEnergy == false) &&
      (isComputeForces == false) &&
      (isComputeProcess_dEdr == false) &&
      (isComputeProcess_d2Edr2 == false))
    return ier;

  // initialize energy and forces
  if (isComputeEnergy == true) {
    *energy = 0.0;
  }
  if (isComputeParticleEnergy == true) {
    for (int i = 0; i < Nparticles; ++i) {
      particleEnergy[i] = 0.0;
    }
  }
  if (isComputeForces == true) {
    for (int i = 0; i < Nparticles; ++i) {
      for (int j = 0; j < DIM; ++j)
        forces[i][j] = 0.0;
    }
  }

  // setting up generalzied coords matrix and the derivative matrix
  double** generalizedCoords;
//  double*** dGeneralizedCoords;
  int Ndescriptors = descriptor_->get_num_descriptors();
  AllocateAndInitialize2DArray(generalizedCoords, Ncontrib, Ndescriptors);

  // calculate generalized coordiantes
  //
  // Setup loop over contributing particles
  int ii = 0;
  int numnei = 0;
  int* n1atom = 0;
  double* pRij = 0;
  int const baseConvert = baseconvert_;
	double const* const* const  constCutoffsSq2D = cutoffsSq2D_;

	for (Iter iterator(pkim, get_neigh, baseConvert, Ncontrib, &ii, &numnei,
                     &n1atom, &pRij);
       iterator.done() == false;
       iterator.next(&ii, &numnei, &n1atom, &pRij))
  {
    int const numNei = numnei;
    int const * const n1Atom = n1atom;
    int const i = ii;
    int const iSpecies = particleSpecies[i];

    // Setup loop over neighbors of current particle
    for (int jj = 0; jj < numNei; ++jj)
    {
      // adjust index of particle neighbor
      int const j = n1Atom[jj] + baseConvert;
      int const jSpecies = particleSpecies[j];
      double rij[DIM];

			// Compute rij
			for (int dim = 0; dim < DIM; ++dim) {
				rij[dim] = coordinates[j][dim] - coordinates[i][dim];
			}

      // compute distance squared
      double const rijmag = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
      double const rcutij = sqrt(constCutoffsSq2D[iSpecies][jSpecies]);

      // if particles i and j not interact
      if (rijmag > rcutij) continue;

      // two-body descriptors
      for (size_t p=0; p<descriptor_->name.size(); p++) {

        if (descriptor_->name[p] != "g1" &&
            descriptor_->name[p] != "g2" &&
            descriptor_->name[p] != "g3") {
          continue;
        }
        int idx = descriptor_->starting_index[p];

        for(int q=0; q<descriptor_->num_param_sets[p]; q++) {

          double gc;
          if (descriptor_->name[p] == "g1") {
            descriptor_->sym_g1(rijmag, rcutij, gc);
          }
          else if (descriptor_->name[p] == "g2") {
            double eta = descriptor_->params[p][q][0];
            double Rs = descriptor_->params[p][q][1];
            descriptor_->sym_g2(eta, Rs, rijmag, rcutij, gc);
          }
          else if (descriptor_->name[p] == "g3") {
            double kappa = descriptor_->params[p][q][0];
            descriptor_->sym_g3(kappa, rijmag, rcutij, gc);
          }

          generalizedCoords[i][idx] += gc;
          idx += 1;

        } // loop over same descriptor but different parameter set
      } // loop over descriptors


      // three-body descriptors
      if (descriptor_->has_three_body == false) continue;

      for (int kk = jj+1; kk < numNei; ++kk) {

        // adjust index of particle neighbor
        int const k = n1Atom[kk] + baseConvert;
        int const kSpecies = particleSpecies[k];

        // Compute rik, rjk and their squares
        double rik[DIM];
        double rjk[DIM];
        for (int dim = 0; dim < DIM; ++dim) {
          rik[dim] = coordinates[k][dim] - coordinates[i][dim];
          rjk[dim] = coordinates[k][dim] - coordinates[j][dim];
        }
        double const rikmag = sqrt(rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2]);
        double const rjkmag = sqrt(rjk[0]*rjk[0] + rjk[1]*rjk[1] + rjk[2]*rjk[2]);
        double const rcutik = sqrt(constCutoffsSq2D[iSpecies][kSpecies]);
        double const rcutjk = sqrt(constCutoffsSq2D[jSpecies][kSpecies]);

        double const rvec[3] = {rijmag, rikmag, rjkmag};
        double const rcutvec[3] = {rcutij, rcutik, rcutjk};

        if (rikmag > rcutik) continue; // three-dody not interacting

        for (size_t p=0; p<descriptor_->name.size(); p++) {

          if (descriptor_->name[p] != "g4" &&
              descriptor_->name[p] != "g5") {
            continue;
          }
          int idx = descriptor_->starting_index[p];

          for(int q=0; q<descriptor_->num_param_sets[p]; q++) {

            double gc;
            if (descriptor_->name[p] == "g4") {
              double zeta = descriptor_->params[p][q][0];
              double lambda = descriptor_->params[p][q][1];
              double eta = descriptor_->params[p][q][2];
              descriptor_->sym_g4(zeta, lambda, eta, rvec, rcutvec, gc);
            }
            else if (descriptor_->name[p] == "g5") {
              double zeta = descriptor_->params[p][q][0];
              double lambda = descriptor_->params[p][q][1];
              double eta = descriptor_->params[p][q][2];
              descriptor_->sym_g5(zeta, lambda, eta, rvec, rcutvec, gc);
            }

            generalizedCoords[i][idx] += gc;
            idx += 1;

          } // loop over same descriptor but different parameter set
        }  // loop over descriptors
      }  // loop over kk (three body neighbors)
    }  // end of first neighbor loop
  }  // end of loop over contributing particles


  // centering and normalization
  if (descriptor_->center_and_normalize) {
    for (int i=0; i<Ncontrib; i++) {
      for (int j=0; j<Ndescriptors; j++) {
        generalizedCoords[i][j] = (generalizedCoords[i][j] -
            descriptor_->features_mean[j]) / descriptor_->features_std[j];
      }
    }
  }


  // NN feedforward
  network_->forward(generalizedCoords[0], Ncontrib, Ndescriptors);
  // NN backpropagation to compute derivative of energy w.r.t generalized coords
  network_->backward();

  // get access to derivatives of energy w.r.t generalized coords
  double** dEdGeneralizedCoords;
  int extentZero = Ncontrib;
  int extentOne = Ndescriptors;
  dEdGeneralizedCoords = new double*[extentZero];
  dEdGeneralizedCoords[0] = network_->get_grad_input();
  for (int i = 1; i < extentZero; ++i) {
    dEdGeneralizedCoords[i] = dEdGeneralizedCoords[i-1] + extentOne;
  }


  // Contribution to energy
  if (isComputeEnergy == true) {
    *energy = network_->get_sum_output();
  }

  // Contribution to particle energy
  if (isComputeParticleEnergy == true) {
    double* Epart;
    Epart = network_->get_output();
    for (int i=0; i<Ncontrib; i++) {
      particleEnergy[i] = Epart[i];
    }
  }


  // Compute derivative of energy w.r.t coords
  if ((isComputeProcess_dEdr == true) || (isComputeForces == true))
  {

    for (Iter iterator(pkim, get_neigh, baseConvert, Ncontrib, &ii, &numnei,
                     &n1atom, &pRij);
        iterator.done() == false;
        iterator.next(&ii, &numnei, &n1atom, &pRij))
    {
      int const numNei = numnei;
      int const * const n1Atom = n1atom;
      int const i = ii;
      int const iSpecies = particleSpecies[i];

      // Setup loop over neighbors of current particle
      for (int jj = 0; jj < numNei; ++jj)
      {
        // adjust index of particle neighbor
        int const j = n1Atom[jj] + baseConvert;
        int const jSpecies = particleSpecies[j];
        double rij[DIM];

        // Compute rij
        for (int dim = 0; dim < DIM; ++dim) {
          rij[dim] = coordinates[j][dim] - coordinates[i][dim];
        }

        // compute distance squared
        double const rijmag = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        double const rcutij = sqrt(constCutoffsSq2D[iSpecies][jSpecies]);

        // if particles i and j not interact
        if (rijmag > rcutij) continue;

        // two-body descriptors
        for (size_t p=0; p<descriptor_->name.size(); p++) {

          if (descriptor_->name[p] != "g1" &&
              descriptor_->name[p] != "g2" &&
              descriptor_->name[p] != "g3") {
            continue;
          }
          int idx = descriptor_->starting_index[p];

          for(int q=0; q<descriptor_->num_param_sets[p]; q++) {

            double gc;
            double dgcdr_two;
            if (descriptor_->name[p] == "g1") {
              descriptor_->sym_d_g1(rijmag, rcutij, gc, dgcdr_two);
            }
            else if (descriptor_->name[p] == "g2") {
              double eta = descriptor_->params[p][q][0];
              double Rs = descriptor_->params[p][q][1];
              descriptor_->sym_d_g2(eta, Rs, rijmag, rcutij, gc, dgcdr_two);
            }
            else if (descriptor_->name[p] == "g3") {
              double kappa = descriptor_->params[p][q][0];
              descriptor_->sym_d_g3(kappa, rijmag, rcutij, gc, dgcdr_two);
            }

            // centering and normalization
            if (descriptor_->center_and_normalize) {
              dgcdr_two /= descriptor_->features_std[idx];
            }

            for (int kdim = 0; kdim < DIM; ++kdim) {
              double phi = dEdGeneralizedCoords[i][idx]*dgcdr_two*rij[kdim]/rijmag;
              forces[i][kdim] += phi;
              forces[j][kdim] -= phi;
            }
            idx += 1;

          } // loop over same descriptor but different parameter set
        } // loop over descriptors


        // three-body descriptors
        if (descriptor_->has_three_body == false) continue;

        for (int kk = jj+1; kk < numNei; ++kk) {

          // adjust index of particle neighbor
          int const k = n1Atom[kk] + baseConvert;
          int const kSpecies = particleSpecies[k];

          // Compute rik, rjk and their squares
          double rik[DIM];
          double rjk[DIM];
          for (int dim = 0; dim < DIM; ++dim) {
            rik[dim] = coordinates[k][dim] - coordinates[i][dim];
            rjk[dim] = coordinates[k][dim] - coordinates[j][dim];
          }
          double const rikmag = sqrt(rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2]);
          double const rjkmag = sqrt(rjk[0]*rjk[0] + rjk[1]*rjk[1] + rjk[2]*rjk[2]);
          double const rcutik = sqrt(constCutoffsSq2D[iSpecies][kSpecies]);
          double const rcutjk = sqrt(constCutoffsSq2D[jSpecies][kSpecies]);

          double const rvec[3] = {rijmag, rikmag, rjkmag};
          double const rcutvec[3] = {rcutij, rcutik, rcutjk};

          if (rikmag > rcutik) continue; // three-dody not interacting

          for (size_t p=0; p<descriptor_->name.size(); p++) {

            if (descriptor_->name[p] != "g4" &&
                descriptor_->name[p] != "g5") {
              continue;
            }
            int idx = descriptor_->starting_index[p];

            for(int q=0; q<descriptor_->num_param_sets[p]; q++) {

              double gc;
              double dgcdr_three[3];
              if (descriptor_->name[p] == "g4") {
                double zeta = descriptor_->params[p][q][0];
                double lambda = descriptor_->params[p][q][1];
                double eta = descriptor_->params[p][q][2];
                descriptor_->sym_d_g4(zeta, lambda, eta, rvec, rcutvec, gc, dgcdr_three);
              }
              else if (descriptor_->name[p] == "g5") {
                double zeta = descriptor_->params[p][q][0];
                double lambda = descriptor_->params[p][q][1];
                double eta = descriptor_->params[p][q][2];
                descriptor_->sym_d_g5(zeta, lambda, eta, rvec, rcutvec, gc, dgcdr_three);
              }

              // centering and normalization
              if (descriptor_->center_and_normalize) {
                dgcdr_three[0] /= descriptor_->features_std[idx];
                dgcdr_three[1] /= descriptor_->features_std[idx];
                dgcdr_three[2] /= descriptor_->features_std[idx];
              }

              for (int kdim = 0; kdim < DIM; ++kdim) {
                double phi_ij = dEdGeneralizedCoords[i][idx]*
                  dgcdr_three[0]*rij[kdim]/rijmag;
                double phi_ik = dEdGeneralizedCoords[i][idx]*
                  dgcdr_three[1]*rik[kdim]/rikmag;
                double phi_jk = dEdGeneralizedCoords[i][idx]*
                  dgcdr_three[2]*rjk[kdim]/rjkmag;
                forces[i][kdim] += phi_ij + phi_ik;
                forces[j][kdim] += -phi_ij + phi_jk;
                forces[k][kdim] += -phi_ik - phi_jk;
              }
              idx += 1;

            } // loop over same descriptor but different parameter set
          }  // loop over descriptors
        }  // loop over kk (three body neighbors)
      }  // loop over first neighbor
    }  // loop over i atoms
  } // compute force



  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

#endif  // ANN_IMPLEMENTATION_HPP_
