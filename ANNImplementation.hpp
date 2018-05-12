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

#include <iomanip>
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
  double* cutoffs_samelayer_;

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
  double** cutoffsSq2D_samelayer_;

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


  // configurations of last computation
  int numberOfParticles_last_call_;
  std::vector<int> particleSpecies_last_call_;
  std::vector<double> coordinates_last_call_;

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

  if ((isComputeEnergy == false) &&
      (isComputeParticleEnergy == false) &&
      (isComputeForces == false) &&
      (isComputeProcess_dEdr == false) &&
      (isComputeProcess_d2Edr2 == false))
    return ier;

  if (isComputeProcess_d2Edr2 == true)
    std::cerr<<"KIM potential: not supported ComputeProcess_d2Edr2"<<std::endl;

  bool need_dE = (isComputeProcess_dEdr == true) || (isComputeForces == true);

  // ANNImplementation: values that does not change
  const int Nparticles = cachedNumberOfParticles_;
  const int Ncontrib = cachedNumberContributingParticles_;


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


  // Allocate memory for precompute values of sym g4
  size_t n_lambda = descriptor_->g4_distinct_lambda.size();
  size_t n_zeta = descriptor_->g4_distinct_zeta.size();
  size_t n_eta = descriptor_->g4_distinct_eta.size();
  double** costerm;
  double*** dcosterm_dr;
  double* eterm;
  double** determ_dr;
  AllocateAndInitialize2DArray(costerm, n_lambda, n_zeta);
  AllocateAndInitialize3DArray(dcosterm_dr, n_lambda, n_zeta, 3);
  AllocateAndInitialize1DArray(eterm, n_eta);
  AllocateAndInitialize2DArray(determ_dr, n_eta, 3);


  // number of descriptors
  const int Ndescriptors = descriptor_->get_num_descriptors();
  const int Ndescriptors_two = descriptor_->get_num_descriptors_two_body();
  const int Ndescriptors_three = descriptor_->get_num_descriptors_three_body();
#ifdef DEBUG
  std::cout<<"@Ndescriptors = " << Ndescriptors<<std::endl;
  std::cout<<"@Ndescriptors_two = " << Ndescriptors_two<<std::endl;
  std::cout<<"@Ndescriptors_three = " << Ndescriptors_three<<std::endl;
#endif


  // index map between 1D two-body, three-body descriptors and global 1D descriptor
  int* map_t_desc_two = new int[Ndescriptors_two];
  int* map_t_desc_three = new int[Ndescriptors_three];
  int t_two = 0;
  int t_three = 0;
  for (size_t p=0; p<descriptor_->name.size(); p++) {
    for(int q=0; q<descriptor_->num_param_sets[p]; q++) {

      if (strcmp(descriptor_->name[p], "g1") == 0 ||
          strcmp(descriptor_->name[p], "g2") == 0 ||
          strcmp(descriptor_->name[p], "g3") == 0) {
        map_t_desc_two[t_two] = descriptor_->get_global_1D_index(p, q);
        t_two += 1;
      }
      else if (strcmp(descriptor_->name[p], "g4") == 0 ||
          strcmp(descriptor_->name[p], "g5") == 0) {
        map_t_desc_three[t_three] = descriptor_->get_global_1D_index(p, q);
        t_three += 1;
      }

    }
  }


  // allocate memory based on approxinate number of neighbors
  // memory will be relallocated if numNei is larger than approx_numNei
  double** dGCdr_two;
  double*** dGCdr_three;
  int approx_numNei = 100;
  int Npairs_two = approx_numNei;
  int Npairs_three = approx_numNei*(approx_numNei-1)/2;
  AllocateAndInitialize2DArray(dGCdr_two, Npairs_two, Ndescriptors_two);
  AllocateAndInitialize3DArray(dGCdr_three, Npairs_three, Ndescriptors_three, 3);



  // calculate generalized coordinates
  //
  // Setup loop over contributing particles
  for (int ii=0; ii<Ncontrib; ii++) {

    int const i = ii;
    int const iSpecies = particleSpecies[i];

    // get neighbors of atom i
    int one = 1;
    int dummy;
    int numnei = 0;
    int* n1atom = 0;
    double* pRij = 0;
    int const baseConvert = baseconvert_;
    int request = i - baseConvert;
    get_neigh( reinterpret_cast<void**>(const_cast<KIM_API_model**>(&pkim)),
        &one, &request, &dummy, &numnei, &n1atom, &pRij);

    int const numNei = numnei;
    int const * const n1Atom = n1atom;


    // generalized coords of atom i and its derivatives w.r.t. pair distances
    double* GC;
    AllocateAndInitialize1DArray(GC, Ndescriptors);

    const int Npairs_two = numNei;
    const int Npairs_three = numNei*(numNei-1)/2;
    // realloate memory is numNei is larger than approx_numNei
    if (numNei > approx_numNei) {
      Deallocate2DArray(dGCdr_two);
      Deallocate3DArray(dGCdr_three);
      AllocateAndInitialize2DArray(dGCdr_two, Npairs_two, Ndescriptors_two);
      AllocateAndInitialize3DArray(dGCdr_three, Npairs_three, Ndescriptors_three, 3);
      approx_numNei = numNei;
    }


    // Setup loop over neighbors of current particle
    for (int jj = 0; jj < numNei; ++jj)
    {
      // adjust index of particle neighbor
      int const j = n1Atom[jj] + baseConvert;
      int const jSpecies = particleSpecies[j];

      // cutoff between ij
      double rcutij = sqrt(cutoffsSq2D_[iSpecies][jSpecies]);

      // Compute rij
      double rij[DIM];
      for (int dim = 0; dim < DIM; ++dim) {
        rij[dim] = coordinates[j][dim] - coordinates[i][dim];
      }
      double const rijsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
      double const rijmag = sqrt(rijsq);

      // if particles i and j not interact
      if (rijmag > rcutij) continue;

      // pre-compute two-body cut function
      double fcij = descriptor_->cutoff(rijmag, rcutij);
      double dfcij = descriptor_->d_cutoff(rijmag, rcutij);


      int s_two = jj;   // row index of dGCdr_two
      int t_two = 0;   // column index of dGCdr_two
      for (size_t p=0; p<descriptor_->name.size(); p++) {

        if (strcmp(descriptor_->name[p], "g1") != 0 &&
            strcmp(descriptor_->name[p], "g2") != 0 &&
            strcmp(descriptor_->name[p], "g3") != 0) {
          continue;
        }


        for(int q=0; q<descriptor_->num_param_sets[p]; q++) {

          double gc;
          double dgcdr_two;

//          if (strcmp(descriptor_->name[p], "g1") == 0) {
//            if (need_dE) {
//              descriptor_->sym_d_g1(rijmag, rcutij, gc, dgcdr_two);
//            } else {
//              descriptor_->sym_g1(rijmag, rcutij, gc);
//            }
//          }
//          else if (strcmp(descriptor_->name[p], "g2") == 0) {
          double eta = descriptor_->params[p][q][0];
          double Rs = descriptor_->params[p][q][1];

//          if (need_dE) {
            descriptor_->sym_d_g2(eta, Rs, rijmag, rcutij, fcij, dfcij, gc, dgcdr_two);
//          } else {
//            descriptor_->sym_g2(eta, Rs, rijmag, rcutij, gc);
//          }
            //          }
//          else if (strcmp(descriptor_->name[p], "g3") == 0) {
//            double kappa = descriptor_->params[p][q][0];
//            if (need_dE) {
//              descriptor_->sym_d_g3(kappa, rijmag, rcutij, gc, dgcdr_two);
//            } else {
//              descriptor_->sym_g3(kappa, rijmag, rcutij, gc);
//            }
//          }
//

          int desc_idx = descriptor_->get_global_1D_index(p, q);
          GC[desc_idx] += gc;
//          if (need_dE) {
            dGCdr_two[s_two][t_two] = dgcdr_two;
            t_two += 1;
//          }

        } // loop over same descriptor but different parameter set
      } // loop over descriptors



      // three-body descriptors
//      if (descriptor_->has_three_body == false) continue;

      for (int kk = jj+1; kk < numNei; ++kk) {

        // adjust index of particle neighbor
        int const k = n1Atom[kk] + baseConvert;
        int const kSpecies = particleSpecies[k];

        // cutoff between ik and jk
        double const rcutik = sqrt(cutoffsSq2D_samelayer_[iSpecies][kSpecies]);
        double const rcutjk = sqrt(cutoffsSq2D_samelayer_[jSpecies][kSpecies]);

        // Compute rik, rjk and their squares
        double rik[DIM];
        double rjk[DIM];
        for (int dim = 0; dim < DIM; ++dim) {
          rik[dim] = coordinates[k][dim] - coordinates[i][dim];
          rjk[dim] = coordinates[k][dim] - coordinates[j][dim];
        }
        double const riksq = rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2];
        double const rjksq = rjk[0]*rjk[0] + rjk[1]*rjk[1] + rjk[2]*rjk[2];
        double const rikmag = sqrt(riksq);
        double const rjkmag = sqrt(rjksq);


        double const rvec[3] = {rijmag, rikmag, rjkmag};
        double const rcutvec[3] = {rcutij, rcutik, rcutjk};

        if (rikmag > rcutik) continue; // three-dody not interacting

        //@TEMP only for g4, should delete this if we have g5
        if (rjkmag > rcutjk) continue; // only for g4, not for g4


        // cutoff term, i.e. the product of fc(rij), fc(rik), and fc(rjk)
        double fcik = cut_cos(rikmag, rcutik);
        double fcjk = cut_cos(rjkmag, rcutjk);
        double dfcik = d_cut_cos(rikmag, rcutik);
        double dfcjk = d_cut_cos(rjkmag, rcutjk);
        double fcprod = fcij*fcik*fcjk;
        double dfcprod_dr[3];    // dfcprod/drij, dfcprod/drik, dfcprod/drjk
        dfcprod_dr[0] = dfcij*fcik*fcjk;
        dfcprod_dr[1] = dfcik*fcij*fcjk;
        dfcprod_dr[2] = dfcjk*fcij*fcik;


        int s_three = (numNei-1 + numNei-jj)*jj/2 + (kk-jj-1); // row index of dGCdr_three
#ifdef DEBUG
        std::cout<<"@numNei="<<numNei<<std::endl;
        std::cout<<"@jj="<<jj<<std::endl;
        std::cout<<"@kk="<<kk<<std::endl;
        std::cout<<"@s_three="<<s_three<<std::endl;
#endif

        int t_three = 0; // column index of dGCdr_three

        for (size_t p=0; p<descriptor_->name.size(); p++) {

          if (strcmp(descriptor_->name[p], "g4") != 0 &&
              strcmp(descriptor_->name[p], "g5") != 0) {
            continue;
          }

          // precompute recurring values in cosine terms and exponential terms
          descriptor_->precompute_g4(rijmag, rikmag, rjkmag, rijsq, riksq, rjksq,
              n_lambda, n_zeta, n_eta, costerm, dcosterm_dr, eterm, determ_dr);


          for(int q=0; q<descriptor_->num_param_sets[p]; q++) {

            double gc;
            double dgcdr_three[3];

//            if (strcmp(descriptor_->name[p], "g4") == 0) {

//              if (need_dE) {

                // get values from precomputed
                int izeta = descriptor_->g4_lookup_zeta[q];
                int ilam = descriptor_->g4_lookup_lambda[q];
                int ieta = descriptor_->g4_lookup_eta[q];

                double ct = costerm[ilam][izeta];
                double dct[3];
                dct[0] = dcosterm_dr[ilam][izeta][0];
                dct[1] = dcosterm_dr[ilam][izeta][1];
                dct[2] = dcosterm_dr[ilam][izeta][2];

                double et = eterm[ieta];
                double det[3];
                det[0] = determ_dr[ieta][0];
                det[1] = determ_dr[ieta][1];
                det[2] = determ_dr[ieta][2];

                descriptor_->sym_d_g4_2(rvec, rcutvec, fcprod, dfcprod_dr, ct, dct, et, det, gc, dgcdr_three);

//              } else {
//                descriptor_->sym_g4(zeta, lambda, eta, rvec, rcutvec, gc);
//              }
//            }
//            else if (strcmp(descriptor_->name[p], "g5") == 0) {
//              double zeta = descriptor_->params[p][q][0];
//              double lambda = descriptor_->params[p][q][1];
//              double eta = descriptor_->params[p][q][2];
//              if (need_dE) {
//                descriptor_->sym_d_g5(zeta, lambda, eta, rvec, rcutvec, gc, dgcdr_three);
//              } else {
//                descriptor_->sym_g5(zeta, lambda, eta, rvec, rcutvec, gc);
//              }
//            }
//

            int desc_idx = descriptor_->get_global_1D_index(p, q);
            GC[desc_idx] += gc;
//            if (need_dE) {
              dGCdr_three[s_three][t_three][0] = dgcdr_three[0];
              dGCdr_three[s_three][t_three][1] = dgcdr_three[1];
              dGCdr_three[s_three][t_three][2] = dgcdr_three[2];
              t_three += 1;
//            }

          } // loop over same descriptor but different parameter set
        }  // loop over descriptors
      }  // loop over kk (three body neighbors)

    }  // loop over jj


    // Note, in KIM-API v1 for full neighbor list, `Ncontrib` is actually the total
    // number of atoms. For non-contributing atoms, its numNei is set to 0.
    // So here we need to continue immediately, otherwise energy of noncontributing
    // atoms will be incorrectly added
    if(numNei == 0) continue;


    /*
    //@DEBUG delete debug (print generalized coords normalized)
    std::cout<<"\n# Debug descriptor values before normalization" << std::endl;
    std::cout<<"# atom id    descriptor values ..." << std::endl;
    std::cout<< ii <<"    ";
    for(int j=0; j<Ndescriptors; j++) {
      printf("%.15f ",GC[j]);
    }
    std::cout<<std::endl;
     */


    // centering and normalization
    if (descriptor_->center_and_normalize) {

      for (int t=0; t<Ndescriptors; t++) {
        GC[t] = (GC[t] - descriptor_->features_mean[t]) / descriptor_->features_std[t];
      }

// We moved this below
/*
      if (need_dE) {
        for (int s=0; s<Npairs_two; s++) {
          for (int t=0; t<Ndescriptors_two; t++) {
            int desc_idx = map_t_desc_two[t];
            dGCdr_two[s][t] /= descriptor_->features_std[desc_idx];
          }
        }

        for (int s=0; s<Npairs_three; s++) {
          for (int t=0; t<Ndescriptors_three; t++) {
            int desc_idx = map_t_desc_three[t];
            dGCdr_three[s][t][0] /= descriptor_->features_std[desc_idx];
            dGCdr_three[s][t][1] /= descriptor_->features_std[desc_idx];
            dGCdr_three[s][t][2] /= descriptor_->features_std[desc_idx];
          }
        }
      }
*/
    }



    /*
    //@DEBUG delete debug (print generalized coords normalized)
    std::cout<<"\n\n# Debug descriptor values after normalization" << std::endl;
    std::cout<<"# atom id    descriptor values ..." << std::endl;
    std::cout<< ii <<"    ";
    for(int j=0; j<Ndescriptors; j++) {
      printf("%.15f ",GC[j]);
    }
    std::cout<<std::endl;
    */



    // NN feedforward
    network_->forward(GC, 1, Ndescriptors);

    // NN backpropagation to compute derivative of energy w.r.t generalized coords
    double* dEdGC;
    if (need_dE) {
      network_->backward();
      dEdGC = network_->get_grad_input();
    }

    double Ei = 0.;
    if (isComputeEnergy == true || isComputeParticleEnergy == true) {
      Ei = network_->get_sum_output();
    }

    // Contribution to energy
    if (isComputeEnergy == true) {
      *energy += Ei;
    }

    // Contribution to particle energy
    if (isComputeParticleEnergy == true) {
      particleEnergy[i] += Ei;
    }


    // Contribution to forces and virial
    if (need_dE) {

      // neighboring atoms of i
      for (int jj = 0; jj < numNei; ++jj) {

        // adjust index of particle neighbor
        int const j = n1Atom[jj] + baseConvert;
        int const jSpecies = particleSpecies[j];

        // cutoff between ij
        double rcutij = sqrt(cutoffsSq2D_[iSpecies][jSpecies]);

        // Compute rij
        double rij[DIM];
        const double* prij = rij;
        for (int dim = 0; dim < DIM; ++dim) {
          rij[dim] = coordinates[j][dim] - coordinates[i][dim];
        }
        double const rijsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
        double const rijmag = sqrt(rijsq);

        // if particles i and j not interact
        if (rijmag > rcutij) continue;


        // two-body descriptors

        int s_two = jj;
        double dEdr_two = 0;
        for (int t=0; t<Ndescriptors_two; t++) {
          int desc_idx = map_t_desc_two[t];
          if (descriptor_->center_and_normalize) {
            dEdr_two += dGCdr_two[s_two][t] * (dEdGC[desc_idx] / descriptor_->features_std[desc_idx]);
          } else {
            dEdr_two += dGCdr_two[s_two][t] * dEdGC[desc_idx];
          }
        }

        // forces
        if (isComputeForces) {
          for (int dim = 0; dim < DIM; ++dim) {
            double pair = dEdr_two*rij[dim]/rijmag;
            forces[i][dim] += pair;  // for i atom
            forces[j][dim] -= pair;  // for neighboring atoms of i
          }
        }

        // process_dEdr
        if (isComputeProcess_dEdr) {
          int ier = pkim->process_dEdr(const_cast<KIM_API_model**>(&pkim),
              const_cast<double*>(&dEdr_two),
              const_cast<double*>(&rijmag),
              const_cast<double**>(&prij),
              const_cast<int*>(&i),
              const_cast<int*>(&j));
          if (ier < KIM_STATUS_OK) {
            pkim->report_error(__LINE__, __FILE__, "process_dEdr", ier);
            return ier;
          }
        }


        // three-body descriptors
        for (int kk = jj+1; kk < numNei; ++kk) {

          // adjust index of particle neighbor
          int const k = n1Atom[kk] + baseConvert;
          int const kSpecies = particleSpecies[k];

          // cutoff between ik and jk
          double const rcutik = sqrt(cutoffsSq2D_samelayer_[iSpecies][kSpecies]);
          double const rcutjk = sqrt(cutoffsSq2D_samelayer_[jSpecies][kSpecies]);

          // Compute rik, rjk and their squares
          double rik[DIM];
          double rjk[DIM];
          const double* prik = rik;
          const double* prjk = rjk;
          for (int dim = 0; dim < DIM; ++dim) {
            rik[dim] = coordinates[k][dim] - coordinates[i][dim];
            rjk[dim] = coordinates[k][dim] - coordinates[j][dim];
          }
          double const riksq = rik[0]*rik[0] + rik[1]*rik[1] + rik[2]*rik[2];
          double const rjksq = rjk[0]*rjk[0] + rjk[1]*rjk[1] + rjk[2]*rjk[2];
          double const rikmag = sqrt(riksq);
          double const rjkmag = sqrt(rjksq);

          if (rikmag > rcutik) continue; // three-dody not interacting

          //@TEMP only for g4
          if (rjkmag > rcutjk) continue; // only for g4, not for g4


          int s_three = (numNei-1 + numNei-jj)*jj/2 + (kk-jj-1); // row index of dGCdr_three
          double dEdr_three[3] = {0, 0, 0};
          for (int t=0; t<Ndescriptors_three; t++) {
            int desc_idx = map_t_desc_three[t];
            if (descriptor_->center_and_normalize) {
              dEdr_three[0] += dGCdr_three[s_three][t][0] * (dEdGC[desc_idx] / descriptor_->features_std[desc_idx]);   // dEdrij
              dEdr_three[1] += dGCdr_three[s_three][t][1] * (dEdGC[desc_idx] / descriptor_->features_std[desc_idx]);   // dEdrik
              dEdr_three[2] += dGCdr_three[s_three][t][2] * (dEdGC[desc_idx] / descriptor_->features_std[desc_idx]);   // dEdrjk

            } else {
              dEdr_three[0] += dGCdr_three[s_three][t][0] * dEdGC[desc_idx];   // dEdrij
              dEdr_three[1] += dGCdr_three[s_three][t][1] * dEdGC[desc_idx];   // dEdrik
              dEdr_three[2] += dGCdr_three[s_three][t][2] * dEdGC[desc_idx];   // dEdrjk
            }
          }


          // forces
          if (isComputeForces) {
            for (int dim = 0; dim < DIM; ++dim) {
              double pair_ij = dEdr_three[0]*rij[dim]/rijmag;
              double pair_ik = dEdr_three[1]*rik[dim]/rikmag;
              double pair_jk = dEdr_three[2]*rjk[dim]/rjkmag;
              forces[i][dim] += pair_ij + pair_ik;    // for i atom
              forces[j][dim] += -pair_ij + pair_jk;    // for neighboring atoms of i
              forces[k][dim] += -pair_ik - pair_jk;    // for neighboring atoms of i
            }
          }

          // process_dEdr
          if (isComputeProcess_dEdr) {
            int ier;
            ier = pkim->process_dEdr(const_cast<KIM_API_model**>(&pkim),
                const_cast<double*>(&dEdr_three[0]),
                const_cast<double*>(&rijmag),
                const_cast<double**>(&prij),
                const_cast<int*>(&i),
                const_cast<int*>(&j));
            if (ier < KIM_STATUS_OK) {
              pkim->report_error(__LINE__, __FILE__, "process_dEdr", ier);
              return ier;
            }

            ier = pkim->process_dEdr(const_cast<KIM_API_model**>(&pkim),
                const_cast<double*>(&dEdr_three[1]),
                const_cast<double*>(&rikmag),
                const_cast<double**>(&prik),
                const_cast<int*>(&i),
                const_cast<int*>(&k));
            if (ier < KIM_STATUS_OK) {
              pkim->report_error(__LINE__, __FILE__, "process_dEdr", ier);
              return ier;
            }

            ier = pkim->process_dEdr(const_cast<KIM_API_model**>(&pkim),
                const_cast<double*>(&dEdr_three[2]),
                const_cast<double*>(&rjkmag),
                const_cast<double**>(&prjk),
                const_cast<int*>(&j),
                const_cast<int*>(&k));
            if (ier < KIM_STATUS_OK) {
              pkim->report_error(__LINE__, __FILE__, "process_dEdr", ier);
              return ier;
            }
          }


        }  // loop over kk
      }  // loop over jj
    }  // need_dE


    Deallocate1DArray(GC);

  }  // loop over ii, i.e. contributing particles


  Deallocate2DArray(costerm);
  Deallocate3DArray(dcosterm_dr);
  Deallocate1DArray(eterm);
  Deallocate2DArray(determ_dr);
  delete [] map_t_desc_two;
  delete [] map_t_desc_three;

  Deallocate2DArray(dGCdr_two);
  Deallocate3DArray(dGCdr_three);

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

#endif  // ANN_IMPLEMENTATION_HPP_
