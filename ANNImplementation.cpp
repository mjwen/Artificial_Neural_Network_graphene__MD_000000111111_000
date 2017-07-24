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
// Copyright (c) 2013--2015, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Mingjian Wen
//


#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include "ANN.hpp"
#include "ANNImplementation.hpp"
#include "KIM_API_status.h"

#define MAXLINE 1024
#define IGNORE_RESULT(fn) if(fn){}


//==============================================================================
//
// Implementation of ANNImplementation public member functions
//
//==============================================================================

//******************************************************************************
ANNImplementation::ANNImplementation(
    KIM_API_model* const pkim,
    char const* const  parameterFileNames,
    int const parameterFileNameLength,
    int const numberParameterFiles,
    int* const ier)
    : numberOfSpeciesIndex_(-1),  // initizlize index, pointer, and cached
      numberOfParticlesIndex_(-1),    // member variables
      particleSpeciesIndex_(-1),
      coordinatesIndex_(-1),
      get_neighIndex_(-1),
      process_dEdrIndex_(-1),
      process_d2Edr2Index_(-1),
      cutoffIndex_(-1),
      energyIndex_(-1),
      forcesIndex_(-1),
      particleEnergyIndex_(-1),
      numberModelSpecies_(0),
      numberUniqueSpeciesPairs_(0),
      cutoffs_(0),
      cachedNumberOfParticles_(0),
      cachedNumberContributingParticles_(0)
			// add potential parameters


{
  *ier = SetConstantValues(pkim);
  if (*ier < KIM_STATUS_OK) return;

  AllocateFreeParameterMemory();

  FILE* parameterFilePointers[MAX_PARAMETER_FILES];
  *ier = OpenParameterFiles(pkim, parameterFileNames, parameterFileNameLength,
                            numberParameterFiles, parameterFilePointers);
  if (*ier < KIM_STATUS_OK) return;

//TODO read file
  *ier = ProcessParameterFiles(pkim, parameterFilePointers, numberParameterFiles);

  CloseParameterFiles(parameterFilePointers, numberParameterFiles);
  if (*ier < KIM_STATUS_OK) return;

//TODO enable later
//  *ier = ConvertUnits(pkim);
//  if (*ier < KIM_STATUS_OK) return;

  *ier = SetReinitMutableValues(pkim);
  if (*ier < KIM_STATUS_OK) return;

  *ier = RegisterKIMParameters(pkim);
  if (*ier < KIM_STATUS_OK) return;

  *ier = RegisterKIMFunctions(pkim);
  if (*ier < KIM_STATUS_OK) return;

  // everything is good
  *ier = KIM_STATUS_OK;
  return;
}

//******************************************************************************
ANNImplementation::~ANNImplementation()
{ // note: it is ok to delete a null pointer and we have ensured that
  // everything is initialized to null

  delete [] cutoffs_;
}

//******************************************************************************
int ANNImplementation::Reinit(KIM_API_model* const pkim)
{
  int ier;

  ier = SetReinitMutableValues(pkim);
  if (ier < KIM_STATUS_OK) return ier;

  // nothing else to do for this case

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int ANNImplementation::Compute(KIM_API_model* const pkim)
{
  int ier;

  // KIM API Model Input compute flags
  bool isComputeProcess_dEdr;
  bool isComputeProcess_d2Edr2;
  //
  // KIM API Model Output compute flags
  bool isComputeEnergy;
  bool isComputeForces;
  bool isComputeParticleEnergy;
  //
  // KIM API Model Input
  int const* particleSpecies = 0;
  GetNeighborFunction * get_neigh = 0;
  VectorOfSizeDIM const* coordinates = 0;
  //
  // KIM API Model Output
  double* energy = 0;
  double* particleEnergy = 0;
  VectorOfSizeDIM* forces = 0;
  ier = SetComputeMutableValues(pkim, isComputeProcess_dEdr,
                                isComputeProcess_d2Edr2, isComputeEnergy,
                                isComputeForces, isComputeParticleEnergy,
                                particleSpecies, get_neigh,
                                coordinates, energy, particleEnergy, forces);
  if (ier < KIM_STATUS_OK) return ier;

  // Skip this check for efficiency
  //
  // ier = CheckParticleSpecies(pkim, particleSpecies);
  // if (ier < KIM_STATUS_OK) return ier;


#include "ANNImplementationComputeDispatch.cpp"
  return ier;
}

//==============================================================================
//
// Implementation of ANNImplementation private member functions
//
//==============================================================================

//******************************************************************************
int ANNImplementation::SetConstantValues(KIM_API_model* const pkim)
{
  int ier = KIM_STATUS_FAIL;

  // get baseconvert value from KIM API object
  baseconvert_ = pkim->get_model_index_shift();

  // obtain indices for various KIM API Object arguments
  pkim->getm_index(
      &ier, 3 * 11,
      "numberOfSpecies",             &numberOfSpeciesIndex_,             1,
      "numberOfParticles",           &numberOfParticlesIndex_,           1,
      "particleSpecies",             &particleSpeciesIndex_,             1,
      "coordinates",                 &coordinatesIndex_,                 1,
      "get_neigh",                   &get_neighIndex_,                   1,
      "process_dEdr",                &process_dEdrIndex_,                1,
      "process_d2Edr2",              &process_d2Edr2Index_,              1,
      "cutoff",                      &cutoffIndex_,                      1,
      "energy",                      &energyIndex_,                      1,
      "forces",                      &forcesIndex_,                      1,
      "particleEnergy",              &particleEnergyIndex_,              1);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "getm_index", ier);
    return ier;
  }

  // set numberModelSpecies & numberUniqueSpeciesPairs
  int dummy;
  ier = pkim->get_num_model_species(&numberModelSpecies_, &dummy);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "get_num_model_species", ier);
    return ier;
  }
  numberUniqueSpeciesPairs_ = ((numberModelSpecies_+1)*numberModelSpecies_)/2;

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int ANNImplementation::OpenParameterFiles(
    KIM_API_model* const pkim,
    char const* const parameterFileNames,
    int const parameterFileNameLength,
    int const numberParameterFiles,
    FILE* parameterFilePointers[MAX_PARAMETER_FILES])
{
  int ier;

  if (numberParameterFiles > MAX_PARAMETER_FILES)
  {
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__, "ANN given too many"
                       " parameter files", ier);
    return ier;
  }

  for (int i = 0; i < numberParameterFiles; ++i)
  {
    parameterFilePointers[i]
        = fopen(&parameterFileNames[i * parameterFileNameLength], "r");
    if (parameterFilePointers[i] == 0)
    {
      char message[MAXLINE];
      sprintf(message,
              "ANN parameter file number %d cannot be opened",
              i);
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__, message, ier);
      for (int j = i - 1; i <= 0; --i)
      {
        fclose(parameterFilePointers[j]);
      }
      return ier;
    }
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int ANNImplementation::ProcessParameterFiles(
    KIM_API_model* const pkim,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
    int const numberParameterFiles)
{
  int N, ier;
  int endOfFileFlag = 0;
  char spec1[MAXLINE], spec2[MAXLINE], nextLine[MAXLINE];
  char *nextLinePtr;
  int iIndex, jIndex , indx, iiIndex, jjIndex;
  double nextCutoff, nextEpsilon, nextSigma;
  nextLinePtr = nextLine;

/*
  getNextDataLine(parameterFilePointers[0], nextLinePtr,
                  MAXLINE, &endOfFileFlag);
  ier = sscanf(nextLine, "%d %d", &N);
  if (ier != 2)
  {
    sprintf(nextLine, "unable to read first line of the parameter file");
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__, nextLine, ier);
    fclose(parameterFilePointers[0]);
    return ier;
  }
  if (N != numberModelSpecies_)
  {
    sprintf(nextLine, "The value for N from the parameter file is inconsistent "
            "with numberModelSpecies_");
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__, nextLine, ier);
    fclose(parameterFilePointers[0]);
    return ier;
  }

  // get and correctly order the particle names
  const char** const particleNames = new const char*[numberModelSpecies_];
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    const char* kimModelParticleSpecies;
    ier = pkim->get_model_species(i, &kimModelParticleSpecies);
    if (ier < KIM_STATUS_OK)
    {
      pkim->report_error(__LINE__, __FILE__, "get_model_species", ier);
      delete [] particleNames;
      return ier;
    }
    int const index = pkim->get_species_code(kimModelParticleSpecies, &ier);
    if (index >= numberModelSpecies_)
    {
      pkim->report_error(__LINE__, __FILE__, "get_species_code",
                         KIM_STATUS_FAIL);
      delete [] particleNames;
      return KIM_STATUS_FAIL;
    }
    particleNames[index] = kimModelParticleSpecies;
  }

  // set all values in the arrays to -1 for mixing later
  for (int i = 0; i < ((N+1)*N/2); i++)
  {
    cutoffs_[i]  = -1;
    epsilons_[i] = -1;
    sigmas_[i] = -1;
  }

  // Read and process data lines
  getNextDataLine(parameterFilePointers[0], nextLinePtr,
                  MAXLINE, &endOfFileFlag);
  while (endOfFileFlag == 0)
  {
    ier = sscanf(nextLine, "%s  %s %lg %lg %lg",
	         spec1, spec2, &nextCutoff, &nextEpsilon, &nextSigma);
    if (ier != 5)
    {
      sprintf(nextLine, "error reading lines of the parameter file");
      pkim->report_error(__LINE__, __FILE__, nextLine, KIM_STATUS_FAIL);
      delete [] particleNames;
      return KIM_STATUS_FAIL;
    }
    iIndex = jIndex = -1;
    for (int i = 0; i <  N; i++)
    {
      if (strcmp(spec1, particleNames[i]) == 0)
      {
        iIndex = i;
      }
      if (strcmp(spec2, particleNames[i]) == 0)
      {
        jIndex = i;
      }
    }
    if ((iIndex == -1) || (jIndex == -1))
    {
      sprintf(nextLine, "Unsupported Species name found in parameter file");
      pkim->report_error(__LINE__, __FILE__, nextLine, KIM_STATUS_FAIL);
      delete [] particleNames;
      return KIM_STATUS_FAIL;
    }
    if (iIndex >= jIndex)
    {
      indx = jIndex*N + iIndex - (jIndex*jIndex + jIndex)/2;
    }
    else
    {
      indx = iIndex*N + jIndex - (iIndex*iIndex + iIndex)/2;
    }
    cutoffs_[indx] = nextCutoff;
    epsilons_[indx] = nextEpsilon;
    sigmas_[indx] = nextSigma;

    getNextDataLine(parameterFilePointers[0], nextLinePtr,
                    MAXLINE, &endOfFileFlag);
  }

  // check that we got all like - like pairs
  sprintf(nextLine, "There are not values for like-like pairs of:");
  for (int i = 0; i < N; i++)
  {
    if (cutoffs_[(i*N + i - (i*i + i)/2)] == -1)
    {
      strcat(nextLine, "  ");
      strcat(nextLine, particleNames[i]);
      ier = -1;
    }
  }
  if (ier == -1)
  {
    pkim->report_error(__LINE__, __FILE__, nextLine, KIM_STATUS_FAIL);
    delete [] particleNames;
    return KIM_STATUS_FAIL;
  }

  // Perform Mixing if nessisary
  for (int jIndex = 0; jIndex < N; jIndex++)
  {
    jjIndex = (jIndex*N + jIndex - (jIndex*jIndex + jIndex)/2);
    for (int iIndex = (jIndex+1) ; iIndex < N; iIndex++)
    {
      indx = jIndex*N + iIndex - (jIndex*jIndex + jIndex)/2;
      if (cutoffs_[indx] == -1)
      {
        iiIndex = (iIndex*N + iIndex - (iIndex*iIndex + iIndex)/2);
        epsilons_[indx] = sqrt(epsilons_[iiIndex]*epsilons_[jjIndex]);
        sigmas_[indx] = (sigmas_[iiIndex] + sigmas_[jjIndex])/2.0;
        cutoffs_[indx] = (cutoffs_[iiIndex] + cutoffs_[jjIndex])/2.0;
      }
    }
  }
  delete [] particleNames;

*/

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
void ANNImplementation::getNextDataLine(
    FILE* const filePtr, char* nextLinePtr, int const maxSize,
    int *endOfFileFlag)
{
  do
  {
    if(fgets(nextLinePtr, maxSize, filePtr) == NULL)
    {
       *endOfFileFlag = 1;
       break;
    }
    while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t') ||
           (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r' ))
    {
      nextLinePtr = (nextLinePtr + 1);
    }
  }
  while ((strncmp("#", nextLinePtr, 1) == 0) || (strlen(nextLinePtr) == 0));
}

//******************************************************************************
void ANNImplementation::CloseParameterFiles(
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
    int const numberParameterFiles)
{
  for (int i = 0; i < numberParameterFiles; ++i)
    fclose(parameterFilePointers[i]);
}

//******************************************************************************
void ANNImplementation::AllocateFreeParameterMemory()
{ // allocate memory for data
  cutoffs_ = new double[numberUniqueSpeciesPairs_];

}

//******************************************************************************
int ANNImplementation::ConvertUnits(KIM_API_model* const pkim)
{
  int ier;

  // define default base units
  char length[] = "A";
  char energy[] = "eV";
  char charge[] = "e";
  char temperature[] = "K";
  char time[] = "ps";

/*
  // changing units of cutoffs and sigmas
  double const convertLength
      = pkim->convert_to_act_unit(length, energy, charge, temperature, time,
                                  1.0, 0.0, 0.0, 0.0, 0.0, &ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "convert_to_act_unit", ier);
    return ier;
  }
  if (convertLength != ONE)
  {
    for (int i = 0; i < numberUniqueSpeciesPairs_; ++i)
    {
      cutoffs_[i] *= cutoffs_[i];  // convert to active units
      sigmas_[i] *= sigmas_[i];  // convert to active units
    }
  }
  // changing units of epsilons
  double const convertEnergy
      = pkim->convert_to_act_unit(length, energy, charge, temperature, time,
                                  0.0, 1.0, 0.0, 0.0, 0.0, &ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "convert_to_act_unit", ier);
    return ier;
  }
  if (convertEnergy != ONE)
  {
    for (int i = 0; i < numberUniqueSpeciesPairs_; ++i)
    {
      epsilons_[i] *= convertEnergy;  // convert to active units
    }
  }

*/

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int ANNImplementation::RegisterKIMParameters(
    KIM_API_model* const pkim) const
{
  int ier;

  // publish parameters
/*  pkim->setm_data(&ier, 1 * 4,
                  //
                  "PARAM_FREE_cutoffs",
                  numberUniqueSpeciesPairs_,
                  (void*) cutoffs_,
                  1);
  if (ier < KIM_STATUS_OK) {
    pkim->report_error(__LINE__, __FILE__, "setm_data", ier);
    return ier;
  }
*/
  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int ANNImplementation::RegisterKIMFunctions(
    KIM_API_model* const pkim)
    const
{
  int ier;

  // register the destroy() and reinit() functions
  pkim->setm_method(&ier, 3 * 4,
                    "destroy", 1, (func_ptr) &(ANN::Destroy), 1,
                    "reinit",  1, (func_ptr) &(ANN::Reinit),  1,
                    "compute", 1, (func_ptr) &(ANN::Compute), 1);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "setm_method", ier);
    return ier;
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int ANNImplementation::SetReinitMutableValues(
    KIM_API_model* const pkim)
{ // use (possibly) new values of free parameters to compute other quantities
  int ier;

  // get cutoff pointer
  double* const cutoff
      = static_cast<double*>(pkim->get_data_by_index(cutoffIndex_, &ier));
  if (ier < KIM_STATUS_OK) {
    pkim->report_error(__LINE__, __FILE__, "get_data_by_index", ier);
    return ier;
  }

  // update cutoff value in KIM API object
  *cutoff = 0;
  int numberSpecies, maxStringLength;
  ier = pkim->get_num_sim_species(&numberSpecies, &maxStringLength);
  if (ier < KIM_STATUS_OK) {
    pkim->report_error(__LINE__, __FILE__, "get_num_sim_species", ier);
    return ier;
  }


  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int ANNImplementation::SetComputeMutableValues(
    KIM_API_model* const pkim,
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
    VectorOfSizeDIM*& forces)
{
  int ier = KIM_STATUS_FAIL;

  // get compute flags
  int compEnergy;
  int compForces;
  int compParticleEnergy;
  int compProcess_dEdr;
  int compProcess_d2Edr2;
  pkim->getm_compute_by_index(&ier, 3 * 5,
                              energyIndex_,         &compEnergy,         1,
                              forcesIndex_,         &compForces,         1,
                              particleEnergyIndex_, &compParticleEnergy, 1,
                              process_dEdrIndex_,   &compProcess_dEdr,   1,
                              process_d2Edr2Index_, &compProcess_d2Edr2, 1);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "getm_compute_by_index", ier);
    return ier;
  }

  isComputeEnergy = (compEnergy == KIM_COMPUTE_TRUE);
  isComputeForces = (compForces == KIM_COMPUTE_TRUE);
  isComputeParticleEnergy = (compParticleEnergy == KIM_COMPUTE_TRUE);
  isComputeProcess_dEdr = (compProcess_dEdr == KIM_COMPUTE_TRUE);
  isComputeProcess_d2Edr2 = (compProcess_d2Edr2 == KIM_COMPUTE_TRUE);

  // extract pointers based on compute flags
  //
  // double const* cutoff;            // currently unused
  // int const* numberOfSpecies;  // currently unused
  int const* numberOfParticles;
  pkim->getm_data_by_index(
      &ier, 3 * 6,
      // cutoffIndex_, &cutoff, 1,
      // numberOfSpeciesIndex_, &numberOfSpecies, 1,
      numberOfParticlesIndex_, &numberOfParticles, 1,
      particleSpeciesIndex_, &particleSpecies, 1,
      coordinatesIndex_, &coordinates, 1,
      energyIndex_, &energy, compEnergy,
      particleEnergyIndex_, &particleEnergy, compParticleEnergy,
      forcesIndex_, &forces, compForces);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "getm_data_by_index", ier);
    return ier;
  }

  // update values
  cachedNumberOfParticles_ = *numberOfParticles;

	// set so that it can be used even with a full neighbor list
	cachedNumberContributingParticles_ = *numberOfParticles;


  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int ANNImplementation::CheckParticleSpecies(
    KIM_API_model* const pkim,
    int const* const particleSpecies)
    const
{
  int ier;
  for (int i = 0; i < cachedNumberOfParticles_; ++i)
  {
    if ((particleSpecies[i] < 0) || (particleSpecies[i] >= numberModelSpecies_))
    {
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__,
                         "unsupported particle species detected", ier);
      return ier;
    }
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int ANNImplementation::GetComputeIndex(
    const bool& isComputeProcess_dEdr,
    const bool& isComputeProcess_d2Edr2,
    const bool& isComputeEnergy,
    const bool& isComputeForces,
    const bool& isComputeParticleEnergy) const
{
  const int processdE = 2;
  const int processd2E = 2;
  const int energy = 2;
  const int force = 2;
  const int particleEnergy = 2;


  int index = 0;

  // processdE
  index += (int(isComputeProcess_dEdr))
      * processd2E * energy * force * particleEnergy;

  // processd2E
  index += (int(isComputeProcess_d2Edr2)) * energy * force * particleEnergy;

  // energy
  index += (int(isComputeEnergy)) * force * particleEnergy;

  // force
  index += (int(isComputeForces)) * particleEnergy;

  // particleEnergy
  index += (int(isComputeParticleEnergy));


  return index;
}


//==============================================================================
//
// Implementation of helper functions
//
//==============================================================================

//******************************************************************************
void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
                                  int const extentOne)
{ // allocate memory and set pointers
  arrayPtr = new double*[extentZero];
  arrayPtr[0] = new double[extentZero * extentOne];
  for (int i = 1; i < extentZero; ++i)
  {
    arrayPtr[i] = arrayPtr[i-1] + extentOne;
  }

  // initialize
  for (int i = 0; i < extentZero; ++i)
  {
    for (int j = 0; j < extentOne; ++j)
    {
      arrayPtr[i][j] = 0.0;
    }
  }
}

//******************************************************************************
void Deallocate2DArray(double**& arrayPtr)
{ // deallocate memory
  if (arrayPtr != 0) delete [] arrayPtr[0];
  delete [] arrayPtr;

  // nullify pointer
  arrayPtr = 0;
}
