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

#include "KIM_API_status.h"
#include "ANN.hpp"
#include "ANNImplementation.hpp"
#include "descriptor.h"
#include "helper.h"

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
			particleStatusIndex_(-1),
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
			cutoffsSq2D_(0),
      cachedNumberOfParticles_(0),
      cachedNumberContributingParticles_(0)
			// add potential parameters


{
	// create descriptor and network classes
	descriptor_ = new Descriptor();
	network_ = new NeuralNetwork();

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
  Deallocate2DArray(cutoffsSq2D_);
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
      &ier, 3 * 12,
      "numberOfSpecies",             &numberOfSpeciesIndex_,             1,
      "numberOfParticles",           &numberOfParticlesIndex_,           1,
      "particleSpecies",             &particleSpeciesIndex_,             1,
      "particleStatus",							 &particleStatusIndex_,              1,
      "coordinates",                 &coordinatesIndex_,                 1,
      "get_neigh",                   &get_neighIndex_,                   1,
      "process_dEdr",                &process_dEdrIndex_,                1,
      "process_d2Edr2",              &process_d2Edr2Index_,              1,
      "cutoff",                      &cutoffIndex_,                      1,
      "energy",                      &energyIndex_,                      1,
      "forces",                      &forcesIndex_,                      1,
      "particleEnergy",              &particleEnergyIndex_,              1);
  if (ier < KIM_STATUS_OK) {
    pkim->report_error(__LINE__, __FILE__, "getm_index", ier);
    return ier;
  }

  // set numberModelSpecies & numberUniqueSpeciesPairs
  int dummy;
  ier = pkim->get_num_model_species(&numberModelSpecies_, &dummy);
  if (ier < KIM_STATUS_OK) {
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
  int ier;
  //int N;
  int endOfFileFlag = 0;
  char nextLine[MAXLINE];
  char errorMsg[MAXLINE];
  char name[MAXLINE];
	double cutoff;

  // descriptor
	int numDescTypes;
	int numDescs;
	int numParams;
	int numParamSets;
	double** descParams;

  // network
  int numLayers;
  int* numPerceptrons;
	double** weight;
	double* bias;

  //char spec1[MAXLINE], spec2[MAXLINE];
  //int iIndex, jIndex , indx, iiIndex, jjIndex;
  //double nextCutoff;

	// cutoff
  getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  ier = sscanf(nextLine, "%s %lf", name, &cutoff);
  if (ier != 2) {
    sprintf(errorMsg, "unable to read cutoff from line:\n");
    strcat(errorMsg, nextLine);
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
    fclose(parameterFilePointers[0]);
    return ier;
  }
	descriptor_->set_cutfunc(name);
	cutoffs_[0] = cutoff;

	// number of descriptor types
  getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  ier = sscanf(nextLine, "%d", &numDescTypes);
  if (ier != 1) {
    sprintf(errorMsg, "unable to read number of descriptor types from line:\n");
    strcat(errorMsg, nextLine);
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
    fclose(parameterFilePointers[0]);
    return ier;
  }

  // descriptor
  for (int i=0; i<numDescTypes; i++) {
    // descriptor name and parameter dimensions
    getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
    ier = sscanf(nextLine, "%s %d %d", name, &numParamSets, &numParams);
    if (ier != 3) {
      sprintf(errorMsg, "unable to read descriptor from line:\n");
      strcat(errorMsg, nextLine);
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
      fclose(parameterFilePointers[0]);
      return ier;
    }

    // check size of params is correct w.r.t its name
    if (strcmp(name, "G2") == 0) {
      if (numParams != 2) {
        sprintf(errorMsg, "number of params for descriptor G2 is incorrect, "
            "expecting 2, but given %d.\n", numParams);
        ier = KIM_STATUS_FAIL;
        pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
        fclose(parameterFilePointers[0]);
        return ier;
      }
    }
    else if (strcmp(name, "G3") == 0) {
      if (numParams != 1) {
        sprintf(errorMsg, "number of params for descriptor G3 is incorrect, "
            "expecting 1, but given %d.\n", numParams);
        ier = KIM_STATUS_FAIL;
        pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
        fclose(parameterFilePointers[0]);
        return ier;
      }
    }
    else {
      sprintf(errorMsg, "unsupported descriptor `%s' from line:\n", name);
      strcat(errorMsg, nextLine);
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
      fclose(parameterFilePointers[0]);
      return ier;
    }

    // read descriptor params
    AllocateAndInitialize2DArray(descParams, numParamSets, numParams);
    for (int j=0; j<numParamSets; j++) {
      getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
      ier = getXdouble(nextLine, numParams, descParams[j]);
      if (ier != KIM_STATUS_OK) {
        sprintf(errorMsg, "unable to read descriptor parameters from line:\n");
        strcat(errorMsg, nextLine);
        ier = KIM_STATUS_FAIL;
        pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
        fclose(parameterFilePointers[0]);
        return ier;
      }
    }

    // copy data to Descriptor
    descriptor_->add_descriptor(name, descParams, numParamSets, numParams);
    Deallocate2DArray(descParams);
  }
//TODO delete
  descriptor_->echo_input();


  // network structure
  // number of layers
  getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  ier = sscanf(nextLine, "%d", &numLayers);
  if (ier != 1) {
    sprintf(errorMsg, "unable to read number of layers from line:\n");
    strcat(errorMsg, nextLine);
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
    fclose(parameterFilePointers[0]);
    return ier;
  }
  // number of perceptrons in each layer
  numPerceptrons = new int[numLayers];
  getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  ier = getXint(nextLine, numLayers, numPerceptrons);
  if (ier != KIM_STATUS_OK) {
    sprintf(errorMsg, "unable to read number of perceptrons from line:\n");
    strcat(errorMsg, nextLine);
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
    fclose(parameterFilePointers[0]);
    return ier;
  }
  // number of descriptors, (size of input)
  numDescs = descriptor_->get_num_descriptors();
  // copy to network class
  network_->set_nn_structure(numDescs, numLayers, numPerceptrons);


  // activation function
  getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  ier = sscanf(nextLine, "%s", name);
  if (ier != 1) {
    sprintf(errorMsg, "unable to read activation function from line:\n");
    strcat(errorMsg, nextLine);
    ier = KIM_STATUS_FAIL;
    pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
    fclose(parameterFilePointers[0]);
    return ier;
  }
  network_->set_activation(name);


  // weights and biases
  for (int i=0; i<numLayers; i++) {

    // weights
    int row;
    int col;
    if (i==0) {
      row = numDescs;
      col = numPerceptrons[i];
    } else {
      row = numPerceptrons[i-1];
      col = numPerceptrons[i];
    }

    AllocateAndInitialize2DArray(weight, row, col);
    for (int j=0; j<row; j++) {
      getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
      ier = getXdouble(nextLine, col, weight[j]);
      if (ier != KIM_STATUS_OK) {
        sprintf(errorMsg, "unable to read weight from line:\n");
        strcat(errorMsg, nextLine);
        ier = KIM_STATUS_FAIL;
        pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
        fclose(parameterFilePointers[0]);
        return ier;
      }
    }

    // bias
    AllocateAndInitialize1DArray(bias, col);
    getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
    ier = getXdouble(nextLine, col, bias);
    if (ier != KIM_STATUS_OK) {
      sprintf(errorMsg, "unable to read bias from line:\n");
      strcat(errorMsg, nextLine);
      ier = KIM_STATUS_FAIL;
      pkim->report_error(__LINE__, __FILE__, errorMsg, ier);
      fclose(parameterFilePointers[0]);
      return ier;
    }

    // copy to network class
    network_->add_weight_bias(weight, bias, i);
    Deallocate2DArray(weight);
    Deallocate1DArray(bias);
  }


  delete numPerceptrons;



  network_->echo_input();
exit(1);



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
  char* pch;

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

  // remove comments starting with `#' in a line
  pch = strchr(nextLinePtr, '#');
  if (pch != NULL) {
    *pch = '\0';
  }

}

//******************************************************************************
int ANNImplementation::getXdouble(char* linePtr, const int N, double* list)
{
  int ier;
  char * pch;
  char line[MAXLINE];
  int i = 0;

  strcpy(line, linePtr);
  pch = strtok(line, " \t");
  while (pch != NULL) {
    ier = sscanf(pch, "%lf", &list[i]);
    if (ier != 1) {
      ier = KIM_STATUS_FAIL;
      return ier;
    }
    pch = strtok(NULL, " \t");
    i += 1;
  }
  if (i != N) {
    ier = KIM_STATUS_FAIL;
    return ier;
  }

  ier = KIM_STATUS_OK;
  return ier;
}


//******************************************************************************
int ANNImplementation::getXint(char* linePtr, const int N, int* list)
{
  int ier;
  char * pch;
  char line[MAXLINE];
  int i = 0;

  strcpy(line, linePtr);
  pch = strtok(line, " \t");
  while (pch != NULL) {
    ier = sscanf(pch, "%d", &list[i]);
    if (ier != 1) {
      ier = KIM_STATUS_FAIL;
      return ier;
    }
    pch = strtok(NULL, " \t");
    i += 1;
  }
  if (i != N) {
    ier = KIM_STATUS_FAIL;
    return ier;
  }

  ier = KIM_STATUS_OK;
  return ier;
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
	AllocateAndInitialize2DArray(cutoffsSq2D_, numberModelSpecies_, numberModelSpecies_);
}

//******************************************************************************
int ANNImplementation::ConvertUnits(KIM_API_model* const pkim)
{
  int ier;

  // define default base units
/*  char length[] = "A";
  char energy[] = "eV";
  char charge[] = "e";
  char temperature[] = "K";
  char time[] = "ps";
*/

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
  if (ier < KIM_STATUS_OK) {
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

	// update cutoffsSq (This requires PECIES_001_NAME_STR needs to have code 0,
	// SPECIES_002_NAME_STR needs to have code 1 ... in .kim file)
	for (int i = 0; i < numberModelSpecies_; ++i) {
		for (int j = 0; j <= i ; ++j) {
			int const index = j*numberModelSpecies_ + i - (j*j + j)/2;
			cutoffsSq2D_[i][j] = cutoffsSq2D_[j][i] = (cutoffs_[index]*cutoffs_[index]);
		}
	}

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

	// find the largest cutoff of a subset of all the supported species of the model
	const char* simSpeciesI;
	const char* simSpeciesJ;
	for (int i = 0; i < numberSpecies; i++) {
		ier = pkim->get_sim_species(i, &simSpeciesI);
		if (ier < KIM_STATUS_OK) {
			pkim->report_error(__LINE__, __FILE__, "get_num_sim_species", ier);
			return ier;
		}
		int const indexI = pkim->get_species_code(simSpeciesI, &ier);
		if (indexI >= numberModelSpecies_ || ier<KIM_STATUS_OK ) {
			pkim->report_error(__LINE__, __FILE__, "get_species_code",
					KIM_STATUS_FAIL);
			return KIM_STATUS_FAIL;
		}

		for (int j = 0; j < numberSpecies; j++) {
			ier = pkim->get_sim_species( j, &simSpeciesJ);
			if (ier < KIM_STATUS_OK) {
				pkim->report_error(__LINE__, __FILE__, "get_num_sim_species", ier);
				return ier;
			}
			int const indexJ = pkim->get_species_code(simSpeciesJ, &ier);
			if (indexJ >= numberModelSpecies_ || ier<KIM_STATUS_OK ) {
				pkim->report_error(__LINE__, __FILE__, "get_species_code",
						KIM_STATUS_FAIL);
				return KIM_STATUS_FAIL;
			}
			if (*cutoff < cutoffsSq2D_[indexI][indexJ]) {
				*cutoff = cutoffsSq2D_[indexI][indexJ];
			}
		}
	}
	*cutoff = sqrt(*cutoff);


//TODO update of other free params shoud go here


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
  if (ier < KIM_STATUS_OK) {
    pkim->report_error(__LINE__, __FILE__, "getm_compute_by_index", ier);
    return ier;
  }

  isComputeEnergy = (compEnergy == KIM_COMPUTE_TRUE);
  isComputeForces = (compForces == KIM_COMPUTE_TRUE);
  isComputeParticleEnergy = (compParticleEnergy == KIM_COMPUTE_TRUE);
  isComputeProcess_dEdr = (compProcess_dEdr == KIM_COMPUTE_TRUE);
  isComputeProcess_d2Edr2 = (compProcess_d2Edr2 == KIM_COMPUTE_TRUE);

  // extract pointers based on compute flags
  int const* numberOfParticles;
  int const* particleStatus = 0;
  pkim->getm_data_by_index(
      &ier, 3 * 7,
      numberOfParticlesIndex_, &numberOfParticles, 1,
      particleSpeciesIndex_,	 &particleSpecies,	 1,
      particleStatusIndex_,		 &particleStatus,		 1,
      coordinatesIndex_,			 &coordinates,			 1,
      energyIndex_,						 &energy,						 compEnergy,
      particleEnergyIndex_,		 &particleEnergy,		 compParticleEnergy,
      forcesIndex_,						 &forces,						 compForces);
  if (ier < KIM_STATUS_OK) {
    pkim->report_error(__LINE__, __FILE__, "getm_data_by_index", ier);
    return ier;
  }

	// get neigh function
	get_neigh = (GetNeighborFunction *) pkim->get_method_by_index(get_neighIndex_, &ier);
	if (ier < KIM_STATUS_OK) {
		pkim->report_error(__LINE__, __FILE__, "get_method_by_index", ier);
		return ier;
	}

  // update values
  cachedNumberOfParticles_ = *numberOfParticles;

	// set so that it can be used even with a full neighbor list
	cachedNumberContributingParticles_ = 0;
	for (int i=0; i<*numberOfParticles; i++) {
		if (particleStatus[i] == 1) {
			cachedNumberContributingParticles_ += 1;
		}
	}

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
  //const int processdE = 2;
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

