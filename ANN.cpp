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
//    Mingjian
//


#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include "ANN.hpp"
#include "ANNImplementation.hpp"
#include "KIM_API_status.h"

//==============================================================================
//
// This is the standard interface to KIM Model Drivers
//
//==============================================================================

//******************************************************************************
int model_driver_init(void* km, char* paramfile_names, int* nmstrlen,
                      int* numparamfiles)
{
  KIM_API_model* const pkim = *static_cast<KIM_API_model**>(km);
  int ier;

  // read input files, convert units if needed, compute
  // interpolation coefficients, set cutoff, and publish parameters
  ANN* const modelObject
      = new ANN(pkim, paramfile_names, *nmstrlen,
                            *numparamfiles, &ier);
  if (ier < KIM_STATUS_OK)
  {
    // constructor already reported the error
    delete modelObject;
    return ier;
  }

  // register pointer to ANN object in KIM object
  pkim->set_model_buffer(static_cast<void*>(modelObject), &ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "set_model_buffer", ier);
    delete modelObject;
    return ier;
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//==============================================================================
//
// Implementation of ANN public wrapper functions
//
//==============================================================================

//******************************************************************************
ANN::ANN(
    KIM_API_model* const pkim,
    char const* const parameterFileNames,
    int const parameterFileNameLength,
    int const numberParameterFiles,
    int* const ier)
{
  implementation_ = new ANNImplementation(
      pkim,
      parameterFileNames,
      parameterFileNameLength,
      numberParameterFiles,
      ier);
}

//******************************************************************************
ANN::~ANN()
{
  delete implementation_;
}

//******************************************************************************
int ANN::Destroy(void* kimmdl)  // static member function
{
  KIM_API_model* const pkim = *static_cast<KIM_API_model**>(kimmdl);
  int ier;
  ANN* const modelObject
      = (ANN*) pkim->get_model_buffer(&ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "get_model_buffer", ier);
    return ier;
  }

  if (modelObject != NULL)
  {
    // delete object itself
    delete modelObject;

    // nullify model buffer
    pkim->set_model_buffer(NULL, &ier);
  }

  // everything is good
  ier = KIM_STATUS_OK;
  return ier;
}

//******************************************************************************
int ANN::Reinit(void* kimmdl)  // static member function
{
  KIM_API_model* const pkim = *static_cast<KIM_API_model**>(kimmdl);
  int ier;
  ANN* const modelObject
      = (ANN*) pkim->get_model_buffer(&ier);
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "get_model_buffer", ier);
    return ier;
  }

  return modelObject->implementation_->Reinit(pkim);
}

//******************************************************************************
int ANN::Compute(void* kimmdl)  // static member function
{
  KIM_API_model* const pkim = *static_cast<KIM_API_model**>(kimmdl);
  int ier;
  ANN* const modelObject
      = static_cast<ANN*>(pkim->get_model_buffer(&ier));
  if (ier < KIM_STATUS_OK)
  {
    pkim->report_error(__LINE__, __FILE__, "get_model_buffer", ier);
    return ier;
  }

  return modelObject->implementation_->Compute(pkim);
}
