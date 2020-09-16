// This file was created by H Brooks on 29/06/2020
#include "DagMCBase.hpp"

using namespace DAGMC;

// *****************************************************************************
// Public methods
// *****************************************************************************

// get the float verision of dagmc version string
float DagMCBase::version(std::string* version_string) {
  if (NULL != version_string)
    *version_string = std::string("DagMC version ") + std::string(DAGMC_VERSION_STRING);
  return DAGMC_VERSION;
}

