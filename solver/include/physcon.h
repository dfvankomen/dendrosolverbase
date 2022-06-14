#ifndef DENDROSOLVER_PHYSCON_H_
#define DENDROSOLVER_PHYSCON_H_

#include <iostream>

#include "derivs.h"
#include "grUtils.h"
#include "parameters.h"

void physical_constraints(double **uZipConVars, const double **uZipVars,
                          const unsigned int &offset, const double *pmin,
                          const double *pmax, const unsigned int *sz,
                          const unsigned int &bflag);

#endif  // DENDROSOLVER_PHYSCON_H_