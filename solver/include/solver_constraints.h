#ifndef DENDROSOLVER_SOLVER_CONSTRAINTS_H_
#define DENDROSOLVER_SOLVER_CONSTRAINTS_H_

#include <iostream>

#include "grUtils.h"
#include "parameters.h"

void enforce_solver_constraints(double **uiVar, const unsigned int offset);

#endif