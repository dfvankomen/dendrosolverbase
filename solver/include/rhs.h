#ifndef DENDROSOLVER_RHS_H_
#define DENDROSOLVER_RHS_H_

#include <cmath>
#include <iostream>
#include <time.h>
#include "derivs.h"
#include "parameters.h"
#include "profile_params.h"
#include "grDef.h"
#include "mathUtils.h"
#include "block.h"

#ifdef DENDROSOLVER_ENABLE_CUDA
#include "rhs_cuda.cuh"
#include "params_ch.h"
#include "profile_gpu.h"
#endif

/**@brief computes complete RHS iteratiing over all the blocks.
 * @param[out] unzipVarsRHS: unzipped variables computed RHS
 * @param[in]  unzipVars: unzipped variables. 
 * @param[in]  blkList: block list. 
 * @param[in]  numBlocks: number of blocks. 
 */
void dendroSolverRHS(double **uzipVarsRHS, const double **uZipVars, const ot::Block *blkList, unsigned int numBlocks);

void dendroSolverRHSUnpacked(double **uzipVarsRHS, const double **uZipVars,
             const unsigned int &offset,
             const double *ptmin, const double *ptmax, const unsigned int *sz,
             const unsigned int &bflag);

// void dendroSolverRHSUnpacked_sep(double **uzipVarsRHS, const double **uZipVars,
//                  const unsigned int &offset,
//                  const double *ptmin, const double *ptmax, const unsigned int *sz,
//                  const unsigned int &bflag);

void dendroSolverBCS(double *f_rhs, const double *f,
              const double *dxf, const double *dyf, const double *dzf,
              const double *pmin, const double *pmax,
              const double f_falloff, const double f_asymptotic,
              const unsigned int *sz, const unsigned int &bflag);

void freeze_bcs(double *f_rhs, const unsigned int *sz, const unsigned int &bflag);

void fake_initial_data(double x, double y, double z, double *u);

void max_spacetime_speeds(double *const lambda1max, double *const lambda2max, double *const lambda3max,
                          const double *const alpha,
                          const double *const beta1, const double *const beta2, const double *const beta3,
                          const double *const gtd11, const double *const gtd12, const double *const gtd13,
                          const double *const gtd22, const double *const gtd23, const double *const gtd33,
                          const double *const chi, const unsigned int *sz);

void call_HAD_rhs();

#endif // DENDROSOLVER_RHS_H_