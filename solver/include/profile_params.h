//
// Created by David
//
/**
 * @author Milinda Fernando / David Van Komen
 * School of Computing, University of Utah
 * @brief Contins parameters for the dsolve timing
 */

#ifndef DENDROSOLVER_PROFILE_PARAMS_H_
#define DENDROSOLVER_PROFILE_PARAMS_H_

#include "profiler.h"

namespace dsolve {
namespace timer {
extern profiler_t total_runtime;

extern profiler_t t_f2o;
extern profiler_t t_cons;
extern profiler_t t_bal;
extern profiler_t t_mesh;

extern profiler_t t_rkSolve;
extern profiler_t t_ghostEx_sync;

extern profiler_t t_unzip_sync;
extern profiler_t t_unzip_async;

extern profiler_t t_deriv;
extern profiler_t t_rhs;

// TODO: profilers for other RHS variables? Generate these if going for separate
// calculations
extern profiler_t t_rhs_a;
extern profiler_t t_rhs_b;
extern profiler_t t_rhs_gt;
extern profiler_t t_rhs_chi;
extern profiler_t t_rhs_At;
extern profiler_t t_rhs_K;
extern profiler_t t_rhs_Gt;
extern profiler_t t_rhs_B;

extern profiler_t t_bdyc;

extern profiler_t t_zip;
extern profiler_t t_rkStep;

extern profiler_t t_isReMesh;
extern profiler_t t_gridTransfer;
extern profiler_t t_ioVtu;
extern profiler_t t_ioCheckPoint;

}  // namespace timer
}  // namespace dsolve

#endif  // DENDROSOLVER_PROFILE_PARAMS_H_
