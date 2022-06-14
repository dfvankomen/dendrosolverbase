/**
 * @author Milinda Fernando / David Van Komen
 * School of Computing, University of Utah
 * @brief Contins utility functions for EMDA simulation
 */

#ifndef DENDROSOLVER_GRUTILS_H_
#define DENDROSOLVER_GRUTILS_H_

#include "block.h"
#include "dendroProfileParams.h"
#include "grDef.h"
#include "json.hpp"
#include "lebedev.h"
#include "mesh.h"
#include "parUtils.h"
#include "parameters.h"
#include "point.h"
#include "profile_params.h"
#include "swsh.h"

using json = nlohmann::json;
namespace dsolve {
/**
 * @brief: Read the parameter file and initialize the variables in parameters.h
 * file.
 * @param[in] fName: file name
 * @param[in] comm: MPI communicator.
 * */
// void readParamFile(const char *fName, MPI_Comm comm);

/**
 * @brief dump the read parameter files.
 *
 * @param sout
 * @param root
 * @param comm
 */
// void dumpParamFile(std::ostream &sout, int root, MPI_Comm comm);

// INITIAL DATA FUNCTIONS
// clang-format off
/*[[[cog
import cog
import sys
import os
import importlib.util
import dendrosym

# get the current working directory, should be root of project
current_path = os.getcwd()
output_path = os.path.join(current_path, "gencode")

# the following lines will import any module directly from
spec = importlib.util.spec_from_file_location("dendroconf", CONFIG_FILE_PATH)
dendroconf = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = dendroconf
spec.loader.exec_module(dendroconf)

cog.outl("// INITIAL DATA FUNCTIONS")
cog.outl(dendroconf.dendroConfigs.generate_initial_data_declaration(var_type="evolution"))
]]]*/
// clang-format on

//[[[end]]]

/**
 * @brief Two puncture intiial data from HAD code/
 *
 * @param xx1 : x coord (octree coord)
 * @param yy1 : y coord (octree coord)
 * @param zz1 : z coord (octree coord)
 * @param var : initialized dsolve variables for the grid points
 */
void punctureData(const double xx1, const double yy1, const double zz1,
                  double *var);

/**
 * @brief  evaluate HAD initial data from the physical coordinates of the
 * domain. Uses domain coords no coord transformation inside.
 * @param xx1 x coord
 * @param yy1 y coord
 * @param zz1 z coord
 * @param var Initialized data eval at (x,y,z) point.
 */
void punctureDataPhysicalCoord(const double xx, const double yy,
                               const double zz, double *var);

/**
 * @brief compute the static Kerr-Schild BH data
 *
 * @param xx1 : x coord
 * @param yy1 : y coord
 * @param zz1 : z coord
 * @param var : initialized dsolve variables for the grid points
 */
void KerrSchildData(const double xx1, const double yy1, const double zz1,
                    double *var);

/**
 * @brief add artificial noise to the initial data.
 * @param xx1 : x coord
 * @param yy1 : y coord
 * @param zz1 : z coord
 * @param var : initialized dsolve variables for the grid points
 */
// void noiseData(const double xx1, const double yy1, const double zz1, double
// *var);

/**
 * @brief fake initial data.
 * @param xx1 : x coord
 * @param yy1 : y coord
 * @param zz1 : z coord
 * @param var : initialized dsolve variables for the grid points
 */
void fake_initial_data(double x, double y, double z, double *u);

// NOTE: the following functions were added when adding dsolve functionality
/**
 * @brief calculate and set the initial data for superposed boosted kerr-sen
 * @param xx1 : x coord, GRIDX format
 * @param yy1 : y coord, GRIDX format
 * @param zz1 : z coord, GRIDX format
 * @param var : initialized dsolve variables for the grid points
 */
void initDataFuncToPhysCoords(double xx1, double yy1, double zz1, double *var);

/**
 * @brief calculate and set the initial data for noise
 * @param xx1 : x coord
 * @param yy1 : y coord
 * @param zz1 : z coord
 * @param var : initialized dsolve variables for the grid points
 */
void noiseInit(double x, double y, double z, double *u);

namespace trumpet_data {

void trumpetData(const double xx1, const double yy1, const double zz1,
                 double *var);
void bndcnd(double h, double &x, double y[], double dydx[]);
void derivs(double x, double y[], double dydx[]);
void hunt(double xx[], int n, double x, int *jlo);
void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
          double yerr[], void (*derivs)(double, double[], double[]));
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
          double yscal[], double *hdid, double *hnext,
          void (*derivs)(double, double[], double[]));
void odeint(double ystart[], int nvar, double x1, double x2, double eps,
            double h1, double hmin, int *nok, int *nbad,
            void (*derivs)(double, double[], double[]),
            void (*rkqs)(double[], double[], int, double *, double, double,
                         double[], double *, double *,
                         void (*)(double, double[], double[])),
            int kount);
double interpolation3(double xp[], double yp[], int np, double xb,
                      int *n_nearest_pt);
double interpolation4(double xp[], double yp[], int np, double xb,
                      int *n_nearest_pt);

}  // end of namespace trumpet_data

/**
 * @brief: Generates block adaptive octree for the given binary blockhole
 * problem.
 * @param[out] tmpNodes: created octree tmpNodes
 * @param[in] pt_min: block min point
 * @param[in] pt_max: block max point
 * @param[in] regLev: regular grid level
 * @param[in] maxDepth: maximum refinement level.
 * @param[in] comm: MPI communicator.
 * */
void blockAdaptiveOctree(std::vector<ot::TreeNode> &tmpNodes,
                         const Point &pt_min, const Point &pt_max,
                         const unsigned int regLev, const unsigned int maxDepth,
                         MPI_Comm comm);

/**
 * @brief Compute the wavelet tolerance as a function of space.
 *
 * @param x : x coord.
 * @param y : y coord
 * @param z : z coord
 * @param tol_min : min. tolerance value.
 * @return double
 */
double computeWTol(double x, double y, double z, double tol_min);

/**
 * @brief Compute the wavelet tolerance as a function of space (uses actual
 * domain coordinates not octree coordinates)
 *
 * @param x : x coord.
 * @param y : y coord
 * @param z : z coord
 * @param hx : resolution in x,y,z
 * @return double
 */
double computeWTolDCoords(double x, double y, double z, double *hx);

/**
 * @breif: Compute L2 constraint norms.
 */
template <typename T>
double computeConstraintL2Norm(const T *constraintVec, const T *maskVec,
                               unsigned int lbegin, unsigned int lend,
                               MPI_Comm comm);

/**
 * @brief Compute L2 constraint norms.
 */
template <typename T>
double computeConstraintL2Norm(const ot::Mesh *mesh, const T *constraintVec,
                               const T *maskVector, T maskthreshoold);

/**
 * @breif write constraints to a file.
 */
template <typename T>
double extractConstraints(const ot::Mesh *mesh, const T **constraintVar,
                          const T *maskVec, double maskthreshoold,
                          unsigned int timestep, double stime);

/**@brief : write a block to binary*/
void writeBLockToBinary(const double **unzipVarsRHS, unsigned int offset,
                        const double *pmin, const double *pmax, double *bxMin,
                        double *bxMax, const unsigned int *sz,
                        unsigned int blkSz, double dxFactor,
                        const char *fprefix);

/**@brief returns the octant weight for LTS timestepping. */
unsigned int getOctantWeight(const ot::TreeNode *pNode);

/**
 * @brief Compute the BH locations based on the integration of the shift
 * vectors.
 *
 * @param in current locations of the BHs
 * @param out : evolved locations of the BH
 * @param zipVars : zip representation of the current evolution variables.
 * (assumes ghost synced)
 * @param dt : time step size.
 */
void computeBHLocations(const ot::Mesh *pMesh, const Point *in, Point *out,
                        double **zipVars, double dt);

/**
 * @brief Allocate the derivative workspace for use in RHS functionality
 *
 * @param pMesh
 * @param s_fac
 */
void allocate_deriv_workspace(const ot::Mesh *pMesh, unsigned int s_fac);

/**
 * @brief Deallocate the derivative workspace for use in the RHS functionality
 *
 */
void deallocate_deriv_workspace();

}  // end of namespace dsolve

namespace dsolve {

namespace timer {

/**@brief initialize all the flop counters. */
void initFlops();

/**@brief clears the snapshot counter for time profiler variables*/
void resetSnapshot();

/**@brief reduce min mean max.
 * @param [in] stat: local time
 * @param [out] stat_g 0-min, 1-mean 2-max
 * */
template <typename T>
void computeOverallStats(T *stat, T *stat_g, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    par::Mpi_Reduce(stat, stat_g, 1, MPI_MIN, 0, comm);
    par::Mpi_Reduce(stat, stat_g + 1, 1, MPI_SUM, 0, comm);
    par::Mpi_Reduce(stat, stat_g + 2, 1, MPI_MAX, 0, comm);
    stat_g[1] /= (npes);
}

/** @breif : printout the profile parameters. */
void profileInfo(const char *filePrefix, const ot::Mesh *pMesh);

/** @breif : printout the profile parameters (intermediate profile information).
 */
void profileInfoIntermediate(const char *filePrefix, const ot::Mesh *pMesh,
                             const unsigned int currentStep);

}  // namespace timer

}  // namespace dsolve

namespace GW {
/**
 * @brief : debug function to write psi4 to interpolated to spheres.
 */
void psi4ShpereDump(const ot::Mesh *mesh, DendroScalar **cVar,
                    unsigned int timestep, double time);

}  // end of namespace GW

#include "grUtils.tcc"

#endif  // DENDROSOLVER_GRUTILS_H_
