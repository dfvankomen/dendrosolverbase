

#ifndef DENDROSOLVER_PARAMETERS_H_
#define DENDROSOLVER_PARAMETERS_H_

// library includes
#include <string.h>

#include <iostream>
// toml needs to be in the path (or included via submodule)
#include <toml.hpp>

// dendro only includes
#include "dendro.h"
#include "memory_pool.h"
#include "parUtils.h"

// project-specific includes
#include "bh.h"
#include "grDef.h"

/*[[[cog

import cog
import dendrosym
import os

paramh_str, paramc_str = dendrosym.params.generate_all_parameter_text("dsolve",
PARAM_SETUP_FILE)

cog.outl(paramh_str)

# when that's done, we also generate the sample file, can't think of another
place to put this with open(os.path.join(os.path.dirname(PARAM_SETUP_FILE),
"solver_parameters.sample.toml"), "w") as f:
    f.write(dendrosym.params.generate_sample_config_file_text("dsolve",
PARAM_SETUP_FILE))

]]]*/
namespace dsolve {
void readParamFile(const char* inFile, MPI_Comm comm);
void dumpParamFile(std::ostream& sout, int root, MPI_Comm comm);

extern mem::memory_pool<double> DENDROSOLVER_MEM_POOL;
extern BH BH1;
extern BH BH2;
extern Point DENDROSOLVER_BH_LOC[2];
// extern RefinementMode DENDROSOLVER_REFINEMENT_MODE;

extern double* DENDROSOLVER_DERIV_WORKSPACE;
// number of derivatives, the greater between the RHS and Constraint
// TODO: this needs to be automated!!!!!!!!! ESPECIALLY WITH ADVANCED
// NOTE: THIS NUMBER IS PROBABLY TOO LARGE AS IS!
const unsigned int DENDROSOLVER_NUM_DERIVATIVES = 200;

}  // namespace dsolve

namespace GW {
extern unsigned int DENDROSOLVER_GW_NUM_RADAII;

extern unsigned int DENDROSOLVER_GW_NUM_LMODES;

extern unsigned int DENDROSOLVER_GW_RADAII[6];

extern unsigned int DENDROSOLVER_GW_L_MODES[3];

static const unsigned int DENDROSOLVER_GW_OUTPUT_PRECISION = 10;

}  // namespace GW

namespace BHLOC {
extern unsigned int EXTRACTION_VAR_ID;

extern double EXTRACTION_TOL;

}  // namespace BHLOC

namespace TPID {
extern std::string FILE_PREFIX;

static const double TP_epsilon = 1e-06;

static const int swap_xz = 0;

static const int use_sources = 0;

static const int rescale_sources = 0;

static const int use_external_initial_guess = 0;

static const int do_residuum_debug_output = 1;

static const int do_initial_debug_output = 1;

static const int multiply_old_lapse = 0;

static const double TP_Tiny = 1e-15;

static const double TP_Extend_Radius = 0.0;

static const int Newton_maxit = 5;

extern double par_b;

extern double par_m_plus;

extern double par_m_minus;

extern double target_M_plus;

extern double target_M_minus;

extern double par_P_plus[3];

extern double par_P_minus[3];

extern double par_S_plus[3];

extern double par_S_minus[3];

extern double center_offset[3];

extern unsigned int give_bare_mass;

extern int initial_lapse;

extern unsigned int grid_setup_method;

extern int solve_momentum_constraint;

extern unsigned int verbose;

extern double adm_tol;

extern double Newton_tol;

extern double initial_lapse_psi_exponent;

extern unsigned int npoints_A;

extern unsigned int npoints_B;

extern unsigned int npoints_phi;

}  // namespace TPID

namespace dsolve {
/** @brief: Dendro version number, usually 5.0 especially for this project */
static const double DENDRO_VERSION = 5.0;

/** @brief: Minimum possible value for all points of alpha */
static const double ALPHA_FLOOR = 0.1;

/** @brief: Element order for the computations */
extern unsigned int DENDROSOLVER_ELE_ORDER;

/** @brief: Padding width for each of the blocks */
extern unsigned int DENDROSOLVER_PADDING_WIDTH;

/** @brief: The number of total variables */
static const unsigned int DENDROSOLVER_NUM_VARS = 36;

/** @brief: Number of constraint variables */
static const unsigned int DENDROSOLVER_CONSTRAINT_NUM_VARS = 6;

/** @brief: Number of RK45 stages that should be performed */
static const unsigned int DENDROSOLVER_RK45_STAGES = 6;

/** @brief: Number of RK4 stages that should be performed */
static const unsigned int DENDROSOLVER_RK4_STAGES = 4;

/** @brief: Number of RK3 stages that should be performed */
static const unsigned int DENDROSOLVER_RK3_STAGES = 3;

/** @brief: Adaptive time step update safety factor */
static const double DENDROSOLVER_SAFETY_FAC = 0.8;

/** @brief: Number of internal variables */
static const unsigned int DENDROSOLVER_NUM_VARS_INTENL =
    (DENDROSOLVER_RK45_STAGES + 1) * DENDROSOLVER_NUM_VARS;

/** @brief: Minimum black hole domain, to be added to the parameter file for
 * running! */
extern double DENDROSOLVER_COMPD_MIN[3];

/** @brief: Maximum black hole domain, to be added to the parameter file for
 * running! */
extern double DENDROSOLVER_COMPD_MAX[3];

/** @brief: Minimum coordinates of the OCTREE */
extern double DENDROSOLVER_OCTREE_MIN[3];

/** @brief: Maximum coordinates of the OCTREE */
extern double DENDROSOLVER_OCTREE_MAX[3];

/** @brief: Output frequency for the solution, for saving to VTU file */
extern unsigned int DENDROSOLVER_IO_OUTPUT_FREQ;

/** @brief: Gravitational wave extraction frequency */
extern unsigned int DENDROSOLVER_GW_EXTRACT_FREQ;

/** @brief: Gravitational wave extraction frequency after the merger */
extern unsigned int DENDROSOLVER_GW_EXTRACT_FREQ_AFTER_MERGER;

/** @brief: Timestep output frequency */
extern unsigned int DENDROSOLVER_TIME_STEP_OUTPUT_FREQ;

/** @brief: Frequency for performing remeshing test based on wavelets */
extern unsigned int DENDROSOLVER_REMESH_TEST_FREQ;

/** @brief: Frequency for performing remeshing test based on wavelets after the
 * merger */
extern unsigned int DENDROSOLVER_REMESH_TEST_FREQ_AFTER_MERGER;

/** @brief: Frequency for checkpoint saving */
extern unsigned int DENDROSOLVER_CHECKPT_FREQ;

/** @brief: Option for restoring from a checkpoint (will restore if set to 1) */
extern unsigned int DENDROSOLVER_RESTORE_SOLVER;

/** @brief: Disable AMR and enable block adaptivity */
extern unsigned int DENDROSOLVER_ENABLE_BLOCK_ADAPTIVITY;

/** @brief: File prefix for the VTU files that will be saved */
extern std::string DENDROSOLVER_VTU_FILE_PREFIX;

/** @brief: File prefix for the checkpoint files */
extern std::string DENDROSOLVER_CHKPT_FILE_PREFIX;

/** @brief: File prefix for the intermediate profile files */
extern std::string DENDROSOLVER_PROFILE_FILE_PREFIX;

/** @brief: Number variables for refinement */
extern unsigned int DENDROSOLVER_NUM_REFINE_VARS;

/** @brief: The IDs for the refinement variables, this will depend on the enum
 * that's generated from the Python */
extern unsigned int DENDROSOLVER_REFINE_VARIABLE_INDICES[36];

/** @brief: The number of evolution variables to put in the output of the files
 */
extern unsigned int DENDROSOLVER_NUM_EVOL_VARS_VTU_OUTPUT;

/** @brief: The number of constraint variables written to VTU files */
extern unsigned int DENDROSOLVER_NUM_CONST_VARS_VTU_OUTPUT;

/** @brief: Evolution variable IDs to be written to the VTU files */
extern unsigned int DENDROSOLVER_VTU_OUTPUT_EVOL_INDICES[36];

/** @brief: Constraint variable IDs to be written to the VTU files */
extern unsigned int DENDROSOLVER_VTU_OUTPUT_CONST_INDICES[6];

/** @brief: Solution output gap (instead of frequency, we can use to output the
 * solution if currentTime > lastIOOutputTime + DENDROSOLVER_IO_OUTPUT_GAP) */
extern unsigned int DENDROSOLVER_IO_OUTPUT_GAP;

/** @brief:  Grain size N/p, Where N number of total octants, p number of active
 * cores */
extern unsigned int DENDROSOLVER_DENDRO_GRAIN_SZ;

/** @brief: Dendro coarsening factor, if computed wavelet tol <
 * DENDROSOLVER_DENDRO_AMR_FAC*DENDROSOLVER_WAVELET_TOL */
extern double DENDROSOLVER_DENDRO_AMR_FAC;

/** @brief: Number of grid iterations untill the grid converges */
extern unsigned int DENDROSOLVER_INIT_GRID_ITER;

/** @brief: Splitter fix value */
extern unsigned int DENDROSOLVER_SPLIT_FIX;

/** @brief: The Courant factor: CFL stability number (specifies how
 * dt=DENDROSOLVER_CFL_FACTOR*dx) */
extern double DENDROSOLVER_CFL_FACTOR;

/** @brief: Simulation time begin */
extern unsigned int DENDROSOLVER_RK_TIME_BEGIN;

/** @brief: Simulation time end */
extern unsigned int DENDROSOLVER_RK_TIME_END;

/** @brief: RK method to use (0 -> RK3 , 1 -> RK4, 2 -> RK45) */
extern unsigned int DENDROSOLVER_RK_TYPE;

/** @brief: Prefered time step size (this is overwritten with the specified CFL
 * factor, not recommended to use this) */
extern double DENDROSOLVER_RK45_TIME_STEP_SIZE;

/** @brief: Desired tolerance value for the RK45 method (with adaptive time
 * stepping), NOT CURRENTLY USED */
extern double DENDROSOLVER_RK45_DESIRED_TOL;

/** @brief: The dissipation type to be used */
extern unsigned int DISSIPATION_TYPE;

/** @brief: The dissipation "NC", note this is only called in a comment for
 * "artificial dissipation" which appears to not be defined anywhere */
extern unsigned int DENDROSOLVER_DISSIPATION_NC;

/** @brief: The dissipation "S", note this is only called in a comment for
 * "artificial dissipation" which appears to not be defined anywhere */
extern unsigned int DENDROSOLVER_DISSIPATION_S;

/** @brief: The TS offset for LTS in EMDA */
extern unsigned int DENDROSOLVER_LTS_TS_OFFSET;

/** @brief: Global parameter to track if a merged checkpoint file is written. */
extern bool DENDROSOLVER_MERGED_CHKPT_WRITTEN;

/** @brief: Tolerance for refinement based on EH */
extern double DENDROSOLVER_EH_REFINE_VAL;

/** @brief: Tolerance for coarsening based on EH */
extern double DENDROSOLVER_EH_COARSEN_VAL;

/** @brief: Whether to output only the z slice in the VTU file */
extern bool DENDROSOLVER_VTU_Z_SLICE_ONLY;

/** @brief: Variable group size for the asynchronous unzip operation. This is an
 * async communication. (Upper bound should be DENDROSOLVER_NUM_VARS) */
extern unsigned int DENDROSOLVER_ASYNC_COMM_K;

/** @brief: Dendro load imbalance tolerance for flexible partitioning */
extern double DENDROSOLVER_LOAD_IMB_TOL;

/** @brief: Dimensionality of the octree, (meshing is supported only for 3D) */
extern unsigned int DENDROSOLVER_DIM;

/** @brief: Maximum and minimum levels of refinement of the mesh */
extern unsigned int DENDROSOLVER_MAXDEPTH;

extern unsigned int DENDROSOLVER_MINDEPTH;

/** @brief: Wavelet tolerance */
extern double DENDROSOLVER_WAVELET_TOL;

/** @brief: Wavelet tolerance for GW extraction (after merger) */
extern double DENDROSOLVER_GW_REFINE_WTOL;

/** @brief: Set wavelet tolerance using a function (default 0) */
extern unsigned int DENDROSOLVER_USE_WAVELET_TOL_FUNCTION;

/** @brief: The maximum value of the wavelet tolerance */
extern double DENDROSOLVER_WAVELET_TOL_MAX;

/** @brief: Radius R0 for the wavelet tolerance function */
extern double DENDROSOLVER_WAVELET_TOL_FUNCTION_R0;

/** @brief: Radius R1 for the wavelet tolerance function */
extern double DENDROSOLVER_WAVELET_TOL_FUNCTION_R1;

/** @brief: Fd intergrid transfer enable or disable */
extern bool DENDROSOLVER_USE_FD_GRID_TRANSFER;

/** @brief: Refinement mode: 0 -> WAMR , 1 -> EH, 2 -> EH_WAMR 3 -> BH_loc based
 */
extern RefinementMode DENDROSOLVER_REFINEMENT_MODE;

extern double DENDROSOLVER_BLK_MIN_X;

extern double DENDROSOLVER_BLK_MIN_Y;

extern double DENDROSOLVER_BLK_MIN_Z;

extern double DENDROSOLVER_BLK_MAX_X;

extern double DENDROSOLVER_BLK_MAX_Y;

extern double DENDROSOLVER_BLK_MAX_Z;

extern double ETA_CONST;

extern double ETA_R0;

extern double ETA_DAMPING;

extern double ETA_DAMPING_EXP;

extern double ANG_PAR;

extern double CHI_FLOOR;

extern double DENDROSOLVER_TRK0;

extern double KO_DISS_SIGMA;

extern double DENDROSOLVER_ETA_R0;

extern unsigned int DENDROSOLVER_ID_TYPE;

extern double DENDROSOLVER_GRID_MIN_X;

extern double DENDROSOLVER_GRID_MAX_X;

extern double DENDROSOLVER_GRID_MIN_Y;

extern double DENDROSOLVER_GRID_MAX_Y;

extern double DENDROSOLVER_GRID_MIN_Z;

extern double DENDROSOLVER_GRID_MAX_Z;

/** @brief: AMR radius for the BH location based refinement for the first black
 * hole. */
extern double DENDROSOLVER_BH1_AMR_R;

/** @brief: Skip grid point distance from bh < DENDROSOLVER_BH1_CONSTRAINT_R
 * when computing the constraint norms for BH1 */
extern unsigned int DENDROSOLVER_BH1_CONSTRAINT_R;

/** @brief: Maximum refinement level for BH1 */
extern unsigned int DENDROSOLVER_BH1_MAX_LEV;

/** @brief: AMR radius for the BH location based refinement for the second black
 * hole. */
extern double DENDROSOLVER_BH2_AMR_R;

/** @brief: Skip grid point distance from bh < DENDROSOLVER_BH2_CONSTRAINT_R
 * when computing the constraint norms  for BH2 */
extern unsigned int DENDROSOLVER_BH2_CONSTRAINT_R;

/** @brief: Maximum refinement level for BH2 */
extern unsigned int DENDROSOLVER_BH2_MAX_LEV;

/** @brief: parameters for the eta_damping function */
extern double DENDROSOLVER_ETA_POWER[2];

extern unsigned int DENDROSOLVER_XI_0;

extern unsigned int DENDROSOLVER_XI_1;

extern unsigned int DENDROSOLVER_XI_2;

extern double DENDROSOLVER_BH2_MASS;

extern double DENDROSOLVER_BH2_X;

extern double DENDROSOLVER_BH2_Y;

extern double DENDROSOLVER_BH2_Z;

extern double DENDROSOLVER_BH2_V_X;

extern double DENDROSOLVER_BH2_V_Y;

extern double DENDROSOLVER_BH2_V_Z;

extern double DENDROSOLVER_BH2_SPIN;

extern double DENDROSOLVER_BH2_SPIN_THETA;

extern double DENDROSOLVER_BH2_SPIN_PHI;

extern double DENDROSOLVER_BH1_MASS;

extern double DENDROSOLVER_BH1_X;

extern double DENDROSOLVER_BH1_Y;

extern double DENDROSOLVER_BH1_Z;

extern double DENDROSOLVER_BH1_V_X;

extern double DENDROSOLVER_BH1_V_Y;

extern double DENDROSOLVER_BH1_V_Z;

extern double DENDROSOLVER_BH1_SPIN;

extern double DENDROSOLVER_BH1_SPIN_THETA;

extern double DENDROSOLVER_BH1_SPIN_PHI;

}  // namespace dsolve

//[[[end]]]

#endif  // DENDROSOLVER_PARAMETERS_H_
