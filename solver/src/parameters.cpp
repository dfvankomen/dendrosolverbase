
#include "parameters.h"

/**
 * Global parameters used across the program
 *
 * NOTE: this will be generated via Python scripts in the future
 *
 */

// clang-format off
/*[[[cog

import cog
import dendrosym
cog.outl("// clang-format on")

paramh_str, paramc_str = dendrosym.params.generate_all_parameter_text("dsolve",
PARAM_SETUP_FILE)

cog.outl(paramc_str)

]]]*/
// clang-format on
namespace dsolve {
mem::memory_pool<double> DENDROSOLVER_MEM_POOL =
    mem::memory_pool<double>(0, 16);
BH BH1;
BH BH2;
// RefinementMode DENDROSOLVER_REFINEMENT_MODE = RefinementMode::WAMR;
Point DENDROSOLVER_BH_LOC[2];

// NECESSARY ALLOCATION/START FOR DERIV WORKSPACE
double* EMDA_DERIV_WORKSPACE = nullptr;
}  // namespace dsolve

namespace GW {
unsigned int DENDROSOLVER_GW_NUM_RADAII = 6;
unsigned int DENDROSOLVER_GW_NUM_LMODES = 3;
unsigned int DENDROSOLVER_GW_RADAII[6] = {50, 60, 70, 80, 90, 100};
unsigned int DENDROSOLVER_GW_L_MODES[3] = {2, 3, 4};
}  // namespace GW

namespace BHLOC {
unsigned int EXTRACTION_VAR_ID = 0;
double EXTRACTION_TOL = 0.3;
}  // namespace BHLOC

namespace TPID {
std::string FILE_PREFIX = "rit_q2";
double par_b = 4.0;
double par_m_plus = 0.3171512710491704;
double par_m_minus = 0.6515011227547999;
double target_M_plus = 0.3171512710491704;
double target_M_minus = 0.6515011227547999;
double par_P_plus[3] = {-0.0017777739959719536, 0.10049503227336423, 0.0};
double par_P_minus[3] = {0.0017777739959719536, -0.10049503227336423, 0.0};
double par_S_plus[3] = {0.0, 0.0, 0.0};
double par_S_minus[3] = {0.0, 0.0, 0.0};
double center_offset[3] = {1.323826943609176, 0.0, 0.0};
unsigned int give_bare_mass = 1;
int initial_lapse = 2;
unsigned int grid_setup_method = 1;
int solve_momentum_constraint = 1;
unsigned int verbose = 1;
double adm_tol = 1e-10;
double Newton_tol = 1e-10;
double initial_lapse_psi_exponent = -2.0;
unsigned int npoints_A = 60;
unsigned int npoints_B = 60;
unsigned int npoints_phi = 60;
}  // namespace TPID

namespace dsolve {
unsigned int DENDROSOLVER_LAMBDA[4] = {1, 1, 1, 1};
double DENDROSOLVER_ETA[2] = {0.4, 0.4};
double DENDROSOLVER_ETADAMP = 2.0;
double DENDROSOLVER_ALPHA_THEORY[2] = {1.0, 1.0};
double DENDROSOLVER_LF[2] = {1.0, 0.0};
double DENDROSOLVER_P_EXPO = -1.0;
unsigned int DENDROSOLVER_ELE_ORDER = 6;
unsigned int DENDROSOLVER_PADDING_WIDTH = DENDROSOLVER_ELE_ORDER >> 1u;
double DENDROSOLVER_COMPD_MIN[3] = {-50.0, -50.0, -50.0};
double DENDROSOLVER_COMPD_MAX[3] = {50.0, 50.0, 50.0};
double DENDROSOLVER_OCTREE_MIN[3] = {0.0, 0.0, 0.0};
double DENDROSOLVER_OCTREE_MAX[3] = {(double)(1u << DENDROSOLVER_MAXDEPTH),
                                     (double)(1u << DENDROSOLVER_MAXDEPTH),
                                     (double)(1u << DENDROSOLVER_MAXDEPTH)};
unsigned int DENDROSOLVER_IO_OUTPUT_FREQ = 1000;
unsigned int DENDROSOLVER_GW_EXTRACT_FREQ = 10;
unsigned int DENDROSOLVER_GW_EXTRACT_FREQ_AFTER_MERGER = 10;
unsigned int DENDROSOLVER_TIME_STEP_OUTPUT_FREQ = 10;
unsigned int DENDROSOLVER_REMESH_TEST_FREQ = 10;
unsigned int DENDROSOLVER_REMESH_TEST_FREQ_AFTER_MERGER = 10;
unsigned int DENDROSOLVER_CHECKPT_FREQ = 5000;
unsigned int DENDROSOLVER_RESTORE_SOLVER = 0;
unsigned int DENDROSOLVER_ENABLE_BLOCK_ADAPTIVITY = 0;
std::string DENDROSOLVER_VTU_FILE_PREFIX = "vtu/solver_gr";
std::string DENDROSOLVER_CHKPT_FILE_PREFIX = "cp/solver_cp";
std::string DENDROSOLVER_PROFILE_FILE_PREFIX = "solver_prof";
unsigned int DENDROSOLVER_NUM_REFINE_VARS = 36;
unsigned int DENDROSOLVER_REFINE_VARIABLE_INDICES[36] = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17,
    18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
unsigned int DENDROSOLVER_NUM_EVOL_VARS_VTU_OUTPUT = 14;
unsigned int DENDROSOLVER_NUM_CONST_VARS_VTU_OUTPUT = 1;
unsigned int DENDROSOLVER_VTU_OUTPUT_EVOL_INDICES[36] = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17,
    18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
unsigned int DENDROSOLVER_VTU_OUTPUT_CONST_INDICES[6] = {0, 1, 2, 3, 4, 5};
unsigned int DENDROSOLVER_IO_OUTPUT_GAP = 1;
unsigned int DENDROSOLVER_DENDRO_GRAIN_SZ = 50;
double DENDROSOLVER_DENDRO_AMR_FAC = 0.1;
unsigned int DENDROSOLVER_INIT_GRID_ITER = 10;
unsigned int DENDROSOLVER_SPLIT_FIX = 2;
double DENDROSOLVER_CFL_FACTOR = 0.25;
unsigned int DENDROSOLVER_RK_TIME_BEGIN = 0;
unsigned int DENDROSOLVER_RK_TIME_END = 800;
unsigned int DENDROSOLVER_RK_TYPE = 1;
double DENDROSOLVER_RK45_TIME_STEP_SIZE = 0.01;
double DENDROSOLVER_RK45_DESIRED_TOL = 0.001;
unsigned int DISSIPATION_TYPE = 0;
unsigned int DENDROSOLVER_DISSIPATION_NC = 0;
unsigned int DENDROSOLVER_DISSIPATION_S = 0;
unsigned int DENDROSOLVER_LTS_TS_OFFSET = 0;
bool DENDROSOLVER_MERGED_CHKPT_WRITTEN = false;
double DENDROSOLVER_EH_REFINE_VAL = 0.4;
double DENDROSOLVER_EH_COARSEN_VAL = 0.6;
bool DENDROSOLVER_VTU_Z_SLICE_ONLY = true;
unsigned int DENDROSOLVER_ASYNC_COMM_K = 4;
double DENDROSOLVER_LOAD_IMB_TOL = 0.1;
unsigned int DENDROSOLVER_DIM = 3;
unsigned int DENDROSOLVER_MAXDEPTH = 16;
unsigned int DENDROSOLVER_MINDEPTH = 3;
double DENDROSOLVER_WAVELET_TOL = 1e-05;
double DENDROSOLVER_GW_REFINE_WTOL = 0.0001;
unsigned int DENDROSOLVER_USE_WAVELET_TOL_FUNCTION = 3;
double DENDROSOLVER_WAVELET_TOL_MAX = 0.001;
double DENDROSOLVER_WAVELET_TOL_FUNCTION_R0 = 30.0;
double DENDROSOLVER_WAVELET_TOL_FUNCTION_R1 = 220.0;
bool DENDROSOLVER_USE_FD_GRID_TRANSFER = false;
RefinementMode DENDROSOLVER_REFINEMENT_MODE = static_cast<RefinementMode>(0);
double DENDROSOLVER_BLK_MIN_X = -6.0;
double DENDROSOLVER_BLK_MIN_Y = -6.0;
double DENDROSOLVER_BLK_MIN_Z = -1.0;
double DENDROSOLVER_BLK_MAX_X = 6.0;
double DENDROSOLVER_BLK_MAX_Y = 6.0;
double DENDROSOLVER_BLK_MAX_Z = 1.0;
double ETA_CONST = 2.0;
double ETA_R0 = 30.0;
double ETA_DAMPING = 2.0;
double ETA_DAMPING_EXP = 2.0;
double ANG_PAR = 0.01;
double CHI_FLOOR = 0.0001;
double DENDROSOLVER_TRK0 = 0.0;
double KO_DISS_SIGMA = 0.4;
double DENDROSOLVER_ETA_R0 = 1.31;
unsigned int DENDROSOLVER_ID_TYPE = 0;
double DENDROSOLVER_GRID_MIN_X = -400.0;
double DENDROSOLVER_GRID_MAX_X = 400.0;
double DENDROSOLVER_GRID_MIN_Y = -400.0;
double DENDROSOLVER_GRID_MAX_Y = 400.0;
double DENDROSOLVER_GRID_MIN_Z = -400.0;
double DENDROSOLVER_GRID_MAX_Z = 400.0;
double DENDROSOLVER_BH1_AMR_R = 0.07;
unsigned int DENDROSOLVER_BH1_CONSTRAINT_R = 5;
unsigned int DENDROSOLVER_BH1_MAX_LEV;
double DENDROSOLVER_BH2_AMR_R = 1.3;
unsigned int DENDROSOLVER_BH2_CONSTRAINT_R = 5;
unsigned int DENDROSOLVER_BH2_MAX_LEV;
double DENDROSOLVER_ETA_POWER[2] = {2.0, 2.0};
unsigned int DENDROSOLVER_XI_0 = 0;
unsigned int DENDROSOLVER_XI_1 = 0;
unsigned int DENDROSOLVER_XI_2 = 0;
double DENDROSOLVER_BH2_MASS = 0.6515011227547999;
double DENDROSOLVER_BH2_X = -2.676173056390826;
double DENDROSOLVER_BH2_Y = 0.0;
double DENDROSOLVER_BH2_Z = 0.0;
double DENDROSOLVER_BH2_V_X = -0.0017777739959719536;
double DENDROSOLVER_BH2_V_Y = -0.10049503227336423;
double DENDROSOLVER_BH2_V_Z = 0.0;
double DENDROSOLVER_BH2_SPIN = 0.0;
double DENDROSOLVER_BH2_SPIN_THETA = 0.0;
double DENDROSOLVER_BH2_SPIN_PHI = 0.0;
double DENDROSOLVER_BH1_MASS = 0.3171512710491704;
double DENDROSOLVER_BH1_X = 5.323826943609176;
double DENDROSOLVER_BH1_Y = 0.0;
double DENDROSOLVER_BH1_Z = 0.0;
double DENDROSOLVER_BH1_V_X = -0.0017777739959719536;
double DENDROSOLVER_BH1_V_Y = 0.10049503227336423;
double DENDROSOLVER_BH1_V_Z = 0.0;
double DENDROSOLVER_BH1_SPIN = 0.0;
double DENDROSOLVER_BH1_SPIN_THETA = 0.0;
double DENDROSOLVER_BH1_SPIN_PHI = 0.0;
}  // namespace dsolve
namespace dsolve {
void readParamFile(const char* inFile, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    auto file = toml::parse(inFile);

    if (!rank) {
        if (file.contains("GW::DENDROSOLVER_GW_NUM_RADAII")) {
            GW::DENDROSOLVER_GW_NUM_RADAII =
                file["GW::DENDROSOLVER_GW_NUM_RADAII"].as_integer();
        }

        if (file.contains("GW::DENDROSOLVER_GW_NUM_LMODES")) {
            GW::DENDROSOLVER_GW_NUM_LMODES =
                file["GW::DENDROSOLVER_GW_NUM_LMODES"].as_integer();
        }

        if (file.contains("GW::DENDROSOLVER_GW_RADAII")) {
            for (int i = 0; i < 6; ++i) {
                GW::DENDROSOLVER_GW_RADAII[i] =
                    file["GW::DENDROSOLVER_GW_RADAII"][i].as_integer();
            }
        }

        if (file.contains("GW::DENDROSOLVER_GW_L_MODES")) {
            for (int i = 0; i < 3; ++i) {
                GW::DENDROSOLVER_GW_L_MODES[i] =
                    file["GW::DENDROSOLVER_GW_L_MODES"][i].as_integer();
            }
        }

        if (file.contains("BHLOC::EXTRACTION_VAR_ID")) {
            BHLOC::EXTRACTION_VAR_ID =
                file["BHLOC::EXTRACTION_VAR_ID"].as_integer();
        }

        if (file.contains("BHLOC::EXTRACTION_TOL")) {
            if (-1.0 > file["BHLOC::EXTRACTION_TOL"].as_floating() ||
                1.0 < file["BHLOC::EXTRACTION_TOL"].as_floating()) {
                std::cerr << R"(Invalid value for "BHLOC::EXTRACTION_TOL")"
                          << std::endl;
                exit(-1);
            }

            BHLOC::EXTRACTION_TOL = file["BHLOC::EXTRACTION_TOL"].as_floating();
        }

        if (file.contains("TPID::FILE_PREFIX")) {
            TPID::FILE_PREFIX = file["TPID::FILE_PREFIX"].as_string();
        }

        if (file.contains("TPID::par_b")) {
            if (0 > file["TPID::par_b"].as_floating() ||
                8.0 < file["TPID::par_b"].as_floating()) {
                std::cerr << R"(Invalid value for "TPID::par_b")" << std::endl;
                exit(-1);
            }

            TPID::par_b = file["TPID::par_b"].as_floating();
        }

        if (file.contains("TPID::par_m_plus")) {
            if (0 > file["TPID::par_m_plus"].as_floating() ||
                2.0 < file["TPID::par_m_plus"].as_floating()) {
                std::cerr << R"(Invalid value for "TPID::par_m_plus")"
                          << std::endl;
                exit(-1);
            }

            TPID::par_m_plus = file["TPID::par_m_plus"].as_floating();
        }

        if (file.contains("TPID::par_m_minus")) {
            if (0 > file["TPID::par_m_minus"].as_floating() ||
                2.0 < file["TPID::par_m_minus"].as_floating()) {
                std::cerr << R"(Invalid value for "TPID::par_m_minus")"
                          << std::endl;
                exit(-1);
            }

            TPID::par_m_minus = file["TPID::par_m_minus"].as_floating();
        }

        if (file.contains("TPID::target_M_plus")) {
            if (0 > file["TPID::target_M_plus"].as_floating() ||
                2.0 < file["TPID::target_M_plus"].as_floating()) {
                std::cerr << R"(Invalid value for "TPID::target_M_plus")"
                          << std::endl;
                exit(-1);
            }

            TPID::target_M_plus = file["TPID::target_M_plus"].as_floating();
        }

        if (file.contains("TPID::target_M_minus")) {
            if (0 > file["TPID::target_M_minus"].as_floating() ||
                2.0 < file["TPID::target_M_minus"].as_floating()) {
                std::cerr << R"(Invalid value for "TPID::target_M_minus")"
                          << std::endl;
                exit(-1);
            }

            TPID::target_M_minus = file["TPID::target_M_minus"].as_floating();
        }

        if (file.contains("TPID::par_P_plus")) {
            for (int i = 0; i < 3; ++i) {
                TPID::par_P_plus[i] = file["TPID::par_P_plus"][i].as_floating();
            }
        }

        if (file.contains("TPID::par_P_minus")) {
            for (int i = 0; i < 3; ++i) {
                TPID::par_P_minus[i] =
                    file["TPID::par_P_minus"][i].as_floating();
            }
        }

        if (file.contains("TPID::par_S_plus")) {
            for (int i = 0; i < 3; ++i) {
                TPID::par_S_plus[i] = file["TPID::par_S_plus"][i].as_floating();
            }
        }

        if (file.contains("TPID::par_S_minus")) {
            for (int i = 0; i < 3; ++i) {
                TPID::par_S_minus[i] =
                    file["TPID::par_S_minus"][i].as_floating();
            }
        }

        if (file.contains("TPID::center_offset")) {
            for (int i = 0; i < 3; ++i) {
                TPID::center_offset[i] =
                    file["TPID::center_offset"][i].as_floating();
            }
        }

        if (file.contains("TPID::give_bare_mass")) {
            TPID::give_bare_mass = file["TPID::give_bare_mass"].as_integer();
        }

        if (file.contains("TPID::initial_lapse")) {
            if (0 > file["TPID::initial_lapse"].as_integer() ||
                3 < file["TPID::initial_lapse"].as_integer()) {
                std::cerr << R"(Invalid value for "TPID::initial_lapse")"
                          << std::endl;
                exit(-1);
            }

            TPID::initial_lapse = file["TPID::initial_lapse"].as_integer();
        }

        if (file.contains("TPID::grid_setup_method")) {
            TPID::grid_setup_method =
                file["TPID::grid_setup_method"].as_integer();
        }

        if (file.contains("TPID::solve_momentum_constraint")) {
            if (0 > file["TPID::solve_momentum_constraint"].as_integer() ||
                2 < file["TPID::solve_momentum_constraint"].as_integer()) {
                std::cerr
                    << R"(Invalid value for "TPID::solve_momentum_constraint")"
                    << std::endl;
                exit(-1);
            }

            TPID::solve_momentum_constraint =
                file["TPID::solve_momentum_constraint"].as_integer();
        }

        if (file.contains("TPID::verbose")) {
            TPID::verbose = file["TPID::verbose"].as_integer();
        }

        if (file.contains("TPID::adm_tol")) {
            if (0 > file["TPID::adm_tol"].as_floating() ||
                2e-10 < file["TPID::adm_tol"].as_floating()) {
                std::cerr << R"(Invalid value for "TPID::adm_tol")"
                          << std::endl;
                exit(-1);
            }

            TPID::adm_tol = file["TPID::adm_tol"].as_floating();
        }

        if (file.contains("TPID::Newton_tol")) {
            if (0 > file["TPID::Newton_tol"].as_floating() ||
                2e-10 < file["TPID::Newton_tol"].as_floating()) {
                std::cerr << R"(Invalid value for "TPID::Newton_tol")"
                          << std::endl;
                exit(-1);
            }

            TPID::Newton_tol = file["TPID::Newton_tol"].as_floating();
        }

        if (file.contains("TPID::initial_lapse_psi_exponent")) {
            if (-4.0 > file["TPID::initial_lapse_psi_exponent"].as_floating() ||
                4.0 < file["TPID::initial_lapse_psi_exponent"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "TPID::initial_lapse_psi_exponent")"
                    << std::endl;
                exit(-1);
            }

            TPID::initial_lapse_psi_exponent =
                file["TPID::initial_lapse_psi_exponent"].as_floating();
        }

        if (file.contains("TPID::npoints_A")) {
            TPID::npoints_A = file["TPID::npoints_A"].as_integer();
        }

        if (file.contains("TPID::npoints_B")) {
            TPID::npoints_B = file["TPID::npoints_B"].as_integer();
        }

        if (file.contains("TPID::npoints_phi")) {
            TPID::npoints_phi = file["TPID::npoints_phi"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_LAMBDA")) {
            for (int i = 0; i < 4; ++i) {
                dsolve::DENDROSOLVER_LAMBDA[i] =
                    file["dsolve::DENDROSOLVER_LAMBDA"][i].as_integer();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_ETA")) {
            for (int i = 0; i < 2; ++i) {
                dsolve::DENDROSOLVER_ETA[i] =
                    file["dsolve::DENDROSOLVER_ETA"][i].as_floating();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_ETADAMP")) {
            if (0.5 > file["dsolve::DENDROSOLVER_ETADAMP"].as_floating() ||
                3.2 < file["dsolve::DENDROSOLVER_ETADAMP"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_ETADAMP")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_ETADAMP =
                file["dsolve::DENDROSOLVER_ETADAMP"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_ALPHA_THEORY")) {
            for (int i = 0; i < 2; ++i) {
                dsolve::DENDROSOLVER_ALPHA_THEORY[i] =
                    file["dsolve::DENDROSOLVER_ALPHA_THEORY"][i].as_floating();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_LF")) {
            for (int i = 0; i < 2; ++i) {
                dsolve::DENDROSOLVER_LF[i] =
                    file["dsolve::DENDROSOLVER_LF"][i].as_floating();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_P_EXPO")) {
            if (-1.0 > file["dsolve::DENDROSOLVER_P_EXPO"].as_floating() ||
                4.0 < file["dsolve::DENDROSOLVER_P_EXPO"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_P_EXPO")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_P_EXPO =
                file["dsolve::DENDROSOLVER_P_EXPO"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_ELE_ORDER")) {
            dsolve::DENDROSOLVER_ELE_ORDER =
                file["dsolve::DENDROSOLVER_ELE_ORDER"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_PADDING_WIDTH")) {
            dsolve::DENDROSOLVER_PADDING_WIDTH =
                file["dsolve::DENDROSOLVER_PADDING_WIDTH"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_COMPD_MIN")) {
            for (int i = 0; i < 3; ++i) {
                dsolve::DENDROSOLVER_COMPD_MIN[i] =
                    file["dsolve::DENDROSOLVER_COMPD_MIN"][i].as_floating();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_COMPD_MAX")) {
            for (int i = 0; i < 3; ++i) {
                dsolve::DENDROSOLVER_COMPD_MAX[i] =
                    file["dsolve::DENDROSOLVER_COMPD_MAX"][i].as_floating();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_OCTREE_MIN")) {
            for (int i = 0; i < 3; ++i) {
                dsolve::DENDROSOLVER_OCTREE_MIN[i] =
                    file["dsolve::DENDROSOLVER_OCTREE_MIN"][i].as_floating();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_OCTREE_MAX")) {
            for (int i = 0; i < 3; ++i) {
                dsolve::DENDROSOLVER_OCTREE_MAX[i] =
                    file["dsolve::DENDROSOLVER_OCTREE_MAX"][i].as_floating();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_IO_OUTPUT_FREQ")) {
            dsolve::DENDROSOLVER_IO_OUTPUT_FREQ =
                file["dsolve::DENDROSOLVER_IO_OUTPUT_FREQ"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_GW_EXTRACT_FREQ")) {
            dsolve::DENDROSOLVER_GW_EXTRACT_FREQ =
                file["dsolve::DENDROSOLVER_GW_EXTRACT_FREQ"].as_integer();
        }

        if (file.contains(
                "dsolve::DENDROSOLVER_GW_EXTRACT_FREQ_AFTER_MERGER")) {
            dsolve::DENDROSOLVER_GW_EXTRACT_FREQ_AFTER_MERGER =
                file["dsolve::DENDROSOLVER_GW_EXTRACT_FREQ_AFTER_MERGER"]
                    .as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_TIME_STEP_OUTPUT_FREQ")) {
            dsolve::DENDROSOLVER_TIME_STEP_OUTPUT_FREQ =
                file["dsolve::DENDROSOLVER_TIME_STEP_OUTPUT_FREQ"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_REMESH_TEST_FREQ")) {
            dsolve::DENDROSOLVER_REMESH_TEST_FREQ =
                file["dsolve::DENDROSOLVER_REMESH_TEST_FREQ"].as_integer();
        }

        if (file.contains(
                "dsolve::DENDROSOLVER_REMESH_TEST_FREQ_AFTER_MERGER")) {
            dsolve::DENDROSOLVER_REMESH_TEST_FREQ_AFTER_MERGER =
                file["dsolve::DENDROSOLVER_REMESH_TEST_FREQ_AFTER_MERGER"]
                    .as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_CHECKPT_FREQ")) {
            dsolve::DENDROSOLVER_CHECKPT_FREQ =
                file["dsolve::DENDROSOLVER_CHECKPT_FREQ"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_RESTORE_SOLVER")) {
            dsolve::DENDROSOLVER_RESTORE_SOLVER =
                file["dsolve::DENDROSOLVER_RESTORE_SOLVER"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_ENABLE_BLOCK_ADAPTIVITY")) {
            dsolve::DENDROSOLVER_ENABLE_BLOCK_ADAPTIVITY =
                file["dsolve::DENDROSOLVER_ENABLE_BLOCK_ADAPTIVITY"]
                    .as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_VTU_FILE_PREFIX")) {
            dsolve::DENDROSOLVER_VTU_FILE_PREFIX =
                file["dsolve::DENDROSOLVER_VTU_FILE_PREFIX"].as_string();
        }

        if (file.contains("dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX")) {
            dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX =
                file["dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX"].as_string();
        }

        if (file.contains("dsolve::DENDROSOLVER_PROFILE_FILE_PREFIX")) {
            dsolve::DENDROSOLVER_PROFILE_FILE_PREFIX =
                file["dsolve::DENDROSOLVER_PROFILE_FILE_PREFIX"].as_string();
        }

        if (file.contains("dsolve::DENDROSOLVER_NUM_REFINE_VARS")) {
            dsolve::DENDROSOLVER_NUM_REFINE_VARS =
                file["dsolve::DENDROSOLVER_NUM_REFINE_VARS"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_REFINE_VARIABLE_INDICES")) {
            for (int i = 0; i < 36; ++i) {
                dsolve::DENDROSOLVER_REFINE_VARIABLE_INDICES[i] =
                    file["dsolve::DENDROSOLVER_REFINE_VARIABLE_INDICES"][i]
                        .as_integer();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_NUM_EVOL_VARS_VTU_OUTPUT")) {
            dsolve::DENDROSOLVER_NUM_EVOL_VARS_VTU_OUTPUT =
                file["dsolve::DENDROSOLVER_NUM_EVOL_VARS_VTU_OUTPUT"]
                    .as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_NUM_CONST_VARS_VTU_OUTPUT")) {
            dsolve::DENDROSOLVER_NUM_CONST_VARS_VTU_OUTPUT =
                file["dsolve::DENDROSOLVER_NUM_CONST_VARS_VTU_OUTPUT"]
                    .as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_VTU_OUTPUT_EVOL_INDICES")) {
            for (int i = 0; i < 36; ++i) {
                dsolve::DENDROSOLVER_VTU_OUTPUT_EVOL_INDICES[i] =
                    file["dsolve::DENDROSOLVER_VTU_OUTPUT_EVOL_INDICES"][i]
                        .as_integer();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_VTU_OUTPUT_CONST_INDICES")) {
            for (int i = 0; i < 6; ++i) {
                dsolve::DENDROSOLVER_VTU_OUTPUT_CONST_INDICES[i] =
                    file["dsolve::DENDROSOLVER_VTU_OUTPUT_CONST_INDICES"][i]
                        .as_integer();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_IO_OUTPUT_GAP")) {
            dsolve::DENDROSOLVER_IO_OUTPUT_GAP =
                file["dsolve::DENDROSOLVER_IO_OUTPUT_GAP"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_DENDRO_GRAIN_SZ")) {
            dsolve::DENDROSOLVER_DENDRO_GRAIN_SZ =
                file["dsolve::DENDROSOLVER_DENDRO_GRAIN_SZ"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_DENDRO_AMR_FAC")) {
            if (0.0 >
                    file["dsolve::DENDROSOLVER_DENDRO_AMR_FAC"].as_floating() ||
                0.2 <
                    file["dsolve::DENDROSOLVER_DENDRO_AMR_FAC"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_DENDRO_AMR_FAC")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_DENDRO_AMR_FAC =
                file["dsolve::DENDROSOLVER_DENDRO_AMR_FAC"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_INIT_GRID_ITER")) {
            dsolve::DENDROSOLVER_INIT_GRID_ITER =
                file["dsolve::DENDROSOLVER_INIT_GRID_ITER"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_SPLIT_FIX")) {
            dsolve::DENDROSOLVER_SPLIT_FIX =
                file["dsolve::DENDROSOLVER_SPLIT_FIX"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_CFL_FACTOR")) {
            if (0.0 > file["dsolve::DENDROSOLVER_CFL_FACTOR"].as_floating() ||
                0.5 < file["dsolve::DENDROSOLVER_CFL_FACTOR"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_CFL_FACTOR")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_CFL_FACTOR =
                file["dsolve::DENDROSOLVER_CFL_FACTOR"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_RK_TIME_BEGIN")) {
            dsolve::DENDROSOLVER_RK_TIME_BEGIN =
                file["dsolve::DENDROSOLVER_RK_TIME_BEGIN"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_RK_TIME_END")) {
            dsolve::DENDROSOLVER_RK_TIME_END =
                file["dsolve::DENDROSOLVER_RK_TIME_END"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_RK_TYPE")) {
            dsolve::DENDROSOLVER_RK_TYPE =
                file["dsolve::DENDROSOLVER_RK_TYPE"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE")) {
            if (0.0 > file["dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE"]
                          .as_floating() ||
                0.02 < file["dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE"]
                           .as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE =
                file["dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_RK45_DESIRED_TOL")) {
            if (0.0 > file["dsolve::DENDROSOLVER_RK45_DESIRED_TOL"]
                          .as_floating() ||
                0.002 < file["dsolve::DENDROSOLVER_RK45_DESIRED_TOL"]
                            .as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_RK45_DESIRED_TOL")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_RK45_DESIRED_TOL =
                file["dsolve::DENDROSOLVER_RK45_DESIRED_TOL"].as_floating();
        }

        if (file.contains("dsolve::DISSIPATION_TYPE")) {
            dsolve::DISSIPATION_TYPE =
                file["dsolve::DISSIPATION_TYPE"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_DISSIPATION_NC")) {
            dsolve::DENDROSOLVER_DISSIPATION_NC =
                file["dsolve::DENDROSOLVER_DISSIPATION_NC"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_DISSIPATION_S")) {
            dsolve::DENDROSOLVER_DISSIPATION_S =
                file["dsolve::DENDROSOLVER_DISSIPATION_S"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_LTS_TS_OFFSET")) {
            dsolve::DENDROSOLVER_LTS_TS_OFFSET =
                file["dsolve::DENDROSOLVER_LTS_TS_OFFSET"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_MERGED_CHKPT_WRITTEN")) {
            dsolve::DENDROSOLVER_MERGED_CHKPT_WRITTEN =
                file["dsolve::DENDROSOLVER_MERGED_CHKPT_WRITTEN"].as_boolean();
        }

        if (file.contains("dsolve::DENDROSOLVER_EH_REFINE_VAL")) {
            if (0.0 >
                    file["dsolve::DENDROSOLVER_EH_REFINE_VAL"].as_floating() ||
                1.0 <
                    file["dsolve::DENDROSOLVER_EH_REFINE_VAL"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_EH_REFINE_VAL")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_EH_REFINE_VAL =
                file["dsolve::DENDROSOLVER_EH_REFINE_VAL"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_EH_COARSEN_VAL")) {
            if (0.0 >
                    file["dsolve::DENDROSOLVER_EH_COARSEN_VAL"].as_floating() ||
                1.2 <
                    file["dsolve::DENDROSOLVER_EH_COARSEN_VAL"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_EH_COARSEN_VAL")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_EH_COARSEN_VAL =
                file["dsolve::DENDROSOLVER_EH_COARSEN_VAL"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_VTU_Z_SLICE_ONLY")) {
            dsolve::DENDROSOLVER_VTU_Z_SLICE_ONLY =
                file["dsolve::DENDROSOLVER_VTU_Z_SLICE_ONLY"].as_boolean();
        }

        if (file.contains("dsolve::DENDROSOLVER_ASYNC_COMM_K")) {
            dsolve::DENDROSOLVER_ASYNC_COMM_K =
                file["dsolve::DENDROSOLVER_ASYNC_COMM_K"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_LOAD_IMB_TOL")) {
            if (0.0 > file["dsolve::DENDROSOLVER_LOAD_IMB_TOL"].as_floating() ||
                0.2 < file["dsolve::DENDROSOLVER_LOAD_IMB_TOL"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_LOAD_IMB_TOL")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_LOAD_IMB_TOL =
                file["dsolve::DENDROSOLVER_LOAD_IMB_TOL"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_DIM")) {
            dsolve::DENDROSOLVER_DIM =
                file["dsolve::DENDROSOLVER_DIM"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_MAXDEPTH")) {
            dsolve::DENDROSOLVER_MAXDEPTH =
                file["dsolve::DENDROSOLVER_MAXDEPTH"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_MINDEPTH")) {
            dsolve::DENDROSOLVER_MINDEPTH =
                file["dsolve::DENDROSOLVER_MINDEPTH"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_WAVELET_TOL")) {
            if (0.0 > file["dsolve::DENDROSOLVER_WAVELET_TOL"].as_floating() ||
                2e-05 <
                    file["dsolve::DENDROSOLVER_WAVELET_TOL"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_WAVELET_TOL")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_WAVELET_TOL =
                file["dsolve::DENDROSOLVER_WAVELET_TOL"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_GW_REFINE_WTOL")) {
            if (0.0 >
                    file["dsolve::DENDROSOLVER_GW_REFINE_WTOL"].as_floating() ||
                0.0002 <
                    file["dsolve::DENDROSOLVER_GW_REFINE_WTOL"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_GW_REFINE_WTOL")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_GW_REFINE_WTOL =
                file["dsolve::DENDROSOLVER_GW_REFINE_WTOL"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_USE_WAVELET_TOL_FUNCTION")) {
            dsolve::DENDROSOLVER_USE_WAVELET_TOL_FUNCTION =
                file["dsolve::DENDROSOLVER_USE_WAVELET_TOL_FUNCTION"]
                    .as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_WAVELET_TOL_MAX")) {
            if (0.0 > file["dsolve::DENDROSOLVER_WAVELET_TOL_MAX"]
                          .as_floating() ||
                0.002 < file["dsolve::DENDROSOLVER_WAVELET_TOL_MAX"]
                            .as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_WAVELET_TOL_MAX")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_WAVELET_TOL_MAX =
                file["dsolve::DENDROSOLVER_WAVELET_TOL_MAX"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R0")) {
            if (0.0 > file["dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R0"]
                          .as_floating() ||
                60.0 < file["dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R0"]
                           .as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R0")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R0 =
                file["dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R0"]
                    .as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R1")) {
            if (0.0 > file["dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R1"]
                          .as_floating() ||
                440.0 < file["dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R1"]
                            .as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R1")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R1 =
                file["dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R1"]
                    .as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_USE_FD_GRID_TRANSFER")) {
            dsolve::DENDROSOLVER_USE_FD_GRID_TRANSFER =
                file["dsolve::DENDROSOLVER_USE_FD_GRID_TRANSFER"].as_boolean();
        }

        if (file.contains("dsolve::DENDROSOLVER_REFINEMENT_MODE")) {
            if (0 > file["dsolve::DENDROSOLVER_REFINEMENT_MODE"].as_integer() ||
                3 < file["dsolve::DENDROSOLVER_REFINEMENT_MODE"].as_integer()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_REFINEMENT_MODE")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_REFINEMENT_MODE = static_cast<RefinementMode>(
                file["dsolve::DENDROSOLVER_REFINEMENT_MODE"].as_integer());
        }

        if (file.contains("dsolve::DENDROSOLVER_BLK_MIN_X")) {
            if (-12.0 > file["dsolve::DENDROSOLVER_BLK_MIN_X"].as_floating() ||
                0.0 < file["dsolve::DENDROSOLVER_BLK_MIN_X"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BLK_MIN_X")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BLK_MIN_X =
                file["dsolve::DENDROSOLVER_BLK_MIN_X"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BLK_MIN_Y")) {
            if (-12.0 > file["dsolve::DENDROSOLVER_BLK_MIN_Y"].as_floating() ||
                0.0 < file["dsolve::DENDROSOLVER_BLK_MIN_Y"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BLK_MIN_Y")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BLK_MIN_Y =
                file["dsolve::DENDROSOLVER_BLK_MIN_Y"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BLK_MIN_Z")) {
            if (-12.0 > file["dsolve::DENDROSOLVER_BLK_MIN_Z"].as_floating() ||
                0.0 < file["dsolve::DENDROSOLVER_BLK_MIN_Z"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BLK_MIN_Z")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BLK_MIN_Z =
                file["dsolve::DENDROSOLVER_BLK_MIN_Z"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BLK_MAX_X")) {
            if (0.0 > file["dsolve::DENDROSOLVER_BLK_MAX_X"].as_floating() ||
                12.0 < file["dsolve::DENDROSOLVER_BLK_MAX_X"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BLK_MAX_X")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BLK_MAX_X =
                file["dsolve::DENDROSOLVER_BLK_MAX_X"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BLK_MAX_Y")) {
            if (0.0 > file["dsolve::DENDROSOLVER_BLK_MAX_Y"].as_floating() ||
                12.0 < file["dsolve::DENDROSOLVER_BLK_MAX_Y"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BLK_MAX_Y")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BLK_MAX_Y =
                file["dsolve::DENDROSOLVER_BLK_MAX_Y"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BLK_MAX_Z")) {
            if (0.0 > file["dsolve::DENDROSOLVER_BLK_MAX_Z"].as_floating() ||
                12.0 < file["dsolve::DENDROSOLVER_BLK_MAX_Z"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BLK_MAX_Z")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BLK_MAX_Z =
                file["dsolve::DENDROSOLVER_BLK_MAX_Z"].as_floating();
        }

        if (file.contains("dsolve::ETA_CONST")) {
            if (0.0 > file["dsolve::ETA_CONST"].as_floating() ||
                4.0 < file["dsolve::ETA_CONST"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::ETA_CONST")"
                          << std::endl;
                exit(-1);
            }

            dsolve::ETA_CONST = file["dsolve::ETA_CONST"].as_floating();
        }

        if (file.contains("dsolve::ETA_R0")) {
            if (0.0 > file["dsolve::ETA_R0"].as_floating() ||
                60.0 < file["dsolve::ETA_R0"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::ETA_R0")"
                          << std::endl;
                exit(-1);
            }

            dsolve::ETA_R0 = file["dsolve::ETA_R0"].as_floating();
        }

        if (file.contains("dsolve::ETA_DAMPING")) {
            if (0.0 > file["dsolve::ETA_DAMPING"].as_floating() ||
                4.0 < file["dsolve::ETA_DAMPING"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::ETA_DAMPING")"
                          << std::endl;
                exit(-1);
            }

            dsolve::ETA_DAMPING = file["dsolve::ETA_DAMPING"].as_floating();
        }

        if (file.contains("dsolve::ETA_DAMPING_EXP")) {
            if (0.0 > file["dsolve::ETA_DAMPING_EXP"].as_floating() ||
                4.0 < file["dsolve::ETA_DAMPING_EXP"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::ETA_DAMPING_EXP")"
                          << std::endl;
                exit(-1);
            }

            dsolve::ETA_DAMPING_EXP =
                file["dsolve::ETA_DAMPING_EXP"].as_floating();
        }

        if (file.contains("dsolve::ANG_PAR")) {
            if (0.0 > file["dsolve::ANG_PAR"].as_floating() ||
                0.02 < file["dsolve::ANG_PAR"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::ANG_PAR")"
                          << std::endl;
                exit(-1);
            }

            dsolve::ANG_PAR = file["dsolve::ANG_PAR"].as_floating();
        }

        if (file.contains("dsolve::CHI_FLOOR")) {
            if (0.0 > file["dsolve::CHI_FLOOR"].as_floating() ||
                0.0002 < file["dsolve::CHI_FLOOR"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::CHI_FLOOR")"
                          << std::endl;
                exit(-1);
            }

            dsolve::CHI_FLOOR = file["dsolve::CHI_FLOOR"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_TRK0")) {
            if (0.0 > file["dsolve::DENDROSOLVER_TRK0"].as_floating() ||
                1.0 < file["dsolve::DENDROSOLVER_TRK0"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::DENDROSOLVER_TRK0")"
                          << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_TRK0 =
                file["dsolve::DENDROSOLVER_TRK0"].as_floating();
        }

        if (file.contains("dsolve::KO_DISS_SIGMA")) {
            if (0.0 > file["dsolve::KO_DISS_SIGMA"].as_floating() ||
                0.8 < file["dsolve::KO_DISS_SIGMA"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::KO_DISS_SIGMA")"
                          << std::endl;
                exit(-1);
            }

            dsolve::KO_DISS_SIGMA = file["dsolve::KO_DISS_SIGMA"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_ETA_R0")) {
            if (0.0 > file["dsolve::DENDROSOLVER_ETA_R0"].as_floating() ||
                3.0 < file["dsolve::DENDROSOLVER_ETA_R0"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_ETA_R0")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_ETA_R0 =
                file["dsolve::DENDROSOLVER_ETA_R0"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_ID_TYPE")) {
            dsolve::DENDROSOLVER_ID_TYPE =
                file["dsolve::DENDROSOLVER_ID_TYPE"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_GRID_MIN_X")) {
            if (-500.0 >
                    file["dsolve::DENDROSOLVER_GRID_MIN_X"].as_floating() ||
                0.0 < file["dsolve::DENDROSOLVER_GRID_MIN_X"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_GRID_MIN_X")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_GRID_MIN_X =
                file["dsolve::DENDROSOLVER_GRID_MIN_X"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_GRID_MAX_X")) {
            if (0.0 > file["dsolve::DENDROSOLVER_GRID_MAX_X"].as_floating() ||
                500.0 < file["dsolve::DENDROSOLVER_GRID_MAX_X"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_GRID_MAX_X")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_GRID_MAX_X =
                file["dsolve::DENDROSOLVER_GRID_MAX_X"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_GRID_MIN_Y")) {
            if (-500.0 >
                    file["dsolve::DENDROSOLVER_GRID_MIN_Y"].as_floating() ||
                0.0 < file["dsolve::DENDROSOLVER_GRID_MIN_Y"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_GRID_MIN_Y")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_GRID_MIN_Y =
                file["dsolve::DENDROSOLVER_GRID_MIN_Y"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_GRID_MAX_Y")) {
            if (0.0 > file["dsolve::DENDROSOLVER_GRID_MAX_Y"].as_floating() ||
                500.0 < file["dsolve::DENDROSOLVER_GRID_MAX_Y"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_GRID_MAX_Y")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_GRID_MAX_Y =
                file["dsolve::DENDROSOLVER_GRID_MAX_Y"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_GRID_MIN_Z")) {
            if (-500.0 >
                    file["dsolve::DENDROSOLVER_GRID_MIN_Z"].as_floating() ||
                0.0 < file["dsolve::DENDROSOLVER_GRID_MIN_Z"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_GRID_MIN_Z")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_GRID_MIN_Z =
                file["dsolve::DENDROSOLVER_GRID_MIN_Z"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_GRID_MAX_Z")) {
            if (0.0 > file["dsolve::DENDROSOLVER_GRID_MAX_Z"].as_floating() ||
                500.0 < file["dsolve::DENDROSOLVER_GRID_MAX_Z"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_GRID_MAX_Z")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_GRID_MAX_Z =
                file["dsolve::DENDROSOLVER_GRID_MAX_Z"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_AMR_R")) {
            if (0.0 > file["dsolve::DENDROSOLVER_BH1_AMR_R"].as_floating() ||
                0.14 < file["dsolve::DENDROSOLVER_BH1_AMR_R"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH1_AMR_R")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH1_AMR_R =
                file["dsolve::DENDROSOLVER_BH1_AMR_R"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_CONSTRAINT_R")) {
            dsolve::DENDROSOLVER_BH1_CONSTRAINT_R =
                file["dsolve::DENDROSOLVER_BH1_CONSTRAINT_R"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_MAX_LEV")) {
            dsolve::DENDROSOLVER_BH1_MAX_LEV =
                file["dsolve::DENDROSOLVER_BH1_MAX_LEV"].as_integer();
        }

        else {
            std::cerr
                << R"(No value for "DENDROSOLVER_BH1_MAX_LEV"; "DENDROSOLVER_BH1_MAX_LEV" must be given a value)"
                << std::endl;
            exit(-1);
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_AMR_R")) {
            if (0.0 > file["dsolve::DENDROSOLVER_BH2_AMR_R"].as_floating() ||
                3.0 < file["dsolve::DENDROSOLVER_BH2_AMR_R"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH2_AMR_R")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH2_AMR_R =
                file["dsolve::DENDROSOLVER_BH2_AMR_R"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_CONSTRAINT_R")) {
            dsolve::DENDROSOLVER_BH2_CONSTRAINT_R =
                file["dsolve::DENDROSOLVER_BH2_CONSTRAINT_R"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_MAX_LEV")) {
            dsolve::DENDROSOLVER_BH2_MAX_LEV =
                file["dsolve::DENDROSOLVER_BH2_MAX_LEV"].as_integer();
        }

        else {
            std::cerr
                << R"(No value for "DENDROSOLVER_BH2_MAX_LEV"; "DENDROSOLVER_BH2_MAX_LEV" must be given a value)"
                << std::endl;
            exit(-1);
        }

        if (file.contains("dsolve::DENDROSOLVER_ETA_POWER")) {
            for (int i = 0; i < 2; ++i) {
                dsolve::DENDROSOLVER_ETA_POWER[i] =
                    file["dsolve::DENDROSOLVER_ETA_POWER"][i].as_floating();
            }
        }

        if (file.contains("dsolve::DENDROSOLVER_XI_0")) {
            dsolve::DENDROSOLVER_XI_0 =
                file["dsolve::DENDROSOLVER_XI_0"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_XI_1")) {
            dsolve::DENDROSOLVER_XI_1 =
                file["dsolve::DENDROSOLVER_XI_1"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_XI_2")) {
            dsolve::DENDROSOLVER_XI_2 =
                file["dsolve::DENDROSOLVER_XI_2"].as_integer();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_MASS")) {
            if (-6.0 > file["dsolve::DENDROSOLVER_BH2_MASS"].as_floating() ||
                6.0 < file["dsolve::DENDROSOLVER_BH2_MASS"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH2_MASS")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH2_MASS =
                file["dsolve::DENDROSOLVER_BH2_MASS"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_X")) {
            if (-6.0 > file["dsolve::DENDROSOLVER_BH2_X"].as_floating() ||
                6.0 < file["dsolve::DENDROSOLVER_BH2_X"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::DENDROSOLVER_BH2_X")"
                          << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH2_X =
                file["dsolve::DENDROSOLVER_BH2_X"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_Y")) {
            if (-6.0 > file["dsolve::DENDROSOLVER_BH2_Y"].as_floating() ||
                6.0 < file["dsolve::DENDROSOLVER_BH2_Y"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::DENDROSOLVER_BH2_Y")"
                          << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH2_Y =
                file["dsolve::DENDROSOLVER_BH2_Y"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_Z")) {
            if (-6.0 > file["dsolve::DENDROSOLVER_BH2_Z"].as_floating() ||
                6.0 < file["dsolve::DENDROSOLVER_BH2_Z"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::DENDROSOLVER_BH2_Z")"
                          << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH2_Z =
                file["dsolve::DENDROSOLVER_BH2_Z"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_V_X")) {
            if (-4.0 > file["dsolve::DENDROSOLVER_BH2_V_X"].as_floating() ||
                0 < file["dsolve::DENDROSOLVER_BH2_V_X"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH2_V_X")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH2_V_X =
                file["dsolve::DENDROSOLVER_BH2_V_X"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_V_Y")) {
            if (-4.0 > file["dsolve::DENDROSOLVER_BH2_V_Y"].as_floating() ||
                4.0 < file["dsolve::DENDROSOLVER_BH2_V_Y"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH2_V_Y")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH2_V_Y =
                file["dsolve::DENDROSOLVER_BH2_V_Y"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_V_Z")) {
            if (-4.0 > file["dsolve::DENDROSOLVER_BH2_V_Z"].as_floating() ||
                4.0 < file["dsolve::DENDROSOLVER_BH2_V_Z"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH2_V_Z")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH2_V_Z =
                file["dsolve::DENDROSOLVER_BH2_V_Z"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_SPIN")) {
            if (0.0 > file["dsolve::DENDROSOLVER_BH2_SPIN"].as_floating() ||
                3.0 < file["dsolve::DENDROSOLVER_BH2_SPIN"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH2_SPIN")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH2_SPIN =
                file["dsolve::DENDROSOLVER_BH2_SPIN"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_SPIN_THETA")) {
            if (0.0 >
                    file["dsolve::DENDROSOLVER_BH2_SPIN_THETA"].as_floating() ||
                3.0 <
                    file["dsolve::DENDROSOLVER_BH2_SPIN_THETA"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH2_SPIN_THETA")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH2_SPIN_THETA =
                file["dsolve::DENDROSOLVER_BH2_SPIN_THETA"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH2_SPIN_PHI")) {
            if (0.0 > file["dsolve::DENDROSOLVER_BH2_SPIN_PHI"].as_floating() ||
                3.0 < file["dsolve::DENDROSOLVER_BH2_SPIN_PHI"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH2_SPIN_PHI")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH2_SPIN_PHI =
                file["dsolve::DENDROSOLVER_BH2_SPIN_PHI"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_MASS")) {
            if (0.0 > file["dsolve::DENDROSOLVER_BH1_MASS"].as_floating() ||
                0.7 < file["dsolve::DENDROSOLVER_BH1_MASS"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH1_MASS")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH1_MASS =
                file["dsolve::DENDROSOLVER_BH1_MASS"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_X")) {
            if (-6.0 > file["dsolve::DENDROSOLVER_BH1_X"].as_floating() ||
                6.0 < file["dsolve::DENDROSOLVER_BH1_X"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::DENDROSOLVER_BH1_X")"
                          << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH1_X =
                file["dsolve::DENDROSOLVER_BH1_X"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_Y")) {
            if (-6.0 > file["dsolve::DENDROSOLVER_BH1_Y"].as_floating() ||
                6.0 < file["dsolve::DENDROSOLVER_BH1_Y"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::DENDROSOLVER_BH1_Y")"
                          << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH1_Y =
                file["dsolve::DENDROSOLVER_BH1_Y"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_Z")) {
            if (-6.0 > file["dsolve::DENDROSOLVER_BH1_Z"].as_floating() ||
                6.0 < file["dsolve::DENDROSOLVER_BH1_Z"].as_floating()) {
                std::cerr << R"(Invalid value for "dsolve::DENDROSOLVER_BH1_Z")"
                          << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH1_Z =
                file["dsolve::DENDROSOLVER_BH1_Z"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_V_X")) {
            if (-4.0 > file["dsolve::DENDROSOLVER_BH1_V_X"].as_floating() ||
                4.0 < file["dsolve::DENDROSOLVER_BH1_V_X"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH1_V_X")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH1_V_X =
                file["dsolve::DENDROSOLVER_BH1_V_X"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_V_Y")) {
            if (-4.0 > file["dsolve::DENDROSOLVER_BH1_V_Y"].as_floating() ||
                4.0 < file["dsolve::DENDROSOLVER_BH1_V_Y"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH1_V_Y")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH1_V_Y =
                file["dsolve::DENDROSOLVER_BH1_V_Y"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_V_Z")) {
            if (-4.0 > file["dsolve::DENDROSOLVER_BH1_V_Z"].as_floating() ||
                4.0 < file["dsolve::DENDROSOLVER_BH1_V_Z"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH1_V_Z")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH1_V_Z =
                file["dsolve::DENDROSOLVER_BH1_V_Z"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_SPIN")) {
            if (0.0 > file["dsolve::DENDROSOLVER_BH1_SPIN"].as_floating() ||
                3.0 < file["dsolve::DENDROSOLVER_BH1_SPIN"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH1_SPIN")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH1_SPIN =
                file["dsolve::DENDROSOLVER_BH1_SPIN"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_SPIN_THETA")) {
            if (0.0 >
                    file["dsolve::DENDROSOLVER_BH1_SPIN_THETA"].as_floating() ||
                3.0 <
                    file["dsolve::DENDROSOLVER_BH1_SPIN_THETA"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH1_SPIN_THETA")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH1_SPIN_THETA =
                file["dsolve::DENDROSOLVER_BH1_SPIN_THETA"].as_floating();
        }

        if (file.contains("dsolve::DENDROSOLVER_BH1_SPIN_PHI")) {
            if (0.0 > file["dsolve::DENDROSOLVER_BH1_SPIN_PHI"].as_floating() ||
                3.0 < file["dsolve::DENDROSOLVER_BH1_SPIN_PHI"].as_floating()) {
                std::cerr
                    << R"(Invalid value for "dsolve::DENDROSOLVER_BH1_SPIN_PHI")"
                    << std::endl;
                exit(-1);
            }

            dsolve::DENDROSOLVER_BH1_SPIN_PHI =
                file["dsolve::DENDROSOLVER_BH1_SPIN_PHI"].as_floating();
        }

        dsolve::BH1 =
            BH(DENDROSOLVER_BH1_MASS, DENDROSOLVER_BH1_X, DENDROSOLVER_BH1_Y,
               DENDROSOLVER_BH1_Z, DENDROSOLVER_BH1_V_X, DENDROSOLVER_BH1_V_Y,
               DENDROSOLVER_BH1_V_Z, DENDROSOLVER_BH1_SPIN,
               DENDROSOLVER_BH1_SPIN_THETA, DENDROSOLVER_BH1_SPIN_PHI);

        dsolve::BH2 =
            BH(DENDROSOLVER_BH2_MASS, DENDROSOLVER_BH2_X, DENDROSOLVER_BH2_Y,
               DENDROSOLVER_BH2_Z, DENDROSOLVER_BH2_V_X, DENDROSOLVER_BH2_V_Y,
               DENDROSOLVER_BH2_V_Z, DENDROSOLVER_BH2_SPIN,
               DENDROSOLVER_BH2_SPIN_THETA, DENDROSOLVER_BH2_SPIN_PHI);
    }

    // Broadcast code to send parameters to other processes
    par::Mpi_Bcast(&(GW::DENDROSOLVER_GW_NUM_RADAII), 1, 0, comm);
    par::Mpi_Bcast(&(GW::DENDROSOLVER_GW_NUM_LMODES), 1, 0, comm);
    MPI_Bcast(&(GW::DENDROSOLVER_GW_RADAII), 6, MPI_INT, 0, comm);
    MPI_Bcast(&(GW::DENDROSOLVER_GW_L_MODES), 3, MPI_INT, 0, comm);
    par::Mpi_Bcast(&(BHLOC::EXTRACTION_VAR_ID), 1, 0, comm);
    par::Mpi_Bcast(&(BHLOC::EXTRACTION_TOL), 1, 0, comm);
    MPI_Bcast(const_cast<char*>(TPID::FILE_PREFIX.c_str()),
              TPID::FILE_PREFIX.size() + 1, MPI_CHAR, 0, comm);
    par::Mpi_Bcast(&(TPID::par_b), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::par_m_plus), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::par_m_minus), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::target_M_plus), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::target_M_minus), 1, 0, comm);
    MPI_Bcast(&(TPID::par_P_plus), 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&(TPID::par_P_minus), 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&(TPID::par_S_plus), 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&(TPID::par_S_minus), 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&(TPID::center_offset), 3, MPI_DOUBLE, 0, comm);
    par::Mpi_Bcast(&(TPID::give_bare_mass), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::initial_lapse), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::grid_setup_method), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::solve_momentum_constraint), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::verbose), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::adm_tol), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::Newton_tol), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::initial_lapse_psi_exponent), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::npoints_A), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::npoints_B), 1, 0, comm);
    par::Mpi_Bcast(&(TPID::npoints_phi), 1, 0, comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_LAMBDA), 4, MPI_INT, 0, comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_ETA), 2, MPI_DOUBLE, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_ETADAMP), 1, 0, comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_ALPHA_THEORY), 2, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_LF), 2, MPI_DOUBLE, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_P_EXPO), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_ELE_ORDER), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_PADDING_WIDTH), 1, 0, comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_COMPD_MIN), 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_COMPD_MAX), 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_OCTREE_MIN), 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_OCTREE_MAX), 3, MPI_DOUBLE, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_IO_OUTPUT_FREQ), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_GW_EXTRACT_FREQ), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_GW_EXTRACT_FREQ_AFTER_MERGER), 1, 0,
                   comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_TIME_STEP_OUTPUT_FREQ), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_REMESH_TEST_FREQ), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_REMESH_TEST_FREQ_AFTER_MERGER), 1, 0,
                   comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_CHECKPT_FREQ), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_RESTORE_SOLVER), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_ENABLE_BLOCK_ADAPTIVITY), 1, 0, comm);
    MPI_Bcast(const_cast<char*>(dsolve::DENDROSOLVER_VTU_FILE_PREFIX.c_str()),
              dsolve::DENDROSOLVER_VTU_FILE_PREFIX.size() + 1, MPI_CHAR, 0,
              comm);
    MPI_Bcast(const_cast<char*>(dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX.c_str()),
              dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX.size() + 1, MPI_CHAR, 0,
              comm);
    MPI_Bcast(
        const_cast<char*>(dsolve::DENDROSOLVER_PROFILE_FILE_PREFIX.c_str()),
        dsolve::DENDROSOLVER_PROFILE_FILE_PREFIX.size() + 1, MPI_CHAR, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_NUM_REFINE_VARS), 1, 0, comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_REFINE_VARIABLE_INDICES), 36, MPI_INT, 0,
              comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_NUM_EVOL_VARS_VTU_OUTPUT), 1, 0,
                   comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_NUM_CONST_VARS_VTU_OUTPUT), 1, 0,
                   comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_VTU_OUTPUT_EVOL_INDICES), 36, MPI_INT, 0,
              comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_VTU_OUTPUT_CONST_INDICES), 6, MPI_INT, 0,
              comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_IO_OUTPUT_GAP), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_DENDRO_GRAIN_SZ), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_DENDRO_AMR_FAC), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_INIT_GRID_ITER), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_SPLIT_FIX), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_CFL_FACTOR), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_RK_TIME_BEGIN), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_RK_TIME_END), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_RK_TYPE), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_RK45_DESIRED_TOL), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DISSIPATION_TYPE), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_DISSIPATION_NC), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_DISSIPATION_S), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_LTS_TS_OFFSET), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_MERGED_CHKPT_WRITTEN), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_EH_REFINE_VAL), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_EH_COARSEN_VAL), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_VTU_Z_SLICE_ONLY), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_ASYNC_COMM_K), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_LOAD_IMB_TOL), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_DIM), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_MAXDEPTH), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_MINDEPTH), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_WAVELET_TOL), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_GW_REFINE_WTOL), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_USE_WAVELET_TOL_FUNCTION), 1, 0,
                   comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_WAVELET_TOL_MAX), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R0), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R1), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_USE_FD_GRID_TRANSFER), 1, 0, comm);
    par::Mpi_Bcast((int*)&dsolve::DENDROSOLVER_REFINEMENT_MODE, 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BLK_MIN_X), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BLK_MIN_Y), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BLK_MIN_Z), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BLK_MAX_X), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BLK_MAX_Y), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BLK_MAX_Z), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::ETA_CONST), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::ETA_R0), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::ETA_DAMPING), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::ETA_DAMPING_EXP), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::ANG_PAR), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::CHI_FLOOR), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_TRK0), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::KO_DISS_SIGMA), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_ETA_R0), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_ID_TYPE), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_GRID_MIN_X), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_GRID_MAX_X), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_GRID_MIN_Y), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_GRID_MAX_Y), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_GRID_MIN_Z), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_GRID_MAX_Z), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_AMR_R), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_CONSTRAINT_R), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_MAX_LEV), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_AMR_R), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_CONSTRAINT_R), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_MAX_LEV), 1, 0, comm);
    MPI_Bcast(&(dsolve::DENDROSOLVER_ETA_POWER), 2, MPI_DOUBLE, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_XI_0), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_XI_1), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_XI_2), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_MASS), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_X), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_Y), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_Z), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_V_X), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_V_Y), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_V_Z), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_SPIN), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_SPIN_THETA), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH2_SPIN_PHI), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_MASS), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_X), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_Y), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_Z), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_V_X), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_V_Y), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_V_Z), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_SPIN), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_SPIN_THETA), 1, 0, comm);
    par::Mpi_Bcast(&(dsolve::DENDROSOLVER_BH1_SPIN_PHI), 1, 0, comm);
}

void dumpParamFile(std::ostream& sout, int root, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if (rank == root) {
        sout << "Parameters read: " << std::endl;
        sout << "\tGW::DENDROSOLVER_GW_NUM_RADAII: "
             << GW::DENDROSOLVER_GW_NUM_RADAII << std::endl;
        sout << "\tGW::DENDROSOLVER_GW_NUM_LMODES: "
             << GW::DENDROSOLVER_GW_NUM_LMODES << std::endl;
        sout << "\tGW::DENDROSOLVER_GW_RADAII: [";
        for (unsigned int i = 0; i < 6; ++i) {
            sout << GW::DENDROSOLVER_GW_RADAII[i] << (i < 6 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tGW::DENDROSOLVER_GW_L_MODES: [";
        for (unsigned int i = 0; i < 3; ++i) {
            sout << GW::DENDROSOLVER_GW_L_MODES[i] << (i < 3 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tGW::DENDROSOLVER_GW_OUTPUT_PRECISION: "
             << GW::DENDROSOLVER_GW_OUTPUT_PRECISION << std::endl;
        sout << "\tBHLOC::EXTRACTION_VAR_ID: " << BHLOC::EXTRACTION_VAR_ID
             << std::endl;
        sout << "\tBHLOC::EXTRACTION_TOL: " << BHLOC::EXTRACTION_TOL
             << std::endl;
        sout << "\tTPID::FILE_PREFIX: " << TPID::FILE_PREFIX << std::endl;
        sout << "\tTPID::TP_epsilon: " << TPID::TP_epsilon << std::endl;
        sout << "\tTPID::swap_xz: " << TPID::swap_xz << std::endl;
        sout << "\tTPID::use_sources: " << TPID::use_sources << std::endl;
        sout << "\tTPID::rescale_sources: " << TPID::rescale_sources
             << std::endl;
        sout << "\tTPID::use_external_initial_guess: "
             << TPID::use_external_initial_guess << std::endl;
        sout << "\tTPID::do_residuum_debug_output: "
             << TPID::do_residuum_debug_output << std::endl;
        sout << "\tTPID::do_initial_debug_output: "
             << TPID::do_initial_debug_output << std::endl;
        sout << "\tTPID::multiply_old_lapse: " << TPID::multiply_old_lapse
             << std::endl;
        sout << "\tTPID::TP_Tiny: " << TPID::TP_Tiny << std::endl;
        sout << "\tTPID::TP_Extend_Radius: " << TPID::TP_Extend_Radius
             << std::endl;
        sout << "\tTPID::Newton_maxit: " << TPID::Newton_maxit << std::endl;
        sout << "\tTPID::par_b: " << TPID::par_b << std::endl;
        sout << "\tTPID::par_m_plus: " << TPID::par_m_plus << std::endl;
        sout << "\tTPID::par_m_minus: " << TPID::par_m_minus << std::endl;
        sout << "\tTPID::target_M_plus: " << TPID::target_M_plus << std::endl;
        sout << "\tTPID::target_M_minus: " << TPID::target_M_minus << std::endl;
        sout << "\tTPID::par_P_plus: [";
        for (unsigned int i = 0; i < 3; ++i) {
            sout << TPID::par_P_plus[i] << (i < 3 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tTPID::par_P_minus: [";
        for (unsigned int i = 0; i < 3; ++i) {
            sout << TPID::par_P_minus[i] << (i < 3 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tTPID::par_S_plus: [";
        for (unsigned int i = 0; i < 3; ++i) {
            sout << TPID::par_S_plus[i] << (i < 3 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tTPID::par_S_minus: [";
        for (unsigned int i = 0; i < 3; ++i) {
            sout << TPID::par_S_minus[i] << (i < 3 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tTPID::center_offset: [";
        for (unsigned int i = 0; i < 3; ++i) {
            sout << TPID::center_offset[i] << (i < 3 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tTPID::give_bare_mass: " << TPID::give_bare_mass << std::endl;
        sout << "\tTPID::initial_lapse: " << TPID::initial_lapse << std::endl;
        sout << "\tTPID::grid_setup_method: " << TPID::grid_setup_method
             << std::endl;
        sout << "\tTPID::solve_momentum_constraint: "
             << TPID::solve_momentum_constraint << std::endl;
        sout << "\tTPID::verbose: " << TPID::verbose << std::endl;
        sout << "\tTPID::adm_tol: " << TPID::adm_tol << std::endl;
        sout << "\tTPID::Newton_tol: " << TPID::Newton_tol << std::endl;
        sout << "\tTPID::initial_lapse_psi_exponent: "
             << TPID::initial_lapse_psi_exponent << std::endl;
        sout << "\tTPID::npoints_A: " << TPID::npoints_A << std::endl;
        sout << "\tTPID::npoints_B: " << TPID::npoints_B << std::endl;
        sout << "\tTPID::npoints_phi: " << TPID::npoints_phi << std::endl;
        sout << "\tdsolve::DENDRO_VERSION: " << dsolve::DENDRO_VERSION
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_LAMBDA: [";
        for (unsigned int i = 0; i < 4; ++i) {
            sout << dsolve::DENDROSOLVER_LAMBDA[i] << (i < 4 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_ETA: [";
        for (unsigned int i = 0; i < 2; ++i) {
            sout << dsolve::DENDROSOLVER_ETA[i] << (i < 2 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_ETADAMP: "
             << dsolve::DENDROSOLVER_ETADAMP << std::endl;
        sout << "\tdsolve::DENDROSOLVER_ALPHA_THEORY: [";
        for (unsigned int i = 0; i < 2; ++i) {
            sout << dsolve::DENDROSOLVER_ALPHA_THEORY[i]
                 << (i < 2 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_LF: [";
        for (unsigned int i = 0; i < 2; ++i) {
            sout << dsolve::DENDROSOLVER_LF[i] << (i < 2 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_P_EXPO: " << dsolve::DENDROSOLVER_P_EXPO
             << std::endl;
        sout << "\tdsolve::ALPHA_FLOOR: " << dsolve::ALPHA_FLOOR << std::endl;
        sout << "\tdsolve::DENDROSOLVER_ELE_ORDER: "
             << dsolve::DENDROSOLVER_ELE_ORDER << std::endl;
        sout << "\tdsolve::DENDROSOLVER_PADDING_WIDTH: "
             << dsolve::DENDROSOLVER_PADDING_WIDTH << std::endl;
        sout << "\tdsolve::DENDROSOLVER_NUM_VARS: "
             << dsolve::DENDROSOLVER_NUM_VARS << std::endl;
        sout << "\tdsolve::DENDROSOLVER_CONSTRAINT_NUM_VARS: "
             << dsolve::DENDROSOLVER_CONSTRAINT_NUM_VARS << std::endl;
        sout << "\tdsolve::DENDROSOLVER_RK45_STAGES: "
             << dsolve::DENDROSOLVER_RK45_STAGES << std::endl;
        sout << "\tdsolve::DENDROSOLVER_RK4_STAGES: "
             << dsolve::DENDROSOLVER_RK4_STAGES << std::endl;
        sout << "\tdsolve::DENDROSOLVER_RK3_STAGES: "
             << dsolve::DENDROSOLVER_RK3_STAGES << std::endl;
        sout << "\tdsolve::DENDROSOLVER_SAFETY_FAC: "
             << dsolve::DENDROSOLVER_SAFETY_FAC << std::endl;
        sout << "\tdsolve::DENDROSOLVER_NUM_VARS_INTENL: "
             << dsolve::DENDROSOLVER_NUM_VARS_INTENL << std::endl;
        sout << "\tdsolve::DENDROSOLVER_COMPD_MIN: [";
        for (unsigned int i = 0; i < 3; ++i) {
            sout << dsolve::DENDROSOLVER_COMPD_MIN[i]
                 << (i < 3 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_COMPD_MAX: [";
        for (unsigned int i = 0; i < 3; ++i) {
            sout << dsolve::DENDROSOLVER_COMPD_MAX[i]
                 << (i < 3 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_OCTREE_MIN: [";
        for (unsigned int i = 0; i < 3; ++i) {
            sout << dsolve::DENDROSOLVER_OCTREE_MIN[i]
                 << (i < 3 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_OCTREE_MAX: [";
        for (unsigned int i = 0; i < 3; ++i) {
            sout << dsolve::DENDROSOLVER_OCTREE_MAX[i]
                 << (i < 3 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_IO_OUTPUT_FREQ: "
             << dsolve::DENDROSOLVER_IO_OUTPUT_FREQ << std::endl;
        sout << "\tdsolve::DENDROSOLVER_GW_EXTRACT_FREQ: "
             << dsolve::DENDROSOLVER_GW_EXTRACT_FREQ << std::endl;
        sout << "\tdsolve::DENDROSOLVER_GW_EXTRACT_FREQ_AFTER_MERGER: "
             << dsolve::DENDROSOLVER_GW_EXTRACT_FREQ_AFTER_MERGER << std::endl;
        sout << "\tdsolve::DENDROSOLVER_TIME_STEP_OUTPUT_FREQ: "
             << dsolve::DENDROSOLVER_TIME_STEP_OUTPUT_FREQ << std::endl;
        sout << "\tdsolve::DENDROSOLVER_REMESH_TEST_FREQ: "
             << dsolve::DENDROSOLVER_REMESH_TEST_FREQ << std::endl;
        sout << "\tdsolve::DENDROSOLVER_REMESH_TEST_FREQ_AFTER_MERGER: "
             << dsolve::DENDROSOLVER_REMESH_TEST_FREQ_AFTER_MERGER << std::endl;
        sout << "\tdsolve::DENDROSOLVER_CHECKPT_FREQ: "
             << dsolve::DENDROSOLVER_CHECKPT_FREQ << std::endl;
        sout << "\tdsolve::DENDROSOLVER_RESTORE_SOLVER: "
             << dsolve::DENDROSOLVER_RESTORE_SOLVER << std::endl;
        sout << "\tdsolve::DENDROSOLVER_ENABLE_BLOCK_ADAPTIVITY: "
             << dsolve::DENDROSOLVER_ENABLE_BLOCK_ADAPTIVITY << std::endl;
        sout << "\tdsolve::DENDROSOLVER_VTU_FILE_PREFIX: "
             << dsolve::DENDROSOLVER_VTU_FILE_PREFIX << std::endl;
        sout << "\tdsolve::DENDROSOLVER_CHKPT_FILE_PREFIX: "
             << dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX << std::endl;
        sout << "\tdsolve::DENDROSOLVER_PROFILE_FILE_PREFIX: "
             << dsolve::DENDROSOLVER_PROFILE_FILE_PREFIX << std::endl;
        sout << "\tdsolve::DENDROSOLVER_NUM_REFINE_VARS: "
             << dsolve::DENDROSOLVER_NUM_REFINE_VARS << std::endl;
        sout << "\tdsolve::DENDROSOLVER_REFINE_VARIABLE_INDICES: [";
        for (unsigned int i = 0; i < 36; ++i) {
            sout << dsolve::DENDROSOLVER_REFINE_VARIABLE_INDICES[i]
                 << (i < 36 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_NUM_EVOL_VARS_VTU_OUTPUT: "
             << dsolve::DENDROSOLVER_NUM_EVOL_VARS_VTU_OUTPUT << std::endl;
        sout << "\tdsolve::DENDROSOLVER_NUM_CONST_VARS_VTU_OUTPUT: "
             << dsolve::DENDROSOLVER_NUM_CONST_VARS_VTU_OUTPUT << std::endl;
        sout << "\tdsolve::DENDROSOLVER_VTU_OUTPUT_EVOL_INDICES: [";
        for (unsigned int i = 0; i < 36; ++i) {
            sout << dsolve::DENDROSOLVER_VTU_OUTPUT_EVOL_INDICES[i]
                 << (i < 36 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_VTU_OUTPUT_CONST_INDICES: [";
        for (unsigned int i = 0; i < 6; ++i) {
            sout << dsolve::DENDROSOLVER_VTU_OUTPUT_CONST_INDICES[i]
                 << (i < 6 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_IO_OUTPUT_GAP: "
             << dsolve::DENDROSOLVER_IO_OUTPUT_GAP << std::endl;
        sout << "\tdsolve::DENDROSOLVER_DENDRO_GRAIN_SZ: "
             << dsolve::DENDROSOLVER_DENDRO_GRAIN_SZ << std::endl;
        sout << "\tdsolve::DENDROSOLVER_DENDRO_AMR_FAC: "
             << dsolve::DENDROSOLVER_DENDRO_AMR_FAC << std::endl;
        sout << "\tdsolve::DENDROSOLVER_INIT_GRID_ITER: "
             << dsolve::DENDROSOLVER_INIT_GRID_ITER << std::endl;
        sout << "\tdsolve::DENDROSOLVER_SPLIT_FIX: "
             << dsolve::DENDROSOLVER_SPLIT_FIX << std::endl;
        sout << "\tdsolve::DENDROSOLVER_CFL_FACTOR: "
             << dsolve::DENDROSOLVER_CFL_FACTOR << std::endl;
        sout << "\tdsolve::DENDROSOLVER_RK_TIME_BEGIN: "
             << dsolve::DENDROSOLVER_RK_TIME_BEGIN << std::endl;
        sout << "\tdsolve::DENDROSOLVER_RK_TIME_END: "
             << dsolve::DENDROSOLVER_RK_TIME_END << std::endl;
        sout << "\tdsolve::DENDROSOLVER_RK_TYPE: "
             << dsolve::DENDROSOLVER_RK_TYPE << std::endl;
        sout << "\tdsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE: "
             << dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE << std::endl;
        sout << "\tdsolve::DENDROSOLVER_RK45_DESIRED_TOL: "
             << dsolve::DENDROSOLVER_RK45_DESIRED_TOL << std::endl;
        sout << "\tdsolve::DISSIPATION_TYPE: " << dsolve::DISSIPATION_TYPE
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_DISSIPATION_NC: "
             << dsolve::DENDROSOLVER_DISSIPATION_NC << std::endl;
        sout << "\tdsolve::DENDROSOLVER_DISSIPATION_S: "
             << dsolve::DENDROSOLVER_DISSIPATION_S << std::endl;
        sout << "\tdsolve::DENDROSOLVER_LTS_TS_OFFSET: "
             << dsolve::DENDROSOLVER_LTS_TS_OFFSET << std::endl;
        sout << "\tdsolve::DENDROSOLVER_MERGED_CHKPT_WRITTEN: "
             << dsolve::DENDROSOLVER_MERGED_CHKPT_WRITTEN << std::endl;
        sout << "\tdsolve::DENDROSOLVER_EH_REFINE_VAL: "
             << dsolve::DENDROSOLVER_EH_REFINE_VAL << std::endl;
        sout << "\tdsolve::DENDROSOLVER_EH_COARSEN_VAL: "
             << dsolve::DENDROSOLVER_EH_COARSEN_VAL << std::endl;
        sout << "\tdsolve::DENDROSOLVER_VTU_Z_SLICE_ONLY: "
             << dsolve::DENDROSOLVER_VTU_Z_SLICE_ONLY << std::endl;
        sout << "\tdsolve::DENDROSOLVER_ASYNC_COMM_K: "
             << dsolve::DENDROSOLVER_ASYNC_COMM_K << std::endl;
        sout << "\tdsolve::DENDROSOLVER_LOAD_IMB_TOL: "
             << dsolve::DENDROSOLVER_LOAD_IMB_TOL << std::endl;
        sout << "\tdsolve::DENDROSOLVER_DIM: " << dsolve::DENDROSOLVER_DIM
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_MAXDEPTH: "
             << dsolve::DENDROSOLVER_MAXDEPTH << std::endl;
        sout << "\tdsolve::DENDROSOLVER_MINDEPTH: "
             << dsolve::DENDROSOLVER_MINDEPTH << std::endl;
        sout << "\tdsolve::DENDROSOLVER_WAVELET_TOL: "
             << dsolve::DENDROSOLVER_WAVELET_TOL << std::endl;
        sout << "\tdsolve::DENDROSOLVER_GW_REFINE_WTOL: "
             << dsolve::DENDROSOLVER_GW_REFINE_WTOL << std::endl;
        sout << "\tdsolve::DENDROSOLVER_USE_WAVELET_TOL_FUNCTION: "
             << dsolve::DENDROSOLVER_USE_WAVELET_TOL_FUNCTION << std::endl;
        sout << "\tdsolve::DENDROSOLVER_WAVELET_TOL_MAX: "
             << dsolve::DENDROSOLVER_WAVELET_TOL_MAX << std::endl;
        sout << "\tdsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R0: "
             << dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R0 << std::endl;
        sout << "\tdsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R1: "
             << dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R1 << std::endl;
        sout << "\tdsolve::DENDROSOLVER_USE_FD_GRID_TRANSFER: "
             << dsolve::DENDROSOLVER_USE_FD_GRID_TRANSFER << std::endl;
        sout << "\tdsolve::DENDROSOLVER_REFINEMENT_MODE: "
             << dsolve::DENDROSOLVER_REFINEMENT_MODE << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BLK_MIN_X: "
             << dsolve::DENDROSOLVER_BLK_MIN_X << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BLK_MIN_Y: "
             << dsolve::DENDROSOLVER_BLK_MIN_Y << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BLK_MIN_Z: "
             << dsolve::DENDROSOLVER_BLK_MIN_Z << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BLK_MAX_X: "
             << dsolve::DENDROSOLVER_BLK_MAX_X << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BLK_MAX_Y: "
             << dsolve::DENDROSOLVER_BLK_MAX_Y << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BLK_MAX_Z: "
             << dsolve::DENDROSOLVER_BLK_MAX_Z << std::endl;
        sout << "\tdsolve::ETA_CONST: " << dsolve::ETA_CONST << std::endl;
        sout << "\tdsolve::ETA_R0: " << dsolve::ETA_R0 << std::endl;
        sout << "\tdsolve::ETA_DAMPING: " << dsolve::ETA_DAMPING << std::endl;
        sout << "\tdsolve::ETA_DAMPING_EXP: " << dsolve::ETA_DAMPING_EXP
             << std::endl;
        sout << "\tdsolve::ANG_PAR: " << dsolve::ANG_PAR << std::endl;
        sout << "\tdsolve::CHI_FLOOR: " << dsolve::CHI_FLOOR << std::endl;
        sout << "\tdsolve::DENDROSOLVER_TRK0: " << dsolve::DENDROSOLVER_TRK0
             << std::endl;
        sout << "\tdsolve::KO_DISS_SIGMA: " << dsolve::KO_DISS_SIGMA
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_ETA_R0: " << dsolve::DENDROSOLVER_ETA_R0
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_ID_TYPE: "
             << dsolve::DENDROSOLVER_ID_TYPE << std::endl;
        sout << "\tdsolve::DENDROSOLVER_GRID_MIN_X: "
             << dsolve::DENDROSOLVER_GRID_MIN_X << std::endl;
        sout << "\tdsolve::DENDROSOLVER_GRID_MAX_X: "
             << dsolve::DENDROSOLVER_GRID_MAX_X << std::endl;
        sout << "\tdsolve::DENDROSOLVER_GRID_MIN_Y: "
             << dsolve::DENDROSOLVER_GRID_MIN_Y << std::endl;
        sout << "\tdsolve::DENDROSOLVER_GRID_MAX_Y: "
             << dsolve::DENDROSOLVER_GRID_MAX_Y << std::endl;
        sout << "\tdsolve::DENDROSOLVER_GRID_MIN_Z: "
             << dsolve::DENDROSOLVER_GRID_MIN_Z << std::endl;
        sout << "\tdsolve::DENDROSOLVER_GRID_MAX_Z: "
             << dsolve::DENDROSOLVER_GRID_MAX_Z << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_AMR_R: "
             << dsolve::DENDROSOLVER_BH1_AMR_R << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_CONSTRAINT_R: "
             << dsolve::DENDROSOLVER_BH1_CONSTRAINT_R << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_MAX_LEV: "
             << dsolve::DENDROSOLVER_BH1_MAX_LEV << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_AMR_R: "
             << dsolve::DENDROSOLVER_BH2_AMR_R << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_CONSTRAINT_R: "
             << dsolve::DENDROSOLVER_BH2_CONSTRAINT_R << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_MAX_LEV: "
             << dsolve::DENDROSOLVER_BH2_MAX_LEV << std::endl;
        sout << "\tdsolve::DENDROSOLVER_ETA_POWER: [";
        for (unsigned int i = 0; i < 2; ++i) {
            sout << dsolve::DENDROSOLVER_ETA_POWER[i]
                 << (i < 2 - 1 ? ',' : ']');
        }
        sout << std::endl;
        sout << "\tdsolve::DENDROSOLVER_XI_0: " << dsolve::DENDROSOLVER_XI_0
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_XI_1: " << dsolve::DENDROSOLVER_XI_1
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_XI_2: " << dsolve::DENDROSOLVER_XI_2
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_MASS: "
             << dsolve::DENDROSOLVER_BH2_MASS << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_X: " << dsolve::DENDROSOLVER_BH2_X
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_Y: " << dsolve::DENDROSOLVER_BH2_Y
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_Z: " << dsolve::DENDROSOLVER_BH2_Z
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_V_X: "
             << dsolve::DENDROSOLVER_BH2_V_X << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_V_Y: "
             << dsolve::DENDROSOLVER_BH2_V_Y << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_V_Z: "
             << dsolve::DENDROSOLVER_BH2_V_Z << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_SPIN: "
             << dsolve::DENDROSOLVER_BH2_SPIN << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_SPIN_THETA: "
             << dsolve::DENDROSOLVER_BH2_SPIN_THETA << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH2_SPIN_PHI: "
             << dsolve::DENDROSOLVER_BH2_SPIN_PHI << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_MASS: "
             << dsolve::DENDROSOLVER_BH1_MASS << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_X: " << dsolve::DENDROSOLVER_BH1_X
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_Y: " << dsolve::DENDROSOLVER_BH1_Y
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_Z: " << dsolve::DENDROSOLVER_BH1_Z
             << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_V_X: "
             << dsolve::DENDROSOLVER_BH1_V_X << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_V_Y: "
             << dsolve::DENDROSOLVER_BH1_V_Y << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_V_Z: "
             << dsolve::DENDROSOLVER_BH1_V_Z << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_SPIN: "
             << dsolve::DENDROSOLVER_BH1_SPIN << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_SPIN_THETA: "
             << dsolve::DENDROSOLVER_BH1_SPIN_THETA << std::endl;
        sout << "\tdsolve::DENDROSOLVER_BH1_SPIN_PHI: "
             << dsolve::DENDROSOLVER_BH1_SPIN_PHI << std::endl;
    }
}
}  // namespace dsolve

//[[[end]]]

// end parameters.cpp file