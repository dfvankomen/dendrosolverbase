/**
 * This is the main file for Custom Formulation of GR Simulation
 *
 * @author David Van Komen, Milinda Fernando
 * @brief Formulation of GR simulations, this houses "main"
 *
 * :::License:::
 */

#include "include/solver_main.h"

#include <iostream>
#include <vector>

#include "TreeNode.h"
#include "include/grUtils.h"
#include "include/rkSolver.h"
#include "mesh.h"
#include "meshUtils.h"
#include "mpi.h"
#include "octUtils.h"

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "No parameter file was given, exiting..." << std::endl;
        std::cout << "Usage: " << argv[0] << " paramFile" << std::endl;
        exit(0);
    }

    // seed the randomness used later
    srand(static_cast<unsigned>(time(0)));

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    // initialize the flops timer
    dsolve::timer::initFlops();

    // begin the imer for total run time
    dsolve::timer::total_runtime.start();

    // CMake options to be printed:::
    if (!rank) {
#ifdef DENDROSOLVER_COMPUTE_CONSTRAINTS
        std::cout << GRN << "  Compiled with DENDROSOLVER_COMPUTE_CONSTRAINTS"
                  << NRM << std::endl;
#else
        std::cout << RED
                  << "  Compiled without DENDROSOLVER_COMPUTE_CONSTRAINTS"
                  << NRM << std::endl;
#endif
#ifdef DENDROSOLVER_ENABLE_VTU_CONSTRAINT_OUTPUT
        std::cout << GRN
                  << "  Compiled with DENDROSOLVER_ENABLE_VTU_CONSTRAINT_OUTPUT"
                  << NRM << std::endl;
#else
        std::cout
            << RED
            << "  Compiled without DENDROSOLVER_ENABLE_VTU_CONSTRAINT_OUTPUT"
            << NRM << std::endl;
#endif
#ifdef DENDROSOLVER_ENABLE_VTU_OUTPUT
        std::cout << GRN << "  Compiled with DENDROSOLVER_ENABLE_VTU_OUTPUT"
                  << NRM << std::endl;
#else
        std::cout << RED << "  Compiled without DENDROSOLVER_ENABLE_VTU_OUTPUT"
                  << NRM << std::endl;
#endif
#ifdef DENDROSOLVER_ETA_FUNCTION
        std::cout << GRN << "  Compiled with  DENDROSOLVER_ETA_FUNCTION" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  DENDROSOLVER_ETA_FUNCTION"
                  << NRM << std::endl;
#endif
#ifdef DENDROSOLVER_EXTRACT_BH_LOCATIONS
        std::cout << GRN << "  Compiled with  DENDROSOLVER_EXTRACT_BH_LOCATIONS"
                  << NRM << std::endl;
#else
        std::cout << RED
                  << "  Compiled without  DENDROSOLVER_EXTRACT_BH_LOCATIONS"
                  << NRM << std::endl;
#endif
#ifdef DENDROSOLVER_EXTRACT_GRAVITATIONAL_WAVES
        std::cout << GRN
                  << "  Compiled with  DENDROSOLVER_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#else
        std::cout
            << RED
            << "  Compiled without  DENDROSOLVER_EXTRACT_GRAVITATIONAL_WAVES"
            << NRM << std::endl;
#endif
#ifdef DENDROSOLVER_GAUGE_ROCHESTER
        std::cout << GRN << "  Compiled with  DENDROSOLVER_GAUGE_ROCHESTER"
                  << NRM << std::endl;
#else
        std::cout << RED << "  Compiled without  DENDROSOLVER_GAUGE_ROCHESTER"
                  << NRM << std::endl;
#endif
#ifdef DENDROSOLVER_KERR_SCHILD_TEST
        std::cout << GRN << "  Compiled with  DENDROSOLVER_KERR_SCHILD_TEST"
                  << NRM << std::endl;
#else
        std::cout << RED << "  Compiled without  DENDROSOLVER_KERR_SCHILD_TEST"
                  << NRM << std::endl;
#endif

#ifdef DENDROSOLVER_REFINE_BASE_EH
        std::cout << GRN << "  Compiled with  DENDROSOLVER_REFINE_BASE_EH"
                  << NRM << std::endl;
#else
        std::cout << RED << "  Compiled without  DENDROSOLVER_REFINE_BASE_EH"
                  << NRM << std::endl;
#endif

#ifdef USE_FD_INTERP_FOR_UNZIP
        std::cout << GRN << "  Compiled with  USE_FD_INTERP_FOR_UNZIP" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  USE_FD_INTERP_FOR_UNZIP" << NRM
                  << std::endl;
#endif

#ifdef DENDROSOLVER_USE_4TH_ORDER_DERIVS
        std::cout << GRN << "  Using 4th order FD stencils" << NRM << std::endl;
#endif

#ifdef DENDROSOLVER_USE_6TH_ORDER_DERIVS
        std::cout << GRN << "  Using 6th order FD stencils" << NRM << std::endl;
#endif

#ifdef DENDROSOLVER_USE_8TH_ORDER_DERIVS
        std::cout << GRN << "  Using 8th order FD stencils" << NRM << std::endl;
#endif
    }

    /**
     * STEP 1
     *
     * Read in the Parameter File and check initialization parameters for
     * problems.
     */
    if (!rank) {
        std::cout << " reading parameter file :" << argv[1] << std::endl;
    }
    dsolve::readParamFile(argv[1], comm);

    int root = std::min(1, npes - 1);
    // dump the parameter file
    dsolve::dumpParamFile(std::cout, root, comm);

    // dump parameter data to individual files for sanity
    // uncomment the following to test the dump for each individual process!
#if 1
    for (int ifile = 0; ifile < npes; ifile++) {
        if (rank == ifile) {
            std::ofstream tempfile;
            tempfile.open("dumped_params_process_" + std::to_string(ifile) +
                          ".txt");
            dsolve::dumpParamFile(tempfile, ifile, comm);
            tempfile.close();
        }
    }
#endif

    _InitializeHcurve(dsolve::DENDROSOLVER_DIM);
    m_uiMaxDepth = dsolve::DENDROSOLVER_MAXDEPTH;

    if (dsolve::DENDROSOLVER_NUM_VARS % dsolve::DENDROSOLVER_ASYNC_COMM_K !=
        0) {
        if (!rank)
            std::cout << "[overlap communication error]: total "
                         "DENDROSOLVER_NUM_VARS: "
                      << dsolve::DENDROSOLVER_NUM_VARS
                      << " is not divisable by DENDROSOLVER_ASYNC_COMM_K: "
                      << dsolve::DENDROSOLVER_ASYNC_COMM_K << std::endl;
        MPI_Abort(comm, 0);
    }

    if (dsolve::DENDROSOLVER_GW_EXTRACT_FREQ >
        dsolve::DENDROSOLVER_IO_OUTPUT_FREQ) {
        if (!rank)
            std::cout << " DENDROSOLVER_GW_EXTRACT_FREQ should be less than "
                         "DENDROSOLVER_IO_OUTPUT_FREQ "
                      << std::endl;
        MPI_Abort(comm, 0);
    }

    /**
     * STEP 2
     *
     * Generate the initial grid
     * Either through adaptive blocking or through adaptive mesh refinement
     */
    std::vector<ot::TreeNode> tmpNodes;
    // NOTE: initial grid function needs to do a conversion from grid format to
    // physical format
    std::function<void(double, double, double, double *)> f_init =
        [](double x, double y, double z, double *var) {
            dsolve::initDataFuncToPhysCoords(x, y, z, var);
        };
    // std::function<double(double, double, double)> f_init_alpha = [](double x,
    // double y, double z) { double var[24]; dsolve::punctureData(x,y,z,var);
    // return var[0]; };
    // std::function<void(double,double,double,double*)> f_init=[](double
    // x,double y,double z,double*var){dsolve::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars = dsolve::DENDROSOLVER_NUM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < dsolve::DENDROSOLVER_NUM_VARS; i++)
        varIndex[i] = i;

    DendroIntL localSz, globalSz;
    double t_stat;
    double t_stat_g[3];

    dsolve::timer::t_f2o.start();

    if (dsolve::DENDROSOLVER_ENABLE_BLOCK_ADAPTIVITY) {
        if (!rank)
            std::cout << YLW << "Using block adaptive mesh. AMR disabled "
                      << NRM << std::endl;
        const Point pt_min(dsolve::DENDROSOLVER_BLK_MIN_X,
                           dsolve::DENDROSOLVER_BLK_MIN_Y,
                           dsolve::DENDROSOLVER_BLK_MIN_Z);
        const Point pt_max(dsolve::DENDROSOLVER_BLK_MAX_X,
                           dsolve::DENDROSOLVER_BLK_MAX_Y,
                           dsolve::DENDROSOLVER_BLK_MAX_Z);

        dsolve::blockAdaptiveOctree(tmpNodes, pt_min, pt_max, m_uiMaxDepth - 2,
                                    m_uiMaxDepth, comm);
    } else {
        if (!rank)
            std::cout << YLW << "Using function2Octree. AMR enabled " << NRM
                      << std::endl;
        const unsigned int f2olmin = std::min(dsolve::DENDROSOLVER_BH1_MAX_LEV,
                                              dsolve::DENDROSOLVER_BH2_MAX_LEV);
        if (f2olmin < MAXDEAPTH_LEVEL_DIFF + 2) {
            if (!rank)
                std::cout << "BH min level should be larger than "
                          << (MAXDEAPTH_LEVEL_DIFF + 2) << std::endl;

            MPI_Abort(comm, 0);
        }
        function2Octree(f_init, dsolve::DENDROSOLVER_NUM_VARS, varIndex,
                        interpVars, tmpNodes,
                        (f2olmin - MAXDEAPTH_LEVEL_DIFF - 2),
                        dsolve::DENDROSOLVER_WAVELET_TOL,
                        dsolve::DENDROSOLVER_ELE_ORDER, comm);
    }

    if (!rank) std::cout << "Now generating mesh" << std::endl;

    /**
     * STEP 3
     *
     * Generate the mesh itself
     * Take the initial grid and set up the mesh to be used throughout
     * computations
     */
    ot::Mesh *mesh = ot::createMesh(
        tmpNodes.data(), tmpNodes.size(), dsolve::DENDROSOLVER_ELE_ORDER, comm,
        1, ot::SM_TYPE::FDM, dsolve::DENDROSOLVER_DENDRO_GRAIN_SZ,
        dsolve::DENDROSOLVER_LOAD_IMB_TOL, dsolve::DENDROSOLVER_SPLIT_FIX);

    if (!rank) {
        std::cout << "Mesh generation finished" << std::endl;
    }

    mesh->setDomainBounds(
        Point(dsolve::DENDROSOLVER_GRID_MIN_X, dsolve::DENDROSOLVER_GRID_MIN_Y,
              dsolve::DENDROSOLVER_GRID_MIN_Z),
        Point(dsolve::DENDROSOLVER_GRID_MAX_X, dsolve::DENDROSOLVER_GRID_MAX_Y,
              dsolve::DENDROSOLVER_GRID_MAX_Z));

    // io::vtk::mesh2vtuFine(mesh,"begin",0,NULL,NULL,0,NULL,NULL,0,NULL,NULL,false);

    if (!rank) {
        std::cout << "Domain bounds set" << std::endl;
    }

    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin, lmax);

    if (!rank) {
        std::cout << "Computed min max level:" << lmin << " " << lmax
                  << std::endl;
    }

    if (!rank) {
        std::cout << "================= Grid Info (Before init grid "
                     "converge):==============================================="
                     "========"
                  << std::endl;
        std::cout << "lmin: " << lmin << " lmax:" << lmax << std::endl;
        std::cout << "dx: "
                  << ((dsolve::DENDROSOLVER_COMPD_MAX[0] -
                       dsolve::DENDROSOLVER_COMPD_MIN[0]) *
                      ((1u << (m_uiMaxDepth - lmax)) /
                       ((double)dsolve::DENDROSOLVER_ELE_ORDER)) /
                      ((double)(1u << (m_uiMaxDepth))))
                  << std::endl;
        std::cout << "dt: "
                  << dsolve::DENDROSOLVER_CFL_FACTOR *
                         ((dsolve::DENDROSOLVER_COMPD_MAX[0] -
                           dsolve::DENDROSOLVER_COMPD_MIN[0]) *
                          ((1u << (m_uiMaxDepth - lmax)) /
                           ((double)dsolve::DENDROSOLVER_ELE_ORDER)) /
                          ((double)(1u << (m_uiMaxDepth))))
                  << std::endl;
        std::cout << "========================================================="
                     "======================================================"
                  << std::endl;
    }

    /**
     * STEP 4
     *
     * Prepare for RK solving
     * Includes prepping the the solver and potentially restoring from a
     * checkpoint if enabled
     */
    dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE =
        dsolve::DENDROSOLVER_CFL_FACTOR *
        ((dsolve::DENDROSOLVER_COMPD_MAX[0] -
          dsolve::DENDROSOLVER_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) /
          ((double)dsolve::DENDROSOLVER_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));

    ode::solver::RK_SOLVER rk_dsolve(mesh, dsolve::DENDROSOLVER_RK_TIME_BEGIN,
                                     dsolve::DENDROSOLVER_RK_TIME_END,
                                     dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE,
                                     (RKType)dsolve::DENDROSOLVER_RK_TYPE);

    if (dsolve::DENDROSOLVER_RESTORE_SOLVER == 1)
        rk_dsolve.restoreCheckPoint(
            dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX.c_str(), comm);

    if (!rank) std::cout << GRN << "Now starting solver!" << NRM << std::endl;

    /**-
     * STEP 5
     *
     * Start the solver!
     * This runs the main loop of the program and begins solving everything
     */
    if (!rank) {
        std::cout << GRN << "Now starting solver!" << NRM << std::endl;
    }
    dsolve::timer::t_rkSolve.start();
    rk_dsolve.rkSolve();
    dsolve::timer::t_rkSolve.stop();

    /**
     * STEP 6
     *
     * Solver finished, finalize everything and prepare for exiting the program
     */
    dsolve::timer::total_runtime.stop();

    double t2 = dsolve::timer::t_rkSolve.seconds;
    double t2_g;
    par::Mpi_Allreduce(&t2, &t2_g, 1, MPI_MAX, comm);
    if (!rank) {
        std::cout << "  SOLVER TIME (max): " << t2_g << std::endl;
    }

    rk_dsolve.freeMesh();

    MPI_Finalize();

    if (!rank) {
        std::cout << GRN << "Solver finished!" << NRM << std::endl;
    }

    return 0;
}
