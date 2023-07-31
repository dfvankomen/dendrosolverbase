/**
 * @file dsolveCtx.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Solver contex file.
 * @version 0.1
 * @date 2019-12-20
 *
 * School of Computing, University of Utah.
 * @copyright Copyright (c) 2019
 *
 */

#include "solverCtx.h"

#include <stdlib.h>

namespace dsolve {
SOLVERCtx::SOLVERCtx(ot::Mesh* pMesh) : Ctx() {
    m_uiMesh = pMesh;

    m_var[VL::CPU_EV].create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                                    ot::DVEC_LOC::HOST, DENDROSOLVER_NUM_VARS,
                                    true);
    m_var[VL::CPU_EV_UZ_IN].create_vector(
        m_uiMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        DENDROSOLVER_NUM_VARS, true);
    m_var[VL::CPU_EV_UZ_OUT].create_vector(
        m_uiMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        DENDROSOLVER_NUM_VARS, true);

    m_var[VL::CPU_CV].create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                                    ot::DVEC_LOC::HOST,
                                    DENDROSOLVER_CONSTRAINT_NUM_VARS, true);
    m_var[VL::CPU_CV_UZ_IN].create_vector(
        m_uiMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        DENDROSOLVER_CONSTRAINT_NUM_VARS, true);

    m_uiTinfo._m_uiStep = 0;
    m_uiTinfo._m_uiT = 0;
    m_uiTinfo._m_uiTb = dsolve::DENDROSOLVER_RK_TIME_BEGIN;
    m_uiTinfo._m_uiTe = dsolve::DENDROSOLVER_RK_TIME_END;
    m_uiTinfo._m_uiTh = dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE;

    m_uiElementOrder = dsolve::DENDROSOLVER_ELE_ORDER;

    m_uiMinPt =
        Point(dsolve::DENDROSOLVER_GRID_MIN_X, dsolve::DENDROSOLVER_GRID_MIN_Y,
              dsolve::DENDROSOLVER_GRID_MIN_Z);
    m_uiMaxPt =
        Point(dsolve::DENDROSOLVER_GRID_MAX_X, dsolve::DENDROSOLVER_GRID_MAX_Y,
              dsolve::DENDROSOLVER_GRID_MAX_Z);

    deallocate_deriv_workspace();
    allocate_deriv_workspace(m_uiMesh, 1);

    ot::dealloc_mpi_ctx<DendroScalar>(
        m_uiMesh, m_mpi_ctx, DENDROSOLVER_NUM_VARS, DENDROSOLVER_ASYNC_COMM_K);
    ot::alloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, DENDROSOLVER_NUM_VARS,
                                    DENDROSOLVER_ASYNC_COMM_K);

    return;
}

SOLVERCtx::~SOLVERCtx() {
    for (unsigned int i = 0; i < VL::END; i++) m_var[i].destroy_vector();

    deallocate_deriv_workspace();
    ot::dealloc_mpi_ctx<DendroScalar>(
        m_uiMesh, m_mpi_ctx, DENDROSOLVER_NUM_VARS, DENDROSOLVER_ASYNC_COMM_K);
}

int SOLVERCtx::initialize() {
    if (dsolve::DENDROSOLVER_RESTORE_SOLVER) {
        this->restore_checkpt();
        return 0;
    }

    this->init_grid();

    bool isRefine = false;
    DendroIntL oldElements, oldElements_g;
    DendroIntL newElements, newElements_g;

    DendroIntL oldGridPoints, oldGridPoints_g;
    DendroIntL newGridPoints, newGridPoints_g;

    unsigned int iterCount = 1;
    const unsigned int max_iter = dsolve::DENDROSOLVER_INIT_GRID_ITER;
    const unsigned int rank_global = m_uiMesh->getMPIRankGlobal();
    MPI_Comm gcomm = m_uiMesh->getMPIGlobalCommunicator();

    DendroScalar* unzipVar[dsolve::DENDROSOLVER_NUM_VARS];
    unsigned int refineVarIds[dsolve::DENDROSOLVER_NUM_REFINE_VARS];

    for (unsigned int vIndex = 0; vIndex < dsolve::DENDROSOLVER_NUM_REFINE_VARS;
         vIndex++)
        refineVarIds[vIndex] =
            dsolve::DENDROSOLVER_REFINE_VARIABLE_INDICES[vIndex];

    double wTol = dsolve::DENDROSOLVER_WAVELET_TOL;
    std::function<double(double, double, double, double* hx)> waveletTolFunc =
        [](double x, double y, double z, double* hx) {
            return dsolve::computeWTolDCoords(x, y, z, hx);
        };

    DVec& m_evar = m_var[VL::CPU_EV];
    DVec& m_evar_unz = m_var[VL::CPU_EV_UZ_IN];

    do {
        this->unzip(m_evar, m_evar_unz, dsolve::DENDROSOLVER_ASYNC_COMM_K);
        m_evar_unz.to_2d(unzipVar);
        // isRefine=this->is_remesh();
        // enforce WMAR refinement based refinement initially.
        // TODO: need to fix this for initial refinement!
        isRefine = dsolve::isReMeshWAMR(
            m_uiMesh, (const double**)unzipVar, refineVarIds,
            dsolve::DENDROSOLVER_NUM_REFINE_VARS, waveletTolFunc,
            dsolve::DENDROSOLVER_DENDRO_AMR_FAC);
        if (isRefine) {
            ot::Mesh* newMesh =
                this->remesh(dsolve::DENDROSOLVER_DENDRO_GRAIN_SZ,
                             dsolve::DENDROSOLVER_LOAD_IMB_TOL,
                             dsolve::DENDROSOLVER_SPLIT_FIX);

            oldElements = m_uiMesh->getNumLocalMeshElements();
            newElements = newMesh->getNumLocalMeshElements();

            oldGridPoints = m_uiMesh->getNumLocalMeshNodes();
            newGridPoints = newMesh->getNumLocalMeshNodes();

            par::Mpi_Allreduce(&oldElements, &oldElements_g, 1, MPI_SUM, gcomm);
            par::Mpi_Allreduce(&newElements, &newElements_g, 1, MPI_SUM, gcomm);

            par::Mpi_Allreduce(&oldGridPoints, &oldGridPoints_g, 1, MPI_SUM,
                               m_uiMesh->getMPIGlobalCommunicator());
            par::Mpi_Allreduce(&newGridPoints, &newGridPoints_g, 1, MPI_SUM,
                               m_uiMesh->getMPIGlobalCommunicator());

            if (!rank_global) {
                std::cout << "[dsolveCtx] iter : " << iterCount
                          << " (Remesh triggered) ->  old mesh : "
                          << oldElements_g << " new mesh : " << newElements_g
                          << std::endl;
                std::cout << "[dsolveCtx] iter : " << iterCount
                          << " (Remesh triggered) ->  old mesh (zip nodes) : "
                          << oldGridPoints_g
                          << " new mesh (zip nodes) : " << newGridPoints_g
                          << std::endl;
            }

            this->grid_transfer(newMesh);

            std::swap(m_uiMesh, newMesh);
            delete newMesh;

#ifdef __CUDACC__
            device::MeshGPU*& dptr_mesh = this->get_meshgpu_device_ptr();
            device::MeshGPU* mesh_gpu = this->get_meshgpu_host_handle();

            mesh_gpu->dealloc_mesh_on_device(dptr_mesh);
            dptr_mesh = mesh_gpu->alloc_mesh_on_device(m_uiMesh);
#endif
        }

        iterCount += 1;

    } while (isRefine &&
             (newElements_g != oldElements_g ||
              newGridPoints_g != oldGridPoints_g) &&
             (iterCount < max_iter));

    this->init_grid();

    // // realloc dsolve deriv space
    deallocate_deriv_workspace();
    allocate_deriv_workspace(m_uiMesh, 1);

    unsigned int lmin, lmax;
    m_uiMesh->computeMinMaxLevel(lmin, lmax);
    dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE =
        dsolve::DENDROSOLVER_CFL_FACTOR *
        ((dsolve::DENDROSOLVER_COMPD_MAX[0] -
          dsolve::DENDROSOLVER_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) /
          ((double)dsolve::DENDROSOLVER_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    m_uiTinfo._m_uiTh = dsolve::DENDROSOLVER_RK45_TIME_STEP_SIZE;

    if (!m_uiMesh->getMPIRankGlobal()) {
        const DendroScalar dx_finest =
            ((dsolve::DENDROSOLVER_COMPD_MAX[0] -
              dsolve::DENDROSOLVER_COMPD_MIN[0]) *
             ((1u << (m_uiMaxDepth - lmax)) /
              ((double)dsolve::DENDROSOLVER_ELE_ORDER)) /
             ((double)(1u << (m_uiMaxDepth))));
        const DendroScalar dt_finest =
            dsolve::DENDROSOLVER_CFL_FACTOR * dx_finest;

        std::cout << "================= Grid Info (After init grid "
                     "converge):==============================================="
                     "========"
                  << std::endl;
        std::cout << "lmin: " << lmin << " lmax:" << lmax << std::endl;
        std::cout << "dx: " << dx_finest << std::endl;
        std::cout << "dt: " << dt_finest << std::endl;
        std::cout << "========================================================="
                     "======================================================"
                  << std::endl;
    }

    return 0;
}

int SOLVERCtx::init_grid() {
    DVec& m_evar = m_var[VL::CPU_EV];
    DVec& m_dptr_evar = m_var[VL::GPU_EV];

    const ot::TreeNode* pNodes = &(*(m_uiMesh->getAllElements().begin()));
    const unsigned int eleOrder = m_uiMesh->getElementOrder();
    const unsigned int* e2n_cg = &(*(m_uiMesh->getE2NMapping().begin()));
    const unsigned int* e2n_dg = &(*(m_uiMesh->getE2NMapping_DG().begin()));
    const unsigned int nPe = m_uiMesh->getNumNodesPerElement();
    const unsigned int nodeLocalBegin = m_uiMesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd = m_uiMesh->getNodeLocalEnd();

    DendroScalar* zipIn[dsolve::DENDROSOLVER_NUM_VARS];
    m_evar.to_2d(zipIn);

    DendroScalar var1[dsolve::DENDROSOLVER_NUM_VARS];

    DendroScalar mp, mm, mp_adm, mm_adm, E, J1, J2, J3;
    // set the TP communicator.
    if (dsolve::DENDROSOLVER_ID_TYPE == 0) {
        TP_MPI_COMM = m_uiMesh->getMPIGlobalCommunicator();
        TwoPunctures((double)0, (double)0, (double)0, var1, &mp, &mm, &mp_adm,
                     &mm_adm, &E, &J1, &J2, &J3);
    }

    for (unsigned int elem = m_uiMesh->getElementLocalBegin();
         elem < m_uiMesh->getElementLocalEnd(); elem++) {
        DendroScalar var[dsolve::DENDROSOLVER_NUM_VARS];
        for (unsigned int k = 0; k < (eleOrder + 1); k++)
            for (unsigned int j = 0; j < (eleOrder + 1); j++)
                for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                    const unsigned int nodeLookUp_CG =
                        e2n_cg[elem * nPe +
                               k * (eleOrder + 1) * (eleOrder + 1) +
                               j * (eleOrder + 1) + i];
                    if (nodeLookUp_CG >= nodeLocalBegin &&
                        nodeLookUp_CG < nodeLocalEnd) {
                        const unsigned int nodeLookUp_DG =
                            e2n_dg[elem * nPe +
                                   k * (eleOrder + 1) * (eleOrder + 1) +
                                   j * (eleOrder + 1) + i];
                        unsigned int ownerID, ii_x, jj_y, kk_z;
                        m_uiMesh->dg2eijk(nodeLookUp_DG, ownerID, ii_x, jj_y,
                                          kk_z);
                        const DendroScalar len =
                            (double)(1u << (m_uiMaxDepth -
                                            pNodes[ownerID].getLevel()));

                        const DendroScalar x =
                            pNodes[ownerID].getX() + ii_x * (len / (eleOrder));
                        const DendroScalar y =
                            pNodes[ownerID].getY() + jj_y * (len / (eleOrder));
                        const DendroScalar z =
                            pNodes[ownerID].getZ() + kk_z * (len / (eleOrder));

                        if (dsolve::DENDROSOLVER_ID_TYPE == 0) {
                            const DendroScalar xx = GRIDX_TO_X(x);
                            const DendroScalar yy = GRIDY_TO_Y(y);
                            const DendroScalar zz = GRIDZ_TO_Z(z);

                            TwoPunctures((double)xx, (double)yy, (double)zz,
                                         var, &mp, &mm, &mp_adm, &mm_adm, &E,
                                         &J1, &J2, &J3);
                        } else {
                            // all other values are handled in the initial data
                            // wrapper including an error message
                            initialDataFunctionWrapper((double)x, (double)y,
                                                       (double)z, var);
                        }
                        for (unsigned int v = 0;
                             v < dsolve::DENDROSOLVER_NUM_VARS; v++)
                            zipIn[v][nodeLookUp_CG] = var[v];
                    }
                }
    }

    for (unsigned int node = m_uiMesh->getNodeLocalBegin();
         node < m_uiMesh->getNodeLocalEnd(); node++)
        enforce_dsolve_constraints(zipIn, node);

#ifdef DENDROSOLVER_EXTRACT_BH_LOCATIONS
    m_uiBHLoc[0] = Point(dsolve::BH1.getBHCoordX(), dsolve::BH1.getBHCoordY(),
                         dsolve::BH1.getBHCoordZ());
    m_uiBHLoc[1] = Point(dsolve::BH2.getBHCoordX(), dsolve::BH2.getBHCoordY(),
                         dsolve::BH2.getBHCoordZ());
#endif

    return 0;
}

int SOLVERCtx::finalize() { return 0; }

int SOLVERCtx::rhs(DVec* in, DVec* out, unsigned int sz, DendroScalar time) {
    // all the variables should be packed together.
    // assert(sz==1);
    // DendroScalar * sVar[DENDROSOLVER_NUM_VARS];
    // in->to_2d(sVar);

    this->unzip(*in, m_var[VL::CPU_EV_UZ_IN],
                dsolve::DENDROSOLVER_ASYNC_COMM_K);

#ifdef __PROFILE_CTX__
    this->m_uiCtxpt[ts::CTXPROFILE::RHS].start();
#endif

    DendroScalar* unzipIn[DENDROSOLVER_NUM_VARS];
    DendroScalar* unzipOut[DENDROSOLVER_NUM_VARS];

    m_var[CPU_EV_UZ_IN].to_2d(unzipIn);
    m_var[CPU_EV_UZ_OUT].to_2d(unzipOut);

    const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();

    dendroSolverRHS(unzipOut, (const DendroScalar**)unzipIn, blkList,
                    numBlocks);

#ifdef __PROFILE_CTX__
    this->m_uiCtxpt[ts::CTXPROFILE::RHS].stop();
#endif

    this->zip(m_var[CPU_EV_UZ_OUT], *out);

    return 0;
}

int SOLVERCtx::write_vtu() {
    if (!m_uiMesh->isActive()) return 0;

    DVec& m_evar = m_var[VL::CPU_EV];
    DVec& m_evar_unz = m_var[VL::CPU_EV_UZ_IN];
    DVec& m_cvar = m_var[VL::CPU_CV];
    DVec& m_cvar_unz = m_var[VL::CPU_CV_UZ_IN];

    this->unzip(m_evar, m_evar_unz, DENDROSOLVER_ASYNC_COMM_K);

    DendroScalar* consUnzipVar[dsolve::DENDROSOLVER_CONSTRAINT_NUM_VARS];
    DendroScalar* consVar[dsolve::DENDROSOLVER_CONSTRAINT_NUM_VARS];

    DendroScalar* evolUnzipVar[dsolve::DENDROSOLVER_NUM_VARS];
    DendroScalar* evolVar[dsolve::DENDROSOLVER_NUM_VARS];

    m_evar_unz.to_2d(evolUnzipVar);
    m_cvar_unz.to_2d(consUnzipVar);

    m_evar.to_2d(evolVar);
    m_cvar.to_2d(consVar);

#if DENDROSOLVER_COMPUTE_CONSTRAINTS

    const std::vector<ot::Block> blkList = m_uiMesh->getLocalBlockList();

    unsigned int offset;
    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    double dx, dy, dz;
    const Point pt_min(dsolve::DENDROSOLVER_COMPD_MIN[0],
                       dsolve::DENDROSOLVER_COMPD_MIN[1],
                       dsolve::DENDROSOLVER_COMPD_MIN[2]);
    const Point pt_max(dsolve::DENDROSOLVER_COMPD_MAX[0],
                       dsolve::DENDROSOLVER_COMPD_MAX[1],
                       dsolve::DENDROSOLVER_COMPD_MAX[2]);
    const unsigned int PW = dsolve::DENDROSOLVER_PADDING_WIDTH;

    for (unsigned int blk = 0; blk < blkList.size(); blk++) {
        offset = blkList[blk].getOffset();
        sz[0] = blkList[blk].getAllocationSzX();
        sz[1] = blkList[blk].getAllocationSzY();
        sz[2] = blkList[blk].getAllocationSzZ();

        bflag = blkList[blk].getBlkNodeFlag();

        dx = blkList[blk].computeDx(pt_min, pt_max);
        dy = blkList[blk].computeDy(pt_min, pt_max);
        dz = blkList[blk].computeDz(pt_min, pt_max);

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * dx;
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * dy;
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * dz;

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * dx;
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * dy;
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * dz;

        physical_constraints(consUnzipVar, (const DendroScalar**)evolUnzipVar,
                             offset, ptmin, ptmax, sz, bflag);
    }

    /*double consVecMin[dsolve::DENDROSOLVER_CONSTRAINT_NUM_VARS];
    double consVecMax[dsolve::DENDROSOLVER_CONSTRAINT_NUM_VARS];*/
    double constraintMaskedL2[dsolve::DENDROSOLVER_CONSTRAINT_NUM_VARS];
    this->zip(m_cvar_unz, m_cvar);
    m_uiMesh->readFromGhostBegin(m_cvar.get_vec_ptr(), m_cvar.get_dof());
    m_uiMesh->readFromGhostEnd(m_cvar.get_vec_ptr(), m_cvar.get_dof());

    dsolve::extractConstraints(m_uiMesh, (const DendroScalar**)consVar,
                               evolVar[BHLOC::EXTRACTION_VAR_ID],
                               BHLOC::EXTRACTION_TOL, m_uiTinfo._m_uiStep,
                               m_uiTinfo._m_uiT);
#ifndef DENDROSOLVER_KERR_SCHILD_TEST
#ifdef DENDROSOLVER_EXTRACT_GRAVITATIONAL_WAVES
    GW::extractFarFieldPsi4(m_uiMesh, (const DendroScalar**)consVar,
                            m_uiTinfo._m_uiStep, m_uiTinfo._m_uiT);
#endif
#endif

#endif

#ifdef DENDROSOLVER_ENABLE_VTU_OUTPUT

    if ((m_uiTinfo._m_uiStep % dsolve::DENDROSOLVER_IO_OUTPUT_FREQ) == 0) {
        std::vector<std::string> pDataNames;
        const unsigned int numConstVars =
            dsolve::DENDROSOLVER_NUM_CONST_VARS_VTU_OUTPUT;
        const unsigned int numEvolVars =
            dsolve::DENDROSOLVER_NUM_EVOL_VARS_VTU_OUTPUT;

        double* pData[(numConstVars + numEvolVars)];

        for (unsigned int i = 0; i < numEvolVars; i++) {
            pDataNames.push_back(
                std::string(dsolve::DENDROSOLVER_VAR_NAMES
                                [DENDROSOLVER_VTU_OUTPUT_EVOL_INDICES[i]]));
            pData[i] = evolVar[DENDROSOLVER_VTU_OUTPUT_EVOL_INDICES[i]];
        }

        for (unsigned int i = 0; i < numConstVars; i++) {
            pDataNames.push_back(
                std::string(dsolve::DENDROSOLVER_CONSTRAINT_VAR_NAMES
                                [DENDROSOLVER_VTU_OUTPUT_CONST_INDICES[i]]));
            pData[numEvolVars + i] =
                consVar[DENDROSOLVER_VTU_OUTPUT_CONST_INDICES[i]];
        }

        std::vector<char*> pDataNames_char;
        pDataNames_char.reserve(pDataNames.size());

        for (unsigned int i = 0; i < pDataNames.size(); i++)
            pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

        const char* fDataNames[] = {"Time", "Cycle"};
        const double fData[] = {m_uiTinfo._m_uiT, (double)m_uiTinfo._m_uiStep};

        char fPrefix[256];
        sprintf(fPrefix, "%s_%d", dsolve::DENDROSOLVER_VTU_FILE_PREFIX.c_str(),
                m_uiTinfo._m_uiStep);

        if (dsolve::DENDROSOLVER_VTU_Z_SLICE_ONLY) {
            unsigned int s_val[3] = {1u << (m_uiMaxDepth - 1),
                                     1u << (m_uiMaxDepth - 1),
                                     1u << (m_uiMaxDepth - 1)};
            unsigned int s_norm[3] = {0, 0, 1};
            io::vtk::mesh2vtu_slice(
                m_uiMesh, s_val, s_norm, fPrefix, 2, fDataNames, fData,
                (numEvolVars + numConstVars), (const char**)&pDataNames_char[0],
                (const double**)pData);
        } else
            io::vtk::mesh2vtuFine(m_uiMesh, fPrefix, 2, fDataNames, fData,
                                  (numEvolVars + numConstVars),
                                  (const char**)&pDataNames_char[0],
                                  (const double**)pData);
    }

#endif

#ifdef DENDROSOLVER_EXTRACT_BH_LOCATIONS
    dsolve::writeBHCoordinates((const ot::Mesh*)m_uiMesh,
                               (const Point*)m_uiBHLoc, 2, m_uiTinfo._m_uiStep,
                               m_uiTinfo._m_uiT);
#endif

    return 0;
}

int SOLVERCtx::write_checkpt() {
    if (!m_uiMesh->isActive()) return 0;

    unsigned int cpIndex;
    (m_uiTinfo._m_uiStep % (2 * dsolve::DENDROSOLVER_CHECKPT_FREQ) == 0)
        ? cpIndex = 0
        : cpIndex = 1;  // to support alternate file writing.

    const bool is_merged =
        ((dsolve::DENDROSOLVER_BH_LOC[0] - dsolve::DENDROSOLVER_BH_LOC[1])
             .abs() < 0.1);
    if (is_merged && !dsolve::DENDROSOLVER_MERGED_CHKPT_WRITTEN) {
        cpIndex = 3;
        dsolve::DENDROSOLVER_MERGED_CHKPT_WRITTEN = true;
    }

    unsigned int rank = m_uiMesh->getMPIRank();
    unsigned int npes = m_uiMesh->getMPICommSize();

    DendroScalar* eVar[DENDROSOLVER_NUM_VARS];
    DVec& m_evar = m_var[VL::CPU_EV];
    m_evar.to_2d(eVar);

    char fName[256];
    const ot::TreeNode* pNodes = &(*(m_uiMesh->getAllElements().begin() +
                                     m_uiMesh->getElementLocalBegin()));
    sprintf(fName, "%s_octree_%d_%d.oct",
            dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX.c_str(), cpIndex, rank);
    io::checkpoint::writeOctToFile(fName, pNodes,
                                   m_uiMesh->getNumLocalMeshElements());

    unsigned int numVars = dsolve::DENDROSOLVER_NUM_VARS;
    const char** varNames = dsolve::DENDROSOLVER_VAR_NAMES;

    sprintf(fName, "%s_%d_%d.var",
            dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX.c_str(), cpIndex, rank);
    io::checkpoint::writeVecToFile(fName, m_uiMesh, (const double**)eVar,
                                   dsolve::DENDROSOLVER_NUM_VARS);

    if (!rank) {
        sprintf(fName, "%s_step_%d.cp",
                dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX.c_str(), cpIndex);
        std::cout << "[DENDROSOLVERCtx] \t writing checkpoint file : " << fName
                  << std::endl;
        std::ofstream outfile(fName);
        if (!outfile) {
            std::cout << fName << " file open failed " << std::endl;
            return 0;
        }

        json checkPoint;
        checkPoint["DENDRO_TS_TIME_BEGIN"] = m_uiTinfo._m_uiTb;
        checkPoint["DENDRO_TS_TIME_END"] = m_uiTinfo._m_uiTe;
        checkPoint["DENDRO_TS_ELEMENT_ORDER"] = m_uiElementOrder;

        checkPoint["DENDRO_TS_TIME_CURRENT"] = m_uiTinfo._m_uiT;
        checkPoint["DENDRO_TS_STEP_CURRENT"] = m_uiTinfo._m_uiStep;
        checkPoint["DENDRO_TS_TIME_STEP_SIZE"] = m_uiTinfo._m_uiTh;
        checkPoint["DENDRO_TS_LAST_IO_TIME"] = m_uiTinfo._m_uiT;

        checkPoint["DENDRO_TS_WAVELET_TOLERANCE"] =
            dsolve::DENDROSOLVER_WAVELET_TOL;
        checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"] =
            dsolve::DENDROSOLVER_LOAD_IMB_TOL;
        checkPoint["DENDRO_TS_NUM_VARS"] =
            numVars;  // number of variables to restore.
        checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"] =
            m_uiMesh->getMPICommSize();  // (note that rank 0 is always active).

        checkPoint["DENDRO_BH1_X"] = m_uiBHLoc[0].x();
        checkPoint["DENDRO_BH1_Y"] = m_uiBHLoc[0].y();
        checkPoint["DENDRO_BH1_Z"] = m_uiBHLoc[0].z();

        checkPoint["DENDRO_BH2_X"] = m_uiBHLoc[1].x();
        checkPoint["DENDRO_BH2_Y"] = m_uiBHLoc[1].y();
        checkPoint["DENDRO_BH2_Z"] = m_uiBHLoc[1].z();

        outfile << std::setw(4) << checkPoint << std::endl;
        outfile.close();
    }

    return 0;
}

int SOLVERCtx::restore_checkpt() {
    unsigned int numVars = 0;
    std::vector<ot::TreeNode> octree;
    json checkPoint;

    int rank;
    int npes;
    MPI_Comm comm = m_uiMesh->getMPIGlobalCommunicator();
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    unsigned int activeCommSz;

    char fName[256];
    unsigned int restoreStatus = 0;
    unsigned int restoreStatusGlobal =
        0;  // 0 indicates successfully restorable.

    ot::Mesh* newMesh;
    unsigned int restoreStep[2];
    restoreStep[0] = 0;
    restoreStep[1] = 0;

    unsigned int restoreFileIndex = 0;

    for (unsigned int cpIndex = 0; cpIndex < 2; cpIndex++) {
        restoreStatus = 0;

        if (!rank) {
            sprintf(fName, "%s_step_%d.cp",
                    dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX.c_str(), cpIndex);
            std::ifstream infile(fName);
            if (!infile) {
                std::cout << fName << " file open failed " << std::endl;
                restoreStatus = 1;
            }

            if (restoreStatus == 0) {
                infile >> checkPoint;
                m_uiTinfo._m_uiTb = checkPoint["DENDRO_TS_TIME_BEGIN"];
                m_uiTinfo._m_uiTe = checkPoint["DENDRO_TS_TIME_END"];
                m_uiTinfo._m_uiT = checkPoint["DENDRO_TS_TIME_CURRENT"];
                m_uiTinfo._m_uiStep = checkPoint["DENDRO_TS_STEP_CURRENT"];
                m_uiTinfo._m_uiTh = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
                m_uiElementOrder = checkPoint["DENDRO_TS_ELEMENT_ORDER"];

                dsolve::DENDROSOLVER_WAVELET_TOL =
                    checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                dsolve::DENDROSOLVER_LOAD_IMB_TOL =
                    checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];

                numVars = checkPoint["DENDRO_TS_NUM_VARS"];
                activeCommSz = checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];

                m_uiBHLoc[0] = Point((double)checkPoint["DENDRO_BH1_X"],
                                     (double)checkPoint["DENDRO_BH1_Y"],
                                     (double)checkPoint["DENDRO_BH1_Z"]);
                m_uiBHLoc[1] = Point((double)checkPoint["DENDRO_BH2_X"],
                                     (double)checkPoint["DENDRO_BH2_Y"],
                                     (double)checkPoint["DENDRO_BH2_Z"]);
                restoreStep[cpIndex] = m_uiTinfo._m_uiStep;
            }
        }
    }

    if (!rank) {
        if (restoreStep[0] < restoreStep[1])
            restoreFileIndex = 1;
        else
            restoreFileIndex = 0;
    }

    par::Mpi_Bcast(&restoreFileIndex, 1, 0, comm);

    restoreStatus = 0;
    octree.clear();
    if (!rank)
        std::cout
            << "[DENDROSOLVERCtx] :  Trying to restore from checkpoint index : "
            << restoreFileIndex << std::endl;

    if (!rank) {
        sprintf(fName, "%s_step_%d.cp",
                dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX.c_str(),
                restoreFileIndex);
        std::ifstream infile(fName);
        if (!infile) {
            std::cout << fName << " file open failed " << std::endl;
            restoreStatus = 1;
        }

        if (restoreStatus == 0) {
            infile >> checkPoint;
            m_uiTinfo._m_uiTb = checkPoint["DENDRO_TS_TIME_BEGIN"];
            m_uiTinfo._m_uiTe = checkPoint["DENDRO_TS_TIME_END"];
            m_uiTinfo._m_uiT = checkPoint["DENDRO_TS_TIME_CURRENT"];
            m_uiTinfo._m_uiStep = checkPoint["DENDRO_TS_STEP_CURRENT"];
            m_uiTinfo._m_uiTh = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
            m_uiElementOrder = checkPoint["DENDRO_TS_ELEMENT_ORDER"];

            dsolve::DENDROSOLVER_WAVELET_TOL =
                checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
            dsolve::DENDROSOLVER_LOAD_IMB_TOL =
                checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];

            numVars = checkPoint["DENDRO_TS_NUM_VARS"];
            activeCommSz = checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];

            m_uiBHLoc[0] = Point((double)checkPoint["DENDRO_BH1_X"],
                                 (double)checkPoint["DENDRO_BH1_Y"],
                                 (double)checkPoint["DENDRO_BH1_Z"]);
            m_uiBHLoc[1] = Point((double)checkPoint["DENDRO_BH2_X"],
                                 (double)checkPoint["DENDRO_BH2_Y"],
                                 (double)checkPoint["DENDRO_BH2_Z"]);
            restoreStep[restoreFileIndex] = m_uiTinfo._m_uiStep;
        }
    }

    par::Mpi_Allreduce(&restoreStatus, &restoreStatusGlobal, 1, MPI_MAX, comm);
    if (restoreStatusGlobal == 1) {
        if (!rank)
            std::cout << "[DENDROSOLVERCtx] : Restore step failed, restore "
                         "file corrupted. "
                      << std::endl;
        MPI_Abort(comm, 0);
    }

    MPI_Bcast(&m_uiTinfo, sizeof(ts::TSInfo), MPI_BYTE, 0, comm);
    par::Mpi_Bcast(&dsolve::DENDROSOLVER_WAVELET_TOL, 1, 0, comm);
    par::Mpi_Bcast(&dsolve::DENDROSOLVER_LOAD_IMB_TOL, 1, 0, comm);

    par::Mpi_Bcast(&numVars, 1, 0, comm);
    par::Mpi_Bcast(&m_uiElementOrder, 1, 0, comm);
    par::Mpi_Bcast(&activeCommSz, 1, 0, comm);

    par::Mpi_Bcast(m_uiBHLoc, 2, 0, comm);
    dsolve::DENDROSOLVER_BH_LOC[0] = m_uiBHLoc[0];
    dsolve::DENDROSOLVER_BH_LOC[1] = m_uiBHLoc[1];

    if (activeCommSz > npes) {
        if (!rank)
            std::cout
                << " [DENDROSOLVERCtx] : checkpoint file written from  a "
                   "larger "
                   "communicator than the current global comm. (i.e. "
                   "communicator shrinking not allowed in the restore step. )"
                << std::endl;

        MPI_Abort(comm, 0);
    }

    bool isActive = (rank < activeCommSz);

    MPI_Comm newComm;
    par::splitComm2way(isActive, &newComm, comm);

    if (isActive) {
        int activeRank;
        int activeNpes;

        MPI_Comm_rank(newComm, &activeRank);
        MPI_Comm_size(newComm, &activeNpes);
        assert(activeNpes == activeCommSz);

        sprintf(fName, "%s_octree_%d_%d.oct",
                dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX.c_str(),
                restoreFileIndex, activeRank);
        restoreStatus = io::checkpoint::readOctFromFile(fName, octree);
        assert(par::test::isUniqueAndSorted(octree, newComm));
    }

    par::Mpi_Allreduce(&restoreStatus, &restoreStatusGlobal, 1, MPI_MAX, comm);
    if (restoreStatusGlobal == 1) {
        if (!rank)
            std::cout << "[DENDROSOLVERCtx]: octree (*.oct) restore file is "
                         "corrupted "
                      << std::endl;
        MPI_Abort(comm, 0);
    }

    newMesh = new ot::Mesh(octree, 1, m_uiElementOrder, activeCommSz, comm);
    newMesh->setDomainBounds(
        Point(dsolve::DENDROSOLVER_GRID_MIN_X, dsolve::DENDROSOLVER_GRID_MIN_Y,
              dsolve::DENDROSOLVER_GRID_MIN_Z),
        Point(dsolve::DENDROSOLVER_GRID_MAX_X, dsolve::DENDROSOLVER_GRID_MAX_Y,
              dsolve::DENDROSOLVER_GRID_MAX_Z));
    // no need to transfer data only to resize the contex variables.
    // this->grid_transfer(newMesh);
    for (unsigned int i = 0; i < VL::END; i++) m_var[i].destroy_vector();

    m_var[VL::CPU_EV].create_vector(newMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                                    ot::DVEC_LOC::HOST, DENDROSOLVER_NUM_VARS,
                                    true);
    m_var[VL::CPU_EV_UZ_IN].create_vector(
        newMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        DENDROSOLVER_NUM_VARS, true);
    m_var[VL::CPU_EV_UZ_OUT].create_vector(
        newMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        DENDROSOLVER_NUM_VARS, true);

    m_var[VL::CPU_CV].create_vector(newMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                                    ot::DVEC_LOC::HOST,
                                    DENDROSOLVER_CONSTRAINT_NUM_VARS, true);
    m_var[VL::CPU_CV_UZ_IN].create_vector(
        newMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        DENDROSOLVER_CONSTRAINT_NUM_VARS, true);

    ot::dealloc_mpi_ctx<DendroScalar>(
        m_uiMesh, m_mpi_ctx, DENDROSOLVER_NUM_VARS, DENDROSOLVER_ASYNC_COMM_K);
    ot::alloc_mpi_ctx<DendroScalar>(newMesh, m_mpi_ctx, DENDROSOLVER_NUM_VARS,
                                    DENDROSOLVER_ASYNC_COMM_K);

    // only reads the evolution variables.
    if (isActive) {
        int activeRank;
        int activeNpes;

        DendroScalar* inVec[DENDROSOLVER_NUM_VARS];
        DVec& m_evar = m_var[VL::CPU_EV];
        m_evar.to_2d(inVec);

        MPI_Comm_rank(newComm, &activeRank);
        MPI_Comm_size(newComm, &activeNpes);
        assert(activeNpes == activeCommSz);

        sprintf(fName, "%s_%d_%d.var",
                dsolve::DENDROSOLVER_CHKPT_FILE_PREFIX.c_str(),
                restoreFileIndex, activeRank);
        restoreStatus = io::checkpoint::readVecFromFile(
            fName, newMesh, inVec, dsolve::DENDROSOLVER_NUM_VARS);
    }

    MPI_Comm_free(&newComm);
    par::Mpi_Allreduce(&restoreStatus, &restoreStatusGlobal, 1, MPI_MAX, comm);
    if (restoreStatusGlobal == 1) {
        if (!rank)
            std::cout
                << "[DENDROSOLVERCtx]: varible (*.var) restore file currupted "
                << std::endl;
        MPI_Abort(comm, 0);
    }

    std::swap(m_uiMesh, newMesh);
    delete newMesh;

    // with the mesh restored, we can now allocate the derivative workspace
    deallocate_deriv_workspace();
    allocate_deriv_workspace(m_uiMesh, 1);

    unsigned int localSz = m_uiMesh->getNumLocalMeshElements();
    unsigned int totalElems = 0;
    par::Mpi_Allreduce(&localSz, &totalElems, 1, MPI_SUM, comm);

    if (!rank)
        std::cout << " checkpoint at step : " << m_uiTinfo._m_uiStep
                  << "active Comm. sz: " << activeCommSz
                  << " restore successful: "
                  << " restored mesh size: " << totalElems << std::endl;

    m_uiIsETSSynced = false;
    return 0;
}

int SOLVERCtx::pre_timestep(DVec sIn) { return 0; }

int SOLVERCtx::pre_stage(DVec sIn) { return 0; }

int SOLVERCtx::post_stage(DVec sIn) { return 0; }

int SOLVERCtx::post_timestep(DVec sIn) {
    DendroScalar* evar[DENDROSOLVER_NUM_VARS];
    sIn.to_2d(evar);
    for (unsigned int node = m_uiMesh->getNodeLocalBegin();
         node < m_uiMesh->getNodeLocalEnd(); node++)
        enforce_solver_constraints(evar, node);

    return 0;
}

bool SOLVERCtx::is_remesh() {
    bool isRefine = false;
    if (dsolve::DENDROSOLVER_ENABLE_BLOCK_ADAPTIVITY) return false;

    if (dsolve::DENDROSOLVER_REFINEMENT_MODE ==
        dsolve::RefinementMode::REFINE_MODE_NONE)
        return false;

    MPI_Comm comm = m_uiMesh->getMPIGlobalCommunicator();

    DVec& m_evar = m_var[VL::CPU_EV];
    DVec& m_evar_unz = m_var[VL::CPU_EV_UZ_IN];

    this->unzip(m_evar, m_evar_unz, dsolve::DENDROSOLVER_ASYNC_COMM_K);

    DendroScalar* unzipVar[DENDROSOLVER_NUM_VARS];
    m_evar_unz.to_2d(unzipVar);

    unsigned int refineVarIds[dsolve::DENDROSOLVER_NUM_REFINE_VARS];
    for (unsigned int vIndex = 0; vIndex < dsolve::DENDROSOLVER_NUM_REFINE_VARS;
         vIndex++)
        refineVarIds[vIndex] =
            dsolve::DENDROSOLVER_REFINE_VARIABLE_INDICES[vIndex];

    double wTol = dsolve::DENDROSOLVER_WAVELET_TOL;
    std::function<double(double, double, double, double* hx)> waveletTolFunc =
        [](double x, double y, double z, double* hx) {
            return dsolve::computeWTolDCoords(x, y, z, hx);
        };

    // TODO: can remove the different refinement modes if we don't want BH tied
    // to it
    if (dsolve::DENDROSOLVER_REFINEMENT_MODE == dsolve::RefinementMode::WAMR) {
        isRefine = dsolve::isReMeshWAMR(
            m_uiMesh, (const double**)unzipVar, refineVarIds,
            dsolve::DENDROSOLVER_NUM_REFINE_VARS, waveletTolFunc,
            dsolve::DENDROSOLVER_DENDRO_AMR_FAC);

    } else if (dsolve::DENDROSOLVER_REFINEMENT_MODE ==
               dsolve::RefinementMode::EH) {
        isRefine = dsolve::isRemeshEH(
            m_uiMesh, (const double**)unzipVar, dsolve::VAR::U_ALPHA,
            dsolve::DENDROSOLVER_EH_REFINE_VAL,
            dsolve::DENDROSOLVER_EH_COARSEN_VAL, true);

    } else if (dsolve::DENDROSOLVER_REFINEMENT_MODE ==
               dsolve::RefinementMode::EH_WAMR) {
        const bool isR1 = dsolve::isReMeshWAMR(
            m_uiMesh, (const double**)unzipVar, refineVarIds,
            dsolve::DENDROSOLVER_NUM_REFINE_VARS, waveletTolFunc,
            dsolve::DENDROSOLVER_DENDRO_AMR_FAC);
        const bool isR2 = dsolve::isRemeshEH(
            m_uiMesh, (const double**)unzipVar, dsolve::VAR::U_ALPHA,
            dsolve::DENDROSOLVER_EH_REFINE_VAL,
            dsolve::DENDROSOLVER_EH_COARSEN_VAL, false);

        isRefine = (isR1 || isR2);
    } else if (dsolve::DENDROSOLVER_REFINEMENT_MODE ==
               dsolve::RefinementMode::BH_LOC) {
        isRefine = dsolve::isRemeshBH(m_uiMesh, m_uiBHLoc);
    }

    return isRefine;
}

int SOLVERCtx::update_app_vars() {
    m_uiEVar = m_uiEvolutionVar[0];
    m_uiCVar = m_uiConstrainedVar[0];

    m_uiEUnzip[0] = m_uiEvolutionUnzipVar[0];
    m_uiEUnzip[1] = m_uiEvolutionUnzipVar[1];

    m_uiCUnzip[0] = m_uiConstraintUnzipVar[0];

    return 0;
}

DVec& SOLVERCtx::get_evolution_vars() { return m_var[CPU_EV]; }

DVec& SOLVERCtx::get_constraint_vars() { return m_var[CPU_CV]; }

// DVec SOLVERCtx::get_primitive_vars() { return m_var[CPU_PV]; }

int SOLVERCtx::terminal_output() {
    if (m_uiMesh->isActive()) {
        DendroScalar min = 0, max = 0;
        DVec& m_evar = m_var[VL::CPU_EV];
        min = vecMin(m_uiMesh, m_evar.get_vec_ptr(), ot::VEC_TYPE::CG_NODAL,
                     true);
        max = vecMax(m_uiMesh, m_evar.get_vec_ptr(), ot::VEC_TYPE::CG_NODAL,
                     true);

        // TODO: maybe adjust this one?
        if (!(m_uiMesh->getMPIRank())) {
            std::cout << "[DENDROSOLVERCtx]:  "
                      << dsolve::DENDROSOLVER_VAR_NAMES[0]
                      << " (min,max) : \t ( " << min << ", " << max << " ) "
                      << std::endl;
            if (std::isnan(min) || std::isnan(max)) {
                std::cout << "[Error]: NAN detected " << std::endl;
                MPI_Abort(m_uiMesh->getMPICommunicator(), 0);
            }
        }
    }

    return 0;
}

int SOLVERCtx::grid_transfer(const ot::Mesh* m_new) {
#ifdef __PROFILE_CTX__
    m_uiCtxpt[ts::CTXPROFILE::GRID_TRASFER].start();
#endif
    DVec& m_evar = m_var[VL::CPU_EV];
    DVec::grid_transfer(m_uiMesh, m_new, m_evar);
    // printf("igt ended\n");

    m_var[VL::CPU_CV].destroy_vector();
    m_var[VL::CPU_CV_UZ_IN].destroy_vector();

    m_var[VL::CPU_EV_UZ_IN].destroy_vector();
    m_var[VL::CPU_EV_UZ_OUT].destroy_vector();

    m_var[VL::CPU_CV].create_vector(m_new, ot::DVEC_TYPE::OCT_SHARED_NODES,
                                    ot::DVEC_LOC::HOST,
                                    DENDROSOLVER_CONSTRAINT_NUM_VARS, true);
    m_var[VL::CPU_CV_UZ_IN].create_vector(
        m_new, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        DENDROSOLVER_CONSTRAINT_NUM_VARS, true);

    m_var[VL::CPU_EV_UZ_IN].create_vector(
        m_new, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        DENDROSOLVER_NUM_VARS, true);
    m_var[VL::CPU_EV_UZ_OUT].create_vector(
        m_new, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        DENDROSOLVER_NUM_VARS, true);

    ot::dealloc_mpi_ctx<DendroScalar>(
        m_uiMesh, m_mpi_ctx, DENDROSOLVER_NUM_VARS, DENDROSOLVER_ASYNC_COMM_K);
    ot::alloc_mpi_ctx<DendroScalar>(m_new, m_mpi_ctx, DENDROSOLVER_NUM_VARS,
                                    DENDROSOLVER_ASYNC_COMM_K);

    m_uiIsETSSynced = false;

#ifdef __PROFILE_CTX__
    m_uiCtxpt[ts::CTXPROFILE::GRID_TRASFER].stop();
#endif
    return 0;
}

unsigned int SOLVERCtx::compute_lts_ts_offset() {
    const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
    const ot::TreeNode* pNodes = m_uiMesh->getAllElements().data();

    const unsigned int ldiff = 0;

    unsigned int lmin, lmax;
    m_uiMesh->computeMinMaxLevel(lmin, lmax);
    DENDROSOLVER_LTS_TS_OFFSET = 0;

    const double dt_min = dsolve::DENDROSOLVER_CFL_FACTOR *
                          ((dsolve::DENDROSOLVER_COMPD_MAX[0] -
                            dsolve::DENDROSOLVER_COMPD_MIN[0]) *
                           ((1u << (m_uiMaxDepth - lmax)) /
                            ((double)dsolve::DENDROSOLVER_ELE_ORDER)) /
                           ((double)(1u << (m_uiMaxDepth))));
    const double dt_eta_fac = dt_min * ETA_CONST;

    for (unsigned int blk = 0; blk < numBlocks; blk++) {
        const unsigned int blev =
            pNodes[blkList[blk].getLocalElementBegin()].getLevel();
        ot::TreeNode blkNode = blkList[blk].getBlockNode();
        const unsigned int szby2 =
            (1u << (m_uiMaxDepth - blkNode.getLevel() - 1));

        Point pt_oct(blkNode.minX() + szby2, blkNode.minY() + szby2,
                     blkNode.minZ() + szby2);
        Point pt_domain;
        m_uiMesh->octCoordToDomainCoord(pt_oct, pt_domain);
        const double r_coord = pt_domain.abs();

        double eta;
        if (dsolve::RIT_ETA_FUNCTION == 0) {
            // HAD eta function
            eta = ETA_CONST;
            if (r_coord >= ETA_R0) {
                eta *= pow((ETA_R0 / r_coord), ETA_DAMPING_EXP);
            }
        } else {
            // RIT eta function
            double w = r_coord / dsolve::RIT_ETA_WIDTH;
            double arg = -w * w * w * w;
            eta = (dsolve::RIT_ETA_CENTRAL - dsolve::RIT_ETA_OUTER) * exp(arg) +
                  dsolve::RIT_ETA_OUTER;
        }

        const double dt_eta = dt_eta_fac * (1 / eta);
        const double dt_cfl = (1u << (lmax - blev)) * dt_min;

        const double dt_feasible = std::min(dt_eta, dt_cfl);

        if (dt_feasible > dt_min) {
            unsigned int lts_offset =
                lmax - blev - std::floor(std::log2(dt_feasible / dt_min));
            if (DENDROSOLVER_LTS_TS_OFFSET < lts_offset)
                DENDROSOLVER_LTS_TS_OFFSET = lts_offset;
        }
    }

    unsigned int lts_offset_max = 0;
    par::Mpi_Allreduce(&DENDROSOLVER_LTS_TS_OFFSET, &lts_offset_max, 1, MPI_MAX,
                       m_uiMesh->getMPIGlobalCommunicator());
    DENDROSOLVER_LTS_TS_OFFSET = lts_offset_max;

    if (m_uiMesh->isActive() && (!(m_uiMesh->getMPIRank())))
        std::cout << "LTS offset : " << DENDROSOLVER_LTS_TS_OFFSET << std::endl;

    return DENDROSOLVER_LTS_TS_OFFSET;
}

unsigned int SOLVERCtx::getBlkTimestepFac(unsigned int blev, unsigned int lmin,
                                          unsigned int lmax) {
    const unsigned int ldiff = DENDROSOLVER_LTS_TS_OFFSET;
    if ((lmax - blev) <= ldiff)
        return 1;
    else {
        return 1u << (lmax - blev - ldiff);
    }
}

void SOLVERCtx::evolve_bh_loc(DVec sIn, double dt) {
#ifdef DENDROSOLVER_EXTRACT_BH_LOCATIONS

    m_uiMesh->readFromGhostBegin(
        sIn.GetVecArray() + VAR::U_BETA0 * m_uiMesh->getDegOfFreedom(), 3);
    m_uiMesh->readFromGhostEnd(
        sIn.GetVecArray() + VAR::U_BETA0 * m_uiMesh->getDegOfFreedom(), 3);
    Point bhLoc[2];
    DendroScalar** evar = new DendroScalar*[dsolve::DENDROSOLVER_NUM_VARS];
    sIn.Get2DArray(evar, true);
    dsolve::computeBHLocations((const ot::Mesh*)m_uiMesh, m_uiBHLoc, bhLoc,
                               evar, dt);
    // if(!m_uiMesh->getMPIRankGlobal())
    // {
    //     std::cout<<"bh0 "<<bhLoc[0]<<std::endl;
    //     std::cout<<"bh1 "<<bhLoc[1]<<std::endl;

    // }
    m_uiBHLoc[0] = bhLoc[0];
    m_uiBHLoc[1] = bhLoc[1];

    dsolve::DENDROSOLVER_BH_LOC[0] = m_uiBHLoc[0];
    dsolve::DENDROSOLVER_BH_LOC[1] = m_uiBHLoc[1];

    delete[] evar;

// old bh location extractor.
#if 0
                DendroScalar ** evar = new DendroScalar*[dsolve::DENDROSOLVER_NUM_VARS];
                Point bhLoc[2];
                sIn.Get2DArray(evar,true);
                dsolve::extractBHCoords((const ot::Mesh *)m_uiMesh,(const DendroScalar*)evar[BHLOC::EXTRACTION_VAR_ID],BHLOC::EXTRACTION_TOL,(const Point *) m_uiBHLoc,2,(Point*)bhLoc);
                
                m_uiBHLoc[0] = bhLoc[0];
                m_uiBHLoc[1] = bhLoc[1];

                delete [] evar;
#endif

#endif

    return;
}

void SOLVERCtx::lts_smooth(DVec sIn, LTS_SMOOTH_MODE mode) {
    // NOTE: this function seems entirely unused
    // TODO: remove it or update with it being removed
    // the last version of this function (from what I could tell) was
    // to apply smoothing across the different variables for like KO
    // dissipation
}

}  // end of namespace dsolve.