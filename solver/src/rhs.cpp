#include "rhs.h"

#include "hadrhs.h"
#include "solver_main.h"

using namespace std;
using namespace dsolve;

void dendroSolverRHS(double **uzipVarsRHS, const double **uZipVars,
                     const ot::Block *blkList, unsigned int numBlocks) {
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

#ifdef DENDROSOLVER_ENABLE_CUDA

    // TODO: generate CUDA-based stuff
    cuda::SOLVERComputeParams dsolveParams;
    dsolveParams.DENDROSOLVER_LAMBDA[0] = dsolve::DENDROSOLVER_LAMBDA[0];
    dsolveParams.DENDROSOLVER_LAMBDA[1] = dsolve::DENDROSOLVER_LAMBDA[1];
    dsolveParams.DENDROSOLVER_LAMBDA[2] = dsolve::DENDROSOLVER_LAMBDA[2];
    dsolveParams.DENDROSOLVER_LAMBDA[3] = dsolve::DENDROSOLVER_LAMBDA[3];

    dsolveParams.DENDROSOLVER_LAMBDA_F[0] = dsolve::DENDROSOLVER_LAMBDA_F[0];
    dsolveParams.DENDROSOLVER_LAMBDA_F[1] = dsolve::DENDROSOLVER_LAMBDA_F[1];

    dsolveParams.DENDROSOLVER_ETA_POWER[0] = dsolve::DENDROSOLVER_ETA_POWER[0];
    dsolveParams.DENDROSOLVER_ETA_POWER[1] = dsolve::DENDROSOLVER_ETA_POWER[1];

    dsolveParams.ETA_R0 = dsolve::ETA_R0;
    dsolveParams.ETA_CONST = dsolve::ETA_CONST;
    dsolveParams.ETA_DAMPING = dsolve::ETA_DAMPING;
    dsolveParams.ETA_DAMPING_EXP = dsolve::ETA_DAMPING_EXP;
    dsolveParams.KO_DISS_SIGMA = dsolve::KO_DISS_SIGMA;

    dim3 threadBlock(16, 16, 1);
    cuda::computeRHS(uzipVarsRHS, (const double **)uZipVars, blkList, numBlocks,
                     (const cuda::SOLVERComputeParams *)&dsolveParams,
                     threadBlock, pt_min, pt_max, 1);
#else

    for (unsigned int blk = 0; blk < numBlocks; blk++) {
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

#ifdef DENDROSOLVER_RHS_STAGED_COMP
        dendroSolverRHSUnpacked_sep(uzipVarsRHS, (const double **)uZipVars,
                                    offset, ptmin, ptmax, sz, bflag);
#else
        dendroSolverRHSUnpacked(uzipVarsRHS, (const double **)uZipVars, offset,
                                ptmin, ptmax, sz, bflag);
#endif
    }
#endif
}

/*----------------------------------------------------------------------;
 *
 * vector form of RHS
 *
 *----------------------------------------------------------------------*/
void dendroSolverRHSUnpacked(double **unzipVarsRHS, const double **uZipVars,
                             const unsigned int &offset, const double *pmin,
                             const double *pmax, const unsigned int *sz,
                             const unsigned int &bflag) {
    // clang-format off
    /*[[[cog
    import cog
    import sys
    import os
    import importlib.util
    import dendrosym

    cog.outl('// clang-format on')

    # get the current working directory, should be root of project
    output_folder_from_root = "generated-files"
    output_path = output_folder_from_root

    # the following lines will import any module directly from
    spec = importlib.util.spec_from_file_location("dendroconf", CONFIG_FILE_PATH)
    dendroconf = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = dendroconf
    spec.loader.exec_module(dendroconf)

    cog.outl("// EVOLUTION VARIABLE EXTRACTION NOT RHS")
    cog.outl(dendroconf.dendroConfigs.generate_variable_extraction("evolution",
    use_const=True))

    cog.outl("// EVOLUTION VARIABLE EXTRACTION RHS")
    cog.outl(dendroconf.dendroConfigs.generate_rhs_var_extraction("evolution",
    zip_var_name="unzipVarsRHS"))

    ]]]*/
    // clang-format on

    //[[[end]]]

    mem::memory_pool<double> *__mem_pool = &DENDROSOLVER_MEM_POOL;

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const double hx = (pmax[0] - pmin[0]) / (nx - 1);
    const double hy = (pmax[1] - pmin[1]) / (ny - 1);
    const double hz = (pmax[2] - pmin[2]) / (nz - 1);

    // get derivative workspace
    double *const deriv_base = dsolve::DENDROSOLVER_DERIV_WORKSPACE;

    // clang-format off
    /*[[[cog
    cog.outl('// clang-format on')
    cog.outl("// PARAMETER EXTRACTION FOR EVOLUTION")

    cog.outl(dendroconf.dendroConfigs.gen_parameter_code("evolution"))

    ]]]*/
    // clang-format on

    //[[[end]]]

    int idx[3];
    const unsigned int PW = dsolve::DENDROSOLVER_PADDING_WIDTH;
    unsigned int n = sz[0] * sz[1] * sz[2];
    unsigned int BLK_SZ = n;

    // declare the size of bytes for memory allocation down the line
    const unsigned int bytes = n * sizeof(double);

// Create the necessary pre-derivatives
// clang-format off
    /*[[[cog
    cog.outl('// clang-format on')

    cog.outl("//GENERATED ADVANCED DERIVATIVE EQUATIONS")

    print("Now generating advanced derivatves", file=sys.stderr)
    
    # note that we need to store the deallocation string as well for later down the line!
    (intermediate_grad_str, 
     deallocate_intermediate_grad_str) = dendroconf.dendroConfigs.generate_pre_necessary_derivatives(
         "evolution", dtype="double", include_byte_declaration=False
     )

    print("Finished generating advanced derivatves", file=sys.stderr)

    intermediate_filename = "bssn_rhs_intermediate_grad.cpp.inc"

    with open(os.path.join(output_path, intermediate_filename), "w") as f:
        f.write(intermediate_grad_str)

    print("Saved them to file", file=sys.stderr)

    cog.outl(f'#include "../{output_path}/{intermediate_filename}"')

    ]]]*/
    // clang-format on
    //GENERATED ADVANCED DERIVATIVE EQUATIONS
    #include "../generated-files/solver_rhs_intermediate_grad.cpp.inc"
    //[[[end]]]

    dsolve::timer::t_deriv.start();

// create the files that have the derivative memory allocations and
// calculations
// clang-format off
    /*[[[cog
    

    deriv_alloc, deriv_calc, deriv_dealloc = dendroconf.dendroConfigs.generate_deriv_allocation_and_calc("evolution", include_byte_declaration=False)

    print("Generated derivative allocation, calculation, and deallocation code for Evolution", file=sys.stderr)

    alloc_filename = "solver_rhs_deriv_memalloc.cpp.inc"

    with open(os.path.join(output_path, alloc_filename), "w") as f:
        f.write(deriv_alloc)
    
    cog.outl(f'#include "../{output_folder_from_root}/{alloc_filename}"')

    calc_filename = "solver_rhs_deriv_calc.cpp.inc"

    with open(os.path.join(output_path, calc_filename), "w") as f:
        f.write(deriv_calc)
    
    cog.outl(f'#include "../{output_folder_from_root}/{calc_filename}"')

    dealloc_filename = "solver_rhs_deriv_memdealloc.cpp.inc"

    with open(os.path.join(output_path, dealloc_filename), "w") as f:
        f.write(deriv_dealloc)

    cog.outl('// clang-format on')
    ]]]*/
    #include "../generated-files/solver_rhs_deriv_memalloc.cpp.inc"
    #include "../generated-files/solver_rhs_deriv_calc.cpp.inc"
    // clang-format on
    //[[[end]]]

    dsolve::timer::t_deriv.stop();

    // loop dep. removed allowing compiler to optmize for vectorization.
    // cout << "begin loop" << endl;
    dsolve::timer::t_rhs.start();
    for (unsigned int k = PW; k < nz - PW; k++) {
        for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef DENDROSOLVER_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = PW; i < nx - PW; i++) {
                const double x = pmin[0] + i * hx;
                const double y = pmin[1] + j * hy;
                const double z = pmin[2] + k * hz;

                const unsigned int pp = i + nx * (j + ny * k);
                const double r_coord = sqrt(x * x + y * y + z * z);

// TODO: this was here before! Is it something we need in
// BSSN??? double eta = ETA_CONST;

// if (r_coord >= ETA_R0)
// {
//     eta *= pow((ETA_R0 / r_coord), ETA_DAMPING_EXP);
// }

// clang-format off
                /*[[[cog
                cog.outl('// clang-format on')

                evolution_rhs_code = dendroconf.dendroConfigs.generate_rhs_code("evolution")
                evolution_filename = "solver_rhs_eqns.cpp.inc"

                with open(os.path.join(output_path, evolution_filename), "w") as f:
                    f.write(evolution_rhs_code)
                
                cog.outl(f'#include "../{output_folder_from_root}/{evolution_filename}"')
                
                ]]]*/
                // clang-format on
                #include "../generated-files/solver_rhs_eqns.cpp.inc"
                //[[[end]]]

                // /* debugging */
                // unsigned int qi = 46 - 1;
                // unsigned int qj = 10 - 1;
                // unsigned int qk = 60 - 1;
                // unsigned int qidx = qi + nx*(qj + ny*qk);
                // if (0 && qidx == pp) {
                //     std::cout << ".... end OPTIMIZED debug stuff..." <<
                //     std::endl;
                // }
            }
        }
    }
    dsolve::timer::t_rhs.stop();

    // Deallocate the pre-derivatives
    // TODO: is this the best place to put this? or should it reside at the end
    // with the rest of the freeing?
    /*[[[cog

    cog.outl("//GENERATED DEALLOCATION OF INTERMEDIATE GRAD CALCULATIONS")

    cog.outl(deallocate_intermediate_grad_str)

    ]]]*/

    //[[[end]]]

    if (bflag != 0) {
        dsolve::timer::t_bdyc.start();

        // clang-format off
        /*[[[cog
        cog.outl('// clang-format on')
        cog.outl(dendroconf.dendroConfigs.generate_bcs_calculations("evolution"))
        ]]]*/
        // clang-format on

        //[[[end]]]

        dsolve::timer::t_bdyc.stop();
    }

    dsolve::timer::t_deriv.start();
// TODO: include more types of build options
#include "../../generated-files/solver_rhs_ko_deriv_calc.cpp.inc"
    dsolve::timer::t_deriv.stop();

    dsolve::timer::t_rhs.start();

    const double sigma = KO_DISS_SIGMA;

    for (unsigned int k = PW; k < nz - PW; k++) {
        for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef DENDROSOLVER_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = PW; i < nx - PW; i++) {
                const unsigned int pp = i + nx * (j + ny * k);

                // clang-format off
                /*[[[cog
                cog.outl('// clang-format on')

                cog.outl("// GENERATED KO DISSIPATION CALCULATIONS")
                cog.outl(dendroconf.dendroConfigs.generate_ko_calculations("evolution"))

                ]]]*/
                // clang-format on

                //[[[end]]]
            }
        }
    }

    dsolve::timer::t_rhs.stop();

    dsolve::timer::t_deriv.start();
// clang-format off
    /*[[[cog
    cog.outl('// clang-format on')

    cog.outl(f'#include "../{output_folder_from_root}/{dealloc_filename}"')

    ]]]*/
    // clang-format on
    #include "../generated-files/solver_rhs_deriv_memdealloc.cpp.inc"
    //[[[end]]]

    dsolve::timer::t_deriv.stop();

#if 0
        for (unsigned int m = 0; m < 24; m++) {
            std::cout<<"  || dtu("<<m<<")|| = "<<normLInfty(unzipVarsRHS[m] + offset, n)<<std::endl;
        }
#endif
}

// TODO: this is where the sep function went

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void dendroSolverBCS(double *f_rhs, const double *f, const double *dxf,
                     const double *dyf, const double *dzf, const double *pmin,
                     const double *pmax, const double f_falloff,
                     const double f_asymptotic, const unsigned int *sz,
                     const unsigned int &bflag) {
    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const double hx = (pmax[0] - pmin[0]) / (nx - 1);
    const double hy = (pmax[1] - pmin[1]) / (ny - 1);
    const double hz = (pmax[2] - pmin[2]) / (nz - 1);

    const unsigned int PW = dsolve::DENDROSOLVER_PADDING_WIDTH;

    unsigned int ib = PW;
    unsigned int jb = PW;
    unsigned int kb = PW;
    unsigned int ie = sz[0] - PW;
    unsigned int je = sz[1] - PW;
    unsigned int ke = sz[2] - PW;

    double x, y, z;
    unsigned int pp;
    double inv_r;

    // std::cout<<"boundary dsolverhs: size [ "<<nx<<", "<<ny<<", "<<nz<<"
    // ]"<<std::endl; std::cout<<"boundary dsolverhs: pmin [ "<<pmin[0]<<",
    // "<<pmin[1]<<", "<<pmin[2]<<" ]"<<std::endl; std::cout<<"boundary
    // dsolverhs: pmax [ "<<pmax[0]<<", "<<pmax[1]<<", "<<pmax[2]<<"
    // ]"<<std::endl;

    if (bflag & (1u << OCT_DIR_LEFT)) {
        double x = pmin[0] + ib * hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int j = jb; j < je; j++) {
                y = pmin[1] + j * hy;
                pp = IDX(ib, j, k);
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }

    if (bflag & (1u << OCT_DIR_RIGHT)) {
        x = pmin[0] + (ie - 1) * hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int j = jb; j < je; j++) {
                y = pmin[1] + j * hy;
                pp = IDX((ie - 1), j, k);
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }

    if (bflag & (1u << OCT_DIR_DOWN)) {
        y = pmin[1] + jb * hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i * hx;
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);
                pp = IDX(i, jb, k);

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }

    if (bflag & (1u << OCT_DIR_UP)) {
        y = pmin[1] + (je - 1) * hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i * hx;
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);
                pp = IDX(i, (je - 1), k);

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }

    if (bflag & (1u << OCT_DIR_BACK)) {
        z = pmin[2] + kb * hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j * hy;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i * hx;
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);
                pp = IDX(i, j, kb);

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }

    if (bflag & (1u << OCT_DIR_FRONT)) {
        z = pmin[2] + (ke - 1) * hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j * hy;
            for (unsigned int i = ib; i < ie; i++) {
                x = pmin[0] + i * hx;
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);
                pp = IDX(i, j, (ke - 1));

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void max_spacetime_speeds(double *const lambda1max, double *const lambda2max,
                          double *const lambda3max, const double *const alpha,
                          const double *const beta1, const double *const beta2,
                          const double *const beta3, const double *const gtd11,
                          const double *const gtd12, const double *const gtd13,
                          const double *const gtd22, const double *const gtd23,
                          const double *const gtd33, const double *const chi,
                          const unsigned int *sz) {
    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const unsigned int PW = dsolve::DENDROSOLVER_PADDING_WIDTH;
    unsigned int ib = PW;
    unsigned int jb = PW;
    unsigned int kb = PW;
    unsigned int ie = sz[0] - PW;
    unsigned int je = sz[1] - PW;
    unsigned int ke = sz[2] - PW;

    for (unsigned int k = kb; k < ke; k++) {
        for (unsigned int j = jb; j < je; j++) {
            for (unsigned int i = ib; i < ie; i++) {
                unsigned int pp = IDX(i, j, k);
                /* note: gtu is the inverse tilde metric. It should have detgtd
                 * = 1. So, for the purposes of
                 * calculating wavespeeds, I simple set detgtd = 1. */
                double gtu11 = gtd22[pp] * gtd33[pp] - gtd23[pp] * gtd23[pp];
                double gtu22 = gtd11[pp] * gtd33[pp] - gtd13[pp] * gtd13[pp];
                double gtu33 = gtd11[pp] * gtd22[pp] - gtd12[pp] * gtd12[pp];
                if (gtu11 < 0.0 || gtu22 < 0.0 || gtu33 < 0.0) {
                    std::cout << "Problem computing spacetime characteristics"
                              << std::endl;
                    std::cout << "gtu11 = " << gtu11 << ", gtu22 = " << gtu22
                              << ", gtu33 = " << gtu33 << std::endl;
                    gtu11 = 1.0;
                    gtu22 = 1.0;
                    gtu33 = 1.0;
                }
                double t1 = alpha[pp] * sqrt(gtu11 * chi[pp]);
                double t2 = alpha[pp] * sqrt(gtu22 * chi[pp]);
                double t3 = alpha[pp] * sqrt(gtu33 * chi[pp]);
                lambda1max[pp] =
                    std::max(abs(-beta1[pp] + t1), abs(-beta1[pp] - t1));
                lambda2max[pp] =
                    std::max(abs(-beta2[pp] + t2), abs(-beta2[pp] - t2));
                lambda3max[pp] =
                    std::max(abs(-beta3[pp] + t3), abs(-beta3[pp] - t3));
            }
        }
    }
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void freeze_bcs(double *f_rhs, const unsigned int *sz,
                const unsigned int &bflag) {
    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const unsigned int PW = dsolve::DENDROSOLVER_PADDING_WIDTH;
    unsigned int ib = PW;
    unsigned int jb = PW;
    unsigned int kb = PW;
    unsigned int ie = sz[0] - PW;
    unsigned int je = sz[1] - PW;
    unsigned int ke = sz[2] - PW;

    unsigned int pp;

    if (bflag & (1u << OCT_DIR_LEFT)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int j = jb; j < je; j++) {
                pp = IDX(ib, j, k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u << OCT_DIR_RIGHT)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int j = jb; j < je; j++) {
                pp = IDX((ie - 1), j, k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u << OCT_DIR_DOWN)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i, jb, k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u << OCT_DIR_UP)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i, (je - 1), k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u << OCT_DIR_BACK)) {
        for (unsigned int j = jb; j < je; j++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i, j, kb);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u << OCT_DIR_FRONT)) {
        for (unsigned int j = jb; j < je; j++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp = IDX(i, j, (ke - 1));
                f_rhs[pp] = 0.0;
            }
        }
    }
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/

void dendroSolverRHSUnpacked_sep(double **unzipVarsRHS, const double **uZipVars,
                                 const unsigned int &offset, const double *pmin,
                                 const double *pmax, const unsigned int *sz,
                                 const unsigned int &bflag) {
    // TODO: generate the separate version of the code
}
