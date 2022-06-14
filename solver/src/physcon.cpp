#include "physcon.h"

#include "solver_main.h"

using namespace dsolve;

/*----------------------------------------------------------------------
 *
 * vector form of RHS
 *
 *----------------------------------------------------------------------*/
void physical_constraints(double **uZipConVars, const double **uZipVars,
                          const unsigned int &offset, const double *pmin,
                          const double *pmax, const unsigned int *sz,
                          const unsigned int &bflag) {
    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];
    const unsigned int n = nx * ny * nz;

    // declare the size of bytes for memory allocation down the line
    const unsigned int bytes = n * sizeof(double);

    const double hx = (pmax[0] - pmin[0]) / (nx - 1);
    const double hy = (pmax[1] - pmin[1]) / (ny - 1);
    const double hz = (pmax[2] - pmin[2]) / (nz - 1);

    // HAM is hamiltonian constraint
    // mom is 3 component vector - momentum constraints

    // clang-format off
    /*[[[cog
    import cog
    import sys
    import os
    import importlib.util
    import dendrosym

    cog.outl('// clang-format on')

    # get the current working directory, should be root of project
    current_path = os.getcwd()
    output_path = os.path.join(current_path, "generated_files")

    # the following lines will import any module directly from
    spec = importlib.util.spec_from_file_location("dendroconf",
    CONFIG_FILE_PATH) dendroconf = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = dendroconf
    spec.loader.exec_module(dendroconf)

    cog.outl('//// GENERATED CONSTRAINT VARIABLES')
    cog.outl(dendroconf.dendroConfigs.generate_variable_extraction(
        "constraint", use_const=False, zip_var_name="uZipConVars"
    ))

    cog.outl('/// GENERATED CONSTANT EVOLUTION VARIABLES')
    cog.outl(dendroconf.dendroConfigs.generate_variable_extraction(
        "evolution", use_const=True, zip_var_name="uZipVars"
    ))
    ]]]*/
    // clang-format on

    //[[[end]]]

    const unsigned int PW = dsolve::DENDROSOLVER_PADDING_WIDTH;

    // clang-format off
    /*[[[cog
    cog.outl('// clang-format on')
    cog.outl("// PARAMETER EXTRACTION FOR CONSTRAINTS")

    cog.outl(dendroconf.dendroConfigs.gen_parameter_code("constraint"))

    ]]]*/
    // clang-format on

    //[[[end]]]

    mem::memory_pool<double> *__mem_pool = &DENDROSOLVER_MEM_POOL;

    // get derivative workspace
    double *const deriv_base = emda::EMDA_DERIV_WORKSPACE;

    // clang-format off
    /*[[[cog
    cog.outl('// clang-format on')
    
    cog.outl("//GENERATED ADVANCED DERIVATIVE EQUATIONS")

    print("Now generating advanced derivatves for physcon", file=sys.stderr)
    
    # note that we need to store the deallocation string as well for later down the line!
    (intermediate_grad_str, 
     deallocate_intermediate_grad_str) = dendroconf.dendroConfigs.generate_pre_necessary_derivatives(
         "constraint", dtype="double"
     )

    print("Finished generating advanced derivatves", file=sys.stderr)

    intermediate_filename = "emda_physcon_intermediate_grad.cpp.inc"

    with open(os.path.join(output_path, intermediate_filename), "w") as f:
        f.write(intermediate_grad_str)

    print("Saved them to file", file=sys.stderr)

    cog.outl(f'#include "../gencode/{intermediate_filename}"')
    
    ]]]*/
    // clang-format on

    //[[[end]]]

    // create the files that have the derivative memory allocations and
    // calculations
    // clang-format off
    /*[[[cog
    cog.outl('// clang-format on')

    deriv_alloc, deriv_calc, deriv_dealloc = dendroconf.dendroConfigs.generate_deriv_allocation_and_calc("constraint")

    print("Generated derivative allocation, calculation, and deallocation code for Evolution", file=sys.stderr)

    alloc_filename = "emda_physcon_deriv_memalloc.cpp.inc"

    with open(os.path.join(output_path, alloc_filename), "w") as f:
        f.write(deriv_alloc)

    cog.outl(f'#include "../gencode/{alloc_filename}"')

    calc_filename = "emda_physcon_deriv_calc.cpp.inc"

    with open(os.path.join(output_path, calc_filename), "w") as f:
        f.write(deriv_calc)

    cog.outl(f'#include "../gencode/{calc_filename}"')

    dealloc_filename = "emda_physcon_deriv_memdealloc.cpp.inc"

    with open(os.path.join(output_path, dealloc_filename), "w") as f:
        f.write(deriv_dealloc)

    ]]]*/
    // clang-format on

    //[[[end]]]

    // enforce hamiltonian and momentum constraints
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

                // clang-format off
                /*[[[cog

                cog.outl('// clang-format on')

                physcon_rhs_code = dendroconf.dendroConfigs.generate_rhs_code("constraint", include_rhs_in_name=False)
                physcon_filename = "emda_physcon_eqns.cpp.inc"

                with open(os.path.join(output_path, physcon_filename), "w") as f:
                    f.write(physcon_rhs_code)
                
                cog.outl(f'#include "../gencode/{physcon_filename}"')
                
                ]]]*/
                // clang-format on

                //[[[end]]]

                // TODO: represent this somehow with Python
                if (fabs(x) <= 1e-7 && fabs(y) <= 1e-7) {
                    std::cerr << "ABS OF X AND Y <= 1e-7" << std::endl;
                    psi4_real[pp] = 0.0;
                    psi4_imag[pp] = 0.0;
                }
            }
        }
    }

    // clang-format off
    /*[[[cog

    cog.outl('// clang-format on')

    cog.outl(deallocate_intermediate_grad_str)

    cog.outl(f'#include "../gencode/{dealloc_filename}"')
    

    ]]]*/
    // clang-format on

    //[[[end]]]
}
