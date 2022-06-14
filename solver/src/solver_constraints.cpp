#include "solver_constraints.h"

#include <limits>

using namespace dsolve;

/*----------------------------------------------------------------------;
 *
 * enforce physical constraints on variables:
 *
 *
 *----------------------------------------------------------------------*/
void enforce_solver_constraints(double **uiVar, const unsigned int node) {
    const double one_third = 1.0 / 3.0;

    // clang-format off
    /*[[[cog

    import cog
    import sys
    import importlib.util
    import dendrosym

    cog.outl('// clang-format on')

    # the following lines will import any module directly from
    spec = importlib.util.spec_from_file_location("dendroconf", CONFIG_FILE_PATH)
    dendroconf = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = dendroconf
    spec.loader.exec_module(dendroconf)

    cog.outl('//// EMDA CONSTRAINTS')
    cog.outl(dendroconf.dendroConfigs.generate_evolution_constraints())

    ]]]*/
    // clang-format on

    //[[[end]]]
}
