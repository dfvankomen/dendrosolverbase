//
// Created by milinda on 10/8/18.
//

#ifndef DENDROSOLVER_GRDEF_H_
#define DENDROSOLVER_GRDEF_H_

// TODO: potentially move these to inline function calls, macros aren't
// necessary
#define Rx \
    (dsolve::DENDROSOLVER_COMPD_MAX[0] - dsolve::DENDROSOLVER_COMPD_MIN[0])
#define Ry \
    (dsolve::DENDROSOLVER_COMPD_MAX[1] - dsolve::DENDROSOLVER_COMPD_MIN[1])
#define Rz \
    (dsolve::DENDROSOLVER_COMPD_MAX[2] - dsolve::DENDROSOLVER_COMPD_MIN[2])

#define RgX \
    (dsolve::DENDROSOLVER_OCTREE_MAX[0] - dsolve::DENDROSOLVER_OCTREE_MIN[0])
#define RgY \
    (dsolve::DENDROSOLVER_OCTREE_MAX[1] - dsolve::DENDROSOLVER_OCTREE_MIN[1])
#define RgZ \
    (dsolve::DENDROSOLVER_OCTREE_MAX[2] - dsolve::DENDROSOLVER_OCTREE_MIN[2])

#define GRIDX_TO_X(xg)                                          \
    (((Rx / RgX) * (xg - dsolve::DENDROSOLVER_OCTREE_MIN[0])) + \
     dsolve::DENDROSOLVER_COMPD_MIN[0])
#define GRIDY_TO_Y(yg)                                          \
    (((Ry / RgY) * (yg - dsolve::DENDROSOLVER_OCTREE_MIN[1])) + \
     dsolve::DENDROSOLVER_COMPD_MIN[1])
#define GRIDZ_TO_Z(zg)                                          \
    (((Rz / RgZ) * (zg - dsolve::DENDROSOLVER_OCTREE_MIN[2])) + \
     dsolve::DENDROSOLVER_COMPD_MIN[2])

#define X_TO_GRIDX(xc)                                         \
    (((RgX / Rx) * (xc - dsolve::DENDROSOLVER_COMPD_MIN[0])) + \
     dsolve::DENDROSOLVER_OCTREE_MIN[0])
#define Y_TO_GRIDY(yc)                                         \
    (((RgY / Ry) * (yc - dsolve::DENDROSOLVER_COMPD_MIN[1])) + \
     dsolve::DENDROSOLVER_OCTREE_MIN[1])
#define Z_TO_GRIDZ(zc)                                         \
    (((RgZ / Rz) * (zc - dsolve::DENDROSOLVER_COMPD_MIN[2])) + \
     dsolve::DENDROSOLVER_OCTREE_MIN[2])

// type of the rk method.
enum RKType { RK3, RK4, RK45 };

namespace dsolve {

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

cog.outl(dendroconf.dendroConfigs.gen_enum_code("evolution"))

cog.outl(dendroconf.dendroConfigs.gen_enum_code("constraint", enum_name="VAR_CONSTRAINT"))

]]]*/
// clang-format on

//[[[end]]]

// clang-format off
/*[[[cog
cog.outl('// clang-format on')
cog.outl(dendroconf.dendroConfigs.gen_enum_names("evolution"))
cog.outl(dendroconf.dendroConfigs.gen_enum_names("constraint", enum_name="VAR_CONSTRAINT"))
cog.outl(dendroconf.dendroConfigs.gen_enum_iterable_list("evolution"))
cog.outl(dendroconf.dendroConfigs.gen_enum_iterable_list("constraint", enum_name="VAR_CONSTRAINT"))
]]]*/
// clang-format on

//[[[end]]]

/**
 * @brief Refinement mode types.
 * WAMR : Wavelet based refinement.
 * EH : black hole event horizon based refinement.
 * EH_WAMR: both even horizon as well as WAMR based refinement.
 * BH_LOC BH location based refinement, if turned on track the bh locations.
 */
enum RefinementMode {
    WAMR = 0,
    EH,
    EH_WAMR,
    BH_LOC
};  // TODO: generate the different types of refinement, for now they're these

}  // end of namespace dsolve

#endif  // DENDROSOLVER_GRDEF_H_
