# ##############################################################################
# CMAKE CONFIGURATION FILE FOR GR SOLVER
#
# This file in particular contains all of the options that apply specificially
# to the GR Solver binary as it is compiled.
#
# Note that all options and compiler flags set in the ../CMakeLists.txt file
# (the one in the root directory of the project) *are carried through* here.
# Meaning, at this point, the C++ standard is set, OpenMP, MPI, GSL, and other
# libraries are found and linked. Also, general Dendro options, such as ordering
# options, dimensionality, CMake's build type, and others will also propagate
# through to here.
# ##############################################################################

# The SOLVER solver now requires the toml11 libary, now included in /lib this line
# tells CMAKE to link them
set(EXTERNAL_LIB_HEADERS ${CMAKE_SOURCE_DIR}/lib/toml11)

# ===========================
# == SOLVER SPECIFIC OPTIONS == These are options *specificaly* used for compiling
# and modifying the project
# ===========================
option(COG_SOURCE_FILES "Enable cogging of source files" OFF)
option(DENDROSOLVER_PROFILE_HUMAN_READABLE "Enable human readable profile outputs" ON)
option(
  DENDROSOLVER_CONSEC_COMM_SELECT
  "Rank selection is performed using consecutive ranks (works with an arbitrary number of processes). If off, the ranks are selected in binary tree fashion (only with powers of two processes)"
  ON)
option(DENDROSOLVER_ENABLE_VTU_OUTPUT "Enables VTU Output" ON)
option(DENDROSOLVER_COMPUTE_CONSTRAINTS "Enables constraint equation computation" ON)
option(DENDROSOLVER_ENABLE_VTU_CONSTRAINT_OUTPUT "Writes constraint output to VTU file"
       ON)
option(DENDROSOLVER_ETA_FUNCTION "Use function as ETA damping" OFF)
option(DENDROSOLVER_KERR_SCHILD_TEST
       "Compare with the analytical Kerr-Schild solution for a single BH" OFF)
option(DENDROSOLVER_ENABLE_CUDA "Enable RHS computation with GPU acceleration" OFF)
option(DENDROSOLVER_EXTRACT_BH_LOCATIONS
       "Compute the locations of the BH as time evolves" OFF)
option(DENDROSOLVER_REFINE_BASE_EH "Enable refinements based on BH event horizon" OFF)
option(DENDROSOLVER_EXTRACT_GRAVITATIONAL_WAVES "Extract the gravitational waves" OFF)
option(DENDROSOLVER_USE_4TH_ORDER_DERIVS
       "Use 4th order derivative stencil calculations" OFF)
option(DENDROSOLVER_USE_6TH_ORDER_DERIVS
       "Use 6th order derivative stencil calculations" ON)
option(DENDROSOLVER_USE_8TH_ORDER_DERIVS
       "Use 8th order derivative stencil calculations" OFF)
option(
  DENDROSOLVER_ENABLE_AVX
  "Use vectorized computations for the RHS calculations (intel compilers ONLY)"
  OFF)

# setup for Intel-based vectorized computations
if(DENDROSOLVER_ENABLE_AVX)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    add_definitions(-DDENDROSOLVER_ENABLE_AVX)
  else()
    message(
      "WARNING: AVX vectorization not enabled due to detection of non-Intel compiler. Use `export CC=icc` and `export CXX=icpc` to use an Intel compiler."
    )
  endif()
endif(DENDROSOLVER_ENABLE_AVX)

# NOTE: adjustments to the derivative orders based on priority. If 4th order is
# enabled, it'll ignore the settings for later options and so on
if(DENDROSOLVER_USE_4TH_ORDER_DERIVS)
  set(DENDROSOLVER_USE_6TH_ORDER_DERIVS OFF)
  set(DENDROSOLVER_USE_8TH_ORDER_DERIVS OFF)
elseif(DENDROSOLVER_USE_6TH_ORDER_DERIVS)
  set(DENDROSOLVER_USE_4TH_ORDER_DERIVS OFF)
  set(DENDROSOLVER_USE_8TH_ORDER_DERIVS OFF)
elseif(DENDROSOLVER_USE_8TH_ORDER_DERIVS)
  set(DENDROSOLVER_USE_4TH_ORDER_DERIVS OFF)
  set(DENDROSOLVER_USE_6TH_ORDER_DERIVS OFF)
endif()

# now add the definitions for the compiler so that these options can be accessed
# by the code
if(DENDROSOLVER_USE_4TH_ORDER_DERIVS)
  add_definitions(-DDENDROSOLVER_USE_4TH_ORDER_DERIVS)
endif()

if(DENDROSOLVER_USE_6TH_ORDER_DERIVS)
  add_definitions(-DDENDROSOLVER_USE_6TH_ORDER_DERIVS)
endif()

if(DENDROSOLVER_USE_8TH_ORDER_DERIVS)
  add_definitions(-DDENDROSOLVER_USE_8TH_ORDER_DERIVS)
endif()

if(DENDROSOLVER_PROFILE_HUMAN_READABLE)
  add_definitions(-DDENDROSOLVER_PROFILE_HUMAN_READABLE)
endif()

if(DENDROSOLVER_CONSEC_COMM_SELECT)
  add_definitions(-DDENDROSOLVER_CONSEC_COMM_SELECT)
endif()

if(DENDROSOLVER_ENABLE_VTU_OUTPUT)
  add_definitions(-DDENDROSOLVER_ENABLE_VTU_OUTPUT)
endif()

if(DENDROSOLVER_COMPUTE_CONSTRAINTS)
  add_definitions(-DDENDROSOLVER_COMPUTE_CONSTRAINTS)
endif()

if(DENDROSOLVER_ENABLE_VTU_CONSTRAINT_OUTPUT)
  add_definitions(-DDENDROSOLVER_ENABLE_VTU_CONSTRAINT_OUTPUT)
endif()

if(DENDROSOLVER_ETA_FUNCTION)
  add_definitions(-DUSE_ETA_FUNC)
endif()

if(DENDROSOLVER_KERR_SCHILD_TEST)
  add_definitions(-DDENDROSOLVER_KERR_SCHILD_TEST)
endif()

if(DENDROSOLVER_ENABLE_CUDA)
  add_definitions(-DDENDROSOLVER_ENABLE_CUDA)
endif()

# TODO: is this necessary? if(DENDROSOLVER_RHS_STAGED_COMP)
# add_definitions(-DDENDROSOLVER_RHS_STAGED_COMP) endif()

if(DENDROSOLVER_EXTRACT_BH_LOCATIONS)
  add_definitions(-DDENDROSOLVER_EXTRACT_BH_LOCATIONS)
endif()

if(DENDROSOLVER_REFINE_BASE_EH)
  add_definitions(-DDENDROSOLVER_REFINE_BASE_EH)
endif()

# if(DENDROSOLVER_GAUGE_ROCHESTER) add_definitions(-DUSE_ROCHESTER_GAUGE) endif()

set(DENDROSOLVER_INC
    ${CMAKE_SOURCE_DIR}/solver/include/solver_main.h
    ${CMAKE_SOURCE_DIR}/solver/include/parameters.h
    ${CMAKE_SOURCE_DIR}/solver/include/bh.h
    ${CMAKE_SOURCE_DIR}/solver/include/grUtils.h
    ${CMAKE_SOURCE_DIR}/solver/include/grUtils.tcc
    ${CMAKE_SOURCE_DIR}/solver/include/rhs.h
    ${CMAKE_SOURCE_DIR}/solver/include/derivs.h
    ${CMAKE_SOURCE_DIR}/solver/include/physcon.h
    ${CMAKE_SOURCE_DIR}/solver/include/profile_params.h
    ${CMAKE_SOURCE_DIR}/solver/include/solver_constraints.h
    ${CMAKE_SOURCE_DIR}/solver/include/TPUtilities.h
    ${CMAKE_SOURCE_DIR}/solver/include/TwoPunctures.h
    ${CMAKE_SOURCE_DIR}/solver/include/dataUtils.h
    # ${CMAKE_SOURCE_DIR}/include/aeh_solver.h
    ${CMAKE_SOURCE_DIR}/solver/include/solverCtx.h)

set(DENDROSOLVER_SRC
    src/rkSolver.cpp
    src/parameters.cpp
    src/grUtils.cpp
    src/rhs.cpp
    src/derivs.cpp
    src/physcon.cpp
    src/profile_params.cpp
    src/solver_constraints.cpp
    src/TPCoordTransf.cpp
    src/TPEquations.cpp
    src/TPFuncAndJacobian.cpp
    src/TPNewton.cpp
    src/TPUtilities.cpp
    src/TwoPunctures.cpp
    src/dataUtils.cpp
    src/solverCtx.cpp)

set(SOURCE_FILES solver_main.cpp
                 ${DENDROSOLVER_INC} ${DENDROSOLVER_SRC})
add_executable(grSolver ${SOURCE_FILES})

# source directory is root of this project since this file is only called from
# the outer CMake call
target_include_directories(grSolver PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(grSolver PRIVATE ${CMAKE_SOURCE_DIR}/solver)
target_include_directories(grSolver PRIVATE ${CMAKE_SOURCE_DIR}/solver/include)
target_include_directories(grSolver PRIVATE ${CMAKE_SOURCE_DIR}/scripts)
target_include_directories(grSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendro)
target_include_directories(grSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendro/include)
target_include_directories(grSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendro/test/include)
target_include_directories(grSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendro/examples/include)
target_include_directories(grSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendro/ODE/include)
target_include_directories(grSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendro/FEM/include)
target_include_directories(grSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendro/LinAlg/include)
target_include_directories(grSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendro/IO/include)
target_include_directories(grSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendro/IO/vtk/include)
target_include_directories(grSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendro/LinAlg/include)
target_include_directories(grSolver
                           PRIVATE ${CMAKE_DOURCE_DIR}/dendro/IO/zlib/inc)

# and don't forget to include the MPI and GSL paths
target_include_directories(grSolver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(grSolver PRIVATE ${GSL_INCLUDE_DIRS})

# also the toml library with the git submodule
target_include_directories(grSolver PRIVATE ${EXTERNAL_LIB_HEADERS})

if(DENDROSOLVER_EXTRACT_GRAVITATIONAL_WAVES)
  add_definitions(-DDENDROSOLVER_EXTRACT_GRAVITATIONAL_WAVES)
endif(DENDROSOLVER_EXTRACT_GRAVITATIONAL_WAVES)

# then don't forget to link the libraries!!!!
target_link_libraries(grSolver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES}
                      ${GSL_LIBRARIES} m)

# then add the dependency add custom command to generate the source files from
# cog note that this doesn't actually *output* files, it actually will read
# existing source and update them if there are changes
if(COG_SOURCE_FILES)
  # NOTE: the OUTPUT is set to rhs.cpp so that it'll always generate because
  # rhs.cpp is necessary for the project

  if(DEFINED $ENV{DENDROSOLVER_SOURCE_PYTHON})
    set(BASE_PYTHON_SCRIPT ${CMAKE_SOURCE_DIR}/$ENV{DENDROSOLVER_SOURCE_PYTHON})
  else()
    # sets the default python script to be the sample included
    set(BASE_PYTHON_SCRIPT ${CMAKE_SOURCE_DIR}/solver_eqns_configs.py)
  endif()

  message("Will be generating code with the Python file: ${BASE_PYTHON_SCRIPT}")

  # then we also need to define the base parameters file
  if(DEFINED $ENV{DENDROSOLVER_SOURCE_ENV_FILE})
    set(BASE_PARAMETER_FILE ${CMAKE_SOURCE_DIR}/$ENV{DENDROSOLVER_SOURCE_ENV_FILE})
  else()
    set(BASE_PARAMETER_FILE ${CMAKE_SOURCE_DIR}/solver_param_setup.toml)
  endif()

  message("Will be generating parameter code with the file: ${BASE_PARAMETER_FILE}")

  add_custom_command(
    OUTPUT ${CMAKE_SOURCE_DIR}/src/dummy.cpp
    DEPENDS ${CMAKE_SOURCE_DIR}/solver_eqns_configs.py
    COMMAND
      python3 -m cogapp -D
      CONFIG_FILE_PATH=${BASE_PYTHON_SCRIPT} -D
      PARAM_SETUP_FILE=${BASE_PARAMETER_FILE} -r
      @${CMAKE_SOURCE_DIR}/files_to_generate.txt
    COMMENT
      "Generating code based on the configurations in solver_eqns_configs.py. This may take a while as it optimizes."
    VERBATIM USES_TERMINAL
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  add_custom_target(COG_GENERATION DEPENDS ${CMAKE_SOURCE_DIR}/src/dummy.cpp)

  add_dependencies(grSolver COG_GENERATION)

endif()
