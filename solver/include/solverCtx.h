/**
 * @file dsolveCtx.h
 * @author Milinda Fernando
 * @author David Van Komen
 * @brief Application context class for solving the Einstein equations a
 * formulation.
 * @version 0.1
 * @date 2019-12-20
 *
 * @copyright Copyright (c) 2019, University of Utah.
 *
 */

#ifndef DENDROSOLVER_SOLVERCTX_H_
#define DENDROSOLVER_SOLVERCTX_H_

#include <iostream>

#include "TwoPunctures.h"
#include "checkPoint.h"
#include "ctx.h"
#include "dataUtils.h"
#include "grDef.h"
#include "grUtils.h"
#include "gwExtract.h"
#include "oct2vtk.h"
#include "parUtils.h"
#include "parameters.h"
#include "physcon.h"
#include "rhs.h"
#include "solver_constraints.h"

namespace dsolve {

/**@brief smoothing modes avail for LTS recomended for LTS time stepping. */
enum LTS_SMOOTH_MODE { KO = 0, WEIGHT_FUNC };
enum VL {
    CPU_EV = 0,
    CPU_CV,
    CPU_EV_UZ_IN,
    CPU_EV_UZ_OUT,
    CPU_CV_UZ_IN,
    GPU_EV,
    GPU_EV_UZ_IN,
    GPU_EV_UZ_OUT,
    END
};
typedef ot::DVector<DendroScalar, unsigned int> DVec;

class SOLVERCtx : public ts::Ctx<SOLVERCtx, DendroScalar, unsigned int> {
   protected:
            /**@brief: evolution var (zip)*/
            DVec m_var[VL::END];

            Point m_uiBHLoc[2];

   public:
    /**@brief: default constructor*/
    SOLVERCtx(ot::Mesh *pMesh);

    /**@brief: default deconstructor*/
    ~SOLVERCtx();

    /**
     * @brief sets time adaptive offset
     * @param tadapoffst
     */
    void set_time_adap_offset(unsigned int tadapoffst) {
        DENDROSOLVER_LTS_TS_OFFSET = tadapoffst;
    }

    /**@brief : returns the time adaptive offset value*/
    unsigned int get_time_adap_offset() { return DENDROSOLVER_LTS_TS_OFFSET; }

    /**@brief: initial solution and grid convergence calls init_grid()*/
    int initialize();

    /**@brief: initialize the grid, solution. */
    int init_grid();

    /**
     * @brief computes the Right Hand Side
     *
     * @param in : zipped input
     * @param out : zipped output
     * @param sz  : number of variables.
     * @param time : current time.
     * @return int : status. (0) on success.
     */
    int rhs(DVec *in, DVec *out, unsigned int sz, DendroScalar time);

    /**@brief: function execute before each stage
     * @param sIn: stage var in.
     */
    int pre_stage(DVec& sIn);

    /**@brief: function execute after each stage
     * @param sIn: stage var in.
     */
    int post_stage(DVec& sIn);

    /**@brief: function execute before each step*/
    int pre_timestep(DVec& sIn);

    /**@brief: function execute after each step*/
    int post_timestep(DVec& sIn);

    /**@brief: function execute after each step*/
    bool is_remesh();

    /**@brief: write to vtu. */
    int write_vtu();

    /**@brief: writes checkpoint*/
    int write_checkpt();

    /**@brief: restore from check point*/
    int restore_checkpt();

    /**@brief: should be called for free up the contex memory. */
    int finalize();

    /**@brief: pack and returns the evolution variables to one DVector*/
    DVec& get_evolution_vars();

    /**@brief: pack and returns the constraint variables to one DVector*/
    DVec& get_constraint_vars();

    /**@brief: pack and returns the primitive variables to one DVector*/
    DVec& get_primitive_vars();

    /**@brief: prints any messages to the terminal output. */
    int terminal_output();

    /**@brief: returns the async communication batch size. */
    unsigned int get_async_batch_sz() {
        return dsolve::DENDROSOLVER_ASYNC_COMM_K;
    }

    /**@brief: returns the number of variables considered when performing
     * refinement*/
    unsigned int get_num_refine_vars() { return DENDROSOLVER_NUM_REFINE_VARS; }

    /**@brief: return the pointer for containing evolution refinement variable
     * ids*/
    const unsigned int *get_refine_var_ids() {
        return DENDROSOLVER_REFINE_VARIABLE_INDICES;
    }

    /**@brief return the wavelet tolerance function / value*/
    std::function<double(double, double, double, double *hx)>
    get_wtol_function() {
        double wtol = DENDROSOLVER_WAVELET_TOL;
        // TODO: basic WTolDCoords (currently set to black hole stuff)
        std::function<double(double, double, double, double *)> waveletTolFunc =
            [](double x, double y, double z, double *hx) {
                return dsolve::computeWTolDCoords(x, y, z, hx);
            };
        return waveletTolFunc;
    }

    /**@brief computes the LTS TS offset based on the eta damping.*/
    unsigned int compute_lts_ts_offset();

    /**@brief : blk time step factor. */
    static unsigned int getBlkTimestepFac(unsigned int blev, unsigned int lmin,
                                          unsigned int lmax);

    /**@biref: evolve bh locations. */
    void evolve_bh_loc(DVec sIn, double dt);

    /**
     * @brief LTS smooth mode.
     *
     * @param sIn : time synced evolution vector.
     * @param mode : smoothing mode.
     */
    void lts_smooth(DVec sIn, LTS_SMOOTH_MODE mode);

    /**@brief: return true if the BH are merged. */
    bool is_bh_merged(double tol) const {
        return (
            (dsolve::DENDROSOLVER_BH_LOC[0] - dsolve::DENDROSOLVER_BH_LOC[1])
                .abs() < tol);
    };

    int grid_transfer(const ot::Mesh* m_new);
};

}  // end of namespace dsolve

#endif
