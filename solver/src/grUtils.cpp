//
// Created by milinda on 7/26/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief Contains utility functions for EMDA simulation.
 */
//

// TODO: multiple things to do here, but the big idea is that the puncture needs
// to be defined *somehow* for all GR types I worry that if we change var names
// then we'll have to edit this stuff by hand

// TODO: read param file needs to be modified whenever there are new params as
// well

#include "grUtils.h"

namespace dsolve {

// NOTE: the read param file and dump param file are now included in
// parameters.cpp

void initDataFuncToPhysCoords(const double xx1, const double yy1,
                              const double zz1, double *var) {
    // convert "grid" values to physical values
    const double xx = GRIDX_TO_X(xx1);
    const double yy = GRIDY_TO_Y(yy1);
    const double zz = GRIDZ_TO_Z(zz1);

    // TODO: perform a "switch" here for different initialization functions?
    minkowskiInit(xx, yy, zz, var);
}

// clang-format off
/*[[[cog
import cog
import sys
import os
import importlib.util
import dendrosym
cog.outl("// clang-format off")

# get the current working directory, should be root of project
current_path = os.getcwd()
output_path = os.path.join(current_path, "gencode")

# the following lines will import any module directly from
spec = importlib.util.spec_from_file_location("dendroconf", CONFIG_FILE_PATH)
dendroconf = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = dendroconf
spec.loader.exec_module(dendroconf)

cog.outl("// INITIAL DATA FUNCTIONS")
cog.outl(dendroconf.dendroConfigs.generate_initial_data_code(var_type="evolution"))
]]]*/
// clang-format off

//[[[end]]]

void punctureDataPhysicalCoord(const double xx, const double yy,
                               const double zz, double *var) {
    // NOTE: this is currently defined for BSSN and will likely not work with
    // other formulations!
#if 0
        /* Define the Levi-Cevita pseudo-tensor and Kroneckar delta */
        double epijk[3][3][3];
        int i, j, k;
        for (k = 0; k < 3; k++)
        {
            for (j = 0; j < 3; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    epijk[k][j][i] = 0.0;
                }
            }
        }
        epijk[0][1][2] = 1.0;
        epijk[1][2][0] = 1.0;
        epijk[2][0][1] = 1.0;
        epijk[0][2][1] = -1.0;
        epijk[2][1][0] = -1.0;
        epijk[1][0][2] = -1.0;

        double deltaij[3][3];
        for (j = 0; j < 3; j++)
        {
            for (i = 0; i < 3; i++)
            {
                deltaij[j][i] = 0.0;
            }
        }

        deltaij[0][0] = 1.0;
        deltaij[1][1] = 1.0;
        deltaij[2][2] = 1.0;

        double x1, y1, z1, rv1;
        double x2, y2, z2, rv2;
        double vn1[3], vn2[3];

        double vpsibl;
        double v_u_corr, amp_capj, amp_capr, l_r, u0_j, u2_j, mu_j, p2_mu_j, v_u_j1;
        double v1, v2, v3, v4, vt1, vt2;

        int i1, i2, i3, i4;
        double amp_capp, u0_p, u2_p, mu_p, p2_mu_p;
        double v_u_p1, v_u_c1, v_u_j2, v_u_p2;
        double v_u_c2, vpsibl_u, vpsibl_u2;

        // bh 1
        double mass1 = BH1.getBHMass();
        double bh1x = BH1.getBHCoordX();
        double bh1y = BH1.getBHCoordY();
        double bh1z = BH1.getBHCoordZ();

        double vp1[3];
        vp1[0] = BH1.getVx();
        vp1[1] = BH1.getVy();
        vp1[2] = BH1.getVz();

        double vp1tot = sqrt(vp1[0] * vp1[0] + vp1[1] * vp1[1] + vp1[2] * vp1[2]);
        double spin1 = BH1.getBHSpin();
        double spin1_th = BH1.getBHSpinTheta();
        double spin1_phi = BH1.getBHSpinPhi();
        double vs1[3];

        vs1[0] = spin1 * sin(spin1_th) * cos(spin1_phi);
        vs1[1] = spin1 * sin(spin1_th) * sin(spin1_phi);
        vs1[2] = spin1 * cos(spin1_th);

        // bh 2
        double mass2 = BH2.getBHMass();
        double bh2x = BH2.getBHCoordX();
        double bh2y = BH2.getBHCoordY();
        double bh2z = BH2.getBHCoordZ();

        double vp2[3];
        vp2[0] = BH2.getVx();
        vp2[1] = BH2.getVy();
        vp2[2] = BH2.getVz();

        double vp2tot = sqrt(vp2[0] * vp2[0] + vp2[1] * vp2[1] + vp2[2] * vp2[2]);
        double spin2 = BH2.getBHSpin();
        double spin2_th = BH2.getBHSpinTheta();
        double spin2_phi = BH2.getBHSpinPhi();

        double vs2[3];
        vs2[0] = spin2 * sin(spin2_th) * cos(spin2_phi);
        vs2[1] = spin2 * sin(spin2_th) * sin(spin2_phi);
        vs2[2] = spin2 * cos(spin2_th);

        // coordinates with respect to center of bh1
        x1 = xx - bh1x;
        y1 = yy - bh1y;
        z1 = zz - bh1z;

        //locating as a radial form
        rv1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
        vn1[0] = x1 / rv1;
        vn1[1] = y1 / rv1;
        vn1[2] = z1 / rv1;

        //same as BH2
        x2 = xx - bh2x;
        y2 = yy - bh2y;
        z2 = zz - bh2z;

        rv2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
        vn2[0] = x2 / rv2;
        vn2[1] = y2 / rv2;
        vn2[2] = z2 / rv2;

        //Initial data is related with the paper: http://arxiv.org/abs/0711.1165
        //Brill-Lindquist conformal factor
        vpsibl = 1.0 + mass1 / (2.0 * rv1);
        vpsibl = vpsibl + mass2 / (2.0 * rv2);

        v_u_corr = 0.0;
        // bh 1

        //For spinning puncture
        if (fabs(spin1) > 1.e-6)
        {
            amp_capj = 4.0 * spin1 / (mass1 * mass1);
            amp_capr = 2.0 * rv1 / mass1;
            l_r = 1.0 / (1.0 + amp_capr);
            u0_j = (l_r + l_r * l_r + l_r * l_r * l_r - 4.0 * l_r * l_r * l_r * l_r + 2.0 * l_r * l_r * l_r * l_r * l_r) / 40.0;
            u2_j = -pow(l_r, 5) / 20.0;
            mu_j = vn1[0] * vs1[0];
            mu_j = mu_j + vn1[1] * vs1[1];
            mu_j = (mu_j + vn1[2] * vs1[2]) / fabs(spin1);
            p2_mu_j = (3.0 * mu_j * mu_j - 1.0) / 2.0;
            v_u_j1 = amp_capj * amp_capj * (u0_j + u2_j * amp_capr * amp_capr * p2_mu_j);
            v_u_corr = v_u_corr + v_u_j1;
        }
        //For boosting puncture
        if (vp1tot > 1.e-6)
        {
            amp_capp = 2.0 * vp1tot / mass1;
            amp_capr = 2.0 * rv1 / mass1;
            l_r = 1.0 / (1.0 + amp_capr);
            u0_p = l_r - 2.0 * l_r * l_r + 2.0 * pow(l_r, 3);
            u0_p = (u0_p - pow(l_r, 4) + 0.20 * pow(l_r, 5)) * (5.0 / 32.0);
            u2_p = 15.0 * l_r + 132.0 * l_r * l_r + 53.0 * pow(l_r, 3);
            u2_p = u2_p + 96.0 * pow(l_r, 4) + 82.0 * pow(l_r, 5);
            u2_p = u2_p + (84.0 / amp_capr) * (pow(l_r, 5) + log(l_r) / amp_capr);
            u2_p = (u2_p) / (80.0 * amp_capr);
            mu_p = vn1[0] * vp1[0] / vp1tot;
            mu_p = mu_p + vn1[1] * vp1[1] / vp1tot;
            mu_p = mu_p + vn1[2] * vp1[2] / vp1tot;
            p2_mu_p = (3.0 * pow(mu_p, 2) - 1.0) / 2.0;
            v_u_p1 = pow(amp_capp, 2) * (u0_p + u2_p * p2_mu_p);
            v_u_corr = v_u_corr + v_u_p1;
        }
        //For spinning boosted pucture
        if (vp1tot > 1.e-6 && fabs(spin1) > 1.e-6)
        {
            v1 = (vp1[1] * vs1[2] - vp1[2] * vs1[1]) * vn1[0];
            v1 = v1 + (vp1[2] * vs1[0] - vp1[0] * vs1[2]) * vn1[1];
            v1 = v1 + (vp1[0] * vs1[1] - vp1[1] * vs1[0]) * vn1[2];
            v1 = v1 * (16.0 / pow(mass1, 4)) * rv1;

            amp_capr = 2.0 * rv1 / mass1;
            l_r = 1.0 / (1.0 + amp_capr);

            v2 = 1.0 + 5.0 * amp_capr + 10.0 * pow(amp_capr, 2);

            v_u_c1 = (v1 * v2 * pow(l_r, 5)) / 80.0;
            v_u_corr = v_u_corr + v_u_c1;
        }
        // bh 2 same puncture as bh 1
        if (fabs(spin2) > 1.e-6)
        {
            amp_capj = 4.0 * spin2 / (mass2 * mass2);
            amp_capr = 2.0 * rv2 / mass2;
            l_r = 1.0 / (1.0 + amp_capr);
            u0_j = (l_r + l_r * l_r + l_r * l_r * l_r - 4.0 * l_r * l_r * l_r * l_r + 2.0 * l_r * l_r * l_r * l_r * l_r) / 40.0;
            u2_j = -pow(l_r, 5) / 20.0;
            mu_j = vn2[0] * vs2[0];
            mu_j = mu_j + vn2[1] * vs2[1];
            mu_j = (mu_j + vn2[2] * vs2[2]) / fabs(spin2);
            p2_mu_j = (3.0 * mu_j * mu_j - 1.0) / 2.0;
            v_u_j2 = amp_capj * amp_capj * (u0_j + u2_j * amp_capr * amp_capr * p2_mu_j);
            v_u_corr = v_u_corr + v_u_j2;
        }

        if (vp2tot > 1.e-6)
        {
            amp_capp = 2.0 * vp2tot / mass2;
            amp_capr = 2.0 * rv2 / mass2;
            l_r = 1.0 / (1.0 + amp_capr);
            u0_p = l_r - 2.0 * l_r * l_r + 2.0 * pow(l_r, 3);
            u0_p = (u0_p - pow(l_r, 4) + 0.20 * pow(l_r, 5)) * (5.0 / 32.0);
            u2_p = 15.0 * l_r + 132.0 * l_r * l_r + 53.0 * pow(l_r, 3);
            u2_p = u2_p + 96.0 * pow(l_r, 4) + 82.0 * pow(l_r, 5);
            u2_p = u2_p + (84.0 / amp_capr) * (pow(l_r, 5) + log(l_r) / amp_capr);
            u2_p = (u2_p) / (80.0 * amp_capr);
            mu_p = vn2[0] * vp2[0] / vp2tot;
            mu_p = mu_p + vn2[1] * vp2[1] / vp2tot;
            mu_p = mu_p + vn2[2] * vp2[2] / vp2tot;
            p2_mu_p = (3.0 * pow(mu_p, 2) - 1.0) / 2.0;
            v_u_p2 = pow(amp_capp, 2) * (u0_p + u2_p * p2_mu_p);
            v_u_corr = v_u_corr + v_u_p2;
        }

        if (vp2tot > 1.e-6 && fabs(spin2) > 1.e-6)
        {
            v1 = (vp2[1] * vs2[2] - vp2[2] * vs2[1]) * vn2[0];
            v1 = v1 + (vp2[2] * vs2[0] - vp2[0] * vs2[2]) * vn2[1];
            v1 = v1 + (vp2[0] * vs2[1] - vp2[1] * vs2[0]) * vn2[2];
            v1 = v1 * (16.0 / pow(mass2, 4)) * rv2;

            amp_capr = 2.0 * rv2 / mass2;
            l_r = 1.0 / (1.0 + amp_capr);

            v2 = 1.0 + 5.0 * amp_capr + 10.0 * pow(amp_capr, 2);

            v_u_c2 = (v1 * v2 * pow(l_r, 5)) / 80.0;
            v_u_corr = v_u_corr + v_u_c2;
        }

        // vpsibl_u will be used for the conformal factor,
        vpsibl_u = vpsibl + v_u_corr;
        // vpsibl_u2 is for the Aij terms...
        // ! since the corrections are first order...
        // ! adding half of the correction seems to give the best results...
        // ! update - do a fit for spin = 0.6...
        vpsibl_u2 = vpsibl + v_u_corr;

        //////// TODO: HUGE: MAKE SURE THIS WORKS!
        // TODO: generate these for EMDA

        var[VAR::U_ALPHA] = 1.0 / (vpsibl_u * vpsibl_u);
        //std::cout<<"Alpha: "<<u[U_ALPHA]<<" vpsibl_u: "<< vpsibl_u<<std::endl;
        var[VAR::U_ALPHA] = std::max(var[VAR::U_ALPHA], PHI_FLOOR);

        v2 = 1.0 / pow(vpsibl_u, 4);
        var[VAR::U_PHI] = v2;

        if (var[VAR::U_PHI] < PHI_FLOOR)
            var[VAR::U_PHI] = PHI_FLOOR;

        var[VAR::U_K] = 0.0;

        var[VAR::U_BETA0] = 0.0;
        var[VAR::U_BETA1] = 0.0;
        var[VAR::U_BETA2] = 0.0;

        var[VAR::U_GAMMAHAT0] = 0.0;
        var[VAR::U_GAMMAHAT1] = 0.0;
        var[VAR::U_GAMMAHAT2] = 0.0;

        var[VAR::U_B0] = 0.0;
        var[VAR::U_B1] = 0.0;
        var[VAR::U_B2] = 0.0;

        var[VAR::U_GAMMAT00] = 1.0; //XX
        var[VAR::U_GAMMAT01] = 0.0; //XY
        var[VAR::U_GAMMAT02] = 0.0; //XZ
        var[VAR::U_GAMMAT11] = 1.0; //YY
        var[VAR::U_GAMMAT12] = 0.0; //YZ
        var[VAR::U_GAMMAT22] = 1.0; //ZZ

        for (i1 = 0; i1 < 3; i1++)
        {
            for (i2 = 0; i2 < 3; i2++)
            {
                // first BH
                v2 = 0.0;
                for (i3 = 0; i3 < 3; i3++)
                {
                    for (i4 = 0; i4 < 3; i4++)
                    {
                        vt1 = epijk[i1][i3][i4] * vs1[i3] * vn1[i4] * vn1[i2];
                        vt2 = epijk[i2][i3][i4] * vs1[i3] * vn1[i4] * vn1[i1];
                        v2 = v2 + vt1 + vt2;
                    }
                }

                v3 = vp1[i1] * vn1[i2] + vp1[i2] * vn1[i1];
                vt1 = 0.0;
                for (i3 = 0; i3 < 3; i3++)
                {
                    vt1 = vt1 + vp1[i3] * vn1[i3];
                }
                vt1 = vt1 * (vn1[i1] * vn1[i2] - deltaij[i1][i2]);
                v3 = v3 + vt1;

                v1 = 3.0 / (pow(vpsibl_u2, 6) * pow(rv1, 3));
                v4 = v1 * (v2 + (rv1 / 2.0) * v3);

                // second BH
                v2 = 0.0;
                for (i3 = 0; i3 < 3; i3++)
                {
                    for (i4 = 0; i4 < 3; i4++)
                    {
                        vt1 = epijk[i1][i3][i4] * vs2[i3] * vn2[i4] * vn2[i2];
                        vt2 = epijk[i2][i3][i4] * vs2[i3] * vn2[i4] * vn2[i1];
                        v2 = v2 + vt1 + vt2;
                    }
                }

                v3 = vp2[i1] * vn2[i2] + vp2[i2] * vn2[i1];
                vt1 = 0.0;
                for (i3 = 0; i3 < 3; i3++)
                {
                    vt1 = vt1 + vp2[i3] * vn2[i3];
                }
                vt1 = vt1 * (vn2[i1] * vn2[i2] - deltaij[i1][i2]);
                v3 = v3 + vt1;

                v1 = 3.0 / (pow(vpsibl_u2, 6) * pow(rv2, 3));
                v4 = v4 + v1 * (v2 + (rv2 / 2.0) * v3);

                if (i1 == 0 && i2 == 0)
                {
                    var[VAR::U_AT00] = v4; //XX
                }
                else if (i1 == 0 && i2 == 1)
                {
                    var[VAR::U_AT01] = v4; //XY
                }
                else if (i1 == 0 && i2 == 2)
                {
                    var[VAR::U_AT02] = v4; //XZ
                }
                else if (i1 == 1 && i2 == 1)
                {
                    var[VAR::U_AT11] = v4; //YY
                }
                else if (i1 == 1 && i2 == 2)
                {
                    var[VAR::U_AT12] = v4; //YZ
                }
                else if (i1 == 2 && i2 == 2)
                {
                    var[VAR::U_AT22] = v4; //ZZ
                }
            }
        }
#endif
}

void punctureData(const double xx1, const double yy1, const double zz1,
                  double *var) {
    const double xx = GRIDX_TO_X(xx1);
    const double yy = GRIDY_TO_Y(yy1);
    const double zz = GRIDZ_TO_Z(zz1);

    punctureDataPhysicalCoord(xx, yy, zz, var);
}

// TODO: remove? trumpet data is not called in the BSSN code anywhere (it is
// stored in that repo)
#if 0
    namespace trumpet_data
    {

        const unsigned int KMAX = 5000;
        double dxsav;
        static double *xp;
        static double *yp[2];

        void trumpetData(const double xx1, const double yy1, const double zz1, double *var)
        {

            const double xx = GRIDX_TO_X(xx1);
            const double yy = GRIDY_TO_Y(yy1);
            const double zz = GRIDZ_TO_Z(zz1);

            double eps = 1.e-10;  // tolerance on the integrator
            double h1 = -1.e-5;   // starting step off the critical point
            double hmin = 1.e-25; // min allowed stepsize for adaptive integrator

            double bh_mass = 7.5;
            double alpha_c = 0.16227766016837933200;
            double C_sq = 1.5543095902183302040;
            double bigR_c = 1.5405694150420948330 * bh_mass;
            double r_c = 0.30405997036 * bh_mass;

            static bool firstcall = true;

            const int neq = 2;
            static double *alpha0;
            static double *bigR0;
            static double *r_iso;

            static int np;

            const double third = 1.0 / 3.0;

            if (firstcall == true)
            {
                // solve ODE system
                std::cout << "trumpetData:  first call. Solve the ODE" << std::endl;

                xp = new double[KMAX];
                yp[0] = new double[KMAX];
                yp[1] = new double[KMAX];
                alpha0 = new double[KMAX];
                bigR0 = new double[KMAX];
                r_iso = new double[KMAX];

                // integrate inward from r_C to r=0.
                double rmin = 0.0;    // min value for inward r integration
                double rmax = 1000.0; // max value for outward r integration
                double rstart, ystart[2], dydr[2];
                bndcnd(h1, rstart, ystart, dydr);
                double rend = rmin + fabs(h1);
                int nok, nbad;
                int kount;
                odeint(ystart, neq, rstart, rend, eps, h1, hmin, &nok, &nbad, derivs, rkqs, kount);

                int kountin = kount;
                for (unsigned int i = 0; i < kountin; i++)
                {
                    r_iso[i] = xp[kountin - 1 - i];
                    alpha0[i] = yp[0][kountin - 1 - i];
                    bigR0[i] = yp[1][kountin - 1 - i];
                    //std::cout<<"<in> xp = "<<xp[i]<<", i="<<i<<", alpha0 = "<<yp[0][i]<<", bigR0 = "<<yp[1][i]<<std::endl;
                }

                // integrate outwards from r_c to large r
                h1 = fabs(h1);
                rend = rmax;

                std::cout << "integrating outwards ... " << std::endl;
                bndcnd(h1, rstart, ystart, dydr);
                odeint(ystart, neq, rstart, rend, eps, h1, hmin, &nok, &nbad, derivs, rkqs, kount);

                int kountout = kount;
                for (unsigned int i = 0; i < kountout; i++)
                {
                    r_iso[i + kountin] = xp[i];
                    alpha0[i + kountin] = yp[0][i];
                    bigR0[i + kountin] = yp[1][i];
                    //std::cout<<"<out> xp = "<<xp[i]<<", alpha0 = "<<yp[0][i]<<", bigR0 = "<<yp[1][i]<<std::endl;
                }
                np = kountin + kountout;

                firstcall = false;
            }

            int nrby_indx = np / 2;
            double ax_eps = 1.0e-5;

            double tenh1 = 10.0 * fabs(h1);
            double rbar = sqrt(xx * xx + yy * yy + zz * zz);
            if (fabs(xx) < tenh1 && fabs(yy) < tenh1 && fabs(zz) < tenh1)
            {
                rbar = tenh1;
            }
            double alpha = interpolation4(r_iso, alpha0, np, rbar, &nrby_indx);
            double bigR = interpolation4(r_iso, bigR0, np, rbar, &nrby_indx);

            if (fabs(xx) < ax_eps && fabs(yy) < ax_eps && fabs(zz) < ax_eps)
            {
                rbar = sqrt(xx * xx + yy * yy + zz * zz + ax_eps * ax_eps);
            }

            double f0 = 1.0 - 2.0 * bh_mass / bigR;

            double Rsq_dalpha_dR = 4.0 * (alpha * alpha - f0 - 0.5 * bh_mass / bigR) / (alpha * alpha - 2.0 * alpha - f0) * bigR;

            double tmp_sqrt = sqrt(alpha * alpha - f0);

            double tmp_chi = (rbar * rbar) / (bigR * bigR);

            double tmp_trK = 2.0 * tmp_sqrt / bigR + (alpha * Rsq_dalpha_dR - bh_mass) / (tmp_sqrt * bigR * bigR);

            double tmp_beta = tmp_sqrt / bigR;

            double tmp_Atilde = ((alpha * Rsq_dalpha_dR - bh_mass) / tmp_sqrt - bigR * tmp_sqrt) / bigR / bigR;

            // TODO: fix this for EMDA?
            var[VAR::U_ALPHA] = alpha;
            var[VAR::U_ALPHA] = std::max(var[VAR::U_ALPHA], PHI_FLOOR);

            var[VAR::U_PHI] = tmp_chi;

            if (var[VAR::U_PHI] < PHI_FLOOR)
                var[VAR::U_PHI] = PHI_FLOOR;

            var[VAR::U_K] = tmp_trK;

            var[VAR::U_BETA0] = xx * tmp_beta;
            var[VAR::U_BETA1] = yy * tmp_beta;
            var[VAR::U_BETA2] = zz * tmp_beta;

            var[VAR::U_GAMMAHAT0] = 0.0;
            var[VAR::U_GAMMAHAT1] = 0.0;
            var[VAR::U_GAMMAHAT2] = 0.0;

            var[VAR::U_B0] = 0.0;
            var[VAR::U_B1] = 0.0;
            var[VAR::U_B2] = 0.0;

            var[VAR::U_GAMMAT00] = 1.0; //XX
            var[VAR::U_GAMMAT01] = 0.0; //XY
            var[VAR::U_GAMMAT02] = 0.0; //XZ
            var[VAR::U_GAMMAT11] = 1.0; //YY
            var[VAR::U_GAMMAT12] = 0.0; //YZ
            var[VAR::U_GAMMAT22] = 1.0; //ZZ

            var[VAR::U_AT00] = tmp_Atilde * ((xx / rbar) * (xx / rbar) - third);
            var[VAR::U_AT01] = tmp_Atilde * xx * yy / (rbar * rbar);
            var[VAR::U_AT02] = tmp_Atilde * xx * zz / (rbar * rbar);
            var[VAR::U_AT11] = tmp_Atilde * ((yy / rbar) * (yy / rbar) - third);
            var[VAR::U_AT12] = tmp_Atilde * yy * zz / (rbar * rbar);
            var[VAR::U_AT22] = tmp_Atilde * ((zz / rbar) * (zz / rbar) - third);
        }

        /*---------------------------------------------------------------------
    *
    *
    *
    *---------------------------------------------------------------------*/
        void bndcnd(double h, double &x, double y[], double dydx[])
        {

            double bh_mass = 7.5;
            double alpha_c = 0.16227766016837933200;
            double bigR_c = 1.5405694150420948330 * bh_mass;
            double r_c = 0.30405997036 * bh_mass;

            double alpha = alpha_c;
            double bigR = bigR_c;
            double r = r_c;

            double bigRsq = bigR * bigR;
            double alphasq = alpha * alpha;

            double dbigR_dr = alpha_c * bigR_c / r_c;

            double df0_dr = 2.0 * bh_mass / bigRsq * dbigR_dr;

            double tmp1 = 1.5 * bh_mass / bigRsq * dbigR_dr;

            double dalpha_dr = 0.25 * (r * df0_dr + 8.0 * alphasq) / (alpha - 1.0) - 0.25 * sqrt(pow((r * df0_dr + 8.0 * alphasq), 2) + 32.0 * alpha * (1.0 - alpha) * r * tmp1) / (alpha - 1.0);

            dalpha_dr = dalpha_dr / r;

            double ddbigR_drdr = (dalpha_dr * bigR + alpha * dbigR_dr) / r - alpha * bigR / (r * r);

            double ddf0_drdr = -4.0 * bh_mass * dbigR_dr * dbigR_dr / (bigR * bigRsq) + 2.0 * bh_mass * ddbigR_drdr / bigRsq;

            double tmp2 = -3.00 * bh_mass * dbigR_dr * dbigR_dr / (bigR * bigRsq) + 1.50 * bh_mass * ddbigR_drdr / bigRsq;

            double ddalpha_drdr = -2.0 * pow((r * dalpha_dr), 3) - 3.0 * r * dalpha_dr * (2.0 * r * dalpha_dr * (alpha - 1.0) - r * df0_dr) + (r * dalpha_dr) * r * r * ddf0_drdr + (8.0 * r * dalpha_dr + 4.0 * alpha) * (2.0 * alpha * r * dalpha_dr - r * tmp1) + 4.0 * alpha * (2.0 * pow((r * dalpha_dr), 2) - r * r * tmp2);

            ddalpha_drdr = ddalpha_drdr / (6.0 * pow(r, 3) * (alpha - 1.0) * dalpha_dr - 2.0 * pow(r, 3) * df0_dr - 8.0 * r * r * alphasq);

            x = r + h;

            y[0] = alpha + h * dalpha_dr + 0.5 * h * h * ddalpha_drdr;
            y[1] = bigR + h * dbigR_dr + 0.5 * h * h * ddbigR_drdr;

            dydx[0] = dalpha_dr + h * ddalpha_drdr;
            dydx[1] = dbigR_dr + h * ddbigR_drdr;
        }

        /*---------------------------------------------------------------------
    *
    *
    *
    *---------------------------------------------------------------------*/
        void derivs(double x, double y[], double dydx[])
        {

            double alpha = y[0];
            double bigR = y[1];
            double r = x;

            double bh_mass = 7.5;

            double alpha_sq = alpha * alpha;

            double alpha_sq_minus_1 = alpha_sq - 1.0;
            double M_over_R = bh_mass / bigR;

            double dalpha_dr = 4.0 * alpha * (alpha_sq_minus_1 + 1.5 * M_over_R) / (alpha_sq_minus_1 - 2.0 * alpha + 2.0 * M_over_R) / r;

            double dbigR_dr = alpha * bigR / r;

            dydx[0] = dalpha_dr;
            dydx[1] = dbigR_dr;
        }

        /*---------------------------------------------------------------------
    *
    *
    *
    *---------------------------------------------------------------------*/
        void hunt(double xx[], int n, double x, int *jlo)
        {
            unsigned long jm, jhi, inc;
            int ascnd;

            ascnd = (xx[n - 1] > xx[0]);
            //if (*jlo < 0 || *jlo > n)
            if (*jlo < 0 || *jlo > n - 1)
            {
                //*jlo=0;
                *jlo = -1;
                //jhi=n-1;
                jhi = n;
            }
            else
            {
                inc = 1;
                if (x >= xx[*jlo] == ascnd)
                {
                    if (*jlo == n - 1)
                        return;
                    jhi = (*jlo) + 1;
                    while (x >= xx[jhi] == ascnd)
                    {
                        *jlo = jhi;
                        inc += inc;
                        jhi = (*jlo) + inc;
                        if (jhi > n - 1)
                        {
                            jhi = n;
                            break;
                        }
                    }
                }
                else
                {
                    if (*jlo == 0)
                    {
                        *jlo = -1;
                        return;
                    }
                    jhi = (*jlo)--;
                    while (x < xx[*jlo] == ascnd)
                    {
                        jhi = (*jlo);
                        inc <<= 1;
                        if (inc >= jhi)
                        {
                            *jlo = -1;
                            break;
                        }
                        else
                            *jlo = jhi - inc;
                    }
                }
            }
            while (jhi - (*jlo) != 1)
            {
                jm = (jhi + (*jlo)) >> 1;
                if (x > xx[jm] == ascnd)
                    *jlo = jm;
                else
                    jhi = jm;
            }
        } /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

        /*---------------------------------------------------------------------
    *
    *
    *
    *---------------------------------------------------------------------*/
        void rkck(double y[], double dydx[], int n, double x, double h,
                  double yout[], double yerr[],
                  void (*derivs)(double, double[], double[]))
        {
            int i;
            static double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875, b21 = 0.2,
                          b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9, b43 = 1.2,
                          b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0, b54 = 35.0 / 27.0,
                          b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0,
                          b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0, c1 = 37.0 / 378.0,
                          c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
                          dc5 = -277.0 / 14336.0;
            double dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
                   dc4 = c4 - 13525.0 / 55296.0, dc6 = c6 - 0.25;

            double *ak2 = new double[n];
            double *ak3 = new double[n];
            double *ak4 = new double[n];
            double *ak5 = new double[n];
            double *ak6 = new double[n];
            double *ytemp = new double[n];

            for (i = 0; i < n; i++)
                ytemp[i] = y[i] + b21 * h * dydx[i];
            (*derivs)(x + a2 * h, ytemp, ak2);
            for (i = 0; i < n; i++)
                ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
            (*derivs)(x + a3 * h, ytemp, ak3);
            for (i = 0; i < n; i++)
                ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
            (*derivs)(x + a4 * h, ytemp, ak4);
            for (i = 0; i < n; i++)
                ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
            (*derivs)(x + a5 * h, ytemp, ak5);
            for (i = 0; i < n; i++)
                ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
            (*derivs)(x + a6 * h, ytemp, ak6);
            for (i = 0; i < n; i++)
                yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
            for (i = 0; i < n; i++)
                yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
            delete[] ytemp;
            delete[] ak6;
            delete[] ak5;
            delete[] ak4;
            delete[] ak3;
            delete[] ak2;
        } /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

        /*---------------------------------------------------------------------
    *
    *
    *
    *---------------------------------------------------------------------*/
        void rkqs(double y[], double dydx[], int n, double *x, double htry,
                  double eps, double yscal[], double *hdid, double *hnext,
                  void (*derivs)(double, double[], double[]))
        {
            const double SAFETY = 0.9;
            const double PGROW = -0.2;
            const double PSHRNK = -0.25;
            const double ERRCON = 1.89e-4;

            void rkck(double y[], double dydx[], int n, double x, double h,
                      double yout[], double yerr[],
                      void (*derivs)(double, double[], double[]));
            int i;
            double errmax, h, xnew;

            double *yerr = new double[n];
            double *ytemp = new double[n];
            h = htry;
            for (;;)
            {
                rkck(y, dydx, n, *x, h, ytemp, yerr, derivs);
                errmax = 0.0;
                for (i = 0; i < n; i++)
                    errmax = std::max(errmax, fabs(yerr[i] / yscal[i]));
                errmax /= eps;
                if (errmax > 1.0)
                {
                    h = SAFETY * h * pow(errmax, PSHRNK);
                    if (h < 0.1 * h)
                        h *= 0.1;
                    xnew = (*x) + h;
                    if (xnew == *x)
                        std::cerr << "stepsize underflow in rkqs" << std::endl;
                    continue;
                }
                else
                {
                    if (errmax > ERRCON)
                        *hnext = SAFETY * h * pow(errmax, PGROW);
                    else
                        *hnext = 5.0 * h;
                    *x += (*hdid = h);
                    for (i = 0; i < n; i++)
                        y[i] = ytemp[i];
                    break;
                }
            }
            delete[] ytemp;
            delete[] yerr;
        } /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

        /*---------------------------------------------------------------------
    *
    *
    *
    *---------------------------------------------------------------------*/
        void odeint(double ystart[], int nvar, double x1, double x2,
                    double eps, double h1, double hmin, int *nok, int *nbad,
                    void (*derivs)(double, double[], double[]),
                    void (*rkqs)(double[], double[], int, double *,
                                 double, double, double[], double *,
                                 double *,
                                 void (*)(double, double[], double[])),
                    int kount)
        {
            int nstp, i;
            double xsav, x, hnext, hdid, h;
            const int MAXSTP = 10000;
            const double TINY = 1.0e-30;

            double *yscal = new double[nvar];
            double *y = new double[nvar];
            double *dydx = new double[nvar];
            x = x1;
            h = copysign(h1, x2 - x1);
            *nok = (*nbad) = kount = 0;
            for (i = 0; i < nvar; i++)
                y[i] = ystart[i];
            if (KMAX > 0)
                xsav = x - dxsav * 2.0;
            for (nstp = 0; nstp < MAXSTP; nstp++)
            {
                //std::cout<<"odeint: nstp="<<nstp<<", kount="<<kount<<std::endl;
                (*derivs)(x, y, dydx);
                for (i = 0; i < nvar; i++)
                    yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
                if (KMAX > 0 && kount < KMAX - 1 && fabs(x - xsav) > fabs(dxsav))
                {
                    xp[kount] = x;
                    for (i = 0; i < nvar; i++)
                        yp[i][kount] = y[i];
                    xsav = x;
                    kount++;
                }
                if ((x + h - x2) * (x + h - x1) > 0.0)
                    h = x2 - x;
                (*rkqs)(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);
                if (hdid == h)
                    ++(nok);
                else
                    ++(nbad);
                if ((x - x2) * (x2 - x1) >= 0.0)
                {
                    for (i = 0; i < nvar; i++)
                        ystart[i] = y[i];
                    if (KMAX)
                    {
                        xp[kount] = x;
                        for (i = 0; i < nvar; i++)
                            yp[i][kount] = y[i];
                        kount++;
                    }
                    delete[] dydx;
                    delete[] y;
                    delete[] yscal;
                    return;
                }
                if (fabs(hnext) <= hmin)
                    std::cerr << "Step size too small in odeint" << std::endl;
                h = hnext;
            }
            std::cerr << "Too many steps in routine odeint" << std::endl;
        }
        /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

        /*---------------------------------------------------------------------
    *
    *
    *
    *---------------------------------------------------------------------*/
        double interpolation3(double xp[], double yp[], int np, double xb,
                              int *n_nearest_pt)
        {
            int k;     /* index of 1st point */
            int m = 4; /* degree of interpolation */
            double y;  /* intermediate value */

            hunt(xp, np, xb, n_nearest_pt);

            k = std::min(std::max((*n_nearest_pt) - (m - 1) / 2,
                                  1),
                         np + 1 - m);

            double DBL_EPSILON = 1.e-12;

            if (xb == xp[k] || xb == xp[k + 1] || xb == xp[k + 2] || xb == xp[k + 3])
            {
                xb += 3.0 * DBL_EPSILON;
            }

            y = (xb - xp[k + 1]) * (xb - xp[k + 2]) * (xb - xp[k + 3]) * yp[k] / ((xp[k] - xp[k + 1]) * (xp[k] - xp[k + 2]) * (xp[k] - xp[k + 3]))

                + (xb - xp[k]) * (xb - xp[k + 2]) * (xb - xp[k + 3]) * yp[k + 1] / ((xp[k + 1] - xp[k]) * (xp[k + 1] - xp[k + 2]) * (xp[k + 1] - xp[k + 3]))

                + (xb - xp[k]) * (xb - xp[k + 1]) * (xb - xp[k + 3]) * yp[k + 2] / ((xp[k + 2] - xp[k]) * (xp[k + 2] - xp[k + 1]) * (xp[k + 2] - xp[k + 3]))

                + (xb - xp[k]) * (xb - xp[k + 1]) * (xb - xp[k + 2]) * yp[k + 3] / ((xp[k + 3] - xp[k]) * (xp[k + 3] - xp[k + 1]) * (xp[k + 3] - xp[k + 2]));

            return (y);
        }

        /*---------------------------------------------------------------------
    *
    *
    *
    *---------------------------------------------------------------------*/
        double interpolation4(double xx[], double yy[], int np, double xb,
                              int *n_nearest_pt)
        {
            int k;     /* index of 1st point */
            int m = 5; /* degree of interpolation */
            double y;  /* intermediate value */

            hunt(xx, np, xb, n_nearest_pt);

            k = std::min(std::max((*n_nearest_pt) - (m - 1) / 2,
                                  1),
                         np + 1 - m);

#if 0
        if ( fabs( xb - 1.0 ) < 1.e-13 )
        { // xb is 1.0 so just return with the corresponding y value -- no interp necessary
            y = yy[np] ;
            return(y) ;
        }
        if ( fabs( xb ) < 1.e-13 )
        { // xb is zero so just return with the corresponding y value -- no interp necessary
            y = yy[0] ;
            return(y) ;
        }
#endif

            double DBL_EPSILON = 1.e-12;

            if (xb == xx[k] || xb == xx[k + 1] || xb == xx[k + 2] || xb == xx[k + 3] || xb == xx[k + 4])
            {
                xb += 3.0 * DBL_EPSILON;
            }

            double xtmp0 = xb - xx[k];
            double xtmp1 = xb - xx[k + 1];
            double xtmp2 = xb - xx[k + 2];
            double xtmp3 = xb - xx[k + 3];
            double xtmp4 = xb - xx[k + 4];
            double xdiff01 = xx[k] - xx[k + 1];
            double xdiff02 = xx[k] - xx[k + 2];
            double xdiff03 = xx[k] - xx[k + 3];
            double xdiff04 = xx[k] - xx[k + 4];
            double xdiff12 = xx[k + 1] - xx[k + 2];
            double xdiff13 = xx[k + 1] - xx[k + 3];
            double xdiff14 = xx[k + 1] - xx[k + 4];
            double xdiff23 = xx[k + 2] - xx[k + 3];
            double xdiff24 = xx[k + 2] - xx[k + 4];
            double xdiff34 = xx[k + 3] - xx[k + 4];

            y = xtmp1 * xtmp2 * xtmp3 * xtmp4 * yy[k] / (xdiff01 * xdiff02 * xdiff03 * xdiff04)

                - xtmp0 * xtmp2 * xtmp3 * xtmp4 * yy[k + 1] / (xdiff01 * xdiff12 * xdiff13 * xdiff14)

                + xtmp0 * xtmp1 * xtmp3 * xtmp4 * yy[k + 2] / (xdiff02 * xdiff12 * xdiff23 * xdiff24)

                - xtmp0 * xtmp1 * xtmp2 * xtmp4 * yy[k + 3] / (xdiff03 * xdiff13 * xdiff23 * xdiff34) + xtmp0 * xtmp1 * xtmp2 * xtmp3 * yy[k + 4] / (xdiff04 * xdiff14 * xdiff24 * xdiff34);

            return (y);
        }

    } // end of namespace trumpet_data

#endif

void blockAdaptiveOctree(std::vector<ot::TreeNode> &tmpNodes,
                         const Point &pt_min, const Point &pt_max,
                         const unsigned int regLev, const unsigned int maxDepth,
                         MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    double pt_g_min[3];
    double pt_g_max[3];

    pt_g_min[0] = X_TO_GRIDX(pt_min.x());
    pt_g_min[1] = Y_TO_GRIDY(pt_min.y());
    pt_g_min[2] = Z_TO_GRIDZ(pt_min.z());

    pt_g_max[0] = X_TO_GRIDX(pt_max.x());
    pt_g_max[1] = Y_TO_GRIDY(pt_max.y());
    pt_g_max[2] = Z_TO_GRIDZ(pt_max.z());

    std::cout << pt_g_min[0] << " " << maxDepth << " " << (1u << maxDepth)
              << std::endl;
    std::cout << pt_g_max[0] << " " << maxDepth << " " << (1u << maxDepth)
              << std::endl;

    assert(pt_g_min[0] >= 0 && pt_g_min[0] <= (1u << maxDepth));
    assert(pt_g_min[1] >= 0 && pt_g_min[1] <= (1u << maxDepth));
    assert(pt_g_min[2] >= 0 && pt_g_min[2] <= (1u << maxDepth));

    assert(pt_g_max[0] >= 0 && pt_g_max[0] <= (1u << maxDepth));
    assert(pt_g_max[1] >= 0 && pt_g_max[1] <= (1u << maxDepth));
    assert(pt_g_max[2] >= 0 && pt_g_max[2] <= (1u << maxDepth));

    unsigned int xRange_b, xRange_e;
    unsigned int yRange_b = pt_g_min[1], yRange_e = pt_g_max[1];
    unsigned int zRange_b = pt_g_min[2], zRange_e = pt_g_max[2];

    xRange_b =
        pt_g_min[0];  //(rank*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];
    xRange_e =
        pt_g_max[1];  //((rank+1)*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];

    unsigned int stepSz = 1u << (maxDepth - regLev);

    /* std::cout<<" x min: "<<xRange_b<<" x_max: "<<xRange_e<<std::endl;
    std::cout<<" y min: "<<yRange_b<<" y_max: "<<yRange_e<<std::endl;
    std::cout<<" z min: "<<zRange_b<<" z_max: "<<zRange_e<<std::endl;*/

    std::cout << "Now adjusting tmpNodes" << std::endl;

    for (unsigned int x = xRange_b; x < xRange_e; x += stepSz)
        for (unsigned int y = yRange_b; y < yRange_e; y += stepSz)
            for (unsigned int z = zRange_b; z < zRange_e; z += stepSz) {
                if (x >= (1u << maxDepth)) x = x - 1;
                if (y >= (1u << maxDepth)) y = y - 1;
                if (z >= (1u << maxDepth)) z = z - 1;

                tmpNodes.push_back(
                    ot::TreeNode(x, y, z, regLev, m_uiDim, maxDepth));
            }

    return;
}

double computeWTol(double x, double y, double z, double tolMin) {
    double origin[3];
    origin[0] = (double)(1u << dsolve::DENDROSOLVER_MAXDEPTH - 1);
    origin[1] = (double)(1u << dsolve::DENDROSOLVER_MAXDEPTH - 1);
    origin[2] = (double)(1u << dsolve::DENDROSOLVER_MAXDEPTH - 1);

    double r =
        sqrt(GRIDX_TO_X(x) * GRIDX_TO_X(x) + GRIDY_TO_Y(y) * GRIDY_TO_Y(y) +
             GRIDZ_TO_Z(z) * GRIDZ_TO_Z(z));

    const double tolMax = dsolve::DENDROSOLVER_WAVELET_TOL_MAX;
    const double R0 = dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R0;
    const double R1 = dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R1;

    if (dsolve::DENDROSOLVER_USE_WAVELET_TOL_FUNCTION == 1) {
        return std::min(
            tolMax, std::max(tolMin, (tolMax - tolMin) / (R1 - R0) * (r - R0) +
                                         tolMin));
    } else {
        return tolMin;
    }
}

double computeWTolDCoords(double x, double y, double z, double *hx) {
    const double tolMax = dsolve::DENDROSOLVER_WAVELET_TOL_MAX;
    const double tolMin = dsolve::DENDROSOLVER_WAVELET_TOL;
    const unsigned int eleOrder = dsolve::DENDROSOLVER_ELE_ORDER;

    if (dsolve::DENDROSOLVER_USE_WAVELET_TOL_FUNCTION == 1) {
        const double R0 = dsolve::DENDROSOLVER_BH1_AMR_R;
        const double R1 = dsolve::DENDROSOLVER_BH2_AMR_R;
        const double dbh =
            (dsolve::DENDROSOLVER_BH_LOC[0] - dsolve::DENDROSOLVER_BH_LOC[1])
                .abs();

        // R_Max is defined based on the initial separation.
        const double R_MAX =
            (dsolve::BH1.getBHCoord() - dsolve::BH2.getBHCoord()).abs() + R0 +
            R1;

        Point grid_p(x, y, z);
        const double dbh0 = (grid_p - dsolve::DENDROSOLVER_BH_LOC[0]).abs();
        const double dbh1 = (grid_p - dsolve::DENDROSOLVER_BH_LOC[1]).abs();

#ifdef DENDROSOLVER_EXTRACT_GRAVITATIONAL_WAVES
        if (dbh < 0.1) {
            if ((dbh0 > R_MAX) && (dbh1 > R_MAX)) {
                const double dr = sqrt(x * x + y * y + z * z);
                if (dr <
                    (GW::DENDROSOLVER_GW_RADAII[GW::DENDROSOLVER_GW_NUM_RADAII -
                                                1] +
                     10))
                    return DENDROSOLVER_GW_REFINE_WTOL;
                else
                    return tolMax;
            } else {
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                            const double xx = x + i * hx[0];
                            const double yy = y + j * hx[1];
                            const double zz = z + k * hx[2];

                            const Point grid_pp(xx, yy, zz);

                            const double dd0 =
                                (grid_pp - dsolve::DENDROSOLVER_BH_LOC[0])
                                    .abs();
                            const double dd1 =
                                (grid_pp - dsolve::DENDROSOLVER_BH_LOC[1])
                                    .abs();

                            // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                            // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                            // "<<hx[0]<<std::endl;

                            if (dd0 < R0 || dd1 < R1) {
                                // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                                // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                                // "<<hx[0]<<std::endl;
                                return tolMin;
                            }
                        }

                const double dr = sqrt(x * x + y * y + z * z);
                if (dr <
                    (GW::DENDROSOLVER_GW_RADAII[GW::DENDROSOLVER_GW_NUM_RADAII -
                                                1] +
                     10))
                    return DENDROSOLVER_GW_REFINE_WTOL;
                else
                    return tolMax;
            }
        } else {
            if ((dbh0 > R_MAX) &&
                (dbh1 >
                 R_MAX))  // no need to check individual points in the element
                return tolMax;
            else {
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                            const double xx = x + i * hx[0];
                            const double yy = y + j * hx[1];
                            const double zz = z + k * hx[2];

                            const Point grid_pp(xx, yy, zz);

                            const double dd0 =
                                (grid_pp - dsolve::DENDROSOLVER_BH_LOC[0])
                                    .abs();
                            const double dd1 =
                                (grid_pp - dsolve::DENDROSOLVER_BH_LOC[1])
                                    .abs();

                            // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                            // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                            // "<<hx[0]<<std::endl;

                            if (dd0 < R0 || dd1 < R1) {
                                // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                                // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                                // "<<hx[0]<<std::endl;
                                return tolMin;
                            }
                        }

                //@milinda 21/11/2020 - smooth transition of the wtol.
                if (dbh0 < dbh1)
                    return std::min(
                        tolMax,
                        std::max(tolMin, (tolMax - tolMin) / (R_MAX - R0) *
                                                 (dbh0 - R0) +
                                             tolMin));
                else
                    return std::min(
                        tolMax,
                        std::max(tolMin, (tolMax - tolMin) / (R_MAX - R1) *
                                                 (dbh1 - R1) +
                                             tolMin));
            }
        }

#else
        if ((dbh0 > R_MAX) &&
            (dbh1 >
             R_MAX))  // no need to check individual points in the element
            return tolMax;
        else {
            for (unsigned int k = 0; k < (eleOrder + 1); k++)
                for (unsigned int j = 0; j < (eleOrder + 1); j++)
                    for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                        const double xx = x + i * hx[0];
                        const double yy = y + j * hx[1];
                        const double zz = z + k * hx[2];

                        const Point grid_pp(xx, yy, zz);

                        const double dd0 =
                            (grid_pp - dsolve::DENDROSOLVER_BH_LOC[0]).abs();
                        const double dd1 =
                            (grid_pp - dsolve::DENDROSOLVER_BH_LOC[1]).abs();

                        // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<" dd0:
                        // "<<dd0<<" dd1: "<<dd1<<" hx: "<<hx[0]<<std::endl;

                        if (dd0 < R0 || dd1 < R1) {
                            // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                            // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                            // "<<hx[0]<<std::endl;
                            return tolMin;
                        }
                    }

            //@milinda 21/11/2020 - smooth transition of the wtol.
            if (dbh0 < dbh1)
                return std::min(
                    tolMax, std::max(tolMin, (tolMax - tolMin) / (R_MAX - R0) *
                                                     (dbh0 - R0) +
                                                 tolMin));
            else
                return std::min(
                    tolMax, std::max(tolMin, (tolMax - tolMin) / (R_MAX - R1) *
                                                     (dbh1 - R1) +
                                                 tolMin));
        }
#endif
    } else if (dsolve::DENDROSOLVER_USE_WAVELET_TOL_FUNCTION == 2) {
        const double r = sqrt(x * x + y * y + z * z);
        const double R0 = dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R0;
        const double R1 = dsolve::DENDROSOLVER_WAVELET_TOL_FUNCTION_R1;
        return std::min(
            tolMax, std::max(tolMin, (tolMax - tolMin) / (R1 - R0) * (r - R0) +
                                         tolMin));
    } else {
        return tolMin;
    }
}

void writeBLockToBinary(const double **unzipVarsRHS, unsigned int offset,
                        const double *pmin, const double *pmax, double *bxMin,
                        double *bxMax, const unsigned int *sz,
                        unsigned int blkSz, double dxFactor,
                        const char *fprefix) {
    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const unsigned int ib = 3;
    const unsigned int jb = 3;
    const unsigned int kb = 3;

    const unsigned int ie = nx - 3;
    const unsigned int je = ny - 3;
    const unsigned int ke = nz - 3;

    const unsigned int blkInlSz = (nx - 3) * (ny - 3) * (nz - 3);

    double hx = (pmax[0] - pmin[0]) / (nx - 1);
    double hy = (pmax[1] - pmin[1]) / (ny - 1);
    double hz = (pmax[2] - pmin[2]) / (nz - 1);

    const double dx = (dsolve::DENDROSOLVER_COMPD_MAX[0] -
                       dsolve::DENDROSOLVER_COMPD_MIN[0]) *
                      (1.0 / (double)(1u << dsolve::DENDROSOLVER_MAXDEPTH));
    unsigned int level =
        dsolve::DENDROSOLVER_MAXDEPTH - ((unsigned int)(hx / dx) - 1);

    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);
    // std::cout<<"ranl: "<<rank<<"npes: "<<npes<<std::endl;

    // std::cout<<"nx: "<<nx<<" level: "<<level<<" hx: "<<hx<<" dx:
    // "<<dx<<std::endl;

    if ((hx > (dxFactor * dx)) ||
        (pmin[0] < bxMin[0] || pmin[1] < bxMin[1] || pmin[2] < bxMin[2]) ||
        (pmax[0] > bxMax[0] || pmax[1] > bxMax[1] || pmax[2] > bxMax[2]))
        return;

    double *blkInternal = new double[blkInlSz];
    for (unsigned int var = 0; var < dsolve::DENDROSOLVER_NUM_VARS; var++) {
        char fName[256];
        sprintf(fName, "%s_%s_n_%d_r_%d_p_%d.bin", fprefix,
                dsolve::DENDROSOLVER_VAR_NAMES[var], nx, rank, npes);
        FILE *outfile = fopen(fName, "w");
        if (outfile == NULL) {
            std::cout << fName << " file open failed " << std::endl;
        }

        for (unsigned int k = kb; k < ke; k++)
            for (unsigned int j = jb; j < je; j++)
                for (unsigned int i = ib; i < ie; i++)
                    blkInternal[k * (ny - 3) * (nx - 3) + j * (nx - 3) + i] =
                        unzipVarsRHS[var]
                                    [offset + k * (ny * nx) + j * (ny) + i];

        fwrite(blkInternal, sizeof(double), blkInlSz,
               outfile);  // write out the number of elements.
        fclose(outfile);
    }

    delete[] blkInternal;
}

unsigned int getOctantWeight(const ot::TreeNode *pNode) {
    return (1u << (3 * pNode->getLevel())) * 1;
}

void computeBHLocations(const ot::Mesh *pMesh, const Point *in, Point *out,
                        double **zipVars, double dt) {
    MPI_Comm commActive = pMesh->getMPICommunicator();

    Point grid_limits[2];
    Point domain_limits[2];

    grid_limits[0] = Point(dsolve::DENDROSOLVER_OCTREE_MIN[0],
                           dsolve::DENDROSOLVER_OCTREE_MIN[1],
                           dsolve::DENDROSOLVER_OCTREE_MIN[2]);
    grid_limits[1] = Point(dsolve::DENDROSOLVER_OCTREE_MAX[0],
                           dsolve::DENDROSOLVER_OCTREE_MAX[1],
                           dsolve::DENDROSOLVER_OCTREE_MAX[2]);

    domain_limits[0] = Point(dsolve::DENDROSOLVER_COMPD_MIN[0],
                             dsolve::DENDROSOLVER_COMPD_MIN[1],
                             dsolve::DENDROSOLVER_COMPD_MIN[2]);
    domain_limits[1] = Point(dsolve::DENDROSOLVER_COMPD_MAX[0],
                             dsolve::DENDROSOLVER_COMPD_MAX[1],
                             dsolve::DENDROSOLVER_COMPD_MAX[2]);

    double beta0[2];
    double beta1[2];
    double beta2[2];

    std::vector<unsigned int> validIndex_beta0;
    std::vector<unsigned int> validIndex_beta1;
    std::vector<unsigned int> validIndex_beta2;

    double beta3vec[6] = {0, 0, 0, 0, 0, 0};
    double bh_pts[6] = {0, 0, 0, 0, 0, 0};

    bh_pts[0] = in[0].x();
    bh_pts[1] = in[0].y();
    bh_pts[2] = in[0].z();

    bh_pts[3] = in[1].x();
    bh_pts[4] = in[1].y();
    bh_pts[5] = in[1].z();

    if (pMesh->isActive()) {
        unsigned int activeRank = pMesh->getMPIRank();

        ot::da::interpolateToCoords(pMesh, zipVars[VAR::U_BETA0], bh_pts, 6,
                                    grid_limits, domain_limits, beta0,
                                    validIndex_beta0);
        ot::da::interpolateToCoords(pMesh, zipVars[VAR::U_BETA1], bh_pts, 6,
                                    grid_limits, domain_limits, beta1,
                                    validIndex_beta1);
        ot::da::interpolateToCoords(pMesh, zipVars[VAR::U_BETA2], bh_pts, 6,
                                    grid_limits, domain_limits, beta2,
                                    validIndex_beta2);

        assert(validIndex_beta0.size() == validIndex_beta1.size());
        assert(validIndex_beta1.size() == validIndex_beta2.size());
    }

    unsigned int red_ranks[2] = {0, 0};
    unsigned int red_ranks_g[2] = {0, 0};
    // global bcast
    for (unsigned int ind = 0; ind < validIndex_beta0.size(); ind++) {
        assert(validIndex_beta0[ind] == validIndex_beta1[ind]);
        assert(validIndex_beta0[ind] == validIndex_beta2[ind]);
        const unsigned int gRank = pMesh->getMPIRankGlobal();

        beta3vec[validIndex_beta0[ind] * 3 + 0] = beta0[validIndex_beta0[ind]];
        beta3vec[validIndex_beta1[ind] * 3 + 1] = beta1[validIndex_beta1[ind]];
        beta3vec[validIndex_beta2[ind] * 3 + 2] = beta2[validIndex_beta2[ind]];

        // std::cout<<"rank: "<<gRank<<"beta["<<(validIndex_beta0[ind]*3)<<"]: (
        // "<<beta3vec[validIndex_beta0[ind]*3 + 0]<<",
        // "<<beta3vec[validIndex_beta0[ind]*3 + 1]<<",
        // "<<beta3vec[validIndex_beta0[ind]*3 + 2]<<")"<<std::endl;
        red_ranks[validIndex_beta2[ind]] = gRank;
    }

    par::Mpi_Allreduce(red_ranks, red_ranks_g, 2, MPI_MAX,
                       pMesh->getMPIGlobalCommunicator());
    MPI_Bcast(&beta3vec[0], 3, MPI_DOUBLE, red_ranks_g[0],
              pMesh->getMPIGlobalCommunicator());
    MPI_Bcast(&beta3vec[3], 3, MPI_DOUBLE, red_ranks_g[1],
              pMesh->getMPIGlobalCommunicator());

    // if(!pMesh->getMPIRankGlobal())
    //     std::cout<<"beta bh0: ( "<<beta3vec[0]<<", "<<beta3vec[1]<<",
    //     "<<beta3vec[2]<<") :  beta 1 ( "<<beta3vec[3]<<", "<<beta3vec[4]<<",
    //     "<<beta3vec[5]<<") "<<std::endl;

    double x[2], y[2], z[2];
    for (unsigned int bh = 0; bh < 2; bh++) {
        x[bh] = in[bh].x() - beta3vec[bh * 3 + 0] * dt;
        y[bh] = in[bh].y() - beta3vec[bh * 3 + 1] * dt;
        z[bh] = in[bh].z() - beta3vec[bh * 3 + 2] * dt;

        out[bh] = Point(x[bh], y[bh], z[bh]);
    }

    return;
}

void allocate_deriv_workspace(const ot::Mesh *pMesh, unsigned int s_fac) {

    // start with deallocation of the workspace, needed when remeshing
    deallocate_deriv_workspace();

    if(!pMesh->isActive())
        return;

    // then get the largest block size from the mesh
    const std::vector<ot::Block> & blkList = pMesh->getLocalBlockList();
    unsigned int max_blk_sz=0;
    for (unsigned int i = 0; i < blkList.size(); i++)
    {
        unsigned int blk_sz = blkList[i].getAllocationSzX() * blkList[i].getAllocationSzY() * blkList[i].getAllocationSzZ();
        if (blk_sz > max_blk_sz)
            max_blk_sz = blk_sz;
    }

    // make sure the derivatives are deallocated? seems unnecessary since it's done earlier?
    deallocate_deriv_workspace();

    // allocate the new memory
    emda::EMDA_DERIV_WORKSPACE = new double[s_fac * max_blk_sz * emda::EMDA_NUM_DERIVATIVES];

}

void deallocate_deriv_workspace() {
    // if the pointer isn't already null, delete the allocated memory 
    if (emda::EMDA_DERIV_WORKSPACE != nullptr) {
        delete[] emda::EMDA_DERIV_WORKSPACE;
        emda::EMDA_DERIV_WORKSPACE = nullptr;
    }
}

}  // end of namespace dsolve

namespace dsolve {

namespace timer {
void initFlops() {
    total_runtime.start();
    t_f2o.start();
    t_cons.start();
    t_bal.start();
    t_mesh.start();
    t_rkSolve.start();
    t_ghostEx_sync.start();
    t_unzip_sync.start();

    for (unsigned int i = 0; i < NUM_FACES; i++)
        dendro::timer::t_unzip_sync_face[i].start();

    dendro::timer::t_unzip_async_internal.start();
    dendro::timer::t_unzip_sync_edge.start();
    dendro::timer::t_unzip_sync_vtex.start();
    dendro::timer::t_unzip_p2c.start();
    dendro::timer::t_unzip_sync_nodalval.start();
    dendro::timer::t_unzip_sync_cpy.start();
    dendro::timer::t_unzip_sync_f_c1.start();
    dendro::timer::t_unzip_sync_f_c2.start();
    dendro::timer::t_unzip_sync_f_c3.start();

    t_unzip_async.start();
    dendro::timer::t_unzip_async_comm.start();

    dendro::timer::t_unzip_async_internal.start();
    dendro::timer::t_unzip_async_external.start();
    dendro::timer::t_unzip_async_comm.start();
    t_deriv.start();
    t_rhs.start();

    t_rhs_a.start();
    t_rhs_b.start();
    t_rhs_gt.start();
    t_rhs_chi.start();
    t_rhs_At.start();
    t_rhs_K.start();
    t_rhs_Gt.start();
    t_rhs_B.start();

    t_bdyc.start();

    t_zip.start();
    t_rkStep.start();
    t_isReMesh.start();
    t_gridTransfer.start();
    t_ioVtu.start();
    t_ioCheckPoint.start();
}

void resetSnapshot() {
    total_runtime.snapreset();
    t_f2o.snapreset();
    t_cons.snapreset();
    t_bal.snapreset();
    t_mesh.snapreset();
    t_rkSolve.snapreset();
    t_ghostEx_sync.snapreset();
    t_unzip_sync.snapreset();

    for (unsigned int i = 0; i < NUM_FACES; i++)
        dendro::timer::t_unzip_sync_face[i].snapreset();

    dendro::timer::t_unzip_sync_internal.snapreset();
    dendro::timer::t_unzip_sync_edge.snapreset();
    dendro::timer::t_unzip_sync_vtex.snapreset();
    dendro::timer::t_unzip_p2c.snapreset();
    dendro::timer::t_unzip_sync_nodalval.snapreset();
    dendro::timer::t_unzip_sync_cpy.snapreset();

    dendro::timer::t_unzip_sync_f_c1.snapreset();
    dendro::timer::t_unzip_sync_f_c2.snapreset();
    dendro::timer::t_unzip_sync_f_c3.snapreset();

    t_unzip_async.snapreset();
    dendro::timer::t_unzip_async_internal.snapreset();
    dendro::timer::t_unzip_async_external.snapreset();
    dendro::timer::t_unzip_async_comm.snapreset();

    t_deriv.snapreset();
    t_rhs.snapreset();

    t_rhs_a.snapreset();
    t_rhs_b.snapreset();
    t_rhs_gt.snapreset();
    t_rhs_chi.snapreset();
    t_rhs_At.snapreset();
    t_rhs_K.snapreset();
    t_rhs_Gt.snapreset();
    t_rhs_B.snapreset();

    t_bdyc.snapreset();

    t_zip.snapreset();
    t_rkStep.snapreset();
    t_isReMesh.snapreset();
    t_gridTransfer.snapreset();
    t_ioVtu.snapreset();
    t_ioCheckPoint.snapreset();
}

void profileInfo(const char *filePrefix, const ot::Mesh *pMesh) {
    int activeRank, activeNpes, globalRank, globalNpes;

    MPI_Comm commActive;
    MPI_Comm commGlobal;

    if (pMesh->isActive()) {
        commActive = pMesh->getMPICommunicator();
        activeRank = pMesh->getMPIRank();
        activeNpes = pMesh->getMPICommSize();
    }

    globalRank = pMesh->getMPIRankGlobal();
    globalNpes = pMesh->getMPICommSizeGlobal();
    commGlobal = pMesh->getMPIGlobalCommunicator();

    double t_stat;
    double t_stat_g[3];

    const char separator = ' ';
    const int nameWidth = 30;
    const int numWidth = 10;

    char fName[256];
    std::ofstream outfile;

    DendroIntL localSz, globalSz;

    if (!activeRank) {
        sprintf(fName, "%s_final.prof", filePrefix);
        outfile.open(fName);
        if (outfile.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            return;
        }

        outfile << "active npes : " << activeNpes << std::endl;
        outfile << "global npes : " << globalNpes << std::endl;
        outfile << "partition tol : " << dsolve::DENDROSOLVER_LOAD_IMB_TOL
                << std::endl;
        outfile << "wavelet tol : " << dsolve::DENDROSOLVER_WAVELET_TOL
                << std::endl;
        outfile << "maxdepth : " << dsolve::DENDROSOLVER_MAXDEPTH << std::endl;
    }

    MPI_Comm comm = commActive;
    unsigned int rank = activeRank;

    localSz = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "Elements : " << globalSz << std::endl;

    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(zip) : " << globalSz << std::endl;

    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(unzip) : " << globalSz << std::endl;

    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "step";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "min(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "mean(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "max(s)" << std::endl;

    t_stat = total_runtime.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "+runtime(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_f2o.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << " ++f2o";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_cons.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << " ++construction";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rkSolve.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << " ++rkSolve";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_bal.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --2:1 balance";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_mesh.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --mesh";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rkStep.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --rkstep";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ghostEx_sync.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --ghostExchge.";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_sync.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_async.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++unzip_async";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

#ifdef ENABLE_DENDRO_PROFILE_COUNTERS
    t_stat = dendro::timer::t_unzip_async_internal.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_internal";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_async_external.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_external";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_async_comm.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_comm (comm) ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;
#endif

    t_stat = t_deriv.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --deriv ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_bdyc.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --boundary con ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_zip.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --zip";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioVtu.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --vtu";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioCheckPoint.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --checkpoint";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    if (!rank) outfile.close();
}

void profileInfoIntermediate(const char *filePrefix, const ot::Mesh *pMesh,
                             const unsigned int currentStep) {
    int activeRank, activeNpes, globalRank, globalNpes;

    MPI_Comm commActive;
    MPI_Comm commGlobal;

    if (pMesh->isActive()) {
        commActive = pMesh->getMPICommunicator();
        activeRank = pMesh->getMPIRank();
        activeNpes = pMesh->getMPICommSize();
    }

    globalRank = pMesh->getMPIRankGlobal();
    globalNpes = pMesh->getMPICommSizeGlobal();
    commGlobal = pMesh->getMPIGlobalCommunicator();

    double t_stat;
    double t_stat_g[3];

    const char separator = ' ';
    const int nameWidth = 30;
    const int numWidth = 10;

    char fName[256];
    std::ofstream outfile;

    DendroIntL localSz, globalSz;

    DendroIntL ghostElements;
    DendroIntL localElements;

    DendroIntL ghostNodes;
    DendroIntL localNodes;

    DendroIntL totalSendNode;
    DendroIntL totalRecvNode;

    DendroIntL numCalls;

#ifdef DENDROSOLVER_PROFILE_HUMAN_READABLE
    if (!activeRank) {
        sprintf(fName, "%s_im.prof", filePrefix);
        outfile.open(fName, std::fstream::app);
        if (outfile.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            return;
        }

        outfile << "active npes : " << activeNpes << std::endl;
        outfile << "global npes : " << globalNpes << std::endl;
        outfile << "current step : " << currentStep << std::endl;
        outfile << "partition tol : " << dsolve::DENDROSOLVER_LOAD_IMB_TOL
                << std::endl;
        outfile << "wavelet tol : " << dsolve::DENDROSOLVER_WAVELET_TOL
                << std::endl;
        outfile << "maxdepth : " << dsolve::DENDROSOLVER_MAXDEPTH << std::endl;
    }

    MPI_Comm comm = commActive;
    unsigned int rank = activeRank;

    localSz = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "Elements : " << globalSz << std::endl;

    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(zip) : " << globalSz << std::endl;

    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(unzip) : " << globalSz << std::endl;

    ghostElements =
        pMesh->getNumPreGhostElements() + pMesh->getNumPostGhostElements();
    localElements = pMesh->getNumLocalMeshElements();

    ghostNodes = pMesh->getNumPreMeshNodes() + pMesh->getNumPostMeshNodes();
    localNodes = pMesh->getNumLocalMeshNodes();

    if (!rank)
        outfile << "========================= MESH "
                   "==========================================================="
                   "============ "
                << std::endl;

    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "step";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "min(#)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "mean(#)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "max(#)" << std::endl;

    t_stat = ghostElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "ghost Elements";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = localElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "local Elements";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = ghostNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "ghost Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = localNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "local Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = pMesh->getGhostExcgTotalSendNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "send Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = pMesh->getGhostExcgTotalRecvNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "recv Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    if (!rank)
        outfile << "========================= RUNTIME "
                   "==========================================================="
                   "======== "
                << std::endl;
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "step";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "min(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "mean(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "max(s)" << std::endl;

    /* t_stat=total_runtime.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"+runtime(s)"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_f2o.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++f2o"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_cons.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++construction"; if(!rank)outfile << std::left
    << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_rkSolve.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++rkSolve"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;*/

    t_stat = t_bal.snap;
    // numCalls=t_bal.num_calls;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++2:1 balance";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_mesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++mesh";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rkStep.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++rkstep";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ghostEx_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++ghostExchge.";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++unzip_sync";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_async.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++unzip_async";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

#ifdef ENABLE_DENDRO_PROFILE_COUNTERS

    t_stat = dendro::timer::t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_comm_wait (comm) ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_nodalval.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_nodalVal";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c1.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c1";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c2.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c2";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c3.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c3";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_cpy.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_cpy";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_internal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_internal";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[0].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_left";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[1].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_right";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[2].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_down";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[3].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_up";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[4].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_back";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[5].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_front";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_edge.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_edge";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_vtex.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_vtex";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_p2c.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_p2c";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;
#endif

    /*
    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
    t_stat=t_unzip_async_internal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_internal"; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

    t_stat=t_unzip_async_external.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_external"; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_comm (comm) "; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
    #endif
    */
    t_stat = t_isReMesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++isReMesh";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_gridTransfer.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++gridTransfer";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_deriv.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++deriv ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++compute_rhs ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_a.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_a ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_b.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_b ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_gt.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_gt ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_chi.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_chi ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_At.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_At ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_K.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_K ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_Gt.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_Gt ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_B.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_B ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_bdyc.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++boundary con ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_zip.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++zip";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioVtu.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++vtu";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioCheckPoint.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++checkpoint";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    if (!rank) outfile.close();
#else

    if (!activeRank) {
        sprintf(fName, "%s_im.prof", filePrefix);
        outfile.open(fName, std::fstream::app);
        if (outfile.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            return;
        }

        // writes the header
        if (currentStep == 0)
            outfile << "step\t act_npes\t glb_npes\t part_tol\t wave_tol\t "
                       "maxdepth\t numOcts\t dof_zip\t dof_unzip\t"
                    << "element_ghost_min\t element_ghost_mean\t "
                       "element_ghost_max\t"
                    << "element_local_min\t element_local_mean\t "
                       "element_local_max\t"
                    << "nodes_local_min\t nodes_local_mean\t nodes_local|max\t"
                    << "send_nodes_min\t send_nodes_mean\t send_nodes_max\t"
                    << "recv_nodes_min\t recv_nodes_mean\t recv_nodes_max\t"
                    << "bal_min\t bal_mean\t bal_max\t"
                    << "mesh_min\t mesh_mean\t mesh_max\t"
                    << "rkstep_min\t rkstep_mean\t rkstep_max\t"
                    << "ghostEx_min\t ghostEx_mean\t ghostEx_max\t"
                    << "unzip_sync_min\t unzip_sync_mean\t unzip_sync_max\t"
                    << "unzip_async_min\t unzip_async_mean\t unzip_async_max\t"
                    << "unzip_async_wait_min\t unzip_async_wait_mean\t "
                       "unzip_async_wait_max\t"
                    << "isRemesh_min\t isRemesh_mean\t isRemesh_max\t"
                    << "GT_min\t GT_mean\t GT_max\t"
                    << "deriv_min\t deriv_mean\t deriv_max\t"
                    << "rhs_min\t rhs_mean\t rhs_max\t" << std::endl;
    }

    MPI_Comm comm = commActive;
    unsigned int rank = activeRank;

    if (!rank) outfile << currentStep << "\t ";
    if (!rank) outfile << activeNpes << "\t ";
    if (!rank) outfile << globalNpes << "\t ";
    if (!rank) outfile << dsolve::DENDROSOLVER_LOAD_IMB_TOL << "\t ";
    if (!rank) outfile << dsolve::DENDROSOLVER_WAVELET_TOL << "\t ";
    if (!rank) outfile << dsolve::DENDROSOLVER_MAXDEPTH << "\t ";

    localSz = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    ghostElements =
        pMesh->getNumPreGhostElements() + pMesh->getNumPostGhostElements();
    localElements = pMesh->getNumLocalMeshElements();

    t_stat = ghostElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = localElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    ghostNodes = pMesh->getNumPreMeshNodes() + pMesh->getNumPostMeshNodes();
    localNodes = pMesh->getNumLocalMeshNodes();

    /*t_stat=ghostNodes;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t
    ";*/

    t_stat = localNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = pMesh->getGhostExcgTotalSendNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = pMesh->getGhostExcgTotalRecvNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_bal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_mesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_rkStep.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_ghostEx_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_unzip_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_unzip_async.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = dendro::timer::t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_isReMesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_gridTransfer.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_deriv.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_rhs.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    if (!rank) outfile << std::endl;
    if (!rank) outfile.close();
#endif
}

}  // namespace timer

}  // namespace dsolve

namespace GW {
void psi4ShpereDump(const ot::Mesh *mesh, DendroScalar **cVar,
                    unsigned int timestep, double time, double dtheta,
                    double dphi) {
    unsigned int rankGlobal = mesh->getMPIRankGlobal();
    unsigned int npesGlobal = mesh->getMPICommSizeGlobal();
    MPI_Comm commGlobal = mesh->getMPIGlobalCommunicator();

    const unsigned int nTheta = (M_PI) / dtheta;
    const unsigned int nPhi = (2 * M_PI) / dphi;
    const unsigned int numPts = nTheta * nPhi;

    unsigned int totalModes = 0;
    for (unsigned int l = 0; l < DENDROSOLVER_GW_NUM_LMODES; l++)
        totalModes += 2 * DENDROSOLVER_GW_L_MODES[l] + 1;

    const unsigned int TOTAL_MODES = totalModes;

    DendroComplex *swsh_coeff =
        new DendroComplex[DENDROSOLVER_GW_NUM_RADAII * TOTAL_MODES];
    DendroComplex *swsh_coeff_g =
        new DendroComplex[DENDROSOLVER_GW_NUM_RADAII * TOTAL_MODES];

    std::vector<unsigned int> lmCounts;
    std::vector<unsigned int> lmOffset;

    lmCounts.resize(DENDROSOLVER_GW_NUM_LMODES);
    lmOffset.resize(DENDROSOLVER_GW_NUM_LMODES);

    for (unsigned int l = 0; l < DENDROSOLVER_GW_NUM_LMODES; l++)
        lmCounts[l] = 2 * DENDROSOLVER_GW_L_MODES[l] + 1;

    lmOffset[0] = 0;
    omp_par::scan(&(*(lmCounts.begin())), &(*(lmOffset.begin())),
                  DENDROSOLVER_GW_NUM_LMODES);

    if (mesh->isActive()) {
        const unsigned int rankActive = mesh->getMPIRank();
        const unsigned int npesActive = mesh->getMPICommSize();

        std::vector<double> coords;
        coords.reserve(3 * numPts);

        std::vector<double> psi4_real;
        psi4_real.resize(numPts);

        std::vector<double> psi4_imag;
        psi4_imag.resize(numPts);

        Point grid_limits[2];
        Point domain_limits[2];

        grid_limits[0] = Point(dsolve::DENDROSOLVER_OCTREE_MIN[0],
                               dsolve::DENDROSOLVER_OCTREE_MIN[1],
                               dsolve::DENDROSOLVER_OCTREE_MIN[2]);
        grid_limits[1] = Point(dsolve::DENDROSOLVER_OCTREE_MAX[0],
                               dsolve::DENDROSOLVER_OCTREE_MAX[1],
                               dsolve::DENDROSOLVER_OCTREE_MAX[2]);

        domain_limits[0] = Point(dsolve::DENDROSOLVER_COMPD_MIN[0],
                                 dsolve::DENDROSOLVER_COMPD_MIN[1],
                                 dsolve::DENDROSOLVER_COMPD_MIN[2]);
        domain_limits[1] = Point(dsolve::DENDROSOLVER_COMPD_MAX[0],
                                 dsolve::DENDROSOLVER_COMPD_MAX[1],
                                 dsolve::DENDROSOLVER_COMPD_MAX[2]);

        std::vector<unsigned int> validIndex;

        for (unsigned int k = 0; k < DENDROSOLVER_GW_NUM_RADAII; k++) {
            for (unsigned int i = 0; i < nTheta; i++)
                for (unsigned int j = 0; j < nPhi; j++) {
                    double x = DENDROSOLVER_GW_RADAII[k] * sin(j * dtheta) *
                               cos(i * dphi);
                    double y = DENDROSOLVER_GW_RADAII[k] * sin(j * dtheta) *
                               sin(i * dphi);
                    double z = DENDROSOLVER_GW_RADAII[k] * cos(j * dtheta);

                    coords.push_back(x);
                    coords.push_back(y);
                    coords.push_back(z);
                }

            validIndex.clear();
            ot::da::interpolateToCoords(
                mesh, cVar[dsolve::VAR_CONSTRAINT::C_PSI4_REAL],
                &(*(coords.begin())), coords.size(), grid_limits, domain_limits,
                &(*(psi4_real.begin())), validIndex);

            validIndex.clear();
            ot::da::interpolateToCoords(
                mesh, cVar[dsolve::VAR_CONSTRAINT::C_PSI4_IMAG],
                &(*(coords.begin())), coords.size(), grid_limits, domain_limits,
                &(*(psi4_imag.begin())), validIndex);
        }
    }
}

}  // namespace GW
