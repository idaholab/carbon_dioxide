///////////////////////////////////////////////////////////////////////////
// LibSBTL_vu_CO2 - SBTL library for carbon dioxide based on:
//
// Span, R. and Wagner, W.:
//
//   "A New Equation of State for Carbon Dioxide Covering the Fluid Region
//    from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa"
//   J. Phys. Chem. Ref. Data, 25(6):1509-1596, 1996.
//
//
// Laesecke, A. and Muzny, C.D.:
//
//   "Reference Correlation for the Viscosity of Carbon Dioxide"
//   J. Phys. Chem. Ref. Data, 46, 013107, 2017.
//
//
// Huber, M.L., Sykioti, E.A., Assael, M.J., and Perkins, R.A.:
//
//   "Reference Correlation of the Thermal Conductivity of Carbon Dioxide
//    from the Triple Point to 1100 K and up to 200 MPa"
//   J. Phys. Chem. Ref. Data, 45, 013102, 2016.
//
//
// Copyright (C) Idaho National Laboratory.
// All rights reserved.
//
// Disclaimer:
// The Idaho National Laboratory (INL) uses its best efforts to deliver a high-quality software and to verify that the computed information is correct.
// However, INL makes no warranties to that effect, and INL shall not be liable for any damage that may result from errors or omissions in the software.
//
// Version: 0.9.0
//
// U_VP_L_CO2
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_CO2.h"
#include "SBTL_call_conv.h"
#include "SBTL_def.h"
//
extern "C" double __stdcall U_VP_L_INI_CO2(double v, double p);
//
extern "C" void __stdcall DIFF_P_VU_L_CO2(double v, double u, double& p, double& dpdv, double& dpdu, double& dudv);
//
SBTLAPI double __stdcall U_VP_L_CO2(double x1_val, double x2_val) throw()
{
    static const double df=1.e-10;       //absolute deviation in f in newtons method
    double u;
    double f=1.e6;
    double p_calc, dpdv_u, dpdu_v, dudv_p;
    int icount=0;
//
    u=U_VP_L_INI_CO2(x1_val,x2_val);
//
//calculate inverse spline by iteration
    while(fabs(f)>df) {
        DIFF_P_VU_L_CO2(x1_val, u, p_calc, dpdv_u, dpdu_v, dudv_p);
        f=p_calc-x2_val;
        u=u-f/dpdu_v;
        if(u<40.)  u=40.;
        if(u>320.) u=320.;
        if(icount++>20) {
            return ERR_VAL;
        }
    }
    return u;
}
//
SBTLAPI void __stdcall DIFF_U_VP_L_CO2(double x1_val, double x2_val, double &u, double& dudv_p, double& dudp_v, double& dpdv_u) throw()
{
    static const double df=1.e-10;       //absolute deviation in f in newtons method
    double f=1.e6;
    double p_calc, dpdu_v;
    int icount=0;
//
    u=U_VP_L_INI_CO2(x1_val,x2_val);
//
//calculate inverse spline by iteration
    while(fabs(f)>df) {
        DIFF_P_VU_L_CO2(x1_val, u, p_calc, dpdv_u, dpdu_v, dudv_p);
        f=p_calc-x2_val;
        u=u-f/dpdu_v;
        if(u<40.)  u=40.;
        if(u>320.) u=320.;
        if(icount++>20) {
            u     =ERR_VAL;
            dudv_p=ERR_VAL;
            dudp_v=ERR_VAL;
            dpdv_u=ERR_VAL;
            return;
        }
    }
    DIFF_P_VU_L_CO2(x1_val, u, p_calc, dpdv_u, dpdu_v, dudv_p);
    dudp_v=1./dpdu_v;
}
