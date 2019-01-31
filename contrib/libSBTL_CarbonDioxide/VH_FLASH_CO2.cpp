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
// HV_FLASH
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_CO2.h"
#include "SBTL_call_conv.h"
//
#define ITMAX 10
//
//initial guess from auxiliary splines
extern "C" double __stdcall U_VH_L_INI_CO2(double v, double h);
extern "C" double __stdcall U_VH_G_INI_CO2_T(double vt, double h);
//
//scaling functions
extern "C" double __stdcall V_U_PMAX_CO2(double u);
extern "C" void   __stdcall DIFF_V_U_PMAX_CO2(double u, double& v, double& dvdu);
extern "C" double __stdcall V1_U_SPL_CO2(double u);
extern "C" void   __stdcall DIFF_V1_U_SPL_CO2(double u, double& v, double& dvdu);
//
//forward functions with derivatives
extern "C" void __stdcall DIFF_P_VU_L_CO2_T(double vs, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_L_CO2_SC(double vs, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_G_CO2_T(double vt, double v, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_G_CO2_TT(double vt, double u, double& p, double& dpdv, double& dpdu, double& dudv);
//
SBTLAPI int __stdcall VH_FLASH_L_CO2(double v, double h, double& u) throw()
{
    double x1tmax,v_u_pmax,dvdu_pmax,dvdu_max;
    double vs,x1tmin,dvdu_min,dvdu_vt,K;
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);
    static const double pc=7.37729837321;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double hc=uc+pc*vc*1.e3;

    static const double df_h=1.e-8;  //-8   //abs. deviation in h

    double hx,px;
    double dpdv_u, dpdu_v, dudv_p;
    double dhdu_v;

    if(fabs(v-vc)/vc<1.e-5 && fabs(h-hc)<0.2) {
        u=uc;
        return I_OK;
    }

    //calculate initial guess
    if(h>hc) {
        u=U_VH_G_INI_CO2_T(log(v),h);
    } else {
        u=U_VH_L_INI_CO2(v,h);
    }

    //newtons method
    double f_h=-1.;
    int icount=0;
    while(fabs(f_h)>df_h) {
        DIFF_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax);
        DIFF_V1_U_SPL_CO2(u, x1tmax, dvdu_max);
        x1tmin=v_u_pmax;
        dvdu_min=dvdu_pmax;
        vs=(v-x1tmin)/(x1tmax-x1tmin)*(x1zmax-x1zmin)+x1zmin;
        dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
        K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
        DIFF_P_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, px, dpdv_u, dpdu_v, dudv_p);
        hx=u+px*v*1.e3;
        dhdu_v=1.+dpdu_v*v*1.e3;
        f_h=hx-h;
        u=u-f_h/dhdu_v;
        if(u<-50.)  u=-50.;
        if(u>1950.) u=1950.;
        if(icount++>ITMAX) {
            u=ERR_VAL;
            return I_ERR;
        }
    }
    return I_OK;
}
//
SBTLAPI void __stdcall VH_FLASH_DERIV_L_CO2(double v, double u, double& dudv_h, double& dudh_v, double& dvdh_u) throw()
{
    double x1tmax,v_u_pmax,dvdu_pmax,dvdu_max,p_;
    double vs,x1tmin,dvdu_min,dvdu_vt,K;
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);
    double dpdv_u, dpdu_v, dudv_p;
    double dhdv_u, dhdu_v;

    DIFF_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax);
    DIFF_V1_U_SPL_CO2(u, x1tmax, dvdu_max);
    x1tmin=v_u_pmax;
    dvdu_min=dvdu_pmax;
    vs=(v-x1tmin)/(x1tmax-x1tmin)*(x1zmax-x1zmin)+x1zmin;
    dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
    K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
    DIFF_P_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, p_, dpdv_u, dpdu_v, dudv_p);
    dhdv_u=(dpdv_u*v+p_)*1.e3;
    dhdu_v=1.+dpdu_v*v*1.e3;
    dudv_h=-dhdv_u/dhdu_v;
    //
    dudh_v=1./dhdu_v;
    dvdh_u=1./dhdv_u;
}
//
SBTLAPI int __stdcall VH_FLASH_G_CO2_T(double v, double vt, double h, double& u) throw()
{
    //static const double pc=7.37729837321;
    //static const double vc=1./467.60000128174;
    //static const double uc=316.468709888;
    //static const double hc=uc+pc*vc*1.e3;

    static const double df_h=1.e-8;  //-8   //abs. deviation in h

    double hx,px;
    double dhdu_v;
    double dpdv_u, dpdu_v, dudv_p;

    //if(fabs(v-vc)/vc<1.e-5 && fabs(h-hc)<0.2) {
    //    return uc;
    //}

    //calculate initial guess
    u=U_VH_G_INI_CO2_T(vt,h);

    //newtons method
    double f_h=-1.;
    int icount=0;
    while(fabs(f_h)>df_h) {
        DIFF_P_VU_G_CO2_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        hx=u+px*v*1.e3;
        dhdu_v=1.+dpdu_v*v*1.e3;
        f_h=hx-h;
        u=u-f_h/dhdu_v;
        if(icount++>ITMAX) {
            u=ERR_VAL;
            return I_ERR;
        }
    }
    return I_OK;
}
//
SBTLAPI void __stdcall VH_FLASH_DERIV_G_CO2(double v, double vt, double u, double& dudv_h, double& dudh_v, double& dvdh_u) throw()
{
    double dpdv_u, dpdu_v, dudv_p;
    double dhdv_u, dhdu_v;
    double p_;

    //derivatives
    DIFF_P_VU_G_CO2_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
    dhdv_u=(dpdv_u*v+p_)*1.e3;
    dhdu_v=1.+dpdu_v*v*1.e3;
    dudv_h=-dhdv_u/dhdu_v;
    //
    dudh_v=1./dhdu_v;
    dvdh_u=1./dhdv_u;
}
//
