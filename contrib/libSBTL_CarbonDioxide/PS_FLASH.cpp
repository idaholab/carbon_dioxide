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
// PH_FLASH
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_CO2.h"
#include "SBTL_call_conv.h"
//
#define ITMAX 20
//
//initial guess from auxiliary splines
extern "C" void  __stdcall VU_SP_L_INI_CO2(double s, double p, double& v, double& u);
extern "C" void  __stdcall VU_SP_G_INI_CO2(double s, double p, double& vt, double& u);
//
//scaling functions
extern "C" double __stdcall V_U_PMAX_CO2(double u);
extern "C" void   __stdcall DIFF_V_U_PMAX_CO2(double u, double& v, double& dvdu);
extern "C" double __stdcall V1_U_SPL_CO2(double u);
extern "C" void   __stdcall DIFF_V1_U_SPL_CO2(double u, double& v, double& dvdu);
//
//forward functions with derivatives
extern "C" void __stdcall DIFF_P_VU_L_CO2_T(double vs, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_L_CO2_T(double vs, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_L_CO2_SC(double vs, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_L_CO2_SC(double vs, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_G_CO2_T(double vt, double v, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_G_CO2_T(double vt, double v, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_G_CO2_TT(double vt, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_G_CO2_TT(double vt, double u, double& p, double& dpdv, double& dpdu, double& dudv);
//
SBTLAPI int __stdcall PS_FLASH_L(double p, double s, double& v, double& u) throw()
{
    double vs,x1tmin,x1tmax,v_u_pmax,dvdu_pmax,dvdu_max;
    double dvdu_min,dvdu_vt,K;
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);
    static const double pc=7.37729837321;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double sc=1.43362534304;

    static const double df_p=1.e-10; //-8   //rel. deviation in p
    static const double df_s=1.e-10; //-8   //abs. deviation in s

    double sx,px,den;
    double dsdv_u, dsdu_v, dudv_s;
    double dpdv_u, dpdu_v, dudv_p;

    if(fabs(p-pc)<0.1 && fabs(s-sc)<0.2) {
        v=vc;
        u=uc;
        return I_OK;
    }

    //calculate initial guesses
    VU_SP_L_INI_CO2(s, p, v, u);

    //newtons method
    double f_p=-1.,f_s=-1.,p_inv=1./p;
    int icount=0;
    while(fabs(f_p*p_inv)>df_p || fabs(f_s)>df_s) {
        DIFF_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax);
        DIFF_V1_U_SPL_CO2(u, x1tmax, dvdu_max);
        x1tmin=v_u_pmax;
        dvdu_min=dvdu_pmax;
        vs=(v-x1tmin)/(x1tmax-x1tmin)*(x1zmax-x1zmin)+x1zmin;
        dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
        K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
        DIFF_S_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, sx, dsdv_u, dsdu_v, dudv_s);      // sx, derivatives
        DIFF_P_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, px, dpdv_u, dpdu_v, dudv_p);      // px, derivatives
        f_p=px-p;
        f_s=sx-s;
        den=dsdu_v*dpdv_u-dsdv_u*dpdu_v;
        v=v+(-dsdu_v*f_p+f_s*dpdu_v)/den;
        u=u+(-f_s*dpdv_u+dsdv_u*f_p)/den;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    return I_OK;
}
//
SBTLAPI void __stdcall PS_FLASH_DERIV_L(double v, double u, double& dvdp_s, double& dvds_p, double& dpds_v, double& dudp_s, double& duds_p, double& dpds_u) throw()
{
    double vs,x1tmin,x1tmax,v_u_pmax,dvdu_pmax,dvdu_max;
    double dvdu_min,dvdu_vt,K,p_,s_;
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);
    double dsdv_u, dsdu_v, dudv_s;
    double dpdv_u, dpdu_v, dudv_p;

    DIFF_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax);
    DIFF_V1_U_SPL_CO2(u, x1tmax, dvdu_max);
    x1tmin=v_u_pmax;
    dvdu_min=dvdu_pmax;
    vs=(v-x1tmin)/(x1tmax-x1tmin)*(x1zmax-x1zmin)+x1zmin;
    dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
    K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
    DIFF_S_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, s_, dsdv_u, dsdu_v, dudv_s);
    DIFF_P_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, p_, dpdv_u, dpdu_v, dudv_p);
    //
    dvdp_s=1./(dpdv_u+dpdu_v*dudv_s);
    dvds_p=1./(dsdv_u+dsdu_v*dudv_p);
    dpds_v=-dvds_p/dvdp_s;
    //
    dudp_s=1./(dpdu_v+dpdv_u/dudv_s);
    duds_p=1./(dsdu_v+dsdv_u/dudv_p);
    dpds_u=-duds_p/dudp_s;
}
//
SBTLAPI int __stdcall PS_FLASH_G(double p, double s, double& v, double& vt, double& u) throw()
{
    static const double pc=7.37729837321;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double sc=1.43362534304;

    static const double df_p=1.e-10; //-8   //rel. deviation in p
    static const double df_s=1.e-10; //-8   //abs. deviation in t

    double sx,px,den;
    double dsdv_u, dsdu_v, dudv_s;
    double dpdv_u, dpdu_v, dudv_p;

    if(fabs(p-pc)<0.1 && fabs(s-sc)<0.2) {
        v=vc;
        vt=log(vc);
        u=uc;
        return I_OK;
    }

    //calculate initial guesses
    VU_SP_G_INI_CO2(s, p, vt, u);

    //newtons method
    double f_p=-1.,f_s=-1.,p_inv=1./p;
    int icount=0;
    while(fabs(f_p*p_inv)>df_p || fabs(f_s)>df_s) {
        DIFF_S_VU_G_CO2_TT(vt, u, sx, dsdv_u, dsdu_v, dudv_s);      // sx, transformed derivatives
        DIFF_P_VU_G_CO2_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        f_p=px-p;
        f_s=sx-s;
        den=dsdu_v*dpdv_u-dsdv_u*dpdu_v;
        vt=vt+(-dsdu_v*f_p+f_s*dpdu_v)/den;
        u =u +(-f_s*dpdv_u+dsdv_u*f_p)/den;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    v=exp(vt);
    return I_OK;
}
//
SBTLAPI void __stdcall PS_FLASH_DERIV_G(double v, double vt, double u, double& dvdp_s, double& dvds_p, double& dpds_v, double& dudp_s, double& duds_p, double& dpds_u) throw()
{
    double dsdv_u, dsdu_v, dudv_s;
    double dpdv_u, dpdu_v, dudv_p;
    double p_,s_;

    //derivatives
    DIFF_P_VU_G_CO2_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
    DIFF_S_VU_G_CO2_T(vt, v, u, s_, dsdv_u, dsdu_v, dudv_s);
    //
    dvdp_s=1./(dpdv_u+dpdu_v*dudv_s);
    dvds_p=1./(dsdv_u+dsdu_v*dudv_p);
    dpds_v=-dvds_p/dvdp_s;
    //
    dudp_s=1./(dpdu_v+dpdv_u/dudv_s);
    duds_p=1./(dsdu_v+dsdv_u/dudv_p);
    dpds_u=-duds_p/dudp_s;
}
//
SBTLAPI int PS_FLASH_G_T(double p, double s, double& vt, double& u) throw()
{
    static const double pc=7.37729837321;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double sc=1.43362534304;

    static const double df_p=1.e-10; //-8   //rel. deviation in p
    static const double df_s=1.e-10; //-8   //abs. deviation in s

    double sx,px,den;
    double dsdv_u, dsdu_v, dudv_s;
    double dpdv_u, dpdu_v, dudv_p;

    if(fabs(p-pc)<0.1 && fabs(s-sc)<0.2) {
        vt=log(vc);
        u=uc;
        return I_OK;
    }

    //calculate initial guesses
    VU_SP_G_INI_CO2(s, p, vt, u);

    //newtons method
    double f_p=-1.,f_s=-1.,p_inv=1./p;
    int icount=0;
    while(fabs(f_p*p_inv)>df_p || fabs(f_s)>df_s) {
        DIFF_S_VU_G_CO2_TT(vt, u, sx, dsdv_u, dsdu_v, dudv_s);      // sx, transformed derivatives
        DIFF_P_VU_G_CO2_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        f_p=px-p;
        f_s=sx-s;
        den=dsdu_v*dpdv_u-dsdv_u*dpdu_v;
        vt=vt+(-dsdu_v*f_p+f_s*dpdu_v)/den;
        u =u +(-f_s*dpdv_u+dsdv_u*f_p)/den;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    return I_OK;
}
