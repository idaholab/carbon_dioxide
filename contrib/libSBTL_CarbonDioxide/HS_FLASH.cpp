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
// HS_FLASH
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
extern "C" void  __stdcall VU_SH_L_INI_CO2(double s, double h, double& v, double& u);
extern "C" void  __stdcall VU_SH_G_INI_CO2(double s, double h, double& vt, double& u);
//
// auxiliary functions
extern "C" double __stdcall H_S_UC_AUX_CO2(double s);
//
//scaling functions
extern "C" double __stdcall V_U_PMAX_CO2(double u);
extern "C" void   __stdcall DIFF_V_U_PMAX_CO2(double u, double& v, double& dvdu);
extern "C" double __stdcall V1_U_SPL_CO2(double u);
extern "C" void   __stdcall DIFF_V1_U_SPL_CO2(double u, double& v, double& dvdu);
//
//forward functions with derivatives
extern "C" void __stdcall DIFF_P_VU_L_CO2_T(double vs, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_T_VU_L_CO2_T(double vs, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& t, double& dtdv, double& dtdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_L_CO2_T(double vs, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_L_CO2_SC(double vs, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_L_CO2_SC(double vs, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_G_CO2_T(double vt, double v, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_T_VU_G_CO2_T(double vt, double v, double u, double& t, double& dtdv, double& dtdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_G_CO2_T(double vt, double v, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_G_CO2_TT(double vt, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_G_CO2_TT(double vt, double u, double& p, double& dpdv, double& dpdu, double& dudv);
//
SBTLAPI int __stdcall HS_FLASH_L(double h, double s, double& v, double& u) throw()
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
    static const double sc=1.43362534304;

    static const double df_h=1.e-8;  //-8   //abs. deviation in h
    static const double df_s=1.e-10; //-8   //abs. deviation in s

    double hx,sx,px,den;
    double dhdv_u, dhdu_v;
    double dpdv_u, dpdu_v, dudv_p;
    double dsdv_u, dsdu_v, dudv_s;

    if(fabs(h-hc)<0.2 && fabs(s-sc)<0.2) {
        v=vc;
        u=uc;
        return I_OK;
    }

    //calculate initial guesses
    VU_SH_L_INI_CO2(s, h, v, u);

    //newtons method
    double f_h=-1.,f_s=-1.;
    int icount=0;
    while(fabs(f_h)>df_h || fabs(f_s)>df_s) {
        DIFF_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax);  
        DIFF_V1_U_SPL_CO2(u, x1tmax, dvdu_max);
        x1tmin=v_u_pmax;
        dvdu_min=dvdu_pmax;
        vs=(v-x1tmin)/(x1tmax-x1tmin)*(x1zmax-x1zmin)+x1zmin;
        dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
        K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
        DIFF_P_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, px, dpdv_u, dpdu_v, dudv_p);      // px, derivatives
        DIFF_S_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, sx, dsdv_u, dsdu_v, dudv_s);      // sx, derivatives
        hx=u+px*v*1.e3;
        dhdv_u=(dpdv_u*v+px)*1.e3;
        dhdu_v=1.+dpdu_v*v*1.e3;
        f_h=hx-h;
        f_s=sx-s;
        den=dhdv_u*dsdu_v-dhdu_v*dsdv_u;
        v=v+(-dsdu_v*f_h+f_s*dhdu_v)/den;
        u=u+(-f_s*dhdv_u+dsdv_u*f_h)/den;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    return I_OK;
}
//
SBTLAPI void __stdcall HS_FLASH_DERIV_L(double v, double u, double& dvdh_s, double& dvds_h, double& dhds_v, double& dudh_s, double& duds_h, double& dhds_u) throw()
{
    double x1tmax,v_u_pmax,dvdu_pmax,dvdu_max,p_,s_;
    double vs,x1tmin,dvdu_min,dvdu_vt,K;
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);
    double dhdv_u, dhdu_v, dudv_h;
    double dpdv_u, dpdu_v, dudv_p;
    double dsdv_u, dsdu_v, dudv_s;

    DIFF_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax);  
    DIFF_V1_U_SPL_CO2(u, x1tmax, dvdu_max);
    x1tmin=v_u_pmax;
    dvdu_min=dvdu_pmax;
    vs=(v-x1tmin)/(x1tmax-x1tmin)*(x1zmax-x1zmin)+x1zmin;
    dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
    K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
    DIFF_P_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, p_, dpdv_u, dpdu_v, dudv_p);
    DIFF_S_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, s_, dsdv_u, dsdu_v, dudv_s);
    dhdv_u=(dpdv_u*v+p_)*1.e3;
    dhdu_v=1.+dpdu_v*v*1.e3;
    dudv_h=-dhdv_u/dhdu_v;
    //
    dvdh_s=1./(dhdv_u+dhdu_v*dudv_s);
    dvds_h=1./(dsdv_u+dsdu_v*dudv_h);
    dhds_v=-dvds_h/dvdh_s;
    //
    dudh_s=1./(dhdu_v+dhdv_u/dudv_s);
    duds_h=1./(dsdu_v+dsdv_u/dudv_h);
    dhds_u=-duds_h/dudh_s;
}
//
SBTLAPI void __stdcall HS_PT_FLASH_DERIV_L(double v, double u, double& p, double& dpdh_s, double& dpds_h, double& dhds_p, double& t, double& dtdh_s, double& dtds_h, double& dhds_t) throw()
{
    double x1tmax,v_u_pmax,dvdu_pmax,dvdu_max,s_;
    double vs,x1tmin,dvdu_min,dvdu_vt,K;
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);
    double dhdv_u, dhdu_v;
    double dpdv_u, dpdu_v, dudv_p;
    double dtdv_u, dtdu_v, dudv_t;
    double dsdv_u, dsdu_v, dudv_s;

    DIFF_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax);  
    DIFF_V1_U_SPL_CO2(u, x1tmax, dvdu_max);
    x1tmin=v_u_pmax;
    dvdu_min=dvdu_pmax;
    vs=(v-x1tmin)/(x1tmax-x1tmin)*(x1zmax-x1zmin)+x1zmin;
    dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
    K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
    DIFF_P_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, p, dpdv_u, dpdu_v, dudv_p);
    DIFF_T_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, t, dtdv_u, dtdu_v, dudv_t);
    DIFF_S_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, s_, dsdv_u, dsdu_v, dudv_s);
    dhdv_u=(dpdv_u*v+p)*1.e3;
    dhdu_v=1.+dpdu_v*v*1.e3;
    //
    dpdh_s=(dpdv_u*dsdu_v-dpdu_v*dsdv_u)/(dhdv_u*dsdu_v-dhdu_v*dsdv_u);
    dpds_h=(dpdv_u*dhdu_v-dpdu_v*dhdv_u)/(dsdv_u*dhdu_v-dsdu_v*dhdv_u);
    dhds_p=-dpds_h/dpdh_s;
    //
    dtdh_s=(dtdv_u*dsdu_v-dtdu_v*dsdv_u)/(dhdv_u*dsdu_v-dhdu_v*dsdv_u);
    dtds_h=(dtdv_u*dhdu_v-dtdu_v*dhdv_u)/(dsdv_u*dhdu_v-dsdu_v*dhdv_u);
    dhds_t=-dtds_h/dtdh_s;
}
//
SBTLAPI int __stdcall HS_FLASH_G(double h, double s, double& v, double& vt, double& u) throw()
{
    static const double pc=7.37729837321;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double hc=uc+pc*vc*1.e3;
    static const double sc=1.43362534304;

    static const double df_h=1.e-8;  //-8   //abs. deviation in h
    static const double df_s=1.e-10; //-8   //abs. deviation in s

    double hg;
    double hx,sx,px,den;
    double dhdv_u, dhdu_v;
    double dpdv_u, dpdu_v, dudv_p;
    double dsdv_u, dsdu_v, dudv_s;

    if(fabs(h-hc)<0.2 && fabs(s-sc)<0.2) {
        v=vc;
        u=uc;
        return I_OK;
    }

    //calculate initial guesses
    if(s<sc) {
        hg=H_S_UC_AUX_CO2(s);
        if(h>hg) {
            VU_SH_G_INI_CO2(s, h, vt, u);
            v=exp(vt);
        } else {
            VU_SH_L_INI_CO2(s, h, v, u);
            vt=log(v);
        }
    } else {
        VU_SH_G_INI_CO2(s, h, vt, u);
        v=exp(vt);
    }

    //newtons method
    double f_h=-1.,f_s=-1.;
    int icount=0;
    while(fabs(f_h)>df_h || fabs(f_s)>df_s) {
        DIFF_P_VU_G_CO2_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        DIFF_S_VU_G_CO2_TT(vt, u, sx, dsdv_u, dsdu_v, dudv_s);      // sx, transformed derivatives
        hx=u+px*v*1.e3;
        dhdv_u=(dpdv_u*v+px*v)*1.e3;
        dhdu_v=1.+dpdu_v*v*1.e3;
        f_h=hx-h;
        f_s=sx-s;
        den=dhdv_u*dsdu_v-dhdu_v*dsdv_u;
        vt=vt+(-dsdu_v*f_h+f_s*dhdu_v)/den;
        u =u +(-f_s*dhdv_u+dsdv_u*f_h)/den;
        v=exp(vt);
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    return I_OK;
}
//
SBTLAPI void __stdcall HS_FLASH_DERIV_G(double v, double vt, double u, double& dvdh_s, double& dvds_h, double& dhds_v, double& dudh_s, double& duds_h, double& dhds_u) throw()
{
    double dhdv_u, dhdu_v, dudv_h;
    double dpdv_u, dpdu_v, dudv_p;
    double dsdv_u, dsdu_v, dudv_s;
    double p_,s_;

    //derivatives
    DIFF_P_VU_G_CO2_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
    DIFF_S_VU_G_CO2_T(vt, v, u, s_, dsdv_u, dsdu_v, dudv_s);
    dhdv_u=(dpdv_u*v+p_)*1.e3;
    dhdu_v=1.+dpdu_v*v*1.e3;
    dudv_h=-dhdv_u/dhdu_v;
    //
    dvdh_s=1./(dhdv_u+dhdu_v*dudv_s);
    dvds_h=1./(dsdv_u+dsdu_v*dudv_h);
    dhds_v=-dvds_h/dvdh_s;
    //
    dudh_s=1./(dhdu_v+dhdv_u/dudv_s);
    duds_h=1./(dsdu_v+dsdv_u/dudv_h);
    dhds_u=-duds_h/dudh_s;
}
//
SBTLAPI void __stdcall HS_PT_FLASH_DERIV_G(double v, double vt, double u, double& p, double& dpdh_s, double& dpds_h, double& dhds_p, double& t, double& dtdh_s, double& dtds_h, double& dhds_t) throw()
{
    double dhdv_u, dhdu_v;
    double dpdv_u, dpdu_v, dudv_p;
    double dtdv_u, dtdu_v, dudv_t;
    double dsdv_u, dsdu_v, dudv_s;
    double s_;

    //derivatives
    DIFF_P_VU_G_CO2_T(vt, v, u, p, dpdv_u, dpdu_v, dudv_p);
    DIFF_T_VU_G_CO2_T(vt, v, u, t, dtdv_u, dtdu_v, dudv_t);
    DIFF_S_VU_G_CO2_T(vt, v, u, s_, dsdv_u, dsdu_v, dudv_s);
    dhdv_u=(dpdv_u*v+p)*1.e3;
    dhdu_v=1.+dpdu_v*v*1.e3;
    //
    dpdh_s=(dpdv_u*dsdu_v-dpdu_v*dsdv_u)/(dhdv_u*dsdu_v-dhdu_v*dsdv_u);
    dpds_h=(dpdv_u*dhdu_v-dpdu_v*dhdv_u)/(dsdv_u*dhdu_v-dsdu_v*dhdv_u);
    dhds_p=-dpds_h/dpdh_s;
    //
    dtdh_s=(dtdv_u*dsdu_v-dtdu_v*dsdv_u)/(dhdv_u*dsdu_v-dhdu_v*dsdv_u);
    dtds_h=(dtdv_u*dhdu_v-dtdu_v*dhdv_u)/(dsdv_u*dhdu_v-dsdu_v*dhdv_u);
    dhds_t=-dtds_h/dtdh_s;
}