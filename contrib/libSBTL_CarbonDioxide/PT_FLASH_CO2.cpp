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
// PT_FLASH
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
extern "C" void  __stdcall VU_TP_L_INI_CO2(double t, double p, double& v, double& u);
extern "C" void  __stdcall VU_TP_G_INI_CO2(double t, double p, double& vt, double& u);
//extern "C" double __stdcall V2_P_AUX_T(double pt);
//extern "C" double __stdcall U2_T_AUX(double t);
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
extern "C" void __stdcall DIFF_T_VU_L_CO2_SC(double vs, double u, double& t, double& dtdv, double& dtdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_L_CO2_SC(double vs, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_G_CO2_T(double vt, double v, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_T_VU_G_CO2_T(double vt, double v, double u, double& t, double& dtdv, double& dtdu, double& dudv);
extern "C" void __stdcall DIFF_T_VU_G_CO2_TT(double vt, double u, double& t, double& dtdv, double& dtdu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_G_CO2_TT(double vt, double u, double& p, double& dpdv, double& dpdu, double& dudv);
//
SBTLAPI double __stdcall V1_T_AUX_CO2(double t);
SBTLAPI double __stdcall U1_T_AUX_CO2(double t);
SBTLAPI double __stdcall PS_T_INV_CO2(double t);
SBTLAPI double __stdcall V2_P_AUX_CO2_T(double pt);
SBTLAPI double __stdcall U2_T_AUX_CO2(double t);
//
SBTLAPI int __stdcall PT_FLASH_L_CO2(double p, double t, double& v, double& u) throw()
{
    double vs,x1tmin,x1tmax;
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double psxx=7.;

    static const double df_p=1.e-10; //-8   //rel. deviation in p
    static const double df_t=1.e-10; //-8   //abs. deviation in t

    double tx,px,den;
    double dtdv_u, dtdu_v, dudv_t;
    double dpdv_u, dpdu_v, dudv_p;

    if(fabs(p-pc)<0.1 && fabs(t-tc)<0.2) {
        v=vc;
        u=uc;
        return I_OK;
    }

    //calculate initial guesses
    VU_TP_L_INI_CO2(t, p, v, u);
    if(p>psxx && p<pc && t<tc) {
        double v1,u1;
        v1=V1_T_AUX_CO2(t);
        u1=U1_T_AUX_CO2(t);
        if(u>u1) u=u1;
        if(v>v1) v=v1;
    }

    //scale v (iteration in transformed coordinates)
    x1tmin=V_U_PMAX_CO2(u);
    x1tmax=V1_U_SPL_CO2(u);
    vs=(v-x1tmin)/(x1tmax-x1tmin)*(x1zmax-x1zmin)+x1zmin;

    //newtons method
    double f_p=-1.,f_t=-1.,p_inv=1./p;
    int icount=0;
    while(fabs(f_p*p_inv)>df_p || fabs(f_t)>df_t) {
        DIFF_T_VU_L_CO2_SC(vs, u, tx, dtdv_u, dtdu_v, dudv_t);      // tx, scaled derivatives
        DIFF_P_VU_L_CO2_SC(vs, u, px, dpdv_u, dpdu_v, dudv_p);      // px, scaled derivatives
        f_p=px-p;
        f_t=tx-t;
        den=dtdu_v*dpdv_u-dtdv_u*dpdu_v;
        vs=vs+(-dtdu_v*f_p+f_t*dpdu_v)/den;
        u =u +(-f_t*dpdv_u+dtdv_u*f_p)/den;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    x1tmin=V_U_PMAX_CO2(u);
    x1tmax=V1_U_SPL_CO2(u);
    v=(vs-x1zmin)/(x1zmax-x1zmin)*(x1tmax-x1tmin)+x1tmin;
    return I_OK;
}
//
SBTLAPI int __stdcall PT_FLASH_DERIV_L_CO2(double p, double t, double& v, double& dvdp_t, double& dvdt_p, double& dpdt_v, double& u, double& dudp_t, double& dudt_p, double& dpdt_u) throw()
{
    double vs,x1tmin,x1tmax,p_,t_;
    double v_u_pmax,dvdu_pmax,v_u_spndl,dvdu_spndl,dvdu_min,dvdu_max,K,dvdu_vt;
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double psxx=7.;
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);

    static const double df_p=1.e-10; //-8  //rel. deviation in p
    static const double df_t=1.e-10; //-8  //abs. deviation in t

    double tx,px,den;
    double dtdv_u, dtdu_v, dudv_t;
    double dpdv_u, dpdu_v, dudv_p;

    if(fabs(p-pc)<0.1 && fabs(t-tc)<0.2) {
        v=vc;
        u=uc;
        return I_OK;
    }

    //calculate initial guesses
    VU_TP_L_INI_CO2(t, p, v, u);
    if(p>psxx && p<pc && t<tc) {
        double v1,u1;
        v1=V1_T_AUX_CO2(t);
        u1=U1_T_AUX_CO2(t);
        if(u>u1) u=u1;
        if(v>v1) v=v1;
    }

    //scale v (iteration in transformed coordinates)
    x1tmin=V_U_PMAX_CO2(u);
    x1tmax=V1_U_SPL_CO2(u);
    vs=(v-x1tmin)/(x1tmax-x1tmin)*(x1zmax-x1zmin)+x1zmin;

    //newtons method
    double f_p=-1.,f_t=-1.,p_inv=1./p;
    int icount=0;
    while(fabs(f_p*p_inv)>df_p || fabs(f_t)>df_t) {
        DIFF_T_VU_L_CO2_SC(vs, u, tx, dtdv_u, dtdu_v, dudv_t);      // tx, scaled derivatives
        DIFF_P_VU_L_CO2_SC(vs, u, px, dpdv_u, dpdu_v, dudv_p);      // px, scaled derivatives
        f_p=px-p;
        f_t=tx-t;
        den=dtdu_v*dpdv_u-dtdv_u*dpdu_v;
        vs=vs+(-dtdu_v*f_p+f_t*dpdu_v)/den;
        u =u +(-f_t*dpdv_u+dtdv_u*f_p)/den;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    DIFF_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax);
    DIFF_V1_U_SPL_CO2(u, v_u_spndl, dvdu_spndl);
    x1tmin=v_u_pmax;
    x1tmax=v_u_spndl;
    dvdu_min=dvdu_pmax;
    dvdu_max=dvdu_spndl;
    K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
    v=(vs-x1zmin)/K+x1tmin;
    //derivatives
    dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
    DIFF_P_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, p_, dpdv_u, dpdu_v, dudv_p);
    DIFF_T_VU_L_CO2_T(vs, x1tmin, x1tmax, dvdu_vt, K, u, t_, dtdv_u, dtdu_v, dudv_t);
    //
    dvdp_t=1./(dpdv_u+dpdu_v*dudv_t);
    dvdt_p=1./(dtdv_u+dtdu_v*dudv_p);
    dpdt_v=-dvdt_p/dvdp_t;
    //
    dudp_t=1./(dpdu_v+dpdv_u/dudv_t);
    dudt_p=1./(dtdu_v+dtdv_u/dudv_p);
    dpdt_u=-dudt_p/dudp_t;
    return I_OK;
}
//
SBTLAPI int __stdcall PT_FLASH_G_CO2(double p, double t, double& v, double& vt, double& u) throw()
{
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double psxx=7.;

    static const double df_p=1.e-10; //-8   //rel. deviation in p
    static const double df_t=1.e-10; //-8   //abs. deviation in t

    double tx,px,den;
    double dtdv_u, dtdu_v, dudv_t;
    double dpdv_u, dpdu_v, dudv_p;

    if(fabs(p-pc)<0.1 && fabs(t-tc)<0.2) {
        v=vc;
        vt=log(vc);
        u=uc;
        return I_OK;
    }

    //calculate initial guesses
    VU_TP_G_INI_CO2(t, p, vt, u);
    if(p>psxx && p<pc && t<tc) {
        double ps,pst,v2t,u2;
        ps=PS_T_INV_CO2(t);
        pst=log(ps);
        v2t=V2_P_AUX_CO2_T(pst);
        u2=U2_T_AUX_CO2(t);
        if(u<u2) u=u2;
        if(vt<v2t) vt=v2t;
    }

    //newtons method
    double f_p=-1.,f_t=-1.,p_inv=1./p;
    int icount=0;
    while(fabs(f_p*p_inv)>df_p || fabs(f_t)>df_t) {
        DIFF_T_VU_G_CO2_TT(vt, u, tx, dtdv_u, dtdu_v, dudv_t);      // tx, transformed derivatives
        DIFF_P_VU_G_CO2_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        f_p=px-p;
        f_t=tx-t;
        den=dtdu_v*dpdv_u-dtdv_u*dpdu_v;
        vt=vt+(-dtdu_v*f_p+f_t*dpdu_v)/den;
        u =u +(-f_t*dpdv_u+dtdv_u*f_p)/den;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    v=exp(vt);
    return I_OK;
}

SBTLAPI int __stdcall PT_FLASH_DERIV_G_CO2(double p, double t, double& v, double& vt, double& dvdp_t, double& dvdt_p, double& dpdt_v, double& u, double& dudp_t, double& dudt_p, double& dpdt_u) throw()
{
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double psxx=7.;

    static const double df_p=1.e-10; //-8   //rel. deviation in p
    static const double df_t=1.e-10; //-8   //abs. deviation in t

    double tx,px,den;
    double dtdv_u, dtdu_v, dudv_t;
    double dpdv_u, dpdu_v, dudv_p;
    double p_,t_;

    if(fabs(p-pc)<0.1 && fabs(t-tc)<0.2) {
        v=vc;
        vt=log(vc);
        u=uc;
        return I_OK;
    }

    //calculate initial guesses
    VU_TP_G_INI_CO2(t, p, vt, u);
    if(p>psxx && p<pc && t<tc) {
        double ps,pst,v2t,u2;
        ps=PS_T_INV_CO2(t);
        pst=log(ps);
        v2t=V2_P_AUX_CO2_T(pst);
        u2=U2_T_AUX_CO2(t);
        if(u<u2) u=u2;
        if(vt<v2t) vt=v2t;
    }

    //newtons method
    double f_p=-1.,f_t=-1.,p_inv=1./p;
    int icount=0;
    while(fabs(f_p*p_inv)>df_p || fabs(f_t)>df_t) {
        DIFF_T_VU_G_CO2_TT(vt, u, tx, dtdv_u, dtdu_v, dudv_t);      // tx, transformed derivatives
        DIFF_P_VU_G_CO2_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        f_p=px-p;
        f_t=tx-t;
        den=dtdu_v*dpdv_u-dtdv_u*dpdu_v;
        vt=vt+(-dtdu_v*f_p+f_t*dpdu_v)/den;
        u =u +(-f_t*dpdv_u+dtdv_u*f_p)/den;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    v=exp(vt);
    //derivatives
    DIFF_P_VU_G_CO2_T(vt, v, u, p_, dpdv_u, dpdu_v, dudv_p);
    DIFF_T_VU_G_CO2_T(vt, v, u, t_, dtdv_u, dtdu_v, dudv_t);
    //
    dvdp_t=1./(dpdv_u+dpdu_v*dudv_t);
    dvdt_p=1./(dtdv_u+dtdu_v*dudv_p);
    dpdt_v=-dvdt_p/dvdp_t;
    //
    dudp_t=1./(dpdu_v+dpdv_u/dudv_t);
    dudt_p=1./(dtdu_v+dtdv_u/dudv_p);
    dpdt_u=-dudt_p/dudp_t;
    return I_OK;
}
//
SBTLAPI int PT_FLASH_G_CO2_T(double p, double t, double& vt, double& u) throw()
{
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double psxx=7.;

    static const double df_p=1.e-10; //-8   //rel. deviation in p
    static const double df_t=1.e-10; //-8   //abs. deviation in t

    double /*v,*/tx,px,den;
    double dtdv_u, dtdu_v, dudv_t;
    double dpdv_u, dpdu_v, dudv_p;

    if(fabs(p-pc)<0.1 && fabs(t-tc)<0.2) {
        vt=log(vc);
        u=uc;
        return I_OK;
    }

    //calculate initial guesses
    VU_TP_G_INI_CO2(t, p, vt, u);
    if(p>psxx && p<pc && t<tc) {
        double ps,pst,v2t,u2;
        ps=PS_T_INV_CO2(t);
        pst=log(ps);
        v2t=V2_P_AUX_CO2_T(pst);
        u2=U2_T_AUX_CO2(t);
        if(u<u2) u=u2;
        if(vt<v2t) vt=v2t;
    }

    //newtons method
    double f_p=-1.,f_t=-1.,p_inv=1./p;
    int icount=0;
    while(fabs(f_p*p_inv)>df_p || fabs(f_t)>df_t) {
        DIFF_T_VU_G_CO2_TT(vt, u, tx, dtdv_u, dtdu_v, dudv_t);      // tx, transformed derivatives
        DIFF_P_VU_G_CO2_TT(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        f_p=px-p;
        f_t=tx-t;
        den=dtdu_v*dpdv_u-dtdv_u*dpdu_v;
        vt=vt+(-dtdu_v*f_p+f_t*dpdu_v)/den;
        u =u +(-f_t*dpdv_u+dtdv_u*f_p)/den;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    return I_OK;
}
