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
// SAT_VU_CO2
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_CO2.h"
#include "SBTL_call_conv.h"

#ifdef __INTEL_COMPILER
#include "mkl_lapack.h"
#else
//#include "lapack.h"
extern "C" void dgetrf_( int* m, int* n, double* a, int* lda, int* ipiv, int* info );
extern "C" void dgetrs_( char* trans, int* n, int* nrhs, 
             double* a, int* lda, int* ipiv, 
             double* b, int* ldb, int* info );
#endif

#define ITMAX 20

//initial guesses from auxiliary splines
extern "C" double __stdcall PS_VU_INI_CO2(double v, double u);
extern "C" double __stdcall V1_T_AUX_CO2(double t);
extern "C" double __stdcall V2_P_AUX_CO2_T(double pt);
extern "C" double __stdcall U1_T_AUX_CO2(double t);
extern "C" double __stdcall U2_T_AUX_CO2(double t);
extern "C" double __stdcall V1_H_AUX_CO2(double h);
extern "C" double __stdcall U2_V_AUX_CO2_T(double vt);
extern "C" double __stdcall PS_U_AUX_CO2_T(double u);
extern "C" double __stdcall PS_V_AUX_CO2_T(double vt);
extern "C" double __stdcall PS_H_AUX_CO2_T(double h);
extern "C" double __stdcall PS_S_AUX_CO2_T(double s);
extern "C" double __stdcall PS_SH_INI_CO2_T(double h, double s);
extern "C" double __stdcall PS_VH_INI_CO2_T(double vt, double h);

//
//perform iteration in transformed coordinates (ln(v), vt(u) , and sqrt(p))
extern "C" double __stdcall TS_P_CO2_T(double pt);
extern "C" void   __stdcall DIFF_TS_P_CO2_T(double pt, double& ts, double& dtsdpt);
extern "C" void   __stdcall DIFF_P_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void   __stdcall DIFF_T_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& t, double& dtdv, double& dtdu, double& dudv);
extern "C" double __stdcall S_VU_L_CO2_T(double vt, double u);
extern "C" void   __stdcall DIFF_S_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void   __stdcall DIFF_P_VU_G_CO2_T(double vt, double v, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void   __stdcall DIFF_T_VU_G_CO2_T(double vt, double v, double u, double& t, double& dtdv, double& dtdu, double& dudv);
extern "C" void   __stdcall DIFF_S_VU_G_CO2_T(double vt, double v, double u, double& s, double& dsdv, double& dsdu, double& dudv);
//
extern "C" void __stdcall DIFF_P_VU_L_CO2_SC(double vt, double u, double& p, double& dpdv_u, double& dpdu_v, double& dudv_p);
extern "C" void __stdcall DIFF_T_VU_L_CO2_SC(double vt, double u, double& t, double& dtdv_u, double& dtdu_v, double& dudv_t);
extern "C" void __stdcall DIFF_P_VU_G_CO2_TT(double vt, double uv, double& p, double& dpdv_u, double& dpdu_v, double& dudv_p);
extern "C" void __stdcall DIFF_T_VU_G_CO2_TT(double vt, double uv, double& t, double& dtdv_u, double& dtdu_v, double& dudv_t);
extern "C" void __stdcall DIFF_S_VU_G_CO2_TT(double vt, double uv, double& s, double& dsdv_u, double& dsdu_v, double& dudv_s);
//
extern "C" double __stdcall V_U_PMAX_CO2(double u);
extern "C" void   __stdcall DIFF_V_U_PMAX_CO2(double u, double& v, double& dvdu);
extern "C" double __stdcall V1_U_SPL_CO2(double u);
extern "C" void   __stdcall DIFF_V1_U_SPL_CO2(double u, double& v, double& dvdu);
//
extern "C" void   __stdcall DIFF_TS_P_CO2(double p, double& ts, double& dtsdp);
//
SBTLAPI int __stdcall SAT_VU_SPL(double v, double u, double& ps, double& ts, double& xvap, double& vl, double& vv, double& vvt, double& ul, double& uv) throw()
{
    //declare variables and constants
    static const double tol_ps=1.e-10; //rel. tolerance in ps
    static const double tol_ts=1.e-8;  //abs. tolerance in ts
    static const double tol_x =1.e-8;  //abs. tolerance in x
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    //static const double hc=uc+pc*vc*1.e3;
    //static const double sc=1.43362534304;
    static const double pcx=7.33;
    //static const double tcx=303.84873305000252;
    static const double v1x=0.0018221754024023962;
    static const double v2x=0.0021592622018924488;
    static const double u1x=301.03552539317360;
    static const double u2x=316.77324527264307;
    //static const double h1x=u1x+pcx*v1x*1.e3;
    //static const double h2x=u2x+pcx*v2x*1.e3;
    //static const double s1x=1.3752027569977439;
    //static const double s2x=1.4351295210653687;
    static const double mx=(u2x-u1x)/(v2x-v1x);
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);
    double x1tmin, x1tmax, vlt;
    double dvdu_min, dvdu_max, K, dvdu_vt;
    double pl, dpdv_u_l, dpdu_v_l, dudv_p_l;
    double pv, dpdv_u_v, dpdu_v_v, dudv_p_v;
    double tl, dtdv_u_l, dtdu_v_l, dudv_t_l;
    double tv, dtdv_u_v, dtdu_v_v, dudv_t_v;
    double dtsdpt;
    double du, du_inv, du_inv2;
    double dv, dv_inv, dv_inv2;
    double dpl,dpv,dtl,dtv,dx;

    //declare arrays
    double J[5][5], F[5];

    //declare and initialize flags and return values
    char TRANS='N';         //dgetrs does not transpose the matrix for LU-decomposition
    int INFO=0;             //return value (0-OK)
    int  LDA=5;             //leading dimension of J
    int  LDB=5;             //leading dimension of F (columns)
    int    N=5;             //number of rows in F
    int NRHS=1;             //number of right hand sides
    int IPIV[5];            //indices of pivot elements

    xvap=ERR_VAL;
 
    //vicinity of the critical point
    if(u>(mx*(v-v1x)+u1x)) {
        ps=pc;
        ts=tc;
        xvap=0.5;
        vl=vc;
        vv=vc;
        vvt=log(vc);
        ul=uc;
        uv=uc;
        return I_OK;
    }

    //initialize values, use the transformed variable pst during the iteration
    double pst;
    ps =PS_VU_INI_CO2(v,u);
    if(ps>pcx) ps=pcx;
    pst=sqrt(ps);
    ts =TS_P_CO2_T(pst);
    vl=V1_T_AUX_CO2(ts);
    double lnp=log(ps);
    vvt=V2_P_AUX_CO2_T(lnp);
    vv=exp(vvt);
    ul=U1_T_AUX_CO2(ts);
    uv=U2_T_AUX_CO2(ts);

    int icount=0;
    do {
        //populate right hand side and jacobian matrix - use dpt derivatives, rather than dp 
        DIFF_V_U_PMAX_CO2(ul, x1tmin, dvdu_min);
        DIFF_V1_U_SPL_CO2(ul, x1tmax, dvdu_max);
        K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
        vlt=(vl-x1tmin)*K+x1zmin;
        dvdu_vt=(vlt-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
        DIFF_P_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, pl, dpdv_u_l, dpdu_v_l, dudv_p_l);
        DIFF_P_VU_G_CO2_T(vvt, vv, uv, pv, dpdv_u_v, dpdu_v_v, dudv_p_v);
        DIFF_TS_P_CO2_T(pst, ts, dtsdpt);
        DIFF_T_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, tl, dtdv_u_l, dtdu_v_l, dudv_t_l);
        DIFF_T_VU_G_CO2_T(vvt, vv, uv, tv, dtdv_u_v, dtdu_v_v, dudv_t_v);
        du=uv-ul;
        du_inv=1./du;
        du_inv2=du_inv*du_inv;
        dv=vv-vl;
        dv_inv=1./dv;
        dv_inv2=dv_inv*dv_inv;
        dpl=pl-ps;
        dpv=pv-ps;
        dtl=tl-ts;
        dtv=tv-ts;
        dx =(v-vl)*dv_inv-(u-ul)*du_inv;

        F[0] = dpl;
        F[1] = dpv;
        F[2] = dtl;
        F[3] = dtv;
        F[4] = dx;

        //matrix is stored in column major order (transposed for Fortran-like interfaces of dgetrf and dgetrs)
        J[0][0]=-2.*pst;    J[1][0]=dpdv_u_l;               J[2][0]=dpdu_v_l;               J[3][0]=0.;                 J[4][0]=0.;
        J[0][1]=-2.*pst;    J[1][1]=0.;                     J[2][1]=0.;                     J[3][1]=dpdv_u_v;           J[4][1]=dpdu_v_v;
        J[0][2]=-dtsdpt;    J[1][2]=dtdv_u_l;               J[2][2]=dtdu_v_l;               J[3][2]=0.;                 J[4][2]=0.;
        J[0][3]=-dtsdpt;    J[1][3]=0.;                     J[2][3]=0.;                     J[3][3]=dtdv_u_v;           J[4][3]=dtdu_v_v;
        J[0][4]=0.;         J[1][4]= ((v-vl)-dv)*dv_inv2;   J[2][4]=-((u-ul)-du)*du_inv2;   J[3][4]=-(v-vl)*dv_inv2;    J[4][4]= (u-ul)*du_inv2;

        dgetrf_(&N, &N, &J[0][0], &LDA, IPIV, &INFO);
        if(INFO) return I_ERR;
        else {
            dgetrs_(&TRANS, &N, &NRHS, &J[0][0], &LDA, &IPIV[0], &F[0], &LDB, &INFO);
            if(INFO) return I_ERR;
        }
        //assign new values
        pst=pst-F[0];
        vl =vl -F[1];
        ul =ul -F[2];
        vv =vv -F[3];
        uv =uv -F[4];
        ps=pst*pst;
        vvt=log(vv);
        if(icount++ > ITMAX) {
            return I_ERR;
        }
    } while (fabs(dpl/ps)>tol_ps || fabs(dpv/ps)>tol_ps ||
             fabs(dtl   )>tol_ts || fabs(dtv   )>tol_ts ||
             fabs(dx    )>tol_x);
    xvap=(u-ul)/(uv-ul);
    return I_OK;
}
/*
//Determine the two-phase equilibrium from pressure-, temperature-, Gibbs-energy-equilibria, rather than using a function ts(p).
//Note: * In order to do this consistently ts(p), ps(t), v'(u), u"(v) must be implemented this way as well.
//       (This is computationally more intensive than using ts(p).)
//      * g is computed from g=u+p*v-t*s in this routine.
SBTLAPI int __stdcall SAT_VU_SPL(double v, double u, double& ps, double& ts, double& xvap, double& vl, double& vv, double& vvt, double& ul, double& uv) throw()
{
    //declare variables and constants
    static const double tol_ps=1.e-10; //rel. tolerance in ps
    static const double tol_ts=1.e-8;  //abs. tolerance in ts
    static const double tol_g=1.e-8;
    static const double tol_x =1.e-8;  //abs. tolerance in x
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double pcx=7.33;
    //static const double tcx=303.84873305000252;
    static const double v1x=0.0018221754024023962;
    static const double v2x=0.0021592622018924488;
    static const double u1x=301.03552539317360;
    static const double u2x=316.77324527264307;
    //static const double h1x=u1x+pcx*v1x*1.e3;
    //static const double h2x=u2x+pcx*v2x*1.e3;
    //static const double s1x=1.3752027569977439;
    //static const double s2x=1.4351295210653687;
    static const double mx=(u2x-u1x)/(v2x-v1x);
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);
    double x1tmin, x1tmax, vlt;
    double dvdu_min, dvdu_max, K, dvdu_vt;
    double pl, dpdv_u_l, dpdu_v_l, dudv_p_l;
    double pv, dpdv_u_v, dpdu_v_v, dudv_p_v;
    double tl, dtdv_u_l, dtdu_v_l, dudv_t_l;
    double tv, dtdv_u_v, dtdu_v_v, dudv_t_v;
    double sl, dsdv_u_l, dsdu_v_l, dudv_s_l;
    double sv, dsdv_u_v, dsdu_v_v, dudv_s_v;
    double gl,gv,dgdv_u_l,dgdu_v_l,dgdv_u_v,dgdu_v_v;
    double dp,dt,dg;
    double du, du_inv, du_inv2;
    double dv, dv_inv, dv_inv2;
    double dx;

    //declare arrays
    double J[4][4], F[4];

    //declare and initialize flags and return values
    char TRANS='N';     //dgetrs does not transpose the matrix for LU-decomposition
    int INFO=0;         //return value (0-OK)
    int  LDA=4;         //leading dimension of J
    int  LDB=4;         //leading dimension of F (columns)
    int    N=4;         //number of rows in F
    int NRHS=1;         //number of right hand sides
    int IPIV[4];        //indices of pivot elements

    xvap=ERR_VAL;
 
    //vicinity of the critical point
    if(u>(mx*(v-v1x)+u1x)) {
        ps=pc;
        ts=tc;
        xvap=0.5;
        vl=vc;
        vv=vc;
        vvt=log(vc);
        ul=uc;
        uv=uc;
        return I_OK;
    }

    //initialize values
    double pst;
    ps =PS_VU_INI(v,u);
    pst=sqrt(ps);
    ts =TS_P_95_T(pst);
    vl=V1_T_AUX_CO2(ts);
    double lnp=log(ps);
    vvt=V2_P_AUX_CO2_T(lnp);
    vv=exp(vvt);
    ul=U1_T_AUX_CO2(ts);
    uv=U2_T_AUX_CO2(ts);

    int icount=0;
    do {
        //populate right hand side and jacobian matrix 
        DIFF_V_U_PMAX_CO2(ul, x1tmin, dvdu_min);
        DIFF_V1_U_SPL_CO2(ul, x1tmax, dvdu_max);
        K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
        vlt=(vl-x1tmin)*K+x1zmin;
        dvdu_vt=(vlt-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
        DIFF_P_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, pl, dpdv_u_l, dpdu_v_l, dudv_p_l);
        DIFF_P_VU_G_CO2_T(vvt, vv, uv, pv, dpdv_u_v, dpdu_v_v, dudv_p_v);
        DIFF_T_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, tl, dtdv_u_l, dtdu_v_l, dudv_t_l);
        DIFF_T_VU_G_CO2_T(vvt, vv, uv, tv, dtdv_u_v, dtdu_v_v, dudv_t_v);
        DIFF_S_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, sl, dsdv_u_l, dsdu_v_l, dudv_s_l);
        DIFF_S_VU_G_CO2_T(vvt, vv, uv, sv, dsdv_u_v, dsdu_v_v, dudv_s_v);
        gl=ul+pl*vl*1.e3-tl*sl;
        dgdv_u_l=(pl+dpdv_u_l*vl)*1.e3-tl*dsdv_u_l-dtdv_u_l*sl;
        dgdu_v_l=1.+dpdu_v_l*vl*1.e3-tl*dsdu_v_l-dtdu_v_l*sl;
        gv=uv+pv*vv*1.e3-tv*sv;
        dgdv_u_v=(pv+dpdv_u_v*vv)*1.e3-tv*dsdv_u_v-dtdv_u_v*sv;
        dgdu_v_v=1.+dpdu_v_v*vv*1.e3-tv*dsdu_v_v-dtdu_v_v*sv;
        du=uv-ul;
        du_inv=1./du;
        du_inv2=du_inv*du_inv;
        dv=vv-vl;
        dv_inv=1./dv;
        dv_inv2=dv_inv*dv_inv;
        dp=pl-pv;
        dt=tl-tv;
        dg=gl-gv;
        dx =(v-vl)*dv_inv-(u-ul)*du_inv;

        F[0] = dp;
        F[1] = dt;
        F[2] = dg;
        F[3] = dx;

        //matrix is stored in column major order (transposed for Fortran-like interfaces of dgetrf and dgetrs)
        J[0][0]=dpdv_u_l;               J[1][0]=dpdu_v_l;               J[2][0]=-dpdv_u_v;              J[3][0]=-dpdu_v_v;
        J[0][1]=dtdv_u_l;               J[1][1]=dtdu_v_l;               J[2][1]=-dtdv_u_v;              J[3][1]=-dtdu_v_v;
        J[0][2]=dgdv_u_l;               J[1][2]=dgdu_v_l;               J[2][2]=-dgdv_u_v;              J[3][2]=-dgdu_v_v;
        J[0][3]= ((v-vl)-dv)*dv_inv2;   J[1][3]=-((u-ul)-du)*du_inv2;   J[2][3]=-(v-vl)*dv_inv2;        J[3][3]= (u-ul)*du_inv2;

        dgetrf_(&N, &N, &J[0][0], &LDA, IPIV, &INFO);
        if(INFO) return I_ERR;
        else {
            dgetrs_(&TRANS, &N, &NRHS, &J[0][0], &LDA, &IPIV[0], &F[0], &LDB, &INFO);
            if(INFO) return I_ERR;
        }
        //assign new values
        vl =vl -F[0];
        ul =ul -F[1];
        vv =vv -F[2];
        uv =uv -F[3];
        vvt=log(vv);
        if(icount++ > ITMAX) {
            return I_ERR;
        }
    } while (fabs(dp/pl)>tol_ps || fabs(dt)>tol_ts ||
             fabs(dg   )>tol_g  || fabs(dx    )>tol_x);
    ps=pl;
    ts=tl;
    xvap=(u-ul)/(uv-ul);
    return I_OK;
}
*/
SBTLAPI void __stdcall DIFF_SAT_VU_SPL(double ps, double xvap, double vl, double vlt, double x1tmin, double x1tmax, double dvdu_vt, double K, double vv, double vvt, double ul, double uv, DERIV_TP& d_tp) throw()
{
    double p_,t_;
    double dpdv_u_l, dpdu_v_l, dudv_p_l, dpdv_u_v, dpdu_v_v, dudv_p_v;
    double dtdv_u_l, dtdu_v_l, dudv_t_l, dtdv_u_v, dtdu_v_v, dudv_t_v;

    //set check values
    d_tp.ps_=ps;
    d_tp.x_ =xvap;

    //calculation of dudv_p=dudv_t=dudv
    d_tp.dudv_pt=(uv-ul)/(vv-vl);

    //compute dtsdp
    DIFF_TS_P_CO2(ps, t_, d_tp.dtsdp);

    //compute dpdv_u, dpdu_v, dudv_p in both phases
    DIFF_P_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, p_, dpdv_u_l, dpdu_v_l, dudv_p_l);
    DIFF_P_VU_G_CO2_T(vvt, vv, uv, p_, dpdv_u_v, dpdu_v_v, dudv_p_v);

    //compute dtdv_u, dtdu_v, dudv_t in both phases
    DIFF_T_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, t_, dtdv_u_l, dtdu_v_l, dudv_t_l);
    DIFF_T_VU_G_CO2_T(vvt, vv, uv, t_, dtdv_u_v, dtdu_v_v, dudv_t_v);

    //compute dvldp and duldp
    d_tp.dvldp=(d_tp.dtsdp-dtdu_v_l/dpdu_v_l)/(dtdv_u_l+dtdu_v_l*dudv_p_l);
    d_tp.duldp=(d_tp.dtsdp-dtdv_u_l*d_tp.dvldp)/dtdu_v_l;

    //compute dvvdp and duvdp
    d_tp.dvvdp=(d_tp.dtsdp-dtdu_v_v/dpdu_v_v)/(dtdv_u_v+dtdu_v_v*dudv_p_v);
    d_tp.duvdp=(d_tp.dtsdp-dtdv_u_v*d_tp.dvvdp)/dtdu_v_v;

    //compute dxdp_u and dxdp_v
    d_tp.dxdp_u=(-d_tp.duldp-xvap*(d_tp.duvdp-d_tp.duldp))/(uv-ul);
    d_tp.dxdp_v=(-d_tp.dvldp-xvap*(d_tp.dvvdp-d_tp.dvldp))/(vv-vl);

    //compute dpdv_u
    d_tp.dpdv_u=1./(d_tp.dvldp+d_tp.dxdp_u*(vv-vl)+xvap*(d_tp.dvvdp-d_tp.dvldp));

    //calculate remaining differential dpdu_v
    d_tp.dpdu_v=-d_tp.dpdv_u/d_tp.dudv_pt;

    //calculate temperature derivatives
    d_tp.dtdv_u=d_tp.dpdv_u*d_tp.dtsdp;
    d_tp.dtdu_v=-d_tp.dtdv_u/d_tp.dudv_pt;
}
//
SBTLAPI int __stdcall SAT_U1_SPL(double u, double& ps, double& ts, double& vl) throw()
{
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double x1zmin=1.;
    static const double x1zmax=100.;

    static const double df_p=1.e-8;   //rel. deviation in p
    static const double df_t=1.e-8;   //abs. deviation in t

    double vt,x1tmin,x1tmax;
    double pst;
    double tx, dtdv_u, dtdu_v, dudv_t;
    double px, dpdv_u, dpdu_v, dudv_p;
    double dtsdpt;
    double den;

    if(fabs(u-uc)<1.e-3) {
       ps=pc;
       ts=tc;
       vl=vc;
       return I_OK;
    }

    //calculate initial guesses
    pst=PS_U_AUX_CO2_T(u);
    ps=pst*pst;
    //vl=V1_U_SPL_CO2(u);

    //scale v (iteration in transformed coordinates)
    x1tmin=V_U_PMAX_CO2(u);
    //x1tmax=V_U_SPNDL_L_CO2(u);
    x1tmax=V1_U_SPL_CO2(u);
    //vt=(vl-x1tmin)/(x1tmax-x1tmin)*(x1zmax-x1zmin)+x1zmin;
    vt=x1zmax;

    //newtons method
    double f_p=-1.,ps_inv=1.,f_t=-1.;
    int icount=0;
    while(fabs(f_p*ps_inv)>df_p || fabs(f_t)>df_t) {
        DIFF_P_VU_L_CO2_SC(vt, u, px, dpdv_u, dpdu_v, dudv_p);      // px, scaled derivatives
        DIFF_T_VU_L_CO2_SC(vt, u, tx, dtdv_u, dtdu_v, dudv_t);      // tx, scaled derivatives
        DIFF_TS_P_CO2_T(pst, ts, dtsdpt);
        f_p=px-ps;
        f_t=tx-ts;
        den=-2.*pst*dtdv_u+dpdv_u*dtsdpt;
        pst=pst+(-dtdv_u*f_p+f_t*dpdv_u)/den;
        vt =vt +(f_t*2.*pst-dtsdpt*f_p)/den;
        ps=pst*pst;
        ps_inv=1./ps;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    vl=(vt-x1zmin)/(x1zmax-x1zmin)*(x1tmax-x1tmin)+x1tmin;
    return I_OK;
}
//
SBTLAPI int __stdcall SAT_H1_SPL(double h, double& ps, double& ts, double& vl, double& ul, double& sl) throw()
{
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double hc=uc+pc*vc*1.e3;
    static const double sc=1.43362534304;
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);

    static const double df_p=1.e-8;   //rel. deviation in p
    static const double df_t=1.e-8;   //abs. deviation in t

    double pst;
    double v_u_pmax, dvdu_pmax, v_u_spndl, dvdu_spndl, K;
    double dvdu_vt;
    double dvdu_min,dvdu_max;
    double vt,x1tmin,x1tmax;
    double tx, dtdv_u, dtdu_v, dudv_t;
    double px, dpdv_u, dpdu_v, dudv_p;
    double dtsdpt;
    double dudv_pst,dudpst_v,dfpdv_pst,dfpdpst_v,dftdv_pst,dftdpst_v;
    double den;

    if(fabs(h-hc)<1.e-3) {
       ps=pc;
       ts=tc;
       vl=vc;
       ul=uc;
       sl=sc;
       return I_OK;
    }

    //calculate initial guesses
    pst=PS_H_AUX_CO2_T(h);
    ps=pst*pst;
    vl=V1_H_AUX_CO2(h);

    //newtons method
    double f_p=-1.,ps_inv=1.,f_t=-1.;
    int icount=0;
    while(fabs(f_p*ps_inv)>df_p || fabs(f_t)>df_t) {
        ul=h-ps*vl*1.e3;
        //scale v with scaling parameters for p and t (iteration in transformed coordinates)
        DIFF_V_U_PMAX_CO2(ul, v_u_pmax, dvdu_pmax);
        DIFF_V1_U_SPL_CO2(ul, v_u_spndl, dvdu_spndl);
        //scaling parameters for p and t
        x1tmin=v_u_pmax;
        x1tmax=v_u_spndl;
        dvdu_min=dvdu_pmax;
        dvdu_max=dvdu_spndl;
        K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
        vt=(vl-x1tmin)*K+x1zmin;
        dvdu_vt=(vt-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
        DIFF_P_VU_L_CO2_T(vt, x1tmin, x1tmax, dvdu_vt, K, ul, px, dpdv_u, dpdu_v, dudv_p);      // px, derivatives
        DIFF_T_VU_L_CO2_T(vt, x1tmin, x1tmax, dvdu_vt, K, ul, tx, dtdv_u, dtdu_v, dudv_t);      // tx, derivatives
        DIFF_TS_P_CO2_T(pst, ts, dtsdpt);
        dudv_pst=-ps*1.e3;
        dudpst_v=-vl*1.e3*2.*pst;
        //
        f_p=px-ps;
        f_t=tx-ts;
        //
        dfpdv_pst=dpdv_u+dpdu_v*dudv_pst;
        dfpdpst_v=dpdu_v*dudpst_v-2.*pst;
        dftdv_pst=dtdv_u+dtdu_v*dudv_pst;
        dftdpst_v=dtdu_v*dudpst_v-dtsdpt;
        //
        den=dfpdv_pst*dftdpst_v-dfpdpst_v*dftdv_pst;
        vl = vl+(-dftdpst_v*f_p+f_t*dfpdpst_v)/den;
        pst=pst+( f_p*dftdv_pst-dfpdv_pst*f_t)/den;
        ps=pst*pst;
        ps_inv=1./ps;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    sl=S_VU_L_CO2_T(vt, ul);

    return I_OK;
}
//
SBTLAPI int __stdcall SAT_V2_SPL_T(double vt, double& ps, double& ts, double& uv) throw()
{
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    //static const double vc=1./467.60000128174;
    static const double vtc=-6.1476132323390168;
    static const double uc=316.468709888;

    static const double df_p=1.e-8;   //rel. deviation in p
    static const double df_t=1.e-8;   //abs. deviation in t

    double pst;
    double tx, dtdv_u, dtdu_v, dudv_t;
    double px, dpdv_u, dpdu_v, dudv_p;
    double dtsdpt;
    double den;

    if(fabs(vt-vtc)<1.e-6) {
        ps=pc;
        ts=tc;
        uv=uc;
        return I_OK;
    }

    //calculate initial guesses
    pst=PS_V_AUX_CO2_T(vt);
    ps=pst*pst;
    uv=U2_V_AUX_CO2_T(vt);

    //newtons method
    double f_p=-1.,ps_inv=1.,f_t=-1.;
    int icount=0;
    while(fabs(f_p*ps_inv)>df_p || fabs(f_t)>df_t) {
        DIFF_P_VU_G_CO2_TT(vt, uv, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        DIFF_T_VU_G_CO2_TT(vt, uv, tx, dtdv_u, dtdu_v, dudv_t);      // tx, transformed derivatives
        DIFF_TS_P_CO2_T(pst, ts, dtsdpt);
        f_p=px-ps;
        f_t=tx-ts;
        den=-2.*pst*dtdu_v+dpdu_v*dtsdpt;
        pst=pst+(-dtdu_v*f_p+f_t*dpdu_v)/den;
        uv =uv +(f_t*2.*pst-dtsdpt*f_p)/den;
        ps=pst*pst;
        ps_inv=1./ps;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    return I_OK;
}
//
SBTLAPI int __stdcall SAT_S2_SPL(double s, double& ps, double& ts, double& vv, double& vvt, double& uv) throw()
{
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    //static const double vc=1./467.60000128174;
    static const double vtc=-6.1476132323390168;
    static const double uc=316.468709888;
    static const double sc=1.43362534304;

    static const double df_p=1.e-8;   //rel. deviation in p
    static const double df_t=1.e-8;   //abs. deviation in t
    static const double df_s=1.e-8;   //abs. deviation in s

    double pst;
    double tx, dtdv_u, dtdu_v, dudv_t;
    double px, dpdv_u, dpdu_v, dudv_p;
    double sx, dsdv_u, dsdu_v, dudv_s;
    double dtsdpt;
    double dfpdpst,dfpdvvt,dfpduv,dftdpst,dftdvvt,dftduv,dfsdvvt,dfsduv;
    double den;

    if(fabs(s-sc)<0.2) {
        ps=pc;
        ts=tc;
        vv=exp(vtc);
        vvt=vtc;
        uv=uc;
        return I_OK;
    }

    //calculate initial guesses
    pst=PS_S_AUX_CO2_T(s);
    ps=pst*pst;
    double lnp=log(ps);
    vvt=V2_P_AUX_CO2_T(lnp);
    uv=U2_V_AUX_CO2_T(vvt);

    //newtons method
    double f_p=-1.,ps_inv=1.,f_t=-1.,f_s=-1.;
    int icount=0;
    while(fabs(f_p*ps_inv)>df_p || fabs(f_t)>df_t || fabs(f_s)>df_s) {
        DIFF_P_VU_G_CO2_TT(vvt, uv, px, dpdv_u, dpdu_v, dudv_p);      // px, transformed derivatives
        DIFF_T_VU_G_CO2_TT(vvt, uv, tx, dtdv_u, dtdu_v, dudv_t);      // tx, transformed derivatives
        DIFF_S_VU_G_CO2_TT(vvt, uv, sx, dsdv_u, dsdu_v, dudv_s);      // sx, transformed derivatives
        DIFF_TS_P_CO2_T(pst, ts, dtsdpt);
        f_p=px-ps;
        f_t=tx-ts;
        f_s=sx-s;
        dfpdpst=-2.*pst;
        dfpdvvt=dpdv_u;
        dfpduv =dpdu_v;
        dftdpst=-dtsdpt;
        dftdvvt=dtdv_u;
        dftduv =dtdu_v;
        dfsdvvt=dsdv_u;
        dfsduv =dsdu_v;
        den=(dfpdpst*dftdvvt-dfpdvvt*dftdpst)*dfsduv+(dfpduv*dftdpst-dfpdpst*dftduv)*dfsdvvt;
        pst=pst+(-f_p*(dftdvvt*dfsduv-dftduv*dfsdvvt)-f_t*(dfpduv*dfsdvvt-dfpdvvt*dfsduv)-f_s*(dfpdvvt*dftduv-dfpduv*dftdvvt))/den;
        vvt=vvt+( f_p*dftdpst*dfsduv-f_t*dfpdpst*dfsduv-f_s*(dfpduv*dftdpst-dfpdpst*dftduv))/den;
        uv =uv +(-f_p*(dftdpst*dfsdvvt)+f_t*dfpdpst*dfsdvvt-f_s*(dfpdpst*dftdvvt-dfpdvvt*dftdpst))/den;
        ps=pst*pst;
        ps_inv=1./ps;
        if(icount++>ITMAX) {
            return I_ERR;
        }
    }
    vv=exp(vvt);
    return I_OK;
}
//
SBTLAPI int __stdcall SAT_HS_SPL(double h, double s, double& ps, double& ts, double& xvap, double& vl, double& vv, double& vvt, double& ul, double& uv) throw()
{
    //declare variables and constants
    static const double tol_ps=1.e-10; //rel. tolerance in ps
    static const double tol_ts=1.e-8;  //abs. tolerance in ts
    static const double tol_x =1.e-8;  //abs. tolerance in x
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double pcx=7.33;
    //static const double tcx=303.84873305000252;
    static const double v1x=0.0018221754024023962;
    static const double v2x=0.0021592622018924488;
    static const double u1x=301.03552539317360;
    static const double u2x=316.77324527264307;
    static const double h1x=u1x+pcx*v1x*1.e3;
    static const double h2x=u2x+pcx*v2x*1.e3;
    static const double s1x=1.3752027569977439;
    static const double s2x=1.4351295210653687;
    static const double mx=(h2x-h1x)/(s2x-s1x);
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);
    double x1tmin, x1tmax, vlt;
    double dvdu_min, dvdu_max, K, dvdu_vt;
    double pl, dpdv_u_l, dpdu_v_l, dudv_p_l;
    double pv, dpdv_u_v, dpdu_v_v, dudv_p_v;
    double tl, dtdv_u_l, dtdu_v_l, dudv_t_l;
    double tv, dtdv_u_v, dtdu_v_v, dudv_t_v;
    double sl, dsdv_u_l, dsdu_v_l, dudv_s_l;
    double sv, dsdv_u_v, dsdu_v_v, dudv_s_v;
    double hl,hv,dhldvl,dhldul,dhvdvv,dhvduv;
    double dtsdpt;
    double dh, dh_inv, dh_inv2;
    double ds, ds_inv, ds_inv2;
    double dpl,dpv,dtl,dtv,dx;
    double ddxdvl,ddxdvv,ddxdul,ddxduv;

    //declare arrays
    double J[5][5], F[5];

    //declare and initialize flags and return values
    char TRANS='N';         //dgetrs does not transpose the matrix for LU-decomposition
    int INFO=0;             //return value (0-OK)
    int  LDA=5;             //leading dimension of J
    int  LDB=5;             //leading dimension of F (columns)
    int    N=5;             //number of rows in F
    int NRHS=1;             //number of right hand sides
    int IPIV[5];            //indices of pivot elements

    xvap=ERR_VAL;
 
    //vicinity of the critical point
    if(h>(mx*(s-s1x)+h1x)) {
        ps=pc;
        ts=tc;
        xvap=0.5;
        vl=vc;
        vv=vc;
        vvt=log(vc);
        ul=uc;
        uv=uc;
        return I_OK;
    }

    //initialize values, use the transformed variable pst during the iteration
    double pst;
    pst=PS_SH_INI_CO2_T(s,h);
    ps=pst*pst;
    ts =TS_P_CO2_T(pst);
    vl=V1_T_AUX_CO2(ts);
    double lnp=log(ps);
    vvt=V2_P_AUX_CO2_T(lnp);
    vv=exp(vvt);
    ul=U1_T_AUX_CO2(ts);
    uv=U2_T_AUX_CO2(ts);

    int icount=0;
    do {
        //populate right hand side and jacobian matrix - use dpt derivatives, rather than dp 
        DIFF_V_U_PMAX_CO2(ul, x1tmin, dvdu_min);
        //scaling of s in region L is treated differently from p and t
        DIFF_V1_U_SPL_CO2(ul, x1tmax, dvdu_max);
        K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
        vlt=(vl-x1tmin)*K+x1zmin;
        dvdu_vt=(vlt-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
        DIFF_P_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, pl, dpdv_u_l, dpdu_v_l, dudv_p_l);
        DIFF_P_VU_G_CO2_T(vvt, vv, uv, pv, dpdv_u_v, dpdu_v_v, dudv_p_v);
        DIFF_TS_P_CO2_T(pst, ts, dtsdpt);
        DIFF_T_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, tl, dtdv_u_l, dtdu_v_l, dudv_t_l);
        DIFF_T_VU_G_CO2_T(vvt, vv, uv, tv, dtdv_u_v, dtdu_v_v, dudv_t_v);
        DIFF_S_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, sl, dsdv_u_l, dsdu_v_l, dudv_s_l);
        DIFF_S_VU_G_CO2_T(vvt, vv, uv, sv, dsdv_u_v, dsdu_v_v, dudv_s_v);
        //
        hl=ul+pl*vl*1.e3;
        hv=uv+pv*vv*1.e3;
        dhldvl=(dpdv_u_l*vl+pl)*1.e3;
        dhldul=1.+dpdu_v_l*vl*1.e3;
        dhvdvv=(dpdv_u_v*vv+pv)*1.e3;
        dhvduv=1.+dpdu_v_v*vv*1.e3;
        dh=hv-hl;
        dh_inv=1./dh;
        dh_inv2=dh_inv*dh_inv;
        ds=sv-sl;
        ds_inv=1./ds;
        ds_inv2=ds_inv*ds_inv;
        dpl=pl-ps;
        dpv=pv-ps;
        dtl=tl-ts;
        dtv=tv-ts;
        dx =(h-hl)*dh_inv-(s-sl)*ds_inv;

        ddxdvl=dhldvl*(h-hv)*dh_inv2-dsdv_u_l*(s-sv)*ds_inv2;
        ddxdvv=dhvdvv*(hl-h)*dh_inv2-dsdv_u_v*(sl-s)*ds_inv2;
        ddxdul=dhldul*(h-hv)*dh_inv2-dsdu_v_l*(s-sv)*ds_inv2;
        ddxduv=dhvduv*(hl-h)*dh_inv2-dsdu_v_v*(sl-s)*ds_inv2;

        F[0] = dpl;
        F[1] = dpv;
        F[2] = dtl;
        F[3] = dtv;
        F[4] = dx;

        //matrix is stored in column major order (transposed for Fortran-like interfaces of dgetrf and dgetrs)
        J[0][0]=-2.*pst;    J[1][0]=dpdv_u_l;               J[2][0]=dpdu_v_l;               J[3][0]=0.;                 J[4][0]=0.;
        J[0][1]=-2.*pst;    J[1][1]=0.;                     J[2][1]=0.;                     J[3][1]=dpdv_u_v;           J[4][1]=dpdu_v_v;
        J[0][2]=-dtsdpt;    J[1][2]=dtdv_u_l;               J[2][2]=dtdu_v_l;               J[3][2]=0.;                 J[4][2]=0.;
        J[0][3]=-dtsdpt;    J[1][3]=0.;                     J[2][3]=0.;                     J[3][3]=dtdv_u_v;           J[4][3]=dtdu_v_v;
        J[0][4]=0.;         J[1][4]=ddxdvl;                 J[2][4]=ddxdul;                 J[3][4]=ddxdvv;             J[4][4]=ddxduv;

        dgetrf_(&N, &N, &J[0][0], &LDA, IPIV, &INFO);
        if(INFO) return I_ERR;
        else {
            dgetrs_(&TRANS, &N, &NRHS, &J[0][0], &LDA, &IPIV[0], &F[0], &LDB, &INFO);
            if(INFO) return I_ERR;
        }
        //assign new values
        pst=pst-F[0];
        vl =vl -F[1];
        ul =ul -F[2];
        vv =vv -F[3];
        uv =uv -F[4];
        ps=pst*pst;
        vvt=log(vv);
        if(icount++ > ITMAX) {
            return I_ERR;
        }
    } while (fabs(dpl/ps)>tol_ps || fabs(dpv/ps)>tol_ps ||
             fabs(dtl   )>tol_ts || fabs(dtv   )>tol_ts ||
             fabs(dx    )>tol_x);
    xvap=(h-hl)/(hv-hl);
    return I_OK;
}
//
SBTLAPI int __stdcall SAT_VH_SPL(double v, double h, double& ps, double& ts, double& xvap, double& vl, double& vv, double& vvt, double& ul, double& uv) throw()
{
    //declare variables and constants
    static const double tol_ps=1.e-10; //rel. tolerance in ps
    static const double tol_ts=1.e-8;  //abs. tolerance in ts
    static const double tol_x =1.e-8;  //abs. tolerance in x
    static const double pc=7.37729837321;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double pcx=7.33;
    //static const double tcx=303.84873305000252;
    static const double v1x=0.0018221754024023962;
    static const double v2x=0.0021592622018924488;
    static const double u1x=301.03552539317360;
    static const double u2x=316.77324527264307;
    static const double h1x=u1x+pcx*v1x*1.e3;
    static const double h2x=u2x+pcx*v2x*1.e3;
    //static const double s1x=1.3752027569977439;
    //static const double s2x=1.4351295210653687;
    static const double mx=(h2x-h1x)/(v2x-v1x);
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);
    double x1tmin, x1tmax, vlt;
    double dvdu_min, dvdu_max, K, dvdu_vt;
    double pl, dpdv_u_l, dpdu_v_l, dudv_p_l;
    double pv, dpdv_u_v, dpdu_v_v, dudv_p_v;
    double tl, dtdv_u_l, dtdu_v_l, dudv_t_l;
    double tv, dtdv_u_v, dtdu_v_v, dudv_t_v;
    double hl,hv,dhldvl,dhldul,dhvdvv,dhvduv;
    double dtsdpt;
    double dv, dv_inv, dv_inv2;
    double dh, dh_inv, dh_inv2;
    double dpl,dpv,dtl,dtv,dx;
    double ddxdvl,ddxdvv,ddxdul,ddxduv;

    //declare arrays
    double J[5][5], F[5];

    //declare and initialize flags and return values
    char TRANS='N';         //dgetrs does not transpose the matrix for LU-decomposition
    int INFO=0;             //return value (0-OK)
    int  LDA=5;             //leading dimension of J
    int  LDB=5;             //leading dimension of F (columns)
    int    N=5;             //number of rows in F
    int NRHS=1;             //number of right hand sides
    int IPIV[5];            //indices of pivot elements

    xvap=ERR_VAL;
 
    //vicinity of the critical point
    if(h>(mx*(v-v1x)+h1x)) {
        ps=pc;
        ts=tc;
        xvap=0.5;
        vl=vc;
        vv=vc;
        vvt=log(vc);
        ul=uc;
        uv=uc;
        return I_OK;
    }

    //initialize values, use the transformed variable pst during the iteration
    double pst;
    double vt=log(v);
    pst=PS_VH_INI_CO2_T(vt,h);
    pst=pst*pst;
    ps=pst*pst;
    ts =TS_P_CO2_T(pst);
    vl=V1_T_AUX_CO2(ts);
    double lnp=log(ps);
    vvt=V2_P_AUX_CO2_T(lnp);
    vv=exp(vvt);
    ul=U1_T_AUX_CO2(ts);
    uv=U2_T_AUX_CO2(ts);

    int icount=0;
    do {
        //populate right hand side and jacobian matrix - use dpt derivatives, rather than dp 
        DIFF_V_U_PMAX_CO2(ul, x1tmin, dvdu_min);
        DIFF_V1_U_SPL_CO2(ul, x1tmax, dvdu_max);
        K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
        vlt=(vl-x1tmin)*K+x1zmin;
        dvdu_vt=(vlt-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
        DIFF_P_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, pl, dpdv_u_l, dpdu_v_l, dudv_p_l);
        DIFF_P_VU_G_CO2_T(vvt, vv, uv, pv, dpdv_u_v, dpdu_v_v, dudv_p_v);
        DIFF_TS_P_CO2_T(pst, ts, dtsdpt);
        DIFF_T_VU_L_CO2_T(vlt, x1tmin, x1tmax, dvdu_vt, K, ul, tl, dtdv_u_l, dtdu_v_l, dudv_t_l);
        DIFF_T_VU_G_CO2_T(vvt, vv, uv, tv, dtdv_u_v, dtdu_v_v, dudv_t_v);
        hl=ul+pl*vl*1.e3;
        hv=uv+pv*vv*1.e3;
        dhldvl=(dpdv_u_l*vl+pl)*1.e3;
        dhldul=1.+dpdu_v_l*vl*1.e3;
        dhvdvv=(dpdv_u_v*vv+pv)*1.e3;
        dhvduv=1.+dpdu_v_v*vv*1.e3;
        dh=hv-hl;
        dh_inv=1./dh;
        dh_inv2=dh_inv*dh_inv;
        dv=vv-vl;
        dv_inv=1./dv;
        dv_inv2=dv_inv*dv_inv;
        dpl=pl-ps;
        dpv=pv-ps;
        dtl=tl-ts;
        dtv=tv-ts;
        dx =(v-vl)*dv_inv-(h-hl)*dh_inv;

        ddxdvl=((v-vl)-dv)*dv_inv2-dhldvl*(h-hv)*dh_inv2;
        ddxdvv=-(v-vl)*dv_inv2+dhvdvv*(h-hl)*dh_inv2;
        ddxdul=dhldul*(hv-h)*dh_inv2;
        ddxduv=dhvduv*(h-hl)*dh_inv2;

        F[0] = dpl;
        F[1] = dpv;
        F[2] = dtl;
        F[3] = dtv;
        F[4] = dx;

        //matrix is stored in column major order (transposed for Fortran-like interfaces of dgetrf and dgetrs)
        J[0][0]=-2.*pst;    J[1][0]=dpdv_u_l;               J[2][0]=dpdu_v_l;               J[3][0]=0.;                 J[4][0]=0.;
        J[0][1]=-2.*pst;    J[1][1]=0.;                     J[2][1]=0.;                     J[3][1]=dpdv_u_v;           J[4][1]=dpdu_v_v;
        J[0][2]=-dtsdpt;    J[1][2]=dtdv_u_l;               J[2][2]=dtdu_v_l;               J[3][2]=0.;                 J[4][2]=0.;
        J[0][3]=-dtsdpt;    J[1][3]=0.;                     J[2][3]=0.;                     J[3][3]=dtdv_u_v;           J[4][3]=dtdu_v_v;
        J[0][4]=0.;         J[1][4]=ddxdvl;                 J[2][4]=ddxdul;                 J[3][4]=ddxdvv;             J[4][4]=ddxduv;

        dgetrf_(&N, &N, &J[0][0], &LDA, IPIV, &INFO);
        if(INFO) return I_ERR;
        else {
            dgetrs_(&TRANS, &N, &NRHS, &J[0][0], &LDA, &IPIV[0], &F[0], &LDB, &INFO);
            if(INFO) return I_ERR;
        }
        //assign new values
        pst=pst-F[0];
        vl =vl -F[1];
        ul =ul -F[2];
        vv =vv -F[3];
        uv =uv -F[4];
        ps=pst*pst;
        vvt=log(vv);
        if(icount++ > ITMAX) {
            return I_ERR;
        }
    } while (fabs(dpl/ps)>tol_ps || fabs(dpv/ps)>tol_ps ||
             fabs(dtl   )>tol_ts || fabs(dtv   )>tol_ts ||
             fabs(dx    )>tol_x);
    xvap=(h-hl)/(hv-hl);
    return I_OK;
}
//
SBTLAPI double __stdcall PS_VU_SPL(double v, double u) throw()
{
    double ps, ts, xvap, vl, vv, vvt, ul, uv;

    int ierr=SAT_VU_SPL(v,  u, ps, ts, xvap, vl, vv, vvt, ul, uv);
    if(ierr==0) {
        return ps;
    } else {
        return ERR_VAL;
    }
}
//
SBTLAPI double __stdcall TS_VU_SPL(double v, double u) throw()
{
    double ps, ts, xvap, vl, vv, vvt, ul, uv;

    int ierr=SAT_VU_SPL(v,  u, ps, ts, xvap, vl, vv, vvt, ul, uv);
    if(ierr==0) {
        return ts;
    } else {
        return ERR_VAL;
    }
}
//
SBTLAPI double __stdcall X_VU_SPL(double v, double u) throw()
{
    double ps, ts, xvap, vl, vv, vvt, ul, uv;

    int ierr=SAT_VU_SPL(v,  u, ps, ts, xvap, vl, vv, vvt, ul, uv);
    if(ierr==0) {
        return xvap;
    } else {
        return ERR_VAL;
    }
}