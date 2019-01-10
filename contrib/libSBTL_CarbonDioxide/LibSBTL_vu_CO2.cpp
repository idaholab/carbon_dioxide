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
// LibSBTL_vu_CO2
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_CO2.h"
#include "SBTL_call_conv.h"
#include "SBTL_func.h"
//
SBTLAPI void __stdcall ireg_vu_SBTLCO2(double v, double u, STR_vu_SBTL_CO2& r) throw()
{
    int ierr;
    double umax,vmin,vmax,ps_,ts_,u2_;
    double v_u_pmax,v1_u,x1tmin,x1tmax;
    double dvdu_pmax,dvdu_v1;
    double dvdu_min,dvdu_max,K,vs;
#ifndef SBTL_USE_C_AUX
    double d2vdu2_pmax,d2vdu2_v1;
    double d2vdu2_min,d2vdu2_max;
    double Z,N;
#endif
//
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double u2max= 398.81; // 398.805136516975;
    static const double upmaxtmax=1354.05666708;
    static const double vpmaxtmax=0.00311821748722;
    static const double upmintmax=1400.83548841;
    static const double vpmintmax=491.203124314;
    static const double vgmin=1.07029e-003;
    static const double vgmax=491.203;

    //saturation states at triple point (exactly)
    static const double u1tr=79.5960006450;
    static const double u2tr=392.775831067;
    static const double v1tr=0.000848563173213;
    static const double v2tr=0.0726697446683;
    static const double qtr=1./(u2tr-u1tr);
//
    static const double x1zmin=1.;
    static const double x1zmax=100.;
    static const double K2=1./(x1zmax-x1zmin);
//
    r.v_=v;
    r.u_=u;
    r.p_=ERR_VAL;   //reset for given (p,v), (p,h), or (p,s)
    r.h_=ERR_VAL;   //reset for given (p,h) or (h,s)
    r.s_=ERR_VAL;   //reset for given (p,s) or (h,s)

    r.ireg  =IREG_ERR;
    r.vls   =ERR_VAL;
    r.vt    =ERR_VAL;
    r.ps    =ERR_VAL;
    r.ts    =ERR_VAL;
    r.x     =ERR_VAL;
    r.v1    =ERR_VAL;
    r.v1s   =ERR_VAL;
    r.v2    =ERR_VAL;
    r.v2t   =ERR_VAL;
    r.u1    =ERR_VAL;
    r.u2    =ERR_VAL;
//
    if(v>vc) {
        if(u>uc) {
            if(u>u2max) {
                if(u<upmaxtmax) {
                    if(v>vpmaxtmax) {
                        if(v<v2tr) {
                            r.vt=log(v);
                            r.ireg=IREG_G;
                        } else {
                            if(v<=vpmintmax) {
                                r.vt=log(v);
                                r.ireg=IREG_G;
                            }
                        }
                    } else {
                        if(v>=vgmin) {
                            r.vt=log(v);
                            r.ireg=IREG_G;
                        }
                    }
                } else if(u<=upmintmax) {
                    r.vt=log(v);
                    umax=U_V_TMAX_AUX_CO2_T(r.vt);
                    if(u<=umax) {
                        r.ireg=IREG_G;
                    }
                }
            } else {
                if(v<=v2tr) {
                    r.vt=log(v);
                    u2_=U2_V_AUX_CO2_T(r.vt);
                    if(u>(u2_+0.1)) {
                        r.ireg=IREG_G;
                    } else if(u<(u2_-0.1)) {
                        vmax=v1tr+(u-u1tr)*qtr*(v2tr-v1tr);
                        if(v<=vmax) {
                            ierr=SAT_VU_SPL(v,u,r.ps,r.ts,r.x,r.v1,r.v2,r.v2t,r.u1,r.u2);
                            if(ierr==I_OK) {
                                r.ireg=IREG_TP;
#if defined(SBTL_USE_C_AUX)
                                //scaling curves:
                                DIFF_V_U_PMAX_CO2(r.u1, v_u_pmax, dvdu_pmax);   //dvdu_pmax: (dv/du)_p_max
                                DIFF_V1_U_SPL_CO2(r.u1, v1_u, dvdu_v1);    //dvdu_v1:   (dv/du)_v1
                                x1tmin=v_u_pmax;
                                x1tmax=v1_u;
                                r.dz_1.x1tmin=x1tmin;
                                r.dz_1.x1tmax=x1tmax;
                                dvdu_min=dvdu_pmax;
                                dvdu_max=dvdu_v1;
                                K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                                r.dz_1.K=K;
                                vs=(r.v1-x1tmin)*K+x1zmin;
                                r.v1s=vs;
                                r.dz_1.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
#else
                                //derivatives for second order partial derivatives
                                DIFF2_V_U_PMAX_CO2(r.u1, v_u_pmax, dvdu_pmax, d2vdu2_pmax);
                                DIFF2_V1_U_SPL_CO2(r.u1, v1_u, dvdu_v1, d2vdu2_v1);
                                //scaling curves:
                                x1tmin=v_u_pmax;
                                x1tmax=v1_u;
                                r.dz_1.x1tmin=x1tmin;
                                r.dz_1.x1tmax=x1tmax;
                                dvdu_min=dvdu_pmax;
                                d2vdu2_min=d2vdu2_pmax;
                                dvdu_max=dvdu_v1;
                                d2vdu2_max=d2vdu2_v1;
                                r.dz_1.dvdu_min=dvdu_min;
                                r.dz_1.d2vdu2_min=d2vdu2_min;
                                r.dz_1.dvdu_max=dvdu_max;
                                r.dz_1.d2vdu2_max=d2vdu2_max;
                                r.dz_1.d2vtdvdu=-(x1zmax-x1zmin)*(dvdu_max-dvdu_min)/((x1tmax-x1tmin)*(x1tmax-x1tmin));
                                K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                                r.dz_1.K=K;
                                vs=(r.v1-x1tmin)*K+x1zmin;
                                r.v1s=vs;
                                r.dz_1.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
                                Z=(-dvdu_min*(x1tmax-x1tmin)-(r.v1-x1tmin)*(dvdu_max-dvdu_min));
                                N=(x1tmax-x1tmin)*(x1tmax-x1tmin);
                                r.dz_1.dvtdu_v=(x1zmax-x1zmin)*Z/N;
                                r.dz_1.d2vtdu2_v=(x1zmax-x1zmin)*((-d2vdu2_min*(x1tmax-x1tmin)-dvdu_min*(dvdu_max-dvdu_min)-(-dvdu_min*(dvdu_max-dvdu_min)+(r.v1-x1tmin)*(d2vdu2_max-d2vdu2_min)))*N-Z*2.*(x1tmax-x1tmin)*(dvdu_max-dvdu_min))/(N*N);
#endif
                            }
                        }
                    } else {
                        ierr=SAT_V2_SPL_T(r.vt, ps_, ts_, u2_);
                        if(u>=u2_) {
                            r.ireg=IREG_G;
                        } else {
                            vmax=v1tr+(u-u1tr)*qtr*(v2tr-v1tr);
                            if(v<=vmax) {
                                ierr=SAT_VU_SPL(v,u,r.ps,r.ts,r.x,r.v1,r.v2,r.v2t,r.u1,r.u2);
                                if(ierr==I_OK) {
                                    r.ireg=IREG_TP;
#if defined(SBTL_USE_C_AUX)
                                    //scaling curves:
                                    DIFF_V_U_PMAX_CO2(r.u1, v_u_pmax, dvdu_pmax);   //dvdu_pmax:    (dv/du)_p_max
                                    DIFF_V1_U_SPL_CO2(r.u1, v1_u, dvdu_v1);         //dvdu_v1:      (dv/du)_v1
                                    x1tmin=v_u_pmax;
                                    x1tmax=v1_u;
                                    r.dz_1.x1tmin=x1tmin;
                                    r.dz_1.x1tmax=x1tmax;
                                    dvdu_min=dvdu_pmax;
                                    dvdu_max=dvdu_v1;
                                    K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                                    r.dz_1.K=K;
                                    vs=(r.v1-x1tmin)*K+x1zmin;
                                    r.v1s=vs;
                                    r.dz_1.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
#else
                                    //derivatives for second order partial derivatives
                                    DIFF2_V_U_PMAX_CO2(r.u1, v_u_pmax, dvdu_pmax, d2vdu2_pmax);
                                    DIFF2_V1_U_SPL_CO2(r.u1, v1_u, dvdu_v1, d2vdu2_v1);
                                    //scaling curves:
                                    x1tmin=v_u_pmax;
                                    x1tmax=v1_u;
                                    r.dz_1.x1tmin=x1tmin;
                                    r.dz_1.x1tmax=x1tmax;
                                    dvdu_min=dvdu_pmax;
                                    d2vdu2_min=d2vdu2_pmax;
                                    dvdu_max=dvdu_v1;
                                    d2vdu2_max=d2vdu2_v1;
                                    r.dz_1.dvdu_min=dvdu_min;
                                    r.dz_1.d2vdu2_min=d2vdu2_min;
                                    r.dz_1.dvdu_max=dvdu_max;
                                    r.dz_1.d2vdu2_max=d2vdu2_max;
                                    r.dz_1.d2vtdvdu=-(x1zmax-x1zmin)*(dvdu_max-dvdu_min)/((x1tmax-x1tmin)*(x1tmax-x1tmin));
                                    K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                                    r.dz_1.K=K;
                                    vs=(r.v1-x1tmin)*K+x1zmin;
                                    r.v1s=vs;
                                    r.dz_1.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
                                    Z=(-dvdu_min*(x1tmax-x1tmin)-(r.v1-x1tmin)*(dvdu_max-dvdu_min));
                                    N=(x1tmax-x1tmin)*(x1tmax-x1tmin);
                                    r.dz_1.dvtdu_v=(x1zmax-x1zmin)*Z/N;
                                    r.dz_1.d2vtdu2_v=(x1zmax-x1zmin)*((-d2vdu2_min*(x1tmax-x1tmin)-dvdu_min*(dvdu_max-dvdu_min)-(-dvdu_min*(dvdu_max-dvdu_min)+(r.v1-x1tmin)*(d2vdu2_max-d2vdu2_min)))*N-Z*2.*(x1tmax-x1tmin)*(dvdu_max-dvdu_min))/(N*N);
#endif
                                }
                            }
                        }
                    }
                } else {
                    if(u>u2tr) {
                        if(v<=vgmax) {
                            r.vt=log(v);
                            r.ireg=IREG_G;
                        }
                    }
                }
            }
        } else {
            vmax=v1tr+(u-u1tr)*qtr*(v2tr-v1tr);
            if(v<=vmax) {
                ierr=SAT_VU_SPL(v,u,r.ps,r.ts,r.x,r.v1,r.v2,r.v2t,r.u1,r.u2);
                if(ierr==I_OK) {
                    r.ireg=IREG_TP;
#if defined(SBTL_USE_C_AUX)
                    //scaling curves:
                    DIFF_V_U_PMAX_CO2(r.u1, v_u_pmax, dvdu_pmax);       //dvdu_pmax:    (dv/du)_p_max
                    DIFF_V1_U_SPL_CO2(r.u1, v1_u, dvdu_v1);             //dvdu_v1:      (dv/du)_v1
                    x1tmin=v_u_pmax;
                    x1tmax=v1_u;
                    r.dz_1.x1tmin=x1tmin;
                    r.dz_1.x1tmax=x1tmax;
                    dvdu_min=dvdu_pmax;
                    dvdu_max=dvdu_v1;
                    K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                    r.dz_1.K=K;
                    vs=(r.v1-x1tmin)*K+x1zmin;
                    r.v1s=vs;
                    r.dz_1.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
#else
                    //derivatives for second order partial derivatives
                    DIFF2_V_U_PMAX_CO2(r.u1, v_u_pmax, dvdu_pmax, d2vdu2_pmax);
                    DIFF2_V1_U_SPL_CO2(r.u1, v1_u, dvdu_v1, d2vdu2_v1);
                    //scaling curves:
                    x1tmin=v_u_pmax;
                    x1tmax=v1_u;
                    r.dz_1.x1tmin=x1tmin;
                    r.dz_1.x1tmax=x1tmax;
                    dvdu_min=dvdu_pmax;
                    d2vdu2_min=d2vdu2_pmax;
                    dvdu_max=dvdu_v1;
                    d2vdu2_max=d2vdu2_v1;
                    r.dz_1.dvdu_min=dvdu_min;
                    r.dz_1.d2vdu2_min=d2vdu2_min;
                    r.dz_1.dvdu_max=dvdu_max;
                    r.dz_1.d2vdu2_max=d2vdu2_max;
                    r.dz_1.d2vtdvdu=-(x1zmax-x1zmin)*(dvdu_max-dvdu_min)/((x1tmax-x1tmin)*(x1tmax-x1tmin));
                    K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                    r.dz_1.K=K;
                    vs=(r.v1-x1tmin)*K+x1zmin;
                    r.v1s=vs;
                    r.dz_1.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
                    Z=(-dvdu_min*(x1tmax-x1tmin)-(r.v1-x1tmin)*(dvdu_max-dvdu_min));
                    N=(x1tmax-x1tmin)*(x1tmax-x1tmin);
                    r.dz_1.dvtdu_v=(x1zmax-x1zmin)*Z/N;
                    r.dz_1.d2vtdu2_v=(x1zmax-x1zmin)*((-d2vdu2_min*(x1tmax-x1tmin)-dvdu_min*(dvdu_max-dvdu_min)-(-dvdu_min*(dvdu_max-dvdu_min)+(r.v1-x1tmin)*(d2vdu2_max-d2vdu2_min)))*N-Z*2.*(x1tmax-x1tmin)*(dvdu_max-dvdu_min))/(N*N);
#endif
                }
            }
        }
    } else {
        vmin=V_U_PMAX_AUX_CO2(u);
        if(v>=vmin*0.999) {
            if(u>=uc) {
                r.vt=log(v);
                r.ireg=IREG_G;
            } else {
                r.v1=V1_U_SPL_CO2(u);
                if(v>(r.v1+0.0005)) {
                    vmax=v1tr+(u-u1tr)*qtr*(v2tr-v1tr);
                    if(v<=vmax) {
                        ierr=SAT_VU_SPL(v,u,r.ps,r.ts,r.x,r.v1,r.v2,r.v2t,r.u1,r.u2);
                        if(ierr==I_OK) {
                            r.ireg=IREG_TP;
#if defined(SBTL_USE_C_AUX)
                            //scaling curves:
                            DIFF_V_U_PMAX_CO2(r.u1, v_u_pmax, dvdu_pmax);   //dvdu_pmax:    (dv/du)_p_max
                            DIFF_V1_U_SPL_CO2(r.u1, v1_u, dvdu_v1);         //dvdu_v1:      (dv/du)_v1
                            x1tmin=v_u_pmax;
                            x1tmax=v1_u;
                            r.dz_1.x1tmin=x1tmin;
                            r.dz_1.x1tmax=x1tmax;
                            dvdu_min=dvdu_pmax;
                            dvdu_max=dvdu_v1;
                            K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                            r.dz_1.K=K;
                            vs=(r.v1-x1tmin)*K+x1zmin;
                            r.v1s=vs;
                            r.dz_1.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
#else
                            //derivatives for second order partial derivatives
                            DIFF2_V_U_PMAX_CO2(r.u1, v_u_pmax, dvdu_pmax, d2vdu2_pmax);
                            DIFF2_V1_U_SPL_CO2(r.u1, v1_u, dvdu_v1, d2vdu2_v1);
                            //scaling curves:
                            x1tmin=v_u_pmax;
                            x1tmax=v1_u;
                            r.dz_1.x1tmin=x1tmin;
                            r.dz_1.x1tmax=x1tmax;
                            dvdu_min=dvdu_pmax;
                            d2vdu2_min=d2vdu2_pmax;
                            dvdu_max=dvdu_v1;
                            d2vdu2_max=d2vdu2_v1;
                            r.dz_1.dvdu_min=dvdu_min;
                            r.dz_1.d2vdu2_min=d2vdu2_min;
                            r.dz_1.dvdu_max=dvdu_max;
                            r.dz_1.d2vdu2_max=d2vdu2_max;
                            r.dz_1.d2vtdvdu=-(x1zmax-x1zmin)*(dvdu_max-dvdu_min)/((x1tmax-x1tmin)*(x1tmax-x1tmin));
                            K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                            r.dz_1.K=K;
                            vs=(r.v1-x1tmin)*K+x1zmin;
                            r.v1s=vs;
                            r.dz_1.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
                            Z=(-dvdu_min*(x1tmax-x1tmin)-(r.v1-x1tmin)*(dvdu_max-dvdu_min));
                            N=(x1tmax-x1tmin)*(x1tmax-x1tmin);
                            r.dz_1.dvtdu_v=(x1zmax-x1zmin)*Z/N;
                            r.dz_1.d2vtdu2_v=(x1zmax-x1zmin)*((-d2vdu2_min*(x1tmax-x1tmin)-dvdu_min*(dvdu_max-dvdu_min)-(-dvdu_min*(dvdu_max-dvdu_min)+(r.v1-x1tmin)*(d2vdu2_max-d2vdu2_min)))*N-Z*2.*(x1tmax-x1tmin)*(dvdu_max-dvdu_min))/(N*N);
#endif
                        }
                    }
                } else if(v<(r.v1-0.0005)) {
                    r.ireg=IREG_L;
#if defined(SBTL_USE_C_AUX)
                    //scaling curves:
                    DIFF_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax);  //dvdu_pmax:    (dv/du)_p_max
                    DIFF_V1_U_SPL_CO2(u, v1_u, dvdu_v1);        //dvdu_v1:      (dv/du)_v1
                    x1tmin=v_u_pmax;
                    x1tmax=v1_u;
                    r.dz_l.x1tmin=x1tmin;
                    r.dz_l.x1tmax=x1tmax;
                    dvdu_min=dvdu_pmax;
                    dvdu_max=dvdu_v1;
                    K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                    r.dz_l.K=K;
                    vs=(v-x1tmin)*K+x1zmin;
                    r.vls=vs;
                    r.dz_l.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
#else
                    //derivatives for second order partial derivatives
                    DIFF2_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax, d2vdu2_pmax);
                    DIFF2_V1_U_SPL_CO2(u, v1_u, dvdu_v1, d2vdu2_v1);
                    //scaling parameters for p and t
                    x1tmin=v_u_pmax;
                    x1tmax=v1_u;
                    r.dz_l.x1tmin=x1tmin;
                    r.dz_l.x1tmax=x1tmax;
                    dvdu_min=dvdu_pmax;
                    d2vdu2_min=d2vdu2_pmax;
                    dvdu_max=dvdu_v1;
                    d2vdu2_max=d2vdu2_v1;
                    r.dz_l.dvdu_min=dvdu_min;
                    r.dz_l.d2vdu2_min=d2vdu2_min;
                    r.dz_l.dvdu_max=dvdu_max;
                    r.dz_l.d2vdu2_max=d2vdu2_max;
                    r.dz_l.d2vtdvdu=-(x1zmax-x1zmin)*(dvdu_max-dvdu_min)/((x1tmax-x1tmin)*(x1tmax-x1tmin));
                    K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                    r.dz_l.K=K;
                    vs=(v-x1tmin)*K+x1zmin;
                    r.vls=vs;
                    r.dz_l.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
                    Z=(-dvdu_min*(x1tmax-x1tmin)-(v-x1tmin)*(dvdu_max-dvdu_min));
                    N=(x1tmax-x1tmin)*(x1tmax-x1tmin);
                    r.dz_l.dvtdu_v=(x1zmax-x1zmin)*Z/N;
                    r.dz_l.d2vtdu2_v=(x1zmax-x1zmin)*((-d2vdu2_min*(x1tmax-x1tmin)-dvdu_min*(dvdu_max-dvdu_min)-(-dvdu_min*(dvdu_max-dvdu_min)+(v-x1tmin)*(d2vdu2_max-d2vdu2_min)))*N-Z*2.*(x1tmax-x1tmin)*(dvdu_max-dvdu_min))/(N*N);
#endif
                } else {
                    ierr=SAT_U1_SPL(u, ps_, ts_, r.v1);
                    if(v>r.v1) {
                        vmax=v1tr+(u-u1tr)*qtr*(v2tr-v1tr);
                        if(v<=vmax) {
                            ierr=SAT_VU_SPL(v,u,r.ps,r.ts,r.x,r.v1,r.v2,r.v2t,r.u1,r.u2);
                            if(ierr==I_OK) {
                                r.ireg=IREG_TP;
#if defined(SBTL_USE_C_AUX)
                                //scaling curves:
                                DIFF_V_U_PMAX_CO2(r.u1, v_u_pmax, dvdu_pmax);   //dvdu_pmax:    (dv/du)_p_max
                                DIFF_V1_U_SPL_CO2(r.u1, v1_u, dvdu_v1);         //dvdu_v1:      (dv/du)_v1
                                x1tmin=v_u_pmax;
                                x1tmax=v1_u;
                                r.dz_1.x1tmin=x1tmin;
                                r.dz_1.x1tmax=x1tmax;
                                dvdu_min=dvdu_pmax;
                                dvdu_max=dvdu_v1;
                                K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                                r.dz_1.K=K;
                                vs=(r.v1-x1tmin)*K+x1zmin;
                                r.v1s=vs;
                                r.dz_1.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
#else
                                //derivatives for second order partial derivatives
                                DIFF2_V_U_PMAX_CO2(r.u1, v_u_pmax, dvdu_pmax, d2vdu2_pmax);
                                DIFF2_V1_U_SPL_CO2(r.u1, v1_u, dvdu_v1, d2vdu2_v1);
                                //scaling curves:
                                x1tmin=v_u_pmax;
                                x1tmax=v1_u;
                                r.dz_1.x1tmin=x1tmin;
                                r.dz_1.x1tmax=x1tmax;
                                dvdu_min=dvdu_pmax;
                                d2vdu2_min=d2vdu2_pmax;
                                dvdu_max=dvdu_v1;
                                d2vdu2_max=d2vdu2_v1;
                                r.dz_1.dvdu_min=dvdu_min;
                                r.dz_1.d2vdu2_min=d2vdu2_min;
                                r.dz_1.dvdu_max=dvdu_max;
                                r.dz_1.d2vdu2_max=d2vdu2_max;
                                r.dz_1.d2vtdvdu=-(x1zmax-x1zmin)*(dvdu_max-dvdu_min)/((x1tmax-x1tmin)*(x1tmax-x1tmin));
                                K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                                r.dz_1.K=K;
                                vs=(r.v1-x1tmin)*K+x1zmin;
                                r.v1s=vs;
                                r.dz_1.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
                                Z=(-dvdu_min*(x1tmax-x1tmin)-(r.v1-x1tmin)*(dvdu_max-dvdu_min));
                                N=(x1tmax-x1tmin)*(x1tmax-x1tmin);
                                r.dz_1.dvtdu_v=(x1zmax-x1zmin)*Z/N;
                                r.dz_1.d2vtdu2_v=(x1zmax-x1zmin)*((-d2vdu2_min*(x1tmax-x1tmin)-dvdu_min*(dvdu_max-dvdu_min)-(-dvdu_min*(dvdu_max-dvdu_min)+(r.v1-x1tmin)*(d2vdu2_max-d2vdu2_min)))*N-Z*2.*(x1tmax-x1tmin)*(dvdu_max-dvdu_min))/(N*N);
#endif
                            }
                        }
                    } else {
                        r.ireg=IREG_L;
#if defined(SBTL_USE_C_AUX)
                        //scaling curves:
                        DIFF_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax);  //dvdu_pmax:    (dv/du)_p_max
                        DIFF_V1_U_SPL_CO2(u, v1_u, dvdu_v1);        //dvdu_v1:      (dv/du)_v1
                        //offsets:
                            //p,t (incl. 1st derivatives)
                        x1tmin=v_u_pmax;
                        x1tmax=v1_u;
                        r.dz_l.x1tmin=x1tmin;
                        r.dz_l.x1tmax=x1tmax;
                        dvdu_min=dvdu_pmax;
                        dvdu_max=dvdu_v1;
                        K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                        r.dz_l.K=K;
                        vs=(v-x1tmin)*K+x1zmin;
                        r.vls=vs;
                        r.dz_l.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
#else
                        //derivatives for second order partial derivatives
                        DIFF2_V_U_PMAX_CO2(u, v_u_pmax, dvdu_pmax, d2vdu2_pmax);
                        DIFF2_V1_U_SPL_CO2(u, v1_u, dvdu_v1, d2vdu2_v1);
                        //scaling curves:
                        x1tmin=v_u_pmax;
                        x1tmax=v1_u;
                        r.dz_l.x1tmin=x1tmin;
                        r.dz_l.x1tmax=x1tmax;
                        dvdu_min=dvdu_pmax;
                        d2vdu2_min=d2vdu2_pmax;
                        dvdu_max=dvdu_v1;
                        d2vdu2_max=d2vdu2_v1;
                        r.dz_l.dvdu_min=dvdu_min;
                        r.dz_l.d2vdu2_min=d2vdu2_min;
                        r.dz_l.dvdu_max=dvdu_max;
                        r.dz_l.d2vdu2_max=d2vdu2_max;
                        r.dz_l.d2vtdvdu=-(x1zmax-x1zmin)*(dvdu_max-dvdu_min)/((x1tmax-x1tmin)*(x1tmax-x1tmin));
                        K=(x1zmax-x1zmin)/(x1tmax-x1tmin);
                        r.dz_l.K=K;
                        vs=(v-x1tmin)*K+x1zmin;
                        r.vls=vs;
                        r.dz_l.dvdu_vt=(vs-x1zmin)*K2*(dvdu_max-dvdu_min)+dvdu_min;
                        Z=(-dvdu_min*(x1tmax-x1tmin)-(v-x1tmin)*(dvdu_max-dvdu_min));
                        N=(x1tmax-x1tmin)*(x1tmax-x1tmin);
                        r.dz_l.dvtdu_v=(x1zmax-x1zmin)*Z/N;
                        r.dz_l.d2vtdu2_v=(x1zmax-x1zmin)*((-d2vdu2_min*(x1tmax-x1tmin)-dvdu_min*(dvdu_max-dvdu_min)-(-dvdu_min*(dvdu_max-dvdu_min)+(v-x1tmin)*(d2vdu2_max-d2vdu2_min)))*N-Z*2.*(x1tmax-x1tmin)*(dvdu_max-dvdu_min))/(N*N);
#endif
                    }
                }
            }
        }
    }
}
//
SBTLAPI int __stdcall ireg_pt_SBTLCO2(double p, double t) throw()
{
    double ts,tg,vtg;
    int ireg;
//
    static const double uc=316.468709888;
    static const double p0=0.0005;
    static const double pc=7.37729837321;
    static const double pmax=100.01;
    static const double t0=216.592;
    static const double tc=304.1282;
    static const double tmax=1305.;
//
    ireg=IREG_ERR;
//
    if(p<p0 || p>pmax || t<t0 || t>tmax) return IREG_ERR;
    if(p<pc) {
        if(t>tc) {
            ireg=IREG_G;
        } else {
            ts=TS_P_CO2(p);
            if(t>ts) {
                ireg=IREG_G;
            } else {
                ireg=IREG_L;
            }
        }
    } else {
        if(t<tc) {
            ireg=IREG_L;
        } else if(t>405.) {
            ireg=IREG_G;
        } else {
            tg=T_P_UC_AUX_CO2(p);
            if(t>(tg+1.e-2)) {
                ireg=IREG_G;
            } else if(t<(tg-1.e-2)) {
                ireg=IREG_L;
            } else {
                vtg=V_P_UC_CO2_T(p);
                tg=T_VU_G_CO2_T(vtg,uc);
                if(t>tg) {
                    ireg=IREG_G;
                } else {
                    ireg=IREG_L;
                }
            }
        }
    }
    return ireg;
}
//
SBTLAPI void __stdcall ireg_pv_SBTLCO2(double p, double v, STR_vu_SBTL_CO2& r) throw()
{
    double ts,vtg,vg,v1,u1,v2,v2t,u2,u;
    int ireg;
//
    static const double p0=0.0005;
    static const double pc=7.37729837321;
    static const double pmax=100.01;
//
    ireg=IREG_ERR;
//
    if(p<p0 || p>pmax) {
        ireg=IREG_ERR;
    } else if(p<pc) {
        ts=TS_P_CO2(p);
        PT_FLASH_L(p, ts, v1, u1);
        if(v<=v1) {
            ireg=IREG_L;
        } else {
            PT_FLASH_G(p, ts, v2, v2t, u2);
            if(v>=v2) {
                ireg=IREG_G;
            } else {
                ireg=IREG_TP;
            }
        }
    } else {
        vtg=V_P_UC_CO2_T(p);
        vg=exp(vtg);
        if(v>vg) {
            ireg=IREG_G;
        } else {
            ireg=IREG_L;
        }
    }

    switch(ireg) {
    case(IREG_L):
        u=U_VP_L_CO2(v, p);
        break;
    case(IREG_G):
        u=U_VP_G_CO2(v, p);
        break;
    case(IREG_TP):
        u=u1+(v-v1)/(v2-v1)*(u2-u1);
        break;
    default:
        r.ireg=IREG_ERR;
        return;
    }

    ireg_vu_SBTLCO2(v, u, r);
    r.p_=p;
}
//
SBTLAPI void __stdcall ireg_ps_SBTLCO2(double p, double s, STR_vu_SBTL_CO2& r) throw()
{
    double ts,vtg,sg,v1,u1,v2,v2t,u2,s1,s2,v,vt,u,x;
    double x1tmin,x1tmax,v1s;
    int ireg;
//
    static const double p0=0.0005;
    static const double pc=7.37729837321;
    static const double uc=316.468709888;
    static const double pmax=100.01;
    static const double x1zmin=1.;
    static const double x1zmax=100.;
//
    ireg=IREG_ERR;
//
    if(p<p0 || p>pmax) {
        ireg=IREG_ERR;
    } else if(p<pc) {
        ts=TS_P_CO2(p);
        PT_FLASH_L(p, ts, v1, u1);
        x1tmin=V_U_PMAX_CO2(u1);
        x1tmax=V1_U_SPL_CO2(u1);
        v1s=(v1-x1tmin)*(x1zmax-x1zmin)/(x1tmax-x1tmin)+x1zmin;
        s1=S_VU_L_CO2_T(v1s,u1);
        if(s<=s1) {
            ireg=IREG_L;
        } else {
            PT_FLASH_G(p, ts, v2, v2t, u2);
            s2=S_VU_G_CO2_T(v2t,u2);
            if(s>=s2) {
                ireg=IREG_G;
            } else {
                ireg=IREG_TP;
            }
        }
    } else {
        vtg=V_P_UC_CO2_T(p);
        sg=S_VU_G_CO2_T(vtg,uc);
        if(s>sg) {
            ireg=IREG_G;
        } else {
            ireg=IREG_L;
        }
    }

    switch(ireg) {
    case(IREG_L):
        PS_FLASH_L(p, s, v, u);
        break;
    case(IREG_G):
        PS_FLASH_G(p, s, v, vt, u);
        break;
    case(IREG_TP):
        x=(s-s1)/(s2-s1);
        v=v1+x*(v2-v1);
        u=u1+x*(u2-u1);
        break;
    default:
        r.ireg=IREG_ERR;
        return;
    }

    ireg_vu_SBTLCO2(v, u, r);
    r.p_=p;
    r.s_=s;
}
//
SBTLAPI void __stdcall ireg_ph_SBTLCO2(double p, double h, STR_vu_SBTL_CO2& r) throw()
{
    double ts,vtg,vg,hg,v1,u1,v2,v2t,u2,h1,h2,v,vt,u,x;
    int ireg;
//
    static const double p0=0.0005;
    static const double pc=7.37729837321;
    static const double uc=316.468709888;
    static const double pmax=100.01;
//
    ireg=IREG_ERR;
//
    if(p<p0 || p>pmax) {
        ireg=IREG_ERR;
    } else if(p<pc) {
        ts=TS_P_CO2(p);
        PT_FLASH_L(p, ts, v1, u1);
        h1=u1+p*v1*1.e3;
        if(h<=h1) {
            ireg=IREG_L;
        } else {
            PT_FLASH_G(p, ts, v2, v2t, u2);
            h2=u2+p*v2*1.e3;
            if(h>=h2) {
                ireg=IREG_G;
            } else {
                ireg=IREG_TP;
            }
        }
    } else {
        vtg=V_P_UC_CO2_T(p);
        vg=exp(vtg);
        hg=uc+p*vg*1.e3;
        if(h>hg) {
            ireg=IREG_G;
        } else {
            ireg=IREG_L;
        }
    }

    switch(ireg) {
    case(IREG_L):
        PH_FLASH_L(p, h, v, u);
        break;
    case(IREG_G):
        PH_FLASH_G(p, h, v, vt, u);
        break;
    case(IREG_TP):
        x=(h-h1)/(h2-h1);
        v=v1+x*(v2-v1);
        u=u1+x*(u2-u1);
        break;
    default:
        r.ireg=IREG_ERR;
        return;
    }

    ireg_vu_SBTLCO2(v, u, r);
    r.p_=p;
    r.h_=h;
}
//
SBTLAPI void __stdcall ireg_hs_SBTLCO2(double h, double s, STR_vu_SBTL_CO2& r) throw()
{
    double smin,smax,hmax,s1,hg;
    double ps,ts,x,v1,u1,v2,v2t,u2,h2;
    double v,vt,u;
    int ireg;
//
    static const double pc=7.37729837321;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double hc=uc+pc*vc*1.e3;
    static const double sc=1.43362534304;
    //static const double hl_min=27.5;      //lower boundary of h in L
    static const double h2_max=437.1;       //maximum of h" (round up)
    static const double h1tr=80.0355261118;
    static const double h2tr=430.416167646;
    static const double s1tr=0.521319785219;
    static const double s2tr=2.13901869162;
    static const double qtr=1./(h2tr-h1tr);
    //static const double upmintmax=1400.83548841;
    static const double hpmintmax=1646.43705056;
    static const double spmintmax=5.33648159436;
    //static const double upmaxtmax=1354.05666708;
    static const double hpmaxtmax=1665.87841580;
    static const double spmaxtmax=2.99766140114;
//
    ireg=IREG_ERR;
//
    if(h>h1tr) {
        smin=S_H_PMAX_AUX_CO2(h);
        if(s>=smin) {
            if(h<hc) {
                s1=S1_H_AUX_CO2(h);   //auxiliary function
                if(s<(s1-1.e-3)) {
                    ireg=IREG_L;
                } else if(s>(s1+1.e-3)) {
                    smax=s1tr+(h-h1tr)*qtr*(h2tr-h1tr);
                    if(s<=smax) {
                        ireg=IREG_TP;
                    }
                } else {
                    SAT_H1_SPL(h, ps, ts, v1, u1, s1);
                    if(s<s1) {
                        ireg=IREG_L;
                    } else {
                        ireg=IREG_TP;
                    }
                }
            } else if(h<h2_max) {
                if(s<sc) {
                    hg=H_S_UC_CO2(s);
                    if(h<hg) {
                        ireg=IREG_L;
                    } else {
                        ireg=IREG_G;
                    }
                } else if(s<s2tr) {
                    h2=H2_S_AUX_CO2(s);
                    if(h<(h2-0.01)) {
                        smax=s1tr+(h-h1tr)*qtr*(h2tr-h1tr);
                        if(s<=smax) {
                            ireg=IREG_TP;
                        }
                    } else if(h>(h2+0.01)) {
                        ireg=IREG_G;
                    } else {
                        SAT_S2_SPL(s, ps, ts, v2, v2t, u2);
                        h2=u2+ps*v2*1.e3;
                        if(h>=h2) {
                            ireg=IREG_G;
                        } else {
                            ireg=IREG_TP;
                        }
                    }
                } else if(s<=spmintmax) {
                    smax=S_H_PMIN_AUX_CO2(h);
                    if(s<=smax) {
                        ireg=IREG_G;
                    }
                }
            } else if(h<hpmaxtmax) {
                smax=S_H_PMIN_AUX_CO2(h);
                if(s<=smax) {
                    ireg=IREG_G;
                }
            } else if(h<hpmintmax) {
                if(s>=spmaxtmax && s<=spmintmax) {
                    hmax=H_S_TMAX_AUX_CO2(s);
                    if(h<=hmax) {
                        ireg=IREG_G;
                    }
                }
            }
        }
    }

    switch(ireg) {
    case(IREG_L):
        HS_FLASH_L(h, s, v, u);
        break;
    case(IREG_G):
        HS_FLASH_G(h, s, v, vt, u);
        break;
    case(IREG_TP):
        SAT_HS_SPL(h, s, ps, ts, x, v1, v2, v2t, u1, u2);
        v=v1+x*(v2-v1);
        u=u1+x*(u2-u1);
        break;
    default:
        r.ireg=IREG_ERR;
        return;
    }

    ireg_vu_SBTLCO2(v, u, r);
    r.h_=h;
    r.s_=s;
}
//
SBTLAPI void __stdcall ireg_vh_SBTLCO2(double v, double h, STR_vu_SBTL_CO2& r) throw()
{
    double vmin,vmax,hmax,hg;
    double s1,h1,h2;
    double ps,ts,x,v1,v2,v2t,u1,u2;
    double vt,u;
    int ireg;
//
    static const double pc=7.37729837321;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;
    static const double hc=uc+pc*vc*1.e3;
    static const double h2_max=437.1;       //maximum of h" (round up)
    static const double h1tr=80.0355261118;
    static const double h2tr=430.416167646;
    static const double v1tr=0.000848563173213;
    static const double v2tr=0.0726697446683;
    static const double qtr=1./(h2tr-h1tr);
    static const double hpmaxtmax=1665.87841580;
    static const double vpmaxtmax=0.00311821748722;
    static const double hpmintmax=1646.43705056;
    static const double vpmintmax=491.203124314;
//
    ireg=IREG_ERR;
    vt=ERR_VAL;
//
    if(h>h1tr) {
        vmin=V_H_PMAX_AUX_CO2(h);
        if(v>=vmin) {
            if(h<hc) {
                v1=V1_H_AUX_CO2(h);   //auxiliary function
                if(v<(v1*0.9999)) {
                    ireg=IREG_L;
                } else if(v>(v1*1.0001)) {
                    vmax=v1tr+(h-h1tr)*qtr*(h2tr-h1tr);
                    if(v<=vmax) {
                        ireg=IREG_TP;
                    }
                } else {
                    SAT_H1_SPL(h, ps, ts, v1, u1, s1);
                    if(v<v1) {
                        ireg=IREG_L;
                    } else {
                        ireg=IREG_TP;
                    }
                }
            } else if(h<h2_max) {
                if(v<vc) {
                    hg=H_V_UC_CO2(v);
                    if(h<hg) {
                        ireg=IREG_L;
                    } else {
                        ireg=IREG_G;
                    }
                } else if(v<v2tr) {
                    vt=log(v);
                    h2=H2_V_AUX_CO2_T(vt);
                    if(h<(h2-0.01)) {
                        vmax=v1tr+(h-h1tr)*qtr*(h2tr-h1tr);
                        if(v<=vmax) {
                            ireg=IREG_TP;
                        }
                    } else if(h>(h2+0.01)) {
                        ireg=IREG_G;
                    } else {
                        SAT_V2_SPL_T(vt, ps, ts, u2);
                        h2=u2+ps*v*1.e3;
                        if(h>=h2) {
                            ireg=IREG_G;
                        } else {
                            ireg=IREG_TP;
                        }
                    }
                } else if(v<=vpmintmax) {
                    vmax=V_H_PMIN_AUX_CO2(h);
                    if(v<=vmax) {
                        ireg=IREG_G;
                    }
                }
            } else if(h<hpmaxtmax) {
                vmax=V_H_PMIN_AUX_CO2(h);
                if(v<=vmax) {
                    ireg=IREG_G;
                }
            } else if(h<hpmintmax) {
                if(v>=vpmaxtmax && v<=vpmintmax) {
                    vt=log(v);
                    hmax=H_V_TMAX_AUX_CO2_T(vt);
                    if(h<=hmax) {
                        ireg=IREG_G;
                    }
                }
            }
        }
    }

    switch(ireg) {
    case(IREG_L):
        VH_FLASH_L(v, h, u);
        break;
    case(IREG_G):
        if(vt==ERR_VAL) vt=log(v);
        VH_FLASH_G_T(v, vt, h, u);
        break;
    case(IREG_TP):
        SAT_VH_SPL(v, h, ps, ts, x, v1, v2, v2t, u1, u2);
        h1=u1+ps*v1*1.e3;
        h2=u2+ps*v2*1.e3;
        x=(h-h1)/(h2-h1);
        u=u1+x*(u2-u1);
        break;
    default:
        r.ireg=IREG_ERR;
        return;
    }

    ireg_vu_SBTLCO2(v, u, r);
    r.v_=v;
    r.h_=h;
}
//
SBTLAPI int __stdcall flash_vu_SBTLCO2(double v, double u, int& ireg, double& p, double& t, double& x, double& alpha_l, double& s, double& cp, double& cv, double& w, double& eta, double& lambda) throw()
{
    STR_vu_SBTL_CO2 str;
    double s1,s2;
    double alpha_v,K1,K2,rho_m,K_m;
    double cp1,cp2,cv1,cv2,w1,w2,eta1,eta2,lambda1,lambda2;
//
    p      =ERR_VAL;
    t      =ERR_VAL;
    x      =ERR_VAL;
    alpha_l=ERR_VAL;
    s      =ERR_VAL;
    cp     =ERR_VAL;
    cv     =ERR_VAL;
    w      =ERR_VAL;
    eta    =ERR_VAL;
    lambda =ERR_VAL;
//
    ireg_vu_SBTLCO2(v, u, str);
    ireg=str.ireg;
//
    switch(ireg) {
    case(IREG_L):
        p     =P_VU_L_CO2_T(str.vls,u);
        t     =T_VU_L_CO2_T(str.vls,u);
        s     =S_VU_L_CO2_T(str.vls ,u);
#ifdef SBTL_USE_C_AUX
        cp    =CP_VU_L_CO2_T(str.vls,u);
        cv    =CV_VU_L_CO2_T(str.vls,u);
#else
        //p,t-scaling here, because cp and cv are derived from those
        cp    =CP_VU_L_CO2_T(str.vls,str.dz_l.x1tmin,str.dz_l.x1tmax,str.dz_l.dvdu_vt,str.dz_l.K,u);
        cv    =CV_VU_L_CO2_T(str.vls,str.dz_l.x1tmin,str.dz_l.x1tmax,str.dz_l.dvdu_vt,str.dz_l.K,u);
#endif
        w     =W_VU_L_CO2_T(str.vls,u);
        eta   =ETA_VU_L_CO2_T(str.vls,u);
        lambda=LAMBDA_VU_L_CO2_T(str.vls,u);
        break;
    case(IREG_G):
        p     =P_VU_G_CO2_T(str.vt,u);
        t     =T_VU_G_CO2_T(str.vt,u);
        s     =S_VU_G_CO2_T(str.vt,u);
#ifdef SBTL_USE_C_AUX
        cp    =CP_VU_G_CO2_T(str.vt,u);
        cv    =CV_VU_G_CO2_T(str.vt,u);
#else
        cp    =CP_VU_G_CO2_T(str.vt,v,u);
        cv    =CV_VU_G_CO2_T(str.vt,v,u);
#endif
        w     =W_VU_G_CO2_T(str.vt,u);
        eta   =ETA_VU_G_CO2_T(str.vt,u);
        lambda=LAMBDA_VU_G_CO2_T(str.vt,u);
        break;
    case(IREG_TP):
        p=str.ps;
        t=str.ts;
        x=str.x;
        alpha_l=(1.-x)*str.v1/v;
        s1    =S_VU_L_CO2_T(str.v1s,str.u1);
        s2    =S_VU_G_CO2_T(str.v2t,str.u2);
        s=s1+x*(s2-s1);
        //cp is infinity in the two phase region (mass average is calculated)
        //cv: the derivative (du/dt)_v is calculable in the two-phase region but its slope is discontinuous at dew and bubble curve (mass average is calculated)
#ifdef SBTL_USE_C_AUX
        cp1   =CP_VU_L_CO2_T(str.v1s,str.u1);
        cv1   =CV_VU_L_CO2_T(str.v1s,str.u1);
        cp2   =CP_VU_G_CO2_T(str.v2t,str.u1);
        cv2   =CV_VU_G_CO2_T(str.v2t,str.u2);
#else
        //p,t-scaling here, because cp and cv are derived from those
        cp1   =CP_VU_L_CO2_T(str.v1s,str.dz_1.x1tmin,str.dz_1.x1tmax,str.dz_1.dvdu_vt,str.dz_1.K,str.u1);
        cv1   =CV_VU_L_CO2_T(str.v1s,str.dz_1.x1tmin,str.dz_1.x1tmax,str.dz_1.dvdu_vt,str.dz_1.K,str.u1);
        cp2   =CP_VU_G_CO2_T(str.v2t,str.v2,str.u1);
        cv2   =CV_VU_G_CO2_T(str.v2t,str.v2,str.u2);
#endif
        cp=cp1+x*(cp2-cp1);
        cv=cv1+x*(cv2-cv1);
        //Wood's equation for speed of sound
        alpha_v=1.-alpha_l;                         //alpha_v=x*v2/v;
        w1     =W_VU_L_CO2_T(str.v1s,str.u1);
        w2     =W_VU_G_CO2_T(str.v2t  ,str.u2);
        K1     =str.v1/(w1*w1);
        K2     =str.v2/(w2*w2);
        rho_m  =alpha_v/str.v2+(1.-alpha_v)/str.v1;   //density of mixture
        K_m    =alpha_v*K2+(1.-alpha_v)*K1;         //comstressibility of mixture
        w      = 1./sqrt(K_m*rho_m);
        //viscosity and thermal conductivity depend on the vapor-volume fraction
        eta1   =ETA_VU_L_CO2_T(str.v1s,str.u1);
        eta2   =ETA_VU_G_CO2_T(str.v2t,str.u2);
        eta    =eta1+alpha_v*(eta2-eta1);
        lambda1=LAMBDA_VU_L_CO2_T(str.v1s,str.u1);
        lambda2=LAMBDA_VU_G_CO2_T(str.v2t  ,str.u2);
        lambda =lambda1+alpha_v*(lambda2-lambda1);
    break;
    default:
    return I_ERR;
    }
    return I_OK;
}
//
SBTLAPI int __stdcall flash_deriv_vu_SBTLCO2(double v, double u, int& ireg, double& p, double& t, double& x, double& alpha_l, double& s, double& cp, double& cv, double& w, double& eta, double& lambda,
                                            double& dpdv_u, double& dpdu_v, double& dudv_p,
                                            double& dtdv_u, double& dtdu_v, double& dudv_t,
                                            double& dsdv_u, double& dsdu_v, double& dudv_s,
                                            double& dcpdv_u, double& dcpdu_v, double& dudv_cp,
                                            double& dcvdv_u, double& dcvdu_v, double& dudv_cv,
                                            double& dwdv_u, double& dwdu_v, double& dudv_w,
                                            double& detadv_u, double& detadu_v, double& dudv_eta,
                                            double& dlambdadv_u, double& dlambdadu_v, double& dudv_lambda) throw()
{
    STR_vu_SBTL_CO2 str;
    double s1,s2;
    double alpha_v,F1,F2,rho_m,F_m;
    double cp1,cp2,cv1,cv2,w1,w2,eta1,eta2,lambda1,lambda2;
    double dsdv_u_1, dsdu_v_1, dudv_s_1, dsdv_u_2, dsdu_v_2, dudv_s_2;
    double ds1dps,ds2dps,dsdps_u,dsdu_ps,dsdv_ps,dsdps_v;
//
    p      =ERR_VAL;
    t      =ERR_VAL;
    x      =ERR_VAL;
    alpha_l=ERR_VAL;
    s      =ERR_VAL;
    cp     =ERR_VAL;
    cv     =ERR_VAL;
    w      =ERR_VAL;
    eta    =ERR_VAL;
    lambda =ERR_VAL;
    dpdv_u=ERR_VAL;         dpdu_v=ERR_VAL;         dudv_p=ERR_VAL;
    dtdv_u=ERR_VAL;         dtdu_v=ERR_VAL;         dudv_t=ERR_VAL;
    dsdv_u=ERR_VAL;         dsdu_v=ERR_VAL;         dudv_s=ERR_VAL;
    dcpdv_u=ERR_VAL;        dcpdu_v=ERR_VAL;        dudv_cp=ERR_VAL;
    dcvdv_u=ERR_VAL;        dcvdu_v=ERR_VAL;        dudv_cv=ERR_VAL;
    dwdv_u=ERR_VAL;         dwdu_v=ERR_VAL;         dudv_w=ERR_VAL;
    detadv_u=ERR_VAL;       detadu_v=ERR_VAL;       dudv_eta=ERR_VAL;
    dlambdadv_u=ERR_VAL;    dlambdadu_v=ERR_VAL;    dudv_lambda=ERR_VAL;
//
    ireg_vu_SBTLCO2(v, u, str);
    ireg=str.ireg;
//
    switch(ireg) {
    case(IREG_L):
        //p, t, and their derivatives
        DIFF_P_VU_L_CO2_T(str.vls, str.dz_l.x1tmin, str.dz_l.x1tmax, str.dz_l.dvdu_vt, str.dz_l.K, u, p, dpdv_u, dpdu_v, dudv_p);
        DIFF_T_VU_L_CO2_T(str.vls, str.dz_l.x1tmin, str.dz_l.x1tmax, str.dz_l.dvdu_vt, str.dz_l.K, u, t, dtdv_u, dtdu_v, dudv_t);
        DIFF_S_VU_L_CO2_T(str.vls, str.dz_l.x1tmin, str.dz_l.x1tmax, str.dz_l.dvdu_vt, str.dz_l.K, u, s, dsdv_u, dsdu_v, dudv_s);
#ifndef SBTL_USE_C_AUX
        DIFF_CV_VU_L_CO2_T(str.vls, str.dz_l.x1tmin, str.dz_l.x1tmax, str.dz_l.dvtdu_v, str.dz_l.dvdu_vt, str.dz_l.K, str.dz_l.d2vtdu2_v, str.dz_l.d2vtdvdu, u, cv, dcvdv_u, dcvdu_v, dudv_cv);
        DIFF_CP_VU_L_CO2_T(str.vls, str.dz_l.x1tmin, str.dz_l.x1tmax, str.dz_l.dvtdu_v, str.dz_l.dvdu_vt, str.dz_l.K, str.dz_l.d2vtdu2_v, str.dz_l.d2vtdvdu, u, cp, dcpdv_u, dcpdu_v, dudv_cp);
#else
        DIFF_CP_VU_L_CO2_T(str.vls, str.dz_l.x1tmin, str.dz_l.x1tmax, str.dz_l.dvdu_vt, str.dz_l.K, u, cp, dcpdv_u, dcpdu_v, dudv_cp);
        DIFF_CV_VU_L_CO2_T(str.vls, str.dz_l.x1tmin, str.dz_l.x1tmax, str.dz_l.dvdu_vt, str.dz_l.K, u, cv, dcvdv_u, dcvdu_v, dudv_cv);
#endif
        DIFF_W_VU_L_CO2_T(str.vls, str.dz_l.x1tmin, str.dz_l.x1tmax, str.dz_l.dvdu_vt, str.dz_l.K, u, w, dwdv_u, dwdu_v, dudv_w);
        DIFF_ETA_VU_L_CO2_T(str.vls, str.dz_l.x1tmin, str.dz_l.x1tmax, str.dz_l.dvdu_vt, str.dz_l.K, u, eta, detadv_u, detadu_v, dudv_eta);
        DIFF_LAMBDA_VU_L_CO2_T(str.vls, str.dz_l.x1tmin, str.dz_l.x1tmax, str.dz_l.dvdu_vt, str.dz_l.K, u, lambda, dlambdadv_u, dlambdadu_v, dudv_lambda);
        break;
    case(IREG_G):
        DIFF_P_VU_G_CO2_T(str.vt, v, u, p, dpdv_u, dpdu_v, dudv_p);
        DIFF_T_VU_G_CO2_T(str.vt, v, u, t, dtdv_u, dtdu_v, dudv_t);
        DIFF_S_VU_G_CO2_T(str.vt, v, u, s, dsdv_u, dsdu_v, dudv_s);
        DIFF_CP_VU_G_CO2_T(str.vt, v, u, cp, dcpdv_u, dcpdu_v, dudv_cp);
        DIFF_CV_VU_G_CO2_T(str.vt, v, u, cv, dcvdv_u, dcvdu_v, dudv_cv);
        DIFF_W_VU_G_CO2_T(str.vt, v, u, w, dwdv_u, dwdu_v, dudv_w);
        DIFF_ETA_VU_G_CO2_T(str.vt, v, u, eta, detadv_u, detadu_v, dudv_eta);
        DIFF_LAMBDA_VU_G_CO2_T(str.vt, v, u, lambda, dlambdadv_u, dlambdadu_v, dudv_lambda);
        break;
    case(IREG_TP):
        //Note: Derivatives of g, s, cp, cv, w, eta, and lambda are not calculated in the two-phase region currently.
        DIFF_SAT_VU_SPL(str.ps, str.x, str.v1, str.v1s, str.dz_1.x1tmin, str.dz_1.x1tmax, str.dz_1.dvdu_vt, str.dz_1.K, str.v2, str.v2t, str.u1, str.u2, str.d_tp);
        dudv_t=dudv_p=str.d_tp.dudv_pt;
        p=str.ps;
        t=str.ts;
        x=str.x;
        alpha_l=(1.-x)*str.v1/v;
        DIFF_S_VU_L_CO2_T(str.v1s, str.dz_1.x1tmin, str.dz_1.x1tmax, str.dz_1.dvdu_vt, str.dz_1.K, str.u1, s1, dsdv_u_1, dsdu_v_1, dudv_s_1);
        DIFF_S_VU_G_CO2_T(str.v2t, str.v2, str.u2, s2, dsdv_u_2, dsdu_v_2, dudv_s_2);
        s=s1+x*(s2-s1);
        ds1dps=dsdu_v_1*str.d_tp.duldp+dsdv_u_1*str.d_tp.dvldp;
        ds2dps=dsdu_v_2*str.d_tp.duvdp+dsdv_u_2*str.d_tp.dvvdp;
        dsdps_u=ds1dps+str.d_tp.dxdp_u*(s2-s1)+x*(ds2dps-ds1dps);
        dsdu_ps=(s2-s1)/(str.u2-str.u1);
        dsdu_v=dsdu_ps+dsdps_u*dpdu_v;
        dsdv_ps=(s2-s1)/(str.v2-str.v1);
        dsdps_v=ds1dps+str.d_tp.dxdp_v*(s2-s1)+x*(ds2dps-ds1dps);
        dsdv_u=dsdv_ps-dsdps_v*dudv_p*dpdu_v;
        dudv_s=-dsdv_u/dsdu_v;
        //For cp, cv, w, eta, and lambda the properties itself are calculated, but not the derivatives.
        //cp is infinity in the two phase region (mass average is calculated)
        //the derivative (du/dt)_v is calculable in the two-phase region but its slope is discontinuous at dew and bubble curve (mass average is calculated)
#ifdef SBTL_USE_C_AUX
        cp1   =CP_VU_L_CO2_T(str.v1s, str.u1);
        cv1   =CV_VU_L_CO2_T(str.v1s, str.u1);
        cp2   =CP_VU_G_CO2_T(str.v2t,str.u1);
        cv2   =CV_VU_G_CO2_T(str.v2t,str.u2);
#else
        cp1   =CP_VU_L_CO2_T(str.v1s,str.dz_1.x1tmin,str.dz_1.x1tmax,str.dz_1.dvdu_vt,str.dz_1.K,str.u1);
        cv1   =CV_VU_L_CO2_T(str.v1s,str.dz_1.x1tmin,str.dz_1.x1tmax,str.dz_1.dvdu_vt,str.dz_1.K,str.u1);
        cp2   =CP_VU_G_CO2_T(str.v2t,str.v2,str.u1);
        cv2   =CV_VU_G_CO2_T(str.v2t,str.v2,str.u2);
#endif
        cp=cp1+x*(cp2-cp1);
        cv=cv1+x*(cv2-cv1);
        //Wood's equation for speed of sound
        alpha_v=1.-alpha_l;                         //alpha_v=x*v2/v;
        w1     =W_VU_L_CO2_T(str.v1s,str.u1);
        w2     =W_VU_G_CO2_T(str.v2t,str.u2);
        F1     =str.v1/(w1*w1);
        F2     =str.v2/(w2*w2);
        rho_m  =alpha_v/str.v2+(1.-alpha_v)/str.v1; //density of mixture
        F_m    =alpha_v*F2+(1.-alpha_v)*F1;         //comstressibility of mixture
        w      = 1./sqrt(F_m*rho_m);
        //viscosity and thermal conductivity depend on the vapor-volume fraction
        eta1   =ETA_VU_L_CO2_T(str.v1s,str.u1);
        eta2   =ETA_VU_G_CO2_T(str.v2t,str.u2);
        eta    =eta1+alpha_v*(eta2-eta1);
        lambda1=LAMBDA_VU_L_CO2_T(str.v1s,str.u1);
        lambda2=LAMBDA_VU_G_CO2_T(str.v2t,str.u2);
        lambda =lambda1+alpha_v*(lambda2-lambda1);
        break;
    default:
        return I_ERR;
    }
    return I_OK;
}
/*
SBTLAPI int __stdcall flash_pt_SBTLCO2(double p, double t, double& v, double& u) throw()
{
    int ireg;
    double vt=ERR_VAL;

    ireg=ireg_pt_SBTLCO2(p, t);

    switch(ireg) {
    case(IREG_L):
        return PT_FLASH_L(p, t, v, u);
    case(IREG_G):
        return PT_FLASH_G(p, t, v, vt, u);
    default:
        return I_ERR;
    }
    return I_OK;
}
//
SBTLAPI int __stdcall flash_deriv_pt_SBTLCO2(double p, double t, double& v, double& dvdp_t, double& dvdt_p, double& dpdt_v, double& u, double& dudp_t, double& dudt_p, double& dpdt_u) throw()
{
    int ireg;
    double vt=ERR_VAL;

    ireg=ireg_pt_SBTLCO2(p, t);

    switch(ireg) {
    case(IREG_L):
        return PT_FLASH_DERIV_L(p, t, v, dvdp_t, dvdt_p, dpdt_v, u, dudp_t, dudt_p, dpdt_u);
    case(IREG_G):
        return PT_FLASH_DERIV_G(p, t, v, vt, dvdp_t, dvdt_p, dpdt_v, u, dudp_t, dudt_p, dpdt_u);
    default:
        return I_ERR;
    }
    return I_OK;
}
//
SBTLAPI int __stdcall flash_px_SBTLCO2(double p, double x, double& v, double& u) throw()
{
    double ts,v1,v2,v2t,u1,u2;

    static const double p0=0.0005;
    static const double pc=7.37729837321;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;

    if(p<p0 || p>pc) {
        v=ERR_VAL;
        u=ERR_VAL;
        return I_ERR;
    }
    if((pc-p)<0.1) {
        v=vc;
        u=uc;
        return I_OK;
    } else {
        ts=TS_P_CO2(p);
        PT_FLASH_L(p, ts, v1, u1);
        PT_FLASH_G(p, ts, v2, v2t, u2);
        v=v1+x*(v2-v1);
        u=u1+x*(u2-u1);
        return I_OK;
    }
}
//
SBTLAPI int __stdcall flash_tx_SBTLCO2(double t, double x, double& v, double& u) throw()
{
    double ps,v1,v2,v2t,u1,u2;

    static const double t0=216.592;
    static const double tc=304.1282;
    static const double vc=1./467.60000128174;
    static const double uc=316.468709888;

    if(t<t0 || t>tc) {
        v=ERR_VAL;
        u=ERR_VAL;
        return I_ERR;
    }
    if((tc-t)<0.2) {
        v=vc;
        u=uc;
        return I_OK;
    } else {
        ps=PS_T_INV_CO2(t);
        PT_FLASH_L(ps, t, v1, u1);
        PT_FLASH_G(ps, t, v2, v2t, u2);
        v=v1+x*(v2-v1);
        u=u1+x*(u2-u1);
        return I_OK;
    }
}
//
SBTLAPI double __stdcall ts_p_SBTLCO2(double p) throw()
{
    static const double p0=0.0005;
    static const double pc=7.37729837321;

    if(p<p0 || p>pc) return ERR_VAL;
    return TS_P_CO2(p);
}
//
SBTLAPI double __stdcall ps_t_SBTLCO2(double t) throw()
{
    static const double t0=216.592;
    static const double tc=304.1282;

    if(t<t0 || t>tc) return ERR_VAL;
    return PS_T_INV_CO2(t);
}
//
SBTLAPI double __stdcall hvap_t_SBTLCO2(double t) throw()
{
    static const double t0=216.592;
    static const double tc=304.1282;

    double ps,vl,ul,vv,vvt,uv,hl,hv;

    if(t<t0 || t>tc) return ERR_VAL;
    if((tc-t)<0.2) return 0.;
    else {
        ps=PS_T_INV_CO2(t);
        if(PT_FLASH_L(ps,t,vl,ul)) return ERR_VAL;
        if(PT_FLASH_G(ps,t,vv,vvt,uv)) return ERR_VAL;
        hl=ul+ps*vl*1.e3;
        hv=uv+ps*vv*1.e3;
        return hv-hl;
    }
}
//
SBTLAPI int __stdcall flash_ps_SBTLCO2(double p, double s, double& v, double& u) throw()
{
    STR_vu_SBTL_CO2 str;
//
    v      =ERR_VAL;
    u      =ERR_VAL;
//
    //region determination: returns transformed value vt of vapor phase and saturation properties ps, ts, x, v1, v2, v2t, u1, u2 (if applicable)
    int str_state=str.GetStatePS(p,s);
    if(str_state<STR_PDP) ireg_ps_SBTLCO2(p, s, str);
//
    switch(str.ireg) {
    case(IREG_L):
        v     =str.v_;
        u     =str.u_;
        break;
    case(IREG_G):
        v     =str.v_;
        u     =str.u_;
        break;
    case(IREG_TP):
        v     =str.v_;
        u     =str.u_;
        break;
    default:
        return I_ERR;
    }
    return I_OK;
}
//
SBTLAPI int __stdcall flash_deriv_ps_SBTLCO2(double p, double s, double& v, double& dvdp_s, double& dvds_p, double& dpds_v, double& u, double& dudp_s, double& duds_p, double& dpds_u) throw()
{
    STR_vu_SBTL_CO2 str;
//
    double u;
    double vt;
    double dpdv_u,dpdu_v,dudv_p;
    double s1, dsdv_u_1, dsdu_v_1, dudv_s_1, s2, dsdv_u_2, dsdu_v_2, dudv_s_2;
    double ds1dps,ds2dps,dsdps_u,dsdu_ps,dsdu_v,dsdv_ps,dsdps_v,dsdv_u,dudv_s;

    int str_state=str.GetStatePS(p,s);
    if(str_state<STR_PDP) ireg_ps_SBTLCO2(p, s, str);

    switch(str.ireg) {
    case(IREG_L):
        v=str.v_;
        u=str.u_;
        PS_FLASH_DERIV_L(p, s, v, u, dvdp_s, dvds_p, dpds_v, dudp_s, duds_p, dpds_u);
        return I_OK;
    case(IREG_G):
        v=str.v_;
        vt=str.vt;
        u=str.u_;
        PS_FLASH_DERIV_G(p, s, v, vt, u, dvdp_s, dvds_p, dpds_v, dudp_s, duds_p, dpds_u);
        return I_OK;
    case(IREG_TP):
        v=str.v_;
        u=str.u_;
        if(str_state<STR_DTP) DIFF_SAT_VU_SPL(str.ps, str.x, str.v1, str.v1s, str.dz_1.x1tmin, str.dz_1.x1tmax, str.dz_1.dvdu_vt, str.dz_1.K, str.v2, str.v2t, str.u1, str.u2, str.d_tp);
        //p derivatives
        dpdv_u=str.d_tp.dpdv_u;
        dpdu_v=str.d_tp.dpdu_v;
        dudv_p=str.d_tp.dudv_pt;
        //s derivatives
        DIFF_S_VU_L_CO2_T(str.v1s, str.dz_1.x1tmin, str.dz_1.x1tmax, str.dz_1.dvdu_vt, str.dz_1.K, str.u1, s1, dsdv_u_1, dsdu_v_1, dudv_s_1);
        DIFF_S_VU_G_CO2_T(str.v2t, str.v2, str.u2, s2, dsdv_u_2, dsdu_v_2, dudv_s_2);
        ds1dps=dsdu_v_1*str.d_tp.duldp+dsdv_u_1*str.d_tp.dvldp;
        ds2dps=dsdu_v_2*str.d_tp.duvdp+dsdv_u_2*str.d_tp.dvvdp;
        dsdps_u=ds1dps+str.d_tp.dxdp_u*(s2-s1)+str.x*(ds2dps-ds1dps);
        dsdu_ps=(s2-s1)/(str.u2-str.u1);
        dsdu_v=dsdu_ps+dsdps_u*str.d_tp.dpdu_v;
        dsdv_ps=(s2-s1)/(str.v2-str.v1);
        dsdps_v=ds1dps+str.d_tp.dxdp_v*(s2-s1)+str.x*(ds2dps-ds1dps);
        dsdv_u=dsdv_ps-dsdps_v*str.d_tp.dudv_pt*str.d_tp.dpdu_v;
        dudv_s=-dsdv_u/dsdu_v;
        //p,v,s derivatives
        dpds_v=dpdu_v/dsdu_v;
        dvds_p=(str.v2-str.v1)/(s2-s1);
        dvdp_s=-dvds_p/dpds_v;
        return I_OK;
    default:
        v=ERR_VAL;
        u=ERR_VAL;
        dvdp_s=ERR_VAL;
        dvds_p=ERR_VAL;
        dpds_v=ERR_VAL;
        dudp_s=ERR_VAL;
        duds_p=ERR_VAL;
        dpds_u=ERR_VAL;
        return I_ERR;
    }
}
*/
