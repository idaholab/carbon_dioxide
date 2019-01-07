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
// SBTL_func
//
#pragma once
//
#include "SBTL_CO2.h"
#include "SBTL_call_conv.h"
//
//-----------------------------------------------------------------------------
//SBTL functions
//-----------------------------------------------------------------------------
//
// saturation curve
extern "C" double __stdcall TS_P_CO2(double p);
extern "C" double __stdcall PS_T_INV_CO2(double t);
extern "C" void   __stdcall DIFF_TS_P_CO2(double p, double& ts, double& dtsdp);
//
// (p,t)-flash routines
extern "C" int    __stdcall PT_FLASH_L(double p, double t, double& v, double& u);
extern "C" int    __stdcall PT_FLASH_DERIV_L(double p, double t, double& v, double& dvdp_t, double& dvdt_p, double& dpdt_v, double& u, double& dudp_t, double& dudt_p, double& dpdt_u);
extern "C" int    __stdcall PT_FLASH_G(double p, double t, double& v, double& vt, double& u);
extern "C" int    __stdcall PT_FLASH_DERIV_G(double p, double t, double& v, double& vt, double& dvdp_t, double& dvdt_p, double& dpdt_v, double& u, double& dudp_t, double& dudt_p, double& dpdt_u);
//
// (p,v)-inverse spline functions
extern "C" double __stdcall U_VP_L_CO2(double v, double p);
extern "C" void   __stdcall DIFF_U_VP_L_CO2(double v, double p, double &u, double& dudv_p, double& dudp_v, double& dpdv_u);
extern "C" double __stdcall U_VP_G_CO2(double v, double p);
extern "C" void   __stdcall DIFF_U_VP_G_CO2(double v, double p, double &u, double& dudv_p, double& dudp_v, double& dpdv_u);
//
// (p,s)-flash routines
extern "C" int    __stdcall PS_FLASH_L(double p, double s, double& v, double& u);
extern "C" void   __stdcall PS_FLASH_DERIV_L(double v, double u, double& dvdp_s, double& dvds_p, double& dpds_v, double& dudp_s, double& duds_p, double& dpds_u);
extern "C" int    __stdcall PS_FLASH_G(double p, double s, double& v, double& vt, double& u);
extern "C" void   __stdcall PS_FLASH_DERIV_G(double v, double vt, double u, double& dvdp_s, double& dvds_p, double& dpds_v, double& dudp_s, double& duds_p, double& dpds_u);
//
// (p,h)-flash routines
extern "C" int    __stdcall PH_FLASH_L(double p, double h, double& v, double& u);
extern "C" void   __stdcall PH_FLASH_DERIV_L(double p, double v, double u, double& dvdp_h, double& dvdh_p, double& dpdh_v, double& dudp_h, double& dudh_p, double& dpdh_u);
extern "C" void   __stdcall PH_T_FLASH_DERIV_L(double p, double v, double u, double& t, double& dtdp_h, double& dtdh_p, double& dpdh_t);
extern "C" int    __stdcall PH_FLASH_G(double p, double h, double& v, double& vt, double& u);
extern "C" void   __stdcall PH_FLASH_DERIV_G(double p, double v, double vt, double u, double& dvdp_h, double& dvdh_p, double& dpdh_v, double& dudp_h, double& dudh_p, double& dpdh_u);
extern "C" void   __stdcall PH_T_FLASH_DERIV_G(double p, double v, double vt, double u, double& t, double& dtdp_h, double& dtdh_p, double& dpdh_t);
//
// (h,s)-flash routines
extern "C" int    __stdcall HS_FLASH_L(double h, double s, double& v, double& u);
extern "C" void   __stdcall HS_FLASH_DERIV_L(double v, double u, double& dvdh_s, double& dvds_h, double& dhds_v, double& dudh_s, double& duds_h, double& dhds_u);
extern "C" void   __stdcall HS_PT_FLASH_DERIV_L(double v, double u, double& p, double& dpdh_s, double& dpds_h, double& dhds_p, double& t, double& dtdh_s, double& dtds_h, double& dhds_t);
extern "C" int    __stdcall HS_FLASH_G(double h, double s, double& v, double& vt, double& u);
extern "C" void   __stdcall HS_FLASH_DERIV_G(double v, double vt, double u, double& dvdh_s, double& dvds_h, double& dhds_v, double& dudh_s, double& duds_h, double& dhds_u);
extern "C" void   __stdcall HS_PT_FLASH_DERIV_G(double v, double vt, double u, double& p, double& dpdh_s, double& dpds_h, double& dhds_p, double& t, double& dtdh_s, double& dtds_h, double& dhds_t);
extern "C" int    __stdcall SAT_HS_SPL(double h, double s, double& ps, double& ts, double& xvap, double& vl, double& vv, double& vvt, double& ul, double& uv);
//
// (v,h)-flash routines
extern "C" int    __stdcall SAT_VH_SPL(double v, double h, double& ps, double& ts, double& xvap, double& vl, double& vv, double& vvt, double& ul, double& uv);
extern "C" int    __stdcall VH_FLASH_L(double v, double h, double& u);
extern "C" int    __stdcall VH_FLASH_G_T(double v, double vt, double h, double& u);
//
// (v,u)-functions with transformed inputs
extern "C" double __stdcall P_VU_L_CO2_T(double vs, double u);
extern "C" double __stdcall P_VU_G_CO2_T(double vt, double u);
extern "C" double __stdcall T_VU_L_CO2_T(double vs, double u);
extern "C" double __stdcall T_VU_G_CO2_T(double vt, double u);
extern "C" double __stdcall S_VU_L_CO2_T(double vs, double u);
extern "C" double __stdcall S_VU_G_CO2_T(double vt, double u);
#ifdef SBTL_USE_C_AUX
extern "C" double __stdcall CP_VU_L_CO2_T(double vs, double u);
extern "C" double __stdcall CP_VU_G_CO2_T(double vt, double u);
extern "C" double __stdcall CV_VU_L_CO2_T(double vs, double u);
extern "C" double __stdcall CV_VU_G_CO2_T(double vt, double u);
#else
extern "C" double __stdcall CP_VU_L_CO2_T(double vs, double x1tmin, double x1tmax, double dvdu_vt, double K, double u);
extern "C" double __stdcall CP_VU_G_CO2_T(double vt, double v, double u);
extern "C" double __stdcall CV_VU_L_CO2_T(double vs, double x1tmin, double x1tmax, double dvdu_vt, double K, double u);
extern "C" double __stdcall CV_VU_G_CO2_T(double vt, double v, double u);
#endif
extern "C" double __stdcall W_VU_L_CO2_T(double vs, double u);
extern "C" double __stdcall W_VU_G_CO2_T(double vt, double u);
extern "C" double __stdcall ETA_VU_L_CO2_T(double vs, double u);
extern "C" double __stdcall ETA_VU_G_CO2_T(double vt, double u);
extern "C" double __stdcall LAMBDA_VU_L_CO2_T(double vs, double u);
extern "C" double __stdcall LAMBDA_VU_G_CO2_T(double vt, double u);
//
// scaling curves for the liquid region
extern "C" double __stdcall V_U_PMAX_CO2(double u);
extern "C" void   __stdcall DIFF_V_U_PMAX_CO2(double u, double& v, double& dvdu);
extern "C" void   __stdcall DIFF2_V_U_PMAX_CO2(double u, double& v, double& dvdu, double& d2vdu2);
extern "C" double __stdcall V1_U_SPL_CO2(double u);
extern "C" void   __stdcall DIFF_V1_U_SPL_CO2(double u, double& v, double& dvdu);
extern "C" void   __stdcall DIFF2_V1_U_SPL_CO2(double u, double& v, double& dvdu, double& d2vdu2);
//
// (v,u)-functions with transformed inputs and derivatives
extern "C" void __stdcall DIFF_P_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_T_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& t, double& dtdv, double& dtdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& s, double& dsdv, double& dsdu, double& dudv);
#ifdef SBTL_USE_C_AUX
extern "C" void __stdcall DIFF_CP_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& cp, double& dcpdv, double& dcpdu, double& dudv);
extern "C" void __stdcall DIFF_CV_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& cv, double& dcvdv, double& dcvdu, double& dudv);
#else
extern "C" void __stdcall DIFF_CP_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvtdu_v, double dvdu_vt, double K, double d2vtdu2_v, double d2vtdvdu, double u, double& cp, double& dcpdv, double& dcpdu, double& dudv);
extern "C" void __stdcall DIFF_CV_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvtdu_v, double dvdu_vt, double K, double d2vtdu2_v, double d2vtdvdu, double u, double& cv, double& dcvdv, double& dcvdu, double& dudv);
#endif
extern "C" void __stdcall DIFF_W_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& w, double& dwdv, double& dwdu, double& dudv);
extern "C" void __stdcall DIFF_ETA_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& eta, double& detadv, double& detadu, double& dudv);
extern "C" void __stdcall DIFF_LAMBDA_VU_L_CO2_T(double vt, double x1tmin, double x1tmax, double dvdu_vt, double K, double u, double& lambda, double& dlambdadv, double& dlambdadu, double& dudv);
extern "C" void __stdcall DIFF_P_VU_G_CO2_T(double vt, double v, double u, double& p, double& dpdv, double& dpdu, double& dudv);
extern "C" void __stdcall DIFF_T_VU_G_CO2_T(double vt, double v, double u, double& t, double& dtdv, double& dtdu, double& dudv);
extern "C" void __stdcall DIFF_S_VU_G_CO2_T(double vt, double v, double u, double& s, double& dsdv, double& dsdu, double& dudv);
extern "C" void __stdcall DIFF_CP_VU_G_CO2_T(double vt, double v, double u, double& cp, double& dcpdv, double& dcpdu, double& dudv);
extern "C" void __stdcall DIFF_CV_VU_G_CO2_T(double vt, double v, double u, double& cv, double& dcvdv, double& dcvdu, double& dudv);
extern "C" void __stdcall DIFF_W_VU_G_CO2_T(double vt, double v, double u, double& w, double& dwdv, double& dwdu, double& dudv);
extern "C" void __stdcall DIFF_ETA_VU_G_CO2_T(double vt, double v, double u, double& eta, double& detadv, double& detadu, double& dudv);
extern "C" void __stdcall DIFF_LAMBDA_VU_G_CO2_T(double vt, double v, double u, double& lambda, double& dlambdadv, double& dlambdadu, double& dudv);
extern "C" void __stdcall DIFF_SAT_VU_SPL(double ps, double x, double vl, double vlt, double x1tmin, double x1tmax, double dvdu_vt, double K, double vv, double vvt, double ul, double uv, DERIV_TP& d_tp);
//
// functions to be used in ireg_xx_SBTLCO2, ... (region boundaries)
extern "C" double __stdcall U_V_TMAX_AUX_CO2_T(double vt);
extern "C" double __stdcall U2_V_AUX_CO2_T(double vt);
extern "C" double __stdcall V_U_PMAX_AUX_CO2(double u);
extern "C" double __stdcall V1_U_SPL_CO2(double u);
extern "C" int    __stdcall SAT_U1_SPL(double u, double& ps, double& ts, double& vl);
extern "C" int    __stdcall SAT_V2_SPL_T(double vt, double& ps, double& ts, double& uv);
extern "C" int    __stdcall SAT_VU_SPL(double v, double u, double& ps, double& ts, double& x, double& vl, double& vv, double& vvt, double& ul, double& uv);
extern "C" double __stdcall T_P_UC_AUX_CO2(double p);
extern "C" double __stdcall V_P_UC_CO2_T(double p);
extern "C" double __stdcall S_H_PMAX_AUX_CO2(double h);
extern "C" double __stdcall S_H_PMIN_AUX_CO2(double h);
extern "C" double __stdcall H_S_TMAX_AUX_CO2(double s);
extern "C" double __stdcall H2_S_AUX_CO2(double s);
extern "C" double __stdcall S1_H_AUX_CO2(double h);
extern "C" int    __stdcall SAT_H1_SPL(double h, double& ps, double& ts, double& vl, double& ul, double& sl);
extern "C" int    __stdcall SAT_S2_SPL(double s, double& ps, double& ts, double& vv, double& vvt, double& uv);
extern "C" double __stdcall H_S_UC_CO2(double s);
extern "C" double __stdcall H_V_UC_CO2(double v);
extern "C" double __stdcall V1_H_AUX_CO2(double h);
extern "C" double __stdcall V_H_PMAX_AUX_CO2(double h);
extern "C" double __stdcall H2_V_AUX_CO2_T(double vt);
extern "C" double __stdcall V_H_PMIN_AUX_CO2(double h);
extern "C" double __stdcall H_V_TMAX_AUX_CO2_T(double vt);
//
