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
#pragma once
//
#include "SBTL_CO2.h"
#include "SBTL_call_conv.h"
//
//-----------------------------------------------------------------------------
// forward declarations of ireg-functions:
//-----------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////////
// ireg_deriv_vu_SBTLCO2:
//      This function is to be called before a sequence of calls to:
//          - p_vu_SBTLCO2(double v, double u, STR_vu_SBTLCO2 r);
//          - t_vu_SBTLCO2(double v, double u, STR_vu_SBTLCO2 r);
//          - ...
//          - p_deriv_vu_SBTLCO2(double v, double u, STR_vu_SBTLCO2 r, double& p, double& dpdv_u, double& dpdu_v, double& dudv_p);
//          - t_deriv_vu_SBTLCO2(double v, double u, STR_vu_SBTLCO2 r, double& t, double& dtdv_u, double& dtdu_v, double& dudv_t);
//          - ...
///////////////////////////////////////////////////////////////////////////////
//
SBTLAPI void __stdcall ireg_vu_SBTLCO2(double v, double u, STR_vu_SBTL_CO2& r);
//
///////////////////////////////////////////////////////////////////////////////
// ireg_pt_SBTLCO2
///////////////////////////////////////////////////////////////////////////////
//
SBTLAPI int __stdcall ireg_pt_SBTLCO2(double p, double t);
//
///////////////////////////////////////////////////////////////////////////////
// ireg_pv_SBTLCO2
///////////////////////////////////////////////////////////////////////////////
//
SBTLAPI void __stdcall ireg_pv_SBTLCO2(double p, double v, STR_vu_SBTL_CO2& r);
//
///////////////////////////////////////////////////////////////////////////////
// ireg_ps_SBTLCO2
///////////////////////////////////////////////////////////////////////////////
//
SBTLAPI void __stdcall ireg_ps_SBTLCO2(double p, double s, STR_vu_SBTL_CO2& r);
//
///////////////////////////////////////////////////////////////////////////////
// ireg_ph_SBTLCO2
///////////////////////////////////////////////////////////////////////////////
//
SBTLAPI void __stdcall ireg_ph_SBTLCO2(double p, double h, STR_vu_SBTL_CO2& r);
//
///////////////////////////////////////////////////////////////////////////////
// ireg_hs_SBTLCO2
///////////////////////////////////////////////////////////////////////////////
//
SBTLAPI void __stdcall ireg_hs_SBTLCO2(double h, double s, STR_vu_SBTL_CO2& r);
//
///////////////////////////////////////////////////////////////////////////////
// ireg_vh_SBTLCO2
///////////////////////////////////////////////////////////////////////////////
//
SBTLAPI void __stdcall ireg_vh_SBTLCO2(double v,double h, STR_vu_SBTL_CO2& r);
//