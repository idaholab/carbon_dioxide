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
// H_V_UC_CO2
//
///////////////////////////////////////////////////////////////////////////
//
#include "math.h"
#include "SBTL_call_conv.h"
#include "SBTL_def.h"
//
extern const double x1_PVUGCO2[];
extern const double x2_PVUGCO2[];
extern const double data_PVUGCO2[];
//
SBTLAPI double __stdcall H_V_UC_CO2(double x1_val) throw()
{
    unsigned int i=0;
    double x1f;
    double x1t;
    double dx1, p;
//
    static const double uc=316.468709888;
    static const double dx2=-0.50098708170003192;
//
    static const double x1_sub_RS_0=-6.8476110345088;
    static const double x1_sub_RS_1=-5.3061027618644;
    static const double ZS_1=-5.2694350176765;
    static const double dist_x1_inv_0=64.222814601034;
    static const double dist_x1_inv_1=17.311611400605;
    //static const double x2_sub_RS_0=309.41919191919;
    //static const double x2_sub_RS_1=424.41919191919;
    //static const double TS_1=429.92848484848;
    //static const double dist_x2_inv_0=0.86086956521739;
    //static const double dist_x2_inv_1=0.10145105755042;
//
// transformations
    x1t=log(x1_val);
//
    if(x1t>ZS_1) {
        x1f=(x1t-ZS_1)*dist_x1_inv_1;
        i=IROUND(x1f)+100;
        if(i>298) i=298;
    } else if(x1t<x1_sub_RS_1) {
        x1f=(x1t-x1_sub_RS_0)*dist_x1_inv_0;
        if(x1f>0.) i=IROUND(x1f);
        else i=0;
    } else {
        i=99;
    }
//
    //if(uc>TS_1) {
    //    x2f=(uc-TS_1)*dist_x2_inv_1;
    //    j=IROUND(x2f)+100;
    //    if(j>198) j=198;
    //} else if(uc<x2_sub_RS_1) {
    //    x2f=(uc-x2_sub_RS_0)*dist_x2_inv_0;
    //    if(x2f>0.) j=IROUND(x2f);
    //    else j=0;
    //} else {
    //    j=99;
    //}
//
    //const double *val=&data_PVUGCO2[9*(j*299+i)];
    const double *val=&data_PVUGCO2[9*(1794+i)];
//
    dx1=x1t-x1_PVUGCO2[i];
    //dx2=uc-x2_PVUGCO2[j];
//
    p= val[0]
        + dx2 * ( val[1] + dx2 *   val[2]  )
        + dx1 * ( val[3] + dx2 * ( val[4] + dx2 *  val[5]  ) + dx1 * ( val[6]
        + dx2 * ( val[7] + dx2 *   val[8]
          ) ) );
//
    return uc+p*x1_val*1.e3;
}
