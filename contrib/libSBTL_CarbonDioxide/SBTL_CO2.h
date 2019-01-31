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
// SBTL_CO2.h
//
#pragma once
//
//-----------------------------------------------------------------------------
// return values (error flags)
//-----------------------------------------------------------------------------
//
#define ERR_VAL -1.e3
#define I_OK  0
#define I_ERR 1
//
//-----------------------------------------------------------------------------
// region flags
//-----------------------------------------------------------------------------
//
#define IREG_ERR 0
#define IREG_L   1
#define IREG_G   2
#define IREG_TP  3
//
//-----------------------------------------------------------------------------
// struct states
//-----------------------------------------------------------------------------
//
                    //known properties:
#define STR_ERR 0   // - none (needs to be initialized with ireg-call)
#define STR_PDP 1   // - region
                    // - auxiliary variables for property calculations (single-phase region & two-phase region)
                    // - auxiliary variables for derivatives in the single-phase region
#define STR_DTP 2   // - all above
                    // - auxiliary variables for derivatives in the two-phase region (calculated if required, preserved for multiple use)
//
//-----------------------------------------------------------------------------
// structs to handle preliminary results and auxiliary variables:
//-----------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////////
// - auxiliary variables of first derivatives (liquid)
///////////////////////////////////////////////////////////////////////////////
//
typedef struct _DERIV_1 {
    double x1tmin;  //v(u)_min - 100 MPa isobar
    double x1tmax;  //v(u)_max - saturated liquid curve
    double dvdu_vt; //(dv/du)_vt
    double K;       //slope of scaling function: (x1zmax-x1zmin)/(x1tmax-x1tmin)
// constructor
    _DERIV_1() { reset();}
// reset
    void reset() {
        x1tmin      =ERR_VAL;
        x1tmax      =ERR_VAL;
        dvdu_vt     =ERR_VAL;
        K           =ERR_VAL;
    }
} DERIV_1;
//
///////////////////////////////////////////////////////////////////////////////
// - auxiliary variables of first and second derivatives (liquid)
///////////////////////////////////////////////////////////////////////////////
//
typedef struct _DERIV_2 {
    double x1tmin;      //v(u)_min
    double dvdu_min;    //(dv/du)_min
    double d2vdu2_min;  //(d2v/du2)_min
    double x1tmax;      //v(u)_max
    double dvdu_max;    //(dv/du)_max
    double d2vdu2_max;  //(d2v/du2)_max
    double dvtdu_v;     //(dvt/du)_v
    double dvdu_vt;     //(dv/du)_vt
    double K;           //(x1zmax-x1zmin)/(x1tmax-x1tmin)
    double d2vtdu2_v;   //(d2vt/du2)_v
    double d2vtdvdu;    //(d2vt/dvdu)
// constructor
    _DERIV_2() { reset();}
// reset
    void reset() {
        x1tmin      =ERR_VAL;
        dvdu_min    =ERR_VAL;
        d2vdu2_min  =ERR_VAL;
        x1tmax      =ERR_VAL;
        dvdu_max    =ERR_VAL;
        d2vdu2_max  =ERR_VAL;
        dvtdu_v     =ERR_VAL;
        dvdu_vt     =ERR_VAL;
        K           =ERR_VAL;
        d2vtdu2_v   =ERR_VAL;
        d2vtdvdu    =ERR_VAL;
    }
} DERIV_2;
//
///////////////////////////////////////////////////////////////////////////////
// - derivatives in the two-phase region
///////////////////////////////////////////////////////////////////////////////
//
typedef struct _DERIV_TP {
    double ps_;     //values below computed for this
    double x_;      //state point

    double dtsdp;   //dts/dp
    double dvldp;   //dv'/dp
    double dvvdp;   //dv"/dp
    double duldp;   //du'/dp
    double duvdp;   //du"/dp
    double dxdp_v;  //(dx/dp)_v
    double dxdp_u;  //(dx/dp)_u
    double dpdv_u;  //(dp/dv)_u
    double dpdu_v;  //(dp/du)_v
    double dtdv_u;  //(dt/dv)_u
    double dtdu_v;  //(dt/du)_v
    double dudv_pt; //(du/dv)_p = (du/dv)_t
// constructor
    _DERIV_TP() { reset();}
// reset
    void reset() {
        ps_    =ERR_VAL;
        x_     =ERR_VAL;

        dtsdp  =ERR_VAL;
        dvldp  =ERR_VAL;
        dvvdp  =ERR_VAL;
        duldp  =ERR_VAL;
        duvdp  =ERR_VAL;
        dxdp_v =ERR_VAL;
        dxdp_u =ERR_VAL;
        dpdv_u =ERR_VAL;
        dpdu_v =ERR_VAL;
        dtdv_u =ERR_VAL;
        dtdu_v =ERR_VAL;
        dudv_pt=ERR_VAL;
    }
// 
    bool IsCalculated(double ps, double x) {
        //if(ps_>0. && x_>0. && ps_==ps && x_==x) {
        if(ps_==ps && x_==x) {
            return true;
        } else {
            reset();
            return false;
        }
    }

} DERIV_TP;
//
///////////////////////////////////////////////////////////////////////////////
// - struct to be used by ireg_vu_SBTL_CO2 (DERIV_TP d_tp is used by DIFF_SAT_VU_SPL_CO2 directly)
// - this is also used by ireg_pv_SBTL_CO2 and ireg_ps_SBTL_CO2
///////////////////////////////////////////////////////////////////////////////
//
typedef struct _STR_vu_SBTL_CO2 {
//
    double v_;      //values below computed for this
    double u_;      //state point
    double p_;      //check p_,v_ for given (p,v), p_,s_ for given (p,s), and p_,h_ for given (p,h)
    double h_;      //check p_,h_ for given (p,h) and h_,s_ for given (h,s)
    double s_;      //check p_,s_ for given (p,s) and h_,s_ for given (h,s)
//
    int    ireg;    //region flag (see above)
    double vls;     //scaled volume             (liquid phase only)
#ifndef SBTL_USE_C_AUX
    DERIV_2 dz_l;   //1st and 2nd derivatives of p and t are required to calculate cp, cv, and their derivatives (liquid phase only)
#else
    DERIV_1 dz_l;   //1st derivatives of z (liquid phase only)
#endif
    double vt;      //transformed volume        (gas phase only)
    double ps;      //saturation pressure       (two-phase region only)
    double ts;      //saturation temperature    (two-phase region only)
    double x;       //vapor fraction            (two-phase region only)
    double v1;      //sat. liquid volume        (two-phase region only)
    double v1s;     //scaled sat. liquid volume (two-phase region only)
#ifndef SBTL_USE_C_AUX
    DERIV_2 dz_1;   //1st and 2nd derivatives of p and t are required to calculate cp1, cv1, and their derivatives (two-phase region only)
#else
    DERIV_1 dz_1;   //sat. liquid: derivatives of z (two-phase region only)
#endif
    double v2;      //sat. vapor volume         (two-phase region only)
    double v2t;     //transf. sat. vapor volume (two-phase region only)
    double u1;      //sat. liquid int. energy   (two-phase region only)
    double u2;      //sat. vapor int. energy    (two-phase region only)
    DERIV_TP d_tp;  //derivatives in two-phase region
// constructor
    _STR_vu_SBTL_CO2() { reset();}
// reset
    void reset() {
        v_      =ERR_VAL;
        u_      =ERR_VAL;
        p_      =ERR_VAL;
        h_      =ERR_VAL;
        s_      =ERR_VAL;
        //
        ireg    =IREG_ERR;
        vls     =ERR_VAL;
        vt      =ERR_VAL;
        ps      =ERR_VAL;
        ts      =ERR_VAL;
        x       =ERR_VAL;
        v1      =ERR_VAL;
        v1s     =ERR_VAL;
        v2      =ERR_VAL;
        v2t     =ERR_VAL;
        u1      =ERR_VAL;
        u2      =ERR_VAL;
    }
//
    int GetState(double v, double u) {
        if(v_==v && u_==u) {
            if(ireg==IREG_TP) {
                if(d_tp.IsCalculated(ps,x)) return STR_DTP;
                else                        return STR_PDP;
            } else {
                return STR_PDP;
            }
        } else {
            reset();
            return STR_ERR;
        }
    }
//
    int GetStatePV(double p, double v) {
        if(p_==p && v_==v) {
            if(ireg==IREG_TP) {
                if(d_tp.IsCalculated(ps,x)) return STR_DTP;
                else                        return STR_PDP;
            } else {
                return STR_PDP;
            }
        } else {
            reset();
            return STR_ERR;
        }
    }
//
    int GetStatePS(double p, double s) {
        if(p_==p && s_==s) {
            if(ireg==IREG_TP) {
                if(d_tp.IsCalculated(ps,x)) return STR_DTP;
                else                        return STR_PDP;
            } else {
                return STR_PDP;
            }
        } else {
            reset();
            return STR_ERR;
        }
    }
//
    int GetStatePH(double p, double h) {
        if(p_==p && h_==h) {
            if(ireg==IREG_TP) {
                if(d_tp.IsCalculated(ps,x)) return STR_DTP;
                else                        return STR_PDP;
            } else {
                return STR_PDP;
            }
        } else {
            reset();
            return STR_ERR;
        }
    }
//
    int GetStateHS(double h, double s) {
        if(h_==h && s_==s) {
            if(ireg==IREG_TP) {
                if(d_tp.IsCalculated(ps,x)) return STR_DTP;
                else                        return STR_PDP;
            } else {
                return STR_PDP;
            }
        } else {
            reset();
            return STR_ERR;
        }
    }
//
    int GetStateVH(double v, double h) {
        if(v_==v && h_==h) {
            if(ireg==IREG_TP) {
                if(d_tp.IsCalculated(ps,x)) return STR_DTP;
                else                        return STR_PDP;
            } else {
                return STR_PDP;
            }
        } else {
            reset();
            return STR_ERR;
        }
    }
} STR_vu_SBTL_CO2;
////
/////////////////////////////////////////////////////////////////////////////////
//// - struct to be used by ireg_vu_SBTL_CO2M
//// - this is also used by ireg_pv_SBTL_CO2M and ireg_ps_SBTL_CO2M
/////////////////////////////////////////////////////////////////////////////////
////
//typedef struct _STR_vu_SBTL_CO2M {
////
//    double v_;      //values below computed for this
//    double u_;      //state point
//    double p_;      //check p_,v_ for given (p,v), p_,s_ for given (p,s), and p_,h_ for given (p,h)
//    double h_;      //check p_,h_ for given (p,h) and h_,s_ for given (h,s)
//    double s_;      //check p_,s_ for given (p,s) and h_,s_ for given (h,s)
////
//    int    ireg;    //region flag (see above)
//    double vls;     //scaled volume             (liquid phase only)
//#ifndef SBTL_USE_C_AUX
//    DERIV_2 dz_l;   //1st and 2nd derivatives of p and t are required to calculate cp, cv, and their derivatives (liquid phase only)
//#else
//    DERIV_1 dz_l;   //derivatives of z (liquid phase only)
//#endif
//    double vt;      //transformed volume        (gas phase only)
//// constructor
//    _STR_vu_SBTL_CO2M() { reset();}
//// reset
//    void reset() {
//        v_      =ERR_VAL;
//        u_      =ERR_VAL;
//        p_      =ERR_VAL;
//        h_      =ERR_VAL;
//        s_      =ERR_VAL;
//        //
//        ireg    =IREG_ERR;
//        vls     =ERR_VAL;
//        vt      =ERR_VAL;
//    }
////
//    int GetState(double v, double u) {
//        if(v_==v && u_==u) {
//            return STR_PDP;
//        } else {
//            return STR_ERR;
//        }
//    }
////
//    int GetStatePV(double p, double v) {
//        if(p_==p && v_==v) {
//            return STR_PDP;
//        } else {
//            reset();
//            return STR_ERR;
//        }
//    }
////
//    int GetStatePS(double p, double s) {
//        if(p_==p && s_==s) {
//            return STR_PDP;
//        } else {
//            reset();
//            return STR_ERR;
//        }
//    }
////
//    int GetStatePH(double p, double h) {
//        if(p_==p && h_==h) {
//            return STR_PDP;
//        } else {
//            reset();
//            return STR_ERR;
//        }
//    }
////
//    int GetStateHS(double h, double s) {
//        if(h_==h && s_==s) {
//            return STR_PDP;
//        } else {
//            reset();
//            return STR_ERR;
//        }
//    }
////
//    int GetStateVH(double v, double h) {
//        if(v_==v && h_==h) {
//            return STR_PDP;
//        } else {
//            reset();
//            return STR_ERR;
//        }
//    }
//} STR_vu_SBTL_CO2M;
////
