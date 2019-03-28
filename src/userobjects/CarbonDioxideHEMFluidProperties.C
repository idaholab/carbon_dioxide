#include "CarbonDioxideHEMFluidProperties.h"
#include "contrib/libSBTL_CarbonDioxide/LibSBTL_vu_CO2.h"
#include "contrib/libSBTL_CarbonDioxide/SBTL_CO2.h"
#include "contrib/libSBTL_CarbonDioxide/SBTL_def.h"
#include "contrib/libSBTL_CarbonDioxide/SBTL_func.h"

registerMooseObject("CarbonDioxideApp", CarbonDioxideHEMFluidProperties);

const Real CarbonDioxideHEMFluidProperties::_P_critical = 7.37729837321E+6;

extern "C" double __stdcall SIGMA_TS_CO2(double t);
extern "C" void __stdcall DIFF_TS_P_CO2(double p, double & ts, double & dtsdp);

template <>
InputParameters
validParams<CarbonDioxideHEMFluidProperties>()
{
  MooseEnum ss_flag("WOOD FROZEN THERMAL_EQ", "WOOD");
  InputParameters params = validParams<HEMFluidProperties>();
  params += validParams<NaNInterface>();
  params.addClassDescription(
      "Fluid properties of carbon dioxide for the homogeneous equilibrium model.");
  params.addParam<MooseEnum>(
      "ss_flag", ss_flag, "flag for type of material for sound speed method");
  return params;
}

CarbonDioxideHEMFluidProperties::CarbonDioxideHEMFluidProperties(const InputParameters & parameters)
  : HEMFluidProperties(parameters),
    NaNInterface(this),
    _to_MPa(1e-6),
    _to_Pa(1e6),
    _to_kJ(1e-3),
    _to_J(1e3)
{
  MooseEnum ss_flag(getParam<MooseEnum>("ss_flag"));
  if (ss_flag == "WOOD")
    _ss_flag = SS_WOOD;
  else if (ss_flag == "FROZEN")
    _ss_flag = SS_FROZEN;
  else if (ss_flag == "THERMAL_EQ")
    _ss_flag = SS_THERMAL_EQ;
}

Real
CarbonDioxideHEMFluidProperties::pressure(Real v, Real e) const
{
  e *= _to_kJ;
  if (_td_props.GetState(v, e) < STR_PDP)
    ireg_vu_SBTLCO2(v, e, _td_props);
  switch (_td_props.ireg)
  {
    case (IREG_L):
      return P_VU_L_CO2_T(_td_props.vls, e) * _to_Pa;
    case (IREG_G):
      return P_VU_G_CO2_T(_td_props.vt, e) * _to_Pa;
    case (IREG_TP):
      return _td_props.ps * _to_Pa;
    default:
      return getNaN();
  }
}

Real
CarbonDioxideHEMFluidProperties::temperature(Real v, Real e) const
{
  if (_td_props.GetState(v, e * _to_kJ) < STR_PDP)
    ireg_vu_SBTLCO2(v, e * _to_kJ, _td_props);
  switch (_td_props.ireg)
  {
    case (IREG_L):
      return T_VU_L_CO2_T(_td_props.vls, e * _to_kJ);
    case (IREG_G):
      return T_VU_G_CO2_T(_td_props.vt, e * _to_kJ);
    case (IREG_TP):
      return _td_props.ts;
    default:
      return getNaN();
  }
}

Real
CarbonDioxideHEMFluidProperties::h(Real pressure, Real temperature, Real quality) const
{
  int ierr;
  Real v, vt, u;
  Real p = pressure * _to_MPa;
  Real t = temperature;

  if (quality <= 0.0)
    ierr = PT_FLASH_L_CO2(p, t, v, u);

  else if (quality >= 1.0)
    ierr = PT_FLASH_G_CO2(p, t, v, vt, u);

  else
  {
    Real hf, hg, hfg;
    h_lat(temperature, hf, hg, hfg);
    return hf + quality * (hg - hf);
  }

  if (ierr == I_OK)
    return (u * _to_J + pressure * v);
  else
    return getNaN();
}

Real
CarbonDioxideHEMFluidProperties::quality(Real v, Real e) const
{
  int state = _td_props.GetState(v, e * _to_kJ);
  if (state < STR_PDP)
    ireg_vu_SBTLCO2(v, e * _to_kJ, _td_props);

  return _td_props.x;
}

Real
CarbonDioxideHEMFluidProperties::quality_Tsat_h(Real Tsat, Real h) const
{
  Real hf, hg, hfg;
  h_lat(Tsat, hf, hg, hfg);
  return (h - hf) / (hg - hf);
}

Real
CarbonDioxideHEMFluidProperties::alpha_vapor(Real v, Real e) const
{
  if (_td_props.GetState(v, e * _to_kJ) < STR_PDP)
    ireg_vu_SBTLCO2(v, e * _to_kJ, _td_props);
  switch (_td_props.ireg)
  {
    case (IREG_L):
      return 0;
    case (IREG_G):
      return 1;
    case (IREG_TP):
      // x - quality, m - mass, v - specific volume, vol - volume
      // _t - total, _l - liquid, _v - vapor
      //   1 - (1 - x)*v_l/v_t
      // = 1 - (1 - m_v/m_t)*(vol_l/m_l)/(vol_t/m_t)
      // = 1 - m_l/m_t * m_t/m_l * vol_l/vol_t
      // = 1 - vol_l/vol_t
      // = vol_v/vol_t
      return 1.0 - (1.0 - _td_props.x) * _td_props.v1 / v;
    default:
      return getNaN();
  }
}

Real
CarbonDioxideHEMFluidProperties::c(Real v, Real e) const
{
  Real alpha_l, alpha_v, w1, w2, K1, K2, rho_m, K_m, L, cp1, ts, dt_dp;

  if (_td_props.GetState(v, e * _to_kJ) < STR_PDP)
    ireg_vu_SBTLCO2(v, e * _to_kJ, _td_props);

  switch (_td_props.ireg)
  {
    case (IREG_L):
      return W_VU_L_CO2_T(_td_props.vls, e * _to_kJ);
    case (IREG_G):
      return W_VU_G_CO2_T(_td_props.vt, e * _to_kJ);
    case (IREG_TP):
      alpha_v = _td_props.x * _td_props.v2 / v;
      alpha_l = 1. - alpha_v;
      w1 = W_VU_L_CO2_T(_td_props.v1s, _td_props.u1);
      w2 = W_VU_G_CO2_T(_td_props.v2t, _td_props.u2);
      rho_m = alpha_l / _td_props.v1 + alpha_v / _td_props.v2; // density of mixture
      switch (_ss_flag)
      {
        case (SS_WOOD):
          K1 = _td_props.v1 / (w1 * w1);
          K2 = _td_props.v2 / (w2 * w2);
          K_m = alpha_l * K1 + alpha_v * K2;
          return 1. / sqrt(K_m * rho_m);
        case (SS_FROZEN):
          K1 = alpha_l / (_td_props.v1 * rho_m);
          K2 = alpha_v / (_td_props.v2 * rho_m);
          return sqrt(K1 * w1 * w1 + K2 * w2 * w2);
        case (SS_THERMAL_EQ):
          K1 = _td_props.v1 / (w1 * w1);
          K2 = _td_props.v2 / (w2 * w2);
          L = (_td_props.u2 - _td_props.u1) + _td_props.ps * (_td_props.v2 - _td_props.v1);
#ifdef SBTL_USE_C_AUX
          cp1 = CP_VU_L_CO2_T(_td_props.v1s, _td_props.u1);
#else
          cp1 = CP_VU_L_CO2_T(_td_props.v1s,
                              _td_props.dz_1.x1tmin,
                              _td_props.dz_1.x1tmax,
                              _td_props.dz_1.dvdu_vt,
                              _td_props.dz_1.K,
                              _td_props.u1);
#endif
          DIFF_TS_P_CO2(_td_props.ps, ts, dt_dp);
          K_m = alpha_l * K1 + alpha_v * K2 +
                alpha_l * _td_props.v2 * cp1 * dt_dp / (_td_props.v1 * L * 1.e6);
          return 1. / sqrt(K_m * rho_m);
        default:
          mooseError("Unknown 2phase flag");
      }
    default:
      return getNaN();
  }
}

Real
CarbonDioxideHEMFluidProperties::cp(Real v, Real e) const
{
  Real cp1, cp2;

  if (_td_props.GetState(v, e * _to_kJ) < STR_PDP)
    ireg_vu_SBTLCO2(v, e * _to_kJ, _td_props);

  switch (_td_props.ireg)
  {
    case (IREG_L):
#ifdef SBTL_USE_C_AUX
      return CP_VU_L_CO2_T(_td_props.vls, e * _to_kJ) * _to_J;
#else
      return CP_VU_L_CO2_T(_td_props.vls,
                           _td_props.dz_l.x1tmin,
                           _td_props.dz_l.x1tmax,
                           _td_props.dz_l.dvdu_vt,
                           _td_props.dz_l.K,
                           e * _to_kJ) *
             _to_J;
#endif
    case (IREG_G):
#ifdef SBTL_USE_C_AUX
      return CP_VU_G_CO2_T(_td_props.vt, e * _to_kJ) * _to_J;
#else
      return CP_VU_G_CO2_T(_td_props.vt, v, e * _to_kJ) * _to_J;
#endif
    case (IREG_TP):
// cp is infinity in the two phase region (mass average is calculated)
#ifdef SBTL_USE_C_AUX
      cp1 = CP_VU_L_CO2_T(_td_props.v1s, _td_props.u1);
      cp2 = CP_VU_G_CO2_T(_td_props.v2t, _td_props.u2);
#else
      cp1 = CP_VU_L_CO2_T(_td_props.v1s,
                          _td_props.dz_1.x1tmin,
                          _td_props.dz_1.x1tmax,
                          _td_props.dz_1.dvdu_vt,
                          _td_props.dz_1.K,
                          _td_props.u1);
      cp2 = CP_VU_G_CO2_T(_td_props.v2t, _td_props.v2, _td_props.u2);
#endif
      return (cp1 + _td_props.x * (cp2 - cp1)) * _to_J;
    default:
      return getNaN();
  }
}

Real
CarbonDioxideHEMFluidProperties::cv(Real v, Real e) const
{
  Real cv1, cv2;

  if (_td_props.GetState(v, e * _to_kJ) < STR_PDP)
    ireg_vu_SBTLCO2(v, e * _to_kJ, _td_props);

  switch (_td_props.ireg)
  {
    case (IREG_L):
#ifdef SBTL_USE_C_AUX
      return CV_VU_L_CO2_T(_td_props.vls, e * _to_kJ) * _to_J;
#else
      return CV_VU_L_CO2_T(_td_props.vls,
                           _td_props.dz_l.x1tmin,
                           _td_props.dz_l.x1tmax,
                           _td_props.dz_l.dvdu_vt,
                           _td_props.dz_l.K,
                           e * _to_kJ) *
             _to_J;
#endif
      break;
    case (IREG_G):
#ifdef SBTL_USE_C_AUX
      return CV_VU_G_CO2_T(_td_props.vt, e * _to_kJ) * _to_J;
#else
      return CV_VU_G_CO2_T(_td_props.vt, v, e * _to_kJ) * _to_J;
#endif
    case (IREG_TP):
// cv: the derivative (du/dt)_v is calculable in the two-phase region but its
// slope is discontinuous at dew and bubble curve (mass average is calculated)
#ifdef SBTL_USE_C_AUX
      cv1 = CV_VU_L_CO2_T(_td_props.v1s, _td_props.u1);
      cv2 = CV_VU_G_CO2_T(_td_props.v2t, _td_props.u2);
#else
      cv1 = CV_VU_L_CO2_T(_td_props.v1s,
                          _td_props.dz_1.x1tmin,
                          _td_props.dz_1.x1tmax,
                          _td_props.dz_1.dvdu_vt,
                          _td_props.dz_1.K,
                          _td_props.u1);
      cv2 = CV_VU_G_CO2_T(_td_props.v2t, _td_props.v2, _td_props.u2);
#endif
      return (cv1 + _td_props.x * (cv2 - cv1)) * _to_J;
    default:
      return getNaN();
  }
}

Real
CarbonDioxideHEMFluidProperties::gamma(Real v, Real e) const
{
  return cp(v, e) / cv(v, e);
}

Real
CarbonDioxideHEMFluidProperties::mu(Real v, Real e) const
{
  Real alpha_v, eta1, eta2;

  if (_td_props.GetState(v, e * _to_kJ) < STR_PDP)
    ireg_vu_SBTLCO2(v, e * _to_kJ, _td_props);

  switch (_td_props.ireg)
  {
    case (IREG_L):
      return ETA_VU_L_CO2_T(_td_props.vls, e * _to_kJ);
    case (IREG_G):
      return ETA_VU_G_CO2_T(_td_props.vt, e * _to_kJ);
    case (IREG_TP):
      alpha_v = _td_props.x * _td_props.v2 / v;
      eta1 = ETA_VU_L_CO2_T(_td_props.v1s, _td_props.u1);
      eta2 = ETA_VU_G_CO2_T(_td_props.v2t, _td_props.u2);
      return eta1 + alpha_v * (eta2 - eta1);
    default:
      return getNaN();
  }
}

Real
CarbonDioxideHEMFluidProperties::k(Real v, Real e) const
{
  Real alpha_v, lambda1, lambda2;

  if (_td_props.GetState(v, e * _to_kJ) < STR_PDP)
    ireg_vu_SBTLCO2(v, e * _to_kJ, _td_props);

  switch (_td_props.ireg)
  {
    case (IREG_L):
      return LAMBDA_VU_L_CO2_T(_td_props.vls, e * _to_kJ);
    case (IREG_G):
      return LAMBDA_VU_G_CO2_T(_td_props.vt, e * _to_kJ);
    case (IREG_TP):
      alpha_v = _td_props.x * _td_props.v2 / v;
      lambda1 = LAMBDA_VU_L_CO2_T(_td_props.v1s, _td_props.u1);
      lambda2 = LAMBDA_VU_G_CO2_T(_td_props.v2t, _td_props.u2);
      return lambda1 + alpha_v * (lambda2 - lambda1);
    default:
      return getNaN();
  }
}

Real
CarbonDioxideHEMFluidProperties::rho(Real pressure, Real temperature, Real quality) const
{
  int ierr;
  Real v, vt, u, vl, ul, density;
  Real p = pressure * _to_MPa;
  Real t = temperature;

  if (quality <= 0.0)
  {
    ierr = PT_FLASH_L_CO2(p, t, vl, ul);
    density = 1.0 / vl;
  }
  else if (quality >= 1.0)
  {
    ierr = PT_FLASH_G_CO2(p, t, v, vt, u);
    density = 1.0 / v;
  }
  else
  {
    ierr = PT_FLASH_L_CO2(p, t, vl, ul);
    ierr += PT_FLASH_G_CO2(p, t, v, vt, u);
    density = 1.0 / ((1.0 - quality) * vl + quality * v);
  }

  if (ierr == I_OK)
    return density;
  else
    return getNaN();
}

void
CarbonDioxideHEMFluidProperties::rho_dpT(Real /*pressure*/,
                                         Real /*temperature*/,
                                         Real & /*rho*/,
                                         Real & /*drho_dp*/,
                                         Real & /*drho_dT*/) const
{
  // LibSBTL does not provide a way to compute derivatives of rho wrt p and T
  mooseError("Not Implemented.");
}

Real
CarbonDioxideHEMFluidProperties::beta(Real pressure, Real temperature) const
{
  Real rho, drho_dp, drho_dT;
  rho_dpT(pressure * _to_MPa, temperature, rho, drho_dp, drho_dT);
  return -drho_dT / rho;
}

void
CarbonDioxideHEMFluidProperties::rho_e(Real pressure, Real temperature, Real & rho, Real & e) const
{
  Real v;
  const unsigned int ierr = PT_FLASH_L_CO2(pressure * _to_MPa, temperature, v, e);
  if (ierr != I_OK)
  {
    rho = getNaN();
    e = getNaN();
    return;
  }

  e *= _to_J; // [J/kg]
  rho = 1 / v;
}

Real
CarbonDioxideHEMFluidProperties::e(Real pressure, Real rho) const
{
  Real v = 1. / rho;
  Real p = pressure;
  if (_td_props.GetStatePV(p * _to_MPa, v) < STR_PDP)
    ireg_pv_SBTLCO2(p * _to_MPa, v, _td_props);

  return _td_props.u_ * _to_J;
}

void
CarbonDioxideHEMFluidProperties::dp_duv(
    Real v, Real e, Real & p, Real & dpdv_u, Real & dpdu_v, Real & dudv_p) const
{
  int str_state = _td_props.GetState(v, e * _to_kJ);
  if (str_state < STR_PDP)
    ireg_vu_SBTLCO2(v, e * _to_kJ, _td_props);

  switch (_td_props.ireg)
  {
    case (IREG_L):
      DIFF_P_VU_L_CO2_T(_td_props.vls,
                        _td_props.dz_l.x1tmin,
                        _td_props.dz_l.x1tmax,
                        _td_props.dz_l.dvdu_vt,
                        _td_props.dz_l.K,
                        e * _to_kJ,
                        p,
                        dpdv_u,
                        dpdu_v,
                        dudv_p);
      p *= _to_Pa;
      dpdv_u *= _to_Pa;
      dpdu_v *= (_to_Pa / _to_J);
      dudv_p *= _to_J;
      break;

    case (IREG_G):
      DIFF_P_VU_G_CO2_T(_td_props.vt, v, e * _to_kJ, p, dpdv_u, dpdu_v, dudv_p);
      p *= _to_Pa;
      dpdv_u *= _to_Pa;
      dpdu_v *= (_to_Pa / _to_J);
      dudv_p *= _to_J;
      break;
    case (IREG_TP):
      if (str_state < STR_DTP)
        DIFF_SAT_VU_SPL_CO2(_td_props.ps,
                            _td_props.x,
                            _td_props.v1,
                            _td_props.v1s,
                            _td_props.dz_1.x1tmin,
                            _td_props.dz_1.x1tmax,
                            _td_props.dz_1.dvdu_vt,
                            _td_props.dz_1.K,
                            _td_props.v2,
                            _td_props.v2t,
                            _td_props.u1,
                            _td_props.u2,
                            _td_props.d_tp);
      p = _td_props.ps * _to_Pa;
      dpdv_u = _td_props.d_tp.dpdv_u * _to_Pa;
      dpdu_v = _td_props.d_tp.dpdu_v * (_to_Pa / _to_J);
      dudv_p = _td_props.d_tp.dudv_pt * _to_J;
      break;
  }
}

void
CarbonDioxideHEMFluidProperties::h_dpT(
    Real /*pressure*/, Real /*temperature*/, Real & /*h*/, Real & /*dh_dp*/, Real & /*dh_dT*/) const
{
  // LibSBTL does not provide a way to compute derivatives of h wrt p and T
  mooseError("Not Implemented.");
}

Real
CarbonDioxideHEMFluidProperties::saturation_T(Real pressure) const
{
  return TS_P_CO2(pressure * _to_MPa);
}

Real
CarbonDioxideHEMFluidProperties::saturation_P(Real temperature) const
{
  return PS_T_INV_CO2(temperature) * _to_Pa;
}

Real
CarbonDioxideHEMFluidProperties::surfaceTension(Real temperature) const
{
  return SIGMA_TS_CO2(temperature) * 1e-3;
}

Real
CarbonDioxideHEMFluidProperties::dT_dP_saturation(Real pressure) const
{
  double temperature, dtsdp;
  DIFF_TS_P_CO2(pressure * 1e-6, temperature, dtsdp);
  return dtsdp * 1e-6;
}

Real
CarbonDioxideHEMFluidProperties::p_critical() const
{
  return _P_critical;
}

void
CarbonDioxideHEMFluidProperties::h_lat(Real temperature, Real & hf, Real & hg, Real & hfg) const
{
  static const Real t0 = 216.592;
  static const Real tc = 304.1282;
  static const Real pc = 7.37729837321;
  static const Real vc = 1. / 467.60000128174;
  static const Real uc = 316.468709888;
  static const Real hc = uc + pc * vc * 1.e3;
  Real t = temperature;
  Real ps, vl, ul, vv, vvt, uv, hl, hv;

  if (t < t0 || t > tc)
    mooseError("not in temperature range");
  if ((tc - t) < 0.2)
  {
    hf = hc * _to_J;
    hg = hc * _to_J;
    hfg = 0.;
  }
  else
  {
    ps = PS_T_INV_CO2(t);
    if (PT_FLASH_L_CO2(ps, t, vl, ul))
      mooseError("The temperature provided is not in the valid range of the library");
    if (PT_FLASH_G_CO2(ps, t, vv, vvt, uv))
      mooseError("The temperature provided is not in the valid range of the library");

    hl = (ul + ps * vl * 1.e3) * _to_J;
    hv = (uv + ps * vv * 1.e3) * _to_J;
    hf = hl;
    hg = hv;
    hfg = hg - hf;
  }
}

Real
CarbonDioxideHEMFluidProperties::v_ph(Real pressure, Real enthalpy) const
{
  pressure *= _to_MPa;
  enthalpy *= _to_kJ;
  int str_state = _td_props.GetStatePH(pressure, enthalpy);
  if (str_state < STR_PDP)
    ireg_ph_SBTLCO2(pressure, enthalpy, _td_props);

  switch (_td_props.ireg)
  {
    case (IREG_L):
      return _td_props.v_;
    case (IREG_G):
      return _td_props.v_;
    case (IREG_TP):
      return _td_props.v_;
    default:
      return getNaN();
  }
}

Real
CarbonDioxideHEMFluidProperties::molarMass() const
{
  return 0.0440098;
}

Real
CarbonDioxideHEMFluidProperties::criticalTemperature() const
{
  return 304.1282;
}

Real
CarbonDioxideHEMFluidProperties::criticalDensity() const
{
  return 467.60000128174;
}

Real
CarbonDioxideHEMFluidProperties::criticalInternalEnergy() const
{
  return 316.468709888 * _to_J;
}
