#include "CarbonDioxideLiquidFluidProperties.h"
#include "contrib/libSBTL_CarbonDioxide/SBTL_CO2.h"

extern "C" double P_VU_L_CO2(double v, double u);
extern "C" double T_VU_L_CO2(double v, double u);
extern "C" int PT_FLASH_L_CO2(double p, double t, double & v, double & u);
extern "C" int PT_FLASH_DERIV_L_CO2(double p,
                                    double t,
                                    double & v,
                                    double & dvdp_t,
                                    double & dvdt_p,
                                    double & dpdt_v,
                                    double & u,
                                    double & dudp_t,
                                    double & dudt_p,
                                    double & dpdt_u);
extern "C" int PS_FLASH_L_CO2(double p, double s, double & v, double & u);
extern "C" void PS_FLASH_DERIV_L_CO2(double v,
                                     double u,
                                     double & dvdp_s,
                                     double & dvds_p,
                                     double & dpds_v,
                                     double & dudp_s,
                                     double & duds_p,
                                     double & dpds_u);
extern "C" double W_VU_L_CO2(double v, double u);
extern "C" double CP_VU_L_CO2(double v, double u);
extern "C" double CV_VU_L_CO2(double v, double u);
extern "C" double ETA_VU_L_CO2(double v, double u);
extern "C" double LAMBDA_VU_L_CO2(double v, double u);
extern "C" double S_VU_L_CO2(double v, double u);
extern "C" double G_VU_L_CO2(double v, double e);
extern "C" double U_VP_L_CO2(double v, double p);
extern "C" int HS_FLASH_L_CO2(double h, double s, double & v, double & u);
extern "C" void HS_FLASH_DERIV_L_CO2(double v,
                                     double u,
                                     double & dvdh_s,
                                     double & dvds_h,
                                     double & dhds_v,
                                     double & dudh_s,
                                     double & duds_h,
                                     double & dhds_u);
extern "C" void
DIFF_U_VP_L_CO2(double v, double p, double & u, double & dudv_p, double & dudp_v, double & dpdv_u);
extern "C" double SIGMA_TS_CO2(double t);

// SBTL functions with derivatives
extern "C" void
DIFF_P_VU_L_CO2(double v, double u, double & p, double & dpdv, double & dpdu, double & dudv);
extern "C" void
DIFF_T_VU_L_CO2(double v, double u, double & t, double & dtdv, double & dtdu, double & dudv);
extern "C" void
DIFF_S_VU_L_CO2(double v, double u, double & s, double & dsdv, double & dsdu, double & dudv);
extern "C" void
DIFF_W_VU_L_CO2(double v, double u, double & w, double & dwdv, double & dwdu, double & dudv);

extern "C" int VH_FLASH_L_CO2(double v, double h, double & u);
extern "C" int PH_FLASH_L_CO2(double p, double h, double & v, double & u);

registerMooseObject("CarbonDioxideApp", CarbonDioxideLiquidFluidProperties);

template <>
InputParameters
validParams<CarbonDioxideLiquidFluidProperties>()
{
  InputParameters params = validParams<SinglePhaseFluidProperties>();
  params += validParams<NaNInterface>();
  params.addClassDescription("Fluid properties of liquid carbon dioxide (meta stable).");
  return params;
}

CarbonDioxideLiquidFluidProperties::CarbonDioxideLiquidFluidProperties(
    const InputParameters & parameters)
  : SinglePhaseFluidProperties(parameters),
    LiquidFluidPropertiesInterface(),
    NaNInterface(this),
    _to_MPa(1e-6),
    _to_Pa(1e6),
    _to_kJ(1e-3),
    _to_J(1e3),
    _to_N(1e-3)
{
}

Real
CarbonDioxideLiquidFluidProperties::p_from_v_e(Real v, Real e) const
{
  return P_VU_L_CO2(v, e * _to_kJ) * _to_Pa;
}

void
CarbonDioxideLiquidFluidProperties::p_from_v_e(
    Real v, Real e, Real & p, Real & dp_dv, Real & dp_de) const
{
  double de_dv_p;

  e *= _to_kJ;
  DIFF_P_VU_L_CO2(v, e, p, dp_dv, dp_de, de_dv_p);
  p *= _to_Pa;
  dp_dv *= _to_Pa;
  dp_de *= _to_Pa / _to_J;
}

Real
CarbonDioxideLiquidFluidProperties::T_from_v_e(Real v, Real e) const
{
  return T_VU_L_CO2(v, e * _to_kJ);
}

void
CarbonDioxideLiquidFluidProperties::T_from_v_e(
    Real v, Real e, Real & T, Real & dT_dv, Real & dT_de) const
{
  double de_dv_T;

  e *= _to_kJ;
  DIFF_T_VU_L_CO2(v, e, T, dT_dv, dT_de, de_dv_T);
  dT_de /= _to_J;
}

Real
CarbonDioxideLiquidFluidProperties::c_from_v_e(Real v, Real e) const
{
  return W_VU_L_CO2(v, e * _to_kJ);
}

void
CarbonDioxideLiquidFluidProperties::c_from_v_e(
    Real v, Real e, Real & c, Real & dc_dv, Real & dc_de) const
{
  double de_dv_c;

  e *= _to_kJ;
  DIFF_W_VU_L_CO2(v, e, c, dc_dv, dc_de, de_dv_c);
  dc_de /= _to_J;
}

Real
CarbonDioxideLiquidFluidProperties::e_from_v_h(Real v, Real h) const
{
  double e;
  VH_FLASH_L_CO2(v, h * _to_kJ, e);
  return e * _to_J;
}

void
CarbonDioxideLiquidFluidProperties::e_from_v_h(
    Real v, Real h, Real & e, Real & de_dv, Real & de_dh) const
{
  e = e_from_v_h(v, h);

  Real p, dp_dv, dp_de;
  p_from_v_e(v, e, p, dp_dv, dp_de);

  Real dv_dh = 1. / (p + dp_dv * v);
  de_dh = 1. / (1. + dp_de * v);
  de_dv = -de_dh / dv_dh;
}

Real
CarbonDioxideLiquidFluidProperties::cp_from_v_e(Real v, Real e) const
{
  return CP_VU_L_CO2(v, e * _to_kJ) * _to_J;
}

void
CarbonDioxideLiquidFluidProperties::cp_from_v_e(
    Real v, Real e, Real & cp, Real & dcp_dv, Real & dcp_de) const
{
  cp = cp_from_v_e(v, e);
  dcp_dv = 0;
  dcp_de = 0;
}

Real
CarbonDioxideLiquidFluidProperties::cv_from_v_e(Real v, Real e) const
{
  return CV_VU_L_CO2(v, e * _to_kJ) * _to_J;
}

Real
CarbonDioxideLiquidFluidProperties::mu_from_v_e(Real v, Real e) const
{
  return ETA_VU_L_CO2(v, e * _to_kJ);
}

void
CarbonDioxideLiquidFluidProperties::mu_from_v_e(
    Real v, Real e, Real & mu, Real & dmu_dv, Real & dmu_de) const
{
  mu = mu_from_v_e(v, e);
  // currently, there is no API for the derivatives in the SBTL package
  dmu_de = 0;
  dmu_dv = 0;
}

Real
CarbonDioxideLiquidFluidProperties::k_from_v_e(Real v, Real e) const
{
  return LAMBDA_VU_L_CO2(v, e * _to_kJ);
}

Real
CarbonDioxideLiquidFluidProperties::s_from_v_e(Real v, Real e) const
{
  return S_VU_L_CO2(v, e * _to_kJ) * _to_J;
}

void
CarbonDioxideLiquidFluidProperties::s_from_v_e(
    Real v, Real e, Real & s, Real & ds_dv, Real & ds_de) const
{
  double de_dv_s;

  e *= _to_kJ;
  DIFF_S_VU_L_CO2(v, e, s, ds_dv, ds_de, de_dv_s);
  s *= _to_J;
  ds_dv *= _to_J;
  ds_de *= _to_J / _to_J;
}

Real
CarbonDioxideLiquidFluidProperties::s_from_h_p(Real h, Real p) const
{
  double v, e;
  PH_FLASH_L_CO2(p * _to_MPa, h * _to_kJ, v, e);
  return s_from_v_e(v, e * _to_J);
}

void
CarbonDioxideLiquidFluidProperties::s_from_h_p(
    Real h, Real p, Real & s, Real & ds_dh, Real & ds_dp) const
{
  double v, e;
  PH_FLASH_L_CO2(p * _to_MPa, h * _to_kJ, v, e);
  e *= _to_J;
  Real dummy_p, dp_dv, dp_de;
  p_from_v_e(v, e, dummy_p, dp_dv, dp_de);
  Real dh_dv = (p + dp_dv * v);
  Real dh_de = 1. + dp_de * v;
  Real ds_dv, ds_de;
  s_from_v_e(v, e, s, ds_dv, ds_de);
  ds_dp = (ds_dv * dh_de - ds_de * dh_dv) / (dp_dv * dh_de - dp_de * dh_dv);
  ds_dh = (ds_dv * dp_de - ds_de * dp_dv) / (dh_dv * dp_de - dh_de * dp_dv);
}

Real
CarbonDioxideLiquidFluidProperties::beta_from_p_T(Real p, Real T) const
{
  double rho, drho_dp, drho_dT;
  rho_from_p_T(p, T, rho, drho_dp, drho_dT);
  return -drho_dT / rho;
}

void
CarbonDioxideLiquidFluidProperties::beta_from_p_T(
    Real p, Real T, Real & beta, Real & dbeta_dp, Real & dbeta_dT) const
{
  beta = beta_from_p_T(p, T);
  dbeta_dT = 0;
  dbeta_dp = 0;
}

Real
CarbonDioxideLiquidFluidProperties::rho_from_p_T(Real p, Real T) const
{
  double v, e;
  const unsigned int ierr = PT_FLASH_L_CO2(p * _to_MPa, T, v, e);
  if (ierr != I_OK)
    return getNaN();

  return 1 / v;
}

void
CarbonDioxideLiquidFluidProperties::rho_from_p_T(
    Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  double v, dv_dp, dv_dT, dp_dT_v;
  double e, de_dp, de_dT, dp_dT_e;
  const unsigned int ierr =
      PT_FLASH_DERIV_L_CO2(p * _to_MPa, T, v, dv_dp, dv_dT, dp_dT_v, e, de_dp, de_dT, dp_dT_e);
  if (ierr != I_OK)
  {
    rho = getNaN();
    drho_dp = getNaN();
    drho_dT = getNaN();
    return;
  }

  rho = 1. / v;
  const double drho_dv = -1. / v / v;
  drho_dp = drho_dv * dv_dp / _to_Pa;
  drho_dT = drho_dv * dv_dT;
}

Real
CarbonDioxideLiquidFluidProperties::e_from_p_rho(Real p, Real rho) const
{
  double v = 1.0 / rho;
  return U_VP_L_CO2(v, p * _to_MPa) * _to_J;
}

void
CarbonDioxideLiquidFluidProperties::e_from_p_rho(
    Real p, Real rho, Real & e, Real & de_dp, Real & de_drho) const
{
  const double v = 1.0 / rho;
  const double dv_drho = -1.0 / rho / rho;

  double de_dv, dp_dv_e;
  DIFF_U_VP_L_CO2(v, p * _to_MPa, e, de_dv, de_dp, dp_dv_e);

  e *= _to_J;
  de_dp *= _to_J / _to_Pa;
  de_drho = de_dv * _to_J * dv_drho;
}

Real
CarbonDioxideLiquidFluidProperties::h_from_p_T(Real p, Real T) const
{
  double v, e;
  const unsigned int ierr = PT_FLASH_L_CO2(p * _to_MPa, T, v, e);
  if (ierr != I_OK)
    return getNaN();

  return e * _to_J + p * v;
}

void
CarbonDioxideLiquidFluidProperties::h_from_p_T(
    Real p, Real T, Real & h, Real & dh_dp, Real & dh_dT) const
{
  double v, dvdp_T, dvdT_p, dpdT_v;
  double e, dedp_T, dedT_p, dpdT_e;
  const unsigned int ierr =
      PT_FLASH_DERIV_L_CO2(p * _to_MPa, T, v, dvdp_T, dvdT_p, dpdT_v, e, dedp_T, dedT_p, dpdT_e);
  if (ierr != I_OK)
  {
    h = getNaN();
    dh_dp = getNaN();
    dh_dT = getNaN();
    return;
  }

  h = e * _to_J + p * v;
  dh_dp = dedp_T * _to_J / _to_Pa + (v + p * _to_MPa * dvdp_T);
  dh_dT = dedT_p * _to_J + p * dvdT_p;
}

Real
CarbonDioxideLiquidFluidProperties::p_from_h_s(Real h, Real s) const
{
  double v, e;
  HS_FLASH_L_CO2(h * _to_kJ, s * _to_kJ, v, e);
  return P_VU_L_CO2(v, e) * _to_Pa;
}

void
CarbonDioxideLiquidFluidProperties::p_from_h_s(
    Real h, Real s, Real & p, Real & dp_dh, Real & dp_ds) const
{
  double v, e;
  HS_FLASH_L_CO2(h * _to_kJ, s * _to_kJ, v, e);

  double dv_dh, dv_ds, dh_ds_v, de_dh, de_ds, dh_ds_e;
  HS_FLASH_DERIV_L_CO2(v, e, dv_dh, dv_ds, dh_ds_v, de_dh, de_ds, dh_ds_e);

  double dp_dv, dp_de;
  p_from_v_e(v, e * _to_J, p, dp_dv, dp_de);

  dp_dh = dp_dv * dv_dh / _to_J + dp_de * de_dh;
  dp_ds = dp_dv * dv_ds / _to_J + dp_de * de_ds;
}

Real
CarbonDioxideLiquidFluidProperties::g_from_v_e(Real v, Real e) const
{
  return G_VU_L_CO2(v, e * _to_kJ) * _to_J;
}

Real
CarbonDioxideLiquidFluidProperties::rho_from_p_s(Real p, Real s) const
{
  double v, e;
  const unsigned int ierr = PS_FLASH_L_CO2(p * _to_MPa, s * _to_kJ, v, e);
  if (ierr != I_OK)
    return getNaN();

  return 1.0 / v;
}

void
CarbonDioxideLiquidFluidProperties::rho_from_p_s(
    Real p, Real s, Real & rho, Real & drho_dp, Real & drho_ds) const
{
  double v, e;
  const unsigned int ierr = PS_FLASH_L_CO2(p * _to_MPa, s * _to_kJ, v, e);
  if (ierr != I_OK)
  {
    rho = getNaN();
    drho_dp = getNaN();
    drho_ds = getNaN();
    return;
  }

  double dv_dp, dv_ds, dp_ds_v, de_dp, de_ds, dp_ds_e;
  PS_FLASH_DERIV_L_CO2(v, e, dv_dp, dv_ds, dp_ds_v, de_dp, de_ds, dp_ds_e);

  rho = 1. / v;
  const double drho_dv = -1. / v / v;
  drho_dp = drho_dv * dv_dp / _to_Pa;
  drho_ds = drho_dv * dv_ds / _to_J;
}

Real
CarbonDioxideLiquidFluidProperties::sigma_from_p_T(Real, Real T) const
{
  return SIGMA_TS_CO2(T) * _to_N;
}

Real
CarbonDioxideLiquidFluidProperties::molarMass() const
{
  return 0.0440098;
}

Real
CarbonDioxideLiquidFluidProperties::criticalTemperature() const
{
  return 304.1282;
}

Real
CarbonDioxideLiquidFluidProperties::criticalDensity() const
{
  return 467.60000128174;
}

Real
CarbonDioxideLiquidFluidProperties::criticalInternalEnergy() const
{
  return 316.468709888 * _to_J;
}
