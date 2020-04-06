#pragma once

#include "HEMFluidProperties.h"
#include "NaNInterface.h"
#include "contrib/libSBTL_CarbonDioxide/SBTL_CO2.h"
#include "contrib/libSBTL_CarbonDioxide/SBTL_flags.h"

/**
 * Carbon Dioxide HEM fluid properties
 *
 * Note on units used in libSBTL_CarbonDioxide:
 * - pressure: MPa
 * - internal energy: kJ/kg
 * - density: kg/m3
 * - speed of sound: m/s
 */
class CarbonDioxideHEMFluidProperties : public HEMFluidProperties, public NaNInterface
{

public:
  CarbonDioxideHEMFluidProperties(const InputParameters & parameters);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

  virtual Real pressure(Real v, Real e) const override;
  virtual Real temperature(Real v, Real e) const override;
  virtual Real quality(Real v, Real e) const override;
  virtual Real quality_Tsat_h(Real Tsat, Real h) const override;
  virtual Real alpha_vapor(Real v, Real e) const override;
  virtual Real c(Real v, Real e) const override;
  virtual Real cp(Real v, Real e) const override;
  virtual Real cv(Real v, Real e) const override;
  virtual Real gamma(Real v, Real e) const override;
  virtual Real mu(Real v, Real e) const override;
  virtual Real k(Real v, Real e) const override;
  virtual Real beta(Real pressure, Real temperature) const override;
  virtual Real rho(Real pressure, Real temperature, Real quality) const override;
  virtual void rho_dpT(
      Real pressure, Real temperature, Real & rho, Real & drho_dp, Real & drho_dT) const override;
  virtual void rho_e(Real pressure, Real temperature, Real & rho, Real & e) const override;
  virtual Real e(Real pressure, Real rho) const override;
  virtual void
  dp_duv(Real v, Real u, Real & p, Real & dpdv_u, Real & dpdu_v, Real & dudv_p) const override;
  virtual Real saturation_T(Real pressure) const override;
  virtual Real saturation_P(Real temperature) const override;
  virtual Real surfaceTension(Real temperature) const override;
  virtual Real dT_dP_saturation(Real pressure) const override;
  virtual Real h(Real pressure, Real temperature, Real quality) const override;
  virtual void
  h_dpT(Real pressure, Real temperature, Real & h, Real & dh_dp, Real & dh_dT) const override;
  virtual Real p_critical() const override;
  virtual void h_lat(Real temperature, Real & hf, Real & hg, Real & hfg) const override;
  virtual Real v_ph(Real pressure, Real enthalpy) const override;
  virtual Real molarMass() const override;
  virtual Real criticalTemperature() const override;
  virtual Real criticalDensity() const override;
  virtual Real criticalInternalEnergy() const override;

#pragma GCC diagnostic pop

private:
  /// Thermodynamic Properties
  mutable STR_vu_SBTL_CO2 _td_props;
  /// Sound Speed Material Flag
  mutable SOUND_2PHASE _ss_flag;

protected:
  /// Conversion factor from Pa to MPa
  const Real _to_MPa;
  /// Conversion factor from MPa to Pa
  const Real _to_Pa;
  /// Conversion factor from J to kJ
  const Real _to_kJ;
  /// Conversion factor from kJ to J
  const Real _to_J;

protected:
  // Critical pressure
  static const Real _P_critical;

public:
  static InputParameters validParams();
};
