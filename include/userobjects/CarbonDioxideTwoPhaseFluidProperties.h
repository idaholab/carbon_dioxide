#pragma once

#include "TwoPhaseFluidProperties.h"
#include "NaNInterface.h"

class CarbonDioxideTwoPhaseFluidProperties;
class SinglePhaseFluidProperties;

template <>
InputParameters validParams<CarbonDioxideTwoPhaseFluidProperties>();

/**
 * Two-phase carbon dioxide fluid properties
 *
 * Range of validity:
 *   0.0005 MPa <= p <= 100 MPa
 *   T_triple (216.59 K) <= T <= 1300 K
 */
class CarbonDioxideTwoPhaseFluidProperties : public TwoPhaseFluidProperties, public NaNInterface
{
public:
  CarbonDioxideTwoPhaseFluidProperties(const InputParameters & parameters);

  virtual Real p_critical() const override;
  virtual Real T_sat(Real pressure) const override;
  virtual Real p_sat(Real temperature) const override;
  virtual Real dT_sat_dp(Real pressure) const override;

  virtual bool supportsPhaseChange() const override { return true; }

protected:
  // Critical pressure
  static const Real _P_critical;
};