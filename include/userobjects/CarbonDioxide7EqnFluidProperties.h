#ifndef CARBONDIOXIDE7EQNFLUIDPROPERTIES_H
#define CARBONDIOXIDE7EQNFLUIDPROPERTIES_H

#include "TwoPhaseFluidProperties.h"
#include "NaNInterface.h"

class CarbonDioxide7EqnFluidProperties;
class SinglePhaseFluidProperties;

template <>
InputParameters validParams<CarbonDioxide7EqnFluidProperties>();

/**
 * CarbonDioxide interface for 7-eqn model
 */
class CarbonDioxide7EqnFluidProperties : public TwoPhaseFluidProperties, public NaNInterface
{
public:
  CarbonDioxide7EqnFluidProperties(const InputParameters & parameters);

  virtual Real p_critical() const override;
  virtual Real T_sat(Real pressure) const override;
  virtual Real p_sat(Real temperature) const override;
  virtual Real dT_sat_dp(Real pressure) const override;

  virtual bool supportsPhaseChange() const override { return true; }

protected:
  // Critical pressure
  static const Real _P_critical;
};

#endif /* CARBONDIOXIDE7EQNFLUIDPROPERTIES_H */
