#include "CarbonDioxideTwoPhaseFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"

TEST_F(CarbonDioxideTwoPhaseFluidPropertiesTest, test)
{
  const Real relative_perturbation = 1e-6;

  Real T = 20. + 273.15; // K
  Real p = 1.e6;         // Pa

  // Tsat + derivatives
  REL_TEST(_fp->T_sat(p), 233.02825123991201, REL_TOL_SAVED_VALUE);
  {
    Real dT_dPsat = _fp->dT_sat_dp(p);

    Real dp = relative_perturbation * p;
    Real dT_dPsat_fd = (_fp->T_sat(p + dp) - _fp->T_sat(p - dp)) / (2 * dp);

    REL_TEST(dT_dPsat, dT_dPsat_fd, REL_TOL_DERIVATIVE);
  }

  // Psat
  REL_TEST(_fp->p_sat(T), 5729052.6071695453, REL_TOL_SAVED_VALUE);
}
