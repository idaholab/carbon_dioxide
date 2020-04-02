#include "CarbonDioxideLiquidFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"

TEST_F(CarbonDioxideLiquidFluidPropertiesTest, test)
{
  const Real T = 273.15;
  const Real p = 5.e6;

  const Real rho_from_p_T = _fp->rho_from_p_T(p, T);
  const Real rho = rho_from_p_T;

  const Real h_from_p_T = _fp->h_from_p_T(p, T);
  const Real h = h_from_p_T;

  const Real e_from_p_rho = _fp->e_from_p_rho(p, rho);
  const Real e = e_from_p_rho;

  const Real v = 1 / rho;

  const Real s_from_v_e = _fp->s_from_v_e(v, e);
  const Real s = s_from_v_e;

  // p
  REL_TEST(_fp->p_from_v_e(v, e), p, REL_TOL_CONSISTENCY);
  REL_TEST(_fp->p_from_h_s(h, s), p, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->p_from_v_e, v, e, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->p_from_h_s, h, s, REL_TOL_DERIVATIVE);

  // T
  REL_TEST(_fp->T_from_v_e(v, e), T, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->T_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // rho
  // REL_TEST(rho, rho_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(rho_from_p_T, 940.51698090428886, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->rho_from_p_s(p, s), rho_from_p_T, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->rho_from_p_T, p, T, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->rho_from_p_s, p, s, REL_TOL_DERIVATIVE * 10.0);

  // e
  // REL_TEST(e, e_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(e_from_p_rho, 193217.29762783623, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->e_from_v_h(v, h), e, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->e_from_p_rho, p, rho, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->e_from_v_h, v, h, REL_TOL_DERIVATIVE);

  // c
  const Real c = _fp->c_from_v_e(v, e);
  // TODO: REL_TEST(c, c_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(c, 567.7165137694559, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->c_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // cp
  const Real cp = _fp->cp_from_v_e(v, e);
  // TODO: REL_TEST(cp, cp_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(cp, 2417.2744155462779, REL_TOL_SAVED_VALUE);

  // cv
  const Real cv = _fp->cv_from_v_e(v, e);
  // TODO: REL_TEST(cv, cv_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(cv, 935.30996220147495, REL_TOL_SAVED_VALUE);

  // mu
  Real mu = _fp->mu_from_v_e(v, e);
  // TODO: REL_TEST(mu, mu_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(mu, 0.00010422532124202926, REL_TOL_SAVED_VALUE);
  Real dmu_dv, dmu_de;
  _fp->mu_from_v_e(v, e, mu, dmu_dv, dmu_de);
  REL_TEST(mu, 0.00010422532124202926, REL_TOL_SAVED_VALUE);
  REL_TEST(dmu_dv, 0., REL_TOL_SAVED_VALUE);
  REL_TEST(dmu_de, 0., REL_TOL_SAVED_VALUE);

  // k
  const Real k = _fp->k_from_v_e(v, e);
  // TODO: REL_TEST(k, k_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(k, 0.11188889323403757, REL_TOL_SAVED_VALUE);

  // s
  REL_TEST(s, 988.69426684335974, REL_TOL_EXTERNAL_VALUE);
  // TODO: REL_TEST(s, s_saved, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->s_from_h_p(h, p), s, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->s_from_v_e, v, e, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->s_from_h_p, h, p, REL_TOL_DERIVATIVE);

  // g
  const Real g = _fp->g_from_v_e(v, e);
  REL_TEST(g, -71.528316240076441e+03, REL_TOL_EXTERNAL_VALUE);
  // TODO: REL_TEST(g, g_saved, REL_TOL_SAVED_VALUE);

  // h
  // TODO: REL_TEST(h, h_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(h_from_p_T, 198533.52274818713, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->h_from_p_T, p, T, REL_TOL_DERIVATIVE * 10.0);

  // beta
  // const Real beta = _fp->beta_from_p_T(p, T);
  // TODO: REL_TEST(beta, beta_external, REL_TOL_EXTERNAL_VALUE);
  // TODO: REL_TEST(beta, beta_saved, REL_TOL_SAVED_VALUE);

  // sigma
  const Real sigma = _fp->sigma_from_p_T(p, 290.);
  REL_TEST(sigma, 0.0016750845120663005, REL_TOL_SAVED_VALUE);

  // fluid specific constants
  REL_TEST(_fp->molarMass(), 0.0440098, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->criticalTemperature(), 304.1282, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->criticalDensity(), 467.60000128174, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->criticalInternalEnergy(), 316.468709888e3, REL_TOL_SAVED_VALUE);
}
