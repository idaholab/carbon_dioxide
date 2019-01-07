#include "CarbonDioxideVaporFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"

TEST_F(CarbonDioxideVaporFluidPropertiesTest, test)
{
  const Real T = 20.0 + 273.15;
  const Real p = 101325;

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
  // TODO: REL_TEST(rho, rho_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(rho_from_p_T, 1.8393470512522803, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->rho_from_p_s(p, s), rho_from_p_T, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->rho_from_p_T, p, T, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->rho_from_p_s, p, s, REL_TOL_DERIVATIVE);

  // e
  // TODO: REL_TEST(e, e_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(e_from_p_rho, 446511.50544600945, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->e_from_v_h(v, h), e, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->e_from_p_rho, p, rho, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->e_from_v_h, v, h, REL_TOL_DERIVATIVE * 10.0);

  // c
  const Real c = _fp->c_from_v_e(v, e);
  // TODO: REL_TEST(c, c_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(c, 266.55544702510394, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->c_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // cp
  const Real cp = _fp->cp_from_v_e(v, e);
  // TODO: REL_TEST(cp, cp_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(cp, 846.00136881550372, REL_TOL_SAVED_VALUE);

  // cv
  const Real cv = _fp->cv_from_v_e(v, e);
  // TODO: REL_TEST(cv, cv_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(cv, 652.40847669234245, REL_TOL_SAVED_VALUE);

  // mu
  const Real mu = _fp->mu_from_v_e(v, e);
  // TODO: REL_TEST(mu, mu_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(mu, 1.4674772387909222e-5, REL_TOL_SAVED_VALUE);

  // k
  const Real k = _fp->k_from_v_e(v, e);
  // TODO: REL_TEST(k, k_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(k, 0.016250508737563615, REL_TOL_SAVED_VALUE);

  // s
  REL_TEST(s, 2722.4590650036321, REL_TOL_EXTERNAL_VALUE);
  // TODO: REL_TEST(s, s_saved, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->s_from_h_p(h, p), s, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->s_from_v_e, v, e, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->s_from_h_p, h, p, REL_TOL_DERIVATIVE);

  // g
  REL_TEST(_fp->g_from_v_e(v, e), -296489.883538, REL_TOL_EXTERNAL_VALUE);
  // TODO: REL_TEST(g, g_saved, REL_TOL_SAVED_VALUE);

  // h
  // TODO: REL_TEST(h, h_external, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(h_from_p_T, 501598.98876299144, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->h_from_p_T, p, T, REL_TOL_DERIVATIVE);

  // beta
  // const Real beta = _fp->beta_from_p_T(p, T);
  // TODO: REL_TEST(beta, beta_external, REL_TOL_EXTERNAL_VALUE);
  // TODO: REL_TEST(beta, beta_saved, REL_TOL_SAVED_VALUE);

  // fluid specific constants
  REL_TEST(_fp->molarMass(), 0.0440098, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->criticalTemperature(), 304.1282, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->criticalDensity(), 467.60000128174, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->criticalInternalEnergy(), 316.468709888e3, REL_TOL_SAVED_VALUE);
}
