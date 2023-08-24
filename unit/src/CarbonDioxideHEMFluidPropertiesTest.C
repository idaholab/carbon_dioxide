//* This file is part of carbon_dioxide
//* https://github.com/idaholab/carbon_dioxide
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/carbon_dioxide/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CarbonDioxideHEMFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"

TEST_F(CarbonDioxideHEMFluidPropertiesTest, general)
{
  REL_TEST(_fp->molarMass(), 0.0440098, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->criticalTemperature(), 304.1282, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->criticalDensity(), 467.60000128174, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->criticalInternalEnergy(), 316.468709888e3, REL_TOL_SAVED_VALUE);
}

TEST_F(CarbonDioxideHEMFluidPropertiesTest, test_2phase)
{
  const Real relative_perturbation = 1e-6;

  const Real v = 0.00662481817619;      // m3/kg
  const Real u = 290.307769134 * 1000; // J/kg

  Real p = _fp->pressure(v, u);
  REL_TEST(p, 2.9999980569563876 * 1e6, REL_TOL_SAVED_VALUE); // [Pa]
  Real T = _fp->temperature(v, u);
  REL_TEST(T, 267.59784757383147, REL_TOL_SAVED_VALUE); // [K]
  Real quality = _fp->quality(v, u);
  REL_TEST(quality, 0.50000019015471786, REL_TOL_SAVED_VALUE);
  Real quality2 = _fp->quality_Tsat_h(T, _fp->h(p, T, quality));
  REL_TEST(quality2, 0.50000019015471797, REL_TOL_SAVED_VALUE);

  REL_TEST(_fp->c(v, u), 164.8924365638984, REL_TOL_SAVED_VALUE);        // [m/s]
  REL_TEST(_fp->cp(v, u), 2.0199580442214963 * 1e3, REL_TOL_SAVED_VALUE); // [J/(kg-K)]
  REL_TEST(_fp->cv(v, u), 0.88670680333081452 * 1e3, REL_TOL_SAVED_VALUE); // [J/(kg-K)]
  REL_TEST(_fp->gamma(v, u), 2.2780450501042178, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->mu(v, u), 2.1611284488617797e-05, REL_TOL_SAVED_VALUE); // [Pa-s]
  REL_TEST(_fp->k(v, u), 0.026285182362021972, REL_TOL_SAVED_VALUE);     // [W/ (m-K)]

  Real rho = _fp->rho(p, T, quality);
  REL_TEST(rho, 150.94753899716585, REL_TOL_SAVED_VALUE);
  Real rho_l = _fp->rho(p, T, 0.0);
  REL_TEST(rho_l, 959.25259503940947, REL_TOL_SAVED_VALUE);
  Real rho_v = _fp->rho(p, T, 1.0);
  REL_TEST(rho_v, 81.919177444716638, REL_TOL_EXTERNAL_VALUE);

  REL_TEST(_fp->alpha_vapor(v, u), (rho - rho_l) / (rho_v - rho_l), REL_TOL_SAVED_VALUE);
  Real hl = 0.0;
  Real hv = 0.0;
  Real hvap = 0.0;
  _fp->h_lat(T, hl, hv, hvap);
  REL_TEST(hl, 186.753 * 1e3, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(hv, 433.611 * 1e3, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(hvap, 246.858 * 1e3, REL_TOL_EXTERNAL_VALUE);

  REL_TEST(_fp->saturation_T(p), 267.59784757379242, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->saturation_P(T), 2.9999980569596114 * 1e6, REL_TOL_EXTERNAL_VALUE);

  Real e;
  // rho_e gives rho and e at quality = 0.0
  _fp->rho_e(p, T, rho, e);
  REL_TEST(rho, 959.25259503940947, REL_TOL_SAVED_VALUE);
  REL_TEST(e, 183.62621014317023 * 1000, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->e(p, rho), 183.62621014317023 * 1000, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->h(p, T, quality), 310.1822* 1e3, REL_TOL_EXTERNAL_VALUE); // [J/kg]

  {
    Real dp = relative_perturbation * p;
    Real dT_dPsat = _fp->dT_dP_saturation(p);
    Real dT_dPsat_fd = (_fp->saturation_T(p + dp) - _fp->saturation_T(p - dp)) / (2 * dp);
    REL_TEST(dT_dPsat, dT_dPsat_fd, REL_TOL_DERIVATIVE);
  }

  {
    Real pressure = 1e6;

    Real de = relative_perturbation * e;
    Real dv = relative_perturbation * v;

    Real dp_dv_fd = (_fp->pressure(v + dv, e) - _fp->pressure(v - dv, e)) / (2 * dv);
    Real dp_de_fd = (_fp->pressure(v, e + de) - _fp->pressure(v, e - de)) / (2 * de);

    Real dpdv = 0, dpdu = 0, dudvp = 0;
    _fp->dp_duv(v, e, pressure, dpdv, dpdu, dudvp);

    REL_TEST(dpdv, dp_dv_fd, REL_TOL_DERIVATIVE);
    REL_TEST(dpdu, dp_de_fd, REL_TOL_DERIVATIVE);
  }
}

TEST_F(CarbonDioxideHEMFluidPropertiesTest, test_gas)
{
  const Real relative_perturbation = 1e-6;

  const Real v = 0.00609178146204;        // m3/kg
  const Real u = 438.152985295 * 1000; // J/kg

  REL_TEST(_fp->pressure(v, u), 7.9999712023691693 * 1e6, REL_TOL_SAVED_VALUE); // [Pa]
  REL_TEST(_fp->temperature(v, u), 349.99943332804094, REL_TOL_SAVED_VALUE);    // [K]
  REL_TEST(_fp->quality(v, u), -1000.0, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->alpha_vapor(v, u), 1.0, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->c(v, u), 254.0979902013965, REL_TOL_SAVED_VALUE);        // [m/s]
  REL_TEST(_fp->cp(v, u), 1.546540207458901 * 1e3, REL_TOL_SAVED_VALUE);  // [J/(kg-K)]
  REL_TEST(_fp->cv(v, u), 0.83103207277429578 * 1e3, REL_TOL_SAVED_VALUE); // [J/(kg-K)]
  REL_TEST(_fp->gamma(v, u), 1.8609873892062572, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->mu(v, u), 2.0082409716418526e-05, REL_TOL_SAVED_VALUE); // [Pa-s]
  REL_TEST(_fp->k(v, u), 0.028847493483532941, REL_TOL_SAVED_VALUE);     // [W/ (m-K)]

  Real p = _fp->pressure(v, u);
  Real T = _fp->temperature(v, u);
  // By sending in 1.0 as quality, we indicate that this is a gas
  REL_TEST(_fp->h(p, T, 1.0), 486.88706156244641 * 1e3, REL_TOL_SAVED_VALUE); // [J/kg]
  Real rho = _fp->rho(p, T, 1.0);
  REL_TEST(rho, 164.15559327453656, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->e(p, rho), u, REL_TOL_CONSISTENCY);

  T = 300.;
  Real h = 400e3; // J/kg
  Real quality2 = _fp->quality_Tsat_h(T, h);
  REL_TEST(quality2, 1.1245819276180786, REL_TOL_EXTERNAL_VALUE);

  {
    Real e = 438.152985295 * 1e3; // J/kg

    Real de = relative_perturbation * e;
    Real dv = relative_perturbation * v;

    Real dp_dv_fd = (_fp->pressure(v + dv, e) - _fp->pressure(v - dv, e)) / (2 * dv);
    Real dp_de_fd = (_fp->pressure(v, e + de) - _fp->pressure(v, e - de)) / (2 * de);

    Real dpdv = 0, dpdu = 0, dudvp = 0;
    _fp->dp_duv(v, e, p, dpdv, dpdu, dudvp);

    REL_TEST(dpdv, dp_dv_fd, REL_TOL_DERIVATIVE);
    REL_TEST(dpdu, dp_de_fd, REL_TOL_DERIVATIVE);
  }
}

TEST_F(CarbonDioxideHEMFluidPropertiesTest, test_liquid)
{
  const Real relative_perturbation = 1e-6;

  Real v = 0.00132772602326;      // m3/kg
  Real u = 259.336330993 * 1000; // J/kg

  REL_TEST(_fp->pressure(v, u), 7.9999989132258333 * 1e6, REL_TOL_SAVED_VALUE); // [Pa]
  REL_TEST(_fp->temperature(v, u), 300.00000245824401, REL_TOL_SAVED_VALUE);     // [K]
  REL_TEST(_fp->quality(v, u), -1000.0, REL_TOL_SAVED_VALUE);
  ABS_TEST(_fp->alpha_vapor(v, u), 0.0, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->c(v, u), 343.66047365374646, REL_TOL_SAVED_VALUE);        // [m/s]
  REL_TEST(_fp->cp(v, u), 3.9320875425750659 * 1e3, REL_TOL_SAVED_VALUE); // [J/(kg-K)]
  REL_TEST(_fp->cv(v, u), 0.98935273335062123 * 1e3, REL_TOL_SAVED_VALUE); // [J/(kg-K)]
  REL_TEST(_fp->gamma(v, u), 3.974404082615048, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->mu(v, u), 0.00006368640046921546, REL_TOL_SAVED_VALUE); // [Pa-s]
  REL_TEST(_fp->k(v, u), 0.082402329789538362, REL_TOL_SAVED_VALUE);     // [W/ (m-K)]

  Real p = _fp->pressure(v, u);
  Real T = _fp->temperature(v, u);
  // By sending in 0.0 as quality, we indicate that this is a liquid
  REL_TEST(_fp->h(p, T, 0.0), 269.95813773614151 * 1e3, REL_TOL_SAVED_VALUE); // [J/kg]
  Real rho = _fp->rho(p, T, 0.0);
  REL_TEST(rho, 753.16743249836645, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->e(p, rho), u, REL_TOL_CONSISTENCY);
  p = 7.e6;
  REL_TEST(_fp->saturation_T(p), 301.8325152953297, REL_TOL_EXTERNAL_VALUE);

  {
    v = 0.00132772602326;                           // m3/kg
    Real e = 259.336330993 * 1e3;                  // J/kg
    Real pressure = _fp->pressure(v, u); // Pa

    Real de = relative_perturbation * e;
    Real dv = relative_perturbation * v;

    Real dp_dv_fd = (_fp->pressure(v + dv, e) - _fp->pressure(v - dv, e)) / (2 * dv);
    Real dp_de_fd = (_fp->pressure(v, e + de) - _fp->pressure(v, e - de)) / (2 * de);

    Real dpdv = 0, dpdu = 0, dudvp = 0;
    _fp->dp_duv(v, e, pressure, dpdv, dpdu, dudvp);

    REL_TEST(dpdv, dp_dv_fd, REL_TOL_DERIVATIVE);
    REL_TEST(dpdu, dp_de_fd, REL_TOL_DERIVATIVE);
  }
}
