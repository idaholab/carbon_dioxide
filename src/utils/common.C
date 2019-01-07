#include "SBTL_CO2.h"
#include "SBTL_func.h"
#include "common.h"

extern "C" double __stdcall SIGMA_TS_CO2(double t);
extern "C" void __stdcall DIFF_TS_P_CO2(double p, double & ts, double & dtsdp);

Real
T_sat(Real pressure)
{
  return TS_P_CO2(pressure);
}

Real
p_sat(Real temperature)
{
  return PS_T_INV_CO2(temperature);
}

Real
dT_dP_sat(Real pressure)
{
  double temperature, dtsdp;
  DIFF_TS_P_CO2(pressure * 1e-6, temperature, dtsdp);
  return dtsdp * 1e-6;
}

Real
sigma(Real temperature)
{
  return SIGMA_TS_CO2(temperature) * 1e-3;
}
