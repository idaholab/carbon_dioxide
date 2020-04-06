#include "CarbonDioxideTwoPhaseFluidProperties.h"
#include "SinglePhaseFluidProperties.h"

extern "C" double TS_P_CO2(double p);
extern "C" void DIFF_TS_P_CO2(double x1_val, double & ts, double & dtsdp);
extern "C" double PS_T_INV_CO2(double t);

const Real CarbonDioxideTwoPhaseFluidProperties::_P_critical = 7.37729837321E+6;

registerMooseObject("CarbonDioxideApp", CarbonDioxideTwoPhaseFluidProperties);

InputParameters
CarbonDioxideTwoPhaseFluidProperties::validParams()
{
  InputParameters params = TwoPhaseFluidProperties::validParams();
  params += NaNInterface::validParams();
  params.addClassDescription("Two-phase carbon dioxide fluid properties");
  return params;
}

CarbonDioxideTwoPhaseFluidProperties::CarbonDioxideTwoPhaseFluidProperties(
    const InputParameters & parameters)
  : TwoPhaseFluidProperties(parameters), NaNInterface(this)
{
  if (_tid == 0)
  {
    std::string class_name = "CarbonDioxideLiquidFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    _fe_problem.addUserObject(class_name, _liquid_name, params);
  }
  _fp_liquid = &_fe_problem.getUserObjectTempl<SinglePhaseFluidProperties>(_liquid_name, _tid);

  if (_tid == 0)
  {
    std::string class_name = "CarbonDioxideVaporFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    _fe_problem.addUserObject(class_name, _vapor_name, params);
  }
  _fp_vapor = &_fe_problem.getUserObjectTempl<SinglePhaseFluidProperties>(_vapor_name, _tid);
}

Real
CarbonDioxideTwoPhaseFluidProperties::p_critical() const
{
  return _P_critical;
}

Real
CarbonDioxideTwoPhaseFluidProperties::T_sat(Real pressure) const
{
  pressure *= 1.e-6;

  static const double p0 = 0.0005; // extrapolated; triple point would be 0.518 MPa
  static const double pc = 7.37729837321;

  if (p0 <= pressure && pressure <= pc)
    return TS_P_CO2(pressure);
  else
    return getNaN();
}

Real
CarbonDioxideTwoPhaseFluidProperties::p_sat(Real temperature) const
{
  static const double t0 = 216.592;
  static const double tc = 304.1282;

  if (t0 < temperature && temperature < tc)
    return PS_T_INV_CO2(temperature) * 1e6;
  else
    return getNaN();
}

Real
CarbonDioxideTwoPhaseFluidProperties::dT_sat_dp(Real pressure) const
{
  pressure *= 1.e-6;

  double temperature, dtsdp;

  static const double p0 = 0.0005; // extrapolated; triple point would be 0.518 MPa
  static const double pc = 7.37729837321;

  if (p0 <= pressure && pressure <= pc)
  {
    DIFF_TS_P_CO2(pressure, temperature, dtsdp);
    return dtsdp * 1e-6;
  }
  else
    return getNaN();
}
