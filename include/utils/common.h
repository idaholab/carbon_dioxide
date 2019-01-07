#ifndef COMMON_H
#define COMMON_H

#include "FluidProperties.h"

/**
 * Saturation temperature
 * @param pressure pressure
 */
Real T_sat(Real pressure);

/**
 * Saturation pressure
 * @param temperature temperature
 */
Real p_sat(Real temperature);

/**
 * dT/dp along the saturation line
 * @param pressure pressure
 */
Real dT_dP_sat(Real pressure);

/**
 * Surface tension
 * @param temperature temperature
 */
Real sigma(Real temperature);

#endif /* COMMON_H */
