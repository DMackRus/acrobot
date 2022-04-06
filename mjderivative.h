//
// Created by david on 05/04/2022.
//

#ifndef MUJOCO_ACROBOT_CONTROL_MJDERIVATIVE_H
#define MUJOCO_ACROBOT_CONTROL_MJDERIVATIVE_H

#include "mujoco.h"

typedef mjtNum (*stepCostFn_t)(const mjData*);

void calcMJDerivatives(mjModel* m, mjData* dmain, mjtNum* deriv, stepCostFn_t stepCostFn);

#endif //MUJOCO_ACROBOT_CONTROL_MJDERIVATIVE_H
