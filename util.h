//
// Created by david on 05/04/2022.
//

#ifndef MUJOCO_ACROBOT_CONTROL_UTIL_H
#define MUJOCO_ACROBOT_CONTROL_UTIL_H

#include "mujoco.h"

void cpMjData(const mjModel* m, mjData* d_dest, const mjData* d_src);

#endif //MUJOCO_ACROBOT_CONTROL_UTIL_H
