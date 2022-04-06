//
// Created by david on 05/04/2022.
//

#ifndef MUJOCO_ACROBOT_CONTROL_UPDATE_H
#define MUJOCO_ACROBOT_CONTROL_UPDATE_H

#include "mujoco.h"


void forwardStep(mjModel* model, mjData* data);

void forwardFrame(mjModel* model, mjData* data);

#endif //MUJOCO_ACROBOT_CONTROL_UPDATE_H
