//
// Created by davem on 16/02/2022.
//

#ifndef MUJOCO_SANDBOX_ILQRCORE_H
#define MUJOCO_SANDBOX_ILQRCORE_H

#include "../Utility/MujocoController/MujocoController.h"
#include "iLQR_funcs.h"

void differentiateDynamics(m_state *X, m_dof *U, m_state_state *f_x, m_state_dof *f_u, float *l, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu);
bool backwardsPass(m_state_state *f_x, m_state_dof *f_u, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K);
bool backwardsPassTest(m_state_state *f_x, m_state_dof *f_u, float l, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K, m_state *X);
float forwardsPass(m_state *X, m_state *X_new, m_dof *U, m_dof *U_new, m_dof *k, m_dof_state *K);
float forwardsPassTest(m_state *X, m_state *X_best, m_dof *U, m_dof *U_best, m_dof *k, m_dof_state *K);
bool checkForConvergence(float newCost, m_state *X, m_state *X_new, m_dof *U, m_dof *U_new, bool *costImprovement);
void increaseLamda();
bool checkValidAandBMatrices(m_state_state A, m_state_dof B);

#endif //MUJOCO_SANDBOX_ILQRCORE_H
