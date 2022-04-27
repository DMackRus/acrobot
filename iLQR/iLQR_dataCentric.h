//
// Created by david on 13/04/2022.
//

#ifndef MUJOCO_ACROBOT_CONTROL_ILQR_DATACENTRIC_H
#define MUJOCO_ACROBOT_CONTROL_ILQR_DATACENTRIC_H

#include "../Utility/MujocoController/MujocoController.h"
#include "ilqrCore.h"
#include "iLQR_funcs.h"
#include "../util.h"

#define MUJOCO_DT 0.002
#define ILQR_DT 0.01
#define NUM_MJSTEPS_PER_CONTROL 5
#define MUJOCO_HORIZON_LENGTH 2500
#define ILQR_HORIZON_LENGTH 500
#define TEST_LINEARISATION 0

void iLQR(m_state X0, m_dof *U, m_state *X);
void differentiateDynamics(m_state_state *f_x, m_state_dof *f_u, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu);
void scaleLinearisation(Ref<m_state_state> A_scaled, Ref<m_state_dof> B_scaled, Ref<m_state_state> A, Ref<m_state_dof> B, int num_steps_per_dt);
bool backwardsPass_Quu_reg(m_state_state *A, m_state_dof *B, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K);
void backwardsPass_Vxx_reg(m_state_state *A, m_state_dof *B, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K);
float forwardsPass(m_dof *k, m_dof_state *K, float oldCost);
bool checkForConvergence(float newCost, float oldCost);

void lineariseDynamicsParallel(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum);
void lineariseDynamicsSerial(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum);
void lineariseDynamicsSerial_trial(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum, float ilqr_dt);
void stepSimulation(const Ref<m_state> currentState, const Ref<m_dof> U, Ref<m_state> Xnew, Ref<m_state> Xdot, int numSimSteps);
float rollOutTrajectory(const Ref<const m_state> X0, m_state *X, m_dof *U, int numControls);
float calcStateCost(mjData *d);
float calcCostDerivatives(Ref<m_state> l_x, Ref<m_state_state> l_xx, Ref<m_dof> l_u, Ref<m_dof_dof> l_uu, int controlNum);

m_dof returnStateControls(mjData *d);
m_state returnState(mjData *d);
void setState(mjData *data, m_state desiredState);
void setControl(mjData *data, m_dof desiredControl);

m_dof returnNextControl(int controlNum, bool finalControl);

void initCostMatrices();
void initDesiredState();

void testILQR();
void simpleTest();

void saveStates(m_state *X_dyn, m_state *X_lin);
void saveTrajecToCSV(m_dof *U, m_state *X);

#endif //MUJOCO_ACROBOT_CONTROL_ILQR_DATACENTRIC_H
