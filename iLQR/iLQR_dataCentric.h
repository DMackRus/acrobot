//
// Created by david on 13/04/2022.
//

#ifndef MUJOCO_ACROBOT_CONTROL_ILQR_DATACENTRIC_H
#define MUJOCO_ACROBOT_CONTROL_ILQR_DATACENTRIC_H

#include "mujoco.h"
#include "glfw3.h"
#include "../Utility/stdInclude/stdInclude.h"

#define MUJOCO_DT 0.002
#define ILQR_DT 0.01
#define NUM_MJSTEPS_PER_CONTROL 5
#define MUJOCO_HORIZON_LENGTH 2500
#define ILQR_HORIZON_LENGTH 500

;class iLQR
{
    public:
    // constructor - mujoco model, data, initial controls and initial state
    iLQR(mjModel* m, mjData* d, m_state X0);

    /*      Data     */
    // MuJoCo model and data
    mjModel* model;
    mjData* mdata = NULL;

    // Array of mujoco data structure along the trajectory
    mjData* dArray[ILQR_HORIZON_LENGTH+1];
    // Mujoco data for the initial state of the system
    mjData* d_init;

    /**************************************************************************
     *
     *  iLQR Parameters
     *
     *
     */
    float maxLamda = 10000;             // Maximum lambda before canceliing optimisation
    float minLamda = 0.00001;           // Minimum lamda
    float lamdaFactor = 10;             // Lamda multiplicative factor
    float epsConverge = 0.005;          // Satisfactory convergence of cost function

    float controlCost[NUM_DOF] = {0.01, 0.01};
    float stateCosts[NUM_STATES] = {10, 10, 0.01, 0.01};
    //float terminalScalarConstant = 4;
    //int torqueLims[NUM_DOF] = {100};

    m_dof *finalControls = new m_dof[ILQR_HORIZON_LENGTH];
    m_dof *initControls = new m_dof[ILQR_HORIZON_LENGTH];

    float lamda = 0.1;
    int numIterations = 0;
    m_state X_desired;
    m_dof_dof R;
    m_state_state Q;

    void optimise(m_dof *U, m_state *X);

    void getDerivatives(m_state_state *f_x, m_state_dof *f_u, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu);
    void scaleLinearisation(Ref<m_state_state> A_scaled, Ref<m_state_dof> B_scaled, Ref<m_state_state> A, Ref<m_state_dof> B, int num_steps_per_dt);

    bool backwardsPass_Quu_reg(m_state_state *A, m_state_dof *B, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K);
    void backwardsPass_Vxx_reg(m_state_state *A, m_state_dof *B, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K);

    float forwardsPass(m_dof *k, m_dof_state *K, float oldCost);

    bool checkForConvergence(float newCost, float oldCost);

    void lineariseDynamicsSerial(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum);
    void lineariseDynamicsSerial_trial(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum, float dt);

    float calcStateCost(mjData *d);
    float calcCostDerivatives(Ref<m_state> l_x, Ref<m_state_state> l_xx, Ref<m_dof> l_u, Ref<m_dof_dof> l_uu, int controlNum);

    m_dof returnStateControls(mjData *d);
    m_state returnState(mjData *d);
    void setState(mjData *data, m_state desiredState);
    void setControl(mjData *data, m_dof desiredControl);

    void initCostMatrices();
    void initDesiredState();

};

void cpMjData(const mjModel* m, mjData* d_dest, const mjData* d_src);

#endif //MUJOCO_ACROBOT_CONTROL_ILQR_DATACENTRIC_H
