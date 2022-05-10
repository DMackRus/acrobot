//
// Created by david on 13/04/2022.
//

#ifndef MUJOCO_ACROBOT_CONTROL_ILQR_DATACENTRIC_H
#define MUJOCO_ACROBOT_CONTROL_ILQR_DATACENTRIC_H

#include "mujoco.h"
#include "glfw3.h"
#include "../Utility/stdInclude/stdInclude.h"
#include "../model/acrobotModel.h"
#include "../model/boxPushBoxModel.h"
#include "../model/frankaArmAndBox.h"
#include "../Utility/MujocoController/MujocoUI.h"

#define MUJOCO_DT 0.002
#define ILQR_DT 0.002
#define NUM_MJSTEPS_PER_CONTROL 1
#define ILQR_HORIZON_LENGTH 5000

;
//template<int HORIZON_LENGTH>
class iLQR
{
    public:

    // constructor - mujoco model, data, initial controls and initial state
    iLQR(mjModel* m, mjData* d, m_state _X0, frankaModel* _modelTranslator, MujocoController* _mujocoController);

    /*      Data     */
    // MuJoCo model and data
    mjModel* model;
    mjData* mdata = NULL;
    frankaModel *modelTranslator;
    MujocoController *mujocoController;

    // Array of mujoco data structure along the trajectory
    mjData* dArray[ILQR_HORIZON_LENGTH+1];
    // Mujoco data for the initial state of the system
    mjData* d_init;
    m_state X0;

    /**************************************************************************
     *
     *  iLQR Parameters
     *
     *
     */
    float maxLamda = 10000;             // Maximum lambda before canceliing optimisation
    float minLamda = 0.00001;           // Minimum lamda
    float lamdaFactor = 10;             // Lamda multiplicative factor
    float epsConverge = 0.01;          // Satisfactory convergence of cost function
    int maxIterations = 20;

    //int torqueLims[NUm_ctrl] = {100};

    std::vector<m_ctrl> initControls;
    std::vector<m_ctrl> finalControls;

    // Initialise partial differentiation matrices for all timesteps T
    // for linearised dynamics
    std::vector<m_state_state> f_x;
    std::vector<m_state_ctrl> f_u;

    std::vector<m_state_state> A_scaled;
    std::vector<m_state_ctrl> B_scaled;

    // Quadratic cost partial derivatives
    std::vector<m_state> l_x;
    std::vector<m_state_state> l_xx;
    std::vector<m_ctrl> l_u;
    std::vector<m_ctrl_ctrl> l_uu;

    // Initialise state feedback gain matrices
    std::vector<m_ctrl> k;
    std::vector<m_ctrl_state> K;

    // Initialise new controls and states storage for evaluation
    std::vector<m_ctrl> U_new;
    std::vector<m_state> X_new;

    float lamda = 0.1;
    int numIterations = 0;

    void optimise();

    void getDerivatives();
    void scaleLinearisation(Ref<m_state_state> A_scaled, Ref<m_state_ctrl> B_scaled, Ref<m_state_state> A, Ref<m_state_ctrl> B, int num_steps_per_dt);

    bool backwardsPass_Quu_reg();
    void backwardsPass_Vxx_reg();
    bool isMatrixPD(Ref<MatrixXd> matrix);

    float forwardsPass(float oldCost);

    bool checkForConvergence(float newCost, float oldCost);

    void lineariseDynamicsSerial(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum);
    void lineariseDynamicsSerial_trial(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum, float ilqr_dt);

    void lineariseDynamicsSerial_trial_step(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum, float dt);

    m_ctrl returnDesiredControl(int controlIndex, bool finalControl);
    void setInitControls(std::vector<m_ctrl> _initControls);
    void makeDataForOptimisation();


};

void cpMjData(const mjModel* m, mjData* d_dest, const mjData* d_src);

#endif //MUJOCO_ACROBOT_CONTROL_ILQR_DATACENTRIC_H
