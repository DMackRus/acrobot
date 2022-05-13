//
// Created by david on 03/05/22.
//

#ifndef MUJOCO_ACROBOT_CONTROL_FRANKAARMANDBOX_H
#define MUJOCO_ACROBOT_CONTROL_FRANKAARMANDBOX_H

#include "mujoco.h"
#include "../Utility/stdInclude/stdInclude.h"
#include "../Utility/MujocoController/MujocoController.h"

;class frankaModel{
public:
    frankaModel(mjModel *m, m_state _desiredState);

    float controlCost[NUM_CTRL] = {0, 0, 0, 0, 0, 0, 0};

    // State vector is: 7 joint angles, two cube pos (X and Y), cube rot, 7 joint velocities, two cube velocities (X and Y)
    float stateCosts[(2 * DOF)] = {0, 0, 0, 0, 0, 0, 0,
                                    5, 5, 5,
                                    0, 0, 0, 0, 0,0, 0,
                                    0.1, 0.1};

//    float stateCosts[(2 * DOF)] = {0, 0, 0, 0, 0, 0, 0,
//                                   10, 0.5,
//                                   0, 0, 0, 0, 0,0, 0,
//                                   0.1, 0.1};

    int terminalConstant = 10;

    std::vector<std::string> stateNames;
    int stateIndexToStateName[DOF] = {0, 1, 2, 3, 4, 5, 6, 7, 7, 7};

    mjModel* model;
    m_state X_desired;
    m_ctrl_ctrl R;
    m_state_state Q;

    // Given a set of mujoco data, what is the cost of its state and controls
    float costFunction(mjData *d, int controlNum, int totalControls);

    // Given a set of mujoco data, what are its cost derivates with respect to state and control
    void costDerivatives(mjData *d, Ref<m_state> l_x, Ref<m_state_state> l_xx, Ref<m_ctrl> l_u, Ref<m_ctrl_ctrl> l_uu, int controlNum, int totalControls);

    // set the state of a mujoco data object as per this model
    void setState(mjData *d, m_state X);

    // Return the state of a mujoco data model
    m_state returnState(mjData *d);

    // Set the controls of a mujoco data object
    void setControls(mjData *d, m_ctrl U);

    // Return the controls of a mujoco data object
    m_ctrl returnControls(mjData *d);

    m_dof returnVelocities(mjData *d);
    m_dof returnAccelerations(mjData *d);

    void perturbVelocity(mjData *perturbedData, mjData *origData, int stateIndex, double eps);

    void perturbPosition(mjData *perturbedData, mjData *origData, int stateIndex, double eps);


};

#endif //MUJOCO_ACROBOT_CONTROL_FRANKAARMANDBOX_H
