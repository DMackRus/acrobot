//
// Created by david on 03/05/22.
//

#include "frankaArmAndBox.h"

extern MujocoController *globalMujocoController;

frankaModel::frankaModel(mjModel *m, m_state _desiredState){
    X_desired = _desiredState.replicate(1, 1);

    R.setIdentity();
    for(int i = 0; i < NUM_CTRL; i++){
        R(i, i) = controlCost[i];
    }
    Q.setIdentity();
    for(int i = 0; i < DOF; i++){
        Q(i, i) = stateCosts[i];
    }

    model = m;

    // Magic number, 7 joint names and also the cube name
    for(int i = 0; i < 8; i++){
        stateNames.push_back(std::string());

    }
    stateNames[0] = "panda0_link1";
    stateNames[1] = "panda0_link2";
    stateNames[2] = "panda0_link3";
    stateNames[3] = "panda0_link4";
    stateNames[4] = "panda0_link5";
    stateNames[5] = "panda0_link6";
    stateNames[6] = "panda0_link7";
    stateNames[7] = "box_obstacle_1";

}

// Given a set of mujoco data, what is the cost of its state and controls
float frankaModel::costFunction(mjData *d, int controlNum, int totalControls){
    float stateCost;
    m_state X_diff;
    m_state X;
    m_ctrl U;

    X = returnState(d);
    U = returnControls(d);

    VectorXd temp(1);

    // actual - desired
    X_diff = X - X_desired;
    float percentageDone = (float)controlNum / (float)totalControls;
    float terminalScalar = (percentageDone * terminalConstant) + 1;

    temp = (terminalScalar * (X_diff.transpose() * Q * X_diff)) + (U.transpose() * R * U);

    stateCost = temp(0);

    return stateCost;
}

// Given a set of mujoco data, what are its cost derivates with respect to state and control
void frankaModel::costDerivatives(mjData *d, Ref<m_state> l_x, Ref<m_state_state> l_xx, Ref<m_ctrl> l_u, Ref<m_ctrl_ctrl> l_uu, int controlNum, int totalControls){
    m_state X_diff;
    m_state X;
    m_ctrl U;

    X = returnState(d);
    U = returnControls(d);

    // actual - desired
    X_diff = X - X_desired;
    float percentageDone = (float)controlNum / (float)totalControls;
    float terminalScalar = (percentageDone * terminalConstant) + 1;

    l_x = 2 * terminalScalar *  Q * X_diff;
    l_xx = 2 *  terminalScalar * Q;

    l_u = 2 * R * U;
    l_uu = 2 * R;
}

// set the state of a mujoco data object as per this model
void frankaModel::setState(mjData *d, m_state X){

    // Firstly set all the required franka panda arm joints
    for(int i = 0; i < 7; i++){
        int bodyId = mj_name2id(model, mjOBJ_BODY, stateNames[i].c_str());
        globalMujocoController->set_qPosVal(model, d, bodyId, false, 0, X(i));
        globalMujocoController->set_qVelVal(model, d, bodyId, false, 0, X(i + 10));
    }

    // Set the positions of the cube
    int boxId = mj_name2id(model, mjOBJ_BODY, stateNames[7].c_str());
    globalMujocoController->set_qPosVal(model, d, boxId, true, 0, X(7));
    globalMujocoController->set_qPosVal(model, d, boxId, true, 1, X(8));
    globalMujocoController->set_qVelVal(model, d, boxId, true, 0, X(17));
    globalMujocoController->set_qVelVal(model, d, boxId, true, 1, X(18));
    globalMujocoController->set_qVelVal(model, d, boxId, true, 4, X(19));

    m_quat boxQuat = globalMujocoController->returnBodyQuat(model, d, boxId);
    //cout << "returned box quat " << boxQuat << endl;
    m_point axisAngle = globalMujocoController->quat2Axis(boxQuat);
    //cout << "returned box axisAngle " << axisAngle << endl;
    axisAngle(1) = X(9);
    //cout << "new box axisAngle " << axisAngle << endl;
    boxQuat = globalMujocoController->axis2Quat(axisAngle);
    //cout << "new box quat " << boxQuat << endl;
    globalMujocoController->setBodyQuat(model, d, boxId, boxQuat);

}

// Return the state of a mujoco data model
m_state frankaModel::returnState(mjData *d){
    m_state state;

    // Firstly set all the required franka panda arm joints
    for(int i = 0; i < 7; i++){
        int bodyId = mj_name2id(model, mjOBJ_BODY, stateNames[i].c_str());
        state(i) = globalMujocoController->return_qPosVal(model, d, bodyId, false, 0);
        state(i + 10) = globalMujocoController->return_qVelVal(model, d, bodyId, false, 0);
    }

    // Set the positions of the cube
    int boxId = mj_name2id(model, mjOBJ_BODY, stateNames[7].c_str());
    state(7) = globalMujocoController->return_qPosVal(model, d, boxId, true, 0);
    state(8) = globalMujocoController->return_qPosVal(model, d, boxId, true, 1);
    state(17) = globalMujocoController->return_qVelVal(model, d, boxId, true, 0);
    state(18) = globalMujocoController->return_qVelVal(model, d, boxId, true, 1);
    state(19) = globalMujocoController->return_qVelVal(model, d, boxId, true, 4);

    m_quat boxQuat = globalMujocoController->returnBodyQuat(model, d, boxId);
    //cout << "box qaut returned: " << boxQuat << endl;
    m_point axisAngle = globalMujocoController->quat2Axis(boxQuat);
    //cout << "box axis returned: " << boxQuat << endl;
    state(9) = axisAngle(1);


    return state;
}

// Set the controls of a mujoco data object
void frankaModel::setControls(mjData *d, m_ctrl U){
    for(int i = 0; i < NUM_CTRL; i++){
        d->ctrl[i] = U(i);
    }
}

// Return the controls of a mujoco data object
m_ctrl frankaModel::returnControls(mjData *d){
    m_ctrl controls;
    for(int i = 0; i < NUM_CTRL; i++){
        controls(i) = d->ctrl[i];
    }

    return controls;
}

m_dof frankaModel::returnVelocities(mjData *d){
    m_dof velocities;

    for(int i = 0; i < NUM_CTRL; i++){
        int bodyId = mj_name2id(model, mjOBJ_BODY, stateNames[i].c_str());
        velocities(i) = globalMujocoController->return_qVelVal(model, d, bodyId, false, 0);
    }

    int boxId = mj_name2id(model, mjOBJ_BODY, stateNames[7].c_str());
    velocities(NUM_CTRL) = globalMujocoController->return_qVelVal(model, d, boxId, true, 0);
    velocities(NUM_CTRL + 1) = globalMujocoController->return_qVelVal(model, d, boxId, true, 1);
    velocities(NUM_CTRL + 2) = globalMujocoController->return_qVelVal(model, d, boxId, true, 4);

    return velocities;
}

m_dof frankaModel::returnAccelerations(mjData *d){
    m_dof accelerations;

    for(int i = 0; i < NUM_CTRL; i++){
        int bodyId = mj_name2id(model, mjOBJ_BODY, stateNames[i].c_str());
        accelerations(i) = globalMujocoController->return_qAccVal(model, d, bodyId, false, 0);
    }

    int boxId = mj_name2id(model, mjOBJ_BODY, stateNames[7].c_str());
    accelerations(NUM_CTRL) = globalMujocoController->return_qAccVal(model, d, boxId, true, 0);
    accelerations(NUM_CTRL + 1) = globalMujocoController->return_qAccVal(model, d, boxId, true, 1);
    accelerations(NUM_CTRL + 2) = globalMujocoController->return_qAccVal(model, d, boxId, true, 4);


    return accelerations;
}

void frankaModel::perturbVelocity(mjData *perturbedData, mjData *origData, int stateIndex, double eps){
    int stateNameIndex = stateIndexToStateName[stateIndex];
    int bodyId = mj_name2id(model, mjOBJ_BODY, stateNames[stateNameIndex].c_str());
    bool freeJoint;
    int freeJntIndex = 0;

    if(stateIndex <= 6){
        freeJoint = false;
    }
    else{
        freeJoint = true;
        if(stateIndex == 7){
            freeJntIndex = 0;
        }
        else if(stateIndex == 8){
            freeJntIndex = 1;
        }
        else{
            freeJntIndex = 4;
        }
    }

    double origVelocity = globalMujocoController->return_qVelVal(model, origData, bodyId, freeJoint, freeJntIndex);
    double perturbedVel = origVelocity + eps;
    globalMujocoController->set_qVelVal(model, perturbedData, bodyId, freeJoint, freeJntIndex, perturbedVel);

}

void frankaModel::perturbPosition(mjData *perturbedData, mjData *origData, int stateIndex, double eps){
//    int posIndex = stateIndexToPosIndex[stateIndex];
//
//    if(stateIndex == DOF){
//        perturbedData->qpos[posIndex] = origData->qpos[posIndex] + eps;
//    }
//    else{
//        const std::string boxstring = "box_obstacle_1";
//        int boxId = mj_name2id(model, mjOBJ_BODY, boxstring.c_str());
//        int jointIndex = model->body_jntadr[boxId];
//
//        m_quat cubeQuat = globalMujocoController->returnBodyQuat(model, origData, boxId);
//        m_point axisAngle = globalMujocoController->quat2Axis(cubeQuat);
//
//        axisAngle(1) += eps;
//        cubeQuat = globalMujocoController->axis2Quat(axisAngle);
//        globalMujocoController->setBodyQuat(model, origData, boxId, cubeQuat);
//
//    }
}

