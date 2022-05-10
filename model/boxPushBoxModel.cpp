//
// Created by david on 29/04/22.
//

#include "boxPushBoxModel.h"

boxModel::boxModel(m_state _desiredState){
    X_desired = _desiredState.replicate(1, 1);

    R.setIdentity();
    for(int i = 0; i < NUM_CTRL; i++){
        R(i, i) = controlCost[i];
    }
    Q.setIdentity();
    for(int i = 0; i < DOF; i++){
        Q(i, i) = stateCosts[i];
    }
}

// Given a set of mujoco data, what is the cost of its state and controls
float boxModel::costFunction(mjData *d){
    float stateCost;
    m_state X_diff;
    m_state X;
    m_ctrl U;

    X = returnState(d);
    U = returnControls(d);

    VectorXd temp(1);

    // actual - desired
    X_diff = X - X_desired;

    temp = (X_diff.transpose() * Q * X_diff) + (U.transpose() * R * U);

    stateCost = temp(0);

    return stateCost;
}

// Given a set of mujoco data, what are its cost derivates with respect to state and control
void boxModel::costDerivatives(mjData *d, Ref<m_state> l_x, Ref<m_state_state> l_xx, Ref<m_ctrl> l_u, Ref<m_ctrl_ctrl> l_uu){
    m_state X_diff;
    m_state X;
    m_ctrl U;

    X = returnState(d);
    U = returnControls(d);

    // actual - desired
    X_diff = X - X_desired;

    l_x = 2 * Q * X_diff;
    l_xx = 2 * Q;

    l_u = 2 * R * U;
    l_uu = 2 * R;
}

// set the state of a mujoco data object as per this model
void boxModel::setState(mjData *d, m_state X){

    d->qpos[0] = X(0);
    d->qpos[1] = X(1);
    d->qpos[7] = X(2);
    d->qpos[8] = X(3);

    d->qvel[0] = X(4);
    d->qvel[1] = X(5);
    d->qvel[6] = X(6);
    d->qvel[7] = X(7);


//    // Two boxes in the state vector
//    for(int i = 0; i < 2; i++){
//        for(int j = 0; j < 2; j++){
//            d->qpos[i] = X(i);
//            d->qvel[i] = X(i+DOF);
//        }
//    }
//
//    for(int i = 0; i < DOF; i++){
//        d->qpos[i] = X(i);
//        d->qvel[i] = X(i+DOF);
//    }
}

// Return the state of a mujoco data model
m_state boxModel::returnState(mjData *d){
    m_state state;

    state(0) = d->qpos[0];
    state(1) = d->qpos[1];
    state(2) = d->qpos[7];
    state(3) = d->qpos[8];
    state(4) = d->qvel[0];
    state(5) = d->qvel[1];
    state(6) = d->qvel[6];
    state(7) = d->qvel[7];


//    for(int i = 0; i < DOF; i++){
//        state(i) = d->qpos[i];
//        state(i + 2) = d->qvel[i];
//    }

    return state;
}

// Set the controls of a mujoco data object
void boxModel::setControls(mjData *d, m_ctrl U){
    for(int i = 0; i < NUM_CTRL; i++){
        d->ctrl[i] = U(i);
    }
}

// Return the controls of a mujoco data object
m_ctrl boxModel::returnControls(mjData *d){
    m_ctrl controls;
    for(int i = 0; i < NUM_CTRL; i++){
        controls(i) = d->ctrl[i];
    }

    return controls;
}

m_dof boxModel::returnAccelerations(mjData *d){
    m_dof accelerations;

    accelerations(0) = d->qacc[0];
    accelerations(1) = d->qacc[1];
    accelerations(2) = d->qacc[6];
    accelerations(3) = d->qacc[7];

//    for(int i = 0; i < DOF; i++){
//        accelerations(i) = d->qacc[i];
//    }

    return accelerations;
}

void boxModel::perturbVelocity(mjData *perturbedData, mjData *origData, int stateIndex, double eps){
    int velIndex = stateIndexToVelIndex[stateIndex];
    perturbedData->qvel[stateIndex] = origData->qvel[stateIndex] + eps;
}

void boxModel::perturbPosition(mjData *perturbedData, mjData *origData, int stateIndex, double eps){
    int posIndex = stateIndexToPosIndex[stateIndex];
    perturbedData->qpos[stateIndex] = origData->qpos[stateIndex] + eps;
}
