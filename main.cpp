#include "Utility/MujocoController/MujocoUI.h"
#include "iLQR/iLQR_dataCentric.h"
#include "model/acrobotModel.h"
//#include "model/boxPushBoxModel.h"
#include "model/frankaArmAndBox.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include <Eigen/Core>
#include "mujoco.h"

#define RUN_ILQR 1
#define TEST_LINEARISATION 1


extern MujocoController *globalMujocoController;
extern mjModel* model;						// MuJoCo model
extern mjData* mdata;						// MuJoCo data

extern iLQR* optimiser;
frankaModel* modelTranslator;

std::vector<m_state> X_dyn;
std::vector<m_state> X_lin;

m_state X0;
m_state X_desired;

ofstream outputDiffDyn;
std::string diffDynFilename = "diffDyn.csv";
ofstream outputFile;
std::string filename = "finalTrajectory.csv";

extern std::vector<m_ctrl> testInitControls;
extern mjData* d_init_test;

void saveStates();
void saveTrajecToCSV();

void testILQR(m_state X0);

void simpleTest();
void initControls();

void setStateTest();

int main() {
    initMujoco();
    // Franka arm with end effector parallel to ground configuration
//    X0 << 0, 0.275, -0.119, -2.76, 2.97, 0, 0,
//            0.7, 0,
//            0, 0, 0, 0, 0, 0, 0,
//            0, 0;

    X0 << -0.12, -0.147, 0.06, -3, 0, 1.34, 0,
            0.4, 0, 0,
            0, 0, 0, 0, 0, 0, 0,
            0, 0, 0;

    X_desired << 0 ,0.88, 0.0297, -1.13, -0.0593, 0, 0,
            0.8, 0.3, 0,
            0, 0, 0, 0, 0, 0, 0,
            0, 0, 0;


    if(RUN_ILQR){

        modelTranslator = new frankaModel(model, X_desired);
        optimiser = new iLQR(model, mdata, X0, modelTranslator, globalMujocoController);
        initControls();
        optimiser->setInitControls(testInitControls);
        optimiser->makeDataForOptimisation();
        testILQR(X0);
        render();
    }
    else{
        modelTranslator = new frankaModel(model, X_desired);
        simpleTest();
        render_simpleTest();
    }

    return 0;
}

void setStateTest(){
    d_init_test = mj_makeData(model);
    modelTranslator = new frankaModel(model, X_desired);
    cout << "X0 is: " << X0 << endl;
    modelTranslator->setState(mdata, X0);

    for(int i = 0; i < 5; i++){
        mj_step(model, mdata);
    }
    cpMjData(model, d_init_test, mdata);
}

void initControls(){
    d_init_test = mj_makeData(model);
    modelTranslator->setState(mdata, X0);
    for(int i = 0; i < 5; i++){
        mj_step(model, mdata);
    }
    cpMjData(model, d_init_test, mdata);

    const std::string test1 = "panda0_link0";
    const std::string test2 = "panda0_rightfinger";
    const std::string test3 = "box_obstacle_1";
    int test1id = mj_name2id(model, mjOBJ_BODY, test1.c_str());
    int test2id = mj_name2id(model, mjOBJ_BODY, test2.c_str());
    int test3id = mj_name2id(model, mjOBJ_BODY, test3.c_str());
    cout << "id 1: " << test1id << "id 2: " << test2id << " id: 3: " << test3id <<  endl;
    cout << "model nv: " << model->nv << " model nq: " << model->nq << " model nu: " << model->nu << endl;

    const std::string endEffecName = "franka_gripper";
    int endEffecId = mj_name2id(model, mjOBJ_BODY, endEffecName.c_str());

    m_pose startPose = globalMujocoController->returnBodyPose(model, mdata, endEffecId);
    m_quat startQuat = globalMujocoController->returnBodyQuat(model, mdata, endEffecId);

    cout << "start pose is: " << startPose << endl;
    m_pose endPose;
    m_pose direction;
    direction << 0.6, 0, 0, 0, 0, 0;
    float magnitudeDiff = sqrt((pow(direction(0), 2)) + (pow(direction(1), 2)) + (pow(direction(2), 2)));
    float forceMagnitude = 30;
    // normalise vector diff
    direction /= magnitudeDiff;

    m_pose linearInterpolationDesiredForce;
    linearInterpolationDesiredForce = direction * forceMagnitude;
    endPose = startPose + direction;

    for(int i = 0; i <= ILQR_HORIZON_LENGTH; i++){

        m_pose currentEEPose = globalMujocoController->returnBodyPose(model, mdata, endEffecId);
        m_quat currentQuat = globalMujocoController->returnBodyQuat(model, mdata, endEffecId);
        m_quat invCurrentQuat = globalMujocoController->invQuat(currentQuat);

        m_quat quatDiff = globalMujocoController->multQuat(startQuat, invCurrentQuat);
        cout << "quat diff: " << quatDiff << endl;
        cout << "end effec pose: " << currentEEPose << endl;
        MatrixXd Jac = globalMujocoController->calculateJacobian(model, mdata, endEffecId);
        MatrixXd Jac_t = Jac.transpose();

        m_point axisDiff = globalMujocoController->quat2Axis(quatDiff);

        m_ctrl desiredControls;
        m_pose desiredEEForce;
        m_pose diff;
        // cout << "currentEEPoint: " << currentEEPoint << endl;
        diff = (currentEEPose - endPose);
        diff(3) = axisDiff(0);
        diff(4) = axisDiff(1);
        diff(5) = axisDiff(2);
        desiredEEForce = linearInterpolationDesiredForce;

        float zAxisRedFactor = 100 * diff(2);
        float rollAxisRedFactor = 10 * diff(3);
        float pitchAxisRedFactor = 10 * diff(4);
        float yawAxisRedFactor = 10 * diff(5);
        desiredEEForce(2) -= zAxisRedFactor;
        desiredEEForce(3) -= rollAxisRedFactor;
        desiredEEForce(4) -= pitchAxisRedFactor;
        desiredEEForce(5) -= yawAxisRedFactor;

        //cout << "Jac_t: " << Jac_t << endl;
        //cout << "desiredEEForce: " << desiredEEForce << endl;

        desiredControls = Jac_t * desiredEEForce;
        //cout << "desired controls: " << desiredControls << endl;

        for(int k = 0; k < NUM_CTRL; k++){

            testInitControls.push_back(m_ctrl());
            testInitControls[i](k) = desiredControls(k) + mdata->qfrc_bias[k];
            mdata->ctrl[k] = testInitControls[i](k);
        }

        for(int j = 0; j < NUM_MJSTEPS_PER_CONTROL; j++){
            mj_step(model, mdata);
        }
    }
}


void simpleTest(){

    initControls();
//    for(int i=0; i < ILQR_HORIZON_LENGTH; i++) {
//        testInitControls.push_back(m_ctrl());
//    }
//
//    d_init_test = mj_makeData(model);
//    modelTranslator = new frankaModel(model, X_desired);
//    cout << "X0 is: " << X0 << endl;
//    modelTranslator->setState(mdata, X0);
//
//    for(int i = 0; i < 5; i++){
//        mj_step(model, mdata);
//    }
//    cpMjData(model, d_init_test, mdata);
//
//    const std::string test1 = "panda0_link0";
//    const std::string test2 = "panda0_rightfinger";
//    const std::string test3 = "box_obstacle_1";
//    int test1id = mj_name2id(model, mjOBJ_BODY, test1.c_str());
//    int test2id = mj_name2id(model, mjOBJ_BODY, test2.c_str());
//    int test3id = mj_name2id(model, mjOBJ_BODY, test3.c_str());
//    cout << "id 1: " << test1id << "id 2: " << test2id << " id: 3: " << test3id <<  endl;
//    cout << "model nv: " << model->nv << " model nq: " << model->nq << " model nu: " << model->nu << endl;
//
//
//    const std::string endEffecName = "franka_gripper";
//    int endEffecId = mj_name2id(model, mjOBJ_BODY, endEffecName.c_str());
//
//    m_pose startPose = globalMujocoController->returnBodyPose(model, mdata, endEffecId);
//    m_quat startQuat = globalMujocoController->returnBodyQuat(model, mdata, endEffecId);
//
//    cout << "start pose is: " << startPose << endl;
//    m_pose endPose;
//    m_pose direction;
//    direction << 0.6, 0, 0, 0, 0, 0;
//    float magnitudeDiff = sqrt((pow(direction(0), 2)) + (pow(direction(1), 2)) + (pow(direction(2), 2)));
//    float forceMagnitude = 10;
//    // normalise vector diff
//    direction /= magnitudeDiff;
//
//    m_pose linearInterpolationDesiredForce;
//    linearInterpolationDesiredForce = direction * forceMagnitude;
//    endPose = startPose + direction;
//
//    for(int i = 0; i <= ILQR_HORIZON_LENGTH; i++){
//
//        m_pose currentEEPose = globalMujocoController->returnBodyPose(model, mdata, endEffecId);
//        m_quat currentQuat = globalMujocoController->returnBodyQuat(model, mdata, endEffecId);
//        m_quat invCurrentQuat = globalMujocoController->invQuat(currentQuat);
//
//        m_quat quatDiff = globalMujocoController->multQuat(startQuat, invCurrentQuat);
//        cout << "quat diff: " << quatDiff << endl;
//        cout << "end effec pose: " << currentEEPose << endl;
//        MatrixXd Jac = globalMujocoController->calculateJacobian(model, mdata, endEffecId);
//        MatrixXd Jac_t = Jac.transpose();
//
//        m_point axisDiff = globalMujocoController->quat2Axis(quatDiff);
//
//        m_ctrl desiredControls;
//        m_pose desiredEEForce;
//        m_pose diff;
//       // cout << "currentEEPoint: " << currentEEPoint << endl;
//        diff = (currentEEPose - endPose);
//        diff(3) = axisDiff(0);
//        diff(4) = axisDiff(1);
//        diff(5) = axisDiff(2);
//        desiredEEForce = linearInterpolationDesiredForce;
//
//        float zAxisRedFactor = 100 * diff(2);
//        float rollAxisRedFactor = 10 * diff(3);
//        float pitchAxisRedFactor = 10 * diff(4);
//        float yawAxisRedFactor = 10 * diff(5);
//        desiredEEForce(2) -= zAxisRedFactor;
//        desiredEEForce(3) -= rollAxisRedFactor;
//        desiredEEForce(4) -= pitchAxisRedFactor;
//        desiredEEForce(5) -= yawAxisRedFactor;
//
//        //cout << "Jac_t: " << Jac_t << endl;
//        //cout << "desiredEEForce: " << desiredEEForce << endl;
//
//        desiredControls = Jac_t * desiredEEForce;
//        //cout << "desired controls: " << desiredControls << endl;
//
//        for(int k = 0; k < NUM_CTRL; k++){
//
//            testInitControls[i](k) = desiredControls(k) + mdata->qfrc_bias[k];
//            mdata->ctrl[k] = testInitControls[i](k);
//        }
//
//        for(int j = 0; j < NUM_MJSTEPS_PER_CONTROL; j++){
//            mj_step(model, mdata);
//        }
//    }

}

void testILQR(m_state X0){

    outputFile.open(filename);
    outputFile << "Control Number" << ",";

    for(int i = 0; i < NUM_CTRL; i++){
        string title = "T" + i;
        outputFile << title << ",";
    }

    for(int i = 0; i < DOF / 2; i++){
        string title = "P" + i;
        outputFile << title << ",";
    }

    for(int i = 0; i < DOF / 2; i++){
        string title = "V" + i;
        outputFile << title << ",";
    }
    outputFile << endl;

    auto iLQRStart = high_resolution_clock::now();

    for(int i = 0; i < ILQR_HORIZON_LENGTH+1; i++){
        X_dyn.push_back(m_state());
        X_lin.push_back(m_state());
    }



    double cubeVelDiff = 0.0f;

    bool firstXDot = true;
    m_state prevXDot;
    m_ctrl prevControl;
    if(TEST_LINEARISATION){
        for(int i = 0;  i < ILQR_HORIZON_LENGTH; i++){
            MatrixXd A = ArrayXXd::Zero((2 * DOF), (2 * DOF));
            MatrixXd B = ArrayXXd::Zero((2 * DOF), NUM_CTRL);
            MatrixXd A_dt = ArrayXXd::Zero((2 * DOF), (2 * DOF));
            MatrixXd B_dt = ArrayXXd::Zero((2 * DOF), NUM_CTRL);

            //optimiser->lineariseDynamicsSerial_trial(A, B, i, ILQR_DT);

            //optimiser->scaleLinearisation(A_dt, B_dt, A, B, 1);

            if(i == 0){
                // set X0 dyn and lin as the initial state of the system
                X_dyn[0] = modelTranslator->returnState(optimiser->dArray[0]);
                X_lin[0] = modelTranslator->returnState(optimiser->dArray[0]);

                m_dof velocities = modelTranslator->returnVelocities(mdata);
                m_dof accelerations = modelTranslator->returnAccelerations(mdata);
                for(int j = 0; j < DOF; j++){
                    prevXDot(j) = velocities(j);
                    prevXDot(j + DOF) = accelerations(j);
                }

                X_dyn[1] = modelTranslator->returnState(optimiser->dArray[1]);
                X_lin[1] = modelTranslator->returnState(optimiser->dArray[1]);
            }
//            else if(i == 1){
//                // calculate X1 dyn via mj_step, set X1_lin to this alsogetDerivatives as we dont have enough information yet for linearisation
//                X_dyn[i] = modelTranslator->returnState(optimiser->dArray[i]);
//                X_lin[i] = modelTranslator->returnState(optimiser->dArray[i]);
//
//                m_ctrl appliedControl = modelTranslator->returnControls(optimiser->dArray[i]);
//                modelTranslator->setControls(mdata, appliedControl);
//
//                for(int k = 0; k < NUM_MJSTEPS_PER_CONTROL; k++){
//                    mj_step(model, mdata);
//                }
//
//                m_dof velocities = modelTranslator->returnVelocities(mdata);
//                m_dof accelerations = modelTranslator->returnAccelerations(mdata);
//                for(int j = 0; j < DOF; j++){
//                    prevXDot(j) = velocities(j);
//                    prevXDot(j + DOF) = accelerations(j);
//                }
//            }
            else{
                // Calculate X bar and U bar at current iteration by comparing current state and control with last state and control
                m_state currentState_dyn, lastState_dyn, X_bar;
                m_ctrl currentControl, lastControl, U_bar;
                m_state X_bar_dot, X_dot;
                X_dyn[i] = modelTranslator->returnState(optimiser->dArray[i]);

                currentState_dyn = X_dyn[i].replicate(1,1);
                lastState_dyn = X_dyn[i - 1].replicate(1,1);
                X_bar = currentState_dyn - lastState_dyn;

                lastControl = modelTranslator->returnControls(optimiser->dArray[i-1]);
                currentControl = modelTranslator->returnControls(optimiser->dArray[i]);
                U_bar = currentControl - lastControl;

                // Calculate A and B matrices by linearising around previous state
                optimiser->lineariseDynamicsSerial_trial_step(A, B, i-1, MUJOCO_DT);

                if(i == 2500){
                    cout << "A is" << endl << A << endl;
                    cout << "B is" << endl << B << endl;
                }

                // Calculate X bar dot via X(.) = Ax + BU
                X_bar_dot = (A * X_bar) + (B * U_bar);

                // Calculate current X dot via X(.)(t) = X(.)(t-1) + X_bar_dot
                X_dot = prevXDot + X_bar_dot;

                // Calculate next state via linearisation using euler integration of current X_dot
                if(1){
                    X_lin[i + 1] = X_dyn[i] + (X_dot * MUJOCO_DT);
                }
                else{
                    X_lin[i + 1] = X_lin[i] + (X_dot * MUJOCO_DT);
                }

                // update previous x_dot
                if(1){
                    // either use previous x dot as mujoco x dot
                    m_dof velocities = modelTranslator->returnVelocities(optimiser->dArray[i]);
                    m_dof accelerations = modelTranslator->returnAccelerations(optimiser->dArray[i]);
                    for(int j = 0; j < DOF; j++){
                        prevXDot(j) = velocities(j);
                        prevXDot(j + DOF) = accelerations(j);
                    }
                }
                else{
                    // use linealry calculated  dot (NOTE: this will accumulate error but may be a more valid test)
                    prevXDot = X_dot.replicate(1,1);
                }

                cubeVelDiff += pow(((X_lin[i + 1](17) - X_dyn[i + 1](17)) * ILQR_DT),2);

            }

            if(firstXDot){
                firstXDot = false;
                m_dof velocities = modelTranslator->returnVelocities(optimiser->dArray[i]);
                m_dof accelerations = modelTranslator->returnAccelerations(optimiser->dArray[i]);
                for(int j = 0; j < DOF; j++){
                    prevXDot(j) = velocities(j);
                    prevXDot(j + DOF) = accelerations(j);
                }
            }
        }
        saveStates();
    }
    cout << "cube vel diff total: " << cubeVelDiff << endl;

    // Reset data to initial state
    cpMjData(model, mdata, optimiser->d_init);

    optimiser->optimise();

    auto iLQRStop = high_resolution_clock::now();
    auto iLQRDur = duration_cast<microseconds>(iLQRStop - iLQRStart);

    cout << "iLQR took " << iLQRDur.count()/1000000 << " seconds" << endl;

    saveTrajecToCSV();

}

void saveTrajecToCSV(){

    for(int i = 0; i < ILQR_HORIZON_LENGTH; i++){
        outputFile << i << ",";

        for(int j = 0; j < NUM_CTRL; j++){
            outputFile << optimiser->finalControls[i](j) << ",";
        }

        for(int j = 0; j < DOF; j++){
            outputFile << optimiser->X_new[i](j) << ",";
        }

        for(int j = 0; j < DOF; j++){
            outputFile << optimiser->X_new[i](j+DOF) << ",";
        }
        outputFile << endl;
    }

    outputFile.close();
}

void saveStates(){

    cout << "X_dyn[end]: " << X_dyn[ILQR_HORIZON_LENGTH] << endl;
    cout << "X_lin[end]: " << X_lin[ILQR_HORIZON_LENGTH] << endl;

    cout << "X_dyn[0]: " << X_dyn[0] << endl;
    cout << "X_lin[0]: " << X_lin[0] << endl;


    outputDiffDyn.open(diffDynFilename);
    outputDiffDyn << "Joint 0 dyn" << "," << "Joint 0 lin" << "," << "Joint 0 diff" << "," << "Joint 1 dyn" << "," << "Joint 1 lin" << "," << "Joint 1 diff" << ",";
    outputDiffDyn << "Joint 2 dyn" << "," << "Joint 2 lin" << "," << "Joint 2 diff" << "," << "Joint 3 dyn" << "," << "Joint 3 lin" << "," << "Joint 3 diff" << ",";
    outputDiffDyn << "Joint 4 dyn" << "," << "Joint 4 lin" << "," << "Joint 4 diff" << "," << "Joint 5 dyn" << "," << "Joint 5 lin" << "," << "Joint 5 diff" << ",";
    outputDiffDyn << "Joint 6 dyn" << "," << "Joint 6 lin" << "," << "Joint 6 diff" << ",";
    outputDiffDyn << "Cube X dyn" << "," << "Cube X lin" << "," << "Cube X diff" << "," << "Cube Y dyn" << "," << "Cube Y lin" << "," << "Cube Y diff" << "," << "Cube rot dyn" << "," << "Cube rot lin" << "," << "Cube rot diff" << ",";
    outputDiffDyn << "Joint 0 vel dyn" << "," << "Joint 0 vel lin" << "," << "Joint 0 vel diff" << ",";
    outputDiffDyn << "Joint 1 vel dyn" << "," << "Joint 1 vel lin" << "," << "Joint 1 vel diff" << "," << "Joint 2 vel dyn" << "," << "Joint 2 vel lin" << "," << "Joint 2 vel diff" << ",";
    outputDiffDyn << "Joint 3 vel dyn" << "," << "Joint 3 vel lin" << "," << "Joint 3 vel diff" << "," << "Joint 4 vel dyn" << "," << "Joint 4 vel lin" << "," << "Joint 4 vel diff" << ",";
    outputDiffDyn << "Joint 5 vel dyn" << "," << "Joint 5 vel lin" << "," << "Joint 5 vel diff" << "," << "Joint 6 vel dyn" << "," << "Joint 6 vel lin" << "," << "Joint 6 vel diff" << ",";
    outputDiffDyn << "Cube X  vel dyn" << "," << "Cube X vel lin" << "," << "Cube X vel diff" << "," << "Cube Y vel dyn" << "," << "Cube Y vel lin" << "," << "Cube Y vel diff" << "," << "Cube rot vel dyn" << "," << "Cube rot vel lin" << "," << "Cube rot vel diff" << endl;

    for(int i = 0; i < ILQR_HORIZON_LENGTH; i++){
        for(int j = 0; j < (2 * DOF); j++){
            float val;
            val = X_dyn[i](j);
            outputDiffDyn << val << ",";
            val = X_lin[i](j);
            outputDiffDyn << val << ",";
            val = X_lin[i](j) - X_dyn[i](j);
            outputDiffDyn << val << ",";
        }
        outputDiffDyn << endl;
    }

    outputDiffDyn.close();
}
