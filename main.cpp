#include "Utility/MujocoController/MujocoUI.h"
#include "iLQR/iLQR_dataCentric.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include <Eigen/Core>
#include "mujoco.h"

#define RUN_ILQR 1
#define TEST_LINEARISATION 0


extern MujocoController *globalMujocoController;
extern mjModel* model;						// MuJoCo model
extern mjData* mdata;						// MuJoCo data

extern iLQR* optimiser;

ofstream outputDiffDyn;
std::string diffDynFilename = "diffDyn.csv";
ofstream outputFile;
std::string filename = "finalTrajectory.csv";

extern m_dof *finalControls;
extern m_dof *initControls;

void saveStates(m_state *X_dyn, m_state *X_lin);
void saveTrajecToCSV(m_dof *U, m_state *X);

void testILQR();

int main() {
    initMujoco();
    //initialseController();


    if(RUN_ILQR){
        m_state X0;
        X0 << 1.5, 0, 0, 0;
        optimiser = new iLQR(model, mdata, X0);
        testILQR();
        render();
    }
    else{
        //simpleTest();
        //render_simpleTest();
    }

    return 0;
}

void testILQR(){
    m_state X0;
    X0 << 1.5, 0, 0, 0;

    outputFile.open(filename);
    outputFile << "Control Number" << ",";

    for(int i = 0; i < NUM_DOF; i++){
        string title = "T" + i;
        outputFile << title << ",";
    }

    for(int i = 0; i < NUM_STATES / 2; i++){
        string title = "P" + i;
        outputFile << title << ",";
    }

    for(int i = 0; i < NUM_STATES / 2; i++){
        string title = "V" + i;
        outputFile << title << ",";
    }
    outputFile << endl;


    auto iLQRStart = high_resolution_clock::now();
    m_state *X_dyn = new m_state[ILQR_HORIZON_LENGTH + 1];


//    for(int i = 0; i < NUM_DOF; i++){
//        mdata->qpos[i] = X0(i);
//        mdata->qvel[i] = X0(i+2);
//    }

//    m_state *X_lin = new m_state[ILQR_HORIZON_LENGTH + 1];

//    X_lin[0] = X0;
//    X_dyn[0] = X0;
//
//    m_state diff;
//    float sumDiff[NUM_STATES] = {0};
//    MatrixXd A_last = ArrayXXd::Zero(NUM_STATES, NUM_STATES);
//    MatrixXd B_last = ArrayXXd::Zero(NUM_STATES, NUM_DOF);

    if(TEST_LINEARISATION){
        for(int i = 0;  i < ILQR_HORIZON_LENGTH; i++){
            //globalMujocoController->setSystemState(X_dyn[i]);
            MatrixXd A = ArrayXXd::Zero(NUM_STATES, NUM_STATES);
            MatrixXd B = ArrayXXd::Zero(NUM_STATES, NUM_DOF);

            //lineariseDynamicsParallel(A, B, i);

            mj_step(model, mdata);

            //X_dyn[i + 1] = returnState(dArray[i]);


            m_state xdot_lin;
            m_state actual_x_dot;

//            m_state_state A_dt;
//            m_state_dof B_dt;
//            int scale = mujoco_steps_per_dt / linearising_num_sim_steps;
            //scaleLinearisation(A, B, A_dt, B_dt, scale);

            //X_lin[i+1] = (A * X_lin[i]) + (B * initControls);
            //xdot_lin = X_lin[i+1] - X_lin[i];
            //actual_x_dot = X_dyn[i+1] - X_dyn[i];

//            cout << "xdot linearising was: " << endl;
//            cout << xdot_lin << endl;
//            cout << "xdot from dynamics was: " << endl;
//            cout << actual_x_dot << endl;
//            cout << "XLin is: " << endl;
//            cout << X_lin[i+1] << endl;
//            cout << "XDyn is: " << endl;
//            cout << X_dyn[i+1] << endl;

//            m_state A_temp = (A * X_dyn[i]);
//            m_state B_temp = (B * initControls[i]);
//            //cout << "A part is: " << endl << A_temp << endl;
//            //cout << "b part is: " << endl << B_temp << endl;
//            cout << "A" << endl << A_dt << endl;
//            cout << "B" << endl << B_dt << endl;

            //diff = (actual_x_dot - xdot_lin);

//            for(int j = 0; j < NUM_STATES; j++){
//                if(diff(j) > 0.1){
//                    cout << "index: " << j << endl;
//                    int a = 1;
//                }
//                if(diff(j) < -0.1){
//                    cout << "index: " << j << endl;
//                    int a = 1;
//                }
//            }

//            for(int  j = 0; j < NUM_STATES; j++){
//                sumDiff[j] += pow((X_dyn[i](j) - X_lin[i](j)), 2);
//
//            }
        }

//        cout << "sum squared diff at end: " << endl;
//        float totalBadness = 0.0f;
//        for(int i = 0; i < NUM_STATES; i++){
//            cout << sumDiff[i] << " "  << endl;
//            totalBadness += sumDiff[i];
//        }
//        cout << "total badness: " << totalBadness << endl;
//        saveStates(X_dyn, X_lin);
    }

    optimiser->optimise(finalControls, X_dyn);

    auto iLQRStop = high_resolution_clock::now();
    auto iLQRDur = duration_cast<microseconds>(iLQRStop - iLQRStart);

    cout << "iLQR took " << iLQRDur.count()/1000000 << " seconds" << endl;

    saveTrajecToCSV(finalControls, X_dyn);

}

void saveTrajecToCSV(m_dof *U, m_state *X){

    for(int i = 0; i < ILQR_HORIZON_LENGTH; i++){
        outputFile << i << ",";

        for(int j = 0; j < NUM_DOF; j++){
            outputFile << U[i](j) << ",";
        }

        for(int j = 0; j < NUM_STATES / 2; j++){
            outputFile << X[i](j) << ",";
        }

        for(int j = 0; j < NUM_STATES / 2; j++){
            outputFile << X[i](j+2) << ",";
        }
        outputFile << endl;
    }

    outputFile.close();
}

void saveStates(m_state *X_dyn, m_state *X_lin){
    ofstream outputDiffDyn;

    std::string diffDynFilename = "diffDyn.csv";
    outputDiffDyn.open(diffDynFilename);
    outputDiffDyn << "X0 dyn" << "," << "X0 lin" << "," << "X0 diff" << "," << "X1 dyn" << "," << "X1 lin" << "," << "X1 diff" << ",";
    outputDiffDyn << "X2 dyn" << "," << "X2 lin" << "," << "X2 diff" << "," << "X3 dyn" << "," << "X3 lin" << "," << "X3 diff" << ",";
    outputDiffDyn << "X4 dyn" << "," << "X4 lin" << "," << "X4 diff" << "," << "X5 dyn" << "," << "X5 lin" << "," << "X5 diff" << ",";
    outputDiffDyn << "X6 dyn" << "," << "X6 lin" << "," << "X6 diff" << "," << "X7 dyn" << "," << "X7 lin" << "," << "X7 diff" << ",";
    outputDiffDyn << "X8 dyn" << "," << "X8 lin" << "," << "X8 diff" << "," << "X9 dyn" << "," << "X9 lin" << "," << "X9 diff" << ",";
    outputDiffDyn << "X10 dyn" << "," << "X10 lin" << "," << "X10 diff" << "," << "X11 dyn" << "," << "X11 lin" << "," << "X11 diff" << ",";
    outputDiffDyn << "X12 dyn" << "," << "X12 lin" << "," << "X12 diff" << "," << "X13 dyn" << "," << "X13 lin" << "," << "X13 diff" << ",";
    outputDiffDyn << "X14 dyn" << "," << "X14 lin" << "," << "X14 diff" << "," << "X15 dyn" << "," << "X15 lin" << "," << "X15 diff" << "," << endl;

    for(int i = 0; i < ILQR_HORIZON_LENGTH; i++){
        for(int j = 0; j < NUM_STATES; j++){
            outputDiffDyn << X_dyn[i](j) << ",";
            outputDiffDyn << X_lin[i](j) << ",";
            outputDiffDyn << X_lin[i](j) - X_dyn[i](j) << ",";
        }
        outputDiffDyn << endl;
    }
}




