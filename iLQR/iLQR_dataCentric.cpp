//
// Created by david on 13/04/2022.
//

#include "iLQR_dataCentric.h"

/**************************************************************************
 *
 *  iLQR Parameters
 *
 *
 */

int linearising_num_sim_steps = 1;  // How many mujoco steps to use for linearising
bool alphaSearchEnabled = true;     // Whether alpha search is enabled to maximise optimisation at each forwards pass
float maxLamda = 10000;             // Maximum lambda before canceliing optimisation
float minLamda = 0.00001;            // Minimum lamda
float lamdaFactor = 10;             // Lamda multiplicative factor
float epsConverge = 0.005;


//float controlCost[DOF] = {0.0001, 0.0001, 0.0001, 0.0001, 0.00005, 0.00005, 0.00005};
float controlCost[NUM_DOF] = {0.01, 0.01};
float stateCosts[NUM_STATES] = {10, 1, 0.1, 0.1};
float terminalScalarConstant = 4;
int torqueLims[NUM_DOF] = {100};

extern MujocoController *globalMujocoController;
extern mjModel* model;						// MuJoCo model
extern mjData* mdata;						// MuJoCo data
extern mjvCamera cam;                   // abstract camera
extern mjvScene scn;                    // abstract scene
extern mjvOption opt;			        // visualization options
extern mjrContext con;				    // custom GPU context
extern GLFWwindow *window;

mjData* dArray[ILQR_HORIZON_LENGTH+1];
mjData* d_init;

m_state X_desired;
m_dof_dof R;
m_state_state Q;

float lamb = 0.1;
int numIterations = 0;

ofstream outputDiffDyn;

std::string diffDynFilename = "diffDyn.csv";

ofstream outputFile;

std::string filename = "finalTrajectory.csv";

m_dof *finalControls = new m_dof[ILQR_HORIZON_LENGTH];
m_dof *initControls = new m_dof[ILQR_HORIZON_LENGTH];

inline mjtNum stepCost(const mjData* d)
{
    mjtNum cost = d->qpos[0];
    return cost;
}

void simpleTest(){
    m_state X0;
    X0 << 2.8, 0, 0, 0;

    for(int i = 0; i < NUM_DOF; i++){
        mdata->qpos[i] = X0(i);
        mdata->qvel[i] = X0(i+2);
    }

    cout << "number of controls: " << model->nu << endl;

    mdata->ctrl[1] = 10;

//    m_state_state *f_x = new m_state_state[numControls];
//    m_state_dof *f_u = new m_state_dof[numControls];
//
//    float l[numControls];
//    m_state *l_x = new m_state[numControls + 1];
//    m_state_state *l_xx = new m_state_state[numControls + 1];
//    m_dof *l_u = new m_dof[numControls];
//    m_dof_dof *l_uu = new m_dof_dof[numControls];
//
//    // Initialise state feedback gain matrices
//    m_dof *k = new m_dof[numControls];
//    m_dof_state *K = new m_dof_state[numControls];
//
//    l_xx[numControls] << 0.02, 0, 0, 0,
//                         0, 0.02, 0, 0,
//                         0, 0, 0.01, 0,
//                         0, 0, 0, 0.01;
//
//    int bPStart = 2;
//
//    l_xx[numControls - bPStart] << 0.02, 0, 0, 0,
//            0, 0.02, 0, 0,
//            0, 0, 0.01, 0,
//            0, 0, 0, 0.01;
//
//    l_uu[numControls] << 0.002;
//    l_uu[numControls - bPStart] << 0.002;
//
//    f_x[numControls - bPStart] << 0.891, 0.0561, 0.0967, 0.00192,
//                                  0.154, 0.846, 0.00475, 0.0950,
//                                  -1.92, 0.980, 0.910, 0.0498,
//                                  2.68, -2.716, 0.127, 0.868;
//
//    f_u[numControls - bPStart] << 0.540, 0.014, 9.67, 0.474;
//
//    l_x[numControls] << 0.00437, 0.00916, -0.00239, 0.00517;
//
//    l_x[numControls - bPStart] << 0.0051, 0.0056, -0.0006, 0.012;
//
//    l_u[numControls - bPStart] << 0.00000127;
//
//    // ---------------------------------------------------------------- //
//
//    l_xx[numControls - bPStart - 1] << 0.02, 0, 0, 0,
//            0, 0.02, 0, 0,
//            0, 0, 0.01, 0,
//            0, 0, 0, 0.01;
//
//    l_uu[numControls - bPStart - 1] << 0.002;
//
//    f_x[numControls - bPStart - 1] << 0.895, 0.0542, 0.0966, 0.000917,
//                                    0.146, 0.849, 0.00518, 0.0964,
//                                    -1.87, 0.955, 0.910, 0.0316,
//                                    2.56, -2.669, 0.133, 0.894;
//
//    f_u[numControls - bPStart - 1] << 0.540, 0.0156, 9.66, 0.505;
//
//    l_x[numControls - bPStart - 1] << 0.00498, 0.00318, 0.0026, 0.0110;
//
//    l_u[numControls - bPStart - 1] << -0.00000522;
//
//    // ------------------------------------------------------------------ //
//
//    backwardsPass_Quu_reg(f_x, f_u, l_x, l_xx, l_u, l_uu, k, K);
}

void testILQR(){
    m_state X0;
    X0 << 1.5, 0, 0, 0;

    auto iLQRStart = high_resolution_clock::now();

    initCostMatrices();
    initDesiredState();

    for(int i = 0; i < NUM_DOF; i++){
        mdata->qpos[i] = X0(i);
        mdata->qvel[i] = X0(i+2);
    }

    m_state *X_lin = new m_state[ILQR_HORIZON_LENGTH + 1];
    m_state *X_dyn = new m_state[ILQR_HORIZON_LENGTH + 1];
    X_lin[0] = X0;
    X_dyn[0] = X0;

    m_state diff;
    float sumDiff[NUM_STATES] = {0};
    MatrixXd A_last = ArrayXXd::Zero(NUM_STATES, NUM_STATES);
    MatrixXd B_last = ArrayXXd::Zero(NUM_STATES, NUM_DOF);

    d_init = mj_makeData(model);
    cpMjData(model, d_init, mdata);
    dArray[0] = mj_makeData(model);
    cpMjData(model, dArray[0], mdata);

    for(int i = 0; i <= ILQR_HORIZON_LENGTH; i++){
        dArray[i] = mj_makeData(model);
        float vel = mdata->qvel[0];
        float kp = 0.1;
        for(int k = 0; k < NUM_DOF; k++){
            //initControls[i](k) = mdata->qfrc_bias[k];
            initControls[i](k) = -kp * vel;
            mdata->ctrl[k] = initControls[i](k);
        }
        cpMjData(model, dArray[i], mdata);

        for(int i = 0; i < NUM_MJSTEPS_PER_CONTROL; i++){
            mj_step(model, mdata);
        }

        // populate dArray with mujoco data objects from start of trajec until end

        // copy current data into current data array place

        // step simulation with initialised controls
    }
    cpMjData(model, dArray[ILQR_HORIZON_LENGTH], mdata);

    cpMjData(model, mdata, d_init);

    if(TEST_LINEARISATION){
        for(int i = 0;  i < ILQR_HORIZON_LENGTH; i++){
            //globalMujocoController->setSystemState(X_dyn[i]);
            MatrixXd A = ArrayXXd::Zero(NUM_STATES, NUM_STATES);
            MatrixXd B = ArrayXXd::Zero(NUM_STATES, NUM_DOF);

            lineariseDynamicsParallel(A, B, i);

            mj_step(model, mdata);

            X_dyn[i + 1] = returnState(dArray[i]);


            m_state xdot_lin;
            m_state actual_x_dot;

//            m_state_state A_dt;
//            m_state_dof B_dt;
//            int scale = mujoco_steps_per_dt / linearising_num_sim_steps;
            //scaleLinearisation(A, B, A_dt, B_dt, scale);

            //X_lin[i+1] = (A * X_lin[i]) + (B * initControls);
            xdot_lin = X_lin[i+1] - X_lin[i];
            actual_x_dot = X_dyn[i+1] - X_dyn[i];

//            cout << "xdot linearising was: " << endl;
//            cout << xdot_lin << endl;
//            cout << "xdot from dynamics was: " << endl;
//            cout << actual_x_dot << endl;
            cout << "XLin is: " << endl;
            cout << X_lin[i+1] << endl;
            cout << "XDyn is: " << endl;
            cout << X_dyn[i+1] << endl;

//            m_state A_temp = (A * X_dyn[i]);
//            m_state B_temp = (B * initControls[i]);
//            //cout << "A part is: " << endl << A_temp << endl;
//            //cout << "b part is: " << endl << B_temp << endl;
//            cout << "A" << endl << A_dt << endl;
//            cout << "B" << endl << B_dt << endl;

            diff = (actual_x_dot - xdot_lin);

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

        cout << "sum squared diff at end: " << endl;
        float totalBadness = 0.0f;
        for(int i = 0; i < NUM_STATES; i++){
            cout << sumDiff[i] << " "  << endl;
            totalBadness += sumDiff[i];
        }
        cout << "total badness: " << totalBadness << endl;
        saveStates(X_dyn, X_lin);
    }

    iLQR(X0, finalControls, X_dyn);

    auto iLQRStop = high_resolution_clock::now();
    auto iLQRDur = duration_cast<microseconds>(iLQRStop - iLQRStart);

    cout << "iLQR took " << iLQRDur.count()/1000000 << " seconds" << endl;

    cout << "final control: " << finalControls[ILQR_HORIZON_LENGTH - 1] << endl;

    saveTrajecToCSV(finalControls, X_dyn);

}

void iLQR(m_state X0, m_dof *U, m_state *X){
    bool optimisationFinished = false;
    int numIterations = 0;
    float newCost = 0;
    float oldCost = 1000;

    // Initialise partial differentiation matrices for all timesteps T
    // for linearised dynamics
    m_state_state *f_x = new m_state_state[ILQR_HORIZON_LENGTH];
    m_state_dof *f_u = new m_state_dof[ILQR_HORIZON_LENGTH];

    // Quadratic cost partial derivatives
    m_state *l_x = new m_state[ILQR_HORIZON_LENGTH + 1];
    m_state_state *l_xx = new m_state_state[ILQR_HORIZON_LENGTH + 1];
    m_dof *l_u = new m_dof[ILQR_HORIZON_LENGTH];
    m_dof_dof *l_uu = new m_dof_dof[ILQR_HORIZON_LENGTH];

    // Initialise state feedback gain matrices
    m_dof *k = new m_dof[ILQR_HORIZON_LENGTH];
    m_dof_state *K = new m_dof_state[ILQR_HORIZON_LENGTH];

    // Initialise new controls and states storage for evaluation
    m_dof *U_new = new m_dof[ILQR_HORIZON_LENGTH];
    m_state *X_new = new m_state[ILQR_HORIZON_LENGTH + 1];

    //Perform a forward pass with initialised controls - FILL
    //oldCost = rollOutTrajectory(X0, X, U, numControls);
    //cout << "INITIAL COST FROM WARM START TRAJECTORY IS: " << oldCost << endl;


    // iterate until optimisation finished
    while(!optimisationFinished){
        numIterations++;


        auto start = high_resolution_clock::now();

        // Linearise the dynamics and save cost values at each state
        // STEP 1 - Linearise dynamics and calculate cost quadratics at every time step
        differentiateDynamics(f_x, f_u, l_x, l_xx, l_u, l_uu);

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Linearising model: " << duration.count()/1000 << " milliseconds" << endl;

        bool validBackPass = false;
        bool lamdaExit = false;

        // STEP 2 - Backwards pass to compute optimal linear and feedback gain matrices k and K
        // Keep running this step until the backwards pass does NOT encounter a non_PD Q_uu_reg matrix
        while(!validBackPass) {

            validBackPass = backwardsPass_Quu_reg(f_x, f_u, l_x, l_xx, l_u, l_uu, k, K);

            if (!validBackPass) {
                if (lamb < maxLamda) {
                    lamb *= lamdaFactor;
                } else {
                    lamdaExit = true;
                    optimisationFinished = true;
                    break;
                }
            } else {
                if (lamb > minLamda) {
                    lamb /= lamdaFactor;
                }
            }
        }

        if(!lamdaExit){
            // STEP 3 - Forwards pass to calculate new optimal controls - with optional alpha binary search
            newCost = forwardsPass(k, K, oldCost);

            // STEP 4 - Check for convergence
            optimisationFinished = checkForConvergence(newCost, oldCost);

            oldCost = newCost;
        }
    }

    for(int i = 0; i < ILQR_HORIZON_LENGTH; i++){
        U[i] = returnStateControls(dArray[i]);
        X[i] = returnState(dArray[i]);
    }
}

void differentiateDynamics(m_state_state *f_x, m_state_dof *f_u, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu){

    MatrixXd A = ArrayXXd::Zero(NUM_STATES, NUM_STATES);
    MatrixXd B = ArrayXXd::Zero(NUM_STATES, NUM_DOF);

    // Linearise the dynamics along the trajectory
    cout << "------------------------ LINEARISE DYNAMICS -----------------------------" << endl;

    for(int t = 0; t < ILQR_HORIZON_LENGTH; t++){
        // Calculate linearised dynamics for current time step via finite differencing
        lineariseDynamicsSerial(A, B, t);

//        m_state_state A_dt;
//        m_state_dof B_dt;
//        int scaling = mujoco_steps_per_dt / linearising_num_sim_steps;
//        scaleLinearisation(A, B, A_dt, B_dt,  scaling);

        f_x[t] = A;
        f_u[t] = B;

        calcCostDerivatives(l_x[t], l_xx[t], l_u[t], l_uu[t], t);

        l_x[t]  *= ILQR_DT;
        l_xx[t] *= ILQR_DT;
        l_u[t]  *= ILQR_DT;
        l_uu[t] *= ILQR_DT;

    }

    calcCostDerivatives(l_x[ILQR_HORIZON_LENGTH], l_xx[ILQR_HORIZON_LENGTH], l_u[ILQR_HORIZON_LENGTH-1], l_uu[ILQR_HORIZON_LENGTH-1], ILQR_HORIZON_LENGTH);
    l_x [ILQR_HORIZON_LENGTH] *= ILQR_DT;
    l_xx[ILQR_HORIZON_LENGTH] *= ILQR_DT;

    cout << "---------------------- END LINEARISE DYNAMICS -----------------------------" << endl;
}

bool backwardsPass_Quu_reg(m_state_state *A, m_state_dof *B, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K){
    m_state V_x;
    V_x = l_x[ILQR_HORIZON_LENGTH];
    m_state_state V_xx;
    V_xx = l_xx[ILQR_HORIZON_LENGTH];

    for(int t = ILQR_HORIZON_LENGTH - 2; t > -1; t--){
        m_state Q_x;
        m_dof Q_u;
        m_state_state Q_xx;
        m_dof_dof Q_uu;
        m_dof_state Q_ux;

//        cout << "iteration: " << t << endl;
//        cout << "V_x " << V_x << endl;
//        cout << "V_xx " << V_xx << endl;
//        cout << "A[t] " << A[t] << endl;
//        cout << "B[t] " << B[t] << endl;
//        cout << "l_x[t] " << l_x[t] << endl;
//        cout << "l_xx[t] " << l_xx[t] << endl;
//        cout << "l_u[t] " << l_u[t] << endl;
//        cout << "l_uu[t] " << l_uu[t] << endl;

        Q_u = l_u[t] + (B[t].transpose() * V_x);

        Q_x = l_x[t] + (A[t].transpose() * V_x);

        Q_ux = (B[t].transpose() * (V_xx * A[t]));

        Q_uu = l_uu[t] + (B[t].transpose() * (V_xx * B[t]));

        Q_xx = l_xx[t] + (A[t].transpose() * (V_xx * A[t]));

//        cout << "A " << A[t] << endl;
//        cout << "B  " << B[t] << endl;
//
//        cout << "Q_x " << Q_x << endl;
//        cout << "Q_u " << Q_u << endl;
//
//        cout << "Q_xx " << Q_xx << endl;
//        cout << "Q_uu " << Q_uu << endl;
//        cout << "Q_ux " << Q_ux << endl;

        m_dof_dof Q_uu_reg = Q_uu.replicate(1, 1);

        //cout << "Q_uu_reg " << Q_uu_reg << endl;

        for(int i = 0; i < NUM_DOF; i++){
            Q_uu_reg(i, i) += lamb;
        }



        if(Q_uu_reg(0, 0) < 0){
            return false;
        }

        auto temp = (Q_uu_reg).ldlt();
        m_dof_dof I;
        I.setIdentity();
        m_dof_dof Q_uu_inv = temp.solve(I);


        //cout << "Q_uu_inv " << Q_uu_inv << endl;

        // Caluclate Q_uu_inverse via eigen vector regularisation
//        SelfAdjointEigenSolver<m_dof_dof> eigensolver(Q_uu);
//
//        m_dof eigVals = eigensolver.eigenvalues();
//        m_dof_dof eigVecs = eigensolver.eigenvectors();
//
//        for(int i = 0; i < NUM_DOF; i++){
//            if(eigVals(i) < 0) eigVals(i) = 0.0f;
//            eigVals(i) += lamb;
//            eigVals(i) = 1 / eigVals(i);
//        }
//
//        m_dof_dof diagMat;
//
//        diagMat = eigVals.asDiagonal();
//
//        m_dof_dof Q_uu_inv;
//        Q_uu_inv = eigVecs * diagMat * eigVecs.transpose();


        k[t] = -Q_uu_inv * Q_u;
        K[t] = -Q_uu_inv * Q_ux;
//        cout << "k[t] " << k[t] << endl;
//        cout << "K[t] " << K[t] << endl;

        V_x = Q_x + (K[t].transpose() * (Q_uu * k[t])) + (K[t].transpose() * Q_u) + (Q_ux.transpose() * k[t]);
        V_xx = Q_xx + (K[t].transpose() * (Q_uu * K[t])) + (K[t].transpose() * Q_ux) + (Q_ux.transpose() * K[t]);

        V_xx = (V_xx + V_xx.transpose()) / 2;

        int df = 1;

        if(t > ILQR_HORIZON_LENGTH){
            cout << "---------------------------------------" << endl;
            cout << "control number  " << t << endl;
            cout << "k[t] " << k[t] << endl;
            cout << "K[t]" << K[t] << endl;
            cout << "l_x[t] " << endl << l_x[t] << endl;
            cout << "l_xx[t] " << endl  << l_xx[t] << endl;
            cout << "Q_ux" << Q_ux<< endl;
            cout << "V_xx" << V_xx<< endl;
            cout << "B" << B[t] << endl;
            cout << "A" << A[t]<< endl;
            int a = 1;

        }
    }

    return true;
}

void backwardsPass_Vxx_reg(m_state_state *A, m_state_dof *B, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K){
    m_state V_x;
    V_x = l_x[ILQR_HORIZON_LENGTH - 1];
    m_state_state V_xx;
    V_xx = l_xx[ILQR_HORIZON_LENGTH - 1];

    for(int t = ILQR_HORIZON_LENGTH - 1; t > -1; t--){
        m_state Q_x;
        m_dof Q_u;
        m_state_state Q_xx;
        m_dof_dof Q_uu;
        m_dof_state Q_ux;
        m_state_state V_xx_reg;

//        cout << "V_x " << V_x << endl;
//        cout << "V_xx " << V_xx << endl;

        V_xx_reg = V_xx.replicate(1, 1);
        for(int i = 0; i < NUM_STATES; i++){
            V_xx_reg(i, i) += lamb;
        }

        Q_x = l_x[t] + (A[t].transpose() * V_x);

        Q_u = l_u[t] + (B[t].transpose() * V_x);

        Q_xx = l_xx[t] + (A[t].transpose() * (V_xx * A[t]));

        Q_uu = l_uu[t] + (B[t].transpose() * (V_xx * B[t]));

        Q_ux = (B[t].transpose() * (V_xx * A[t]));


//        cout << "A " << A[t] << endl;
//        cout << "B  " << B[t] << endl;
//
//        cout << "Q_x " << Q_x << endl;
//        cout << "Q_u " << Q_u << endl;
//
//        cout << "Q_xx " << Q_xx << endl;
//        cout << "Q_uu " << Q_uu << endl;
//        cout << "Q_ux " << Q_ux << endl;

        m_dof_dof Q_uu_reg;
        m_dof_state Q_ux_reg;

        Q_uu_reg = l_uu[t] + (B[t].transpose() * (V_xx_reg * B[t]));

        Q_ux_reg = (B[t].transpose() * (V_xx_reg * A[t]));

        auto temp = (Q_uu_reg).ldlt();
        m_dof_dof I;
        I.setIdentity();
        m_dof_dof Q_uu_inv = temp.solve(I);
        //cout << "Q_uu_inv " << Q_uu_inv << endl;

        k[t] = -Q_uu_inv * Q_u;
        K[t] = -Q_uu_inv * Q_ux_reg;
//        cout << "k[t] " << k[t] << endl;
//        cout << "K[t] " << K[t] << endl;

        V_x = Q_x + (K[t].transpose() * (Q_uu * k[t])) + (K[t].transpose() * Q_u) + (Q_ux.transpose() * k[t]);

        V_xx = Q_xx + (K[t].transpose() * (Q_uu * K[t])) + (K[t].transpose() * Q_ux) + (Q_ux.transpose() * K[t]);


        //V_x = Q_x - (Q_u * Q_uu_inv * Q_ux).transpose();
        //cout << "V_x " << V_x << endl;
        //V_xx = Q_xx + (Q_ux.transpose() * Q_uu_inv * Q_ux);
        //cout << "V_xx " << V_xx << endl;
        int df = 1;



        if(t > ILQR_HORIZON_LENGTH){
            cout << "---------------------------------------" << endl;
            cout << "control number  " << t << endl;
            cout << "k[t] " << k[t] << endl;
            cout << "K[t]" << K[t] << endl;
            cout << "l_x[t] " << endl << l_x[t] << endl;
            cout << "l_xx[t] " << endl  << l_xx[t] << endl;
            cout << "Q_ux" << Q_ux<< endl;
            cout << "V_xx" << V_xx<< endl;
            cout << "B" << B[t] << endl;
            cout << "A" << A[t]<< endl;
            int a = 1;

        }
    }
}

float forwardsPass(m_dof *k, m_dof_state *K, float oldCost){
    m_dof *U_new = new m_dof[ILQR_HORIZON_LENGTH];
    float alpha = 1.0;
    float newCost = 0;
    bool costReduction = false;

    //cout << "-------------------- START FORWARDS PASS ------------------------" << endl;
    while(!costReduction){
        cpMjData(model, mdata, d_init);
        newCost = 0;
        m_state stateFeedback;
        m_state X;
        m_state X_new;
        m_dof U_last = returnStateControls(dArray[0]);

        // calculate initial trajectory
        U_new[0] = U_last + (alpha * k[0]);

        for(int k = 0; k < NUM_DOF; k++){
            mdata->ctrl[k] = U_new[0](k);
        }

        for(int t = 0; t < ILQR_HORIZON_LENGTH - 1; t++){



            X = returnState(dArray[t]);
            X_new = returnState(mdata);
            U_last = returnStateControls(dArray[t]);

            stateFeedback = X_new - X;
            m_dof feedBackGain = K[t] * stateFeedback;

            U_new[t] = U_last + (alpha * k[t]) + feedBackGain;

                //            cout << "old state was: " << endl << X << endl;
//            cout << "new state is: " << endl << X_new << endl;
//
//            cout << "feedback : " << endl << stateFeedback << endl;
//            cout << "feedback gain: " << endl << feedBackGain << endl;

//                cout << "old control was: " << endl << U_last << endl;
//                cout << "new control is: " << endl << U_new[controlNumber] << endl;

                // constrain new torque within limits
            for(int k = 0; k < NUM_DOF; k++){
                //if(U_new[t](k) > torqueLims[k]) U_new[t](k) = torqueLims[k];
                //if(U_new[t](k) < -torqueLims[k]) U_new[t](k) = -torqueLims[k];

                mdata->ctrl[k] = U_new[t](k);
            }

            // calc current state cost and keep running tally
            float currentCost = calcStateCost(mdata);

            newCost += (currentCost * ILQR_DT);

            for(int i = 0; i < NUM_MJSTEPS_PER_CONTROL; i++){
                mj_step(model, mdata);
            }

        }

        if(newCost < oldCost){
            costReduction = true;
        }
        else{
            alpha = alpha - 0.1;
            if(alpha <= 0){
                newCost = oldCost;
                break;
            }
        }
    }

    // only update array of mujoco data structure if resulting trajectory is better
    if(newCost < oldCost){

        cpMjData(model, mdata, d_init);

        for(int k = 0; k < NUM_DOF; k++){
            mdata->ctrl[k] = U_new[0](k);
        }

        cpMjData(model, dArray[0], mdata);

        for(int i = 0; i < ILQR_HORIZON_LENGTH; i++){

            for(int k = 0; k < NUM_DOF; k++){
                mdata->ctrl[k] = U_new[i](k);
            }
            cpMjData(model, dArray[i], mdata);

            for(int i = 0; i < NUM_MJSTEPS_PER_CONTROL; i++){
                mj_step(model, mdata);
            }

        }
        cpMjData(model, dArray[ILQR_HORIZON_LENGTH], mdata);
    }

    m_state termStateBest = returnState(dArray[ILQR_HORIZON_LENGTH]);
    cout << "terminal state best: " << endl << termStateBest << endl;
    cout << "middle control 250 was: " << endl << returnStateControls(dArray[250]) << endl;
    cout << "best alpha was " << alpha << endl;
    cout << "best cost was " << newCost << endl;
    cout << "-------------------- END FORWARDS PASS ------------------------" << endl;

    return newCost;
}

bool checkForConvergence(float newCost, float oldCost){
    bool convergence = false;
    m_state terminalState = returnState(dArray[ILQR_HORIZON_LENGTH]);

    std::cout << "--------------------------------------------------" <<  std::endl;
    std::cout << "New cost: " << newCost <<  std::endl;
    m_state X_diff = terminalState- X_desired;
    std::cout << "terminal state is " << endl << terminalState << endl;

    numIterations++;
    float costGrad = (oldCost - newCost)/newCost;

    if((numIterations > 5) && costGrad < epsConverge) {
        convergence = true;
        cout << "ilQR converged, num Iterations: " << numIterations << " final cost: " << newCost << endl;
    }

    return convergence;
}

void lineariseDynamicsParallel(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum){
    const int degreesOfFreedom = 2;
    const int controls = 1;

    using x_t = Differentiator<degreesOfFreedom, controls>::x_t;
    using u_t = Differentiator<degreesOfFreedom, controls>::u_t;

    // linearize around set state and control
    stepCostFn_t stepCostFn = stepCost;
//    cout << " control " << controlNum << " is: " << dArray[controlNum]->ctrl[0] << endl;

    if(controlNum == 250){
        m_state termStateBest = returnState(dArray[ILQR_HORIZON_LENGTH]);
        cout << "terminal state best: " << endl << termStateBest << endl;
        cout << "middle control 250 was: " << endl << returnStateControls(dArray[controlNum])(0) << endl;
    }


    Differentiator<degreesOfFreedom, controls>* differentiator = new Differentiator<degreesOfFreedom, controls>(model, dArray[controlNum], stepCostFn);
    differentiator->setMJData(dArray[controlNum]);
    differentiator->updateDerivatives();



    _A = *(differentiator->A);
    _B = *(differentiator->B);

//    if(controlNum % 50 == 0){
//        cpMjData(model, mdata, dArray[controlNum]);
//        mj_forward(model, mdata);
//        // get framebuffer viewport
//        mjrRect viewport = { 0, 0, 0, 0 };
//        glfwGetFramebufferSize(window, &viewport.width, &viewport.height);
//
//        // update scene and render
//        mjv_updateScene(model, mdata, &opt, NULL, &cam, mjCAT_ALL, &scn);
//        mjr_render(viewport, &scn, &con);
//
//        // swap OpenGL buffers (blocking call due to v-sync)
//        glfwSwapBuffers(window);
//
//        // process pending GUI events, call GLFW callbacks
//        glfwPollEvents();
//        cout << "A matrix is: " << _A << endl;
//        int a = 1;
//    }
}

void lineariseDynamicsSerial(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum){

    float eps = 1e-5;

    // calculate A matrix
    m_state X_copy = returnState(dArray[controlNum]);
    mjData *saveData;
    saveData = mj_makeData(model);
    cpMjData(model, saveData, dArray[controlNum]);

    for(int i = 0; i < NUM_STATES; i++){
        m_state decX = X_copy.replicate(1, 1);
        m_state incX = X_copy.replicate(1, 1);

        decX(i) -= eps;
        incX(i) += eps;

        m_state stateInc;
        m_state stateDec;

        setState(saveData, decX);
        mj_step(model, saveData);
        stateDec = returnState(saveData);

        // Undo pertubation
        cpMjData(model, saveData, dArray[controlNum]);

        setState(saveData, incX);
        mj_step(model, saveData);
        stateInc = returnState(saveData);

        // Undo pertubation
        cpMjData(model, saveData, dArray[controlNum]);

        for(int j = 0; j < NUM_STATES; j++){
            _A(j, i) = (stateInc(j) - stateDec(j)) / (2 * eps);
        }

    }

    // calculate B matrix
    for(int i = 0; i < NUM_DOF; i++){
        m_dof decControl = returnStateControls(dArray[controlNum]);
        m_dof incControl = returnStateControls(dArray[controlNum]);

        decControl(i) -= eps;
        incControl(i) += eps;

        m_state stateInc;
        m_state stateDec;

        setControl(saveData, decControl);
        mj_step(model, saveData);
        stateDec = returnState(saveData);

        // Undo pertubation
        cpMjData(model, saveData, dArray[controlNum]);

        setControl(saveData, incControl);
        mj_step(model, saveData);
        stateInc = returnState(saveData);

        // Undo pertubation
        cpMjData(model, saveData, dArray[controlNum]);

        for(int j = 0; j < NUM_STATES; j++){
            _B(j, i) = (stateInc(j) - stateDec(j)) / (2 * eps);
        }
    }

    mj_deleteData(saveData);

    //cout << "A matrix is: " << _A << endl;
    //cout << " B Mtrix is: " << _B << endl;

}

float calcStateCost(mjData *d){
    float stateCost;
    m_state X_diff;
    m_state X;
    m_dof U;

    X = returnState(d);
    U = returnStateControls(d);

    VectorXd temp(1);

    // actual - desired
    X_diff = X - X_desired;
//    float power = controlNum / numControls;
//    float scalar = pow(terminalScalarConstant, power);

    temp = (X_diff.transpose() * Q * X_diff) + (U.transpose() * R * U);

    stateCost = temp(0);

    return stateCost;
}

float calcCostDerivatives(Ref<m_state> l_x, Ref<m_state_state> l_xx, Ref<m_dof> l_u, Ref<m_dof_dof> l_uu, int controlNum){
    m_state X_diff;
    m_state X;
    m_dof U;

    X = returnState(dArray[controlNum]);
    U = returnStateControls(dArray[controlNum]);

    // actual - desired
    X_diff = X - X_desired;

//    float power = (float) controlNum / numControls;
//    float scalar = pow(terminalScalarConstant, power);
    l_x = 2 * Q * X_diff;
    l_xx = 2 * Q;

    l_u = 2 * R * U;
    l_uu = 2 * R;
}

m_dof returnStateControls(mjData *data){
    m_dof controls;
    for(int i = 0; i < NUM_DOF; i++){
        controls(i) = data->ctrl[i];
    }

    return controls;
}

m_state returnState(mjData *data){
    m_state state;

    for(int i = 0; i < model->nv; i++){
        state(i) = data->qpos[i];
        state(i + 2) = data->qvel[i];
    }

    return state;
}

void setState(mjData *data, m_state desiredState){
    for(int i = 0; i < model->nv; i++){
        data->qpos[i] = desiredState(i);
        data->qvel[i] = desiredState(i+2);
    }
}

void setControl(mjData *data, m_dof desiredControl){
    for(int i = 0; i < NUM_DOF; i++){
        data->ctrl[i] = desiredControl(i);
    }
}

void initCostMatrices(){
    R.setIdentity();
    for(int i = 0; i < NUM_DOF; i++){
        R(i, i) = controlCost[i];
    }
    Q.setIdentity();
    for(int i = 0; i < NUM_STATES; i++){
        Q(i, i) = stateCosts[i];
    }

}

void initDesiredState(){
    // 0 , 0 is upright for pendulum
    // Pi, 0 is facing towards the ground
    X_desired << 0, 0, 0, 0;

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

m_dof returnNextControl(int controlNum, bool finalControl){
    m_dof nextControl;

    if(finalControl){
        nextControl = finalControls[controlNum];
    }
    else{
        nextControl = initControls[controlNum];
    }

    return nextControl;
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
