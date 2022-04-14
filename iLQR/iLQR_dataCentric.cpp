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
#define NUM_CONTROLS 2500
#define TEST_LINEARISATION 1

float dt = 0.002; // time between controls changing
int numControls = 2500; // number of controls needed as per horizon length and time step
float horizonLength = numControls * dt; // seconds
int mujoco_steps_per_dt = (dt / MUJOCO_TIMESTEP); // Number of mj_steps needed per control time step
int linearising_num_sim_steps = 1;  // How many mujoco steps to use for linearising
bool alphaSearchEnabled = true;     // Whether alpha search is enabled to maximise optimisation at each forwards pass
float maxLamda = 10000;             // Maximum lambda before canceliing optimisation
float minLamda = 0.00001;            // Minimum lamda
float lamdaFactor = 10;             // Lamda multiplicative factor
float epsConverge = 0.005;
bool costFunctionFD = false;
m_dof nextControlSequence;


//float controlCost[DOF] = {0.0001, 0.0001, 0.0001, 0.0001, 0.00005, 0.00005, 0.00005};
float controlCost[NUM_DOF] = {0.01};
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

mjData* dArray[NUM_CONTROLS+1];
mjData* d_init;

m_state X_desired;
m_dof_dof R;
m_state_state Q;

float lamb = 0.1;
int numIterations = 0;

ofstream outputDiffDyn;

std::string diffDynFilename = "diffDyn.csv";

inline mjtNum stepCost(const mjData* d)
{
    mjtNum cost = d->qpos[0];
    return cost;
}

void simpleTest(){
    m_state X0;
    X0 << 3.14, 0, 0, 0;

    for(int i = 0; i < NUM_DOF; i++){
        mdata->qpos[i] = X0(i);
        mdata->qvel[i] = X0(i+2);
    }

    mdata->ctrl[0] = 10;
}

void testILQR(){
    m_state X0;
    X0 << 2.8, 0, 0, 0;

    initCostMatrices();
    initDesiredState();

    //globalMujocoController->setSystemState(X0);

    for(int i = 0; i < NUM_DOF; i++){
        mdata->qpos[i] = X0(i);
        mdata->qvel[i] = X0(i+2);
    }
    m_dof *finalControls = new m_dof[numControls];

//    initControls = new m_dof[numControls];
//
//
//    warmStartControls(initControls, X0);
//
//    for(int i = 0; i < numControls; i++){
//        finalControls[i] = initControls[i];
//    }

//    for(int i = numControls/2; i < numControls; i++){
//        for(int j = 0; j < DOF; j++){
//            initControls[i](j) = 0;
//        }
//    }

    m_state *X_lin = new m_state[numControls + 1];
    m_state *X_dyn = new m_state[numControls + 1];
    X_lin[0] = X0;
    X_dyn[0] = X0;

    m_state diff;
    float sumDiff[NUM_STATES] = {0};
    MatrixXd A_last = ArrayXXd::Zero(NUM_STATES, NUM_STATES);
    MatrixXd B_last = ArrayXXd::Zero(NUM_STATES, NUM_DOF);

    d_init = mj_makeData(model);
    cpMjData(model, d_init, mdata);
    for(int i = 0; i <= NUM_CONTROLS; i++){
        // populate dArray with mujoco data objects from start of trajec until end
        dArray[i] = mj_makeData(model);
        // copy current data into current data array place
        float vel = mdata->qvel[0];
        float kp = 2;
        float ctrlsignal = -2 * vel;
        mdata->ctrl[0] = 0.1;
        cpMjData(model, dArray[i], mdata);
        // step simulation with initialised controls
        mj_step(model, mdata);


        // get framebuffer viewport
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

    }

    cpMjData(model, mdata, d_init);

    m_dof initControls;
    initControls << 0.1;

    if(TEST_LINEARISATION){
        for(int i = 0;  i < numControls; i++){
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

            X_lin[i+1] = (A * X_lin[i]) + (B * initControls);
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

    auto iLQRStart = high_resolution_clock::now();

    iLQR(X0, finalControls, X_dyn);

    auto iLQRStop = high_resolution_clock::now();
    auto iLQRDur = duration_cast<microseconds>(iLQRStop - iLQRStart);

    cout << "iLQR took " << iLQRDur.count()/1000000 << " seconds" << endl;

    //saveTrajecToCSV(finalControls, X_dyn);

    // reset simulation
    // save control sequence to mujoco controller class
    // untick ilqractive

//    iLQRSetControlSequence(finalControls, numControls);
//    mujocoTimesStepsPerControl = mujoco_steps_per_dt;
//    numControlsPerTrajectory = numControls;
//    controlCounter = 0;
//    mujocoTimeStepCounter = 0;
}

void iLQR(m_state X0, m_dof *U, m_state *X){
    bool optimisationFinished = false;
    int numIterations = 0;
    float newCost = 0;
    float oldCost = 1000;

    // Create a copy of the initial state of the system


    // Initialise partial differentiation matrices for all timesteps T
    // for linearised dynamics
    m_state_state *f_x = new m_state_state[numControls];
    m_state_dof *f_u = new m_state_dof[numControls];

    // Quadratic cost partial derivatives
    float l[numControls];
    m_state *l_x = new m_state[numControls + 1];
    m_state_state *l_xx = new m_state_state[numControls + 1];
    m_dof *l_u = new m_dof[numControls];
    m_dof_dof *l_uu = new m_dof_dof[numControls];

    // Initialise state feedback gain matrices
    m_dof *k = new m_dof[numControls];
    m_dof_state *K = new m_dof_state[numControls];

    // Initialise new controls and states storage for evaluation
    m_dof *U_new = new m_dof[numControls];
    m_state *X_new = new m_state[numControls + 1];

    //Perform a forward pass with initialised controls - FILL
    //oldCost = rollOutTrajectory(X0, X, U, numControls);
    //cout << "INITIAL COST FROM WARM START TRAJECTORY IS: " << oldCost << endl;


    // iterate until optimisation finished
    while(!optimisationFinished){

        // Linearise the dynamics and save cost values at each state
        auto start = high_resolution_clock::now();

        // STEP 1 - Linearise dynamics and calculate cost quadratics at every time step
        differentiateDynamics(f_x, f_u, l_x, l_xx, l_u, l_uu);

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Linearising model: " << duration.count()/1000 << " milliseconds" << endl;

        bool costImprovementMade = false;
        while(!costImprovementMade){

            // STEP 2 - Backwards pass to compute optimal linear and feedback gain matrices k and K
            backwardsPass_Vxx_reg(f_x, f_u, l_x, l_xx, l_u, l_uu, k, K);

//            for(int i = 0; i < numControls; i++){
//                cout << "control " << i << " K[i] " << K[i] << endl;
//            }

            // STEP 3 - Forwards pass to calculate new optimal controls - with optional alpha binary search
            newCost = forwardsPass(U_new, k, K, oldCost);

            if(newCost < oldCost){
                costImprovementMade = true;
                // STEP 4 - Check for convergence
                optimisationFinished = checkForConvergence(newCost, oldCost);
                oldCost = newCost;

                if(lamb > minLamda){
                    lamb /= lamdaFactor;
                }

            }
            else{
                lamb *= lamdaFactor;
                if(lamb > maxLamda){
                    optimisationFinished = true;
                    cout << "ilQR finished due to lamda exceeds lamda max " << endl;
                    break;
                }
            }
        }
    }
}

void differentiateDynamics(m_state_state *f_x, m_state_dof *f_u, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu){

    MatrixXd A = ArrayXXd::Zero(NUM_STATES, NUM_STATES);
    MatrixXd B = ArrayXXd::Zero(NUM_STATES, NUM_DOF);

    // Linearise the dynamics along the trajectory
    cout << "------------------------ LINEARISE DYNAMICS -----------------------------" << endl;

    for(int t = 0; t < numControls; t++){
        // Calculate linearised dynamics for current time step via finite differencing
        lineariseDynamicsParallel(A, B, t);

//        m_state_state A_dt;
//        m_state_dof B_dt;
//        int scaling = mujoco_steps_per_dt / linearising_num_sim_steps;
//        scaleLinearisation(A, B, A_dt, B_dt,  scaling);

        f_x[t] = A;
        f_u[t] = B;

        calcCostDerivatives(l_x[t], l_xx[t], l_u[t], l_uu[t], t);

        l_x[t]  *= dt;
        l_xx[t] *= dt;
        l_u[t]  *= dt;
        l_uu[t] *= dt;

//        cout << "------------------- iteration: " << t << " --------------------" << endl;
//        cout << "l_xx[t]: " << l_xx[t] << endl;
//        cout << "l_x[t]: " << l_x[t] << endl;

//        cout << "iteration " << t << endl;
//        cout << " f_x " << f_x[t] << endl;
//        cout << " f_u " << f_u[t] << endl;
        //int a = 1;
//
//        cout << "cost " << l[t] << endl;
//        cout << "cost dif x" << l_x[t] << endl;
//        cout << "cost dif xx" << l_xx[t] << endl;

    }

    cout << "---------------------- END LINEARISE DYNAMICS -----------------------------" << endl;

    //l[numControls] = terminalCost(l_x[numControls], l_xx[numControls], X[numControls]);

//    l   [numControls] *= dt;
//    l_x [numControls] *= dt;
//    l_xx[numControls] *= dt;
}

bool backwardsPass_Quu_reg(m_state_state *A, m_state_dof *B, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K){
    m_state V_x;
    V_x = l_x[numControls - 1];
    m_state_state V_xx;
    V_xx = l_xx[numControls - 1];

    for(int t = numControls - 1; t > -1; t--){
        m_state Q_x;
        m_dof Q_u;
        m_state_state Q_xx;
        m_dof_dof Q_uu;
        m_dof_state Q_ux;
        m_state_state V_xx_reg;

//        cout << "V_x " << V_x << endl;
//        cout << "V_xx " << V_xx << endl;

        V_xx_reg = V_xx.replicate(1, 1);
//        for(int i = 0; i < NUM_STATES; i++){
//            V_xx_reg(i, i) += lamb;
//        }


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

        m_dof_dof Q_uu_reg = Q_uu.replicate(1, 1);

        for(int i = 0; i < NUM_DOF; i++){
            Q_uu_reg(i, i) += lamb;
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


        //V_x = Q_x - (Q_u * Q_uu_inv * Q_ux).transpose();
        //cout << "V_x " << V_x << endl;
        //V_xx = Q_xx + (Q_ux.transpose() * Q_uu_inv * Q_ux);
        //cout << "V_xx " << V_xx << endl;
        int df = 1;



        if(t > numControls){
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
    V_x = l_x[numControls - 1];
    m_state_state V_xx;
    V_xx = l_xx[numControls - 1];

    for(int t = numControls - 1; t > -1; t--){
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



        if(t > numControls){
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

float forwardsPass(m_dof *U_best, m_dof *k, m_dof_state *K, float oldCost){
    m_dof *U_new = new m_dof[numControls];
    //m_state *X_new = new m_state[numControls + 1];
    float alpha = 0.1;
    float bestAlpha = alpha;
    float newCost = 0;
    float bestCost = 1000;
    int numAlphaChecks = 9;

    cout << "-------------------- START FORWARDS PASS ------------------------" << endl;
    for(int i = 0; i < numAlphaChecks; i++){
        cpMjData(model, mdata, d_init);
        newCost = 0;


        for(int t = 0; t < numControls; t++){
            m_state stateFeedback;
            m_state X = returnState(dArray[t]);
            m_state X_new = returnState(mdata);
            m_dof U_last = returnStateControls(dArray[t]);

            stateFeedback = X_new - X;
            m_dof feedBackGain = K[t] * stateFeedback;
            m_dof linearFeedback = (alpha * k[t]);

            // calculate new optimal control for this timestep
            U_new[t] = U_last + (alpha * k[t]) + feedBackGain;

//            if(i == 1 && t == 1){
//                cout << "old control was: " << U_last << endl;
//                cout << "new control is: " << U_new[t] << endl;
//                int a = 1;
//            }

//            cout << "old state was: " << endl << X << endl;
//            cout << "new state is: " << endl << X_new << endl;
//
//            cout << "feedback : " << endl << stateFeedback << endl;
//            cout << "feedback gain: " << endl << feedBackGain << endl;
////
//            cout << "old U was: " << endl << U_last << endl;
//            cout << "new U is: " << endl << U_new[t] << endl;

            // constrain new torque within limits
            for(int k = 0; k < NUM_DOF; k++){
                if(U_new[t](k) > torqueLims[k]) U_new[t](k) = torqueLims[k];
                if(U_new[t](k) < -torqueLims[k]) U_new[t](k) = -torqueLims[k];

                mdata->ctrl[k] = U_new[t](k);
            }

            // calc current state cost and keep running tally
            float currentCost = calcStateCost(mdata);

            newCost += (currentCost * dt);

            // step model
            mj_step(model, mdata);
        }

        if(newCost < bestCost){
            bestCost = newCost;
            for(int j = 0; j < numControls; j++){
                U_best[j] = U_new[j];
            }
            bestAlpha = alpha;
        }
        alpha += 0.1;
    }

    // only update array of mujoco data structure if resulting trajectory is better
    if(bestCost < oldCost){
        cpMjData(model, mdata, d_init);
        for(int i = 0; i < numControls; i++){
            for(int k = 0; k < NUM_DOF; k++){

                mdata->ctrl[k] = U_best[i](k);
                //cout << "current control in copy: " << mdata->ctrl[k] << endl;
            }

            cpMjData(model, dArray[i], mdata);
            //cout << "current control in copy: " << dArray[i]->ctrl[0] << endl;


            mj_step(model, mdata);
            int a = 1;

            if((i % 10 == 0)){
                mjrRect viewport = { 0, 0, 0, 0 };
                glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

                // update scene and render
                mjv_updateScene(model, mdata, &opt, NULL, &cam, mjCAT_ALL, &scn);
                mjr_render(viewport, &scn, &con);

                // swap OpenGL buffers (blocking call due to v-sync)
                glfwSwapBuffers(window);

                // process pending GUI events, call GLFW callbacks
                glfwPollEvents();
            }

        }
        cpMjData(model, dArray[numControls], mdata);
    }

    m_state termStateBest = returnState(dArray[numControls]);
    cout << "terminal state best: " << endl << termStateBest << endl;
//    for(int t = 0; t < numControls; t++){
//        cout << "control " << t << " : " << U_best[t] << endl;
//    }
    cout << "best alpha was " << bestAlpha << endl;
    cout << "best cost was " << bestCost << endl;

    cout << "-------------------- END FORWARDS PASS ------------------------" << endl;

    return bestCost;
}

bool checkForConvergence(float newCost, float oldCost){
    bool convergence = false;
    m_state terminalState = returnState(dArray[numControls]);

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


    Differentiator<degreesOfFreedom, controls>* differentiator = new Differentiator<degreesOfFreedom, controls>(model, dArray[controlNum], stepCostFn);
    differentiator->setMJData(dArray[controlNum]);
    differentiator->updateDerivatives();

    _A = *(differentiator->A);
    _B = *(differentiator->B);
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

void initCostMatrices(){
    R.setIdentity();
    for(int i = 0; i < NUM_DOF; i++){
        R(i, i) = controlCost[i];
    }
    Q.setIdentity();
    for(int i = 0; i < NUM_STATES; i++){
        Q(i, i) = stateCosts[i];
    }

    if(mujoco_steps_per_dt % linearising_num_sim_steps != 0){
        cout << "invalid choice for linearising sim steps" << endl;
    }
}

void initDesiredState(){
    X_desired << 3.14, 0, 0, 0;
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

    for(int i = 0; i < numControls; i++){
        for(int j = 0; j < NUM_STATES; j++){
            outputDiffDyn << X_dyn[i](j) << ",";
            outputDiffDyn << X_lin[i](j) << ",";
            outputDiffDyn << X_lin[i](j) - X_dyn[i](j) << ",";
        }
        outputDiffDyn << endl;
    }

}
