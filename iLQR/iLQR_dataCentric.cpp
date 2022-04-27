//
// Created by david on 13/04/2022.
//

#include "iLQR_dataCentric.h"

//extern MujocoController *globalMujocoController;
//extern mjModel* model;						// MuJoCo model
//extern mjData* mdata;						// MuJoCo data
//extern mjvCamera cam;                   // abstract camera
//extern mjvScene scn;                    // abstract scene
//extern mjvOption opt;			        // visualization options
//extern mjrContext con;				    // custom GPU context
//extern GLFWwindow *window;

iLQR* optimiser;

iLQR::iLQR(mjModel* m, mjData* d, m_state X0){
    numIterations = 0;
    lamda = 0.1;

    // Initialise internal iLQR model and data
    model = m;
    mdata = mj_makeData(model);
    cpMjData(model, mdata, d);

    initCostMatrices();
    initDesiredState();

    // Set initial state and run mj_step several times to stabilise system
    setState(mdata, X0);
    for(int i = 0; i < 20; i++){
        mj_step(model, mdata);
    }

    d_init = mj_makeData(model);
    cpMjData(model, d_init, mdata);

    for(int i = 0; i <= ILQR_HORIZON_LENGTH; i++){
        // populate dArray with mujoco data objects from start of trajec until end
        dArray[i] = mj_makeData(model);
        float vel = mdata->qvel[0];
        float kp = 0.1;
        for(int k = 0; k < NUM_DOF; k++){
            //initControls[i](k) = mdata->qfrc_bias[k];
            initControls[i](k) = -kp * vel;
            mdata->ctrl[k] = initControls[i](k);
        }
        // copy current data into current data array at correct timestep
        cpMjData(model, dArray[i], mdata);

        // step simulation with initialised controls
        for(int i = 0; i < NUM_MJSTEPS_PER_CONTROL; i++){
            mj_step(model, mdata);
        }
    }
    dArray[ILQR_HORIZON_LENGTH] = mj_makeData(model);
    cpMjData(model, dArray[ILQR_HORIZON_LENGTH], mdata);

    // reset mdata back to initial state
    cpMjData(model, mdata, d_init);

}

void iLQR::optimise(m_dof *U, m_state *X){
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

    // iterate until optimisation finished
    while(!optimisationFinished){
        numIterations++;


        auto start = high_resolution_clock::now();

        // Linearise the dynamics and save cost values at each state
        // STEP 1 - Linearise dynamics and calculate cost quadratics at every time step
        getDerivatives(f_x, f_u, l_x, l_xx, l_u, l_uu);

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
                if (lamda < maxLamda) {
                    lamda *= lamdaFactor;
                } else {
                    lamdaExit = true;
                    optimisationFinished = true;
                    break;
                }
            } else {
                if (lamda > minLamda) {
                    lamda /= lamdaFactor;
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

void iLQR::getDerivatives(m_state_state *f_x, m_state_dof *f_u, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu){

    MatrixXd A = ArrayXXd::Zero(NUM_STATES, NUM_STATES);
    MatrixXd B = ArrayXXd::Zero(NUM_STATES, NUM_DOF);

    MatrixXd A_test = ArrayXXd::Zero(NUM_STATES, NUM_STATES);
    MatrixXd B_test = ArrayXXd::Zero(NUM_STATES, NUM_DOF);

    MatrixXd A_dt = ArrayXXd::Zero(NUM_STATES, NUM_STATES);
    MatrixXd B_dt = ArrayXXd::Zero(NUM_STATES, NUM_DOF);



    int save_iterations = model->opt.iterations;
    mjtNum save_tolerance = model->opt.tolerance;

    model->opt.iterations = 30;
    model->opt.tolerance = 0;

    // Linearise the dynamics along the trajectory
    for(int t = 0; t < ILQR_HORIZON_LENGTH; t++){
        // Calculate linearised dynamics for current time step via finite differencing
        //lineariseDynamicsSerial_trial(A_test, B_test, t);
        lineariseDynamicsSerial_trial(A, B, t, MUJOCO_DT);

//        cout << "A truth is: " << endl << A << endl;
//        cout << "A trial is: " << endl << A_test << endl;
//        cout << "B truth is: " << endl << B << endl;
//        cout << "B trial is: " << endl << B_test << endl;

        scaleLinearisation(A_dt, B_dt, A, B, NUM_MJSTEPS_PER_CONTROL);

        f_x[t] = A_dt;
        f_u[t] = B_dt;

        calcCostDerivatives(l_x[t], l_xx[t], l_u[t], l_uu[t], t);

        l_x[t]  *= ILQR_DT;
        l_xx[t] *= ILQR_DT;
        l_u[t]  *= ILQR_DT;
        l_uu[t] *= ILQR_DT;

    }

    model->opt.iterations = save_iterations;
    model->opt.tolerance = save_tolerance;

    calcCostDerivatives(l_x[ILQR_HORIZON_LENGTH], l_xx[ILQR_HORIZON_LENGTH], l_u[ILQR_HORIZON_LENGTH-1], l_uu[ILQR_HORIZON_LENGTH-1], ILQR_HORIZON_LENGTH);
    l_x [ILQR_HORIZON_LENGTH] *= ILQR_DT;
    l_xx[ILQR_HORIZON_LENGTH] *= ILQR_DT;

}

//void simpleTest(){
//    m_state X0;
//    X0 << 2.8, 0, 0, 0;
//
//    for(int i = 0; i < NUM_DOF; i++){
//        mdata->qpos[i] = X0(i);
//        mdata->qvel[i] = X0(i+2);
//    }
//
//    cout << "number of controls: " << model->nu << endl;
//
//    mdata->ctrl[1] = 10;
//
////    m_state_state *f_x = new m_state_state[numControls];
////    m_state_dof *f_u = new m_state_dof[numControls];
////
////    float l[numControls];
////    m_state *l_x = new m_state[numControls + 1];
////    m_state_state *l_xx = new m_state_state[numControls + 1];
////    m_dof *l_u = new m_dof[numControls];
////    m_dof_dof *l_uu = new m_dof_dof[numControls];
////
////    // Initialise state feedback gain matrices
////    m_dof *k = new m_dof[numControls];
////    m_dof_state *K = new m_dof_state[numControls];
////
////    l_xx[numControls] << 0.02, 0, 0, 0,
////                         0, 0.02, 0, 0,
////                         0, 0, 0.01, 0,
////                         0, 0, 0, 0.01;
////
////    int bPStart = 2;
////
////    l_xx[numControls - bPStart] << 0.02, 0, 0, 0,
////            0, 0.02, 0, 0,
////            0, 0, 0.01, 0,
////            0, 0, 0, 0.01;
////
////    l_uu[numControls] << 0.002;
////    l_uu[numControls - bPStart] << 0.002;
////
////    f_x[numControls - bPStart] << 0.891, 0.0561, 0.0967, 0.00192,
////                                  0.154, 0.846, 0.00475, 0.0950,
////                                  -1.92, 0.980, 0.910, 0.0498,
////                                  2.68, -2.716, 0.127, 0.868;
////
////    f_u[numControls - bPStart] << 0.540, 0.014, 9.67, 0.474;
////
////    l_x[numControls] << 0.00437, 0.00916, -0.00239, 0.00517;
////
////    l_x[numControls - bPStart] << 0.0051, 0.0056, -0.0006, 0.012;
////
////    l_u[numControls - bPStart] << 0.00000127;
////
////    // ---------------------------------------------------------------- //
////
////    l_xx[numControls - bPStart - 1] << 0.02, 0, 0, 0,
////            0, 0.02, 0, 0,
////            0, 0, 0.01, 0,
////            0, 0, 0, 0.01;
////
////    l_uu[numControls - bPStart - 1] << 0.002;
////
////    f_x[numControls - bPStart - 1] << 0.895, 0.0542, 0.0966, 0.000917,
////                                    0.146, 0.849, 0.00518, 0.0964,
////                                    -1.87, 0.955, 0.910, 0.0316,
////                                    2.56, -2.669, 0.133, 0.894;
////
////    f_u[numControls - bPStart - 1] << 0.540, 0.0156, 9.66, 0.505;
////
////    l_x[numControls - bPStart - 1] << 0.00498, 0.00318, 0.0026, 0.0110;
////
////    l_u[numControls - bPStart - 1] << -0.00000522;
////
////    // ------------------------------------------------------------------ //
////
////    backwardsPass_Quu_reg(f_x, f_u, l_x, l_xx, l_u, l_uu, k, K);
//}

void iLQR::scaleLinearisation(Ref<m_state_state> A_scaled, Ref<m_state_dof> B_scaled, Ref<m_state_state> A, Ref<m_state_dof> B, int num_steps_per_dt){

    // TODO look into ways of speeding up matrix to the power of calculation
    A_scaled = A.replicate(1, 1);
    B_scaled = B.replicate(1, 1);
    m_state_dof currentBTerm;

    for(int i = 0; i < num_steps_per_dt - 1; i++){
        A_scaled *= A;
    }

    currentBTerm = B.replicate(1, 1);
    B_scaled += currentBTerm;
    for(int i = 0; i < num_steps_per_dt - 1; i++){
        currentBTerm = A * currentBTerm;
        B_scaled += currentBTerm;
    }
}

bool iLQR::backwardsPass_Quu_reg(m_state_state *A, m_state_dof *B, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K){
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

        Q_u = l_u[t] + (B[t].transpose() * V_x);

        Q_x = l_x[t] + (A[t].transpose() * V_x);

        Q_ux = (B[t].transpose() * (V_xx * A[t]));

        Q_uu = l_uu[t] + (B[t].transpose() * (V_xx * B[t]));

        Q_xx = l_xx[t] + (A[t].transpose() * (V_xx * A[t]));


        m_dof_dof Q_uu_reg = Q_uu.replicate(1, 1);

        for(int i = 0; i < NUM_DOF; i++){
            Q_uu_reg(i, i) += lamda;
        }

        if(Q_uu_reg(0, 0) < 0){
            return false;
        }

        auto temp = (Q_uu_reg).ldlt();
        m_dof_dof I;
        I.setIdentity();
        m_dof_dof Q_uu_inv = temp.solve(I);

        k[t] = -Q_uu_inv * Q_u;
        K[t] = -Q_uu_inv * Q_ux;

        V_x = Q_x + (K[t].transpose() * (Q_uu * k[t])) + (K[t].transpose() * Q_u) + (Q_ux.transpose() * k[t]);
        V_xx = Q_xx + (K[t].transpose() * (Q_uu * K[t])) + (K[t].transpose() * Q_ux) + (Q_ux.transpose() * K[t]);

        V_xx = (V_xx + V_xx.transpose()) / 2;

    }

    return true;
}

void iLQR::backwardsPass_Vxx_reg(m_state_state *A, m_state_dof *B, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K){
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

        V_xx_reg = V_xx.replicate(1, 1);
        for(int i = 0; i < NUM_STATES; i++){
            V_xx_reg(i, i) += lamda;
        }

        Q_x = l_x[t] + (A[t].transpose() * V_x);

        Q_u = l_u[t] + (B[t].transpose() * V_x);

        Q_xx = l_xx[t] + (A[t].transpose() * (V_xx * A[t]));

        Q_uu = l_uu[t] + (B[t].transpose() * (V_xx * B[t]));

        Q_ux = (B[t].transpose() * (V_xx * A[t]));

        m_dof_dof Q_uu_reg;
        m_dof_state Q_ux_reg;

        Q_uu_reg = l_uu[t] + (B[t].transpose() * (V_xx_reg * B[t]));

        Q_ux_reg = (B[t].transpose() * (V_xx_reg * A[t]));

        auto temp = (Q_uu_reg).ldlt();
        m_dof_dof I;
        I.setIdentity();
        m_dof_dof Q_uu_inv = temp.solve(I);

        k[t] = -Q_uu_inv * Q_u;
        K[t] = -Q_uu_inv * Q_ux_reg;

        V_x = Q_x + (K[t].transpose() * (Q_uu * k[t])) + (K[t].transpose() * Q_u) + (Q_ux.transpose() * k[t]);

        V_xx = Q_xx + (K[t].transpose() * (Q_uu * K[t])) + (K[t].transpose() * Q_ux) + (Q_ux.transpose() * K[t]);

    }
}

float iLQR::forwardsPass(m_dof *k, m_dof_state *K, float oldCost){
    m_dof *U_new = new m_dof[ILQR_HORIZON_LENGTH];
    float alpha = 1.0;
    float newCost = 0;
    bool costReduction = false;

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
//    cout << "terminal state best: " << endl << termStateBest << endl;
//    cout << "best alpha was " << alpha << endl;
//    cout << "best cost was " << newCost << endl;
//    cout << "-------------------- END FORWARDS PASS ------------------------" << endl;

    return newCost;
}

bool iLQR::checkForConvergence(float newCost, float oldCost){
    bool convergence = false;
    m_state terminalState = returnState(dArray[ILQR_HORIZON_LENGTH]);

    std::cout << "--------------------------------------------------" <<  std::endl;
    std::cout << "New cost: " << newCost <<  std::endl;
    m_state X_diff = terminalState- X_desired;
    std::cout << "terminal state is " << endl << terminalState << endl;

    numIterations++;
    float costGrad = (oldCost - newCost)/newCost;

    if((numIterations > 2) && costGrad < epsConverge) {
        convergence = true;
        cout << "ilQR converged, num Iterations: " << numIterations << " final cost: " << newCost << endl;
    }

    return convergence;
}

//void iLQR::lineariseDynamicsSerial(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum){
//
//    float eps = 1e-5;
//
//    // calculate A matrix
//    m_state X_copy = returnState(dArray[controlNum]);
//    mjData *saveData;
//    saveData = mj_makeData(model);
//    cpMjData(model, saveData, dArray[controlNum]);
//
//    for(int i = 0; i < NUM_STATES; i++){
//        m_state decX = X_copy.replicate(1, 1);
//        m_state incX = X_copy.replicate(1, 1);
//
//        decX(i) -= eps;
//        incX(i) += eps;
//
//        m_state stateInc;
//        m_state stateDec;
//
//        setState(saveData, decX);
//        mj_step(model, saveData);
//        stateDec = returnState(saveData);
//
//        // Undo pertubation
//        cpMjData(model, saveData, dArray[controlNum]);
//
//        setState(saveData, incX);
//        mj_step(model, saveData);
//        stateInc = returnState(saveData);
//
//        // Undo pertubation
//        cpMjData(model, saveData, dArray[controlNum]);
//
//        for(int j = 0; j < NUM_STATES; j++){
//            _A(j, i) = (stateInc(j) - stateDec(j)) / (2 * eps);
//        }
//
//    }
//
//    // calculate B matrix
//    for(int i = 0; i < NUM_DOF; i++){
//        m_dof decControl = returnStateControls(dArray[controlNum]);
//        m_dof incControl = returnStateControls(dArray[controlNum]);
//
//        decControl(i) -= eps;
//        incControl(i) += eps;
//
//        m_state stateInc;
//        m_state stateDec;
//
//        setControl(saveData, decControl);
//        mj_step(model, saveData);
//        stateDec = returnState(saveData);
//
//        // Undo pertubation
//        cpMjData(model, saveData, dArray[controlNum]);
//
//        setControl(saveData, incControl);
//        mj_step(model, saveData);
//        stateInc = returnState(saveData);
//
//        // Undo pertubation
//        cpMjData(model, saveData, dArray[controlNum]);
//
//        for(int j = 0; j < NUM_STATES; j++){
//            _B(j, i) = (stateInc(j) - stateDec(j)) / (2 * eps);
//        }
//    }
//
//    mj_deleteData(saveData);
//
//    //cout << "A matrix is: " << _A << endl;
//    //cout << " B Mtrix is: " << _B << endl;
//
//}

void iLQR::lineariseDynamicsSerial_trial(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum, float ilqr_dt){

    // Initialise variables
    int nv = model->nv;
    int nu = model->nu;
    static int nwarmup = 3;
    float eps = 1e-5;

    _A.block(0, 0, nv, nv).setIdentity();
    _A.block(0, nv, nv, nv).setIdentity();
    _A.block(0, nv, nv, nv) *= model->opt.timestep;
    _B.setZero();

    // Initialise matrices for forwards dynamics
    m_dof_dof dqaccdq;
    m_dof_dof dqaccdqvel;
    m_dof_dof dqaccdctrl;

    // Create a copy of the current data that we want to differentiate around
    mjData *saveData;
    saveData = mj_makeData(model);
    cpMjData(model, saveData, dArray[controlNum]);

    // Allocate memory for variables
    mjtNum* temp = mj_stackAlloc(saveData, nv);
    mjtNum* warmstart = mj_stackAlloc(saveData, nv);

    // Compute mj_forward once with no skips
    mj_forward(model, saveData);

    // Compute mj_forward a few times to allow optimiser to get a more accurate value for qacc
    // skips position and velocity stages (TODO LOOK INTO IF THIS IS NEEDED FOR MY METHOD)
    for( int rep=1; rep<nwarmup; rep++ )
        mj_forwardSkip(model, saveData, mjSTAGE_VEL, 1);

    mjtNum* output = saveData->qacc;

    // save output for center point and warmstart (needed in forward only)
    mju_copy(warmstart, saveData->qacc_warmstart, nv);

    // CALCULATE dqaccdctrl
    for(int i = 0; i < NUM_DOF; i++){
        saveData->ctrl[i] = dArray[controlNum]->ctrl[i] + eps;

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_VEL, 1);

        // copy and store +perturbation
        mju_copy(temp, output, nv);

        // perturb selected target -
        cpMjData(model, saveData, dArray[controlNum]);
        saveData->ctrl[i] = dArray[controlNum]->ctrl[i] - eps;

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_VEL, 1);

        for( int j=0; j<nv; j++ ){
            dqaccdctrl(j, i) = (temp[j] - output[j])/(2*eps);
        }

        // undo pertubation
        cpMjData(model, saveData, dArray[controlNum]);

    }

    // CALCULATE dqaccdvel
    for(int i = 0; i < NUM_DOF; i++){
        // perturb velocity +
        saveData->qvel[i] = dArray[controlNum]->qvel[i] + eps;

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_POS, 1);

        // copy and store +perturbation
        mju_copy(temp, output, nv);

        // perturb velocity -
        saveData->qvel[i] = dArray[controlNum]->qvel[i] - eps;

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_POS, 1);

        // compute column i of derivative 1
        for( int j=0; j<nv; j++ ){
            dqaccdqvel(j, i) = (temp[j] - output[j])/(2*eps);
        }

        // undo perturbation
        cpMjData(model, saveData, dArray[controlNum]);
    }

    // CALCULATE dqaccdcqpos
    for(int i = 0; i < NUM_DOF; i++){
        // perturb position +
        saveData->qpos[i] = dArray[controlNum]->qpos[i] + eps;
//        cout << "qacc Number 1: ";
//        for(int j = 0; j < NUM_DOF; j++){
//            cout << saveData->qacc[j] << " ";
//        }
//        cout << endl;

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_NONE, 1);

        // copy and store +perturbation
        mju_copy(temp, output, nv);
//        cout << "qacc Number 1: ";
//        for(int j = 0; j < NUM_DOF; j++){
//            cout << saveData->qacc[j] << " ";
//        }
//        cout << endl;

        // perturb position -
        saveData->qpos[i] = dArray[controlNum]->qpos[i] - eps;

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_NONE, 1);
//        cout << "qacc Number 1: ";
//        for(int j = 0; j < NUM_DOF; j++){
//            cout << saveData->qacc[j] << " ";
//        }
//        cout << endl;

        // compute column i of derivative 1
        for( int j=0; j<nv; j++ ){
            dqaccdq(j, i) = (temp[j] - output[j])/(2*eps);
        }

        // undo perturbation
        cpMjData(model, saveData, dArray[controlNum]);
    }

    mj_deleteData(saveData);

//    cout << " dqaccdqis: " << dqaccdq << endl;
//    cout << " dqaccdqvel: " << dqaccdqvel << endl;
//    cout << " dqaccdctrl: " << dqaccdctrl << endl;

    _A.block(nv, 0, nv, nv) = dqaccdq * ilqr_dt;
    _A.block(nv, nv, nv, nv).setIdentity();
    _A.block(nv, nv, nv, nv) += dqaccdqvel * ilqr_dt;
    _B.block(nv, 0, nv, nu) = dqaccdctrl * ilqr_dt;

//    cout << "A matrix is: " << _A << endl;
//    cout << " B Mtrix is: " << _B << endl;

}

float iLQR::calcStateCost(mjData *d){
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

float iLQR::calcCostDerivatives(Ref<m_state> l_x, Ref<m_state_state> l_xx, Ref<m_dof> l_u, Ref<m_dof_dof> l_uu, int controlNum){
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

m_dof iLQR::returnStateControls(mjData *data){
    m_dof controls;
    for(int i = 0; i < NUM_DOF; i++){
        controls(i) = data->ctrl[i];
    }

    return controls;
}

m_state iLQR::returnState(mjData *data){
    m_state state;

    for(int i = 0; i < model->nv; i++){
        state(i) = data->qpos[i];
        state(i + 2) = data->qvel[i];
    }

    return state;
}

void iLQR::setState(mjData *data, m_state desiredState){
    for(int i = 0; i < model->nv; i++){
        data->qpos[i] = desiredState(i);
        data->qvel[i] = desiredState(i+2);
    }
}

void iLQR::setControl(mjData *data, m_dof desiredControl){
    for(int i = 0; i < NUM_DOF; i++){
        data->ctrl[i] = desiredControl(i);
    }
}



void iLQR::initCostMatrices(){
    R.setIdentity();
    for(int i = 0; i < NUM_DOF; i++){
        R(i, i) = controlCost[i];
    }
    Q.setIdentity();
    for(int i = 0; i < NUM_STATES; i++){
        Q(i, i) = stateCosts[i];
    }

}

void iLQR::initDesiredState(){
    // 0 , 0 is upright for pendulum
    // Pi, 0 is facing towards the ground
    X_desired << 0, 0, 0, 0;

}

void cpMjData(const mjModel* m, mjData* d_dest, const mjData* d_src){
    d_dest->time = d_src->time;
    mju_copy(d_dest->qpos, d_src->qpos, m->nq);
    mju_copy(d_dest->qvel, d_src->qvel, m->nv);
    mju_copy(d_dest->qacc, d_src->qacc, m->nv);
    mju_copy(d_dest->qacc_warmstart, d_src->qacc_warmstart, m->nv);
    mju_copy(d_dest->qfrc_applied, d_src->qfrc_applied, m->nv);
    mju_copy(d_dest->xfrc_applied, d_src->xfrc_applied, 6*m->nbody);
    mju_copy(d_dest->ctrl, d_src->ctrl, m->nu);
}
