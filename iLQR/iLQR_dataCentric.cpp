//
// Created by david on 13/04/2022.
//

#include "iLQR_dataCentric.h"


iLQR* optimiser;

iLQR::iLQR(mjModel* m, mjData* d, m_state _X0, frankaModel* _modelTranslator, MujocoController* _mujocoController){
    numIterations = 0;
    lamda = 0.1;

    for(int i=0; i < ILQR_HORIZON_LENGTH; i++) {
        initControls.push_back(m_ctrl());
        finalControls.push_back(m_ctrl());
    }

    for(int i = 0; i < ILQR_HORIZON_LENGTH; i++){
        f_x.push_back(m_state_state());
        f_u.push_back(m_state_ctrl());

        A_scaled.push_back(m_state_state());
        B_scaled.push_back(m_state_ctrl());

        l_x.push_back(m_state());
        l_xx.push_back(m_state_state());
        l_u.push_back(m_ctrl());
        l_uu.push_back(m_ctrl_ctrl());

        k.push_back(m_ctrl());
        K.push_back(m_ctrl_state());

        U_new.push_back(m_ctrl());
        X_new.push_back(m_state());
    }

    // Extra as one more state than controls
    l_x.push_back(m_state());
    l_xx.push_back(m_state_state());

    // Initialise internal iLQR model and data
    model = m;
    mujocoController = _mujocoController;
    modelTranslator = _modelTranslator;
    mdata = mj_makeData(model);
    cpMjData(model, mdata, d);
    X0 = _X0.replicate(1,1);


}

void iLQR::optimise(){
    bool optimisationFinished = false;
    int numIterations = 0;
    float newCost = 0;
    float oldCost = 1000;

    // iterate until optimisation finished, convergence or if lamda > maxLamda
    for(int i = 0; i < maxIterations; i++){
        numIterations++;

        auto start = high_resolution_clock::now();

        // Linearise the dynamics and save cost values at each state
        // STEP 1 - Linearise dynamics and calculate cost quadratics at every time step
        getDerivatives();

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Linearising model: " << duration.count()/1000 << " milliseconds" << endl;

        bool validBackPass = false;
        bool lamdaExit = false;

        // STEP 2 - Backwards pass to compute optimal linear and feedback gain matrices k and K
        // Keep running this step until the backwards pass does NOT encounter a non_PD Q_uu_reg matrix
        while(!validBackPass) {

            validBackPass = backwardsPass_Quu_reg();

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
            // STEP 3 - Forwards pass to calculate new optimal controls - with optional alpha backtracking line search
            newCost = forwardsPass(oldCost);

            // STEP 4 - Check for convergence
            optimisationFinished = checkForConvergence(newCost, oldCost);
            if(optimisationFinished){
                break;
            }

            oldCost = newCost;
        }
        else{
            cout << "optimisation exited after lamda exceed lamda max, iteration: " << i << endl;
            break;
        }
    }

    for(int i = 0; i < ILQR_HORIZON_LENGTH; i++){
        finalControls[i] = modelTranslator->returnControls(dArray[i]);
        X_new[i] = modelTranslator->returnState(dArray[i]);
    }

}

void iLQR::getDerivatives(){

    MatrixXd A = ArrayXXd::Zero((2 * DOF), (2 * DOF));
    MatrixXd B = ArrayXXd::Zero((2 * DOF), NUM_CTRL);

//    MatrixXd A_dt = ArrayXXd::Zero((2 * DOF), (2 * DOF));
//    MatrixXd B_dt = ArrayXXd::Zero((2 * DOF), NUM_CTRL);

    int save_iterations = model->opt.iterations;
    mjtNum save_tolerance = model->opt.tolerance;

    model->opt.iterations = 30;
    model->opt.tolerance = 0;

    // Linearise the dynamics along the trajectory
    for(int t = 0; t < ILQR_HORIZON_LENGTH; t++){
        // Calculate linearised dynamics for current time step via finite differencing
        lineariseDynamicsSerial_trial_step(A, B, t, MUJOCO_DT);
        //lineariseDynamicsSerial_trial_step(A, B, t, MUJOCO_DT);

        //cout << "A: " << endl << A << endl;
        //cout << "B: " << endl << B << endl;
        scaleLinearisation(A_scaled[t], B_scaled[t], A, B, 1);

        if(t == 1000){
            int a  = 1;
        }
//        cout << "A+_scale: " << endl << A_scaled[t] << endl;
//        cout << "B_scale: " << endl << B_scaled[t] << endl;

        f_x[t] = A_scaled[t].replicate(1,1);
        f_u[t] = B_scaled[t].replicate(1,1);

        modelTranslator->costDerivatives(dArray[t], l_x[t], l_xx[t], l_u[t], l_uu[t], t, ILQR_HORIZON_LENGTH);

        if(t % 50 == 0){
            int b = 5;
        }

        l_x[t]  *= ILQR_DT;
        l_xx[t] *= ILQR_DT;
        l_u[t]  *= ILQR_DT;
        l_uu[t] *= ILQR_DT;

    }

    model->opt.iterations = save_iterations;
    model->opt.tolerance = save_tolerance;

    //TODO FIX FACT THAT THERE SHOULD BE NO CONTROL COST AT END OF TRAJECTORY
    modelTranslator->costDerivatives(dArray[ILQR_HORIZON_LENGTH], l_x[ILQR_HORIZON_LENGTH], l_xx[ILQR_HORIZON_LENGTH], l_u[ILQR_HORIZON_LENGTH-1], l_uu[ILQR_HORIZON_LENGTH-1], ILQR_HORIZON_LENGTH, ILQR_HORIZON_LENGTH);
    l_x [ILQR_HORIZON_LENGTH] *= ILQR_DT;
    l_xx[ILQR_HORIZON_LENGTH] *= ILQR_DT;

}


void iLQR::scaleLinearisation(Ref<m_state_state> A_scaled, Ref<m_state_ctrl> B_scaled, Ref<m_state_state> A, Ref<m_state_ctrl> B, int num_steps_per_dt){

    // TODO look into ways of speeding up matrix to the power of calculation
    A_scaled = A.replicate(1, 1);
    B_scaled = B.replicate(1, 1);
    m_state_ctrl currentBTerm;

    for(int i = 0; i < num_steps_per_dt - 1; i++){
        A_scaled *= A;
    }

    currentBTerm = B.replicate(1, 1);
    for(int i = 0; i < num_steps_per_dt - 1; i++){
        currentBTerm = A * currentBTerm;
        B_scaled += currentBTerm;
    }
}

bool iLQR::backwardsPass_Quu_reg(){
    m_state V_x;
    V_x = l_x[ILQR_HORIZON_LENGTH];
    m_state_state V_xx;
    V_xx = l_xx[ILQR_HORIZON_LENGTH];

    for(int t = ILQR_HORIZON_LENGTH - 2; t > -1; t--){
        m_state Q_x;
        m_ctrl Q_u;
        m_state_state Q_xx;
        m_ctrl_ctrl Q_uu;
        m_ctrl_state Q_ux;

        Q_u = l_u[t] + (B_scaled[t].transpose() * V_x);

        Q_x = l_x[t] + (A_scaled[t].transpose() * V_x);

        Q_ux = (B_scaled[t].transpose() * (V_xx * A_scaled[t]));

        Q_uu = l_uu[t] + (B_scaled[t].transpose() * (V_xx * B_scaled[t]));

        Q_xx = l_xx[t] + (A_scaled[t].transpose() * (V_xx * A_scaled[t]));


        m_ctrl_ctrl Q_uu_reg = Q_uu.replicate(1, 1);

        for(int i = 0; i < NUM_CTRL; i++){
            Q_uu_reg(i, i) += lamda;
        }

        if(!isMatrixPD(Q_uu_reg)){
            return false;
        }

        auto temp = (Q_uu_reg).ldlt();
        m_ctrl_ctrl I;
        I.setIdentity();
        m_ctrl_ctrl Q_uu_inv = temp.solve(I);

        k[t] = -Q_uu_inv * Q_u;
        K[t] = -Q_uu_inv * Q_ux;

        V_x = Q_x + (K[t].transpose() * (Q_uu * k[t])) + (K[t].transpose() * Q_u) + (Q_ux.transpose() * k[t]);
        V_xx = Q_xx + (K[t].transpose() * (Q_uu * K[t])) + (K[t].transpose() * Q_ux) + (Q_ux.transpose() * K[t]);

        V_xx = (V_xx + V_xx.transpose()) / 2;

    }

    return true;
}

void iLQR::backwardsPass_Vxx_reg(){
    m_state V_x;
    V_x = l_x[ILQR_HORIZON_LENGTH - 1];
    m_state_state V_xx;
    V_xx = l_xx[ILQR_HORIZON_LENGTH - 1];

    for(int t = ILQR_HORIZON_LENGTH - 1; t > -1; t--){
        m_state Q_x;
        m_ctrl Q_u;
        m_state_state Q_xx;
        m_ctrl_ctrl Q_uu;
        m_ctrl_state Q_ux;
        m_state_state V_xx_reg;

        V_xx_reg = V_xx.replicate(1, 1);
        for(int i = 0; i < (2 * DOF); i++){
            V_xx_reg(i, i) += lamda;
        }

        Q_x = l_x[t] + (A_scaled[t].transpose() * V_x);

        Q_u = l_u[t] + (B_scaled[t].transpose() * V_x);

        Q_xx = l_xx[t] + (A_scaled[t].transpose() * (V_xx * A_scaled[t]));

        Q_uu = l_uu[t] + (B_scaled[t].transpose() * (V_xx * B_scaled[t]));

        Q_ux = (B_scaled[t].transpose() * (V_xx * A_scaled[t]));

        m_ctrl_ctrl Q_uu_reg;
        m_ctrl_state Q_ux_reg;

        Q_uu_reg = l_uu[t] + (B_scaled[t].transpose() * (V_xx_reg * B_scaled[t]));

        Q_ux_reg = (B_scaled[t].transpose() * (V_xx_reg * A_scaled[t]));

        auto temp = (Q_uu_reg).ldlt();
        m_ctrl_ctrl I;
        I.setIdentity();
        m_ctrl_ctrl Q_uu_inv = temp.solve(I);

        k[t] = -Q_uu_inv * Q_u;
        K[t] = -Q_uu_inv * Q_ux_reg;

        V_x = Q_x + (K[t].transpose() * (Q_uu * k[t])) + (K[t].transpose() * Q_u) + (Q_ux.transpose() * k[t]);

        V_xx = Q_xx + (K[t].transpose() * (Q_uu * K[t])) + (K[t].transpose() * Q_ux) + (Q_ux.transpose() * K[t]);

        V_xx = (V_xx + V_xx.transpose()) / 2;

    }
}

bool iLQR::isMatrixPD(Ref<MatrixXd> matrix){
    bool matrixPD = true;
    //TODO implement cholesky decomp for PD check and maybe use result for inverse Q_uu

    Eigen::LLT<Eigen::MatrixXd> lltOfA(matrix); // compute the Cholesky decomposition of the matrix
    if(lltOfA.info() == Eigen::NumericalIssue)
    {
        matrixPD = false;
        //throw std::runtime_error("Possibly non semi-positive definitie matrix!");
    }

    return matrixPD;
}

float iLQR::forwardsPass(float oldCost){
    // TODO check if ths needs to be changed to a standard vector?
    m_ctrl *U_new = new m_ctrl[ILQR_HORIZON_LENGTH];
    float alpha = 1.0;
    float newCost = 0;
    bool costReduction = false;

    while(!costReduction){
        cpMjData(model, mdata, d_init);
        newCost = 0;
        m_state stateFeedback;
        m_state X;
        m_state X_new;
        m_ctrl U_last = modelTranslator->returnControls(dArray[0]);

        // calculate initial trajectory
        U_new[0] = U_last + (alpha * k[0]);

        for(int k = 0; k < NUM_CTRL; k++){
            mdata->ctrl[k] = U_new[0](k);
        }

        for(int t = 0; t < ILQR_HORIZON_LENGTH - 1; t++){

            X = modelTranslator->returnState(dArray[t]);
            X_new = modelTranslator->returnState(mdata);
            U_last = modelTranslator->returnControls(dArray[t]);

            stateFeedback = X_new - X;
            m_ctrl feedBackGain = K[t] * stateFeedback;

            U_new[t] = U_last + (alpha * k[t]) + feedBackGain;
//            cout << "k " << k[t] << endl;
//            cout << "K " << K[t] << endl;
//            cout << "U_last " << U_last << endl;
//            cout << "New U: " << U_new[t] << endl;

            // constrain new torque within limits
            for(int k = 0; k < NUM_CTRL; k++){
                //if(U_new[t](k) > torqueLims[k]) U_new[t](k) = torqueLims[k];
                //if(U_new[t](k) < -torqueLims[k]) U_new[t](k) = -torqueLims[k];

                mdata->ctrl[k] = U_new[t](k);
            }

            // calc current state cost and keep running tally
            float currentCost = modelTranslator->costFunction(mdata, t, ILQR_HORIZON_LENGTH);

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

        for(int k = 0; k < NUM_CTRL; k++){
            mdata->ctrl[k] = U_new[0](k);
        }

        cpMjData(model, dArray[0], mdata);

        for(int i = 0; i < ILQR_HORIZON_LENGTH; i++){

            for(int k = 0; k < NUM_CTRL; k++){
                mdata->ctrl[k] = U_new[i](k);
            }
            cpMjData(model, dArray[i], mdata);

            for(int i = 0; i < NUM_MJSTEPS_PER_CONTROL; i++){
                mj_step(model, mdata);
            }

        }
        cpMjData(model, dArray[ILQR_HORIZON_LENGTH], mdata);
    }

    m_state termStateBest = modelTranslator->returnState(dArray[ILQR_HORIZON_LENGTH]);
//    cout << "terminal state best: " << endl << termStateBest << endl;
//    cout << "best alpha was " << alpha << endl;
//    cout << "best cost was " << newCost << endl;
//    cout << "-------------------- END FORWARDS PASS ------------------------" << endl;

    return newCost;
}

bool iLQR::checkForConvergence(float newCost, float oldCost){
    bool convergence = false;
    m_state terminalState = modelTranslator->returnState(dArray[ILQR_HORIZON_LENGTH]);

    std::cout << "--------------------------------------------------" <<  std::endl;
    std::cout << "New cost: " << newCost <<  std::endl;
    double cubeX = terminalState(7);
    double cubeY = terminalState(8);
    std::cout << "terminal state is, cube X: " << cubeX << " cube Y: " << cubeY << endl;

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
//    for(int i = 0; i < (2 * DOF); i++){
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
//        for(int j = 0; j < (2 * DOF); j++){
//            _A(j, i) = (stateInc(j) - stateDec(j)) / (2 * eps);
//        }
//
//    }
//
//    // calculate B matrix
//    for(int i = 0; i < NUM_CTRL; i++){
//        m_ctrl decControl = returnStateControls(dArray[controlNum]);
//        m_ctrl incControl = returnStateControls(dArray[controlNum]);
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
//        for(int j = 0; j < (2 * DOF); j++){
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
    static int nwarmup = 3;
    float eps = 1e-2;

    _A.block(0, 0, DOF, DOF).setIdentity();
    _A.block(0, DOF, DOF, DOF).setIdentity();
    _A.block(0, DOF, DOF, DOF) *= model->opt.timestep;
    _B.setZero();

    // Initialise matrices for forwards dynamics
    m_dof_dof dqaccdq;
    m_dof_dof dqaccdqvel;
    m_dof_ctrl dqaccdctrl;

    m_dof acellDec;
    m_dof acellInc;

    // Create a copy of the current data that we want to differentiate around
    mjData *saveData;
    saveData = mj_makeData(model);
    cpMjData(model, saveData, dArray[controlNum]);

    // Allocate memory for variables
    mjtNum* warmstart = mj_stackAlloc(saveData, DOF);

    // Compute mj_forward once with no skips
    mj_forward(model, saveData);

    // Compute mj_forward a few times to allow optimiser to get a more accurate value for qacc
    // skips position and velocity stages (TODO LOOK INTO IF THIS IS NEEDED FOR MY METHOD)
    for( int rep=1; rep<nwarmup; rep++ )
        mj_forwardSkip(model, saveData, mjSTAGE_VEL, 1);


    // save output for center point and warmstart (needed in forward only)
    mju_copy(warmstart, saveData->qacc_warmstart, DOF);

    // CALCULATE dqaccdctrl
    for(int i = 0; i < NUM_CTRL; i++){
        saveData->ctrl[i] = dArray[controlNum]->ctrl[i] + eps;

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_VEL, 1);

        // copy and store +perturbation
        acellInc = modelTranslator->returnAccelerations(saveData);

        // perturb selected target -
        cpMjData(model, saveData, dArray[controlNum]);
        saveData->ctrl[i] = dArray[controlNum]->ctrl[i] - eps;

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_VEL, 1);

        acellDec = modelTranslator->returnAccelerations(saveData);

        for(int j = 0; j < DOF; j++){
            dqaccdctrl(j, i) = (acellInc(j) - acellDec(j))/(2*eps);
        }

        // undo pertubation
        cpMjData(model, saveData, dArray[controlNum]);

    }

    // CALCULATE dqaccdvel
    for(int i = 0; i < DOF; i++){
        // perturb velocity +

        modelTranslator->perturbVelocity(saveData, dArray[controlNum], i, eps);

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_POS, 1);

        // copy and store +perturbation
        acellInc = modelTranslator->returnAccelerations(saveData);
        //cout << "acellInc " << endl << acellInc << endl;

        // perturb velocity -
        modelTranslator->perturbVelocity(saveData, dArray[controlNum], i, -eps);

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_POS, 1);

        acellDec = modelTranslator->returnAccelerations(saveData);
        //cout << "acellDec " << endl << acellDec << endl;

        // compute column i of derivative 1
        for(int j = 0; j < DOF; j++){
            dqaccdqvel(j, i) = (acellInc(j) - acellDec(j))/(2*eps);
        }
        //cout << "dq/dvel " << endl << dqaccdqvel << endl;

        // undo perturbation
        cpMjData(model, saveData, dArray[controlNum]);
    }

    // CALCULATE dqaccdcqpos
    for(int i = 0; i < DOF; i++){
        // perturb position +
        modelTranslator->perturbPosition(saveData, dArray[controlNum], i, eps);
//        cout << "qacc Number 1: ";
//        for(int j = 0; j < NUM_CTRL; j++){
//            cout << saveData->qacc[j] << " ";
//        }
//        cout << endl;

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_NONE, 1);
        acellInc = modelTranslator->returnAccelerations(saveData);
        //cout << "accel inc is: " << endl << acellInc << endl;
//        for( int rep=1; rep<5; rep++ ){
//            mj_forwardSkip(model, saveData, mjSTAGE_NONE, 1);
//            acellInc = modelTranslator->returnAccelerations(saveData);
//            cout << "accel inc is: " << endl << acellInc << endl;
//        }


        // copy and store +perturbation
        acellInc = modelTranslator->returnAccelerations(saveData);
//        cout << "qacc Number 1: ";
//        for(int j = 0; j < NUM_CTRL; j++){
//            cout << saveData->qacc[j] << " ";
//        }
//        cout << endl;

        // perturb position -
        modelTranslator->perturbPosition(saveData, dArray[controlNum], i, -eps);

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_forwardSkip(model, saveData, mjSTAGE_NONE, 1);

        // Additional mj_forward steps when computing dqac/dqpos
        for( int rep=1; rep<nwarmup; rep++ )
            mj_forwardSkip(model, saveData, mjSTAGE_NONE, 1);

        acellDec = modelTranslator->returnAccelerations(saveData);
        //cout << "accel inc is: " << endl << acellInc << endl;
//        cout << "qacc from output ";
//        for(int j = 0; j < NUM_CTRL; j++){
//            cout << output[j] << " ";
//        }
//        cout << endl;
//


        // compute column i of derivative 1
        for(int j = 0; j < DOF; j++){
            dqaccdq(j, i) = (acellInc(j) - acellDec(j))/(2*eps);
        }

        // undo perturbation
        cpMjData(model, saveData, dArray[controlNum]);
    }

    mj_deleteData(saveData);

//    cout << " dqaccdqis: " << dqaccdq << endl;
//    cout << " dqaccdqvel: " << dqaccdqvel << endl;
    //cout << " dqaccdctrl: " << dqaccdctrl << endl;

    _A.block(DOF, 0, DOF, DOF) = dqaccdq * ilqr_dt;
    _A.block(DOF, DOF, DOF, DOF).setIdentity();
    _A.block(DOF, DOF, DOF, DOF) += dqaccdqvel * ilqr_dt;
    _B.block(DOF, 0, DOF, NUM_CTRL) = dqaccdctrl * ilqr_dt;

    //cout << "A matrix is: " << _A << endl;
//    cout << " B Mtrix is: " << _B << endl;

}

void iLQR::lineariseDynamicsSerial_trial_step(Ref<MatrixXd> _A, Ref<MatrixXd> _B, int controlNum, float dt){
    // Initialise variables
    static int nwarmup = 3;
    // best found so far, 1e-1 eps, 1 stpe linearisation, 0.002 control size

    //float epsControls[NUM_CTRL] = {1e-1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    //float epsVelocities[DOF] = {1e-1, 1e-1, 1e-1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1e-1};
    float eps = 1e-2;
    //int numStepsLinearisation = NUM_MJSTEPS_PER_CONTROL * 1;
    int numStepsLinearisation = 1;

    _A.block(0, 0, DOF, DOF).setIdentity();
    _A.block(0, DOF, DOF, DOF).setIdentity();
    _A.block(0, DOF, DOF, DOF) *= model->opt.timestep;
    _B.setZero();

    // Initialise matrices for forwards dynamics
    m_dof_dof dqaccdq;
    m_dof_dof dqaccdqvel;
    m_dof_ctrl dqaccdctrl;

    m_dof velDec;
    m_dof velInc;

    // Create a copy of the current data that we want to differentiate around
    mjData *saveData;
    saveData = mj_makeData(model);
    cpMjData(model, saveData, dArray[controlNum]);

    // Allocate memory for variables
    mjtNum* warmstart = mj_stackAlloc(saveData, DOF);

    //cout << "accel before: " << saveData->qacc[0] << endl;
    // Compute mj_forward once with no skips
    mj_forward(model, saveData);
    //cout << "accel before: " << saveData->qacc[0] << endl;

    // Compute mj_forward a few times to allow optimiser to get a more accurate value for qacc
    // skips position and velocity stages (TODO LOOK INTO IF THIS IS NEEDED FOR MY METHOD)
    for( int rep=1; rep<nwarmup; rep++ )
        mj_forwardSkip(model, saveData, mjSTAGE_VEL, 1);


    // save output for center point and warmstart (needed in forward only)
    mju_copy(warmstart, saveData->qacc_warmstart, DOF);

    // CALCULATE dqaccdctrl
    for(int i = 0; i < NUM_CTRL; i++){
        saveData->ctrl[i] = dArray[controlNum]->ctrl[i] + eps;

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        for(int j = 0; j < numStepsLinearisation; j++){
            mj_step(model, saveData);
        }


        // copy and store +perturbation
        velInc = modelTranslator->returnVelocities(saveData);
        //cout << "velInc " << endl << velInc << endl;

        // perturb selected target -
        cpMjData(model, saveData, dArray[controlNum]);
        saveData->ctrl[i] = dArray[controlNum]->ctrl[i] - eps;

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        for(int j = 0; j < numStepsLinearisation; j++){
            mj_step(model, saveData);
        }

        velDec = modelTranslator->returnVelocities(saveData);
        //cout << "veldec " << endl << velDec << endl;

        for(int j = 0; j < DOF; j++){
            //double diffScaled = (velInc(j) - velDec(j)) / (numStepsLinearisation / NUM_MJSTEPS_PER_CONTROL);
            double diffScaled = (velInc(j) - velDec(j));
            dqaccdctrl(j, i) = diffScaled/(2*eps);
        }
        //cout << "dqaccdctrl " << endl << dqaccdctrl << endl;


        // undo pertubation
        cpMjData(model, saveData, dArray[controlNum]);

    }

//    for(int i = 0; i < NUM_CTRL; i++){
//        saveData->ctrl[i] = dArray[controlNum]->ctrl[i] + eps;
//
//        // evaluate dynamics, with center warmstart
//        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
//        mj_forwardSkip(model, saveData, mjSTAGE_VEL, 1);
//
//        // copy and store +perturbation
//        velInc = modelTranslator->returnAccelerations(saveData);
//        //cout << "velInc " << endl << velInc << endl;
//
//        // perturb selected target -
//        cpMjData(model, saveData, dArray[controlNum]);
//        saveData->ctrl[i] = dArray[controlNum]->ctrl[i] - eps;
//
//        // evaluate dynamics, with center warmstart
//        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
//        mj_forwardSkip(model, saveData, mjSTAGE_VEL, 1);
//
//        velDec = modelTranslator->returnAccelerations(saveData);
//        //cout << "veldec " << endl << velDec << endl;
//
//        for(int j = 0; j < DOF; j++){
//            dqaccdctrl(j, i) = (velInc(j) - velDec(j))/(2*eps);
//        }
//        //cout << "dqaccdctrl " << endl << dqaccdctrl << endl;
//
//
//        // undo pertubation
//        cpMjData(model, saveData, dArray[controlNum]);
//
//    }

    // CALCULATE dqaccdvel
    for(int i = 0; i < DOF; i++){
        // perturb velocity +
        modelTranslator->perturbVelocity(saveData, dArray[controlNum], i, eps);

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        for(int j = 0; j < numStepsLinearisation; j++){
            mj_step(model, saveData);
        }

        // copy and store +perturbation
        velInc = modelTranslator->returnVelocities(saveData);
        //cout << "acellInc " << endl << velInc << endl;

        // perturb velocity -
        modelTranslator->perturbVelocity(saveData, dArray[controlNum], i, -eps);

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        for(int j = 0; j < numStepsLinearisation; j++){
            mj_step(model, saveData);
        }

        velDec = modelTranslator->returnVelocities(saveData);
        //cout << "acellDec " << endl << velDec << endl;

        // compute column i of derivative 1
        for(int j = 0; j < DOF; j++){
            //double diffScaled = (velInc(j) - velDec(j)) / (numStepsLinearisation / NUM_MJSTEPS_PER_CONTROL);
            double diffScaled = (velInc(j) - velDec(j));
            dqaccdqvel(j, i) = diffScaled/(2*eps);
        }
        //cout << "dqaccdqvel " << endl << dqaccdqvel << endl;

        // undo perturbation
        cpMjData(model, saveData, dArray[controlNum]);
    }

     //CALCULATE dqaccdcqpos
    for(int i = 0; i < DOF; i++){
        // perturb position +
        modelTranslator->perturbPosition(saveData, dArray[controlNum], i, eps);

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_step(model, saveData);
        velInc = modelTranslator->returnVelocities(saveData);

        // perturb position -
        modelTranslator->perturbPosition(saveData, dArray[controlNum], i, -eps);

        // evaluate dynamics, with center warmstart
        mju_copy(saveData->qacc_warmstart, warmstart, model->nv);
        mj_step(model, saveData);

        // Additional mj_forward steps when computing dqac/dqpos
//        for( int rep=1; rep<nwarmup; rep++ )
//            mj_forwardSkip(model, saveData, mjSTAGE_NONE, 1);

        velDec = modelTranslator->returnVelocities(saveData);

        // compute column i of derivative 1
        for(int j = 0; j < DOF; j++){
            dqaccdq(j, i) = (velInc(j) - velDec(j))/(2*eps);
        }

        // undo perturbation
        cpMjData(model, saveData, dArray[controlNum]);
    }

    mj_deleteData(saveData);

//    cout << " dqaccdqis: " << dqaccdq << endl;
//    cout << " dqaccdqvel: " << dqaccdqvel << endl;
    //cout << " dqaccdctrl: " << dqaccdctrl << endl;

    _A.block(DOF, 0, DOF, DOF) = dqaccdq;
    _A.block(DOF, DOF, DOF, DOF).setIdentity();
    _A.block(DOF, DOF, DOF, DOF) = dqaccdqvel;
    //_B.block(DOF, 0, DOF, NUM_CTRL) = dqaccdctrl * dt;
    _B.block(DOF, 0, DOF, NUM_CTRL) = dqaccdctrl;

    //cout << "A matrix is: " << _A << endl;
    //cout << " B Mtrix is: " << _B << endl;
    int a = 1;

}

m_ctrl iLQR::returnDesiredControl(int controlIndex, bool finalControl){
    if(finalControl){
        return initControls[controlIndex];
    }
    else{
        return finalControls[controlIndex];
    }
}

void iLQR::setInitControls(std::vector<m_ctrl> _initControls){
    for(int i = 0; i < ILQR_HORIZON_LENGTH; i++){
        initControls[i] = _initControls[i].replicate(1,1);
    }
}

void iLQR::makeDataForOptimisation(){
    // Set initial state and run mj_step several times to stabilise system
    modelTranslator->setState(mdata, X0);
    for(int i = 0; i < 5; i++){
        mj_step(model, mdata);
    }

    d_init = mj_makeData(model);
    cpMjData(model, d_init, mdata);

//    const std::string endEffecName = "franka_gripper";
//    int endEffecId = mj_name2id(model, mjOBJ_BODY, endEffecName.c_str());
//    cout << "id end: " << endEffecId<< endl;
//    m_point startPoint = mujocoController->returnBodyPoint(endEffecId);
//    cout << "start poinmt end: " << startPoint << endl;
//    m_point endPoint;
//    m_point direction;
//    direction << 0.6, 0, 0;
//    endPoint = startPoint + direction;

    for(int i = 0; i <= ILQR_HORIZON_LENGTH; i++){
        // populate dArray with mujoco data objects from start of trajec until end
        dArray[i] = mj_makeData(model);
//        m_point currentEEPoint = mujocoController->returnBodyPoint(endEffecId);
//        MatrixXd Jac = mujocoController->calculateJacobian(endEffecId);
//        cout << "jac is: " << Jac << endl;
//
//        m_ctrl desiredControls;
//        m_point desiredEEForce;
//        m_point diff;
//
//        diff = endPoint - currentEEPoint;
//        desiredEEForce(0) = diff(0) * 40;
//        desiredEEForce(1) = diff(1) * 10;
//        desiredEEForce(2) = diff(2) * 200;
//
//        desiredControls = Jac.transpose() * desiredEEForce;
//        cout << "desiredControls " << desiredControls << endl;

        for(int k = 0; k < NUM_CTRL; k++){
            //initControls[i](k) = desiredControls(k) + mdata->qfrc_bias[k];
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


