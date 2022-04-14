//
// Created by davem on 16/02/2022.
//

#include "ilqrCore.h"

extern MujocoController *globalMujocoController;
extern m_state X_desired;
extern float dt;
extern int mujoco_steps_per_dt;

extern int initStateIndex;
extern int baseStateIndex;
extern float horizonLength;
extern int linearising_num_sim_steps;
extern bool alphaSearchEnabled;
extern float maxLamda;
extern float minLamda;
extern float lamdaFactor;
extern float epsConverge;
extern int numControls;
extern int torqueLims[NUM_DOF];

extern m_dof_dof R;
extern m_state_state Q;
extern m_state q;
extern m_dof r;

float lamb = 0.1;
int numIterations = 0;
float oldCost;

void differentiateDynamics(m_state *X, m_dof *U, m_state_state *f_x, m_state_dof *f_u, float *l, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu){

    MatrixXd I(NUM_STATES, NUM_STATES);
    I.setIdentity(NUM_STATES, NUM_STATES);
    MatrixXd A = ArrayXXd::Zero(NUM_STATES, NUM_STATES);
    MatrixXd B = ArrayXXd::Zero(NUM_STATES, NUM_DOF);
    // Linearise the dynamics
    for(int t = 0; t < numControls; t++){
        // Calculate linearised dynamics for current time step via finite differencing
//        globalMujocoController->setSystemState(X[t]);
//        globalMujocoController->deleteLastMujocoState();
//        globalMujocoController->saveMujocoState();
        //lineariseDynamics(X[t], U[t], A, B);
        lineariseDynamicsParallel(X[t], U[t], A, B);

        m_state_state A_dt;
        m_state_dof B_dt;
        int scaling = mujoco_steps_per_dt / linearising_num_sim_steps;
        scaleLinearisation(A, B, A_dt, B_dt,  scaling);

        f_x[t] = A_dt;
        f_u[t] = B_dt;

        l[t] = immediateCostAndDerivitives(l_x[t], l_xx[t], l_u[t], l_uu[t], X[t], U[t], t);

        l[t]    *= dt;
        l_x[t]  *= dt;
        l_xx[t] *= dt;
        l_u[t]  *= dt;
        l_uu[t] *= dt;

//        cout << "------------------- iteration: " << t << " --------------------" << endl;
//        cout << "l_xx[t]: " << l_xx[t] << endl;
//        cout << "l_x[t]: " << l_x[t] << endl;

        m_state _;
        stepSimulation(X[t], U[t], X[t+1], _, mujoco_steps_per_dt);

//        cout << "iteration " << t << endl;
//        cout << " f_x " << f_x[t] << endl;
//        cout << " f_u " << f_u[t] << endl;
//        int a = 1;
//
//        cout << "cost " << l[t] << endl;
//        cout << "cost dif x" << l_x[t] << endl;
//        cout << "cost dif xx" << l_xx[t] << endl;

    }

    //l[numControls] = terminalCost(l_x[numControls], l_xx[numControls], X[numControls]);

//    l   [numControls] *= dt;
//    l_x [numControls] *= dt;
//    l_xx[numControls] *= dt;
}

bool backwardsPass(m_state_state *A, m_state_dof *B, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K){
    m_state V_x;
    V_x = l_x[numControls - 1];
    m_state_state V_xx;
    V_xx = l_xx[numControls - 1];
    bool validBackwardsPass = true;

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
        int a = 1;

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
            cout << "B" << B<< endl;
            cout << "A" << A[t]<< endl;
            int a = 1;

        }
    }

    return validBackwardsPass;
}

bool backwardsPassTest(m_state_state *A, m_state_dof *B, float l, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K, m_state *X){
    float V = l;
    m_state V_x;
    V_x = l_x[numControls];;
    m_state_state V_xx;
    V_xx = l_xx[numControls];
    bool validBackwardsPass = true;

    for(int t = numControls - 1; t >= 1; t--){

        V_xx = (V_xx + V_xx.transpose()) / 2;

        m_state c = X[t] - X[t + 1];
        V_xx.diagonal().array() += lamb;
//        if(t > numControls - 5){
//            cout << "V_xx" << V_xx << endl;
//        }
        auto temp = (-2*B[t].transpose()*(V_xx)*B[t] - (2*R)).ldlt();

        K[t].noalias() = temp.solve(2*B[t].transpose()*(V_xx)*A[t]);

        m_dof temp1_k = B[t].transpose() * V_x;
        m_dof temp2_k = 2 * B[t].transpose() * V_xx * c;
        k[t].noalias() = temp.solve(temp1_k + temp2_k + r.transpose());

//        cout << "before" << endl;
//        cout << V_x << endl;
        //
        m_state temp1_1 = ((k[t].transpose()*B[t].transpose()) + c.transpose());

        m_state temp1 = 2*temp1_1.transpose()*V_xx*(A[t]+(B[t]*K[t]));
        m_state temp2 = (V_x.transpose())*(A[t]+(B[t]*K[t]));
        m_state temp3 = temp2 + q;
        m_state temp4 = (2*k[t].transpose()*R*K[t]);
//        cout << "temp1" << endl;
//        cout << temp1 << endl;
//        cout << "temp2" << endl;
//        cout << temp2 << endl;
//        cout << "temp3" << endl;
//        cout << temp3 << endl;
//        cout << "temp4" << endl;
//        cout << temp4 << endl;
        V_x = temp1 + temp3 + temp4;
//        cout << "New V_x" << endl;
//        cout << V_x << endl;

        //vn=2(kTnBT+cT)Vn−1(A+BKn)+vn−1(A+BKn)+q+2kTnRKn

        V_xx = (A[t]+B[t]*K[t]).transpose()*(V_xx)*(A[t]+B[t]*K[t])+Q+K[t].transpose()*R*K[t];

//        cout << "k[t] " << endl << k[t] << endl;
//        cout << "K[t]" << endl << K[t] << endl;
//        cout << "V_xx " << endl << V_xx << endl;
//        cout << "V_x" << endl << V_x << endl;
        int a = 1;
    }


    return true;
}

float forwardsPassTest(m_state *X, m_state *X_best, m_dof *U, m_dof *U_best, m_dof *k, m_dof_state *K){
    m_dof *U_new = new m_dof[numControls];
    m_state *X_new = new m_state[numControls + 1];
    float alpha = 1;
    float bestAlpha = alpha;
    float newCost = 0;
    float bestCost = 1000;
    int numAlphaChecks = 1;

    for(int i = 0; i < numAlphaChecks; i++){
        X_new[0] = X[0];
        newCost = 0;
        globalMujocoController->loadSimulationState(initStateIndex);
        m_state systemState = globalMujocoController->returnSystemState();
        cout << "system init state" << endl << systemState << endl;

        for(int t = 0; t < numControls; t++){
            m_state stateFeedback;
            stateFeedback = X_new[t] - X[t];
            m_dof feedBackGain = K[t] * stateFeedback;
            m_dof linearFeedback = (alpha * k[t]);

            // calculate new optimal control for this timestep
            U_new[t] = U[t] + (alpha * k[t]) + feedBackGain;

            cout << "old state was: " << endl << X[t] << endl;
            cout << "new state is: " << endl << X_new[t] << endl;

            cout << "feedback : " << endl << stateFeedback << endl;
            cout << "feedback gain: " << endl << feedBackGain << endl;
//
            cout << "old U was: " << endl << U[t] << endl;
            cout << "new U is: " << endl << U_new[t] << endl;

            // constrain new torque within limits
            for(int k = 0; k < NUM_DOF; k++){
                if(U_new[t](k) > torqueLims[k]) U_new[t](k) = torqueLims[k];
                if(U_new[t](k) < -torqueLims[k]) U_new[t](k) = -torqueLims[k];
            }

            // step simulation wit new control and repeat
            m_state _;
            stepSimulation(X_new[t], U_new[t], X_new[t+1], _, mujoco_steps_per_dt);

            // calc current state cost and keep running tally
            float currentCost = immediateCost(X_new[t], U_new[t], t);
            newCost += (currentCost * dt);

        }

        if(newCost < bestCost){
            bestCost = newCost;
            for(int j = 0; j < numControls; j++){
                X_best[j + 1] = X_new[j + 1];
                U_best[j] = U_new[j];
            }
            bestAlpha = alpha;
        }
        alpha += 0.01;
    }

    cout << "terminal state best: " << endl << X_best[numControls] << endl;
//    for(int t = 0; t < numControls; t++){
//        cout << "control " << t << " : " << U_best[t] << endl;
//    }
    cout << "best alpha was " << bestAlpha << endl;
    cout << "best cost was " << bestCost << endl;


    return bestCost;
}

float forwardsPass(m_state *X, m_state *X_new, m_dof *U, m_dof *U_new, m_dof *k, m_dof_state *K){
    auto forwardPassStart = high_resolution_clock::now();
    bool alphaLineSearch = true;
    float alpha = 1.0f;
    float bestAlphaCost = 0.0f;
    float newAlphaCost = 0.0f;
    bool firstAlpha = true;
    int alphaSearchCount = 0;
    float newCost;
    float initialAlphaCost = 0.0f;
    m_dof *U_best = new m_dof[numControls];
    m_state *X_best = new m_state[numControls + 1];
    int maxAlphaSearchDepth = 5;
    X_best[0] = X[0];

    while(alphaLineSearch){
        m_state xnew(NUM_STATES);
        xnew = X[0];
        globalMujocoController->loadSimulationState(initStateIndex);
        // Calculate new controls using optimal gain matrices
        for(int t = 0; t < numControls; t++){
            m_state stateFeedback;
            stateFeedback = xnew - X[t];
            m_dof feedBackGain = K[t] * stateFeedback;
            m_dof linearFeedback = (alpha * k[t]);

            U_new[t] = U[t] + (alpha * k[t]) + feedBackGain;
            //U_new[t] = (alpha * k[t]) + feedBackGain;

            if(t < numControls){
//                cout << "---------------------------------------" << endl;
//                cout << "x new " << xnew(0) << endl;
//                cout << "stateFeedback" << stateFeedback(0) << endl;
//                cout << "control number  " << t << endl;
//                cout << "feedbackgain " << feedBackGain(0) << endl;
//                cout << "linearFeedback " << linearFeedback(0) << endl;
//                cout << "old control" << U[t](0) << endl;
//                cout << "new control" << U_new[t](0) << endl;
//                cout << "little k " << endl << k[t] << endl;
//                cout << "Big k " << endl << K[t] << endl;
//                int a = 1;
            }


            // simulate new optimal control and calculate new predicted state.
            m_state _x = xnew;
            m_state _;
            for(int k = 0; k < NUM_DOF; k++){
                if(U_new[t](k) > torqueLims[k]) U_new[t](k) = torqueLims[k];
                if(U_new[t](k) < -torqueLims[k]) U_new[t](k) = -torqueLims[k];
            }
            stepSimulation(_x, U_new[t], xnew, _, mujoco_steps_per_dt);
        }

        cout << "new control index 1: " << U_new[1] << endl;
        newAlphaCost = rollOutTrajectory(X[0], X_new, U_new, numControls);


        if(firstAlpha){
            firstAlpha = false;

            bestAlphaCost = newAlphaCost;
            for(int i = 0; i < numControls; i++){
                U_best[i] = U_new[i];
                X_best[i] = X_new[i];
            }

            initialAlphaCost = newAlphaCost;
            if(!alphaSearchEnabled){
                break;
            }
            alpha /= 2;
            alphaSearchCount++;
        }
        else{
            float costGrad = abs(newAlphaCost - bestAlphaCost);
            if(costGrad < 0.1){
                alphaLineSearch = false;
            }
            if(alphaSearchCount > maxAlphaSearchDepth){
                alphaLineSearch = false;
            }
            if(newAlphaCost < bestAlphaCost){
                alpha = alpha - pow(0.5, alphaSearchCount + 1);
                bestAlphaCost = newAlphaCost;
                for(int i = 0; i < numControls; i++){
                    U_best[i] = U_new[i];
                    X_best[i] = X_new[i];
                }

            }
            else{
                alpha = alpha + pow(0.5, alphaSearchCount + 1);
            }

        }
        alphaSearchCount++;
    }

    newCost = bestAlphaCost;
    for(int i = 0; i < numControls; i++){
        U_new[i] = U_best[i];
        X_new[i] = X_best[i];
    }

    auto forwardPassStop = high_resolution_clock::now();
    auto forwardPassDur = duration_cast<microseconds>(forwardPassStop - forwardPassStart);
    float costImprovement = initialAlphaCost - newCost;

    cout << "Forwards pass: " << forwardPassDur.count()/1000 << " milliseconds" << "num iterations: " << alphaSearchCount << " inital cost: " << initialAlphaCost << ", reduced cost: " << newCost  << ", cost improvement: " << costImprovement << endl;
    return newCost;
}

bool checkForConvergence(float newCost, m_state *X, m_state *X_new, m_dof *U, m_dof *U_new, bool *costImprovement){
    bool convergence = false;
    if(newCost < oldCost){
        *costImprovement = true;
        if(lamb > minLamda){
            lamb /= lamdaFactor;
        }

        // replace old controls with new controls
        // replace old states with new states
        for(int k = 0; k < numControls; k++){
            U[k] = U_new[k];
//            cout << "new control timestep:  " << k << U[k] << endl;
            X[k + 1] = X_new[k + 1];
        }

        std::cout << "--------------------------------------------------" <<  std::endl;
        std::cout << "New cost: " << newCost <<  std::endl;
        m_state X_diff = X[numControls] - X_desired;
        std::cout << "terminal state is " << endl << X[numControls] << endl;

        numIterations++;
        float costGrad = (oldCost - newCost)/newCost;

        if((numIterations > 5) && costGrad < epsConverge){
            convergence = true;
            cout << "ilQR converged, num Iterations: " << numIterations << " final cost: " << newCost << endl;
        }
        oldCost = newCost;
    }
    else{
        // if cost improvement was not made, regulate lamda and run backwards pass again
        lamb *= lamdaFactor;
        if(lamb > maxLamda){
            convergence = true;
            cout << "ilQR finished due to lamda exceeds lamda max " << endl;
        }
    }

    return convergence;
}

void increaseLamda(){
    lamb *= lamdaFactor;
}

bool checkValidAandBMatrices(m_state_state A, m_state_dof B){
    bool validMatrices = true;
    int rowCounter = 0;
    int colCounter = 0;

    // Check the 'A' matrix
    while(rowCounter < NUM_STATES){
        while(colCounter < NUM_STATES){
            if(A(rowCounter, colCounter) > 5 || A(rowCounter, colCounter) < -5){
                validMatrices = false;
                break;
            }
            colCounter++;
        }
        colCounter = 0;
        rowCounter++;
    }

    colCounter = 0;
    rowCounter = 0;
    // Check the 'B' matrix
    while(rowCounter < NUM_STATES){
        while(colCounter < NUM_DOF){
            if(B(rowCounter, colCounter) > 5 || B(rowCounter, colCounter) < -5){
                validMatrices = false;
                break;
            }
            colCounter++;
        }
        colCounter = 0;
        rowCounter++;
    }


    return validMatrices;
}
