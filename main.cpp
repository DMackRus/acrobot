#include "Utility/MujocoController/MujocoUI.h"
#include "iLQR/iLQR.h"
#include "iLQG.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include <Eigen/Core>
#include "mujoco.h"


extern MujocoController *globalMujocoController;
extern mjModel* model;						// MuJoCo model
extern mjData* mdata;						// MuJoCo data

inline mjtNum stepCost(const mjData* d)
{
    mjtNum cost;
    cost =
            1.0  * d->qpos[0] * d->qpos[0] +
            10.0 * d->qpos[1] * d->qpos[1] +
            1.0  * d->qvel[0] * d->qvel[0] +
            10.0 * d->qvel[1] * d->qvel[1] +
            1.0  * d->ctrl[0] * d->ctrl[0];
    return cost;
}

int main() {
    initMujoco();
    initialseController();

    //simpleTest();
    testILQR();



    render();

    return 0;
}

//const int degreesOfFreedom = 2;
//const int controls = 1;
//
//using x_t = Differentiator<degreesOfFreedom, controls>::x_t;
//using u_t = Differentiator<degreesOfFreedom, controls>::u_t;
//
//// bind data to eigen
//Eigen::Map<Eigen::VectorXd> qStar(mdata->qpos, model->nv);
//Eigen::Map<Eigen::VectorXd> qvelStar(mdata->qvel, model->nv);
//Eigen::Map<Eigen::VectorXd> ctrlStar(mdata->ctrl, model->nu);
//
//// change controls
//ctrlStar = ctrlStar.array() - 0.1;
//
//// save first state and contorls
//cout << "qStar" << qStar << endl;
//cout << "qvelStar" << qvelStar << endl;
//x_t x_star;
//cout << "x_star" << x_star << endl;
//x_star << qStar, qvelStar;
//u_t u_star; u_star << ctrlStar;
//
//// linearize around set state and contorl
//stepCostFn_t stepCostFn = stepCost;
//Differentiator<degreesOfFreedom, controls>* differentiator = new Differentiator<degreesOfFreedom, controls>(model, mdata, stepCostFn);
//differentiator->updateDerivatives();
//
//// copy data to d_prime
//mjData* d = mj_makeData(model);
//cpMjData(model, d, mdata);
//
//// advance simulation on dStar and save next state
//mj_step(model, mdata);
//x_t xStarNext; xStarNext << qStar, qvelStar;
//
//// bind data to eigen
//Eigen::Map<Eigen::VectorXd> q(d->qpos, model->nv);
//Eigen::Map<Eigen::VectorXd> qvel(d->qvel, model->nv);
//Eigen::Map<Eigen::VectorXd> ctrl(d->ctrl, model->nu);
//
//// perturb around x* and save x and u
//mjtNum eps = 1e-6;
//q = q.array() + eps;
//qvel = qvel.array() + eps;
//ctrl = ctrl.array() + eps;
//x_t x; x << q, qvel;
//u_t u; u << ctrl;
//
//// advance simulation on perturbed data
//mj_step(model, d);
//x_t xNext; xNext << q, qvel;
//u_t uNext; uNext << ctrl;
//
//cout << "Next A matrix is: " << endl << *(differentiator->A) << endl;
//cout << " next B matrix " << endl << *(differentiator->B) << endl;
//
//// compare
//x_t xNextPrediction = *(differentiator->A) * (x - x_star) + *(differentiator->B) * (u - u_star) + xStarNext;
//std::cout << "xNextPrediction - xNext:" << '\n';
//std::cout << (xNextPrediction - xNext).transpose() << '\n';
//std::cout << "--------------------" << '\n';
//std::cout << "xNext - xStarNext:" << '\n';
//std::cout << (xNext - xStarNext).transpose() << '\n';
//std::cout << "--------------------" << '\n';
//std::cout << "xNext:" << '\n';
//std::cout << xNext.transpose() << '\n';
