//
// Created by davem on 20/01/2022.
//
#include "MujocoController.h"

mjModel* model;						// MuJoCo model
mjData* mdata;						// MuJoCo data
mjvCamera cam;						// abstract camera
mjvOption opt;						// visualization options
mjvScene scn;						// abstract scene
mjrContext con;						// custom GPU context
mjvPerturb pert;
GLFWwindow* window;
MujocoController* globalMujocoController;



ofstream saveControlsFile;

std::string saveControlsFilename = "lastControls.csv";

MujocoController::MujocoController(mjModel* m, mjData* d){
    _model = m;
    _data = d;

    saveControlsFile.open(saveControlsFilename);
}

void MujocoController::step(){
    mj_step(_model, _data);
}

mjModel* MujocoController::returnModel() {
    return _model;
}

//mjData* MujocoController::returnCurrentModelData(){
//
//    return newData;
//}

void MujocoController::saveMujocoState(){
    mjData *newData = mj_makeData(_model);
    mj_copyData(newData, _model, _data);
    _mujocoStates.push_back(newData);
}

void MujocoController::deleteLastMujocoState(){
    mj_deleteData(_mujocoStates[1]);
    _mujocoStates.pop_back();
}

void MujocoController::setSystemState(const Ref<const m_state> systemState){


    _data->qpos[0] = systemState[0];
    _data->qpos[1] = systemState[1];
    _data->qvel[0] = systemState[2];
    _data->qvel[0] = systemState[3];
}

m_state MujocoController::returnSystemState(){
    m_state systemState;

    systemState[0] = _data->qpos[0];
    systemState[1] = _data->qpos[1];
    systemState[2] = _data->qvel[0];
    systemState[3] = _data->qvel[1];

    return systemState;
}

void MujocoController::setBodyState(int bodyId, const Ref<const m_pose> pose, const Ref<const m_pose> velocities){

    // TODO extend this to work with rotations also, whether its quaternions or euler angles


    setBodyPose(bodyId, pose);
    setBodyVel(bodyId, velocities);

}

void MujocoController::setBodyPose(int bodyId, const Ref<const m_pose> pose){
    for(int i = 0; i < 3; i++){
        _data->qpos[((bodyId-1) * 7) + i] = pose(i);
    }
}

void MujocoController::setBodyVel(int  bodyId, const Ref<const m_pose> vel){
    for(int i = 0; i < 3; i++){
        _data->qvel[((bodyId - 1) * 6) + i] = vel(i);
    }
}

void MujocoController::setBodyAccelerations(int bodyId, const Ref<const m_pose> acc){
    for(int i = 0; i < 3; i++){
        _data->qacc[((bodyId) * 6) + i] = acc(i);
    }
}

void MujocoController::setBodyForces(int bodyId, const Ref<const m_pose> frc){
    for(int i = 0; i < 3; i++){
        _data->xfrc_applied[((bodyId) * 6) + i] = frc(i);
    }
}

m_pose MujocoController::returnBodyState(int bodyId){
    m_pose bodyPose;

    // TODO extend this to work with rotations also, whether its quaternions or euler angles
    for(int i = 0; i < 3; i++){
        bodyPose(i) = _data->qpos[((bodyId - 1) * 7) + i];
    }

    return bodyPose;
}

m_pose MujocoController::returnBodyVelocities(int bodyId){
    m_pose bodyVels;

    for(int i = 0; i < 3; i++){
        bodyVels(i) = _data->qvel[((bodyId - 1) * 6) + i];
    }

    return bodyVels;
}

m_pose MujocoController::returnBodyAcceleration(int bodyId){
    m_pose bodyAcc;

    for(int i = 0; i < 3; i++){
        bodyAcc(i) = _data->qacc[((bodyId - 1) * 6) + i];
    }

    return bodyAcc;
}

m_pose MujocoController::returnBodyForces(int bodyId){
    m_pose bodyFrc;

    for(int i = 0; i < 3; i++){
        bodyFrc(i) = _data->xfrc_applied[((bodyId) * 6) + i];
    }

    return bodyFrc;

}

void MujocoController::setRobotConfiguration(const Ref<const m_dof> configuration) {

    for (int i = 0; i < NUM_JOINTS; i++) {
        _data->qpos[i] = configuration(i);
    }
    mj_forward(model, _data);
}

m_dof MujocoController::returnRobotConfiguration(){
    m_dof robotConfig;

    for(int i = 0; i < NUM_JOINTS; i++){
        robotConfig(i) = _data->qpos[i];
    }
    return robotConfig;
}

void MujocoController::setRobotVelocities(const Ref<const m_dof> jointVelocities){
    for (int i = 0; i < NUM_JOINTS; i++) {
        _data->qvel[i] = jointVelocities(i);
    }
}

m_dof MujocoController::returnRobotVelocities(){
    m_dof robotVelocities;
    for(int i = 0; i < NUM_JOINTS; i++){
        robotVelocities(i) = _data->qvel[i];
    }
    return robotVelocities;
}

void MujocoController::setRobotAccelerations(const Ref<const m_dof> jointAccelerations){
    for (int i = 0; i < NUM_JOINTS; i++) {
        _data->qacc[i] = jointAccelerations(i);
    }
}

m_dof MujocoController::returnRobotAccelerations(){
    m_dof jointAccelerations;
    for(int i = 0; i < NUM_JOINTS; i++){
        jointAccelerations(i) = _data->qacc[i];
    }


    return jointAccelerations;
}

bool MujocoController::isConfigInCollision(m_dof configuration) {
    bool collision = false;
    int originalResetValue = 0;

    setRobotConfiguration(configuration);
    mj_step1(_model, _data);

    int numberOfCollisions = getRobotNumberOfCollisions();

    mj_resetData(_model, _data);
    if (_resetLevel > 0) {
        mjData* targetData = _mujocoStates.at(_resetLevel - 1);
        mj_copyData(_data, _model, targetData);
    }

    //Appears not to need this
    //mj_forward(m, d);

    if (numberOfCollisions > 0) {
        collision = true;
    }

    return collision;
}

int MujocoController::getRobotNumberOfCollisions() {

    int numContacts = _data->ncon;
    int numCollisions = 0;
    for (int i = 0; i < numContacts; i++) {
        auto contact = _data->contact[i];

        // Get the ids of the two bodies in contacts
        int bodyInContact1 = _model->body_rootid[_model->geom_bodyid[contact.geom1]];
        int bodyInContact2 = _model->body_rootid[_model->geom_bodyid[contact.geom2]];

        // only consider it a collision if robot - robot
        // or robot - table

        bool contact1Robot = false;
        bool contact1Table = false;
        bool contact2Robot = false;
        bool contact2Table = false;
        for (int j = 0; j < 11; j++) {
            if (bodyInContact1 == robotBodyID[j]) {
                contact1Robot = true;
            }

            if (bodyInContact2 == robotBodyID[j]) {
                contact2Robot = true;
            }
        }

//        if (bodyInContact1 == tableId) {
//            contact1Table = true;
//        }
//        if (bodyInContact2 == tableId) {
//            contact2Table = true;
//        }

        if (contact1Robot) {
            if (contact2Robot || contact2Table) {
                numCollisions++;
            }
        }
        else if(contact2Robot) {
            if (contact1Robot || contact1Table) {
                numCollisions++;
            }
        }
    }

    return numCollisions;
}

struct pose MujocoController::returnEndEffectorPos(){
    pose endEffectorPose;

    // FIX add code to assign rotation to endeffector pose
    endEffectorPose.pos.x = _data->xpos[30];
    endEffectorPose.pos.y = _data->xpos[31];
    endEffectorPose.pos.z = _data->xpos[32];

    return endEffectorPose;
}

void MujocoController::saveSimulationState(){
    mjData *newData = mj_makeData(_model);
    mj_copyData(newData, _model, _data);
    _mujocoStates.push_back(newData);

}

void MujocoController::setResetLevel(int resetLevel){
    _resetLevel = resetLevel;
}

void MujocoController::resetSimulation(){
    mj_resetData(_model, _data);

    if (_resetLevel > 0) {
        mjData *targetData = _mujocoStates.at(_resetLevel);
        mj_copyData(_data, _model, targetData);
    }

    mj_forward(_model, _data);
}

void MujocoController::loadSimulationState(int stateIndex){

    if(0){
        mj_copyData(_data, _model, _mujocoStates[stateIndex]);
    }
    else{
        _data->time = _mujocoStates[stateIndex]->time;
        mju_copy(_data->qpos, _mujocoStates[stateIndex]->qpos, _model->nq);
        mju_copy(_data->qvel, _mujocoStates[stateIndex]->qvel, _model->nv);
        mju_copy(_data->act,  _mujocoStates[stateIndex]->act,  _model->na);

        // copy mocap body pose and userdata
        mju_copy(_data->mocap_pos,  _mujocoStates[stateIndex]->mocap_pos,  3*_model->nmocap);
        mju_copy(_data->mocap_quat, _mujocoStates[stateIndex]->mocap_quat, 4*_model->nmocap);
        mju_copy(_data->userdata, _mujocoStates[stateIndex]->userdata, _model->nuserdata);

        // copy warm-start acceleration
        mju_copy(_data->qacc_warmstart, _mujocoStates[stateIndex]->qacc_warmstart, _model->nv);
    }
}

Eigen::MatrixXd MujocoController::calculateJacobian(int bodyId){
    Eigen::MatrixXd kinematicJacobian(3, 7);

    //mjtNum* J_COMi_temp = mj_stackAlloc(_data, 3*_model->nv);
    Matrix<double, Dynamic, Dynamic, RowMajor> J_p(3, _model->nv);
    Matrix<double, Dynamic, Dynamic, RowMajor> J_r(3, _model->nv);

    mj_jacBody(_model, _data, J_p.data(), NULL, bodyId);

    //J_COMi = Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(J_p, 3, _model->nv);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 7; j++) {
            kinematicJacobian(i, j) = J_p(i, j);
            //cout << kinematicJacobian(i, j) << endl;
        }
    }

    return kinematicJacobian;
}

int MujocoController::returnModelID(const std::string& input){
    return(mj_name2id(model, mjOBJ_BODY, input.c_str()));
}