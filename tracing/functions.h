//
// Created by Fengmei JIN on 2023/2/2.
//

#ifndef PCT_FUNCTIONS_H
#define PCT_FUNCTIONS_H

#include "../spatial/Trip.h"

struct ContactEvent {
    string spreaderID;
    TIMETYPE contactStartTime{0};
    TIMETYPE contactEndTime{0};     // 0 represents to the end of the trip
    int kHop{-1};
    double probability{0};

    ContactEvent(string id, TIMETYPE cst, TIMETYPE cet = 0, int k = 1, double prob = 0):
            spreaderID(std::move(id)), contactStartTime(cst), contactEndTime(cet), kHop(k), probability(prob) {};

    TIMETYPE getStartTime() const {
        return contactStartTime;
    }

    TIMETYPE getEndTime() const {
        return contactEndTime;
    }

    void updateEndTime(TIMETYPE et) {
        contactEndTime = et;
    }
};

static bool sortByCST(pair<Trip *, ContactEvent *> &p1, pair<Trip *, ContactEvent *> &p2) {
    TIMETYPE t1 = p1.second->getStartTime(), t2 = p2.second->getStartTime();
    return t1 == t2 ? p1.second->getEndTime() < p2.second->getEndTime() : (t1 < t2);
}


static bool impossible(const Trip *qTraj, const Trip *pTraj, float contactMaxDist) {
    bool tOverlap = qTraj->temporalOverlap(pTraj);
    if(!tOverlap) {
        bool sOverlap = qTraj->spatialOverlap(pTraj, contactMaxDist);
        if(!sOverlap) {
            return true;
        }
    }
    return false;
}

/**
 * For each query trajectory, any candidate who has
 *      1) no time overlap
 *      2) the spatial distance between boundaries is larger than threshold
 * is impossible (no need to compare in detail)
 * */
static void preFiltering(const string &qid, const Trip *qTraj, float contactMaxDist,
                         const vector<pair<Trip *, ContactEvent *>> &candiTrajectories,
                         unordered_set<string> &impossibleTrajs){
    for(const auto &traj: candiTrajectories) {
        string tid = traj.first->getId();
        if(tid != qid && impossible(qTraj, traj.first, contactMaxDist)) {
            impossibleTrajs.insert(tid);
        }
    }
    impossibleTrajs.insert(qid);    // itself, no computation
}

static void preFiltering(const string &qid, const Trip *qTraj, float contactMaxDist,
                         const map<string, Trip *> &candiTrajectories, unordered_set<string> &impossibleTrajs){
    for(const auto &traj: candiTrajectories) {
        string tid = traj.first;
        if(tid != qid && impossible(qTraj, traj.second, contactMaxDist)) {
            impossibleTrajs.insert(tid);
        }
    }
    impossibleTrajs.insert(qid);    // itself, no computation
}

#endif //PCT_FUNCTIONS_H
