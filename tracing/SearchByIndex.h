//
// Created by Fengmei JIN on 2022/10/5.
//

#ifndef PCT_SEARCHBYINDEX_H
#define PCT_SEARCHBYINDEX_H

#include "../index/rtree/RtreeIndex.h"
#include "../index/intervaltree/ItreeIndex.h"
#include "../utils.h"
#include "functions.h"

struct Subseq {
    Point *startP{nullptr};
    Point *endP{nullptr};
    TIMETYPE contactStartTime{0};
    TIMETYPE contactDuration{0};
    double maxContactProb{0};

    Subseq() = default;

    Subseq(Point *sp, Point *ep, TIMETYPE cst, TIMETYPE cet, double cp) :
        startP(sp), endP(ep), contactStartTime(cst), contactDuration(cet - cst), maxContactProb(cp){}
};

/**
 * Construct index @param itree based on the @param candiTrajectories
 * */
char *initializeForItree(ItreeIndex *itree, IndexPoint *pIndex, const string &searchType,
                         TIMETYPE earliestIndexTime, TIMETYPE contactMinDuration, int deltaT,
                         const vector<pair<Trip *, ContactEvent *>> &candiTrajectories);

/* Filter returned intervals to get qualified trajectory pairs */
int getQualifiedPairs(const string &queryTrajID, bool mbrCompare,
                      IDTYPE queryPID, IndexPoint *pIndex, float distThreshold, const vector<Interval *> &returnedIntervals,
                      unordered_map<string, vector<tuple<Point *,TIMETYPE, TIMETYPE>>> &candiTrajPointPairs,
                      const unordered_set<string> &existingTrajs = unordered_set<string>());

/* Transfer from candidate trajectory pairs to subsequences */
void composeSubseq(Point *qtr, const unordered_map<string, vector<tuple<Point *,TIMETYPE, TIMETYPE>>> &candiTrajPointPairs,
                   unordered_map<string, vector<Subseq>> &traj2subseq, bool enableCP = false,
                   IndexPoint *pIndex = nullptr, float contactMaxDist = 0, bool independentCP = true);

/* query from susceptible version */
ContactEvent *checkConnectivity(const string& queryTrajID, Point *qtr, int k, TIMETYPE contactMinDuration,
                                unordered_map<string, vector<Subseq>> &traj2subseq, vector<TIMETYPE> &currWindowTime,
                                bool enableCP, unordered_map<string, ContactEvent *> &spreader2event);

/**
 * Given a @param susceptibleTraj as query, we search contacts on the @param itree
 * */
ContactEvent *searchDirectContact_tree(const Trip *susceptibleTraj, int k, ItreeIndex *itree, const string &searchType,
                                       bool ground_truth, IndexPoint *pIndex,
                                       float contactMaxDist, TIMETYPE contactMinDuration, float contactMinProb, bool independentCP,
                                       vector<TIMETYPE> &timeCosts);



#endif //PCT_SEARCHBYINDEX_H
