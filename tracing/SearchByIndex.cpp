//
// Created by Fengmei JIN on 2022/10/5.
//

#include "SearchByIndex.h"
#include "integral.h"

void prepareTrajPointPair(Point *p, const string &tid, IndexPoint *pIndex, bool needMBR,
                          TIMETYPE validStartT, TIMETYPE validEndT,
                          map<TIMETYPE, map<TIMETYPE, pair<vector<TrajPointPair>, Rectangle *>>> &st2et2trajPair) {
    TIMETYPE st = p->getArriveTimestamp();
    if (st < validStartT) {
        st = validStartT;
    }
    TIMETYPE et = p->getLeaveTimestamp();
    if (et > validEndT) {
        et = validEndT;
    }
    auto itr = st2et2trajPair.find(st);
    if (itr != st2et2trajPair.end()) {   // the start time exists
        auto itr2 = itr->second.find(et);
        if (itr2 != itr->second.end()) {     // the end time exists
            itr2->second.first.emplace_back(TrajPointPair(tid, p));

            if (needMBR && itr2->second.second == nullptr) {
                itr2->second.second = new Rectangle(pIndex->getPointById(p->getPointId()), p->getRadius());
            }
        }
        else {
            Rectangle *mbr = nullptr;
            vector<TrajPointPair> trajPairV;
            trajPairV.emplace_back(TrajPointPair(tid, p));
            if (needMBR) {
                mbr = new Rectangle(pIndex->getPointById(p->getPointId()), p->getRadius());
            }
            itr->second.insert(make_pair(et, make_pair(trajPairV, mbr)));
        }
    }
    else {
        map<TIMETYPE, pair<vector<TrajPointPair>, Rectangle *>> et2trajPair;
        Rectangle *mbr = nullptr;
        vector<TrajPointPair> trajPairV;
        trajPairV.emplace_back(TrajPointPair(tid, p));

        if (needMBR) {
            mbr = new Rectangle(pIndex->getPointById(p->getPointId()), p->getRadius());
        }
        et2trajPair[et] = make_pair(trajPairV, mbr);
        st2et2trajPair[st] = et2trajPair;
    }
}

unsigned int prepareIntervals(vector<Interval *> &sortedIntervals,
                              const map<TIMETYPE, map<TIMETYPE, pair<vector<TrajPointPair>, Rectangle *>>> &st2et2trajPair) {
    // each distinct [st, et] composes an interval
    unsigned int avgPointNum = 0;
    for (const auto &stPair: st2et2trajPair) {
        TIMETYPE st = stPair.first;

        for (auto &etPair: stPair.second) {
            TIMETYPE et = etPair.first;
            vector<TrajPointPair> t2p = etPair.second.first;
            avgPointNum += t2p.size();

            for (const TrajPointPair &tp: t2p) {
                auto *iv = new Interval(st, et, tp, etPair.second.second);
                sortedIntervals.emplace_back(iv);   // should be sorted in time dimension, if push in this way
            }
        }
    }

    // !!! st2et2trajPair is naturally sorted by map, then no need to sort again
    return avgPointNum;
}

char *initializeForItree(ItreeIndex *itree, IndexPoint *pIndex, const string &searchType,
                         TIMETYPE earliestIndexTime, TIMETYPE contactMinDuration, int deltaT,
                         const vector<pair<Trip *, ContactEvent *>> &candiTrajectories) {
    unsigned int trajPointNum = 0, distinctStNum = 0;
    vector<Interval *> sortedIntervals;
    TIMETYPE LatestEndTime = 0, prepareTime = 0;
    TIMETYPE timer = millisecond();
    {
        map<TIMETYPE, map<TIMETYPE, pair<vector<TrajPointPair>, Rectangle *>>> st2et2trajPair;    // start time -> end time -> traj pairs
        // each map has been sorted by the key value

        bool needMBR = (searchType == "interval");
        for (const auto &tripPair: candiTrajectories) {
            Trip *trip = tripPair.first;
            ContactEvent *event = tripPair.second;
            trajPointNum += trip->getLength();
            string tid = trip->getId();
            Point *ptr = trip->getFirstPoint();
            TIMETYPE validStartT = event == nullptr ? trip->getTimeRange().first : event->getStartTime();
            TIMETYPE validEndT = event == nullptr ? trip->getTimeRange().second : event->getEndTime();
            LatestEndTime = MAX(LatestEndTime, validEndT);
            // scan all points of the current traj
            while (ptr != nullptr && !ptr->isTailPtr()) {
                if (!(ptr->getLeaveTimestamp() <= validStartT || ptr->getArriveTimestamp() >= validEndT)) {
                    prepareTrajPointPair(ptr, tid, pIndex, needMBR, validStartT, validEndT, st2et2trajPair);
                }
                ptr = ptr->getNextPointer();
            }
        }

        // avg point number for each distinct interval
        prepareIntervals(sortedIntervals, st2et2trajPair);    // should be false for TS-index
        prepareTime = millisecond() - timer;

        distinctStNum = st2et2trajPair.size();
        st2et2trajPair.clear();
    }

    // index construction
    timer = millisecond();

    itree->constructEmptyTree(earliestIndexTime, LatestEndTime, contactMinDuration, deltaT);

    itree->constructIndexByIntervals(sortedIntervals, searchType == "interval" ? nullptr : pIndex);

    TIMETYPE constructionTime = millisecond() - timer;     // unit: ms

    char *buffer = new char[200];
    // format: preparation, construction, total, tree height, interval number, distinct start time
    snprintf(buffer, 200, "%lu,%lu,%lu,%d,%lu,%d\n",
            prepareTime, constructionTime, prepareTime + constructionTime, itree->getHeight(),
            sortedIntervals.size(), distinctStNum);

    sortedIntervals.clear();
    sortedIntervals.shrink_to_fit();

    return buffer;
}

bool CompareTupleByTime(tuple<Point *, TIMETYPE, TIMETYPE> &p1, tuple<Point *, TIMETYPE, TIMETYPE> &p2) {
    TIMETYPE t1 = get<1>(p1), t2 = get<1>(p2);
    if (t1 == t2)
        return get<2>(p1) < get<2>(p2);
    return t1 < t2;
}

void
composeSubseq(Point *qtr, const unordered_map<string, vector<tuple<Point *, TIMETYPE, TIMETYPE>>> &candiTrajPointPairs,
              unordered_map<string, vector<Subseq>> &traj2subseq, bool enableCP,
              IndexPoint *pIndex, float contactMaxDist, bool independentCP) {
    TIMETYPE qArriveT = qtr->getArriveTimestamp(), qLeaveT = qtr->getLeaveTimestamp();

    // Step-2: Compose subseqences based on the returned points for each candi trajectory
    for (const auto &t2p: candiTrajPointPairs) {
        string tid = t2p.first;
        vector<tuple<Point *, TIMETYPE, TIMETYPE>> pList = t2p.second;
        if (pList.size() > 1) {
            sort(pList.begin(), pList.end(), CompareTupleByTime);   // sort the current candidate points by time
        }
        set<Point *> test;
        bool firstP = true;
        auto itr = traj2subseq.find(tid);
        for (auto tp: pList) {
            Point *p = get<0>(tp);

            TIMETYPE intervalStart = get<1>(tp);
            auto rtn1 = test.insert(p);

            bool appended = false;
            if (itr !=
                traj2subseq.end()) {  // already exist lastSeq, find a proper position for the current point if possible

                // only need to check the last subseq (closest to the current time)
                Subseq lastSeq = itr->second.back();
                if (lastSeq.endP ==
                    p) {    // the end point of the lastSeq covers the current timestamp, expand the contact length
                    lastSeq.contactDuration += (MIN(p->getLeaveTimestamp(), qLeaveT) - MAX(qArriveT, intervalStart));
                    if (enableCP) {
                        lastSeq.maxContactProb = MAX(lastSeq.maxContactProb,
                                                     computeContactProb(qtr, p, pIndex, contactMaxDist, independentCP));
                    }
                    appended = true;
                }
                else if (lastSeq.endP->getNextPointer() == p) {     // the lastSeq could keep continuous by appending p
                    if (firstP && p->getLeaveTimestamp() <= qArriveT || !firstP) {
                        lastSeq.endP = p;
                        lastSeq.contactDuration += (MIN(p->getLeaveTimestamp(), qLeaveT) -
                                                    MAX(qArriveT, intervalStart));// expand the contact length
                        if (enableCP) {
                            lastSeq.maxContactProb = MAX(lastSeq.maxContactProb,
                                                         computeContactProb(qtr, p, pIndex, contactMaxDist,
                                                                            independentCP));
                        }
                        appended = true;
                    }
                    // Special case: the existing lastSeq covers the previous time window
                    // but the previous end point isn't appended in the candidate intervals discovered at the current timestamp
                    // in this case, instead of appending to the end, this point should create a new lastSeq
                }
                itr->second[itr->second.size() - 1] = lastSeq;
            }
            else {
                auto rtn = traj2subseq.insert(make_pair(tid, vector<Subseq>()));
                itr = rtn.first;
            }
            if (!appended) {
                // create a new subseq
                itr->second.emplace_back(
                        Subseq(p, p, MAX(qArriveT, intervalStart), MIN(qLeaveT, p->getLeaveTimestamp()),
                               enableCP ? computeContactProb(qtr, p, pIndex, contactMaxDist, independentCP) : 0));
            }
            if (firstP) {
                firstP = false;
            }
        }
    }
}

int
getQualifiedPairs(const string &queryTrajID, bool mbrCompare, IDTYPE queryPID, IndexPoint *pIndex, float distThreshold,
                  const vector<Interval *> &returnedIntervals,
                  unordered_map<string, vector<tuple<Point *, TIMETYPE, TIMETYPE>>> &candiTrajPointPairs,
                  const unordered_set<string> &existingTrajs) {
    int intervalFilteredByMBR = 0;
    GeoPoint *queryP = pIndex->getPointById(queryPID);
    for (Interval *interval: returnedIntervals) {   // these trajectory points share the same GeoPoint (aka dist2q)
        TrajPointPair tp = interval->getTrajPointPair();
        if (mbrCompare && !interval->existSpatialOverlap(queryP, distThreshold)) {
            intervalFilteredByMBR++;
            continue;
        }

        string tid = tp._trajID;
        if (tid != queryTrajID &&
            existingTrajs.find(tid) == existingTrajs.end()) {  // no further check for the existing contacts
            // find a qualified id, further check distance threshold
            bool qualified = true;
            if (mbrCompare) {
                double d = pIndex->computeDistance(queryPID, tp._ptr->getPointId());
                double threshold = distThreshold + tp._ptr->getRadius();        // contactMaxDist + r1 + r2
                if (d > threshold) {
                    qualified = false;
                }
            }

            if (qualified) {
                auto itr = candiTrajPointPairs.find(tid);
                if (itr != candiTrajPointPairs.end()) {
                    itr->second.emplace_back(make_tuple(tp._ptr, interval->getStartTime(), interval->getEndTime()));
                }
                else {
                    vector<tuple<Point *, TIMETYPE, TIMETYPE>> points;
                    points.emplace_back(make_tuple(tp._ptr, interval->getStartTime(), interval->getEndTime()));
                    candiTrajPointPairs[tid] = points;
                }
            }
        }
    }
    return intervalFilteredByMBR;
}

ContactEvent *checkConnectivity(const string &queryTrajID, Point *qtr, int k, TIMETYPE contactMinDuration,
                                unordered_map<string, vector<Subseq>> &traj2subseq, vector<TIMETYPE> &currWindowTime,
                                bool enableCP, unordered_map<string, ContactEvent *> &spreader2event) {
    TIMETYPE qArriveT = qtr->getArriveTimestamp(), qLeaveT = qtr->getLeaveTimestamp();
    ContactEvent *earliestEvent = nullptr;

    // Step-3: If the time window exceeds the threshold, check all candidate subsequences
    currWindowTime.emplace_back(qArriveT);
    TIMETYPE currWindowStartTime = currWindowTime.front();   // the earliest "qArriveT" in the candidate traj list
    if (qLeaveT - currWindowStartTime >= contactMinDuration) {
        // process current time window to see if exist contacts satisfying the continuous constraints
        // by scanning each timestamp's candidates (both spatial and temporal)

        unsigned int windowLength = currWindowTime.size();
        bool *existValid = new bool[windowLength];
        memset(existValid, false, windowLength);

        TIMETYPE earliestCST_curr = 0;
        string earliestSpreader;
        double earliestProb = 0;

        auto candiTraj = traj2subseq.begin();
        while (candiTraj != traj2subseq.end()) {
            string tid = candiTraj->first;

            // scan its all subseqs to see if anyone satisfies the temporal threshold
            for (const auto &subseq: candiTraj->second) {
                if (subseq.contactDuration >= contactMinDuration) {
                    if (earliestCST_curr == 0 || subseq.contactStartTime < earliestCST_curr ||
                        (subseq.contactStartTime == earliestCST_curr && subseq.maxContactProb > earliestProb)) {
                        earliestCST_curr = subseq.contactStartTime;
                        earliestSpreader = tid;
                        earliestProb = subseq.maxContactProb;
                    }
                    if (enableCP) {
                        auto itr = spreader2event.find(tid);
                        if (itr != spreader2event.end() && subseq.contactStartTime < itr->second->getStartTime()) {
                            itr->second->contactStartTime = subseq.contactStartTime;
                            itr->second->spreaderID = tid;
                            itr->second->probability = subseq.maxContactProb;
                        }
                        else {
                            spreader2event[tid] = new ContactEvent(tid, subseq.contactStartTime, 0, k + 1,
                                                                   subseq.maxContactProb);
                        }
                    }

                    break;  // end of check the current spreader
                }
            }
            candiTraj++;
        }

        if (!earliestSpreader.empty()) {
            earliestEvent = new ContactEvent(earliestSpreader, earliestCST_curr, 0, k + 1, earliestProb);
            return earliestEvent;
        }
        else {
            candiTraj = traj2subseq.begin();
            while (candiTraj != traj2subseq.end()) {
                // only reserve the last subseq for next timestamp (others don't keep continuous)
                auto subseq = candiTraj->second.back();
                candiTraj->second.clear();
                candiTraj->second.emplace_back(subseq);

                unsigned int j = 0;
                TIMETYPE t = subseq.contactStartTime;
                while (j < windowLength && !existValid[j]) {
                    if (j < windowLength - 1 && currWindowTime[j] <= t && t <= currWindowTime[j + 1] ||
                        t >= currWindowTime[j]) {
                        existValid[j] = true;
                        break;
                    }
                    j++;
                }
                candiTraj++;
            }
        }

        unsigned int j = 0;
        // must remove the first timestamp, the following depends on the validability
        while (j < windowLength && !existValid[j]) {
            // pop the earliest one from the current window
            currWindowTime.erase(currWindowTime.begin());
            j++;
        }
    }

    return earliestEvent;
}

ContactEvent *searchDirectContact_tree(const Trip *susceptibleTraj, int k, ItreeIndex *itree, const string &searchType,
                                       bool ground_truth, IndexPoint *pIndex,
                                       float contactMaxDist, TIMETYPE contactMinDuration, float contactMinProb,
                                       bool independentCP,
                                       vector<TIMETYPE> &timeCosts) {
    Point *qtr = susceptibleTraj->getFirstPoint();
    if (qtr == nullptr || itree == nullptr) {
        printf("Error in SearchByIndex::searchDirectContact_tree: a null pointer!\n");
        return nullptr;
    }
    const string queryTrajID = susceptibleTraj->getId();

    vector<Interval *> returnedIntervals;
    vector<TIMETYPE> currWindowTime;
    unordered_map<string, vector<tuple<Point *, TIMETYPE, TIMETYPE>>> candiTrajPointPairs;       // associated with qualified trajectories
    unordered_map<string, vector<Subseq>> traj2subseq;

    TIMETYPE timer, intervalTreeSearchTime = 0, rtreeSearchTime = 0;
    TIMETYPE filterPairTime = 0, transferSeqTime = 0, connectivityTime = 0;
    TIMETYPE qTrajEndTime = susceptibleTraj->getTimeRange().second;
    unsigned long avgReturnIntervalNum = 0, avgFilterIntervalNum = 0, avgCandiPairNum = 0;
    unsigned long validPointNum = 0, checkedPointNum = 0;

    bool mbrCompare = (searchType == "interval");
    bool enableCP = (!ground_truth && contactMinProb > 0);

    ContactEvent *earliestEvent = nullptr;
    unordered_map<string, ContactEvent *> spreader2event;

    while (qtr != nullptr && !qtr->isTailPtr()) {
        TIMETYPE qArriveT = qtr->getArriveTimestamp();

        TIMETYPE possibleContactStartTime = !currWindowTime.empty() ? currWindowTime.front() : qArriveT;
        if (qTrajEndTime - possibleContactStartTime < contactMinDuration) {
            break;  // the following will not exist possible contact for the given query, terminate
        }

        // 1) find intervals based on temporal overlapping
        returnedIntervals.clear();
        timer = millisecond();
        if (searchType == "interval") {
            itree->overlapSearch(qtr, returnedIntervals);
        }
        else if (searchType == "IR") {  // interval -> Rtree
            rtreeSearchTime += itree->overlapSearch(qtr, returnedIntervals, pIndex, contactMaxDist);
        }
        intervalTreeSearchTime += (millisecond() - timer);
        avgReturnIntervalNum += returnedIntervals.size();

        candiTrajPointPairs.clear();
        if (!returnedIntervals.empty()) {
            timer = millisecond();
            float distThreshold = contactMaxDist + qtr->getRadius();
            avgFilterIntervalNum += getQualifiedPairs(queryTrajID, mbrCompare,
                                                      qtr->getPointId(), pIndex, distThreshold,
                                                      returnedIntervals, candiTrajPointPairs);
            filterPairTime += (millisecond() - timer);
        }

        // will check connectivity for the current timestamp if possible
        if (candiTrajPointPairs.empty()) {
            currWindowTime.clear();
            traj2subseq.clear();        // the continuous interrupt at this moment, reset all
        }
        else {
            avgCandiPairNum += candiTrajPointPairs.size();
            validPointNum++;
            timer = millisecond();
            if (enableCP) {
                composeSubseq(qtr, candiTrajPointPairs, traj2subseq, enableCP, pIndex, contactMaxDist, independentCP);
            }
            else {
                composeSubseq(qtr, candiTrajPointPairs, traj2subseq);
            }
            transferSeqTime += (millisecond() - timer);

            timer = millisecond();
            if (enableCP) {
                ContactEvent *rtn = checkConnectivity(queryTrajID, qtr, k, contactMinDuration, traj2subseq,
                                                       currWindowTime, enableCP, spreader2event);
                if (rtn != nullptr && earliestEvent == nullptr) {
                    earliestEvent = rtn;
                }
            }
            else {
                earliestEvent = checkConnectivity(queryTrajID, qtr, k, contactMinDuration, traj2subseq, currWindowTime,
                                                   enableCP, spreader2event);
            }
            connectivityTime += (millisecond() - timer);

            if (earliestEvent != nullptr) {
                if (!enableCP || (enableCP && earliestEvent->probability >= contactMinProb))
                    break;
            }
        }

        checkedPointNum++;

        qtr = qtr->getNextPointer();
    }

    if (enableCP && earliestEvent != nullptr) {
        double aggregatedNonContactProb = 1;
        double maxContactProb = 0;
        for (const auto &sp: spreader2event) {
            aggregatedNonContactProb *= (1 - (1.0 / (k + 1)) * sp.second->probability);
            maxContactProb = MAX(sp.second->probability, maxContactProb);
        }
        double aggregatedContactProb = 1 - aggregatedNonContactProb;
        if (aggregatedContactProb >= contactMinProb || maxContactProb >= contactMinProb) {
            earliestEvent->probability = aggregatedContactProb;
        }
        else {
            earliestEvent = nullptr;
        }
    }

    if (timeCosts.empty()) {
        int len = 5;
        for (int i = 0; i < len; i++)
            timeCosts.emplace_back(0);  // initialize
    }
    timeCosts[0] += intervalTreeSearchTime - rtreeSearchTime;
    timeCosts[1] += rtreeSearchTime;
    timeCosts[2] += filterPairTime;
    timeCosts[3] += transferSeqTime;
    timeCosts[4] += connectivityTime;

    return earliestEvent;
}