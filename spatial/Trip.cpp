//
// Created by Fengmei Jin on 6/7/2022.
//

#include "Trip.h"
#include "../tracing/integral.h"

// used when computing distribution based on the surrounding pois
#define MAX_RADIUS 5

Trip::Trip(const string &_uid, const vector<tuple<IDTYPE, TIMETYPE, TIMETYPE>> &poiSeq, IndexPoint *pIndex,
           StaticGrid *poiGrid, bool needPoiDistribution) {
    uid = _uid;
    head_ptr = new Point;
    tail_ptr = new Point;
    length = 0;

    Point *point, *prevP = head_ptr;

    float minLng = 200, maxLng = -200;  // longitude: [-180, 180]
    float minLat = 100, maxLat = -100;  // latitude:  [-90, 90]

    const unsigned int poiCategoryNum = needPoiDistribution ? poiGrid->getCategoryNum() : 0;

    for (unsigned int i = 0, len = poiSeq.size(); i < len; i++) {
        auto tp = poiSeq.at(i);
        IDTYPE poiId = get<0>(tp);
        TIMETYPE arvTime = get<1>(tp);
        TIMETYPE levTime = get<2>(tp);

        float *poiDistribution = nullptr;
        if(needPoiDistribution) {
            poiDistribution = new float [poiCategoryNum];
            memset(poiDistribution, 0, poiCategoryNum);
            vector<pair<POI *, double>> list;
            list.emplace_back(make_pair(poiGrid->getPOI(poiId), 0));
            poiGrid->computeAggregatedDistribution(list, poiDistribution, false);
        }

        if (i > 0 && prevP->getPointId() == poiId) {
            if(!needPoiDistribution || prevP->equalPoiDistribution(poiDistribution, poiCategoryNum)) {
                // same as the previous point
                prevP->setLeaveTimestamp(levTime);  // update the previous point's timestamp (instead of creating new one)
                continue;
            }
        }

        GeoPoint *p = pIndex->getPointById(poiId);

        // update the spatial range
        minLng = MIN(minLng, p->getLongitude());
        maxLng = MAX(maxLng, p->getLongitude());
        minLat = MIN(minLat, p->getLatitude());
        maxLat = MAX(maxLat, p->getLatitude());

        // create a new point (not region, radius = 0)
        point = new Point(poiId, nullptr, arvTime, levTime, 0, poiDistribution);

        // set pointers
        prevP->setNextPointer(point);
        point->setPrevPointer(prevP);
        prevP = point;
        length++;
    }
    tail_ptr->setPrevPointer(point);
    point->setNextPointer(tail_ptr);     // this is the last point

    timeRange = make_pair(get<1>(poiSeq.front()), get<2>(poiSeq.back()));
    rec = new Rectangle(minLng, minLat, maxLng, maxLat);
}

pair<IDTYPE, float*> getPointID(float &pLng, float &pLat, IndexPoint *pIndex, StaticGrid *poiGrid = nullptr, float radius = 0,
                                unsigned int poiCategoryNum = 0, int alignRaw2POI = 1){
    IDTYPE currPid = 0;
    float *poiDistribution = nullptr;

    vector<pair<POI *, double>> sortedList;
    if(poiGrid != nullptr) {
        if (poiCategoryNum == 0 // don't need semantic distribution
                || alignRaw2POI == 1) {  // or align a raw point to a single POI
            POI *poi = poiGrid->getNearestPOI(pLng, pLat);
            if(poi != nullptr) {
                sortedList.emplace_back(make_pair(poi, -1));
            }
        }
        else {
            // further compute POI distribution for this point
            // one reason using surrounding distribution is that, one POI may be aligned by different raw points
            sortedList = poiGrid->getSurroundingPOIs(pLng, pLat, radius);     // sort by distance to the point
        }
        if(!sortedList.empty() && poiCategoryNum > 0) {
            poiDistribution = new float[poiCategoryNum];
            memset(poiDistribution, 0, poiCategoryNum);
            poiGrid->computeAggregatedDistribution(sortedList, poiDistribution);
        }
    }

    // assign ID for this point
    if(alignRaw2POI > 0 && !sortedList.empty()) {
        POI *mainPOI = sortedList.front().first;   // the closest POI
        currPid = mainPOI->poiID;
        pLng = mainPOI->point->getLongitude();   // update its coordinates
        pLat = mainPOI->point->getLatitude();
    }

    // finally, if still not available
    if (currPid == 0){
        currPid = pIndex->getVertexId(pLng, pLat);    // a new ID
    }
    return make_pair(currPid, poiDistribution);
}

Trip::Trip(const string &_uid, const vector<IndisRegion> &noisyRegions,
           IndexPoint *pIndex, StaticGrid *poiGrid, bool needPoiDistribution, int alignRaw2POI) {
    uid = _uid;
    head_ptr = new Point;
    tail_ptr = new Point;
    length = 0;

    Point *point, *prevP = head_ptr;
    Point *oriPtr;

    TIMETYPE arvTime, levTime;
    const unsigned int poiCategoryNum = needPoiDistribution ? poiGrid->getCategoryNum() : 0;

    for (unsigned int i = 0, len = noisyRegions.size(); i < len; i++) {
        IndisRegion region = noisyRegions.at(i);
        float radius = region.radius;
        oriPtr = region.oriPtr;
        arvTime = oriPtr->getArriveTimestamp();
        levTime = oriPtr->getLeaveTimestamp();
        if(levTime == 0)    // should not happen!!!
            levTime = (i < len - 1) ? noisyRegions.at(i + 1).oriPtr->getArriveTimestamp() : arvTime;

        // assign pid for each noisy point; assign poi distribution for the region
        pair<IDTYPE, float*> rtn = getPointID(region.lng, region.lat, pIndex, poiGrid, radius, poiCategoryNum);
        IDTYPE currPid = rtn.first;

        if (i > 0 && prevP->getPointId() == currPid) {
            if(!needPoiDistribution || prevP->equalPoiDistribution(rtn.second, poiCategoryNum)) {
                // same as the previous point
                prevP->setLeaveTimestamp(levTime);  // update the previous point's timestamp (instead of creating new one)
                continue;
            }
        }

        if(this->rec == nullptr) {
            this->rec = new Rectangle(pIndex->getPointById(currPid), radius);
        }
        else {
            Rectangle approxCircle(pIndex->getPointById(currPid), radius);
            this->rec->enlarge(&approxCircle);
        }

        point = new Point(currPid, oriPtr, arvTime, levTime, radius, rtn.second);

        // set pointers
        prevP->setNextPointer(point);
        point->setPrevPointer(prevP);
        prevP = point;
        length++;
    }
    tail_ptr->setPrevPointer(point);
    point->setNextPointer(tail_ptr);     // this is the last point

    timeRange = make_pair(noisyRegions.front().oriPtr->getArriveTimestamp(), levTime);
//    rec = new Rectangle(minLng, minLat, maxLng, maxLat);
}

Point *Trip::getFirstPoint() const {
    if (length > 0 && head_ptr != nullptr) {
        return head_ptr->getNextPointer();
    }
    return nullptr;
}

void Trip::clearAllPoints() {
    auto p = head_ptr;
    Point *toDelete;
    while (p != nullptr) {
        toDelete = p;
        p = p->getNextPointer();
        toDelete->clear();
        delete toDelete;
        toDelete = nullptr;
    }
}

bool Trip::spatialOverlap(const Trip *pTrip, float distThreshold) const {
    if (!rec->intersect(*pTrip->rec)) {     // two MBRs don't overlap
        double dist = this->rec->computeMinDist(pTrip->rec);
        return dist <= distThreshold;       // but if their minDist <= threshold, they can be contacted possibly
    }
    return true;
}

bool Trip::temporalOverlap(const Trip *pTrip) const {
    return !(timeRange.first > pTrip->timeRange.second || timeRange.second < pTrip->timeRange.first);
}

bool timeSpanOverlap(const Point *p, const Point *q, TIMETYPE timeThreshold = 0) {
    TIMETYPE pt1 = p->getArriveTimestamp(), pt2 = p->getLeaveTimestamp();
    TIMETYPE qt1 = q->getArriveTimestamp(), qt2 = q->getLeaveTimestamp();

    if (pt1 > qt2 || pt2 < qt1) {
        return false;
    } else {
        if (timeThreshold <= 0) {
            return true;
        } else {
            TIMETYPE overlap = 0;
            if (pt1 <= qt1) {
                overlap = (pt2 <= qt2 ? pt2 : qt2) - qt1;
            } else {
                overlap = (qt2 <= pt2 ? qt2 : pt2) - pt1;
            }
            return overlap >= timeThreshold;
        }
    }
}

/**
 * @param validSpreadStartTime: the start time of the @param spreaderTraj being spreadable
 * @param validSpreadEndTime: the end time of the one being spreader with regard to the current hop length ( = 0 means to the end)
 * @param earliestCST_itself: if > 0, this trajectory will become infectious after this moment, so no need to check hereafter
 * @param earliestCST: the earliest contact found so far, expecting a new contact earlier than this,
 *                      otherwise, no need to search further
 * */
tuple<TIMETYPE, int, double> Trip::contactWith(const Trip *spreaderTraj, TIMETYPE validSpreadStartTime, IndexPoint *pIndex,
                                               float distThreshold, TIMETYPE contactMinDuration, bool enableCP, bool independentCP,
                                               TIMETYPE earliestCST) const {
    TIMETYPE contactDurationSum = 0, contactStartTime = 0;
    int cnt = 1;    // at least one
    Point *q = this->getFirstPoint(), *p = spreaderTraj->getFirstPoint();
    while (q != nullptr && !q->isTailPtr() && p != nullptr && !p->isTailPtr()) {
        TIMETYPE pArriveT = p->getArriveTimestamp(), pLeaveT = p->getLeaveTimestamp();
        TIMETYPE qArriveT = q->getArriveTimestamp(), qLeaveT = q->getLeaveTimestamp();     // query

        // there exists a spreader that contacts with the current one earlier than pArriveT (or qArriveT)
        if(earliestCST > 0 && (pArriveT > earliestCST || qArriveT > earliestCST)) {
            break;
        }

        // the spreader isn't infectious until "validSpreadStartTime"
        if(validSpreadStartTime > 0 && pLeaveT <= validSpreadStartTime) {
            p = p->getNextPointer();    // before validSpreadStartTime, the spreaderTraj is not spreadable
            cnt++;
            continue;
        }

        pArriveT = MAX(pArriveT, validSpreadStartTime);
//        pLeaveT = validSpreadEndTime == 0 ? pLeaveT : MIN(pLeaveT, validSpreadEndTime);

        // no need to compare spatial proximity until two points have temporal overlapping
        if (qArriveT >= pLeaveT || pArriveT >= qLeaveT) {
            if (qArriveT >= pLeaveT) {
                p = p->getNextPointer();    // p delays, move ahead
            }
            else if (pArriveT >= qLeaveT) {
                q = q->getNextPointer();    // q delays, move ahead
            }
            contactDurationSum = 0;     // reset duration counter
            contactStartTime = 0;
            cnt++;
        }
        else {
            // exists time overlapping between p and q
            contactStartTime = contactStartTime == 0 ? MAX(pArriveT, qArriveT) : contactStartTime;
            double maxProb = 0;

            double dist_pq = pIndex->computeDistance(p->getPointId(), q->getPointId());
            double distMax = distThreshold + q->getRadius();    // maxR;
            if (dist_pq > distMax) {      // this pair of points cannot satisfy spatial threshold
                contactDurationSum = 0;
                contactStartTime = 0;   // reset
            }
            else {
                // compute contact probability
                if(enableCP){
                    double prob = computeContactProb(q, p, pIndex, distMax, independentCP);
                    maxProb = MAX(maxProb, prob);
                }
                contactDurationSum += MIN(pLeaveT, qLeaveT) - MAX(pArriveT, qArriveT);
                if (contactDurationSum >= contactMinDuration){ // && maxProb >= contactMinProb) {
                    return make_tuple(contactStartTime, cnt / 2, maxProb);
                }
            }

            if (qLeaveT <= pLeaveT) {
                q = q->getNextPointer();    // q moves ahead, p stays for next comparison
                cnt++;
            }
            if (qLeaveT >= pLeaveT) {
                p = p->getNextPointer();    // p moves ahead, q stays
                cnt++;
            }
        }
    }
    return make_tuple(0, cnt / 2, 0);
}

bool sortByStartTime(Trip *t1, Trip *t2) {
    TIMETYPE st1 = t1->getTimeRange().first, st2 = t2->getTimeRange().first;
    if(st1 == st2) {
        return t1->getTimeRange().second < t2->getTimeRange().second;
    }
    return st1 < st2;
}