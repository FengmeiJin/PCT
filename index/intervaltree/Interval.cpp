//
// Created by Fengmei JIN on 2022/9/26.
//

#include "Interval.h"

bool Interval::existTemporalOverlap(TIMETYPE t1, TIMETYPE t2) const {
    return !(this->_startTime >= t2 || this->_endTime <= t1);
}

bool Interval::CompareInterval(Interval *i1, Interval *i2) {
    TIMETYPE t1 = i1->_startTime;
    TIMETYPE t2 = i2->_startTime;
    if (t1 == t2) {
        return i1->_endTime < i2->_endTime;
    }
    return t1 < t2;
}

bool Interval::existSpatialOverlap(GeoPoint *queryP, float distThreshold) const {
    if(this->_mbr != nullptr) {
        if(distThreshold > 0) {
            double mDist = this->_mbr->computeMinDist(queryP);
            return mDist <= distThreshold;
        }
        else {
            return this->_mbr->contain(queryP->getLongitude(), queryP->getLatitude());
        }
    }
    return true;
}

TrajPointPair Interval::getTrajPointPair() {
    return _tp;
}

IDTYPE Interval::getPointId() const {
    return _tp._ptr->getPointId();
}

float Interval::getPointRadius() const {
    return _tp._ptr->getRadius();
}

void Interval::releaseSpace() {
    if(_mbr != nullptr)
        _mbr->releaseSpace();
}
