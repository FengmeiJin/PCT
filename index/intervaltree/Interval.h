//
// Created by Fengmei JIN on 2022/9/26.
//

#ifndef PCT_INTERVAL_H
#define PCT_INTERVAL_H

#include <utility>

#include "../../spatial/Point.h"
#include "../../spatial/IndexPoint.h"

struct TrajPointPair{
    string _trajID;
    Point *_ptr;     // for a distinct interval [st, et], each trajectory can only have one point at most

    TrajPointPair(string tid, Point *p): _trajID(std::move(tid)), _ptr(p) {}
};

// Structure to represent an interval [low, high]
class Interval {

private:
    TIMETYPE _startTime;
    TIMETYPE _endTime;
    TrajPointPair _tp;

    Rectangle *_mbr{nullptr};

public:

    // Constructor
    Interval(TIMETYPE st, TIMETYPE et, TrajPointPair tp, Rectangle *mbr):
        _startTime(st), _endTime(et), _tp(std::move(tp)), _mbr(mbr){};

    static bool CompareInterval(Interval *i1, Interval *i2);    // for sorting

    bool existTemporalOverlap(TIMETYPE t1, TIMETYPE t2) const;

    bool existSpatialOverlap(GeoPoint *queryP, float distThreshold = 0) const;

    inline TIMETYPE getStartTime() const {
        return _startTime;
    }

    inline TIMETYPE getEndTime() const {
        return _endTime;
    }

    TrajPointPair getTrajPointPair();

    IDTYPE getPointId() const;

    float getPointRadius() const;

    void releaseSpace();
};


#endif //PCT_INTERVAL_H
