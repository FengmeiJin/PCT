
//
// Created by Fengmei Jin on 6/7/2022.
//

#ifndef PCT_TRIP_H
#define PCT_TRIP_H


#include "Point.h"
#include "StaticGrid.h"
#include "Rectangle.h"

struct IndisRegion {
    Point *oriPtr;
    float lng, lat;     // these two will define a new PID
    float radius;       // the distance from the new point to the original one

    IndisRegion(Point *ptr, float _lng, float _lat, float _r) : oriPtr(ptr), lng(_lng), lat(_lat), radius(_r) {};
};

class Trip {
    string uid;

    Point* head_ptr{nullptr};    // just pointer, no data
    Point* tail_ptr{nullptr};
    unsigned int length{0};     // exclude the head/tail

    Rectangle *rec{nullptr};     // spatial range
    pair<TIMETYPE, TIMETYPE> timeRange;   // temporal range, [0] min, [1] max

public:

    // Constructor for simulated check-in sequences (aka POI sequences)
    Trip(const string &_uid, const vector<tuple<IDTYPE, TIMETYPE, TIMETYPE>> &poiSeq, IndexPoint *pIndex,
         StaticGrid *poiGrid, bool needPoiDistribution);

    // Constructor for noisy trip (each tuple with a radius)
    Trip(const string &_uid, const vector<IndisRegion> &noisyRegions, IndexPoint *pIndex,
         StaticGrid *poiGrid = nullptr, bool needPoiDistribution = false, int alignRaw2POI = 1);

    inline string getId() const {
        return uid;
    }

    inline unsigned int getLength() const {
        return length;
    }

    inline pair<TIMETYPE, TIMETYPE> getTimeRange() const {
        return timeRange;
    }

    [[nodiscard]] Point * getFirstPoint() const;

    void clearAllPoints();

    bool spatialOverlap(const Trip *pTrip, float distThreshold = 0) const;

    bool temporalOverlap(const Trip *pTrip) const;

    tuple<TIMETYPE, int, double> contactWith(const Trip *spreaderTraj, TIMETYPE validSpreadStartTime, IndexPoint *pIndex, float distThreshold,
                                             TIMETYPE contactMinDuration, bool enableCP = false, bool independentCP = true,
                                             TIMETYPE earliestCST = 0) const;
};

bool sortByStartTime(Trip *t1, Trip *t2);


#endif //PCT_TRIP_H
