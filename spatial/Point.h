//
// Created by Fengmei Jin on 6/7/2022.
//

#ifndef PCT_POINT_H
#define PCT_POINT_H


#include "../commonheader.h"

using namespace std;

class Point {
    IDTYPE pointID;         // uniquely represent a GeoPoint, if = 0, then invalid
    Point *ptr_ori{nullptr};

    // temporal information will not be changed by perturbation, unit: second
    TIMETYPE arvTime;       // arrive timestamp
    TIMETYPE levTime;       // leave timestamp
    float radius{0};           // the range centered by the current point, if = 0, then exact point

    Point* prev_ptr{nullptr};
    Point* next_ptr{nullptr};

    float* poiDistribution{nullptr};

public:

    Point(): pointID(0), arvTime(0), levTime(0){};  // head or tail

    Point(IDTYPE id, TIMETYPE t): pointID(id), arvTime(t), levTime(t){};    // original point without further info

    Point(IDTYPE id, Point *p_ori, TIMETYPE t1, TIMETYPE t2,
          float _r = 0, float *distribution = nullptr) : pointID(id), ptr_ori(p_ori), arvTime(t1), levTime(t2),
                                                         radius(_r), poiDistribution(distribution) {};

    void clear();

    [[nodiscard]] IDTYPE getPointId() const;
    [[nodiscard]] Point *getPoint_ori() const;

    [[nodiscard]] TIMETYPE getArriveTimestamp() const;
    [[nodiscard]] TIMETYPE getLeaveTimestamp() const;
    void setLeaveTimestamp(TIMETYPE newLeaveTime);
    TIMETYPE getDuration() const;

    [[nodiscard]] bool isHeadPtr() const;
    [[nodiscard]] bool isTailPtr() const;

    [[nodiscard]] Point *getPrevPointer() const;
    [[nodiscard]] Point *getNextPointer() const;
    void setPrevPointer(Point *ptr);
    void setNextPointer(Point *ptr);

    void insertBefore(Point *currP);
    void insertAfter(Point *currP);
    void unlink() const;

    [[nodiscard]] float getRadius() const;
    [[nodiscard]] float *getPoiDistribution() const;

    static bool existTemporalOverlap(Point *p1, Point *p2);

    static bool CompareByTime (Point *p1, Point *p2) {
        return p1->arvTime < p2->arvTime;
    }

    bool equalPoiDistribution(const float *distribution, unsigned int array_len) const;

};


#endif //PCT_POINT_H
