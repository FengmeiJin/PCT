//
// Created by Fengmei JIN on 2022/8/29.
//

#ifndef PCT_RTREEELEMENT_H
#define PCT_RTREEELEMENT_H

#include "../../spatial/Rectangle.h"
#include "../../spatial/Point.h"
#include "../intervaltree/Interval.h"

class RtreeElement: public Rectangle{

    // ===== anonymize a point to another one
    // !!!the same center point may have various radius
    GeoPoint *centerPoint {nullptr};    // the rectangle depicts the circle region within a radius
    double radius{0};
    Interval *interval {nullptr};

public:
    RtreeElement(GeoPoint *_point, double _radius, Interval *_iv): Rectangle(_point, _radius) {
        // the MBR is defined based on the largest radius
        centerPoint = _point;
        radius = _radius;
        interval = _iv;
    }

    double computeDistToPoint(const GeoPoint *p) const {
        if(centerPoint == nullptr) {
            printf("Error in RtreeElement::computeDistToPoint: a null pointer!\n");
            return -1;
        }

        return centerPoint->computeDistance(p);
    }

    Interval* getInterval() {
        return interval;
    }

    inline double getRadius() const {
        return radius;
    }
};


#endif //PCT_RTREEELEMENT_H
