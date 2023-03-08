//
// Created by Fengmei Jin on 5/7/2022.
//

#ifndef PCT_GEOPOINT_H
#define PCT_GEOPOINT_H

#include "../commonheader.h"

#define DIMENSION 2
#define EARTH_RADIUS 6378.137f   //unit: km

using namespace std;

class GeoPoint {

    float* coordinates {nullptr};     // [0] longitude (x), [1] latitude (y)

public:
    GeoPoint() = default;
    GeoPoint(float lng, float lat): coordinates(new float[DIMENSION]{lng, lat}) {};
    explicit GeoPoint(GeoPoint *p): coordinates(new float[DIMENSION]{p->getLongitude(), p->getLatitude()}) {};

    [[nodiscard]] float getLongitude() const;
    [[nodiscard]] float getLatitude() const;

    static double rad(double degree) ;
    [[nodiscard]] double computeDistance(float lng, float lat) const;
    [[nodiscard]] double computeDistance(const GeoPoint *point) const;
    static double computeDistanceByCoors(float lng1, float lat1, float lng2, float lat2);   // in general use

    [[nodiscard]] string toString() const;

    void clearCoors();

    float getCoordinate(int axis);

    void updateCoordinates(float _lng, float _lat);
};


#endif //PCT_GEOPOINT_H
