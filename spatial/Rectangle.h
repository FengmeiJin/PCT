//
// Created by Fengmei Jin on 11/7/2022.
//

#ifndef PCT_RECTANGLE_H
#define PCT_RECTANGLE_H

#include "GeoPoint.h"

class Rectangle {

protected:
    float *valueList {nullptr};     // describe the range of this rectangle

public:

    // ===== Constructors
    Rectangle() = default;

    explicit Rectangle(const GeoPoint *p, double _radius = 0, int _dim = DIMENSION);

    Rectangle(float minLng, float minLat, float maxLng, float maxLat, int _dim = DIMENSION) : valueList(
            new float[2 * _dim]{minLng, minLat, maxLng, maxLat}) {};

    explicit Rectangle(Rectangle *_rec, int _dim = DIMENSION);

    // ===== Relationship between rectangles
    bool fullContain(const Rectangle& _rec, int _dim = DIMENSION) const;
    bool intersect(const Rectangle& _rec, int _dim = DIMENSION) const;
    bool inside(const Rectangle& _rec, int _dim = DIMENSION) const;
    bool disjoint(const Rectangle& _rec, int _dim = DIMENSION) const;

    // ===== Relationship with point
    bool contain(float lng, float lat) const;

    // ===== Functions
    float getCoordinate(int axis, int _dim = DIMENSION, bool bottomLeft = true);

    float getValue(int axis) const;

    float *getMin() const;

    float *getMax(int _dim = DIMENSION) const;

    void enlarge(const Rectangle *_rec, int _dim = DIMENSION);

    void releaseSpace();

    char* toString() const;

    // ===== Computation
    static double computePerimeter(float minLng, float minLat, float maxLng, float maxLat);

    double computePerimeter() const;

    double perimeterIncrement(const Rectangle& _rec) const;

    int stateWithSphere(const GeoPoint *queryP, double _radius, int _dim) const;

    double computeMinDist(const Rectangle *rec) const;

    double computeMinDist(const GeoPoint *queryP) const;

    double computeMinDist(float lng, float lat) const;

    // offset in km (unit)
    static pair<double, double> getOffsetDegree(double offsetLng, double offsetLat, float pLat);
};


#endif //PCT_RECTANGLE_H
