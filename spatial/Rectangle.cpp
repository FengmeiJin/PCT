//
// Created by Fengmei Jin on 11/7/2022.
//

#include "Rectangle.h"

/**
*  Detect whether the given rectangle @param _rec is fully contained in this rectangle.
*  @return: If YES, return TRUE. Otherwise, return FALSE.
*/
bool Rectangle::fullContain(const Rectangle &_rec, int _dim) const {
    for (int i = 0; i < _dim; i++) {
        if (_rec.valueList[i] < this->valueList[i]) {
            return false;
        }
        if (_rec.valueList[i + _dim] > this->valueList[i + _dim]) {
            return false;
        }
    }
    return true;
}

bool Rectangle::contain(float lng, float lat) const {
    return lng >= valueList[0] && lng <= valueList[2] && lat >= valueList[1] && lat <= valueList[3];
}

/**
*  Detect whether "this" rectangle is fully inside the given rectangle @param _rec.
 *      (opposite with fullContain)
*  @return: If YES, return TRUE. Otherwise, return FALSE.
*/
bool Rectangle::inside(const Rectangle &_rec, int _dim) const {
    for (int i = 0; i < _dim; i++) {
        if (this->valueList[i] < _rec.valueList[i]) {
            return false;
        }
        if (this->valueList[i + _dim] > _rec.valueList[i + _dim]) {
            return false;
        }
    }
    return true;
}

bool Rectangle::disjoint(const Rectangle &_rec, int _dim) const {
    for (int i = 0; i < _dim; i++) {
        if (this->valueList[i] < _rec.valueList[i]) {
            return true;
        }
        if (this->valueList[i + _dim] > _rec.valueList[i + _dim]) {
            return true;
        }
    }
    return false;
}

/**
*  Detect whether this rectangle intersects with the given rectangle @param _rec.
*      (opposite with disjoint)
*  @return: If YES, return TRUE. Otherwise, return FALSE.
*/
bool Rectangle::intersect(const Rectangle &_rec, int _dim) const {
    for (int i = 0; i < _dim; i++) {
        if (this->valueList[i] >= _rec.valueList[i + _dim]) {    // min > max
            return false;
        }
        if (this->valueList[i + _dim] <= _rec.valueList[i]) {    // max < min
            return false;
        }
    }
    return true;
}

double degree(double radian) {
    return radian * 180 / M_PI;
}

// offset in km (unit)
pair<double, double> Rectangle::getOffsetDegree(double offsetLng, double offsetLat, float pLat) {
    // coordinate offsets in radians
    double dLat = offsetLat / EARTH_RADIUS;
    double dLng = offsetLng / (EARTH_RADIUS * cos(GeoPoint::rad(pLat)));

    return make_pair(degree(dLng), degree(dLat));
}

/**
 * @param _radius: offset in km (unit)
 * */
Rectangle::Rectangle(const GeoPoint *p, double _radius, int _dim) {
    float pLng = p->getLongitude(), pLat = p->getLatitude();
    if(_radius <= 0) {
        valueList = new float[2 * _dim]{pLng, pLat, pLng, pLat};    // a point represented by a rectangle
    }
    else {
        pair<double, double> rtn = getOffsetDegree(_radius, _radius, pLat);
        double degreeLng = rtn.first, degreeLat = rtn.second;

        // Offset Position, decimal degrees
        auto minLng = (float) (pLng - degreeLng);
        auto maxLng = (float) (pLng + degreeLng);
        auto minLat = (float) (pLat - degreeLat);
        auto maxLat = (float) (pLat + degreeLat);
        valueList = new float [2 * _dim]{minLng, minLat, maxLng, maxLat};
    }
}

Rectangle::Rectangle(Rectangle *_rec, int _dim) {
    this->valueList = new float[2 * _dim];
    for (int i = 0; i < 2 * _dim; i++) {
        this->valueList[i] = _rec->valueList[i];
    }
}

/*
*  Release the space of this rectangle.
*/
void Rectangle::releaseSpace() {
    if (this->valueList != nullptr) {
        delete[](this->valueList);
        this->valueList = nullptr;
    }
}

float Rectangle::getCoordinate(int axis, int _dim, bool bottomLeft) {
    if (valueList != nullptr && axis < DIMENSION * 2) {
        return bottomLeft ? valueList[axis] : valueList[axis + _dim];
    }
    return INFINITY;
}

float Rectangle::getValue(int axis) const {
    if (valueList != nullptr && axis < DIMENSION * 2) {
        return valueList[axis];
    }
    return INFINITY;
}

/**
*  Find the minimum bounding rectangle of this and the given rectangles.
*  @param _rec:	the given rectangle
*  @return : 	the resulting MBR
*/
void Rectangle::enlarge(const Rectangle *_rec, int _dim) {
    for (int i = 0; i < _dim; i++) {
        this->valueList[i] = min(this->valueList[i], _rec->valueList[i]);    // [0] minLng, [1] minLat
        this->valueList[i + _dim] = max(this->valueList[i + _dim], _rec->valueList[i + _dim]);   // [2] maxLng, [3] maxLat
    }
}

double Rectangle::computePerimeter(float minLng, float minLat, float maxLng, float maxLat) {
    double left = GeoPoint::computeDistanceByCoors(minLng, minLat, minLng, maxLat);
    double top = GeoPoint::computeDistanceByCoors(minLng, maxLat, maxLng, maxLat);
    double right = GeoPoint::computeDistanceByCoors(maxLng, minLat, maxLng, maxLat);
    double bottom = GeoPoint::computeDistanceByCoors(minLng, minLat, maxLng, minLat);
    return left + top + right + bottom;
}

/**
*  Compute the "perimeter" of this rectangle.
*  @return :    the sum of the length in each dimension.
*/
double Rectangle::computePerimeter() const {
    return computePerimeter(valueList[0], valueList[1], valueList[2], valueList[3]);
}

/**
*  Compute the increment of the "perimeter" if the given rectangle is merged with this rectangle.
*  @param _rec:	the given rectangle
*  @return : 	the sum of increments of the length in each dimension.
*/
double Rectangle::perimeterIncrement(const Rectangle &_rec) const {
    if (this->fullContain(_rec)) {
        return 0;
    }
    double perimeter = this->computePerimeter();
    if (this->inside(_rec)) {
        return _rec.computePerimeter() - perimeter;
    }
    auto *tmp = new Rectangle(_rec);
    tmp->enlarge(this);
    return tmp->computePerimeter() - perimeter;
}

float *Rectangle::getMin() const {
    return this->valueList;
}

float *Rectangle::getMax(int _dim) const {
    return valueList + _dim;
}

/*
*  Compute the state of this rectangle with the specified sphere.
*  Return:
*  	1:		this rectangle is fully inside the sphere.
*  	0:		this rectangle intersects with the sphere.
*  	-1:		this rectangle is fully outside the sphere.
*/
int Rectangle::stateWithSphere(const GeoPoint *queryP, double _radius, int _dim) const {

    Rectangle approxCircle(queryP, _radius, _dim);  // a rectangle instead

    int state = -1;     // fully outside
    if(this->inside(approxCircle) || this->fullContain(approxCircle)) {
        state = 1;      // fully inside
    }
    else if(this->intersect(approxCircle)) {
        state = 0;      // intersect
    }

    return state;
}

char* Rectangle::toString() const {
    char *buffer = new char[100];
//    sprintf(buffer, "(%.7f, %.7f), (%.7f, %.7f)", valueList[0], valueList[1], valueList[2], valueList[3]);
    snprintf(buffer, 100, "%.7f,%.7f,%.7f,%.7f", valueList[0], valueList[1], valueList[2], valueList[3]);
    return buffer;
}

double Rectangle::computeMinDist(float lng, float lat) const {
    if(this->contain(lng, lat)) {
        return 0;
    }
    float x1_min = valueList[0], y1_min = valueList[1], x1_max = valueList[2], y1_max = valueList[3];
    float x2 = lng, y2 = lat;

    bool left = x2 < x1_min;        // queryP is on the left side of this
    bool right = x1_max < x2;       // queryP is on the right side of this
    bool bottom = y2 < y1_min;      // queryP is on the bottom side of this
    bool top = y1_max < y2;         // queryP is on the top side of this

    if (top && left) {
        return GeoPoint::computeDistanceByCoors(x2, y2, x1_min, y1_max);
    }
    if (top && right) {
        return GeoPoint::computeDistanceByCoors(x2, y2, x1_max, y1_max);
    }
    if (bottom && left) {
        return GeoPoint::computeDistanceByCoors(x2, y2, x1_min, y1_min);
    }
    if (bottom && right) {
        return GeoPoint::computeDistanceByCoors(x2, y2, x1_max, y1_min);
    }
    if(top) {
        return GeoPoint::computeDistanceByCoors(x2, y2, x2, y1_max);
    }
    if(bottom) {
        return GeoPoint::computeDistanceByCoors(x2, y2, x2, y1_min);
    }
    if(left) {
        return GeoPoint::computeDistanceByCoors(x2, y2, x1_min, y2);
    }
    if(right) {
        return GeoPoint::computeDistanceByCoors(x2, y2, x1_max, y2);
    }
    return 0;
}

double Rectangle::computeMinDist(const GeoPoint *queryP) const {
    return computeMinDist(queryP->getLongitude(), queryP->getLatitude());
}

// compute minimum distance between two rectangles
double Rectangle::computeMinDist(const Rectangle *rec) const {
    if (this->intersect(*rec) || this->fullContain(*rec) || this->inside(*rec)) {
        return 0;
    }

    float x1_min = valueList[0], y1_min = valueList[1], x1_max = valueList[2], y1_max = valueList[3];
    float x2_min = rec->valueList[0], y2_min = rec->valueList[1], x2_max = rec->valueList[2], y2_max = rec->valueList[3];

    bool left = x2_max < x1_min;        // rec is on the left side of this
    bool right = x1_max < x2_min;       // rec is on the right side of this
    bool bottom = y2_max < y1_min;      // rec is on the bottom side of this
    bool top = y1_max < y2_min;         // rec is on the top side of this

    if (top && left){
        return GeoPoint::computeDistanceByCoors(x1_min, y1_max, x2_max, y2_min);
    }
    if (top && right) {
        return GeoPoint::computeDistanceByCoors(x1_max, y1_max, x2_min, y2_min);
    }
    if (bottom && left) {
        return GeoPoint::computeDistanceByCoors(x1_min, y1_min, x2_max, y2_max);
    }
    if (bottom && right) {
        return GeoPoint::computeDistanceByCoors(x1_max, y1_min, x2_min, y2_max);
    }
    if (top) {
        float xa = MAX(x1_min, x2_min);
        float xb = MIN(x1_max, x2_max);
        double d1 = GeoPoint::computeDistanceByCoors(xa, y1_max, xa, y2_min);
        double d2 = GeoPoint::computeDistanceByCoors(xb, y1_max, xb, y2_min);
        return MIN(d1, d2);
    }
    if (bottom) {
        float xa = MAX(x1_min, x2_min);
        float xb = MIN(x1_max, x2_max);
        double d1 = GeoPoint::computeDistanceByCoors(xa, y1_min, xa, y2_max);
        double d2 = GeoPoint::computeDistanceByCoors(xb, y1_min, xb, y2_max);
        return MIN(d1, d2);
    }
    if (left) {
        float ya = MAX(y1_min, y2_min);
        float yb = MIN(y1_max, y2_max);
        double d1 = GeoPoint::computeDistanceByCoors(x1_min, ya, x2_max, ya);
        double d2 = GeoPoint::computeDistanceByCoors(x1_min, yb, x2_max, yb);
        return MIN(d1, d2);
    }
    if (right) {
        float ya = MAX(y1_min, y2_min);
        float yb = MIN(y1_max, y2_max);
        double d1 = GeoPoint::computeDistanceByCoors(x1_max, ya, x2_min, ya);
        double d2 = GeoPoint::computeDistanceByCoors(x1_max, yb, x2_min, yb);
        return MIN(d1, d2);
    }
    return 0;   // should be unreachable
}
