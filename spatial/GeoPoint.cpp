//
// Created by Fengmei Jin on 5/7/2022.
//

#include "GeoPoint.h"

float GeoPoint::getLongitude() const {
    return coordinates[0];
}

float GeoPoint::getLatitude() const {
    return coordinates[1];
}

double GeoPoint::rad(double degree) {
    return degree * M_PI / 180.0;
}

double GeoPoint::computeDistanceByCoors(float lng1, float lat1, float lng2, float lat2){
    double a = rad(lng1) - rad(lng2);
    double radLat1 = rad(lat1);
    double radLat2 = rad(lat2);
    double b = radLat1 - radLat2;
    double dis = EARTH_RADIUS * 2 * asin(sqrt(pow(sin(b / 2), 2) + cos(radLat1) * cos(radLat2) * pow(sin(a/2),2)));
    return round(dis * 100000.0) / 100000;
}

double GeoPoint::computeDistance(float lng, float lat) const {
    if(this->getLongitude() == lng && this->getLatitude() == lat){
        return 0;
    }
    return computeDistanceByCoors(this->getLongitude(), this->getLatitude(), lng, lat);
}

double GeoPoint::computeDistance(const GeoPoint *point) const {
    if(this->getLongitude() == point->getLongitude() && this->getLatitude() == point->getLatitude()){
        return 0;
    }
    return computeDistanceByCoors(this->getLongitude(), this->getLatitude(), point->getLongitude(), point->getLatitude());
}

string GeoPoint::toString() const {
    string str(to_string(coordinates[0]));
    str.append(",");
    str.append(to_string(coordinates[1]));
    return str;
}

void GeoPoint::clearCoors() {
    delete[] coordinates;
    coordinates = nullptr;
}

float GeoPoint::getCoordinate(int axis) {
    if(axis < DIMENSION) {
        return coordinates[axis];
    }
    return INFINITY;
}

void GeoPoint::updateCoordinates(float _lng, float _lat) {
    if(coordinates != nullptr) {
        coordinates[0] = _lng;
        coordinates[1] = _lat;
    }
    else {
        coordinates = new float[DIMENSION]{_lng, _lat};
    }
}
