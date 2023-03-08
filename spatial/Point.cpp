//
// Created by Fengmei Jin on 6/7/2022.
//

#include "Point.h"

IDTYPE Point::getPointId() const {
    return pointID;
}

Point *Point::getPoint_ori() const {
    return ptr_ori;
}

TIMETYPE Point::getArriveTimestamp() const {
    return arvTime;
}

TIMETYPE Point::getLeaveTimestamp() const {
    return levTime;
}

void Point::setLeaveTimestamp(TIMETYPE newLeaveTime) {
    levTime = newLeaveTime;
}

float Point::getRadius() const {
    return radius;
}

bool Point::isHeadPtr() const{
    return prev_ptr == nullptr;     // the head pointer of a trip
}

bool Point::isTailPtr() const{
    return next_ptr == nullptr;     // the tail pointer of a trip
}

Point* Point::getPrevPointer() const {
    return prev_ptr;
}

Point* Point::getNextPointer() const {
    return next_ptr;
}

void Point::setPrevPointer(Point *ptr) {
    prev_ptr = ptr;
}

void Point::setNextPointer(Point *ptr) {
    next_ptr = ptr;
}

void Point::clear() {
    prev_ptr = nullptr;
    next_ptr = nullptr;
}

// insert this point before the *currP
void Point::insertBefore(Point *currP) {
    this->prev_ptr = currP->prev_ptr;
    this->next_ptr = currP;
    currP->prev_ptr->next_ptr = this;
    currP->prev_ptr = this;
}

// insert this point after the *currP
void Point::insertAfter(Point *currP) {
    this->prev_ptr = currP;
    this->next_ptr = currP->next_ptr;
    currP->next_ptr->prev_ptr = this;
    currP->next_ptr = this;
}

void Point::unlink() const {
    this->prev_ptr->next_ptr = this->next_ptr;
    this->next_ptr->prev_ptr = this->prev_ptr;
}

TIMETYPE Point::getDuration() const {
    return levTime - arvTime;
}

bool Point::existTemporalOverlap(Point *p1, Point *p2) {
    if(p1->arvTime >= p2->levTime || p1->levTime <= p2->arvTime) {
        return false;
    }
    return true;
}

bool Point::equalPoiDistribution(const float *distribution, unsigned int array_len) const {
    if(this->poiDistribution == nullptr || distribution == nullptr) {
        return false;
    }
    for(int i = 0; i < array_len; i++) {
        if(distribution[i] != this->poiDistribution[i])
            return false;
    }
    return true;
}

float *Point::getPoiDistribution() const {
    return poiDistribution;
}
