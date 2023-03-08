//
// Created by Fengmei Jin on 12/7/2022.
//

#include "IndexPoint.h"

IDTYPE IndexPoint::getVertexId(float lng, float lat) {
    auto itr = lat2lng2pid.find(lat);
    if(itr != lat2lng2pid.end()) {
        auto itr2 = itr->second.find(lng);
        if(itr2 != itr->second.end()) {
            return itr2->second;    // already exist
        }
    }

    // not exist
    IDTYPE pid = pid2geopoint.size();   // valid ID starting from 1
    lat2lng2pid[lat].insert(make_pair(lng, pid));
    pid2geopoint.emplace_back(new GeoPoint(lng, lat));
    return pid;
}

pair<float, float> computeMean(const vector<tuple<TIMETYPE, float, float>>& pointSeq,
                               unsigned int arv, unsigned int lev){
    float lng = 0, lat = 0;
    for(unsigned int i = arv, len = pointSeq.size(); i <= lev && i < len; i++) {
        auto p= pointSeq.at(i);
        lng += get<1>(p);
        lat += get<2>(p);
    }
    auto num = (float) (lev - arv + 1);
    return make_pair(lng / num, lat / num);
}

// ACM GIS'08: Mining user similarity based on location history
Rectangle* IndexPoint::mineStayPoints(const vector<tuple<TIMETYPE, float, float>>& pointSeq,
                                      vector<tuple<IDTYPE, TIMETYPE, TIMETYPE>> &stayPoints,
                                      float distThreshold, TIMETYPE timeThreshold) {
    unsigned int pointNum = pointSeq.size(), i = 1;

    float minLng = 200, maxLng = -200;  // longitude: [-180, 180]
    float minLat = 100, maxLat = -100;  // latitude: [-90, 90]

    {   // the first point must be a stay point
        auto point = pointSeq.at(0);
        float lng = get<1>(point), lat = get<2>(point);
        IDTYPE id = getVertexId(lng, lat);
        TIMETYPE arvTime = get<0>(point), levTime = arvTime;
        if(pointNum > 1) {
            levTime = get<0>(pointSeq.at(1));
        }
        stayPoints.emplace_back(make_tuple(id, arvTime, levTime));

        minLng = MIN(minLng, lng);
        maxLng = MAX(maxLng, lng);
        minLat = MIN(minLat, lat);
        maxLat = MAX(maxLat, lat);
    }

    // process the intermediate points
    while (i < pointNum - 1) {
        auto startP = pointSeq.at(i);
        float lng = get<1>(startP), lat = get<2>(startP);
        TIMETYPE arvTime = get<0>(startP);
        unsigned int j = i + 1;
        while (j < pointNum) {
            auto currP = pointSeq.at(j);
            double dist = GeoPoint::computeDistanceByCoors(lng, lat, get<1>(currP), get<2>(currP));
            if(dist > distThreshold || j == pointNum - 1){  // the last point
                TIMETYPE levTime = get<0>(currP);
                TIMETYPE duration = levTime - arvTime;
                if(duration > timeThreshold) {   // stay for a long period
                    pair<float, float> centerPoint = computeMean(pointSeq, i, j);
                    float lngC = centerPoint.first, latC = centerPoint.second;
                    IDTYPE pid = getVertexId(lngC, latC);
#ifdef DEBUG
                    auto tp = stayPoints.back();
                    if(arvTime < get<1>(tp) || arvTime < get<2>(tp)) {
                        printf("[ERROR]\n");
                    }
#endif
                    stayPoints.emplace_back(make_tuple(pid, arvTime, levTime));

                    minLng = MIN(minLng, lngC);
                    maxLng = MAX(maxLng, lngC);
                    minLat = MIN(minLat, latC);
                    maxLat = MAX(maxLat, latC);
                }
                i = j;
                break;
            }
            j++;
        }
    }

    if(pointNum > 1) {
        auto point = pointSeq.at(pointNum - 1);   // the end point must be a stay point
        float lng = get<1>(point), lat = get<2>(point);
        IDTYPE id = getVertexId(lng, lat);
        TIMETYPE arvTime = get<0>(point);
        stayPoints.emplace_back(make_tuple(id, arvTime, arvTime));

        minLng = MIN(minLng, lng);
        maxLng = MAX(maxLng, lng);
        minLat = MIN(minLat, lat);
        maxLat = MAX(maxLat, lat);
    }

    return new Rectangle(minLng, minLat, maxLng, maxLat);
}

double IndexPoint::computeDistance(IDTYPE id1, IDTYPE id2) const {
    if(id1 == id2) {
        return 0;
    }
    if(id1 >= pid2geopoint.size() || id2 >= pid2geopoint.size()) {  // invalid
        return -1;
    }
    GeoPoint *p1 = pid2geopoint[id1];
    GeoPoint *p2 = pid2geopoint[id2];
    double dist = -1;
    if(p1 != nullptr && p2 != nullptr) {
        dist = p1->computeDistance(p2);
    }

    return dist;
}

GeoPoint *IndexPoint::getPointById(IDTYPE id) {
    if(id == 0 || id >= pid2geopoint.size()) {
        return nullptr;
    }
    return pid2geopoint[id];
}
