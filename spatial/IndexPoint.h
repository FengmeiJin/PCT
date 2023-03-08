//
// Created by Fengmei Jin on 12/7/2022.
//

#ifndef PCT_INDEXPOINT_H
#define PCT_INDEXPOINT_H

#include "Rectangle.h"

class IndexPoint {

    unordered_map<float, unordered_map<float, IDTYPE>> lat2lng2pid;
    vector<GeoPoint *> pid2geopoint;

public:

    IndexPoint() {
        pid2geopoint.emplace_back(new GeoPoint);        // to occupy the invalid "0"-th
    }

    IDTYPE getVertexId(float lng, float lat);

    // -- default threshold in the paper
    // if an individual spent more than 30 minutes within a distance of 200 meters,
    // the region is detected as a stay point
    // unit here: km, second
    Rectangle *mineStayPoints(const vector<tuple<TIMETYPE, float, float>> &pointSeq,
                              vector<tuple<IDTYPE, TIMETYPE, TIMETYPE>> &stayPoints,
                              float distThreshold = 0.2, TIMETYPE timeThreshold = 1800);

    double computeDistance(IDTYPE id1, IDTYPE id2) const;

    GeoPoint *getPointById(IDTYPE id);

    inline IDTYPE getPointNum() const {
        return pid2geopoint.size() - 1;  // 0-th invalid
    }
};


#endif //PCT_INDEXPOINT_H
