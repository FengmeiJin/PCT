//
// Created by Fengmei JIN on 2022/9/26.
//

#ifndef PCT_ITREEINDEX_H
#define PCT_ITREEINDEX_H


#include "ItreeNode.h"
#include "../../spatial/Trip.h"

class ItreeIndex {

private:
    ItreeNode *root {nullptr};

public:
    ItreeIndex() = default;
    void releaseSpace() {
        if (root != nullptr) {
            root->releaseSpace();
        }
    }

    void constructIndexByIntervals(const vector<Interval *> &sortedIntervals, IndexPoint *pIndex = nullptr);

    void constructEmptyTree(TIMETYPE earliestIndexTime, TIMETYPE LatestEndTime, TIMETYPE contactMinDuration, int deltaT);

    TIMETYPE overlapSearch(Point *queryP, vector<Interval *> &results,
                           IndexPoint *pIndex = nullptr, float contactMaxDist = 0, bool enableRtreeTimeFilter = false) {
        if (root != nullptr) {
            return root->overlapSearch(queryP, results, pIndex, contactMaxDist, enableRtreeTimeFilter);
        }
        return 0;
    }

    unsigned int getHeight() const {
        return this->root == nullptr ? 0 : this->root->getHeight();
    }
};


#endif //PCT_ITREEINDEX_H
