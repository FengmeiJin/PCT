//
// Created by Fengmei JIN on 2022/8/30.
//

#ifndef PCT_ITREENODE_H
#define PCT_ITREENODE_H

#include "Interval.h"
#include "../rtree/RtreeIndex.h"

// Each node contains a distinct interval and some information about its subtrees
class ItreeNode {
    vector<Interval *> intervals;
    RtreeIndex *rtree{nullptr};

    // Range: [min, max)
    TIMETYPE _minStartTime;     // minimum "low" value in subtree rooted with this node
    TIMETYPE _maxEndTime;       // maximum "high" value in subtree rooted with this node
//    Rectangle *_entireMBR;            // cover all subtrees under this node

    // subtrees rooted with this node
    ItreeNode *_leftChild{nullptr};
    ItreeNode *_rightChild{nullptr};

    // pointers to prev/next nodes on the same level
//    ItreeNode *_prevPtr{nullptr};
    ItreeNode *_nextPtr{nullptr};

    unsigned int _depth{0};            /* Root is 0, its subtrees are 1, ... */

protected:

    unordered_set<ItreeNode *> searchTimeOverlapNodes(TIMETYPE timestamp);
    TIMETYPE searchOverlapIntervals(vector<Interval *> &res, Point *queryP, GeoPoint *gp = nullptr,
                                    float distThreshold = 0, bool enableRtreeTimeFilter = false);

    inline bool fullCover(TIMETYPE t1, TIMETYPE t2 = 0) const;
    inline bool intersect(TIMETYPE t1, TIMETYPE t2) const;

public:

    ItreeNode(TIMETYPE _st, TIMETYPE _et) : _minStartTime(_st), _maxEndTime(_et) {};

    int insertInterval(Interval *iv);

    /**
     * if @param pIndex is valid, then search overlapping on "temporal->spatial"
     *      else, only search temporal overlapping
     * @return the time cost on Rtree if exists
     * */
    TIMETYPE overlapSearch(Point *queryP, vector<Interval *> &results, IndexPoint *pIndex = nullptr,
                           float contactMaxDist = 0, bool enableRtreeTimeFilter = false);

    inline TIMETYPE getMinStartTime() const {
        return this->_minStartTime;
    }

    inline TIMETYPE getMaxEndTime() const {
        return this->_maxEndTime;
    }

    /**
     * if @param pIndex is valid, then intervals are organized on spatial dimension further
     * */
    void organizeIntervals(IndexPoint *pIndex = nullptr);

    void setChildNode(ItreeNode *child, bool left = true);

    void setPointer(ItreeNode *node);

    unsigned int getHeight();

    void releaseSpace();
};


#endif //PCT_ITREENODE_H
