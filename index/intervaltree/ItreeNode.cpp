//
// Created by Fengmei JIN on 2022/8/30.
//

#include "ItreeNode.h"
#include "../../utils.h"

bool ItreeNode::fullCover(TIMETYPE t1, TIMETYPE t2) const {
    return this->_minStartTime <= t1 && this->_maxEndTime > (t2 == 0 ? t1 : t2);    // Range: [min, max)
}

bool ItreeNode::intersect(TIMETYPE t1, TIMETYPE t2) const {
    return !(this->_minStartTime >= t2 || this->_maxEndTime <= t1);     // non-disjoint
}

int ItreeNode::insertInterval(Interval *iv) {
    TIMETYPE iv_st = iv->getStartTime(), iv_et = iv->getEndTime();
    if (!this->fullCover(iv_st, iv_et)) {
        printf("ERROR in ItreeNode::insertInterval "
               "- the given interval <%lu, %lu> cannot be covered by the root!\n", iv_st, iv_et);
        return -1;
    }
    stack<ItreeNode *> candi;
    candi.push(this);
    while (!candi.empty()) {
        ItreeNode *curr = candi.top();
        candi.pop();

        bool leftCover = false, rightCover = false;
        // the curr node fully covers the given interval but its children cannot
        if (curr->_leftChild != nullptr && curr->_leftChild->fullCover(iv_st, iv_et)) {
            candi.push(curr->_leftChild);
            leftCover = true;
        }
        if (curr->_rightChild != nullptr && curr->_rightChild->fullCover(iv_st, iv_et)) {
            candi.push(curr->_rightChild);
            rightCover = true;
        }
        if (!leftCover && !rightCover) {
            curr->intervals.emplace_back(iv);
            return 0;
        }
    }
    return -1;
}

unordered_set<ItreeNode *> ItreeNode::searchTimeOverlapNodes(TIMETYPE timestamp) {
    unordered_set<ItreeNode *> res;
    stack<ItreeNode *> candi;
    candi.push(this);
    while (!candi.empty()) {
        ItreeNode *curr = candi.top();
        candi.pop();

        // compare the duration of this subtree
        if (curr->fullCover(timestamp)) {
            // this interval overlaps with the given timestamp
            res.insert(curr);
            if (curr->_leftChild != nullptr) {
                candi.push(curr->_leftChild);
            }
            if (curr->_rightChild != nullptr) {
                candi.push(curr->_rightChild);
            }
        }
    }
    return res;
}

TIMETYPE ItreeNode::searchOverlapIntervals(vector<Interval *> &res, Point *queryP, GeoPoint *gp,
                                           float distThreshold, bool enableRtreeTimeFilter) {
    TIMETYPE t1 = queryP->getArriveTimestamp(), t2 = queryP->getLeaveTimestamp();
    if (gp != nullptr && this->rtree != nullptr) {
        vector<Interval *> ivs;
        TIMETYPE timer = millisecond(), rtreeSearchTime = 0;
        rtree->rangeQuery(gp, distThreshold, ivs, enableRtreeTimeFilter, t1, t2);
        rtreeSearchTime += (millisecond() - timer);

        if (!enableRtreeTimeFilter) {
            for (Interval *iv: ivs) {
                if (iv->existTemporalOverlap(t1, t2))
                    res.emplace_back(iv);
            }
        }
        else {
            res.insert(res.end(), ivs.begin(), ivs.end());
        }
        return rtreeSearchTime;
    }
    else {
        for (Interval *iv: this->intervals) {
            // intervals with temporal overlapping
            if (iv->existTemporalOverlap(t1, t2))
                res.emplace_back(iv);
        }
        return 0;
    }
}

// given an interval queryP
TIMETYPE ItreeNode::overlapSearch(Point *queryP, vector<Interval *> &res, IndexPoint *pIndex,
                                  float contactMaxDist, bool enableRtreeTimeFilter) {
    // the query time range should be covered by the interval tree
    TIMETYPE startT = queryP->getArriveTimestamp();   startT = MAX(this->_minStartTime, startT);
    TIMETYPE endT = queryP->getLeaveTimestamp();     endT = MIN(this->_maxEndTime, endT);
    unordered_set<ItreeNode *> startNodes = this->searchTimeOverlapNodes(startT);
    unordered_set<ItreeNode *> endNodes = this->searchTimeOverlapNodes(endT);

    // if spatial overlapping is not required, then these two are useless
    GeoPoint *gp = pIndex == nullptr ? nullptr : pIndex->getPointById(queryP->getPointId());
    float distThreshold = queryP->getRadius() + contactMaxDist;
    TIMETYPE rtreeSearchTime = 0;

    if (!startNodes.empty() && !endNodes.empty()) {
        unordered_set<ItreeNode *> checked;
        for (ItreeNode *startN: startNodes) {
            vector<ItreeNode *> path;
            path.emplace_back(startN);

            // check if there exists a valid path from any start node to end node
            ItreeNode *node = startN;
            while (node != nullptr && node->_nextPtr != nullptr && node->_nextPtr->intersect(startT, endT)) {
                node = node->_nextPtr;
                path.emplace_back(node);
            }
            if (node != nullptr && endNodes.find(node) != endNodes.end()) {
                // all nodes on the path overlap with the given point, but not for all intervals
                for (ItreeNode *n: path) {
                    if (!n->intervals.empty() && checked.find(n) == checked.end()) {
                        checked.insert(n);
                        rtreeSearchTime += n->searchOverlapIntervals(res, queryP, gp, distThreshold, enableRtreeTimeFilter);
                    }
                }
            }
        }
    }
    return rtreeSearchTime;
}

unsigned int ItreeNode::getHeight() {
    unsigned int rtn = 0;
    stack<pair<ItreeNode *, unsigned int>> candi;
    candi.push(make_pair(this, this->_depth + 1));   // "depth" starts from 0 (for the root)
    while (!candi.empty()) {
        ItreeNode *curr = candi.top().first;
        unsigned int height = candi.top().second;    // current height
        candi.pop();

        // left child
        if (curr->_leftChild != nullptr)
            candi.push(make_pair(curr->_leftChild, height + 1));
        else
            rtn = MAX(rtn, height);

        // right child
        if (curr->_rightChild != nullptr)
            candi.push(make_pair(curr->_rightChild, height + 1));
        else
            rtn = MAX(rtn, height);
    }
    return rtn;
}

void ItreeNode::setChildNode(ItreeNode *child, bool left) {
    if (left)
        this->_leftChild = child;
    else
        this->_rightChild = child;
}

void ItreeNode::setPointer(ItreeNode *node) {
    this->_nextPtr = node;
}

// prepare RtreeElement for rtree construction
int prepareRtreeElements(const vector<Interval *> &intervals, IndexPoint *pIndex, RtreeElement **elementList){

    // for each "same" GeoPoint, sharing a single RtreeElement
    unordered_map<IDTYPE, vector<Interval*>> gp2intervals;
    {
        for (Interval *iv: intervals) {
            IDTYPE pid = iv->getPointId();
            auto itr = gp2intervals.find(pid);
            if (itr != gp2intervals.end())
                itr->second.emplace_back(iv);
            else {
                vector<Interval *> vec;
                vec.emplace_back(iv);
                gp2intervals[pid] = vec;
            }
        }
    }

    int i = 0;
    for (const auto& gpPair: gp2intervals) {
        GeoPoint *point = pIndex->getPointById(gpPair.first);
        for (Interval *iv: gpPair.second) {
            elementList[i++] = new RtreeElement(point, iv->getPointRadius(), iv);
        }
    }
    return i;
}

void ItreeNode::organizeIntervals(IndexPoint *pIndex) {
    unsigned int avgIntervalNum = 0, validNode = 0, maxIntervalNum = 0;
    unsigned int avgRtreeHeight = 0, validRtree = 0, maxRtreeHeight = 0, minRtreeHeight = INFINITY;

    stack<ItreeNode *> res;
    res.push(this);
    while (!res.empty()) {
        ItreeNode *curr = res.top();
        int intervalNum = (int) curr->intervals.size();
        if (intervalNum > 0) {
            avgIntervalNum += intervalNum;
            validNode++;
            maxIntervalNum = MAX(maxIntervalNum, intervalNum);
            if(pIndex == nullptr) {
                // sorting intervals as organization
                sort(curr->intervals.begin(), curr->intervals.end(), Interval::CompareInterval);
            }
            else {
                // create spatial index for these intervals
                auto **elementList = new RtreeElement *[intervalNum];
                int elemNum = prepareRtreeElements(curr->intervals, pIndex, elementList);
                if (elemNum > 0) {
                    curr->rtree = new RtreeIndex(Dim, BranchFactor, LeafCapacity);
                    curr->rtree->constructRtree_STR(elementList, elemNum,Dim, BranchFactor, LeafCapacity);
                    int h = curr->rtree->getHeight();
                    avgRtreeHeight += h;
                    validRtree++;
                    maxRtreeHeight = MAX(maxRtreeHeight, h);
                    minRtreeHeight = MIN(minRtreeHeight, h);
                }
                elementList = nullptr;
            }
        }
        res.pop();

        if (curr->_leftChild != nullptr) {
            res.push(curr->_leftChild);
        }
        if (curr->_rightChild != nullptr) {
            res.push(curr->_rightChild);
        }
    }

#ifdef DEBUG
    printf("\n\t\t [Temporal Info] # intervals = %u, valid Nodes = %u, AVG intervals per node = %.3f, max interval num = %u\n",
           avgIntervalNum, validNode, avgIntervalNum * 1.0 / validNode, maxIntervalNum);
    if(pIndex != nullptr) {
        printf("\t\t [Spatial Info] valid Rtree = %u, AVG Rtree height = %.3f, max/min Rtree height = %u, %u\n",
               validRtree, avgRtreeHeight * 1.0 / validRtree, maxRtreeHeight, minRtreeHeight);
    }
#endif
}

void ItreeNode::releaseSpace() {
    stack<ItreeNode *> res;
    res.push(this);
    while (!res.empty()) {
        ItreeNode *curr = res.top();
        res.pop();
        if (curr->_leftChild != nullptr) {
            res.push(curr->_leftChild);
        }
        if (curr->_rightChild != nullptr) {
            res.push(curr->_rightChild);
        }
        // Clear this node
        for (Interval *iv: curr->intervals) {
            iv->releaseSpace();
            iv = nullptr;
        }
        curr->intervals.clear();
        curr->intervals.shrink_to_fit();
        if (curr->rtree != nullptr)
            curr->rtree->releaseSpace();
        curr->_leftChild = nullptr;
        curr->_rightChild = nullptr;
        curr->_nextPtr = nullptr;
        delete curr;
        curr = nullptr;
    }
}
