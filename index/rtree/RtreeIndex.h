//
// Created by Fengmei JIN on 2022/8/23.
//

#ifndef PCT_RTREEINDEX_H
#define PCT_RTREEINDEX_H

#include "../../spatial/Trip.h"
#include "RtreeNode.h"
#include "RtreeElement.h"

const int Dim = 2;
const int BranchFactor = 4;
const int LeafCapacity = 4;

// This value indicates the number of child nodes of an internal node
// is in the range of [fanout/FanoutRatio, fanout].
const int FanoutRatio = 4;

/*
*  The tool for sorting points by the specified dimension coordinates.
*/
struct RectangleDimSorter {
    int dim;
    int dimIndex; // in [0, dim-1]
    RectangleDimSorter(int _dim, int _dimLev) : dim(_dim), dimIndex(_dimLev - 1) {};

    bool operator()(const RtreeElement *r1, const RtreeElement *r2) const {
//        return p1->coords[dimIndex] < p2->coords[dimIndex];
//        return p1->getValue(dimIndex) < p2->getValue(dimIndex);

        float coord1 = r1->getValue(dimIndex + dim) + r1->getValue(dimIndex);
        float coord2 = r2->getValue(dimIndex + dim) + r2->getValue(dimIndex);
        return coord1 < coord2;
    }
};

struct NodeDimSorter {
    int dim;
    int dimIndex; // in [0, dim -1]
    NodeDimSorter(int _dim, int _dimLev) : dim(_dim), dimIndex(_dimLev - 1) {};

    bool operator()(const RtreeNode *node1, const RtreeNode *node2) const {
        float coord1 = node1->mbr->getValue(dimIndex + dim) + node1->mbr->getValue(dimIndex);
        float coord2 = node2->mbr->getValue(dimIndex + dim) + node2->mbr->getValue(dimIndex);
        return coord1 < coord2;
    }
};

class RtreeIndex{

private:
    RtreeNode *root;    /* Leaf nodes are at level-0 and elements are at level-(-1). */
    int rootLev;

    // Rtree capacity parameters
    int dim;
    int branchFactor;
    int leafCapacity;

protected:

    void insertNode_SubRoutine(RtreeNode *_subroot, int _subrootLevel, RtreeNode *_toInsNode, int _targetLevel);

    int insertElement_SubRoutine(RtreeNode *_subroot, int _subrootLevel, RtreeElement *_elem);

    void overflowHandler(RtreeNode *_overflowNode, RtreeNode *_toInsNode);

    void underflowHandler(RtreeNode *_underflowNode, int _nodeLevel);

    void releaseSpace_SubRoutine(RtreeNode *_subroot, int _subrootLev);

    int insertNode(RtreeNode *_toInsNode, int _targetLevel);

    int insertElement(RtreeElement *_elem);

    /***************************************************
    *  STR bulk loading related.
    **************************************************/
    /**
    *  Construct leaf level by STR algorithm.
    *  @param _leafLevel:	the target place to store the pointers of the leaf nodes.
    *  @param _dimLev:		which dimension we used for sorting in the current recursion.
    */
    int constructLeafLevel_STR_SubRoutine(RtreeElement **_elementList, int _elemNum, int _dim, int _branchFactor,
                                          int _leafCapacity, vector<RtreeNode *> &_leafLevel, int _dimLev);

    int constructNextLevel_STR_SubRoutine(RtreeNode **_nodeList, int _nodeNum, int _dim, int _branchFactor,
                                          int _leafCapacity, vector<RtreeNode *> &_nextLevel, int _dimLev);

public:
    RtreeIndex();

    RtreeIndex(int _dim, int _branchFactor, int _leafCapacity);

    virtual ~RtreeIndex();

    void releaseSpace();

    /*
    * 	The STR bulk loading algorithm.
    */
    RtreeNode *constructRtree_STR(RtreeElement **_elementList, int _elemNum, int _dim, int _branchFactor, int _leafCapacity);

    /*
    *  Added by jhgan on 2017-01-23.
    *  The implementation of range reporting.
    */
    int rangeQuery(GeoPoint *queryP, double _radius, vector<Interval *> &_targetPlace,
    bool enableRtreeTimeFilter = false, TIMETYPE queryStartT = 0, TIMETYPE queryEndT = 0);

    int getHeight() const;
};


#endif //PCT_RTREEINDEX_H
