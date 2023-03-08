//
// Created by Fengmei JIN on 2022/8/25.
//

#ifndef PCT_RTREENODE_H
#define PCT_RTREENODE_H


#include "RtreeElement.h"

class RtreeNode {

private:
    /*
    *  For an internal node, listPtr is the pointer of the child node list.
    *  For a leaf node, listPtr is the pointer of the element list.
    *      Note that element can be a single point, a circle centered by a point, multiple points
    *      by all means, these geometries can be defined by a rectangle
    */
    int listSize{0};

    char **listPtr{nullptr};

    TIMETYPE minStartTime{0};
    TIMETYPE maxEndTime{0};

    static double computeGroupPerimeter(bool _isLeaf, char** _list, int _start, int _length, int _dim);

public:
    Rectangle *mbr{nullptr};
    RtreeNode *parent{nullptr};

    RtreeNode() = default;

    virtual ~RtreeNode();

    /*
     * element can be a single point, a circle centered by a point, multiple points
     * by all means, these geometries can be defined by a rectangle
     * */
    int addElement(RtreeElement *_elem, int _dim);

    int addChildNode(RtreeNode *_child, int _dim);

    int deleteChildNode(RtreeNode *_child);

    void initialize(bool _isLeaf, int _branchFactor, int _leafCapacity);

    RtreeNode *chooseBestChildNode(Rectangle *_rect);

    RtreeNode *split(bool _isLeaf, char *_toInsPtr, int _dim, int _branchFactor, int _leafCapacity);

    void releaseSpace();

    inline int getListSize() const {
        return listSize;
    }

    inline char **getListPtr() {
        return listPtr;
    }

    bool temporalOverlap(TIMETYPE queryStartT, TIMETYPE queryEndT) const;
};


#endif //PCT_RTREENODE_H
