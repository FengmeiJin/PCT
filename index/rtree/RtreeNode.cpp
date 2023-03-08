//
// Created by Fengmei JIN on 2022/8/25.
//

#include "RtreeNode.h"
#include <algorithm>

/*
*  The largest integer.
*/
//const int MYINFTY = 2147483647;
const double MYINFTY = INFINITY;

/*
*  axisIndex is in the range of [0, DIM - 1].
*/
int axisIndex = 0;

/*
*  The dimensionality of rtree nodes.
*/
int rtreeNode_dim = 2;

/*
*  Compare two given tree node pointers by the specified axisIndex.
*/
bool CompRtreeNodeLess(char* _u1, char* _u2) {
    auto* u1 = (RtreeNode*)_u1;
    auto* u2 = (RtreeNode*)_u2;
    float c1 = u1->mbr->getCoordinate(axisIndex, rtreeNode_dim) + u1->mbr->getCoordinate(axisIndex, rtreeNode_dim, false);
    float c2 = u2->mbr->getCoordinate(axisIndex, rtreeNode_dim) + u2->mbr->getCoordinate(axisIndex, rtreeNode_dim, false);
    return c1 < c2;
}

/*
*  Compare two given element pointers by the specified axisIndex.
*/
bool CompElementLess(char* _u1, char* _u2) {
    auto* m1 = (Rectangle*)_u1;
    auto* m2 = (Rectangle*)_u2;
    float c1 = m1->getCoordinate(axisIndex, rtreeNode_dim) + m1->getCoordinate(axisIndex, rtreeNode_dim, false);
    float c2 = m2->getCoordinate(axisIndex, rtreeNode_dim) + m2->getCoordinate(axisIndex, rtreeNode_dim, false);
    return c1 < c2;
}

RtreeNode::~RtreeNode() {
    this->releaseSpace();
}

// ===== for leaf node

int RtreeNode::addElement(RtreeElement* _elem, int _dim) {
    // Update the MBR for the parent node.
    if (this->listSize == 0) {
        this->mbr = new Rectangle(_elem, _dim);   // initialize the mbr of this node
        this->minStartTime = _elem->getInterval()->getStartTime();
        this->maxEndTime = _elem->getInterval()->getEndTime();
    }
    else if (this->mbr == nullptr) {
        printf("ERROR in RtreeNode::addElement: a null mbr\n");
        return -1;
    }
    else {
        this->mbr->enlarge(_elem, _dim);
        this->minStartTime = MIN(this->minStartTime, _elem->getInterval()->getStartTime());
        this->maxEndTime = MAX(this->maxEndTime, _elem->getInterval()->getEndTime());
    }

    if(this->listPtr == nullptr) {
        printf("ERROR in RtreeNode::addElement: a null listPtr\n");
        return -1;
    }
    this->listPtr[listSize++] = (char*)_elem;

    return listSize;
}

// ===== for non-leaf node

int RtreeNode::addChildNode(RtreeNode* _child, int _dim) {
    // Update the MBR for the parent node.
    if (this->listSize == 0) {
        this->mbr = new Rectangle(_child->mbr, _dim);     // same as its child's mbr
        this->minStartTime = _child->minStartTime;
        this->maxEndTime = _child->maxEndTime;
    }
    else if (this->mbr == nullptr) {
        printf("ERROR in RtreeNode::addChildNode: a null mbr\n");
        return -1;
    }
    else {
        this->mbr->enlarge(_child->mbr, _dim);
        this->minStartTime = MIN(this->minStartTime, _child->minStartTime);
        this->maxEndTime = MAX(this->maxEndTime, _child->maxEndTime);
    }

    if(this->listPtr == nullptr) {
        printf("ERROR in RtreeNode::addChildNode: a null listPtr\n");
        return -1;
    }
    this->listPtr[listSize++] = (char*)_child;
    _child->parent = this;      // set the child's parent node
    return listSize;
}

int RtreeNode::deleteChildNode(RtreeNode* _child) {
    if (this->listSize <= 0 || this->listPtr == nullptr) {
        printf("ERROR in RtreeNode::deleteChildNode: a null node\n");
        return listSize;
    }
    for (int i = 0; i < listSize; i++) {
        if (listPtr[i] == (char*)_child) {
            // Move the last child node to index i and remove the last one.
            listPtr[i] = listPtr[listSize - 1];
            listPtr[listSize - 1] = nullptr;
            listSize--;

            // delete this child node
            _child->releaseSpace();
            delete _child;
            _child = nullptr;
            break;
        }
    }
    return listSize;
}

void RtreeNode::initialize(bool _isLeaf, int _branchFactor, int _leafCapacity) {
    if (listPtr != nullptr) {
        printf("ERROR in RtreeNode::initialize: not an empty node\n");
    }
    listPtr = _isLeaf ? new char*[_leafCapacity] : new char*[_branchFactor];
    listSize = 0;
    // NOTE that the mbr will be initialized when adding the first element/child
}

RtreeNode* RtreeNode::chooseBestChildNode(Rectangle *_rect) {
    if(this->listSize <= 0 || this->listPtr == nullptr) {
        printf("ERROR in RtreeNode::chooseBestChildNode: an empty node\n");
        return nullptr;
    }
    int bestIndex = -1;
    double bestValue = MYINFTY, p;

    int childNodeNum = listSize;
    auto** childList = (RtreeNode**)listPtr;
    for (int i = 0; i < childNodeNum; ++i) {
        p = childList[i]->mbr->perimeterIncrement(*_rect);
        if (bestValue > p) {
            bestValue = p;
            bestIndex = i;
        }
        else if (bestValue == p && bestIndex >= 0) {
            double p1 = childList[i]->mbr->computePerimeter();
            double p2 = childList[bestIndex]->mbr->computePerimeter();
            if (p1 < p2) {
                // If the increment of perimeter are the same, then pick the one whose perimeter is smaller.
                bestValue = p;
                bestIndex = i;
            }
        }
    }
    return childList[bestIndex];
}

RtreeNode* RtreeNode::split(bool _isLeaf, char* _toInsPtr, int _dim, int _branchFactor, int _leafCapacity) {
    if(this->listSize <= 0 || this->listPtr == nullptr) {
        printf("ERROR in RtreeNode::split: an empty node\n");
        return nullptr;
    }

    int childNodeNum = listSize + 1;
    int half = childNodeNum / 2;
    int rest = childNodeNum - half;

    double bestValue = MYINFTY;
    int bestIndex = -1;

    char** tempList = new char*[childNodeNum];
    for (int i = 0; i < listSize; i++) {
        tempList[i] = listPtr[i];
    }
    tempList[listSize] = _toInsPtr;

    // Set the dimensionality of the rtree nodes.
    rtreeNode_dim = _dim;

    // Set the compare function pointer.
    bool(*comp)(char*, char*) = _isLeaf ? CompElementLess : CompRtreeNodeLess;

    // Find the best split axis.
    for (int i = 0; i < _dim; i++) {
        // sort the child nodes on each dimension coordinate.
        axisIndex = i;
        sort(tempList, tempList + childNodeNum, comp);
        double p1 = RtreeNode::computeGroupPerimeter(_isLeaf, tempList, 0, half, _dim);
        double p2 = RtreeNode::computeGroupPerimeter(_isLeaf, tempList, half, rest, _dim);
        if (p1 + p2 < bestValue) {
            // The smaller value the better.
            bestValue = p1 + p2;
            bestIndex = axisIndex;
        }
    }

    // Split this node by bestIndex.
    axisIndex = bestIndex;
    sort(tempList, tempList + childNodeNum, comp);

    // Create a new node u.
    auto* u = new RtreeNode;
    u->initialize(_isLeaf, _branchFactor, _leafCapacity);
    u->parent = this->parent;

    if (!_isLeaf) {
        // Reset this node.
        listSize = 0;
        for (int i = 0; i < half; ++i) {
            this->addChildNode((RtreeNode*)(tempList[i]), _dim);
        }
        for (int i = half; i < childNodeNum; ++i) {
            u->addChildNode((RtreeNode*)(tempList[i]), _dim);
        }
    }
    else {
        // Reset this leaf node.
        listSize = 0;
        for (int i = 0; i < half; ++i) {
            this->addElement((RtreeElement*)(tempList[i]), _dim);
        }
        for (int i = half; i < childNodeNum; ++i) {
            u->addElement((RtreeElement*)(tempList[i]), _dim);
        }
    }
    return u;
}

void RtreeNode::releaseSpace() {
    if(mbr != nullptr)
        mbr->releaseSpace();
    if (listPtr != nullptr) {
        delete[] listPtr;
        listPtr = nullptr;
    }
    listSize = 0;
}

double RtreeNode::computeGroupPerimeter(bool _isLeaf, char** _list, int _start, int _length, int _dim) {
    if (_start < 0 || _length < 1) {
        printf("Error in RtreeNode::computeGroupPerimeter: wrong index given\n");
        return -1;
    }

    auto * vList = new float [_dim * 2];
    float* minList = vList;
    float* maxList = &(vList[_dim]);

    if (!_isLeaf) {
        auto** childList = (RtreeNode**)_list;
        for (int k = 0; k < _dim * 2; k++) {
            vList[k] = childList[_start]->mbr->getValue(k);     // initialize
        }

        int end = _start + _length;
        for (int j = _start + 1; j < end; ++j) {
            for (int k = 0; k < _dim; ++k) {
                minList[k] = MIN(minList[k], childList[j]->mbr->getValue(k));
                maxList[k] = MAX(maxList[k], childList[j]->mbr->getValue(k + _dim));
            }
        }
    }
    else {      // This is a leaf node, containing rectangles
        auto** rectangleList = (Rectangle**)_list;
        for (int k = 0; k < _dim; k++) {
            vList[k] = rectangleList[_start]->getValue(k);     // initialize
        }

        int end = _start + _length;
        for (int j = _start + 1; j < end; ++j) {
            for (int k = 0; k < _dim; ++k) {
                minList[k] = MIN(minList[k], rectangleList[j]->getValue(k));
                maxList[k] = MAX(maxList[k], rectangleList[j]->getValue(k + _dim));
            }
        }
    }

    // Compute the perimeter of the mbr defined by minList and maxList.
    double perimeter = Rectangle::computePerimeter(minList[0], minList[1], maxList[0], maxList[1]);

    delete[] vList;

    return perimeter;
}

bool RtreeNode::temporalOverlap(TIMETYPE queryStartT, TIMETYPE queryEndT) const {
    return !(queryEndT <= minStartTime || queryStartT >= maxEndTime);
}
