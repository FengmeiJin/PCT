//
// Created by Fengmei JIN on 2022/8/23.
//

#include "RtreeIndex.h"

/*
*  The structure used in Rtree query algorithm.
*/
struct QueElement {
    RtreeNode *node{nullptr};
    int lev{0};

    // added by jinfm on 25/8/2022
    QueElement() = default;

    QueElement(RtreeNode *_node, int _lev) : node(_node), lev(_lev) {};

    bool isLeaf() const {
        return lev == 0;
    }
};

/***************************************************************************/

RtreeIndex::RtreeIndex() {
    root = nullptr;
    rootLev = 0;
    dim = -1;
    branchFactor = 12;
    leafCapacity = 12;
}

RtreeIndex::RtreeIndex(int _dim, int _branchFactor, int _leafCapacity) {
    root = nullptr;
    rootLev = 0;
    dim = _dim;
    branchFactor = _branchFactor;
    leafCapacity = _leafCapacity;
}

RtreeIndex::~RtreeIndex() {
    this->releaseSpace();
}

int RtreeIndex::insertNode(RtreeNode *_toInsNode, int _targetLevel) {
    if (root == nullptr) {
        if (dim == -1) {
            printf("Error in RtreeIndex::insertNode: dim = -1! Please initialize the rtree first!\n");
            return 1;
        }
        root = _toInsNode;
        return 0;
    }
    if (rootLev == _targetLevel) {
        // In this case, we need to create a new root.
        auto *newRoot = new RtreeNode;
        newRoot->initialize(false, branchFactor, leafCapacity);
        newRoot->addChildNode(this->root, dim);
        newRoot->addChildNode(_toInsNode, dim);
        this->root = newRoot;
        rootLev = _targetLevel + 1;
    }
    else {
        this->insertNode_SubRoutine(this->root, rootLev, _toInsNode, _targetLevel);
    }
    return 0;
}

int RtreeIndex::insertElement(RtreeElement *_elem) {
    if (root == nullptr) {
        if (dim == -1) {
            printf("Error in RtreeIndex::insertElement: dim = -1! Please initialize the rtree first!\n");
            return 1;
        }

        // Create a new node as the root and insert the element into it
        root = new RtreeNode;
        root->initialize(true, branchFactor, leafCapacity);
        root->addElement(_elem, dim);
        return 0;
    }

    if (rootLev == 0) { // the current root is the leaf
        if (root->getListSize() + 1 > leafCapacity) {
            // the root will overflow after inserting this element.
            RtreeNode *newNode = root->split(true, (char *) _elem, dim, branchFactor, leafCapacity);
            auto *newRoot = new RtreeNode;
            newRoot->initialize(false, branchFactor, leafCapacity);
            newRoot->addChildNode(root, dim);
            newRoot->addChildNode(newNode, dim);
            root = newRoot;
            rootLev++;
        }
        else {
            root->addElement(_elem, dim);
        }
        return 0;
    }

    return this->insertElement_SubRoutine(root, rootLev, _elem);
}

/********************************* Query Functions ******************************************/

/*
*  Added by jhgan on 2017-01-23.
*  The implementation of range reporting.
*/
int RtreeIndex::rangeQuery(GeoPoint *queryP, double _radius, vector<Interval *> &_targetPlace,
                           bool enableRtreeTimeFilter, TIMETYPE queryStartT, TIMETYPE queryEndT) {
    if (root == nullptr) {
        printf("Warning in RtreeIndex::rangeQuery: The root is nullptr!\n");
        return 0;
    }
    // Clear the targetPlace.
    _targetPlace.clear();

    RtreeNode **childList = nullptr;
    int childNum = 0;
    RtreeElement **elemList = nullptr;
    int elemNum = 0;
    int state = 0;

    queue<QueElement> que;
    QueElement temp = QueElement();
    QueElement curElement = QueElement(root, rootLev);
    que.push(curElement);
    while (!que.empty()) {
        curElement = que.front();
        que.pop();

        if (curElement.isLeaf()) {
            // the current node is a leaf node.
            elemList = (RtreeElement **)(curElement.node->getListPtr());
            elemNum = curElement.node->getListSize();
            for (int i = 0; i < elemNum; i++) {
                if (enableRtreeTimeFilter && elemList[i]->getInterval()->existTemporalOverlap(queryStartT, queryEndT)
                    || !enableRtreeTimeFilter) {
                    double dist = elemList[i]->computeDistToPoint(queryP);
                    if (dist <= _radius + elemList[i]->getRadius())
                        _targetPlace.emplace_back(elemList[i]->getInterval());
                }
            }
        }
        else {
            // the current node is an internal node.
            childList = (RtreeNode **) (curElement.node->getListPtr());
            childNum = curElement.node->getListSize();
            for (int i = 0; i < childNum; i++) {
                if (enableRtreeTimeFilter && childList[i]->temporalOverlap(queryStartT, queryEndT)
                    || !enableRtreeTimeFilter) {
                    state = childList[i]->mbr->stateWithSphere(queryP, _radius, dim);
                    if (state != -1) {
                        // intersect or fully inside, then add this node to queue.
                        temp.node = childList[i];
                        temp.lev = curElement.lev - 1;
                        que.push(temp);
                    }
                }
            }
        }
    }
    return 0;
}

void RtreeIndex::releaseSpace() {
    if (root != nullptr) {
        releaseSpace_SubRoutine(root, rootLev);
        delete root;
        root = nullptr;
        dim = -1;
    }
}

/*
*  Private functions.
*/
int RtreeIndex::insertElement_SubRoutine(RtreeNode *_subroot, int _subrootLevel, RtreeElement *_elem) {

    if (_subroot == nullptr || _elem == nullptr) {
        printf("Error in RtreeIndex::insertElement_SubRoutine: _subroot or _elem is nullptr!\n");
        return 1;
    }

    // Version before Jan 8, 2017
    if (_subrootLevel == 0) {
        // _subroot is a leaf node.
        if (_subroot->getListSize() + 1 > leafCapacity) {
            // the _subroot will overflow after inserting this element.
            RtreeNode *newNode = _subroot->split(true, (char *) _elem, dim, branchFactor, leafCapacity);
            this->insertNode_SubRoutine(_subroot->parent, 1, newNode, 0);
        }
        else {
            _subroot->addElement(_elem, dim);
        }
        return 0;
    }

    // subroot is an internal node.
    _subroot->mbr->enlarge(_elem, dim);
    RtreeNode *bestChild = _subroot->chooseBestChildNode(_elem);
    return this->insertElement_SubRoutine(bestChild, _subrootLevel - 1, _elem);
}

void RtreeIndex::insertNode_SubRoutine(RtreeNode *_subroot, int _subrootLevel, RtreeNode *_toInsNode, int _targetLevel) {
    if (_subroot == nullptr) {
        printf("Error in RtreeIndex::InsertNode_SubRoutine: subroot is nullptr!\n");
        return;
    }
    if (_toInsNode == nullptr) {
        printf("Error in RtreeIndex::InsertNode_SubRoutine: toInsNode is nullptr!\n");
        return;
    }

    if (_targetLevel >= _subrootLevel) {
        printf("Error in RtreeIndex::InsertNode_SubRoutine: targetLevel >= subroot->level!\n");
        return;
    }

    // Version before 8 Jan, 2017.
    if (_subrootLevel == _targetLevel + 1) {
        if (_subroot->getListSize() + 1 > branchFactor) {
            // this subroot overflows after this insertion. It needs to split.
            this->overflowHandler(_subroot, _toInsNode);
        }
        else {
            _subroot->addChildNode(_toInsNode, dim);
        }
    }
    else {
        // Enlarge the mbr of subroot.
        _subroot->mbr->enlarge(_toInsNode->mbr, dim);
        RtreeNode *bestChild = _subroot->chooseBestChildNode(_toInsNode->mbr);
        this->insertNode_SubRoutine(bestChild, _subrootLevel - 1, _toInsNode, _targetLevel);
    }
}

void RtreeIndex::overflowHandler(RtreeNode *_overflowNode, RtreeNode *_toInsNode) {
    RtreeNode *newNode = _overflowNode->split(false, (char *) _toInsNode, dim, branchFactor, leafCapacity);

    if (_overflowNode->parent != nullptr) {
        if (_overflowNode->parent->getListSize() + 1 > branchFactor) {
            // its parent node overflows after inserting newNode.
            // recursively handle this case
            this->overflowHandler(_overflowNode->parent, newNode);
        }
        else {
            _overflowNode->parent->addChildNode(newNode, dim);
        }
    }
    else {
        // _overflowNode is the root of this tree.
        // then create a new root for this tree.
        auto *newRoot = new RtreeNode;
        newRoot->initialize(false, branchFactor, leafCapacity);
        newRoot->addChildNode(_overflowNode, dim);
        newRoot->addChildNode(newNode, dim);
        this->root = newRoot;
        rootLev++;
    }
}

void RtreeIndex::underflowHandler(RtreeNode *_underflowNode, int _nodeLevel) {
    if (_underflowNode == nullptr) {
        printf("Error in RtreeIndex::underflowHandler: a null _underflowNode!\n");
        return;
    }
    RtreeNode *parent = _underflowNode->parent;
    if (parent == nullptr) {
        // _underflowNode is the root of this tree.
        int num = _underflowNode->getListSize();
        if (num > 1) {
            // _underflowNode is a legal root.
            return;
        }
        else {
            if (num == 1) {
                if (rootLev > 0) {
                    // _underflowNode is an illegal root.
                    // Let the only child node be the root.
                    this->root = (RtreeNode *) (root->getListPtr()[0]);
                    this->root->parent = nullptr;
                    rootLev = rootLev - 1;
                    _underflowNode->releaseSpace();
                    delete _underflowNode;
                    _underflowNode = nullptr;
                }
            }
            else {
                this->root = nullptr;
                rootLev = 0;
                _underflowNode->releaseSpace();
                delete _underflowNode;
                _underflowNode = nullptr;
            }
        }
    }
    else {
        // _underflowNode is not the root of this tree.

        // Get all the child nodes of _underflowNode.
        int reInsNum = _underflowNode->getListSize();
        char **listPtr = _underflowNode->getListPtr();

        if (reInsNum > 0 && listPtr != nullptr) {
            char **reInsList = new char *[reInsNum];
            for (int i = 0; i < reInsNum; i++)
                reInsList[i] = listPtr[i];

            // Delete underflowNode from its parent.
            if (parent->deleteChildNode(_underflowNode) < branchFactor / FanoutRatio) {
                // parent underflows after deleting underflowNode.
                // Handle the underflowing event of parent.
                this->underflowHandler(parent, _nodeLevel + 1);
            }
            if (_nodeLevel > 0) {
                // Reinsert all the child nodes of underflowNode.
                for (int i = 0; i < reInsNum; i++) {
                    this->insertNode((RtreeNode *) (reInsList[i]), _nodeLevel - 1);
                }
            }
            else {
                // Reinsert all the elements of underflowNode.
                for (int i = 0; i < reInsNum; i++) {
                    this->insertElement((RtreeElement *) (reInsList[i]));
                }
            }

            // added by jinfm on 25/8/2022
            for (int i = 0; i < reInsNum; i++) {
                reInsList[i] = nullptr;
                delete reInsList[i];
            }
            delete[] reInsList;
        }
    }
}

void RtreeIndex::releaseSpace_SubRoutine(RtreeNode *_subroot, int _subrootLev) {
    char **childList = _subroot->getListPtr();
    int childNum = _subroot->getListSize();
    RtreeNode *curChildNode = nullptr;
    if (_subrootLev > 0) {
        for (int i = 0; i < childNum; i++) {
            curChildNode = (RtreeNode *) (childList[i]);
            releaseSpace_SubRoutine(curChildNode, _subrootLev - 1);
            delete curChildNode;
        }
    }
    _subroot->releaseSpace();
}


/***************************************************
*  STR bulk loading related.
**************************************************/
/**
*  Construct leaf level by STR algorithm.
*  @param	_leafLevel:			the target place to store the pointers of the leaf nodes.
*  @param	_dimLev:			which dimension we used for sorting in the current recursion.
*/
int RtreeIndex::constructLeafLevel_STR_SubRoutine(RtreeElement **_elementList, int _elemNum, int _dim,
                                                  int _branchFactor, int _leafCapacity, vector<RtreeNode *> &_leafLevel, int _dimLev) {
    if (_elementList == nullptr) {
        printf("Error in RtreeIndex::constructLeafLevel_STR_SubRoutine: _elementList == nullptr!\n");
        return 1;
    }
    if (_elemNum <= 0) {
        printf("Error in RtreeIndex::constructLeafLevel_STR_SubRoutine: _elemNum <= 0!\n");
        return 1;
    }

    double exp = 1.0 / (_dim - _dimLev + 1);
    int slabNum = ceil(pow((_elemNum * 1.0 / _leafCapacity), exp));
    int slabSize = ceil(_elemNum * 1.0 / slabNum);  // # elements in each slab.

    if (_elemNum <= _leafCapacity) {
        // There is only one node.
        auto *leafNode = new RtreeNode;
        leafNode->initialize(true, _branchFactor, _leafCapacity);
        for (int i = 0; i < _elemNum; ++i) {
            leafNode->addElement(_elementList[i], _dim);
        }
        _leafLevel.emplace_back(leafNode);
        return 0;
    }

    if (slabSize <= _leafCapacity / FanoutRatio || _dimLev == _dim) {
        // There are too few points to make slabs. Or this is the last level in the recursion.
        // Sort all the points by the _dimLev-dim coordinates.
        sort(_elementList, _elementList + _elemNum, RectangleDimSorter(_dim, _dimLev));

        // Scan the sorted list and construct tree nodes.
        int nodeNum = ceil(_elemNum * 1.0 / _leafCapacity); // Note that nodeNum > 1.
        RtreeNode *curNode = nullptr;

        int curIndex = 0;
        int tempNum = _leafCapacity;
        // We construct the first (nodeNum - 2) nodes because we may need to deal with the last two nodes specially.
        for (int i = 0; i < nodeNum - 2; ++i) {
            curNode = new RtreeNode;
            curNode->initialize(true, _branchFactor, _leafCapacity);
            for (int j = 0; j < tempNum; ++j) {
                curNode->addElement(_elementList[curIndex], _dim);
                curIndex++;
            }
            _leafLevel.emplace_back(curNode);
        }
        // We deal with the last two nodes specially.
        tempNum = (_elemNum - curIndex) / 2;
        curNode = new RtreeNode;
        curNode->initialize(true, _branchFactor, _leafCapacity);
        for (int j = 0; j < tempNum; ++j) {
            curNode->addElement(_elementList[curIndex], _dim);
            curIndex++;
        }
        _leafLevel.emplace_back(curNode);

        tempNum = _elemNum - curIndex;
        curNode = new RtreeNode;
        curNode->initialize(true, _branchFactor, _leafCapacity);
        for (int j = 0; j < tempNum; ++j) {
            curNode->addElement(_elementList[curIndex], _dim);
            curIndex++;
        }
        _leafLevel.emplace_back(curNode);
    }
    else {
        // This is not the last level in the recursion and there are enough points in each slab.
        // Make slabs.
        // Note that: There must be at least two slabs because _leafCapacity < _elemNum.
        // Sort all the points by the _dimLev-dim coordinates.
        sort(_elementList, _elementList + _elemNum, RectangleDimSorter(_dim, _dimLev));
        int rtn = 0;
        RtreeElement **curPtr = _elementList;
        for (int i = 0; i < slabNum - 2; ++i) {
            rtn = constructLeafLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _leafLevel,_dimLev + 1);
            if (rtn == 1) {
                printf("Error in RtreeIndex::constructLeafLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
                return 1;
            }
            curPtr += slabSize;
        }

        // Process the last two slabs specially.
        int gapPtrSize = (int) (curPtr - _elementList);
        slabSize = (_elemNum - gapPtrSize) / 2;
        rtn = constructLeafLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _leafLevel,_dimLev + 1);
        if (rtn == 1) {
            printf("Error in RtreeIndex::constructLeafLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
            return 1;
        }
        curPtr += slabSize;

        // Put all the rest points to the last slab.
        gapPtrSize = (int) (curPtr - _elementList);
        slabSize = _elemNum - gapPtrSize;
        rtn = constructLeafLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _leafLevel,_dimLev + 1);
        if (rtn == 1) {
            printf("Error in RtreeIndex::constructLeafLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
            return 1;
        }
        curPtr += slabSize;
    }
    return 0;
}

int RtreeIndex::constructNextLevel_STR_SubRoutine(RtreeNode **_nodeList, int _nodeNum, int _dim,
                                                  int _branchFactor, int _leafCapacity, vector<RtreeNode *> &_nextLevel, int _dimLev) {
    if (_nodeList == nullptr) {
        printf("Error in RtreeIndex::constructNextLevel_STR_SubRoutine: _nodeList == nullptr!\n");
        return 1;
    }
    if (_nodeNum <= 0) {
        printf("Error in RtreeIndex::constructNextLevel_STR_SubRoutine: _nodeNum <= 0!\n");
        return 1;
    }

    double exp = 1.0 / (_dim - _dimLev + 1);
    int slabNum = ceil(pow((_nodeNum * 1.0 / _branchFactor), exp));
    int slabSize = ceil(_nodeNum * 1.0 / slabNum);  // # nodes in each slab.

    if (_nodeNum <= _branchFactor) {
        // There is only one node.
        auto *curNode = new RtreeNode;
        curNode->initialize(false, _branchFactor, _leafCapacity);

        for (int i = 0; i < _nodeNum; ++i) {
            curNode->addChildNode(_nodeList[i], _dim);
        }
        _nextLevel.emplace_back(curNode);
        return 0;
    }

    if (slabSize <= _branchFactor / FanoutRatio || _dimLev == _dim) {
        // There are too few points to make slabs. Or this is the last level in the recursion.
        // Sort all the points by the _dimLev-dim coordinates.
        sort(_nodeList, _nodeList + _nodeNum, NodeDimSorter(_dim, _dimLev));
        // Scan the sorted list and construct tree nodes.
        int nextNodeNum = ceil(_nodeNum * 1.0 / _branchFactor); // Note that nextNodeNum > 1.
        RtreeNode *curNode = nullptr;

        int curIndex = 0;
        int tempNum = _branchFactor;
        // We construct the first (nextNodeNum - 2) nodes because we may need to deal with the last two nodes specially.
        for (int i = 0; i < nextNodeNum - 2; ++i) {
            curNode = new RtreeNode;
            curNode->initialize(false, _branchFactor, _leafCapacity);
            for (int j = 0; j < tempNum; ++j) {
                curNode->addChildNode(_nodeList[curIndex], dim);
                curIndex++;
            }
            _nextLevel.emplace_back(curNode);
        }
        // We deal with the last two nodes specially.
        tempNum = (_nodeNum - curIndex) / 2;
        curNode = new RtreeNode;
        curNode->initialize(false, _branchFactor, _leafCapacity);

        for (int j = 0; j < tempNum; ++j) {
            curNode->addChildNode(_nodeList[curIndex], dim);
            curIndex++;
        }
        _nextLevel.emplace_back(curNode);

        tempNum = _nodeNum - curIndex;
        curNode = new RtreeNode;
        curNode->initialize(false, _branchFactor, _leafCapacity);

        for (int j = 0; j < tempNum; ++j) {
            curNode->addChildNode(_nodeList[curIndex], _dim);
            curIndex++;
        }
        _nextLevel.emplace_back(curNode);
    }
    else {
        // This is not the last level in the recursion and there are enough points in each slab.
        // Make slabs.
        //
        // Note that: There must be at least two slabs because _branchFactor < _nodeNum.
        // Sort all the points by the _dimLev-dim coordinates.
        sort(_nodeList, _nodeList + _nodeNum, NodeDimSorter(_dim, _dimLev));
        int rtn = 0;
        RtreeNode **curPtr = _nodeList;
        for (int i = 0; i < slabNum - 2; ++i) {
            rtn = constructNextLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _nextLevel, _dimLev + 1);
            if (rtn == 1) {
                printf("Error in RtreeIndex::constructNextLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
                return 1;
            }
            curPtr += slabSize;
        }
        // Process the last two slabs specially.
        int gapPtrSize = (int) (curPtr - _nodeList);
        slabSize = (_nodeNum - gapPtrSize) / 2;
        rtn = constructNextLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _nextLevel, _dimLev + 1);
        if (rtn == 1) {
            printf("Error in RtreeIndex::constructNextLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
            return 1;
        }
        curPtr += slabSize;

        // Put all the rest points to the last slab.
        gapPtrSize = (int) (curPtr - _nodeList);
        slabSize = _nodeNum - gapPtrSize;
        rtn = constructNextLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _nextLevel, _dimLev + 1);
        if (rtn == 1) {
            printf("Error in RtreeIndex::constructNextLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
            return 1;
        }
        curPtr += slabSize;
    }
    return 0;
}

/*
* 	The STR bulk loading algorithm.
*/
RtreeNode *RtreeIndex::constructRtree_STR(RtreeElement **_elementList, int _elemNum, int _dim, int _branchFactor, int _leafCapacity) {
    if (root != nullptr) {
        printf("Error in RtreeIndex::constructRtree_STR: The tree is not empty!\n");
        return nullptr;
    }
    if (_branchFactor < FanoutRatio || _leafCapacity < FanoutRatio) {
        printf("Error in RtreeIndex::constructRtree_STR: _branchFactor and _leafCapacity should be >= %d!\n",
               FanoutRatio);
        return nullptr;
    }

    dim = _dim;
    branchFactor = _branchFactor;
    leafCapacity = _leafCapacity;
    rootLev = 0;

    if (_elementList == nullptr) {
        printf("Error in RtreeIndex::constructRtree_STR: _elementList == nullptr!\n");
        return nullptr;
    }
    if (_elemNum <= 0) {
        printf("Error in RtreeIndex::constructRtree_STR: _elemNum <= 0!\n");
        return nullptr;
    }

    int rtn = 0;
    vector<RtreeNode *> curNodeList;
    rtn = constructLeafLevel_STR_SubRoutine(_elementList, _elemNum, _dim, _branchFactor, _leafCapacity, curNodeList, 1);
    if (rtn == 1) {
        printf("Error in RtreeIndex::constructRtree_STR: Leaf level construction failed!\n");
        curNodeList.clear();
        return nullptr;
    }

    int level = 0;
    vector<RtreeNode *> nextNodeList;
    RtreeNode **curNodeListPtr = nullptr;
    int curNodeNum = 0;
    // Start to construct internal levels.
    while (curNodeList.size() > 1) {
        nextNodeList.clear();
        curNodeListPtr = curNodeList.data();
        curNodeNum = (int) curNodeList.size();
        rtn = constructNextLevel_STR_SubRoutine(curNodeListPtr, curNodeNum, _dim,_branchFactor, _leafCapacity,nextNodeList, 1);
        if (rtn == 1) {
            printf("Error in RtreeIndex::constructRtee_STR: Next level construction failed!\n");
            curNodeList.clear();
            nextNodeList.clear();
            return nullptr;
        }
        curNodeList.swap(nextNodeList);
        level++;
    }
    root = curNodeList[0];
    rootLev = level;

    curNodeList.clear();
    nextNodeList.clear();
    return root;
}

int RtreeIndex::getHeight() const {
    return rootLev + 1;
}