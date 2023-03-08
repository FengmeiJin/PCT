//
// Created by Fengmei JIN on 2022/9/26.
//

#include "ItreeIndex.h"

void ItreeIndex::constructIndexByIntervals(const vector<Interval *> &sortedIntervals,
                                           IndexPoint *pIndex) {
    if (this->root == nullptr) {
        printf("ERROR in ItreeIndex::constructIndexByIntervals - null root!\n");
        return;
    }
    for (Interval *iv: sortedIntervals) {
        int rtn = this->root->insertInterval(iv);
        if (rtn != 0) {
            printf("ERROR in ItreeIndex::constructIndexByIntervals - cannot insert the interval (%lu, %lu)!\n", iv->getStartTime(), iv->getEndTime());
        }
    }
    this->root->organizeIntervals(pIndex);
}

void ItreeIndex::constructEmptyTree(TIMETYPE earliestIndexTime, TIMETYPE LatestEndTime, TIMETYPE contactMinDuration, int deltaT) {
    earliestIndexTime = (earliestIndexTime / 100) * 100;
    LatestEndTime = (LatestEndTime / 100 + 1) * 100;
    const TIMETYPE minBinSize = contactMinDuration * deltaT;
    const int minBinNum = (int) ((LatestEndTime - earliestIndexTime) / minBinSize) + 1;
    vector<ItreeNode*> nodeList;
    for (int i = 0; i < minBinNum; i++) {
        auto *node = new ItreeNode(earliestIndexTime + i * minBinSize, earliestIndexTime + (i + 1) * minBinSize);
        nodeList.emplace_back(node);   // leaf nodes
    }

    vector<ItreeNode*> parentList;
    while (nodeList.size() > 1) {
        unsigned int n = nodeList.size();
        // set next pointer for each node
        for (unsigned int i = 0; i < n - 1; i++) {
            nodeList[i]->setPointer(nodeList[i + 1]);           // next pointer
        }
        bool lastRemain = (n % 2 != 0);
        if (lastRemain) {
            n--;
        }
        for (unsigned int i = 0; i < n - 1; i += 2) {
            ItreeNode *node;
            if (n == 2) {
                this->root = new ItreeNode(nodeList[i]->getMinStartTime(), nodeList[i + 1]->getMaxEndTime());
                node = this->root;
            }
            else {
                node = new ItreeNode(nodeList[i]->getMinStartTime(), nodeList[i + 1]->getMaxEndTime());
            }
            node->setChildNode(nodeList[i]);                    // left child
            node->setChildNode(nodeList[i + 1], false);     // right child
            parentList.emplace_back(node);
        }
        if (lastRemain) {
            parentList.emplace_back(nodeList.back());
        }
        nodeList = parentList;
        parentList.clear();
        parentList.shrink_to_fit();
    }
}
