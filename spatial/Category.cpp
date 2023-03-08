//
// Created by Fengmei JIN on 2022/9/20.
//


#include "Category.h"

/**
 * Compute an aggregated distribution for a given poi list
 * @param involveDist: the contribution of a poi to the aggregated distribution
 *                      should also be relevant to the distance between the poi and the ori-point
 * */
void Category::computeAggregatedDistribution(const vector<pair<POI *, double>> &poiList, float *rtn, bool involveDist) const {
    for(int i = 0; i < BeijingCategoryList.size(); i++) {
        rtn[i] = 0;
    }
    // find the minimum non-zero distance
    double minNonZeroDist = poiList.front().second;
    if(involveDist) {
        if(minNonZeroDist < 0)
            involveDist = false;
        if(minNonZeroDist == 0) {
            for (int i = 1; i < poiList.size(); i++) {
                double dist = poiList.at(i).second;
                if (dist > 0) {
                    minNonZeroDist = dist;
                    break;
                }
            }
            if(minNonZeroDist == 0)    // all distance are zero
                involveDist = false;    // no need to consider dist
            else
                minNonZeroDist = MIN(minNonZeroDist, 1);    // no larger than 1
        }
    }


    float sum = 0;
    for(const auto& poiPair : poiList) {
        double dist = poiPair.second > 0 ? poiPair.second : minNonZeroDist;

        auto categoryList = poiPair.first->categoryIDs;
        auto cNum = (float) categoryList.size();
        for(const auto &cpair: categoryList) {
            auto value = (float) (1 / cNum * (involveDist ? 1 / dist : 1));
            rtn[cpair.second] += value;
            sum += value;
        }
    }
    // normalized by the sum
    for(int i = 0; i < BeijingCategoryList.size(); i++) {
        if(rtn[i] > 0) {
            float a = rtn[i];
            float b = rtn[i] / sum;
            if(b <= 1 && b > 0) {
                rtn[i] = b;
            }
            else {
                // ===== debug
                printf("%d, %.3f, %.3f\n", i, a, b);
            }
        }
    }

}

int Category::getCategoryID(const string &subCategoryCode) {
    int cid = -1;
    for (int i = 0; i < BeijingCategoryList.size(); i++) {
        auto category = BeijingCategoryList.at(i);
        for(const auto &prefix: category) {
            if (prefix.length() == 4) {   // it's a whole code of the category instead of prefix
                if (subCategoryCode == prefix) {
                    cid = i;
                }
            }
            else if (subCategoryCode.substr(0, prefix.length()) == prefix) {
                cid = i;
            }
        }
    }
    if(cid < 0) {
        printf("Error in Category::getCategoryID - invalid category code %s\n", subCategoryCode.c_str());
    }
    else {
        BeijingCategoryList[cid].emplace_back(subCategoryCode);     // update the category list
    }
    return cid;
}
