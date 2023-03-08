//
// Created by Fengmei JIN on 2022/9/3.
//

#ifndef PCT_CATEGORY_H
#define PCT_CATEGORY_H


#include "GeoPoint.h"

struct POI {
    IDTYPE poiID{0};
    GeoPoint *point{nullptr};
    set<pair<string, int>> categoryIDs;     // a pair of subcategoryCode and categoryID

    POI() = default;

    POI(IDTYPE _pid, GeoPoint *_p, int _cid, const string &_scode): poiID(_pid), point(_p) {
        categoryIDs.insert(make_pair(_scode, _cid));
    };

    void addCategoryID(int _cid, const string &_scode) {
        categoryIDs.insert(make_pair(_scode, _cid));
    }
};

// describe POI categories
class Category {

    // the prefix of category codes
    vector<string> Food_0{"1"};
    vector<string> Mall_1{"2"};
    vector<string> Medical_2{"2800", "72", "75"};
    vector<string> Vehicle_3{"4", "7100", "7103"};
    vector<string> Entertainment_4{"5", "6"};    // 5 is hotel
    vector<string> Government_5{"70", "7101", "7102"};
    vector<string> Culture_6{"7180", "74", "94", "9600", "9700", "A880"};
    vector<string> Scenery_7{"73", "9080", "9180", "92", "93", "95"};
    vector<string> Business_8{"76", "AF"};
    vector<string> Residence_9{"7700", "A900"};
    vector<string> Life_10{"78", "A080", "A2", "A3", "A4", "A600", "AB"};
    vector<string> Suburbs_11{"808A", "808B", "9800", "BB", "BF"};    //BF00 is bridge
    vector<string> Transport_12{"8"}; // !!!except suburbs
    vector<string> Bank_13{"A1"};
    vector<string> School_14{"A70"};
    vector<string> Company_15{"A980", "A983", "AA"};
    vector<string> Factory_16{"A982"};

    vector<vector<string>> BeijingCategoryList;

public:
    Category(){
        BeijingCategoryList.emplace_back(Food_0);
        BeijingCategoryList.emplace_back(Mall_1);
        BeijingCategoryList.emplace_back(Medical_2);
        BeijingCategoryList.emplace_back(Vehicle_3);
        BeijingCategoryList.emplace_back(Entertainment_4);
        BeijingCategoryList.emplace_back(Government_5);
        BeijingCategoryList.emplace_back(Culture_6);
        BeijingCategoryList.emplace_back(Scenery_7);
        BeijingCategoryList.emplace_back(Business_8);
        BeijingCategoryList.emplace_back(Residence_9);
        BeijingCategoryList.emplace_back(Life_10);
        BeijingCategoryList.emplace_back(Suburbs_11);
        BeijingCategoryList.emplace_back(Transport_12);
        BeijingCategoryList.emplace_back(Bank_13);
        BeijingCategoryList.emplace_back(School_14);
        BeijingCategoryList.emplace_back(Company_15);
        BeijingCategoryList.emplace_back(Factory_16);
    }

    vector<string> getSubCategoryCodes(int cid) const {
        return BeijingCategoryList.at(cid);
    }

    unsigned int getCategoryNum() const {
        return BeijingCategoryList.size();
    }

    int getCategoryID(const string& subCategoryCode);

    void computeAggregatedDistribution(const vector<pair<POI *, double>> &poiList, float *rtn, bool involveDist) const;
};

#endif //PCT_CATEGORY_H
