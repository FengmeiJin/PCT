//
// Created by Fengmei Jin on 11/7/2022.
//

#ifndef PCT_PARAMETERS_H
#define PCT_PARAMETERS_H

#include "../commonheader.h"
#include "FileReader.h"

using namespace std;

class Parameters {
    string homePath, datasetName;
    int data_total{10};

    string privacy_model;
    float epsilon_geo;
    float epsilon_sem;

    bool independent_contact_prob{true};
    float contact_prob_threshold{0};
    TIMETYPE contact_time_threshold{60};    // unit sec
    float identical_poi_dist_threshold{0.1};      // unit km

    float initContactRatio{0.2};
    string searchType;
    int hopLength{1};

public:
    string readParameters(const string &filename);

    string getHomePath();
    string getDatasetName();

    int getObjectNum() const;

    string getPrivacyModel() const;
    float *getPrivacyBudget() const;

    // defining contact
    float getDistThreshold() const;
    TIMETYPE getContactTimeThreshold() const;
    float getContactMinProb() const;
    string getContactProbType() const;

    float getInitSpreaderRatio() const;
    string getSearchType() const;
    int getHopLength() const;
};


#endif //PCT_PARAMETERS_H
