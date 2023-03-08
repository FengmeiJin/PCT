//
// Created by Fengmei Jin on 11/7/2022.
//

#include "Parameters.h"

char *getCurrDirectory() {
    char *currDirectory;
    if ((currDirectory = getcwd(nullptr, 0)) == nullptr) {
        perror("getcwd error");
    }
    else {
        printf("Current Working Directory = %s\n", currDirectory);
    }
    return currDirectory;
}

string printout(const map<string, string> &key2value) {
    string msg = "***********************************************************************\n";
    msg.append("[Parameters]\n");
    for (auto &pair: key2value) {
        msg.append(pair.first);
        msg.append(" = ");
        msg.append(pair.second);
        msg.append("\n");
    }
    msg.append("***********************************************************************\n");
    return msg;
}

string Parameters::readParameters(const string &filename) {
    char *currDirectory = getCurrDirectory();

    string fullpath(currDirectory);
    fullpath = fullpath.substr(0, fullpath.find_last_of('/') + 1);
    fullpath.append(filename);    // fixed

    map<string, string> key2value;
    FileReader::readParameterFile(fullpath, " = ", key2value);

    for (auto &pair: key2value) {
        if (pair.first == "home_path") {
            homePath = pair.second;
        }
        else if (pair.first == "dataset_type") {
            datasetName = pair.second;
        }
        else if (pair.first == "data_total") {
            data_total = stoi(pair.second);
        }
        else if (pair.first == "identical_poi_dist_max") {
            identical_poi_dist_threshold = stof(pair.second);
        }
        else if (pair.first == "contact_min_duration") {
            contact_time_threshold = stol(pair.second);
        }
        else if (pair.first == "contact_min_prob") {
            contact_prob_threshold = stof(pair.second);
        }
        else if (pair.first == "independent_contact_prob") {
            independent_contact_prob = (pair.second == "true");
        }
        else if (pair.first == "privacy_model") {
            privacy_model = pair.second;
        }
        else if (pair.first == "privacy_epsilon_geo") {
            epsilon_geo = stof(pair.second);
        }
        else if (pair.first == "privacy_epsilon_semantic") {
            epsilon_sem = stof(pair.second);
        }
        else if (pair.first == "init_contact_ratio") {
            initContactRatio = stof(pair.second);
        }
        else if (pair.first == "search_type") {
            searchType = pair.second;
        }
        else if (pair.first == "multi_hop") {
            hopLength = stoi(pair.second);
        }
    }
    return printout(key2value);
}

string Parameters::getDatasetName() {
    return datasetName;
}

string Parameters::getHomePath() {
    return homePath;
}

int Parameters::getObjectNum() const {
    return data_total;
}

float Parameters::getDistThreshold() const {
    return identical_poi_dist_threshold;
}

TIMETYPE Parameters::getContactTimeThreshold() const {
    return contact_time_threshold;
}

float* Parameters::getPrivacyBudget() const {
    return new float[2]{epsilon_geo, epsilon_sem};
}

float Parameters::getInitSpreaderRatio() const {
    return initContactRatio;
}

string Parameters::getSearchType() const {
    return searchType;
}

string Parameters::getPrivacyModel() const {
    return privacy_model;
}

int Parameters::getHopLength() const{
    return hopLength;
}

float Parameters::getContactMinProb() const{
    return contact_prob_threshold;
}

string Parameters::getContactProbType() const {
    if (contact_prob_threshold == 0) {
        return "PCT";
    }
    return independent_contact_prob ? "PCTi" : "PCTc";
}
