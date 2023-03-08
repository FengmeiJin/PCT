//
// Created by Fengmei Jin on 11/7/2022.
//

#include "LogFile.h"
#include "../utils.h"

bool LogFile::testOrCreateDir(const string& pathname){
    size_t found = pathname.find_last_of('/');
    string folder = pathname.substr(0, found + 1);
    ifstream fout (folder);
    if(!fout.good()){
        string commend = "mkdir -p " + folder;
        system(commend.c_str());
        fout.open(folder.c_str(), ios_base::in);
        return fout.good();
    }
    return true;
}

LogFile::LogFile(string filepath) {
    logFilepath = std::move(filepath);
    string msg;
    if(!testOrCreateDir(logFilepath)){
        msg = "[ERROR] Cannot create a directory for the log file " + logFilepath + " !\n";
    } else {
        msg = "[INFO] The log file is " + logFilepath + " !\n\n";
    }
    addContent(msg);
}

void LogFile::addContent(string& content, bool screenprint) {
    ofstream outfile;
    outfile.open(logFilepath, ios_base::app | ios_base::out);
    if (outfile.is_open()) {
        outfile << content << endl;
        outfile.close();
    }
    else {
        cout << "[ERROR] Cannot open the log file " << logFilepath << " !\n";
    }
    if(screenprint){
        printf("%s", content.c_str());
        fflush(stdout);
    }
    content.clear();
}

string privacyParameterStr(const string &privacyType, const float *epsilon){
    string str;

    // privacy parameters
    if(privacyType.size() >= 3) {
        if(privacyType == "geo") {
            str.append("_G");
            str.append("_EG" + float2str(epsilon[0]));
        }
        else if (privacyType == "geosemantic") {
            str.append("_GS");
            str.append("_EG" + float2str(epsilon[0]));
            str.append("_ES" + float2str(epsilon[1]));
        }
        else if (privacyType == "semantic") {
            str.append("_S");
            str.append("_ES" + float2str(epsilon[1]));
        }
        else {
            str.append("_" + privacyType);
        }
    }
    return str;
}

string LogFile::composeFilename(const string &datasetType, int dataTotal, int downSample, bool ANONYMIZE_DATA,
                                float identicalPoiDist, TIMETYPE contactMinDuration, float contactMinProb, bool independentCP,
                                int multiHop, const string &searchType, float initContactRatio,
                                const string &privacyType, const float *epsilon, bool FIXSEED) {

    // data parameters
    string filename = datasetType + "_D" + to_string(dataTotal);
    if (datasetType == "Geolife") {
        filename.append("_DS" + to_string(downSample));
    }
    filename.append("_QR" + float2str(initContactRatio));

    // privacy parameters
    filename.append(privacyParameterStr(privacyType, epsilon));

    if(FIXSEED) {
        filename.append("_FS");
    }

    // contact parameters
    filename.append("_MH" + to_string(multiHop));
    if(contactMinProb > 0 && ANONYMIZE_DATA) {
        filename.append(independentCP ? "_I" : "_C");
        filename.append("CP" + float2str(contactMinProb));
    }
    filename.append("_CD" + float2str(identicalPoiDist));
    filename.append("_CT" + to_string(contactMinDuration));

    // search parameters
    filename.append("_" + searchType);

    return filename;
}
