//
// Created by Fengmei Jin on 11/7/2022.
//

#ifndef PCT_LOGFILE_H
#define PCT_LOGFILE_H

#include "../commonheader.h"

using namespace std;

string privacyParameterStr(const string &privacyType, const float *epsilon);

class LogFile {

    string logFilepath;

public:
    static bool testOrCreateDir(const string &pathname);

    explicit LogFile(string filepath);

    void addContent(string &content, bool screenprint = true);

    static string composeFilename(const string &datasetType, int dataTotal, int downSample, bool ANONYMIZE_DATA,
                                  float identicalPoiDist, TIMETYPE contactMinDuration, float contactMinProb, bool independentCP,
                                  int multiHop, const string &searchType, float initContactRatio,
                                  const string &privacyType, const float *epsilon, bool FIXSEED);

};


#endif //PCT_LOGFILE_H
