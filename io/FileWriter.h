//
// Created by Fengmei JIN on 2022/7/18.
//

#ifndef PCT_FILEWRITER_H
#define PCT_FILEWRITER_H

#include "../spatial/Trip.h"
#include "../tracing/functions.h"
#include <sstream>

class FileWriter {

public:

    static void outputContacts(const string &outputFile, int multiHop, bool enableCP,
                               const vector<pair<string, ContactEvent *>> &susceptible2event);

    static void outputNoisyTrips(const string &outputFile, IndexPoint *pIndex, const vector<Trip *> &noisyTrips);
};


#endif //PCT_FILEWRITER_H
