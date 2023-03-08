//
// Created by Fengmei Jin on 6/7/2022.
//

#ifndef PCT_FILEREADER_H
#define PCT_FILEREADER_H

#include "../spatial/StaticGrid.h"
#include "../spatial/Trip.h"
#include "../utils.h"

class FileReader {

public:

    static void readParameterFile(const string &filename, const string& delim, map<string, string>& key2value);

    static string readPOIs(const string &filename, IndexPoint *pIndex, StaticGrid *poiGrid);

    static string readData(const string &filename, const int &dataTotal, IndexPoint *pIndex,
                           StaticGrid *poiGrid, bool needPoiDistribution, vector<Trip *> &trajectories);
};

#endif //PCT_FILEREADER_H
