//
// Created by Fengmei Jin on 6/7/2022.
//

#include "FileReader.h"
#include <sstream>

void FileReader::readParameterFile(const string &filename, const string &delim, map<string, string> &key2value) {
    ifstream fin(filename);
    if(!fin.good()) {
        printf("Error in FileReader::readParameterFile: cannot open %s !\n", filename.c_str());
        return;
    }

    string line, paramName, paramValue;
    while (fin.good()) {
        getline(fin, line);
        if (!line.empty() && line[0] != '#') {
            vector<string> vec;
            split(line, vec, delim);
            paramName = vec.at(0);
            paramValue = vec.at(1);
            key2value[paramName] = paramValue;  // to output all parameters
        }
    }
    fin.close();
}

// format: name,categoryCode,address,,longitude,latitude
string FileReader::readPOIs(const string &filename, IndexPoint *pIndex, StaticGrid *poiGrid) {
    ifstream fin(filename);
    if(!fin.good()) {
        return "[ERROR] cannot open the POI file in " + filename + "\n";
    }

    const int codePos = 1, lngPos = 4, latPos = 5;

    vector<tuple<string, float, float>> poiList;
    string line;
    bool firstLine = true;
    vector<string> tokens;
    while (getline(fin, line)) {     // parse one line indicating one POI
        if(firstLine) {
            firstLine = false;  // the header of the file
            continue;
        }
        tokens.clear();
        split(line, tokens);
        poiList.emplace_back(make_tuple(tokens[codePos], stof(tokens[lngPos]), stof(tokens[latPos])));
    }

    if(poiGrid == nullptr) {
        poiGrid = new StaticGrid();
    }
    poiGrid->initialize(poiList, pIndex);
    char *buffer = new char [100];
    snprintf(buffer, 100, "[INFO] # raw POIs = %lu, # valid POIs in Beijing = %d, # valid grid cells = %d\n",
             poiList.size(), poiGrid->getValidPoiNum(), poiGrid->getValidGridNum());
    return buffer;
}


string FileReader::readData(const string &filename, const int &dataTotal, IndexPoint *pIndex,
                            StaticGrid *poiGrid, bool needPoiDistribution, vector<Trip *> &trajectories) {
    ifstream fin(filename);
    if (!fin.good()) {
        return "[ERROR] cannot open the data file in " + filename + "\n";
    }

    // format: seqID, poiID, start_timestamp, end_timestamp, start_time, end_time, travel_time (s), travel_time (min), duration (s), duration (min)
    const int seqIdPos = 0, poiIdPos = 1, stPos = 2, etPos = 3;

    unsigned int cnt = 0, poiNum = 0;
    vector<tuple<IDTYPE, TIMETYPE, TIMETYPE>> poiSeq;

    printf("[PROGRESS] ");
    fflush(stdout);

    string line;
    bool firstLine = true;
    vector<string> tokens;
    string currId, prevId;
    while (getline(fin, line)) {     // parse one line indicating one point
        if(firstLine) {
            firstLine = false;  // the header of the file
            continue;
        }
        tokens.clear();
        split(line, tokens);
        currId = tokens.at(seqIdPos);
        if (!prevId.empty() && currId != prevId) {
            poiNum += poiSeq.size();

            trajectories.emplace_back(new Trip(prevId, poiSeq, pIndex, poiGrid, needPoiDistribution));

            // if total == 0, then read all data
            if (dataTotal > 0) {
                // print out progress
                if (dataTotal >= 10 && ++cnt % (dataTotal / 10) == 0) {
                    printf("%d%%...", (int) (cnt * 100.0 / dataTotal));
                    fflush(stdout);
                }
                if(trajectories.size() >= dataTotal) {
                    break;
                }
            }
            poiSeq.clear();
            poiSeq.shrink_to_fit();
        }

        prevId = currId;

        IDTYPE poiId = stoi(tokens.at(poiIdPos));
        TIMETYPE startTimestamp = stol(tokens.at(stPos));
        TIMETYPE endTimestamp = stol(tokens.at(etPos));
        poiSeq.emplace_back(make_tuple(poiId, startTimestamp, endTimestamp));
    }
    fin.close();

    // the last one
    if (trajectories.size() < dataTotal) {
        poiNum += poiSeq.size();

        trajectories.emplace_back(new Trip(prevId, poiSeq, pIndex, poiGrid, needPoiDistribution));
    }
    printf("done!\n");
    fflush(stdout);

    poiSeq.clear();
    poiSeq.shrink_to_fit();

    sort(trajectories.begin(), trajectories.end(), sortByStartTime);

    // some statistics message
    string msg = "\t Total Trajectories = ";
    msg.append(to_string(trajectories.size()));
    msg.append(", AVG length of simulated checkin sequences = ");
    msg.append(to_string((float) poiNum / (float) trajectories.size()) + "\n");

    return msg;
}
