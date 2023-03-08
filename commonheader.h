//
// Created by Fengmei Jin on 11/7/2022.
//

#ifndef PCT_COMMONHEADER_H
#define PCT_COMMONHEADER_H

#include <fstream>
#include <iostream>

#include <stdexcept>
#include <unistd.h>
#include <dirent.h>

#include <string>
#include <utility>
#include <cmath>

// container
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <vector>
#include <queue>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

typedef unsigned int IDTYPE;  // the type of an ID (point, traj, ...)
                        // the available ID starts from 1 in case the type is unsigned

typedef unsigned long TIMETYPE;     // the type of timestamp (unit: second)

// beijing range
// latitude [38, 43]; longitude [114, 119]
#define BJ_LNG_MIN 114
#define BJ_LNG_MAX 119
#define BJ_LAT_MIN 38
#define BJ_LAT_MAX 43

const TIMETYPE EarliestStartTime = 1176307200;     // aka 2007-04-12 00:00:00, the earliest start time in Geolife
const int DeltaT = 3;       // for IR-tree index

#endif //PCT_COMMONHEADER_H
