//
// Created by Fengmei Jin on 6/7/2022.
//

#include "StaticGrid.h"
#include <queue>

StaticGrid::StaticGrid() {
    num_lng_col_x = ceil((max_lng - min_lng) / step_lng);
    num_lat_row_y = ceil((max_lat - min_lat) / step_lat);

    for (unsigned int i = 0, n = category.getCategoryNum(); i < n; i++) {
        categoryID_vertices.emplace_back(unordered_set<IDTYPE>());
    }
}

void StaticGrid::initialize(const vector<tuple<string, float, float>> &poiList, IndexPoint *pIndex) {
    vertex_array.emplace_back(nullptr);   // the id-0 is invalid

    map<string, int> code2id;
    int currID = 0;
    for(const auto &poi: poiList) {
        // get category-ID
        string categoryCode = get<0>(poi);
        int categoryID;
        auto itr = code2id.find(categoryCode);
        if(itr == code2id.end()) {
            categoryID = category.getCategoryID(categoryCode);
            code2id.insert(make_pair(categoryCode, categoryID));
        }
        else {
            categoryID = itr->second;
        }
        // get POI-ID
        float lng = get<1>(poi), lat = get<2>(poi);
        IDTYPE pid = pIndex->getVertexId(lng, lat);
        if (pid < vertex_array.size()) {
            vertex_array[pid]->addCategoryID(categoryID, categoryCode);    // one POI may have multiple labels
        }
        else {
            // get grid-ID
            long gid = getGridIdByCoor(lng, lat);
            if (gid >= 0) {
                vertex_array.emplace_back(new POI(++currID, new GeoPoint(lng, lat), categoryID, categoryCode));

                // update the corresponding grid content
                auto itr2 = gridIdx_vertices.find(gid);
                if(itr2 != gridIdx_vertices.end()) {
                    itr2->second.insert(currID);
                }
                else {
                    unordered_set<IDTYPE> ps;   ps.insert(currID);
                    gridIdx_vertices.insert(make_pair(gid, ps));
                }
            }
        }
        // update the corresponding category
        categoryID_vertices[categoryID].insert(pid);
    }
}

bool StaticGrid::gridCover(float lng, float lat) const {
    return !(lng < min_lng || lng > max_lng ||
             lat < min_lat || lat > max_lat);
}

bool StaticGrid::isInvalid(long grid_id) const {
    return grid_id < 0 || grid_id >= num_lng_col_x * num_lat_row_y;
}

long StaticGrid::getGridIdByCoor(float lng, float lat) const {
    if (!gridCover(lng, lat))
        return -1;      // out of range

    long idx_row = floor((max_lat - lat) / step_lat);
    long idx_col = floor((lng - min_lng) / step_lng);
    long grid_id = idx_row * num_lng_col_x + idx_col;

    return isInvalid(grid_id) ? -1 : grid_id;
}

// exclude the centerGridID-self
void StaticGrid::getSurrounding(long centerGridID, vector<long> &results) const {
    if (!isInvalid(centerGridID - 1)) {  // left
        results.emplace_back(centerGridID - 1);
    }
    if (!isInvalid(centerGridID + 1)) {  // right
        results.emplace_back(centerGridID + 1);
    }
    if (!isInvalid(centerGridID - num_lng_col_x)) { // up
        results.emplace_back(centerGridID - num_lng_col_x);
    }
    if (!isInvalid(centerGridID - num_lng_col_x - 1)) {
        results.emplace_back(centerGridID - num_lng_col_x - 1);
    }
    if (!isInvalid(centerGridID - num_lng_col_x + 1)) {
        results.emplace_back(centerGridID - num_lng_col_x + 1);
    }
    if (!isInvalid(centerGridID + num_lng_col_x)) { // down
        results.emplace_back(centerGridID + num_lng_col_x);
    }
    if (!isInvalid(centerGridID + num_lng_col_x - 1)) {
        results.emplace_back(centerGridID + num_lng_col_x - 1);
    }
    if (!isInvalid(centerGridID + num_lng_col_x + 1)) {
        results.emplace_back(centerGridID + num_lng_col_x + 1);
    }
}

bool StaticGrid::getSurrounding(long centerGridID, int iteration, vector<long> &results) const {
    if (iteration == 0) {
        results.emplace_back(centerGridID);
    }
    else if (iteration == 1) {
        getSurrounding(centerGridID, results);
    }
    else if (iteration > 1) {
        long row = floor((float) centerGridID / (float) num_lng_col_x);
        long col = centerGridID - row * num_lng_col_x;
        long rMin = MAX(row - iteration, 0), rMax = MIN(row + iteration, num_lat_row_y - 1);
        long cMin = MAX(col - iteration, 0), cMax = MIN(col + iteration, num_lng_col_x - 1);

        // up and down
        for (long c = cMin; c <= cMax; c++) {
            long grid_id = rMin * num_lng_col_x + c;    // down
            if (!isInvalid(grid_id)) {
                results.emplace_back(grid_id);
            }
            grid_id = rMax * num_lng_col_x + c;         // up
            if (!isInvalid(grid_id)) {
                results.emplace_back(grid_id);
            }
        }
        // left and right
        for (long r = rMin + 1; r < rMax; r++) {
            long grid_id = r * num_lng_col_x + cMin;    // left
            if (!isInvalid(grid_id)) {
                results.emplace_back(grid_id);
            }
            grid_id = r * num_lng_col_x + cMax;         // right
            if (!isInvalid(grid_id)) {
                results.emplace_back(grid_id);
            }
        }
        return rMin == 0 || rMax == num_lat_row_y - 1 || cMax == 0 || cMax == num_lng_col_x - 1;  // the boundary
    }
    return false;
}

bool compareByDist(pair<IDTYPE, double> &p1, pair<IDTYPE, double> &p2) {
    return p1.second < p2.second;
}

POI* StaticGrid::getNearestPOI(float lng, float lat) const {
    IDTYPE vertexID = 0;    // an invalid ID
    int iteration = 0;
    long centerGridID = getGridIdByCoor(lng, lat);
    vector<long> gridIDs;
    while (vertexID == 0 && iteration < MAX_ITERATION) {
        gridIDs.clear();
        bool boundary = getSurrounding(centerGridID, iteration++, gridIDs);

        // find the nearest POI based on the grid index
        double minDist = INFINITY;
        for (long gid: gridIDs) {
            auto pointset = gridIdx_vertices.find(gid);
            if (pointset != gridIdx_vertices.end()) {
                for (IDTYPE poiID: pointset->second) {
                    double distance = vertex_array[poiID]->point->computeDistance(lng, lat);   // raw point and POI
                    if (minDist > distance) {
                        minDist = distance;
                        vertexID = poiID;
                    }
                }
            }
        }
        if(boundary)
            break;
    }
    return vertex_array[vertexID];
}

vector<pair<POI *, double>> StaticGrid::getSurroundingPOIs(float lng, float lat, float maxDist, unsigned int maxNum) const {
    if(maxNum == 0) {
        maxNum = category.getCategoryNum();
    }
    int iteration = 0;
    long centerGridID = getGridIdByCoor(lng, lat);
    vector<long> gridIDs;
    vector<pair<IDTYPE, double>> poiQueue;

    while (iteration < MAX_ITERATION) {
        gridIDs.clear();
        bool boundary = getSurrounding(centerGridID, iteration++, gridIDs);

        // find the nearest POI based on the grid index
        bool existValid = false;
        for (long gid: gridIDs) {
            if(poiQueue.size() > maxNum) {
                Rectangle rec = createBoundaryRec(gid);
                double minDist = rec.computeMinDist(lng, lat);
                if(minDist > poiQueue.back().second || minDist > maxDist) {
                    continue;
                }
            }

            existValid = true;
            auto pointset = gridIdx_vertices.find(gid);
            if (pointset != gridIdx_vertices.end()) {
                for (IDTYPE poiID: pointset->second) {
                    POI *poi = vertex_array[poiID];
                    double distance = poi->point->computeDistance(lng, lat);   // raw point and POI
                    if(distance <= maxDist || poiQueue.empty()) {
                        poiQueue.emplace_back(make_pair(poiID, distance));
//                        if(distance == 0) {
//                            break;  // exact the same coordinates, no further searching
//                        }
                        if(poiQueue.size() >= maxNum) {
                            sort(poiQueue.begin(), poiQueue.end(), compareByDist);  // keep sorting
                        }
                        while (poiQueue.size() > maxNum) {
                            poiQueue.pop_back();
                        }
                    }
                }
            }
        }
        if(boundary || !existValid)
            break;
    }

    vector<pair<POI *, double>> rtn;
    for(const pair<IDTYPE, double> &pair: poiQueue){
        rtn.emplace_back(make_pair(vertex_array[pair.first], pair.second));
        if(rtn.size() == maxNum) {
            break;
        }
    }
    return rtn;
}

Rectangle StaticGrid::createBoundaryRec(long gid) const {
    long row = floor((float) gid / (float) num_lng_col_x);
    long col = gid - row * num_lng_col_x;

    float maxLat = this->max_lat - (float) row * step_lat;
    float minLng = this->min_lng + (float) col * step_lng;
    return {minLng, maxLat - step_lat, minLng + step_lng, maxLat};
}

/**
 * Note that the circular region is approximately represented by a rectangle (larger than itself)
 * IF @param approximate is allowed, then no exact computation for the circle but only its MBR
 */

void StaticGrid::getOverlapGrids(GeoPoint *p, float radius, vector<long> &overlapGrids, bool approximate) const{
    int iteration = 0;
    long centerGridID = getGridIdByCoor(p->getLongitude(), p->getLatitude());
    Rectangle currBoundary = createBoundaryRec(centerGridID);    // will be enlarged with iteration
    string bStr = currBoundary.toString();

    vector<long> gridIDs;
    Rectangle approxCircle(p, radius);  // a rectangle to represent the circular region
    string cStr = approxCircle.toString();
    string pStr = p->toString();

    while (true) {
        gridIDs.clear();
        getSurrounding(centerGridID, iteration++, gridIDs);

        // find the nearest POI based on the grid index
        for (long gid: gridIDs) {
            Rectangle currRec = createBoundaryRec(gid);
            string curStr = currRec.toString();

            if(currRec.intersect(approxCircle)) {   // the grid cell is qualified
                if(approximate) {
                    overlapGrids.emplace_back(gid);
                }
                else {
                    double dist = currRec.computeMinDist(p);
                    if (dist <= radius) {
                        overlapGrids.emplace_back(gid);
                    }
                }
            }
            currBoundary.enlarge(&currRec);
        }

        if(currBoundary.fullContain(approxCircle)) {
            break;            // no need to check next round
        }
    }
}

void StaticGrid::getOverlapPOIs(GeoPoint *queryP, float radius, vector<POI *> &poiList) {
    vector<long> overlapGrids;
    getOverlapGrids(queryP, radius, overlapGrids);

    for(long gid: overlapGrids) {
        auto pointset = gridIdx_vertices.find(gid);
        if(pointset != gridIdx_vertices.end()) {
            for(IDTYPE pid: pointset->second) {
                POI *poi = vertex_array[pid];
                if(poi != nullptr && poi->point != nullptr) {
                    double dist = queryP->computeDistance(poi->point);
                    if(dist <= radius) {
                        poiList.emplace_back(poi);
                    }
                }
            }
        }
    }
}

vector<POI *> StaticGrid::getAllPOIs() {
    return vertex_array;
}

vector<unordered_set<IDTYPE>> StaticGrid::getCategoryPOIs(){
    return categoryID_vertices;
}
