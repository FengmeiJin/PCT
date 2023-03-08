//
// Created by Fengmei Jin on 6/7/2022.
//

#ifndef PCT_STATICGRID_H
#define PCT_STATICGRID_H

#include "../commonheader.h"
#include "IndexPoint.h"
#include "Category.h"

const float step_default = 0.01;
const int MAX_ITERATION = 5;

class StaticGrid {

    // spatial range
    float min_lng{BJ_LNG_MIN};
    float min_lat{BJ_LAT_MIN};
    float max_lng{BJ_LNG_MAX};
    float max_lat{BJ_LAT_MAX};
    float step_lng{step_default};
    float step_lat{step_default};
    IDTYPE num_lng_col_x;
    IDTYPE num_lat_row_y;

    // POI information
    vector<POI *> vertex_array;

    unordered_map<long, unordered_set<IDTYPE>> gridIdx_vertices;   // gridID -> inside POIs

    vector<unordered_set<IDTYPE>> categoryID_vertices;

    Category category;

public:

    StaticGrid();

    void initialize(const vector<tuple<string, float, float>> &poiList, IndexPoint *pIndex);

    // ====== functions based on POI information

    inline unsigned int getValidPoiNum() const {
        return vertex_array.size();
    }
    inline unsigned int getCategoryNum() {
        return category.getCategoryNum();
    }
    inline vector<string> getSubCategoryCodes(int cid) {
        return category.getSubCategoryCodes(cid);
    }

    void computeAggregatedDistribution(const vector<pair<POI *, double>> &poiList, float *rtn, bool involveDist = false) {
        category.computeAggregatedDistribution(poiList, rtn, involveDist);
    }

    // POI with the distance larger than the threshold (e.g., 2km) cannot contribute to the POI distribution of a specific point
    vector<pair<POI *, double>> getSurroundingPOIs(float lng, float lat, float maxDist = 2, unsigned int maxNum = 0) const;

    POI* getNearestPOI(float lng, float lat) const;

    POI* getPOI(IDTYPE pid) const {
        return pid > 0 && pid < vertex_array.size() ? vertex_array[pid] : nullptr;
    }

    void getOverlapPOIs(GeoPoint *queryP, float radius, vector<POI *> &poiList);

    // ====== functions purely based on spatial information

    inline unsigned int getValidGridNum() const {
        return gridIdx_vertices.size();
    }

    long getGridIdByCoor(float lng, float lat) const;

    bool gridCover(float lng, float lat) const;

    bool isInvalid(long grid_id) const;

    void getSurrounding(long centerGridID, vector<long> &results) const;

    bool getSurrounding(long centerGridID, int iteration, vector<long> &results) const;

    void getOverlapGrids(GeoPoint *p, float radius, vector<long> &overlapGrids, bool approximate = true) const;

    Rectangle createBoundaryRec(long gid) const;

    vector<POI *> getAllPOIs();

    vector<unordered_set<IDTYPE>> getCategoryPOIs();
};


#endif //PCT_STATICGRID_H
