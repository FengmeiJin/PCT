//
// Created by Fengmei JIN on 2022/8/10.
//

#ifndef PCT_PERTURBATION_H
#define PCT_PERTURBATION_H

#include "NoiseGenerator.h"

class Perturbation {
private:
    static Trip * applySemanticIndis(float epsilon, IndexPoint *pIndex, StaticGrid *poiGrid, Trip *rawTrip, bool fixedSeed);
    static POI* samplePointFromEntireSpace(StaticGrid *poiGrid, const float *oriDistribution,
                                           float epsilon, bool fixedSeed);

    static Trip *applyGeoSemanticIndis(const float *epsilon, IndexPoint *pIndex, StaticGrid *poiGrid, Trip *rawTrip, bool fixedSeed);
    static double computeJaccardSim(unsigned int cNum, const float *dist1, const float *dist2);
    static double computePoiProbability(float epsilon, const set<pair<string, int>> &categoryIDs, unsigned int cNum, const float *distribution);

    /**
    * @param p: the center point to specify a region, along with the @param radius of the circle
    * @param oriDistribution: the oriDistribution of the original true location
    * */
    static POI *samplePointFromRegion(StaticGrid *poiGrid, GeoPoint *p, float radius, const float *oriDistribution, float epsilon, bool fixedSeed);

    static Trip *applyGeoIndis(float epsilon, IndexPoint *pIndex, Trip *rawTrip, bool fixedSeed);
    static tuple<float, float, float> findGeoIndisRegion(GeoPoint *p, float epsilon, bool fixedSeed);

public:
    static void applyPerturbation(const vector<Trip *> &originalTrajs, vector<Trip *> &anonymizedTrajs,
                                  const float *epsilon, IndexPoint *pIndex, StaticGrid *poiGrid,
                                  bool EnableGeo, bool EnableSemantic, bool EnableGS, bool fixedSeed);
};


#endif //PCT_PERTURBATION_H
