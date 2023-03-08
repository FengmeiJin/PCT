//
// Created by Fengmei JIN on 2022/9/20.
//
#include "Perturbation.h"

void Perturbation::applyPerturbation(const vector<Trip *> &originalTrajs, vector<Trip *> &anonymizedTrajs,
                                     const float *epsilon, IndexPoint *pIndex, StaticGrid *poiGrid,
                                     bool EnableGeo, bool EnableSemantic, bool EnableGS, bool fixedSeed) {
    unsigned int dataTotal = originalTrajs.size(), cnt = 0;

    printf("[PROGRESS] ");
    fflush(stdout);
    for (Trip *rawTrip: originalTrajs) {
        string id = rawTrip->getId();
        if(EnableGeo) {
            Trip *geoTrip = applyGeoIndis(epsilon[0], pIndex, rawTrip, fixedSeed);
            anonymizedTrajs.emplace_back(geoTrip);
        }
        else if(EnableGS) {
            // don't align for geoTrip, since two-phase perturbation
            Trip *gsTrip = applyGeoSemanticIndis(epsilon, pIndex, poiGrid, rawTrip, fixedSeed);
            anonymizedTrajs.emplace_back(gsTrip);
        }
        else if (EnableSemantic){
            Trip *sTrip = applySemanticIndis(epsilon[1], pIndex, poiGrid, rawTrip, fixedSeed);
            anonymizedTrajs.emplace_back(sTrip);
        }

        // print out progress
        if (dataTotal < 10) {
            printf("%d...", ++cnt);
            fflush(stdout);
        }
        else if (++cnt % (dataTotal / 10) == 0) {
            printf("%d%%...", (int) (cnt * 100.0 / dataTotal));
            fflush(stdout);
        }
    }
    printf(" done!\n");
}

// ====================     Geo Indistinguishability     ====================

tuple<float, float, float> Perturbation::findGeoIndisRegion(GeoPoint *p, float epsilon, bool fixedSeed) {
    // sample noise on polar coordinate system
    pair<float, float> noises = NoiseGenerator::fromPlanarLaplaceDistribution(epsilon, fixedSeed);
    float angle = noises.first;
    float radius = noises.second;

    // ===== DEBUG
    if(angle < 0 || radius <= 0) {  // should not happen
        printf("[ALERT] the sampled angle = %.3f, radius = %.3f\n", angle, radius);
    }

    // transfer to the normal coordinates in geographic space
    float lngX = p->getLongitude();
    float latY = p->getLatitude();
    pair<double, double> rtn = Rectangle::getOffsetDegree(cos(angle) * radius, sin(angle) * radius, latY);
    auto degreeLng = (float) rtn.first;
    auto degreeLat = (float) rtn.second;
    if(angle <= 180) {
        latY += degreeLat;
        if (angle <= 90)
            lngX += degreeLng;
        else
            lngX -= degreeLng;
    }
    else {
        latY -= degreeLat;
        if (angle <= 270)
            lngX -= degreeLng;
        else
            lngX += degreeLng;
    }
    return make_tuple(lngX, latY, radius);
}

Trip *Perturbation::applyGeoIndis(float epsilon, IndexPoint *pIndex, Trip *rawTrip, bool fixedSeed)  {
    Point *ptr = rawTrip->getFirstPoint();
    vector<IndisRegion> noisyRegions;
    while (ptr != nullptr && !ptr->isTailPtr()) {
        GeoPoint *p = pIndex->getPointById(ptr->getPointId());
        tuple<float, float, float> noisyRegion = findGeoIndisRegion(p, epsilon, fixedSeed);
        float radius = get<2>(noisyRegion);
        noisyRegions.emplace_back(IndisRegion(ptr, get<0>(noisyRegion), get<1>(noisyRegion), radius));
        ptr = ptr->getNextPointer();
    }

    // we don't do the alignment (from a noisy point to its nearest POI)
    // since it doesn't make sense for a circular region
    // no alignment
    return new Trip(rawTrip->getId(), noisyRegions, pIndex);
}

// ==================== GeoSemantic Indistinguishability ====================
/*
 * implement "semantic distance" by Generalized Jaccard distance = 1 - Jaccard similarity
 * */
double Perturbation::computeJaccardSim(unsigned int cNum, const float *dist1, const float *dist2) {
    if(dist1 == nullptr || dist2 == nullptr)
        return 0;
    double numerator = 0, denominator = 0;
    for(int i = 0; i < cNum; i++) {
        numerator += MIN(dist1[i], dist2[i]);
        denominator += MAX(dist1[i], dist2[i]);
    }
    return numerator / denominator;
}

double Perturbation::computePoiProbability(float epsilon, const set<pair<string, int>>& categoryIDs, unsigned int cNum, const float *distribution) {
    if(categoryIDs.empty()) {
        return 0;
    }
    // distribution for the current candiP
    float *pDistribution = new float [cNum];
    memset(pDistribution, 0, cNum);
    for(const auto &cpair: categoryIDs) {
        pDistribution[cpair.second]++;
    }
    // normalization by the sum
    if(categoryIDs.size() > 1) {
        float sum = (float) categoryIDs.size();
        for (int i = 0; i < cNum; i++) {
            if (pDistribution[i] > 0)
                pDistribution[i] /= sum;
        }
    }
    double jSim = computeJaccardSim(cNum, pDistribution, distribution);      // range [0, 1]
    return exp(epsilon * (1 - jSim) / 2.0);
}

bool compareByProb(pair<double, int> &p1, pair<double, int> &p2) {
    return p1.first < p2.first;
}

/**
 * @param p: the center point of the specific region
 * @param radius: the radius of the circle
 * @param oriDistribution: the oriDistribution of the original true location
 * */
POI* Perturbation::samplePointFromRegion(StaticGrid *poiGrid, GeoPoint *p, float radius, const float *oriDistribution,
                                         float epsilon, bool fixedSeed){
    vector<POI *> overlapPOIs;
    float r = radius, amplifyRatio = 0.5f;
    int amplifyCnt = 0;
    while (overlapPOIs.empty()) {
        poiGrid->getOverlapPOIs(p, r, overlapPOIs);
        r = radius * (1 + amplifyRatio * (float) (++amplifyCnt));  // may result in duplicate computation
    }

    unsigned int i = 0;
    double sum = 0;
    vector<pair<double, int>> prob_array;
    for(POI *poi: overlapPOIs) {
        if(poi != nullptr) {
            // the real prob is proportional to this value
            double prob = computePoiProbability(epsilon, poi->categoryIDs, poiGrid->getCategoryNum(), oriDistribution);
            if(!isnan(prob) && prob > 0) {
                prob_array.emplace_back(make_pair(prob, i++));
                sum += prob;
            }
        }
    }

    unsigned int totalCandi = prob_array.size();
    for(i = 0; i < totalCandi; i++) {
        prob_array[i].first /= sum;   // normalization probability
    }
    sort(prob_array.begin(), prob_array.end(), compareByProb);

    double randomV = NoiseGenerator::doubleFromUniformDistribution(0, 1, fixedSeed);
    sum = 0;
    unsigned int selected = totalCandi;
    for(i = 0; i < totalCandi; i++) {
        if(randomV >= sum && randomV < prob_array[i].first + sum) {
            selected = i;
            break;
        }
        sum += prob_array[i].first;
    }
    if(selected < totalCandi) {
        int idx = prob_array[selected].second;
        return overlapPOIs[idx];
    }
    printf("[ERROR] Fail to sample a POI from the specific region centered by (%s) with the radius of %.3f!\n", p->toString().c_str(), radius);
    return nullptr;
}

// GeoSemantic: the candidate range is not the entire space, but a specified circular region
Trip *Perturbation::applyGeoSemanticIndis(const float *epsilon, IndexPoint *pIndex, StaticGrid *poiGrid,
                                          Trip *rawTrip, bool fixedSeed) {
    string uid = rawTrip->getId();
    Trip *geoTrip = applyGeoIndis(epsilon[0], pIndex, rawTrip, fixedSeed);

    Point *ptr = geoTrip->getFirstPoint();
    vector<IndisRegion> noisyRegions;
    while (ptr != nullptr && !ptr->isTailPtr()) {
        // a center and a radius define a region
        float curRadius = ptr->getRadius();
        GeoPoint *p = pIndex->getPointById(ptr->getPointId());      // noisy point perturbed by geo-indis
        GeoPoint *oriP = pIndex->getPointById(ptr->getPoint_ori()->getPointId());   // raw point

        POI *selectedPOI = samplePointFromRegion(poiGrid, p, curRadius, ptr->getPoiDistribution(), epsilon[1], fixedSeed);
        if(selectedPOI != nullptr) {
            GeoPoint *newP = selectedPOI->point;
            float newRadius = MAX(p->computeDistance(newP), oriP->computeDistance(newP));
            noisyRegions.emplace_back(IndisRegion(ptr->getPoint_ori(), newP->getLongitude(), newP->getLatitude(), newRadius));
        }
        else {
            noisyRegions.emplace_back(IndisRegion(ptr->getPoint_ori(), p->getLongitude(), p->getLatitude(), curRadius));   // keep as the previous noisy point
        }
        ptr = ptr->getNextPointer();
    }

    return new Trip(rawTrip->getId(), noisyRegions, pIndex);
}

// ==================     Semantic Indistinguishability     =================

const int ENTIRE_POIS = -1;

/**
 * @param p: the center point of the specific region
 * @param oriDistribution: the oriDistribution of the original true location
 * @param epsilon: the lower the epsilon, the stronger the privacy is
 * */
POI* Perturbation::samplePointFromEntireSpace(StaticGrid *poiGrid, const float *oriDistribution,
                                              float epsilon, bool fixedSeed){

    vector<POI *> allPOIs = poiGrid->getAllPOIs();  // the first is null as 0 is invalid ID
    unsigned int cNum = poiGrid->getCategoryNum();

    if(ENTIRE_POIS == 0) {
        vector<double> prob_sum_array;  // [k] = sum of 0 ~ k
        prob_sum_array.emplace_back(0);

        double estimatedProb = exp(epsilon / 2);  // upper bound
//    double estimatedProb = 1;   // lower bound
        double normalization = estimatedProb * (float) allPOIs.size();
        double sum = 0;
        double randomV = NoiseGenerator::doubleFromUniformDistribution(0, 1, fixedSeed);
        for(POI *poi: allPOIs) {
            if(poi != nullptr) {
                double prob = computePoiProbability(epsilon, poi->categoryIDs, cNum, oriDistribution);
                normalization = normalization - estimatedProb + prob;
                prob /= normalization;  // a coarse normalization

                if(randomV >= sum && randomV <= prob + sum) {
                    return poi;
                }
                sum += prob;
                prob_sum_array.emplace_back(sum);
            }
        }

        int low = 1, high = (int) prob_sum_array.size() - 1;
        while (low <= high) {
            int mid = (low + high) / 2;
            double v = prob_sum_array[mid];
            double prev = prob_sum_array[mid - 1];
            if(randomV == v || (randomV < v && randomV >= prev)) {
                return allPOIs[mid];
            }
            else if (randomV < v) {
                high = mid - 1;
            }
            else {
                low = mid + 1;
            }
        }
    }
    else {
        unordered_set<IDTYPE> checkedPOIs;
        vector<unordered_set<IDTYPE>> category2pois = poiGrid->getCategoryPOIs();
        vector<pair<double, int>> oriDistribution_sorted;
        for (int i = 0; i < cNum; i++) {
            oriDistribution_sorted.emplace_back(make_pair(oriDistribution[i], i));
        }
        sort(oriDistribution_sorted.begin(), oriDistribution_sorted.end(), compareByProb);
//        for (unsigned int i = cNum - 1; i >= 0; i--) {
        for (unsigned int i = 0; i < cNum; i++) {
            int idx = oriDistribution_sorted[i].second;
            unordered_set<IDTYPE> candidates = category2pois[idx];
            for (IDTYPE pid: candidates) {
                if (checkedPOIs.find(pid) == checkedPOIs.end()) {
                    POI *poi = allPOIs.at(pid);
                    if(poi != nullptr) {
                        // the real prob is proportional to this value
                        double prob = computePoiProbability(epsilon, poi->categoryIDs, cNum, oriDistribution);
                        if(!isnan(prob) && prob > 0) {
                            double normalized_prob = (prob - 1) / (exp(epsilon / 2) - 1);
                            double randomP = NoiseGenerator::doubleFromUniformDistribution(0, 1, fixedSeed);
                            if (randomP >= normalized_prob) {
                                return allPOIs[pid];
                            }
                        }
                    }
                    checkedPOIs.insert(pid);
                }
            }
        }
    }

    int totalCandi = (int) allPOIs.size();
    int randomInt = NoiseGenerator::intFromUniformDistribution(1, totalCandi);
    return allPOIs[randomInt];
}

Trip *Perturbation::applySemanticIndis(float epsilon, IndexPoint *pIndex, StaticGrid *poiGrid,
                                       Trip *rawTrip, bool fixedSeed) {
    Point *ptr = rawTrip->getFirstPoint();
    vector<IndisRegion> noisyRegions;
    while (ptr != nullptr && !ptr->isTailPtr()) {
        POI *selectedPOI = samplePointFromEntireSpace(poiGrid, ptr->getPoiDistribution(), epsilon, fixedSeed);
        if(selectedPOI != nullptr) {
            GeoPoint *p = pIndex->getPointById(ptr->getPointId());
            GeoPoint *newP = selectedPOI->point;
            float newRadius = (float) p->computeDistance(newP);
            noisyRegions.emplace_back(IndisRegion(ptr, newP->getLongitude(), newP->getLatitude(), newRadius));
        }
        else {
            printf("ERROR in Perturbation::applySemanticIndis: return null POI!\n");    // should not happen!
        }
        ptr = ptr->getNextPointer();
    }

    return new Trip(rawTrip->getId(), noisyRegions, pIndex);
}