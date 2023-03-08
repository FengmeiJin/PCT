//
// Created by Fengmei Jin on 5/7/2022.
//

#ifndef PCT_NOISEGENERATOR_H
#define PCT_NOISEGENERATOR_H

#include <random>
#include <set>
#include "../spatial/Trip.h"

class NoiseGenerator {

public:

    static int intFromUniformDistribution(int a, int b, bool fixedSeed = false);

    static double doubleFromUniformDistribution(double a, double b, bool fixedSeed = false);

    static pair<float, float> fromPlanarLaplaceDistribution(float epsilon, bool fixedSeed = false);

    static vector<Trip *> sampleQuery(int expectedNum, const vector<Trip *> &data, bool fixedSeed = false);
};

#endif //PCT_NOISEGENERATOR_H
