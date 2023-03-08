#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <algorithm>
#include <cmath>

#define MY_PI 180
#define MY_PI_RADIAN 3.1415926

// to simulate integral for contact probability
#define SPLIT_NUM pow(10, 3)

// r is the radius of a circle
static double independent_prob_circle(double r, double R){
    return exp(-1 * r) / (1 - exp(-1 * R));
}

static double correlated_prob_rho(double rho) {
    double numerator = exp(-1 * abs(rho));  // lambda is regarded as fixed when integral for rho (outside)
    double denominator = 2 * (1 - exp(-1 * MY_PI));
    return numerator / denominator;
}

static double angle_by_arccos(double v1, double v2, double opposite) {
    double radian = acos((v1 * v1 + v2 * v2 - opposite * opposite) / (2 * v1 * v2));
    return radian * MY_PI / MY_PI_RADIAN;
}

static double angle_by_r(double r, double dist_pq, double delta_s) {
    double theta_angle = angle_by_arccos(r, dist_pq, delta_s);  // "half theta" actually
    double percent = theta_angle / MY_PI;
    if(isnan(percent)) {
        percent = 1.0 / SPLIT_NUM;
    }
    return percent;
}

static double angle_within_PI(double rho) {
    return rho <= MY_PI ? rho : (-1 * (rho - MY_PI));
}

//template<typename FUNC_T, typename RTN_T=double, typename NUMBER_T=unsigned long long>
static double independent_integral(double left, double right, double R,
                            bool partialCircle = false, double pDist = 0, double delta_s = 0) {
    if (left > right) std::swap(left, right);
    double len = (right - left) / SPLIT_NUM;
    double res = 0;
    for (double r = left; r < right; r += len) {
        double v = independent_prob_circle(r, R);
        if(partialCircle) {
            v = v * angle_by_r(r, pDist, delta_s);
        }
        res += v * len;
    }
    return res;
}

/**
 * @param R the radius of a noisy region G_x,z
 * @param dist_pq the distance between two points
 * @param delta_s max distance threshold to define a contact
 * */
static double computeIndependentContactProb(double dist_pq, double R, double delta_s) {
    double prob = 0;

    if(dist_pq >= delta_s + R){   // disjoint
        prob = 0;
    }
    else {
        if(delta_s >= R && dist_pq <= delta_s - R) {  // C(z,R) is fully within C(q,delta_s)
            prob = 1;
        }
        else {
            if (dist_pq >= delta_s) {     // z is out of C(q, delta_s)
                double lb = dist_pq - delta_s;
                double rb = MIN(R, dist_pq + delta_s);
                prob = independent_integral(lb, rb, R, true, dist_pq, delta_s);
            }
            else {  // dist_pq < delta_s,     z is within C(q, delta_s)
                double lb = 0;
                double rb = delta_s - dist_pq;
                prob = independent_integral(lb, rb, R);     // integral for the entire circle within a range

                lb = rb;
                rb = MIN(R, dist_pq + delta_s);
                prob += independent_integral(lb, rb, R, true, dist_pq, delta_s);
            }
        }
    }
    return MIN(prob, 1);    // imperfect due to the split granularity
}

static double correlated_integral(double left, double right, double dist_pq, double dist_p2prev, double dist_q2prev, double R, double delta_s){
    if (left > right) std::swap(left, right);
    double len_r = (right - left) / SPLIT_NUM;
    double res = 0;
    double qzz_angle = angle_by_arccos(dist_q2prev, dist_pq, dist_p2prev);
    for (double r = left; r < right; r += len_r) {
        double theta_half = angle_by_arccos(r, dist_pq, delta_s);
        double rho_left = angle_within_PI(qzz_angle - theta_half);
        double rho_right = angle_within_PI(qzz_angle + theta_half);
        if (rho_left > rho_right) std::swap(rho_left, rho_right);

        // integral for rho first
        double len_rho = (rho_right - rho_left) / SPLIT_NUM;
        double v = 0;
        for(double rho = rho_left; rho < rho_right; rho += len_rho) {
            v += correlated_prob_rho(rho) * len_rho;
        }
        v = MIN(v, 1) * independent_prob_circle(r, R);   // imperfect due to the split granularity
        res += v * len_r;
    }
    return res;
}

/**
 * @param dist_pq the distance between a query point "p" and a candidate point "q"
 * @param dist_p2prev the distance between "p" and "q"'s previous point
 * @param dist_q2prev the distance between "q" and "q"'s previous point
 * @param R the radius of a noisy region G_x,z
 * @param delta_s max distance threshold to define a contact
 * */
static double computeCorrelatedContactProb(double dist_pq, double dist_p2prev, double dist_q2prev, double R, double delta_s) {
    double prob = 0;

    if(dist_pq >= delta_s + R){   // disjoint
        prob = 0;
    }
    else {
        if(delta_s >= R && dist_pq <= delta_s - R) {  // C(z,R) is fully within C(q,delta_s)
            prob = 1;
        }
        else {
            if (dist_pq >= delta_s) {     // z is out of C(q, delta_s)
                double lb = dist_pq - delta_s;
                double rb = MIN(R, dist_pq + delta_s);
                prob = correlated_integral(lb, rb, dist_pq, dist_p2prev, dist_q2prev, R, delta_s);
            }
            else {  // dist_pq < delta_s,     z is within C(q, delta_s)
                double lb = 0;
                double rb = delta_s - dist_pq;
                prob = independent_integral(lb, rb, R);     // here is same as independent version

                lb = rb;
                rb = MIN(R, dist_pq + delta_s);
                prob += correlated_integral(lb, rb, dist_pq, dist_p2prev, dist_q2prev, R, delta_s);
            }
        }
    }
    return MIN(prob, 1);    // imperfect due to the split granularity
}


/** @param qtr: belong to susceptible
 *  @param ptr: belong to spreader
 *  */
static double computeContactProb(Point *qtr, Point *ptr, IndexPoint *pIndex, float distMax, bool independentCP){
    // compute contact probability
    double prob = 0;
    double R = qtr->getRadius();
    if (R > 0) {
        double dist_pq = pIndex->computeDistance(ptr->getPointId(), qtr->getPointId());
        double delta_s = ptr->getRadius() + distMax;       // for a point, it is the threshold;

        // for a region, it sums up its radius
        Point *qPrev = qtr->getPrevPointer();
        if(independentCP || qPrev->isHeadPtr()) {
            prob = computeIndependentContactProb(dist_pq, R, delta_s);
        }
        else {
            double dist_p2prev = pIndex->computeDistance(ptr->getPointId(), qPrev->getPointId());
            double dist_q2prev = pIndex->computeDistance(qtr->getPointId(), qPrev->getPointId());
            prob = computeCorrelatedContactProb(dist_pq, dist_p2prev, dist_q2prev, R, delta_s);
        }
    }
    return prob;
}

#endif