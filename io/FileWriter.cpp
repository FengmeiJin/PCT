//
// Created by Fengmei JIN on 2022/7/18.
//

#include "FileWriter.h"

void FileWriter::outputContacts(const string &outputFile, int multiHop, bool enableCP,
                                const vector<pair<string, ContactEvent *>> &susceptible2event) {
    ofstream fout(outputFile, ios_base::out);
    if(fout.is_open()) {
        stringstream ss;
//        printf("\n[INFO] Printing out the newly found contacts with the contact start time ...\n");
        for(const auto &pair: susceptible2event) {
            ContactEvent *event = pair.second;
            ss << pair.first << ", " << event->spreaderID << ", "     // the corresponding spreaders may vary as long as the cst is same
               << event->contactStartTime;
            if(multiHop > 1) {
                ss << "," << event->kHop;
            }
            if(enableCP) {
                ss << "," << event->probability;
            }
            ss << "\n";
        }
        fout << ss.str();
        fout.close();
        ss.clear();
    }
    else {
        printf("[ERROR] cannot open the file %s for output!\n", outputFile.c_str());
    }
}

// return: ori-distribution, ori-label
string distributionStr(const float *distribution, int dim = 17) {
    string res;
    if(distribution == nullptr) {
        for(int i = 0; i < dim; i++) {
            res.append("0,");
        }
        res.append("0");
        return res;
    }
    char *buffer = new char [10];
    int highestIdx = 0;
    float highestV = distribution[0];
    for(int i = 0; i < dim; i++) {
        float v = distribution[i];
        if(v > highestV) {
            highestV = v;
            highestIdx = i;
        }
        if(v > 0){
            snprintf(buffer, 10, "%.5f,", v);
            res.append(buffer);
            buffer[0] = '\0';
        }
        else {
            res.append("0,");
        }
    }
    res.append(to_string(highestIdx));
    return res;
}

void FileWriter::outputNoisyTrips(const string &outputFile, IndexPoint *pIndex, const vector<Trip *> &noisyTrips) {
    ofstream fout(outputFile, ios_base::out);
    if(!fout.is_open()) {
        printf("[ERROR] cannot open the file %s for output!\n", outputFile.c_str());
    }
    else {
        printf("[INFO] Output the noisy trips to the file: %s...\n", outputFile.c_str());
        fout.close();

        unsigned int cnt = 0, total = noisyTrips.size();
        printf("[PROGRESS] ");
        fflush(stdout);

        int LEN = 1000;
        char *buffer = new char[LEN];
        for(const Trip *traj: noisyTrips) {
            fout.open(outputFile, ios_base::app);

            stringstream ss;
            string id = traj->getId();
            Point *ptr = traj->getFirstPoint();
            while (ptr != nullptr && !ptr->isTailPtr()) {
                GeoPoint *p = pIndex->getPointById(ptr->getPointId());
                Point *ptr_ori = ptr->getPoint_ori();
                GeoPoint *p_ori = nullptr;
                if(ptr_ori != nullptr) {
                    p_ori = pIndex->getPointById(ptr_ori->getPointId());
                }
                else {
                    printf("ERROR in FileWriter::outputNoisyTrips - a null ori ptr for a noisy point!\n");
                    return;
                }
                // one line one point
                // format: id, ori-lng, ori-lat, ori-distribution, ori-label, noisy-lng, noisy-lat, radius, noisy-distribution, noisy-label, arriveT, leaveT
                double dist = p->computeDistance(p_ori);
                snprintf(buffer, LEN, "%s,%.7f,%.7f,%s,%.7f,%.7f,%.7f,%s,%lu,%lu\n", id.c_str(),
                        p_ori->getLongitude(), p_ori->getLatitude(), distributionStr(ptr_ori->getPoiDistribution()).c_str(),
                        p->getLongitude(), p->getLatitude(), dist, distributionStr(ptr->getPoiDistribution()).c_str(),
                        ptr_ori->getArriveTimestamp(), ptr_ori->getLeaveTimestamp());
                ss << buffer;
                buffer[0] = '\0';

                ptr = ptr->getNextPointer();
            }
            fout << ss.str();
            fout.close();
            ss.clear();

            // print out progress
            if (total >= 10 && ++cnt % (total / 10) == 0) {
                printf("%d%%...", (int) (cnt * 100.0 / total));
                fflush(stdout);
            }
        }

        printf("done!\n");
        fflush(stdout);
    }
}