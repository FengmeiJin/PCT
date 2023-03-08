//
// Created by Fengmei JIN on 2023/2/2.
//

#ifndef PCT_INDEX_H
#define PCT_INDEX_H

#include "SearchByIndex.h"
#include "functions.h"

/**
 * For each susceptible @param susceptibleTraj check if he contacted with any spreader in @param spreaderTrips
 * @param spreaderTrips: the spreader trajectory along with his valid contact start time
 *                      !!! sorted by the start time of their spreading !!!
 * @param impossibleTrajs: those which are impossible to contact with @param susceptibleTraj (due to a coarse preFiltering)
 * @param earlyTerminate: 1) if already find a contacted spreader (holding an earliest contact start time so far)
 *                              then terminate searching another spreader if he doesn't have contact event before this time
 *                        2) if find any spreader contacting with the @param susceptibleTraj then terminate (only when @param directContact is true)
 *                              rather than computing all in the @param spreaderTrips
 *
 * @return the spreader with the earliest contact start time to the susceptible query
 * */
ContactEvent* searchDirectContact_linear(const Trip *susceptibleTraj, int qHop, const vector<pair<Trip *, ContactEvent *>> &spreaderTrips,
                                         bool ground_truth, IndexPoint *pIndex, float distThreshold,
                                         TIMETYPE contactMinDuration, float contactMinProb, bool independentCP,
                                         vector<TIMETYPE> &timeCosts) {
    const string queryID = susceptibleTraj->getId();
    const bool enableCP = (!ground_truth && contactMinProb > 0);

    // linearly scan all candidates
    int checkedQueryPoint = 0, checkedCandiNum = 0;
    TIMETYPE earliestCST = 0;               // 0 indicates hasn't found a contacted spreader yet
    string earliestSpreaderID;              // expect to return the one with the earliest contact start time

    double aggregatedNonContactProb = 1;
    for (const auto &spreaderPair: spreaderTrips) {
        Trip *spreaderTraj = spreaderPair.first;
        string spreaderID = spreaderTraj->getId();
        TIMETYPE validSpreadStartTime = spreaderPair.second->getStartTime();

        if (spreaderID != queryID) {
            tuple<TIMETYPE, int, double> rtn = susceptibleTraj->contactWith(spreaderTraj, validSpreadStartTime, pIndex,
                                                                            distThreshold, contactMinDuration, enableCP, independentCP,
                                                                            earliestCST);
            TIMETYPE cst = get<0>(rtn);

            if (cst > 0) {
                if (enableCP) {
                    double prob = get<2>(rtn);
                    aggregatedNonContactProb *= (1 - (1.0 / (qHop + 1)) * prob);
                }

                if(earliestCST == 0 || cst < earliestCST) {
                    earliestCST = cst;
                    earliestSpreaderID = spreaderID;
                }
            }
            checkedQueryPoint += get<1>(rtn);   // statistic
            checkedCandiNum++;
        }
    }
    if (timeCosts.empty()) {
        int len = 2;
        for (int i = 0; i < len; i++)
            timeCosts.emplace_back(0);  // initialize
    }
    timeCosts[0] += checkedQueryPoint / checkedCandiNum;
    timeCosts[1] += checkedCandiNum;

    double aggregatedContactProb = 1 - aggregatedNonContactProb;
    if(earliestSpreaderID.empty() || (enableCP && aggregatedContactProb < contactMinProb)) {
        return nullptr;
    }
    return new ContactEvent(earliestSpreaderID, earliestCST, 0, qHop + 1, aggregatedContactProb);
}

/**
 * Build index for @param initSpreaders
 * Regard each susceptible in @param susceptibleTrajs as query
 * @param susceptible2event: the earliest contact for each susceptible
 * */
string queryFromSusceptible2(const vector<Trip *> &susceptibleTrajs, const vector<Trip *> &initSpreaders, bool ground_truth,
                             IndexPoint *pIndex, const string &searchType, int multiHop, float identicalPoiDist,
                             TIMETYPE contactMinDuration, float contactMinProb, bool independentCP,
                             vector<pair<string, ContactEvent *>> &susceptible2event) {

    vector<pair<Trip *, ContactEvent *>> currSpreaders;   // sorted by the contact start time (record the end of the valid spreading time)
    vector<pair<Trip *, ContactEvent *>> newContacts_k;   // the new contacts returned based on currSpreaders
    for (const auto &traj: initSpreaders) {
        auto *event = new ContactEvent(traj->getId(), traj->getTimeRange().first, traj->getTimeRange().second, 0, 1);
        currSpreaders.emplace_back(traj, event);    // sorted for initial spreader already
    }

    unordered_map<string, Trip *> stillSafeSusceptibles;           // not contacted yet
    for (Trip *traj: susceptibleTrajs)
        stillSafeSusceptibles.insert(make_pair(traj->getId(), traj));


    vector<TIMETYPE> timeCosts;
    TIMETYPE constructionCost = 0, totalSearchCost = 0, timer;

    ItreeIndex *itree = nullptr;
    TIMETYPE earliestIndexTime = EarliestStartTime;

    int LEN = 500;
    char *butter = new char[LEN];
    string msg;
    fflush(stdout);
    int k = 0;

    for (; k < multiHop; k++) {
        printf("\t %d-th hop: ", k + 1);
        fflush(stdout);
        unsigned int cnt = 0, totalQueryNum = stillSafeSusceptibles.size();

        // Create index based on currSpreaders
        if (searchType != "linear") {
            itree = new ItreeIndex;    // interval tree
            timer = millisecond();
            char *rtn = initializeForItree(itree, pIndex, searchType,
                                           earliestIndexTime, contactMinDuration, DeltaT,
                                           currSpreaders);
            constructionCost += (millisecond() - timer);
            msg.append("[Index Construction] ");
            msg.append(rtn);
        }

        timeCosts.clear();
        newContacts_k.clear();
        TIMETYPE earliestSpreaderTime = 0;      // earliestSpreaderTime in newContack_k (found in this round)

        timer = millisecond();

        // Query from susceptibles' view, i.e., the still safe people
        for (const auto &susceptiblePair: stillSafeSusceptibles) {
            string susceptibleID = susceptiblePair.first;
            Trip *susceptibleTraj = susceptiblePair.second;

            ContactEvent *earliestEvent = nullptr;  // a contact found during this round

            if (searchType == "linear") {
                earliestEvent = searchDirectContact_linear(susceptibleTraj, k, currSpreaders, ground_truth,
                                                           pIndex, identicalPoiDist, contactMinDuration,
                                                           contactMinProb, independentCP, timeCosts);
            }
            else if (searchType == "interval" || searchType == "IR") {  // interval tree or IR-tree
                earliestEvent = searchDirectContact_tree(susceptibleTraj, k, itree, searchType, ground_truth, pIndex,
                                                         identicalPoiDist, contactMinDuration, contactMinProb,
                                                         independentCP, timeCosts);
            }

            if (earliestEvent != nullptr) {
                TIMETYPE currCST = earliestEvent->getStartTime();

                // this susceptible is infected this hop, and becomes a new spreader in next round
                newContacts_k.emplace_back(susceptibleTraj, earliestEvent);

                if (earliestSpreaderTime == 0 || earliestSpreaderTime > currCST) {
                    earliestSpreaderTime = currCST;
                }
            }

            // print out progress
            if (++cnt % (totalQueryNum / 10) == 0) {
                printf("%d%%...", (int) (cnt * 100.0 / totalQueryNum));
                fflush(stdout);
            }
        }   // end of scanning all safe susceptibles for k-hop

        printf(" done!\t Found contacts = %lu", newContacts_k.size());

        // Output some statistic information
        TIMETYPE currSearchCost = 0;
        string str;
        if (searchType == "linear") {
            currSearchCost = millisecond() - timer;
            snprintf(butter, LEN,
                     "\t AVG checked points before get contact = %.3f, AVG checked candidates = %.3f, Returned %d-hop contacts = %lu \n",
                     (float) timeCosts[0] / (float) totalQueryNum, (float) timeCosts[1] / (float) totalQueryNum, k + 1,
                     newContacts_k.size());
            str.append(butter);
        }
        else if (searchType == "interval" || searchType == "IR") {
            str = "[TIME COST] " + to_string(k) + ",";
            unsigned int len = 9;
            for (unsigned int i = 0; i < len; i++) {
                if (i < 5) {
                    currSearchCost += timeCosts[i];
                }
                else {
                    timeCosts[i] = timeCosts[i] / totalQueryNum;
                }
                str.append(to_string(timeCosts[i]) + ",");
            }
            str.append(to_string(newContacts_k.size()));

            snprintf(butter, LEN, "\n\t [Search] "
                                  "on I-tree = %lu ms, on R-tree = %lu ms, Filter Pair = %lu ms, Transfer Seq = %lu ms, Check Connect = %lu ms"
                                  "\n\t [TOTAL] Index Construction = %lu ms, Search Total = %lu ms"
                                  "\n\t\t Returned %d-hop contacts = %lu, Remaining susceptibles = %lu, Found contacts total = %lu\n",
                     timeCosts[0], timeCosts[1], timeCosts[2], timeCosts[3], timeCosts[4],
                     constructionCost, currSearchCost, k + 1, newContacts_k.size(), stillSafeSusceptibles.size(),
                     susceptible2event.size());
            str.append(butter);
        }
        totalSearchCost += currSearchCost;

        printf(", %lu ms\n", currSearchCost);

        msg.append(str);

        // Clear the index
        if (itree != nullptr) {
            itree->releaseSpace();
            delete itree;
            itree = nullptr;
        }

        // Reset containers
        butter[0] = '\0';
        currSpreaders.clear();
        currSpreaders.shrink_to_fit();

        if (newContacts_k.empty() || k == multiHop - 1) {
            break;      // no more spreaders or no more susceptibles
        }
        else {
            earliestIndexTime = earliestSpreaderTime;

            // Update the spreader set and susceptible set
            for (pair<Trip *, ContactEvent *> t2e: newContacts_k) {   // to search for next k
                Trip *trip = t2e.first;
                ContactEvent *event = t2e.second;
                if (event->getEndTime() == 0) {
                    event->updateEndTime(trip->getTimeRange().second);
                }
                currSpreaders.emplace_back(trip, event);
                auto itr2 = stillSafeSusceptibles.find(trip->getId());
                if (itr2 != stillSafeSusceptibles.end()) {
                    stillSafeSusceptibles.erase(itr2);
                }
            }
            sort(currSpreaders.begin(), currSpreaders.end(), sortByCST);        // new spreaders are sorted by CST
        }
    }

    msg.append("\n[TIME COST] Contact Search = " + to_string(totalSearchCost) + " ms");
    msg.append("[INFO] total contacts within " + to_string(k + 1) + " hops = " + to_string(susceptible2event.size()));

    return msg;
}

#endif //PCT_INDEX_H
