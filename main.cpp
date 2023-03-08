#include "io/LogFile.h"
#include "io/FileWriter.h"
#include "io/Parameters.h"
#include "privacy/Perturbation.h"
#include "tracing/MultiHopContactSearch.h"

const bool FIXSEED = true;
const bool OUTPUT_Contacts = true;
const bool OUTPUT_NoisyTrips = false;

void runContactTracing(LogFile *logger, const string &output_folder, const string &parameterStr, bool ANONYMIZE_DATA, IndexPoint *pIndex,
                       int multiHop, float identicalPoiDist, TIMETYPE contactMinDuration, float contactMinProb, bool independentCP,
                       const string &searchType, const vector<Trip *> &susceptibleTrajs, const vector<Trip *> &initSpreaders);

int main() {

    // get parameters from file "config.properties"
    auto param = new Parameters;
    string msg = param->readParameters("config.properties");
    const string homePath = param->getHomePath();
    const string poiFilename = homePath + "RoadNetworkInfo/AllPOIs.csv";
    const string inputFolder = homePath + "inputData/";     // a folder
    const string outputNoisyFolder = homePath + "outputs_noisyInd/";    // Output of noisy trips
    const string groundTruthFolder = homePath + "ground_truth/";
    const string outputContactFolder = homePath + "outputs_contacts/";  // Output of contacts

    const string datasetType = param->getDatasetName();         // Geolife or Synthetic
    const int dataTotal = param->getObjectNum();
    const int downSample = 0;   // for Geolife
    const int dayLen = 1;

    const string privacyType = param->getPrivacyModel();
    const float *epsilon = param->getPrivacyBudget();

    const bool EnableGeo = (privacyType == "geo");
    const bool EnableGS = (privacyType == "geosemantic" || privacyType == "GS");
    const bool EnableSemantic = (privacyType == "semantic");
    const bool needPoiDistribution = EnableGS || EnableSemantic;
    const bool ANONYMIZE_DATA = EnableGeo || EnableGS || EnableSemantic;

    // Points that are within the threshold will be regarded as identical
    const float identicalPoiDist = param->getDistThreshold();   // unit: km
    const TIMETYPE contactMinDuration = param->getContactTimeThreshold();  // unit: sec
    const float contactMinProb = param->getContactMinProb();

    // PCT: region-based private contact tracing, deterministic, no probability threshold needed
    // PCTi / PCTc: probabilistic contact tracing with independent/correlated contact probability
    const string contactProbType = param->getContactProbType();
    const bool independentCP = (contactProbType == "PCTi");

    const float initContactRatio = param->getInitSpreaderRatio();
    const string searchType = param->getSearchType();
    const int multiHop = param->getHopLength();

    /* ************************************************************************************************************** */
    // prepare log file
    const string parameterStr = LogFile::composeFilename(datasetType, dataTotal, downSample, ANONYMIZE_DATA,
                                                         identicalPoiDist, contactMinDuration, contactMinProb, independentCP,
                                                         multiHop, searchType, initContactRatio, privacyType, epsilon, FIXSEED);
    const string logfilename = homePath + (ANONYMIZE_DATA ? "logs/" : "ground_truth/") + parameterStr + "_T" +
            to_string(millisecond()) + ".txt";
    auto logger = new LogFile(logfilename);
    logger->addContent(msg);


    /* ************************************************************************************************************** */
    // prepare POI list

    auto *pIndex = new IndexPoint;     // Store the spatial points' coordinates
    auto *poiGrid = new StaticGrid;
    {
        msg = FileReader::readPOIs(poiFilename, pIndex, poiGrid);   // POIs with semantic labels
        logger->addContent(msg);

        msg.append("[INFO] # of distinct points in the pIndex = " + to_string(pIndex->getPointNum()) + "\n");
        logger->addContent(msg);
    }

    /* ************************************************************************************************************** */
    // ===== Step 1: Read original data
    logger->addContent(msg);

    vector<Trip *> originalTrajs;
    {
        string filename = inputFolder + datasetType + "_days" + to_string(dayLen);
        if (datasetType == "Geolife")
            filename.append("_DS" + to_string(downSample));     // Geolife has been preprocessed
        filename.append(".csv");

        msg.append("[INFO] Reading data from " + filename + "\n");
        logger->addContent(msg);

        msg = FileReader::readData(filename, dataTotal, pIndex, poiGrid, needPoiDistribution, originalTrajs);

        // !!! data have been sorted by the start time of a trajectory
        logger->addContent(msg);

        if (originalTrajs.empty()) {
            printf("[ERROR] no data from %s", filename.c_str());
            return -1;
        }
    }

    /* ************************************************************************************************************** */
    // ===== Step 2: Sample the initial spreaders as the start
    vector<Trip *> initSpreaders;
    if(initContactRatio > 0 && initContactRatio < 1){
        msg.append("[INFO] Randomly sampling " + to_string((int) (initContactRatio * 100)) + "% init contact trajectories...\n");
        set<string> initSpreaderIDs;

        const int initSpreaderNum = floor(initContactRatio * (float) dataTotal);
        initSpreaders = NoiseGenerator::sampleQuery(initSpreaderNum, originalTrajs, true);

        for(const Trip *trip: initSpreaders)
            initSpreaderIDs.insert(trip->getId());

        msg.append("\t # of initial spreaders = " + to_string(initSpreaders.size()) + "\n");
        logger->addContent(msg);

        // Remove the initial spreaders from the original dataset
        auto itr = originalTrajs.begin();
        while (itr != originalTrajs.end()) {
            if (initSpreaderIDs.find((*itr)->getId()) != initSpreaderIDs.end())
                itr = originalTrajs.erase(itr);
            else
                itr++;
        }
    }
    else if (initContactRatio == 1) {
        initSpreaders = originalTrajs;
    }
    else {
        printf("[ERROR] Init Ratio should be (0, 1]\n");
        return -1;
    }

    /* ************************************************************************************************************** */
    // ===== Step 3: Apply the corresponding indistinguishability model (i.e., G-Ind, S-Ind, GS-Ind)
    vector<Trip *> susceptibleTrajs;

    if (ANONYMIZE_DATA) {
        msg.append("[INFO] Apply " + privacyType + "-Ind to the data ...\n");
        logger->addContent(msg);

        Perturbation::applyPerturbation(originalTrajs, susceptibleTrajs, epsilon, pIndex, poiGrid,
                                        EnableGeo, EnableSemantic, EnableGS, FIXSEED);

        msg.append("\t # of anonymized trajectories = " + to_string(susceptibleTrajs.size()) + "\n");
        logger->addContent(msg);

        if(OUTPUT_NoisyTrips) {
            msg = "\n------------------------------\n";
            msg.append("[EXTRA] Output Noisy Trips\n");
            msg.append("------------------------------\n");
            logger->addContent(msg);
            string outputFile = outputNoisyFolder + parameterStr + ".csv";
            FileWriter::outputNoisyTrips(outputFile, pIndex, susceptibleTrajs);
        }
    }
    else {
        susceptibleTrajs = originalTrajs;       // no perturbation for data, used for generating ground truth
    }

    /* ************************************************************************************************************** */
    // ===== Step 4: Contact Tracing (Exact/Private/Probabilistic)

    printf("-------------------------------------------------------------\n");

    const string output_folder = !ANONYMIZE_DATA ? groundTruthFolder : outputContactFolder;

    runContactTracing(logger, output_folder, parameterStr, ANONYMIZE_DATA, pIndex,
                      multiHop, identicalPoiDist, contactMinDuration, contactMinProb, independentCP,
                      searchType, susceptibleTrajs, initSpreaders);

    printf("\n\n-------------------------------------------------------------\n");
    printf("Program Complete!\n");

    delete logger;              logger = nullptr;
    delete param;               param = nullptr;
    originalTrajs.clear();      originalTrajs.shrink_to_fit();
    initSpreaders.clear();      initSpreaders.shrink_to_fit();
    susceptibleTrajs.clear();   susceptibleTrajs.shrink_to_fit();
    return 0;
}

void runContactTracing(LogFile *logger, const string &output_folder, const string &parameterStr, bool ANONYMIZE_DATA, IndexPoint *pIndex,
                       int multiHop, float identicalPoiDist, TIMETYPE contactMinDuration, float contactMinProb, bool independentCP,
                       const string &searchType, const vector<Trip *> &susceptibleTrajs, const vector<Trip *> &initSpreaders) {

    vector<pair<string, ContactEvent *>> susceptible2event;

    if (!initSpreaders.empty()) {
        string msg;
        char *buffer = new char[200];
        snprintf(buffer, 200,
                 "[INFO] Thresholds: search-%s, %d-hop [max dist %.3f km], [min duration %lu s], [min prob %.2f]",
                 searchType.c_str(), multiHop, identicalPoiDist, contactMinDuration, contactMinProb);
        msg.append(buffer);
        if (contactMinProb > 0) {
            msg.append(independentCP ? "-i" : "-c");
        }
        msg.append("\n");
        logger->addContent(msg);

        string rtn = queryFromSusceptible2(susceptibleTrajs, initSpreaders, !ANONYMIZE_DATA, pIndex,
                                            searchType, multiHop, identicalPoiDist,
                                            contactMinDuration, contactMinProb, independentCP,
                                            susceptible2event);
        logger->addContent(rtn, false);
    }


    /* ************************************************************************************************************** */
    // ===== Step 5: Output Contacts to the file / screen
    if(OUTPUT_Contacts && !susceptible2event.empty()){

        string outputFile = output_folder + parameterStr
                            + "_MH" + to_string(multiHop) + "_CD" + float2str(identicalPoiDist) + "_CT" + to_string(contactMinDuration);
        if (ANONYMIZE_DATA) {
            outputFile.append("_PCT");
            if (contactMinProb > 0) {
                outputFile.append(independentCP ? "i" : "c");
                outputFile.append(float2str(contactMinProb));
            }
        }
        else {
            outputFile.append("_CT");
        }
        outputFile.append("_" + searchType);
        outputFile.append(".csv");

        FileWriter::outputContacts(outputFile, multiHop, (contactMinProb > 0 && ANONYMIZE_DATA), susceptible2event);
    }
    susceptible2event.clear();
}