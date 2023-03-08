# PCT
Paper Submission for PVLDB_v16 2023 (ID: 1356)

## Title: Efficient Privacy-Preserving Contact Tracing with Semantic-Aware Indistinguishability Guarantee

### Project Structure

    CMakeLists.txt                                -- the version of cmake might be revised accordingly (currently is 3.15)
    config.properties                             -- the program will read parameters here
    cmake-build/                                  -- If needed, clean some files, then re-compile and re-build in your device
        Testing/                                  -- all testing data and some possible outputs
            RoadNetworkInfo/AllPOIs.csv           -- the semantic POIs and their categories (please keep it private)
            inputData/                            -- the input Geolife data after preprocessing (i.e., partition and temporal alignment)
                                                      due to the size limit, only a small sample is given here for testing
            outputs_contacts/                     -- the discovered contacts will be outputted here
    ...                 
    main.cpp                                      -- the entry of the program
    spatial/                                      -- some basic geometry used in this program
    io/                                           -- read/write files
    tracing/                                      -- the functions for privacy-preserving contact tracing
    privacy/                                      -- indistinguishability mechanisms are implemented here
    ...                                           -- other relevant classes and headers (not detailed here)

### Perform the PCT

To compile this project using cmake and c++ compiler, run this script:

    ./buildPCT.sh

If a permission denied error happens, please try to modify the permission of this script file by `chmod 700 buildPCT.sh`

To run the program compiled before, execute the command in the folder of `cmake-build`:

    ./PCT

This command will execute the program in `./main.cpp`.

Again, the test configuration is in `./config.properties` (associated with detailed explanation for each parameter).

### Environment
Tested with CentOS Linux (gcc version 9.2.0) and macOS Monterey (with Apple M1 Pro Chip, Apple clang version 14.0.0).

### Dependency
In particular, [Boost](https://www.boost.org/users/history/version_1_78_0.html) is needed (tested with 1.78.0 version) to support the Laplacian geo-indistinguishability mechanism (details can be found in the [original paper](https://dl.acm.org/doi/10.1145/2508859.2516735)). If the Boost lib has been installed in your device, please update the path in `./CMakeLists.txt` as well.

Thank you!
