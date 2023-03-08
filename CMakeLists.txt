cmake_minimum_required(VERSION 3.22)
project(PCT)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_BUILD_TYPE Release)

set(Boost_INCLUDE_DIR /usr/local/include)
set(Boost_LIBRARY_DIR /usr/local/bin)
find_package(Boost 1.78.0)
include_directories(${Boost_INCLUDE_DIR})

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find boost!")
endif()

add_executable(PCT main.cpp commonheader.h utils.h
        io/FileReader.cpp io/FileReader.h io/Parameters.cpp io/Parameters.h io/LogFile.cpp io/LogFile.h io/FileWriter.cpp io/FileWriter.h
        spatial/GeoPoint.cpp spatial/GeoPoint.h spatial/Point.cpp spatial/Point.h spatial/Trip.cpp spatial/Trip.h spatial/IndexPoint.cpp spatial/IndexPoint.h
        spatial/StaticGrid.cpp spatial/StaticGrid.h spatial/Category.cpp spatial/Category.h
        spatial/Rectangle.cpp spatial/Rectangle.h
        index/rtree/RtreeIndex.cpp index/rtree/RtreeIndex.h
        index/rtree/RtreeNode.cpp index/rtree/RtreeNode.h index/rtree/RtreeElement.h
        index/intervaltree/ItreeIndex.cpp index/intervaltree/ItreeIndex.h
        index/intervaltree/ItreeNode.cpp index/intervaltree/ItreeNode.h
        index/intervaltree/Interval.cpp index/intervaltree/Interval.h
        privacy/NoiseGenerator.cpp privacy/NoiseGenerator.h privacy/Perturbation.cpp privacy/Perturbation.h
        tracing/SearchByIndex.cpp tracing/SearchByIndex.h tracing/integral.h tracing/functions.h
        tracing/MultiHopContactSearch.h)
