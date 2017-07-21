cmake_minimum_required(VERSION 3.3)
project(be_analyzer)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lz -O3")

include_directories(
        ${CMAKE_SOURCE_DIR}/fastq-join
)

set(SOURCE_FILES main_v2.cpp be_fastq_join.cpp be_analyzer.cpp be_fastq-lib.cpp)
add_executable(be_analyzer ${SOURCE_FILES})