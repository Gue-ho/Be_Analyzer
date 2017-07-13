cmake_minimum_required(VERSION 3.3)
project(be_analyzer)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(
        ${CMAKE_SOURCE_DIR}/fastq-join
)

set(SOURCE_FILES main_v2.cpp fastq-join/fastq-lib.cpp)
add_executable(be_analyzer ${SOURCE_FILES})