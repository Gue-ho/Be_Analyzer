cmake_minimum_required(VERSION 3.3)
project(be_analyzer)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(be_analyzer ${SOURCE_FILES})