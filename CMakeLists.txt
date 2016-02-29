cmake_minimum_required (VERSION 2.8.11)
project (CYCLE_HUNT)


add_library(Cycle Cycle.cpp)

add_executable(cycle_removal cycle_removal.cpp)

target_link_libraries(cycle_removal Cycle)

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --std=c++11 -O3 -stdlib=libc++ -g" )
