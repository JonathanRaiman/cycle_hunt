cmake_minimum_required (VERSION 2.8.11)
project (CYCLE_HUNT)

# include the source code under cycle_hunt
add_subdirectory(cycle_hunt)

# connect the Cycle dependency to the executable code.
add_library(Cycle cycle_hunt/Cycle.cpp)
add_executable(cycle_removal cycle_hunt/cycle_removal.cpp)
target_link_libraries(cycle_removal Cycle)

# ensure we are using C++11
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --std=c++11 -O3 -stdlib=libc++ -g" )
