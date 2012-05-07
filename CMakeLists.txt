
cmake_minimum_required(VERSION 2.8)
set(CMAKE_VERBOSE_MAKEFILE OFF)

project(nsga2)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -pedantic -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(NSGA2_SOURCE_FILES
  src/nsga2/global.cpp)

include_directories(src)

add_library(nsga2 SHARED
  ${NSGA2_SOURCE_FILES})

add_executable(test1
  src/nsga2/tests/test1.cpp)

target_link_libraries(test1 nsga2)