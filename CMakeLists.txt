
cmake_minimum_required(VERSION 2.8)
set(CMAKE_VERBOSE_MAKEFILE OFF)

project(nsga2)

# OpenMP
message(STATUS "Checking for OpenMP support")
find_package(OpenMP)
if (OPENMP_FOUND)
  message(STATUS "OpenMP supported. Adding flags ${OpenMP_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  add_definitions(-DUSE_OPENMP)
else (OPENMP_FOUND)
  message(STATUS "OpenMP not supported")
endif (OPENMP_FOUND)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -pedantic -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

set(NSGA2_SOURCE_FILES
  src/nsga2/global.cpp
  src/nsga2/random.cpp
  src/nsga2/NSGA2.cpp
  #src/rand.c # should change to another random generator
  )

set(NSGA2_HEADER_FILES
  src/nsga2/global.h
  src/nsga2/random.h
  src/nsga2/NSGA2.h
  src/nsga2/exception.h
  #src/rand.h
  )

include_directories(src)

add_library(nsga2 SHARED
  ${NSGA2_SOURCE_FILES})

add_library(nsga2-static STATIC
  ${NSGA2_SOURCE_FILES})

# add_executable(test1
#   src/nsga2/tests/test1.cpp)
# target_link_libraries(test1 nsga2)

add_executable(test2
  src/nsga2/tests/test2.cpp)
target_link_libraries(test2 nsga2)

add_executable(test3
  src/nsga2/tests/test3.cpp)
target_link_libraries(test3 nsga2)

add_executable(test4
  src/nsga2/tests/test4.cpp)
target_link_libraries(test4 nsga2)

foreach(HEADER ${NSGA2_HEADER_FILES})
  string(REGEX MATCH "(.*)[/\\]" DIR ${HEADER})
  #message(STATUS "dir is ${DIR}")
  string(REPLACE "src/" "" DIR2 ${DIR})
  #message(STATUS "installing ${PROJECT_SOURCE_DIR}/${HEADER} to include/${DIR2}")
  install(FILES ${PROJECT_SOURCE_DIR}/${HEADER} 
    DESTINATION include/${DIR2} COMPONENT headers)
endforeach(HEADER)

install(TARGETS nsga2 nsga2-static
      RUNTIME DESTINATION bin COMPONENT libraries
      LIBRARY DESTINATION lib COMPONENT libraries
      ARCHIVE DESTINATION lib/static COMPONENT libraries)


