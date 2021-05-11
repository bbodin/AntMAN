cmake_minimum_required (VERSION 2.9)

PROJECT(TestAntMAN)

include(CTest)


set(CMAKE_CXX_FLAGS "--std=c++11 -O3 -Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "--std=c++11 -g -Wall -Wextra -Wno-unused-parameter")

find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

SET(ANTMAN_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/AntMAN/src/")

ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/tests_cpp/ "")

	