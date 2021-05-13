FIND_PACKAGE( Boost 1.33.0 COMPONENTS unit_test_framework  REQUIRED)

FILE(GLOB_RECURSE TEST_FILES *Test.cpp)
MESSAGE(STATUS "prepare tests: ${TEST_FILES}")

ADD_LIBRARY(testutils testutils.cpp)
target_compile_definitions (testutils PRIVATE NO_RCPP)
target_include_directories(testutils PRIVATE ${ANTMAN_INCLUDES})
target_link_libraries(testutils  PRIVATE ${ARMADILLO_LIBRARIES})
#target_include_directories(testutils PRIVATE ${ANTMAN_INCLUDES} ${Boost_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS})
#
	
FOREACH(TEST_FILE  ${TEST_FILES})
	GET_FILENAME_COMPONENT(TEST_NAME ${TEST_FILE} NAME_WE)
	MESSAGE(STATUS "Add ${TEST_NAME} - ${TEST_FILE}")
	ADD_EXECUTABLE(${TEST_NAME} ${TEST_FILE})
	target_include_directories(${TEST_NAME} PRIVATE ${ANTMAN_INCLUDES} ${ARMADILLO_INCLUDE_DIRS})
	target_link_libraries(${TEST_NAME} PRIVATE testutils)
	target_compile_definitions (${TEST_NAME} PRIVATE NO_RCPP)
	ADD_TEST(NAME ${TEST_NAME} COMMAND ${TEST_NAME}  --log_level=all --report_level=detailed  --detect_memory_leak=1)
ENDFOREACH()