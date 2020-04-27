cmake_minimum_required (VERSION 2.9)

PROJECT(TestAntMAN)

include(CTest)


set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wno-unused-parameter")

find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

SET(ANTMAN_INCLUDES "AntMAN/src/")


function(add_antman_test testname)
	add_executable (${testname} tests_cpp/${testname}.cpp)
	target_include_directories(${testname} PRIVATE ${ANTMAN_INCLUDES} ${Boost_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS})
	target_link_libraries (${testname} PRIVATE Boost::program_options)
	target_link_libraries(${testname} PRIVATE ${ARMADILLO_LIBRARIES})
	target_compile_definitions (${testname} PRIVATE NO_RCPP)
	add_test(NAME ${testname}_100  COMMAND $<TARGET_FILE:${testname}> --niter  100 --burnin   0 --thin  1 --verbose 0)
	add_test(NAME ${testname}_500  COMMAND $<TARGET_FILE:${testname}> --niter  500 --burnin  20 --thin  2 --verbose 0)
	add_test(NAME ${testname}_1000 COMMAND $<TARGET_FILE:${testname}> --niter 1000 --burnin 200 --thin  5 --verbose 0)
	add_test(NAME ${testname}_2000 COMMAND $<TARGET_FILE:${testname}> --niter 2000 --burnin 500 --thin 10 --verbose 0)
	
	add_test(NAME ${testname}_v100  COMMAND $<TARGET_FILE:${testname}> --niter  100 --burnin   0 --thin  1 --verbose 3)
	add_test(NAME ${testname}_v500  COMMAND $<TARGET_FILE:${testname}> --niter  500 --burnin  20 --thin  2 --verbose 3)
	add_test(NAME ${testname}_v1000 COMMAND $<TARGET_FILE:${testname}> --niter 1000 --burnin 200 --thin  5 --verbose 3)
	add_test(NAME ${testname}_v2000 COMMAND $<TARGET_FILE:${testname}> --niter 2000 --burnin 500 --thin 10 --verbose 3)
	
                     
	
endfunction(add_antman_test)

add_antman_test("test_mixtures")