C_FILES := $(shell find ../AntMAN/src -name \*.\*pp -not -name RcppExports.cpp)

R_FLAGS= -m64 -std=gnu++11  -DNDEBUG   -fopenmp   -O2 -g -pipe -Wall \
          -fexceptions -fstack-protector-strong -grecord-gcc-switches \
         -mtune=generic -fasynchronous-unwind-tables  -lboost_program_options -lboost_system

TEST_FILES=$(shell ls *.cpp)

all : testAntMAN testAntMANd

testAntMAN : ${TEST_FILES} ../AntMAN/src/verbose.cpp ${C_FILES} 
	g++ ${R_FLAGS}  ${TEST_FILES} ../AntMAN/src/verbose.cpp -DNO_RCPP -I../AntMAN/src -I /usr/include/R/  -larmadillo  -o $@ 

testAntMANd : ${TEST_FILES} ../AntMAN/src/verbose.cpp ${C_FILES} 
	g++ -DVERBOSE_BINARY=1 ${R_FLAGS}  ${TEST_FILES} ../AntMAN/src/verbose.cpp -DNO_RCPP -I../AntMAN/src -I /usr/include/R/  -larmadillo  -o $@ 

clean :
	rm  testAntMAN -f 

