C_FILES := $(shell find ../AntMAN/src -name \*.\*pp -not -name RcppExports.cpp)

R_FLAGS= -m64 -std=gnu++11  -DNDEBUG   -fopenmp   -O3 -g -pipe -Wall \
          -fexceptions -fstack-protector-strong -grecord-gcc-switches \
          -march=native -mtune=native -fasynchronous-unwind-tables  -lboost_program_options -lboost_system

I_FLAGS= -I/usr/local/cuda-10.0/include/  -I../AntMAN/src -I /usr/include/R/
L_FLAGS= -L/usr/local/cuda-10.0/lib64/ -lOpenCL

TEST_FILES=$(shell ls *.cpp)

TARGETS := test_mixtures test_mixtures_debug
all : ${TARGETS}

verbose.o :  ../AntMAN/src/verbose.cpp ${C_FILES} 
	g++ ${R_FLAGS}  ../AntMAN/src/verbose.cpp -DNO_RCPP -I../AntMAN/src -I /usr/include/R/  -larmadillo  -c -o $@ 

test_mixtures : test_mixtures.cpp  verbose.o
	g++ ${R_FLAGS} $^ -DNO_RCPP -I../AntMAN/src -I /usr/include/R/  -larmadillo  -o $@ 

test_mixtures_debug : test_mixtures.cpp  verbose.o
	g++  -DVERBOSE_BINARY=1  ${R_FLAGS}  $^ -DNO_RCPP -I../AntMAN/src -I /usr/include/R/  -larmadillo  -o $@ 

clean :
	rm ${TARGETS} *.o  -f 

