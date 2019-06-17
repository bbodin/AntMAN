
R_FLAGS= -m64 -std=gnu++11  -DNDEBUG   -fopenmp   -O2 -g -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 \
         -Wp,-D_GLIBCXX_ASSERTIONS -fexceptions -fstack-protector-strong -grecord-gcc-switches \
         -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 -specs=/usr/lib/rpm/redhat/redhat-annobin-cc1 -m64 \
         -mtune=generic -fasynchronous-unwind-tables -fstack-clash-protection -fcf-protection 

TEST_FILES=$(shell ls *.cpp)

all : testAntMAN
	./testAntMAN

testAntMAN : ${TEST_FILES} ../src/verbose.cpp
	g++ ${R_FLAGS} $^ -DNO_RCPP -I../src -I /usr/include/R/  -larmadillo  -o $@ 

	
clean :
	rm  testAntMAN 
	