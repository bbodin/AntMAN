C_FILES := $(shell find src -name \*.\*pp)
T_FILES := $(shell find src/tests -name \*.\*pp)
R_FILES := $(shell find R -name \*.R)

R_PATH    := "/usr/include/R"
RLIBS_PATH := "/home/toky/R/x86_64-redhat-linux-gnu-library/3.5/"

all :  AntMan.install/AntMan/libs/antman.so

test-% : src/tests/test-%.cpp
	g++ -m64 -std=gnu++11 -Isrc  -I ${R_PATH}  -DNDEBUG  -I"${RLIBS_PATH}/Rcpp/include" -I"${RLIBS_PATH}/RcppArmadillo/include" -I/usr/local/include  -fopenmp  -fpic  -O2 -g -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wp,-D_GLIBCXX_ASSERTIONS -fexceptions -fstack-protector-strong -grecord-gcc-switches -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 -specs=/usr/lib/rpm/redhat/redhat-annobin-cc1 -m64 -mtune=generic -fasynchronous-unwind-tables -fstack-clash-protection -fcf-protection -c src/tests/test-PriorNegativeBinomial.cpp -o  test1.o
	 g++ -m64 -std=gnu++11  -L/usr/lib64/R/lib -Wl,-z,relro -Wl,-z,now -specs=/usr/lib/rpm/redhat/redhat-hardened-ld test1.o src/utils.cpp -fopenmp -L/usr/lib64/R/lib -lRlapack -L/usr/lib64/R/lib -lRblas -lgfortran -lm -lquadmath -L/usr/lib64/R/lib -lR -o $@



deps :
	echo "To be defined."

AntMan.install/AntMan/libs/antman.so : ${C_FILES} ${R_FILES}
	rm -rf AntMan_1.0.tar.gz AntMan.Rcheck AntMan.install man/AM*.Rd 
	R -e  "Rcpp::compileAttributes(pkgdir = \"./\" , verbose=TRUE);"
	R CMD check . || ( cat  ..Rcheck/00install.out && false )
	R CMD build . 
	mkdir -p AntMan.install
	R CMD INSTALL  -l AntMan.install AntMan_1.0.tar.gz
	R -e  "library(devtools) ; document() ; devtools::build_manual();"



clean : 
	rm -rf current *~ *.install  *_1.0.tar.gz ..Rcheck  *.Rcheck src/*.o src/*.so 	src/*.rds src/RcppExports.cpp  R/RcppExports.R  man/AM*.Rd 

.PHONY: clean

