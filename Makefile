C_FILES := $(shell find src -name *.*pp)
R_FILES := $(shell find R -name *.R)

all : 
	make AntMan.install/AntMan/libs/antman.so

deps :
	echo "To be defined."

AntMan.install/AntMan/libs/antman.so : ${C_FILES} ${R_FILES}
	rm -rf AntMan_1.0.tar.gz AntMan.Rcheck AntMan.install
	R -e  "Rcpp::compileAttributes(pkgdir = \"./\" , verbose=TRUE);"
	R CMD check . || ( cat  AntMan.Rcheck/00install.out && false )
	R CMD build . 
	mkdir -p AntMan.install
	R CMD INSTALL  -l AntMan.install AntMan_1.0.tar.gz
	R -e  "library(devtools) ; document() ; devtools::build_manual();"



clean : 
	rm -rf current *~ *.install  *_1.0.tar.gz *.Rcheck src/*.o src/*.so 	src/*.rds src/RcppExports.cpp  R/RcppExports.R 

.PHONY: clean

