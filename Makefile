C_FILES := $(shell find ./AntMAN/src -name \*.\*pp)
R_FILES := $(shell find ./AntMAN/R -name \*.R)

all :  AntMan.install/AntMan/libs/antman.so

%_1.0.tar.gz : ${C_FILES} ${R_FILES}
	R CMD build ./$*
	
%.Rcheck/ : %_1.0.tar.gz
	R CMD check ./$*

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


clean_object :
	rm -rf AntMAN/src/*.o
clean : 
	rm -rf current *~ *.install  *_1.0.tar.gz ..Rcheck  *.Rcheck src/*.o src/*.so 	src/*.rds src/RcppExports.cpp  R/RcppExports.R  man/AM*.Rd 

.PHONY: clean

