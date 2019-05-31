C_FILES := $(shell find ./AntMAN/src -name \*.\*pp -not -name RcppExports.cpp)
R_FILES := $(shell find ./AntMAN/R ./AntMAN/tests -name \*.R -not -name RcppExports.R)

infos :
	@echo "C_FILES=${C_FILES}"
	@echo "R_FILES=${R_FILES}"
	
all :  AntMan.install/AntMan/libs/antman.so

%/src/RcppExports.cpp  %/R/RcppExports.R : % ${C_FILES} ${R_FILES}
	R -e  "Rcpp::compileAttributes(pkgdir = \"$*\" , verbose=TRUE);"

%_1.0.tar.gz : ${C_FILES} ${R_FILES} %/src/RcppExports.cpp  %/R/RcppExports.R
	rm -rf AntMAN/src/*.o ./AntMAN/src/*.so 
	R CMD build ./$*
	
%.Rcheck/ : ${C_FILES} ${R_FILES} %/src/RcppExports.cpp  %/R/RcppExports.R
	R CMD check ./$*
	
%.Rinstall/ : %_1.0.tar.gz
	mkdir -p $@
	R CMD INSTALL  -l $@ $*_1.0.tar.gz
	
%_1.0.pdf : %
	R -e  "library(devtools) ; pkgbuild::compile_dll(\"$*\");  document(\"$*\"); devtools::build_manual(\"$*\"); "
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
clean : 
	rm -rf current *~ *.install  *_1.0.tar.gz *.Rcheck ./AntMAN/src/*.o ./AntMAN/src/*.so 	./AntMAN/src/*.rds ./AntMAN/src/RcppExports.cpp  ./AntMAN/R/RcppExports.R  ./AntMAN/man/AM*.Rd 

.PHONY: clean

