C_FILES := $(shell find ./AntMAN/src -name \*.\*pp -not -name RcppExports.cpp)
R_FILES := $(shell find ./AntMAN/R ./AntMAN/tests -name \*.R -not -name RcppExports.R)


all :  AntMAN.Rinstall/   AntMAN_1.0.pdf

download : 
	cp ~/Dropbox/AntMan/AntManAPI.R AntMAN/R/AntManAPI.R

upload : 
	cp AntMAN/R/AntManAPI.R ~/Dropbox/AntMan/AntManAPI.R 

infos :
	@echo "C_FILES=${C_FILES}"
	@echo "R_FILES=${R_FILES}"


%/src/RcppExports.cpp  %/R/RcppExports.R : % ${C_FILES} ${R_FILES}
	R -e  "Rcpp::compileAttributes(pkgdir = \"$*\" , verbose=TRUE);"

%_1.0.tar.gz : ${C_FILES} ${R_FILES} %/src/RcppExports.cpp  %/R/RcppExports.R  AntMAN_1.0.pdf
	rm -rf AntMAN/src/*.o ./AntMAN/src/*.so 
	R CMD build ./$*

%.Rcheck/ : ${C_FILES} ${R_FILES} %/src/RcppExports.cpp  %/R/RcppExports.R  AntMAN_1.0.pdf
	R CMD check ./$*

%.Rinstall/ : %_1.0.tar.gz 
	mkdir -p $@
	R CMD INSTALL  -l $@ $*_1.0.tar.gz

%_1.0.pdf : %
	R -e  "library(devtools) ; pkgbuild::compile_dll(\"$*\");  document(\"$*\"); devtools::build_manual(\"$*\"); "

deps :
	echo "To be defined."

clean : 
	rm -rf current *~ *.Rinstall *_1.0.pdf  *_1.0.tar.gz *.Rcheck ./AntMAN/src/*.o ./AntMAN/src/*.so 	./AntMAN/src/*.rds ./AntMAN/src/RcppExports.cpp  ./AntMAN/R/RcppExports.R  ./AntMAN/man/AM*.Rd 

.PHONY: clean

